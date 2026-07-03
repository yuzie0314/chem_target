"""Held-out generalisation evaluation — honest out-of-tuning accuracy.

Motivation
----------
Every conditional rule in ``utils/target_predictor.py`` was tuned against the
curated 20-per-class benchmark (``data/benchmark/curated/compounds.csv``).  The
headline "core-11 190/220 = 86.4%" is therefore a *resubstitution* number
(measured on the tuning set) — it upper-bounds, not estimates, real-world
accuracy.

The held-out set
----------------
``data/benchmark/holdout/compounds.csv`` is a PURE held-out set: it is the
broader ``limit`` download with every tuning compound removed by BOTH ChEMBL id
AND canonical SMILES (zero structural overlap with the curated tuning set).
Rebuild it deterministically from the two frozen inputs with ``--rebuild``.

Scoring the identical model on this set gives an honest out-of-tuning estimate.
It reuses ``run_benchmark.run_predictions`` so FG detection, scoring and class
canonicalisation are byte-identical to the benchmark pipeline.

Usage
-----
    python eval_holdout.py                # score curated + holdout, print the gap
    python eval_holdout.py --rebuild      # re-derive holdout from curated+limit
    python eval_holdout.py --rerun        # force re-prediction on both sets

Fully reproducible offline: curated + limit inputs are frozen in git.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))

from run_benchmark import run_predictions  # noqa: E402

# The 11 mechanistic classes behind the headline 190/220 figure.
CORE11: frozenset[str] = frozenset({
    "mTOR", "tubulin", "carbonic anhydrase", "HDAC", "CYP450", "COX",
    "serine protease", "adenosine receptor", "GPCR", "kinase", "nuclear receptor",
})

_BENCH = _ROOT / "data" / "benchmark"
_OUT = _ROOT / "output" / "benchmark"
_CURATED_CSV = _BENCH / "curated" / "compounds.csv"
_LIMIT_CSV = _BENCH / "limit" / "compounds.csv"
_HOLDOUT_CSV = _BENCH / "holdout" / "compounds.csv"


def _canonical_smiles(smiles: str) -> str | None:
    """Return RDKit canonical SMILES, or None if unparseable/empty."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return Chem.MolToSmiles(mol) if mol else None


def derive_holdout() -> Path:
    """Build the pure held-out set = limit minus curated (id AND canonical SMILES).

    Removes any limit compound that shares a ChEMBL id OR a canonical SMILES with
    the curated tuning set, plus intra-set structural duplicates and unparseable
    rows.  Writes data/benchmark/holdout/compounds.csv and returns its path.
    """
    def load(p: Path) -> list[dict]:
        with open(p, encoding="utf-8") as fh:
            return list(csv.DictReader(fh))

    if not _CURATED_CSV.exists() or not _LIMIT_CSV.exists():
        sys.exit(f"Need both frozen inputs:\n  {_CURATED_CSV}\n  {_LIMIT_CSV}")

    cur, lim = load(_CURATED_CSV), load(_LIMIT_CSV)
    cur_ids = {r["chembl_id"] for r in cur}
    cur_smiles = {c for c in (_canonical_smiles(r["smiles"]) for r in cur) if c}

    kept: list[dict] = []
    seen: set[str] = set()
    d_id = d_smi = d_bad = 0
    for r in lim:
        cs = _canonical_smiles(r["smiles"])
        if not cs:
            d_bad += 1
            continue
        if r["chembl_id"] in cur_ids:
            d_id += 1
            continue
        if cs in cur_smiles:
            d_smi += 1
            continue
        if cs in seen:
            continue
        seen.add(cs)
        kept.append(r)

    _HOLDOUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(_HOLDOUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=lim[0].keys())
        w.writeheader()
        w.writerows(kept)
    print(f"[holdout] built {_HOLDOUT_CSV}: kept {len(kept)}/{len(lim)} "
          f"(dropped id={d_id}, smiles={d_smi}, invalid={d_bad})")
    return _HOLDOUT_CSV


def _ensure_results(mode: str, compounds: Path, rerun: bool) -> Path:
    """Return the results CSV for *mode*, (re)generating it from *compounds*."""
    results = _OUT / f"{mode}_results.csv"
    if not compounds.exists():
        sys.exit(f"Frozen benchmark set missing: {compounds}")
    if rerun or not results.exists():
        _OUT.mkdir(parents=True, exist_ok=True)
        print(f"[{mode}] running predictions on {compounds} ...")
        results = run_predictions(compounds, _OUT, mode)
    return results


def _core11(df: pd.DataFrame) -> pd.DataFrame:
    """Rows whose ground-truth class is one of the core-11 mechanistic classes."""
    return df[df["true_target_class"].isin(CORE11)].copy()


def _acc(df: pd.DataFrame) -> tuple[int, int]:
    """(hits, n) for Top-1 on the given rows."""
    hits = (df["top1_hit"].astype(str) == "1").sum()
    return int(hits), len(df)


def main() -> None:
    """Compute and print the tuning-vs-held-out core-11 Top-1 gap."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rerun", action="store_true",
                    help="Force re-prediction even if results CSVs exist.")
    ap.add_argument("--rebuild", action="store_true",
                    help="Re-derive the holdout set from curated+limit first.")
    args = ap.parse_args()

    if args.rebuild or not _HOLDOUT_CSV.exists():
        derive_holdout()

    cur = pd.read_csv(_ensure_results("curated", _CURATED_CSV, args.rerun), dtype=str)
    hold = pd.read_csv(_ensure_results("holdout", _HOLDOUT_CSV, args.rerun), dtype=str)

    # Sanity: the held-out set must not overlap the tuning set by ChEMBL id.
    overlap = set(cur["chembl_id"]) & set(hold["chembl_id"])
    assert not overlap, f"held-out leaked {len(overlap)} tuning compounds: {overlap}"

    cur_c, hold_c = _core11(cur), _core11(hold)

    print("\n" + "=" * 62)
    print("  Core-11 Top-1 accuracy — tuning set vs PURE held-out")
    print("=" * 62)
    for label, d in [
        ("TUNING set (curated 20/class)", cur_c),
        ("PURE HELD-OUT (limit − curated, zero overlap)  <<<", hold_c),
    ]:
        h, n = _acc(d)
        print(f"  {label:51s} {h:4d}/{n:<4d} = {h / n * 100:5.1f}%")

    h_t, n_t = _acc(cur_c)
    h_h, n_h = _acc(hold_c)
    gap = h_t / n_t * 100 - h_h / n_h * 100
    print("-" * 62)
    print(f"  Overfitting gap (tuning − held-out): {gap:+.1f} percentage points")
    print("=" * 62)

    print("\n  Per-class held-out Top-1 (sorted worst → best):")
    pc: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for _, row in hold_c.iterrows():
        cls = row["true_target_class"]
        pc[cls][1] += 1
        pc[cls][0] += str(row["top1_hit"]) == "1"
    for cls in sorted(pc, key=lambda c: pc[c][0] / pc[c][1]):
        h, n = pc[cls]
        print(f"    {cls:22s} {h:3d}/{n:<3d} = {h / n * 100:5.1f}%")
    print()


if __name__ == "__main__":
    main()
