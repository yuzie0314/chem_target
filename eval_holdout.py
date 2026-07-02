"""Held-out generalisation evaluation — honest out-of-tuning accuracy.

Motivation
----------
Every conditional rule in ``utils/target_predictor.py`` was tuned against the
curated 20-per-class benchmark (``data/benchmark/curated/compounds.csv``).  The
headline "core-11 190/220 = 86.4%" is therefore a *resubstitution* number
(measured on the tuning set) — it upper-bounds, not estimates, real-world
accuracy.

This script quantifies the overfitting gap by scoring the same model against the
broader ``limit`` set (~105 compounds/class) with the curated tuning compounds
removed, i.e. a genuine held-out set.  It reuses ``run_benchmark.run_predictions``
so the FG detection, scoring and class-canonicalisation are byte-identical to the
benchmark pipeline.

Usage
-----
    python eval_holdout.py                # regenerate predictions if missing
    python eval_holdout.py --rerun        # force re-prediction on both sets

Both input sets are frozen in git (curated + limit compounds.csv), so this is
fully reproducible offline.
"""

from __future__ import annotations

import argparse
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


def _ensure_results(mode: str, rerun: bool) -> Path:
    """Return the results CSV for *mode*, (re)generating it from the frozen set."""
    compounds = _BENCH / mode / "compounds.csv"
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
    args = ap.parse_args()

    cur = pd.read_csv(_ensure_results("curated", args.rerun), dtype=str)
    lim = pd.read_csv(_ensure_results("limit", args.rerun), dtype=str)

    tuning_ids = set(cur["chembl_id"])
    cur_c = _core11(cur)
    lim_c = _core11(lim)
    held = lim_c[~lim_c["chembl_id"].isin(tuning_ids)]

    print("\n" + "=" * 62)
    print("  Core-11 Top-1 accuracy — tuning set vs held-out")
    print("=" * 62)
    for label, d in [
        ("TUNING set (curated 20/class)", cur_c),
        ("limit set — all core-11 (incl. tuning overlap)", lim_c),
        ("HELD-OUT (limit minus tuning compounds)  <<<", held),
    ]:
        h, n = _acc(d)
        print(f"  {label:47s} {h:4d}/{n:<4d} = {h / n * 100:5.1f}%")

    h_t, n_t = _acc(cur_c)
    h_h, n_h = _acc(held)
    gap = h_t / n_t * 100 - h_h / n_h * 100
    print("-" * 62)
    print(f"  Overfitting gap (tuning − held-out): {gap:+.1f} percentage points")
    print("=" * 62)

    print("\n  Per-class held-out Top-1 (sorted worst → best):")
    pc: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for _, row in held.iterrows():
        cls = row["true_target_class"]
        pc[cls][1] += 1
        pc[cls][0] += str(row["top1_hit"]) == "1"
    for cls in sorted(pc, key=lambda c: pc[c][0] / pc[c][1]):
        h, n = pc[cls]
        print(f"    {cls:22s} {h:3d}/{n:<3d} = {h / n * 100:5.1f}%")
    print()


if __name__ == "__main__":
    main()
