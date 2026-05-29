"""chem_target — benchmark evaluation pipeline.

Two modes
---------
curated   20 compounds per target class (high-confidence ground truth).
          Compounds are fetched from ChEMBL by pChEMBL value (≥ 6.0 = 1 µM),
          SDF downloaded when available, SMILES used as fallback.

limit     Large-scale batch (N compounds, multi-class, broader coverage).
          Samples across all mapped ChEMBL targets, no per-class cap.

Usage examples
--------------
  # Step 1 — download compound data (writes data/benchmark/)
  python run_benchmark.py download --mode curated
  python run_benchmark.py download --mode limit --limit 2000

  # Step 2 — run predictions (reads data/benchmark/, writes output/benchmark/)
  python run_benchmark.py run --mode curated
  python run_benchmark.py run --mode limit

  # Step 3 — generate summary report
  python run_benchmark.py report --mode curated
  python run_benchmark.py report --mode limit

  # All-in-one (download + run + report)
  python run_benchmark.py all --mode curated
  python run_benchmark.py all --mode limit --limit 2000

Outputs
-------
  data/benchmark/{mode}/compounds.csv         ChEMBL metadata + SMILES + true class
  data/benchmark/{mode}/sdf/*.sdf             3D SDF files (when available)
  output/benchmark/{mode}_results.csv         Per-compound prediction vs ground truth
  output/benchmark/{mode}_summary.csv         Per-class accuracy table
  output/benchmark/{mode}_report.txt          Plain-text summary (publication-ready)

IMPORTANT — Bias disclosure
---------------------------
The residue scoring database (db/fg_residue_table.csv) is derived from
BioLiP 2.0, which aggregates binding events from the RCSB Protein Data Bank
(PDB).  Compounds that have co-crystal structures deposited in the PDB will
have artificially inflated residue scores because their structural context is
directly reflected in the scoring table.  ChEMBL-derived test compounds *may*
overlap with PDB entries.  All benchmark results should be interpreted with
this potential circularity in mind.  See README.md § "Validation bias" for
the full disclosure.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

import requests

# ── Project root on sys.path ───────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))

from utils.target_predictor import predict as _predict

# ── ChEMBL API ─────────────────────────────────────────────────────────────────
_CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_REQUEST_DELAY = 0.4   # seconds between API calls (ChEMBL rate-limit courtesy)

# ── Target class → ChEMBL target ID mapping ────────────────────────────────────
#
# Maps each of our broad "target class" vocabulary words to one or more
# specific ChEMBL single-protein targets.  Compounds from ALL listed targets
# are pooled and the top-N by pChEMBL value are kept.
#
# To add a new class: append an entry here; no other changes needed.
#
TARGET_CLASS_MAP: dict[str, list[str]] = {
    "COX":                  ["CHEMBL230", "CHEMBL220"],           # COX-2, COX-1
    "kinase":               ["CHEMBL203", "CHEMBL5145", "CHEMBL301",
                             "CHEMBL2842"],                       # EGFR, BRAF, CDK2, mTOR(kinase domain)
    "GPCR":                 ["CHEMBL210", "CHEMBL217", "CHEMBL218",
                             "CHEMBL226", "CHEMBL251"],           # β2-AR, D2R, 5-HT2A, A1R, A2AR
    "serine protease":      ["CHEMBL204", "CHEMBL209"],           # thrombin, trypsin
    "cysteine protease":    ["CHEMBL4523", "CHEMBL3227"],         # cathepsin B, cathepsin L
    "nuclear receptor":     ["CHEMBL1871", "CHEMBL206",
                             "CHEMBL2034", "CHEMBL3151"],         # AR, ERα, GR, PR
    "MAO":                  ["CHEMBL2366517", "CHEMBL2828"],      # MAO-A, MAO-B
    "HDAC":                 ["CHEMBL325", "CHEMBL1865",
                             "CHEMBL3192"],                       # HDAC1, HDAC6, HDAC8
    "adenosine receptor":   ["CHEMBL226", "CHEMBL251"],           # A1R, A2AR
    "carbonic anhydrase":   ["CHEMBL205", "CHEMBL3729"],          # CA-II, CA-IX
    "CYP450":               ["CHEMBL340", "CHEMBL1952",
                             "CHEMBL3356"],                       # CYP3A4, CYP2D6, CYP2C9
    "PDE":                  ["CHEMBL1827", "CHEMBL3769"],         # PDE5A, PDE4B
    "mTOR":                 ["CHEMBL2842"],                       # mTOR
    "tubulin":              ["CHEMBL379"],                        # tubulin α1A
    "VKORC1":               ["CHEMBL1953583"],                    # VKORC1
    "topoisomerase":        ["CHEMBL3952", "CHEMBL3191"],         # Topo I, Topo II
    "ribosome":             ["CHEMBL612558"],                     # 50S ribosomal (bacterial)
    "xanthine oxidase":     ["CHEMBL1916"],                       # XO (xanthine dehydrogenase)
    "COMT":                 ["CHEMBL4203"],                       # COMT
}

# Compounds per class (curated mode)
_CURATED_PER_CLASS = 20
# pChEMBL cutoff (6.0 = 1 µM, 7.0 = 100 nM)
_PCHEMBL_MIN = 6.0


# ── Helpers: ChEMBL API ────────────────────────────────────────────────────────

def _get_json(url: str, params: dict | None = None, retries: int = 3) -> dict:
    """GET a JSON endpoint from ChEMBL with retry logic."""
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=20)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as exc:
            if attempt == retries - 1:
                raise
            print(f"  [retry {attempt+1}] {exc}")
            time.sleep(2 ** attempt)
    return {}


def _fetch_activities(
    target_id: str,
    pchembl_min: float = _PCHEMBL_MIN,
    limit: int = 200,
) -> list[dict]:
    """Return activity records for a ChEMBL target sorted by pChEMBL desc."""
    url    = f"{_CHEMBL_BASE}/activity.json"
    params = {
        "target_chembl_id":    target_id,
        "pchembl_value__gte":  str(pchembl_min),
        "assay_type":          "B",          # binding assays only
        "limit":               limit,
        "order_by":            "-pchembl_value",
        "format":              "json",
    }
    data = _get_json(url, params)
    return data.get("activities", [])


def _fetch_sdf(chembl_id: str, sdf_dir: Path) -> Optional[Path]:
    """Download SDF for a ChEMBL compound; return path or None on failure."""
    sdf_path = sdf_dir / f"{chembl_id}.sdf"
    if sdf_path.exists():
        return sdf_path   # cached
    url = f"{_CHEMBL_BASE}/molecule/{chembl_id}.sdf"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 200 and r.text.strip():
            sdf_path.write_text(r.text, encoding="utf-8")
            return sdf_path
    except requests.RequestException:
        pass
    return None


def _smiles_from_sdf(sdf_path: Path) -> Optional[str]:
    """Extract canonical SMILES from a ChEMBL SDF property block."""
    text = sdf_path.read_text(encoding="utf-8", errors="replace")
    # ChEMBL SDF embeds SMILES in a > <chembl_smiles> or <smiles> property tag
    for prop in ("chembl_smiles", "smiles", "canonical_smiles"):
        tag = f"> <{prop}>"
        idx = text.lower().find(tag.lower())
        if idx != -1:
            after = text[idx + len(tag):].strip()
            smi   = after.split("\n")[0].strip()
            if smi and smi != "$$$$":
                return smi
    # Fallback: second line of SDF is often a SMILES in some formats
    lines = [l.rstrip() for l in text.splitlines()]
    for i, line in enumerate(lines):
        if line.strip() == "M  END" and i + 1 < len(lines):
            # Try next non-blank lines for SMILES properties
            pass
    return None


def _smiles_from_chembl(chembl_id: str) -> Optional[str]:
    """Fetch canonical SMILES for a ChEMBL compound via molecule endpoint."""
    url  = f"{_CHEMBL_BASE}/molecule/{chembl_id}.json"
    data = _get_json(url)
    mol  = data.get("molecule_structures") or {}
    return mol.get("canonical_smiles") or mol.get("molfile")


# ── Phase 1: Download ──────────────────────────────────────────────────────────

def download_curated(
    out_dir: Path,
    per_class: int = _CURATED_PER_CLASS,
    pchembl_min: float = _PCHEMBL_MIN,
) -> Path:
    """Download curated test compounds: top-N per target class.

    Returns path to compounds.csv.
    """
    sdf_dir = out_dir / "sdf"
    sdf_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    seen_chembl: set[str] = set()   # deduplicate across targets in same class

    for cls, target_ids in TARGET_CLASS_MAP.items():
        print(f"\n[{cls}] fetching up to {per_class} compounds "
              f"from {len(target_ids)} ChEMBL target(s)...")
        cls_rows: list[dict] = []
        cls_seen: set[str]   = set()

        for tid in target_ids:
            if len(cls_rows) >= per_class:
                break
            activities = _fetch_activities(tid, pchembl_min=pchembl_min, limit=200)
            time.sleep(_REQUEST_DELAY)

            for act in activities:
                if len(cls_rows) >= per_class:
                    break
                cid = act.get("molecule_chembl_id", "")
                if not cid or cid in cls_seen or cid in seen_chembl:
                    continue

                # Prefer SMILES from ChEMBL molecule endpoint
                mol_smiles = act.get("canonical_smiles") or ""
                if not mol_smiles:
                    mol_smiles = _smiles_from_chembl(cid) or ""
                    time.sleep(_REQUEST_DELAY)

                if not mol_smiles:
                    continue   # skip if no SMILES at all

                # Try SDF download
                sdf_path = _fetch_sdf(cid, sdf_dir)
                time.sleep(_REQUEST_DELAY)
                sdf_smiles = _smiles_from_sdf(sdf_path) if sdf_path else None
                has_sdf    = sdf_path is not None and sdf_path.exists()

                # Use SDF-derived SMILES if available; ChEMBL SMILES as fallback
                smiles_used = sdf_smiles or mol_smiles

                cls_seen.add(cid)
                cls_rows.append({
                    "compound_name":    act.get("molecule_pref_name") or cid,
                    "chembl_id":        cid,
                    "smiles":           smiles_used,
                    "true_target_class": cls,
                    "chembl_target_id": tid,
                    "pchembl_value":    act.get("pchembl_value", ""),
                    "standard_type":    act.get("standard_type", ""),
                    "has_sdf":          str(has_sdf),
                    "sdf_file":         sdf_path.name if has_sdf else "",
                })

        seen_chembl.update(cls_seen)
        rows.extend(cls_rows)
        print(f"  → {len(cls_rows)} compounds collected for '{cls}'")

    csv_path = out_dir / "compounds.csv"
    _write_csv(rows, csv_path)
    print(f"\nCurated compounds saved: {csv_path}  ({len(rows)} total)")
    return csv_path


def download_limit(
    out_dir: Path,
    limit: int = 2000,
    pchembl_min: float = _PCHEMBL_MIN,
) -> Path:
    """Download limit-test compounds: broad multi-class sample up to `limit`.

    Returns path to compounds.csv.
    """
    sdf_dir = out_dir / "sdf"
    sdf_dir.mkdir(parents=True, exist_ok=True)

    rows:          list[dict] = []
    seen_chembl:   set[str]   = set()
    per_class_cap = max(10, limit // len(TARGET_CLASS_MAP))

    for cls, target_ids in TARGET_CLASS_MAP.items():
        if len(rows) >= limit:
            break
        print(f"[{cls}] collecting up to {per_class_cap} ...")
        cls_count = 0

        for tid in target_ids:
            if cls_count >= per_class_cap or len(rows) >= limit:
                break
            activities = _fetch_activities(tid, pchembl_min=pchembl_min, limit=300)
            time.sleep(_REQUEST_DELAY)

            for act in activities:
                if cls_count >= per_class_cap or len(rows) >= limit:
                    break
                cid = act.get("molecule_chembl_id", "")
                if not cid or cid in seen_chembl:
                    continue

                mol_smiles = act.get("canonical_smiles") or _smiles_from_chembl(cid) or ""
                time.sleep(_REQUEST_DELAY)
                if not mol_smiles:
                    continue

                # SDF: download only if not already cached (limit test: skip slow downloads)
                sdf_path = (sdf_dir / f"{cid}.sdf") if (sdf_dir / f"{cid}.sdf").exists() else None
                if sdf_path is None:
                    # Quick attempt only (no SDF download in limit mode to save time)
                    sdf_path = None
                has_sdf   = sdf_path is not None and sdf_path.exists()
                sdf_smi   = _smiles_from_sdf(sdf_path) if has_sdf else None

                seen_chembl.add(cid)
                cls_count += 1
                rows.append({
                    "compound_name":     act.get("molecule_pref_name") or cid,
                    "chembl_id":         cid,
                    "smiles":            sdf_smi or mol_smiles,
                    "true_target_class": cls,
                    "chembl_target_id":  tid,
                    "pchembl_value":     act.get("pchembl_value", ""),
                    "standard_type":     act.get("standard_type", ""),
                    "has_sdf":           str(has_sdf),
                    "sdf_file":          sdf_path.name if has_sdf else "",
                })

    csv_path = out_dir / "compounds.csv"
    _write_csv(rows, csv_path)
    print(f"\nLimit-test compounds saved: {csv_path}  ({len(rows)} total)")
    return csv_path


# ── Phase 2: Run predictions ───────────────────────────────────────────────────

def run_predictions(
    compounds_csv: Path,
    results_dir: Path,
    mode: str,
    top_k: int = 10,
) -> Path:
    """Run target prediction on all compounds in compounds.csv.

    Saves per-compound results and returns path to results CSV.
    """
    results_dir.mkdir(parents=True, exist_ok=True)

    compounds = list(_read_csv(compounds_csv))
    n_total   = len(compounds)
    print(f"\nRunning predictions on {n_total} compounds ...")

    rows = []
    for i, row in enumerate(compounds, 1):
        name   = row.get("compound_name", "")
        smiles = row.get("smiles", "").strip()
        true_c = row.get("true_target_class", "")

        if not smiles:
            rows.append(_empty_result(row, "no SMILES"))
            continue

        try:
            pred = _predict(smiles, top_residues=top_k)
        except Exception as exc:
            rows.append(_empty_result(row, str(exc)[:120]))
            continue

        tc  = pred.get("target_class_votes")
        fgs = pred.get("fgs_detected", [])

        if tc is None or tc.empty:
            top1, top2, top3      = "", "", ""
            s1, s2, s3            = "", "", ""
            v1, v2, v3            = "", "", ""
            true_rank             = -1
            top1_hit = top3_hit   = False
        else:
            classes  = list(tc["target_class"])
            scores   = list(tc["score"])
            votes    = list(tc["votes"])
            top1     = classes[0] if len(classes) > 0 else ""
            top2     = classes[1] if len(classes) > 1 else ""
            top3     = classes[2] if len(classes) > 2 else ""
            s1       = round(scores[0], 3) if scores else ""
            s2       = round(scores[1], 3) if len(scores) > 1 else ""
            s3       = round(scores[2], 3) if len(scores) > 2 else ""
            v1       = int(votes[0]) if votes else ""
            v2       = int(votes[1]) if len(votes) > 1 else ""
            v3       = int(votes[2]) if len(votes) > 2 else ""
            # Rank of the true target class (1-based; 0 = not found)
            true_rank = next(
                (r + 1 for r, c in enumerate(classes) if c == true_c), 0
            )
            top1_hit = true_rank == 1
            top3_hit = 1 <= true_rank <= 3

        result = {
            **row,
            "n_fgs":         len(fgs),
            "fgs_detected":  "; ".join(fgs),
            "top1":          top1,
            "top2":          top2,
            "top3":          top3,
            "score_1":       s1,
            "score_2":       s2,
            "score_3":       s3,
            "votes_1":       v1,
            "votes_2":       v2,
            "votes_3":       v3,
            "true_rank":     true_rank,
            "top1_hit":      int(top1_hit),
            "top3_hit":      int(top3_hit),
            "mrr":           round(1 / true_rank, 4) if true_rank > 0 else 0.0,
            "warning":       pred.get("warning") or "",
        }
        rows.append(result)

        if i % 20 == 0 or i == n_total:
            print(f"  {i}/{n_total} done ...")

    results_csv = results_dir / f"{mode}_results.csv"
    _write_csv(rows, results_csv)
    print(f"Results saved: {results_csv}")
    return results_csv


def _empty_result(row: dict, warning: str) -> dict:
    return {
        **row,
        "n_fgs": 0, "fgs_detected": "", "top1": "", "top2": "", "top3": "",
        "score_1": "", "score_2": "", "score_3": "",
        "votes_1": "", "votes_2": "", "votes_3": "",
        "true_rank": -1, "top1_hit": 0, "top3_hit": 0, "mrr": 0.0,
        "warning": warning,
    }


# ── Phase 3: Report ───────────────────────────────────────────────────────────

def generate_report(results_csv: Path, results_dir: Path, mode: str) -> None:
    """Compute accuracy metrics and save per-class summary + plain-text report."""
    import pandas as pd

    df = pd.read_csv(results_csv, dtype=str)
    # Cast numeric columns
    for col in ("top1_hit", "top3_hit", "mrr", "n_fgs", "true_rank", "pchembl_value"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    valid = df[df["warning"].fillna("").eq("") | df["n_fgs"].gt(0)]
    failed = df[df["warning"].fillna("").ne("") & df["n_fgs"].fillna(0).le(0)]

    # ── Per-class table ──────────────────────────────────────────────────────
    cls_rows = []
    for cls, grp in df.groupby("true_target_class"):
        n       = len(grp)
        n_valid = grp["n_fgs"].gt(0).sum()
        top1    = grp["top1_hit"].sum()
        top3    = grp["top3_hit"].sum()
        mrr     = grp["mrr"].mean()
        cls_rows.append({
            "target_class":    cls,
            "n_compounds":     n,
            "n_valid_fgs":     int(n_valid),
            "top1_acc":        round(top1 / n, 3) if n else 0,
            "top3_acc":        round(top3 / n, 3) if n else 0,
            "mean_mrr":        round(mrr, 3),
            "top1_count":      int(top1),
            "top3_count":      int(top3),
        })
    cls_df = pd.DataFrame(cls_rows).sort_values("top3_acc", ascending=False)
    summary_csv = results_dir / f"{mode}_summary.csv"
    cls_df.to_csv(summary_csv, index=False)

    # ── Overall stats ────────────────────────────────────────────────────────
    n_all      = len(df)
    n_valid_fg = int(df["n_fgs"].gt(0).sum())
    top1_all   = df["top1_hit"].sum()
    top3_all   = df["top3_hit"].sum()
    mrr_all    = df["mrr"].mean()

    # ── Plain-text report ────────────────────────────────────────────────────
    lines = [
        "=" * 70,
        f"  chem_target Benchmark Report  ({mode.upper()} mode)",
        "=" * 70,
        "",
        "IMPORTANT — Validation bias disclosure",
        "-" * 40,
        "  The residue scoring layer is trained on BioLiP 2.0 data sourced",
        "  from the RCSB PDB.  ChEMBL test compounds with known crystal",
        "  structures in the PDB may have inflated scores due to structural",
        "  circularity.  These results should be interpreted as an",
        "  upper-bound estimate of prediction performance.",
        "",
        "OVERALL RESULTS",
        "-" * 40,
        f"  Total compounds:         {n_all}",
        f"  With ≥1 FG detected:    {n_valid_fg}  ({100*n_valid_fg/n_all:.1f}%)",
        f"  No FG detected (failed): {n_all - n_valid_fg}",
        f"  Top-1 accuracy:          {top1_all}/{n_all}  ({100*top1_all/n_all:.1f}%)",
        f"  Top-3 accuracy:          {top3_all}/{n_all}  ({100*top3_all/n_all:.1f}%)",
        f"  Mean Reciprocal Rank:    {mrr_all:.3f}",
        "",
        "PER-CLASS BREAKDOWN",
        "-" * 40,
    ]

    header = f"  {'Target class':<22} {'N':>4} {'Top-1':>7} {'Top-3':>7} {'MRR':>6}"
    lines.append(header)
    lines.append("  " + "-" * 50)
    for _, r in cls_df.iterrows():
        lines.append(
            f"  {r['target_class']:<22} "
            f"{int(r['n_compounds']):>4} "
            f"{r['top1_acc']:>7.1%} "
            f"{r['top3_acc']:>7.1%} "
            f"{r['mean_mrr']:>6.3f}"
        )

    lines += [
        "",
        "MOST COMMON FAILURE MODES",
        "-" * 40,
    ]
    if len(failed) > 0:
        lines.append(f"  No FG detected: {len(failed)} compounds")
        warn_counts = Counter(failed["warning"].fillna("unknown"))
        for msg, cnt in warn_counts.most_common(5):
            lines.append(f"    [{cnt}x] {msg[:60]}")

    # FG coverage analysis
    no_fg = df[df["n_fgs"].fillna(0).le(0)]
    if len(no_fg) > 0:
        lines.append(f"\n  True classes with most FG-failures:")
        cls_fail = no_fg["true_target_class"].value_counts().head(5)
        for cls, cnt in cls_fail.items():
            lines.append(f"    {cls}: {cnt} compounds with 0 FGs")

    lines += [
        "",
        "FILES",
        "-" * 40,
        f"  Per-compound results:  {results_csv}",
        f"  Per-class summary:     {summary_csv}",
        "",
        "=" * 70,
    ]

    report_path = results_dir / f"{mode}_report.txt"
    report_text = "\n".join(lines)
    report_path.write_text(report_text, encoding="utf-8")
    print(report_text)
    print(f"\nReport saved: {report_path}")
    print(f"Summary CSV:  {summary_csv}")


# ── I/O helpers ───────────────────────────────────────────────────────────────

def _write_csv(rows: list[dict], path: Path) -> None:
    """Write a list of dicts to CSV, inferring fieldnames from first row."""
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _read_csv(path: Path) -> list[dict]:
    """Read a CSV file into a list of dicts."""
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# ── CLI ───────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_benchmark",
        description="chem_target benchmark evaluation pipeline",
    )
    sub = p.add_subparsers(dest="phase", metavar="<phase>")
    sub.required = True

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--mode", choices=["curated", "limit"], default="curated",
        help="Benchmark mode (default: curated)",
    )

    # download
    dl = sub.add_parser("download", parents=[common],
                        help="Download compounds from ChEMBL")
    dl.add_argument("--per-class", type=int, default=_CURATED_PER_CLASS,
                    help=f"Compounds per class, curated mode (default: {_CURATED_PER_CLASS})")
    dl.add_argument("--limit", type=int, default=2000,
                    help="Total compound limit, limit mode (default: 2000)")
    dl.add_argument("--pchembl-min", type=float, default=_PCHEMBL_MIN,
                    help=f"Minimum pChEMBL value (default: {_PCHEMBL_MIN})")

    # run
    sub.add_parser("run", parents=[common],
                   help="Run target predictions on downloaded compounds")

    # report
    sub.add_parser("report", parents=[common],
                   help="Generate accuracy report from prediction results")

    # all (download + run + report)
    al = sub.add_parser("all", parents=[common],
                        help="Download, run, and report in one step")
    al.add_argument("--per-class", type=int, default=_CURATED_PER_CLASS)
    al.add_argument("--limit", type=int, default=2000)
    al.add_argument("--pchembl-min", type=float, default=_PCHEMBL_MIN)

    return p


def main() -> None:
    """Entry point for the benchmark pipeline."""
    parser  = _build_parser()
    args    = parser.parse_args()
    mode    = args.mode

    bench_dir   = Path("data") / "benchmark" / mode
    results_dir = Path("output") / "benchmark"
    bench_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    compounds_csv = bench_dir / "compounds.csv"
    results_csv   = results_dir / f"{mode}_results.csv"

    if args.phase in ("download", "all"):
        print(f"=== PHASE: DOWNLOAD ({mode}) ===")
        if mode == "curated":
            compounds_csv = download_curated(
                bench_dir,
                per_class=getattr(args, "per_class", _CURATED_PER_CLASS),
                pchembl_min=getattr(args, "pchembl_min", _PCHEMBL_MIN),
            )
        else:
            compounds_csv = download_limit(
                bench_dir,
                limit=getattr(args, "limit", 2000),
                pchembl_min=getattr(args, "pchembl_min", _PCHEMBL_MIN),
            )

    if args.phase in ("run", "all"):
        if not compounds_csv.exists():
            print(f"Error: {compounds_csv} not found. Run 'download' first.",
                  file=sys.stderr)
            sys.exit(1)
        print(f"\n=== PHASE: RUN ({mode}) ===")
        results_csv = run_predictions(compounds_csv, results_dir, mode)

    if args.phase in ("report", "all"):
        if not results_csv.exists():
            print(f"Error: {results_csv} not found. Run 'run' first.",
                  file=sys.stderr)
            sys.exit(1)
        print(f"\n=== PHASE: REPORT ({mode}) ===")
        generate_report(results_csv, results_dir, mode)


if __name__ == "__main__":
    main()
