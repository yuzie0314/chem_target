"""SwissTargetPrediction (STP) comparison benchmark.

Queries SwissTargetPrediction for every compound in the curated benchmark
set and computes Top-1 / Top-3 / MRR against the same ground-truth labels
used for the chem_target internal benchmark.

Usage
-----
  # Run STP on curated set (reads data/benchmark/curated/compounds.csv)
  python run_stp_comparison.py

  # Limit to first N compounds (for quick testing)
  python run_stp_comparison.py --limit 30

  # Skip STP queries; regenerate report from cached stp_raw.csv
  python run_stp_comparison.py --report-only

  # Compare side-by-side (requires both benchmarks already run)
  python run_stp_comparison.py --compare

Outputs
-------
  output/benchmark/stp_raw.csv          Raw STP predictions (one row per
                                         compound×target, all 100 targets)
  output/benchmark/stp_results.csv      Per-compound Top-1/Top-3 results
  output/benchmark/stp_summary.csv      Per-class accuracy table
  output/benchmark/stp_report.txt       Plain-text narrative report
  output/benchmark/comparison_report.txt Side-by-side chem_target vs STP

IMPORTANT — Terms of use
------------------------
SwissTargetPrediction is a free web service for academic use.
(https://www.swisstargetprediction.ch/)
Please use rate-limiting (default 2 s between requests) and do not
distribute the raw STP data commercially.

IMPORTANT — Evaluation methodology
------------------------------------
STP is a fingerprint-similarity tool (FP2 + FP4) trained on ChEMBL.
Our benchmark ChEMBL compounds may appear in STP's training data,
so STP Top-1 / Top-3 numbers here reflect *in-distribution* performance.
chem_target's bias comes from BioLiP/PDB structural data.
Both biases are disclosed in the full report.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Optional

import requests

# ── Force UTF-8 output (Windows cp950 workaround) ─────────────────────────────
# Only rewrap if stdout has a .buffer (raw stream) — avoids double-wrapping
# when stdout is already replaced (e.g. by a log-file wrapper).
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Paths ──────────────────────────────────────────────────────────────────────
_ROOT          = Path(__file__).parent
_CURATED_CSV   = _ROOT / "data" / "benchmark" / "curated" / "compounds.csv"
_OUT_DIR       = _ROOT / "output" / "benchmark"
_STP_RAW_CSV   = _OUT_DIR / "stp_raw.csv"
_STP_RESULTS   = _OUT_DIR / "stp_results.csv"
_STP_SUMMARY   = _OUT_DIR / "stp_summary.csv"
_STP_REPORT    = _OUT_DIR / "stp_report.txt"
_CMP_REPORT    = _OUT_DIR / "comparison_report.txt"
_CHEM_RESULTS  = _OUT_DIR / "curated_results.csv"   # chem_target results

# ── STP API settings ───────────────────────────────────────────────────────────
_STP_PREDICT_URL = "https://www.swisstargetprediction.ch/predict.php"
_STP_ORGANISM    = "Homo_sapiens"
_REQUEST_DELAY   = 2.0   # seconds between STP submissions (rate-limit courtesy)

# ── ChEMBL → benchmark class mapping ──────────────────────────────────────────
# Derived from TARGET_CLASS_MAP in run_benchmark.py.
# When a ChEMBL ID appears in multiple classes (e.g. A1R in GPCR and adenosine
# receptor), ALL matching classes are stored so either can score a hit.

_TARGET_CLASS_MAP: dict[str, list[str]] = {
    "COX":               ["CHEMBL230", "CHEMBL220", "CHEMBL221"],
    "kinase":            ["CHEMBL203", "CHEMBL5145", "CHEMBL301", "CHEMBL2842"],
    "GPCR":              ["CHEMBL210", "CHEMBL217", "CHEMBL218", "CHEMBL226", "CHEMBL251"],
    "serine protease":   ["CHEMBL204", "CHEMBL209"],
    "cysteine protease": ["CHEMBL4523", "CHEMBL3227"],
    "nuclear receptor":  ["CHEMBL1871", "CHEMBL206", "CHEMBL2034", "CHEMBL3151"],
    "MAO":               ["CHEMBL2366517", "CHEMBL2828"],
    "HDAC":              ["CHEMBL325", "CHEMBL1865", "CHEMBL3192"],
    "adenosine receptor":["CHEMBL226", "CHEMBL251"],
    "carbonic anhydrase":["CHEMBL205", "CHEMBL3729"],
    "CYP450":            ["CHEMBL340", "CHEMBL1952", "CHEMBL3356"],
    "PDE":               ["CHEMBL1827", "CHEMBL3769"],
    "mTOR":              ["CHEMBL2842"],
    "tubulin":           ["CHEMBL379"],
    "VKORC1":            ["CHEMBL1953583"],
    "topoisomerase":     ["CHEMBL3952", "CHEMBL3191"],
    "ribosome":          ["CHEMBL612558"],
    "xanthine oxidase":  ["CHEMBL1916"],
    "COMT":              ["CHEMBL4203"],
}

# Invert: ChEMBL ID → list of benchmark classes (one ID may map to several)
_CHEMBL_TO_CLASSES: dict[str, list[str]] = defaultdict(list)
for _cls, _ids in _TARGET_CLASS_MAP.items():
    for _cid in _ids:
        _CHEMBL_TO_CLASSES[_cid].append(_cls)

# ── STP target-class → benchmark class fallback ───────────────────────────────
# Used when STP's predicted ChEMBL ID is not in our target map.
# Maps STP's broad enzyme-class labels to our benchmark classes.
_STP_CLASS_FALLBACK: dict[str, str | None] = {
    "Cytochrome P450":                    "CYP450",
    "Family A G protein-coupled receptor":"GPCR",
    "Family C G protein-coupled receptor":"GPCR",
    "Kinase":                             "kinase",
    "Nuclear receptor":                   "nuclear receptor",
    "Protease":                           "serine protease",   # most common in our benchmark
    "Lyase":                              "carbonic anhydrase",# CA are lyases
    "Hydrolase":                          "serine protease",
    "Eraser":                             "HDAC",
    "Transferase":                        "COMT",
    "Isomerase":                          "topoisomerase",
    "Other cytosolic protein":            "tubulin",
    "Oxidoreductase":                     "COX",               # most common drug target in class
    "Ligand-gated ion channel":           None,
    "Electrochemical transporter":        None,
    "Enzyme":                             None,
    "Phosphatase":                        None,
    "Secreted protein":                   None,
    "Unclassified protein":               None,
    "Fatty acid binding protein family":  None,
    "Ligase":                             None,
}

# ── STP query + parse ─────────────────────────────────────────────────────────

def _query_stp(smiles: str, session: requests.Session,
               retries: int = 3) -> list[dict]:
    """Submit SMILES to STP and return list of prediction dicts.

    Each dict has keys: target, common_name, uniprot_id, chembl_id,
    stp_class, probability (float).
    Returns empty list on failure.
    """
    payload = {"smiles": smiles, "organism": _STP_ORGANISM, "ioi": "2"}

    for attempt in range(retries):
        try:
            r = session.post(_STP_PREDICT_URL, data=payload, timeout=120)
            if r.status_code != 200:
                raise ValueError(f"HTTP {r.status_code}")

            # Extract redirect URL containing job ID
            redirect_matches = re.findall(r'location\.replace\("(.*?)"\)', r.text)
            if not redirect_matches:
                raise ValueError("No redirect URL found in predict.php response")

            result_url = redirect_matches[0]
            r2 = session.get(result_url, timeout=60)
            if r2.status_code != 200:
                raise ValueError(f"result.php HTTP {r2.status_code}")

            return _parse_stp_table(r2.text)

        except Exception as exc:
            if attempt == retries - 1:
                print(f"    [STP error] {exc}", flush=True)
                return []
            wait = 2 ** attempt * 3
            print(f"    [STP retry {attempt+1}] {exc} — waiting {wait}s ...", flush=True)
            time.sleep(wait)

    return []


def _parse_stp_table(html: str) -> list[dict]:
    """Parse the STP result HTML table into a list of prediction dicts."""
    rows = re.findall(r'<tr[^>]*>(.*?)</tr>', html, re.DOTALL | re.IGNORECASE)
    results: list[dict] = []

    for row in rows:
        cells = re.findall(r'<t[hd][^>]*>(.*?)</t[hd]>', row, re.DOTALL | re.IGNORECASE)
        if len(cells) < 6:
            continue
        cleaned = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]
        # Skip header
        if cleaned[0].lower() == 'target':
            continue
        # Parse probability (may be "1.0" or "0.1234...")
        try:
            prob = float(cleaned[5]) if cleaned[5] else 0.0
        except ValueError:
            prob = 0.0

        results.append({
            "target":      cleaned[0],
            "common_name": cleaned[1] if len(cleaned) > 1 else "",
            "uniprot_id":  cleaned[2] if len(cleaned) > 2 else "",
            "chembl_id":   cleaned[3] if len(cleaned) > 3 else "",
            "stp_class":   cleaned[4] if len(cleaned) > 4 else "",
            "probability": prob,
        })

    # Sort by descending probability (STP should already be sorted, but enforce it)
    results.sort(key=lambda x: x["probability"], reverse=True)
    return results


# ── Class assignment ──────────────────────────────────────────────────────────

def _assign_benchmark_classes(stp_pred: dict) -> list[str]:
    """Return the benchmark class(es) for one STP prediction row.

    Priority:
      1. ChEMBL ID exact match to _CHEMBL_TO_CLASSES  (most precise)
      2. STP target class name in _STP_CLASS_FALLBACK  (broad fallback)
      3. [] — no mapping found
    """
    cid = stp_pred.get("chembl_id", "").strip()
    if cid and cid in _CHEMBL_TO_CLASSES:
        return _CHEMBL_TO_CLASSES[cid]

    stpc = stp_pred.get("stp_class", "")
    fallback = _STP_CLASS_FALLBACK.get(stpc)
    return [fallback] if fallback else []


# ── Per-compound evaluation ───────────────────────────────────────────────────

def _evaluate_compound(true_class: str,
                        stp_preds: list[dict]) -> dict:
    """Determine Top-1 / Top-3 hit and rank for one compound.

    Iterates through STP predictions (sorted by probability descending).
    Each prediction is assigned to benchmark class(es); a hit is recorded
    when the true class appears among assigned classes.
    Returns dict with: top1 (bool), top3 (bool), rank (int|None), mrr (float).
    """
    true_rank: int | None = None

    for rank, pred in enumerate(stp_preds, start=1):
        assigned = _assign_benchmark_classes(pred)
        if true_class in assigned:
            true_rank = rank
            break

    return {
        "top1": true_rank == 1,
        "top3": true_rank is not None and true_rank <= 3,
        "rank": true_rank,
        "mrr":  1.0 / true_rank if true_rank else 0.0,
    }


# ── Phase 1: Query STP ────────────────────────────────────────────────────────

def run_stp_queries(compounds: list[dict], delay: float = _REQUEST_DELAY
                    ) -> list[dict]:
    """Query STP for every compound and return annotated result rows.

    Appends stp_preds (list), top1, top3, rank, mrr to each compound dict.
    Raw predictions are also written to stp_raw.csv as they arrive (crash-safe).
    """
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (chem_target benchmark; academic use)",
        "Referer":    "https://www.swisstargetprediction.ch/",
    })

    results: list[dict] = []

    raw_path = _STP_RAW_CSV
    raw_file  = open(raw_path, "w", newline="", encoding="utf-8")
    raw_writer = csv.writer(raw_file)
    raw_writer.writerow(["compound_name", "true_target_class", "smiles",
                          "rank", "target", "common_name", "uniprot_id",
                          "chembl_id", "stp_class", "probability"])

    total = len(compounds)
    for i, cmpd in enumerate(compounds, start=1):
        name       = cmpd.get("compound_name", cmpd.get("name", f"cpd{i}"))
        smiles     = cmpd.get("smiles", "")
        true_class = cmpd.get("true_target_class", cmpd.get("true_class", ""))

        print(f"  [{i}/{total}] {name[:40]} ({true_class}) ...", flush=True)

        preds = _query_stp(smiles, session)

        # Write raw predictions
        for rank, pred in enumerate(preds, start=1):
            raw_writer.writerow([
                name, true_class, smiles,
                rank,
                pred["target"], pred["common_name"],
                pred["uniprot_id"], pred["chembl_id"],
                pred["stp_class"], pred["probability"],
            ])
        raw_file.flush()

        eval_result = _evaluate_compound(true_class, preds)

        results.append({
            "compound_name": name,
            "smiles":        smiles,
            "true_class":    true_class,
            "stp_top1":      eval_result["top1"],
            "stp_top3":      eval_result["top3"],
            "stp_rank":      eval_result["rank"],
            "stp_mrr":       eval_result["mrr"],
            "stp_n_preds":   len(preds),
            "stp_best_class": (
                ", ".join(_assign_benchmark_classes(preds[0])) if preds else ""
            ),
        })

        if i < total:
            time.sleep(delay)

    raw_file.close()
    print(f"\nRaw STP predictions saved: {raw_path}")
    return results


# ── Phase 2: Report generation ────────────────────────────────────────────────

def generate_report(results: list[dict],
                    out_results: Path  = _STP_RESULTS,
                    out_summary: Path  = _STP_SUMMARY,
                    out_report:  Path  = _STP_REPORT) -> dict:
    """Write per-compound CSV, per-class summary CSV, and plain-text report.

    Returns per-class stats dict for use in comparison report.
    """
    _OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Per-compound results ────────────────────────────────────────────────
    with open(out_results, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "compound_name", "smiles", "true_class",
            "stp_top1", "stp_top3", "stp_rank", "stp_mrr",
            "stp_n_preds", "stp_best_class",
        ])
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})

    # ── Per-class summary ───────────────────────────────────────────────────
    class_stats: dict[str, dict] = {}
    for r in results:
        cls = r.get("true_class") or r.get("true_target_class") or "unknown"
        if cls not in class_stats:
            class_stats[cls] = {"n": 0, "top1": 0, "top3": 0, "mrr_sum": 0.0}
        s = class_stats[cls]
        s["n"]       += 1
        s["top1"]    += int(r.get("stp_top1", False))
        s["top3"]    += int(r.get("stp_top3", False))
        s["mrr_sum"] += float(r.get("stp_mrr", 0.0))

    with open(out_summary, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["target_class", "n", "top1", "top1_pct",
                    "top3", "top3_pct", "mrr"])
        for cls, s in sorted(class_stats.items(), key=lambda x: -x[1]["top3"]):
            n = s["n"]
            w.writerow([
                cls, n,
                s["top1"], f'{s["top1"]/n*100:.1f}%',
                s["top3"], f'{s["top3"]/n*100:.1f}%',
                f'{s["mrr_sum"]/n:.3f}',
            ])

    # ── Plain-text report ───────────────────────────────────────────────────
    total       = len(results)
    total_top1  = sum(int(r.get("stp_top1", False)) for r in results)
    total_top3  = sum(int(r.get("stp_top3", False)) for r in results)
    total_mrr   = sum(float(r.get("stp_mrr", 0.0)) for r in results) / max(total, 1)
    no_preds    = sum(1 for r in results if not r.get("stp_n_preds", 0))

    lines: list[str] = [
        "=" * 70,
        "  SwissTargetPrediction (STP) Benchmark Report  (CURATED mode)",
        "=" * 70,
        "",
        "EVALUATION NOTES",
        "-" * 40,
        "  STP uses fingerprint similarity (FP2 + FP4 against ChEMBL).",
        "  Our ChEMBL test compounds may overlap with STP training data —",
        "  results represent an *in-distribution* upper bound for STP.",
        "  chem_target scores via BioLiP residue interactions; its",
        "  structural bias is disclosed separately (curated_report.txt).",
        "  Evaluation: benchmark class assigned via ChEMBL ID match",
        "  (primary) or STP target-class label (fallback).",
        "",
        "OVERALL RESULTS",
        "-" * 40,
        f"  Total compounds:         {total}",
        f"  No STP results:          {no_preds}",
        f"  Top-1 accuracy:          {total_top1}/{total}  ({total_top1/total*100:.1f}%)",
        f"  Top-3 accuracy:          {total_top3}/{total}  ({total_top3/total*100:.1f}%)",
        f"  Mean Reciprocal Rank:    {total_mrr:.3f}",
        "",
        "PER-CLASS BREAKDOWN",
        "-" * 40,
        f"  {'Target class':<24}  {'N':>4}  {'Top-1':>6}  {'Top-3':>6}  {'MRR':>6}",
        "  " + "-" * 54,
    ]
    for cls, s in sorted(class_stats.items(), key=lambda x: -x[1]["mrr_sum"]/x[1]["n"]):
        n = s["n"]
        lines.append(
            f"  {cls:<24}  {n:>4}  {s['top1']/n*100:>5.1f}%  "
            f"{s['top3']/n*100:>5.1f}%  {s['mrr_sum']/n:>6.3f}"
        )
    lines += [
        "",
        "FILES",
        "-" * 40,
        f"  Raw STP predictions:   {_STP_RAW_CSV}",
        f"  Per-compound results:  {out_results}",
        f"  Per-class summary:     {out_summary}",
        "",
        "=" * 70,
        "",
    ]

    report_text = "\n".join(lines)
    out_report.write_text(report_text, encoding="utf-8")
    print(report_text)
    print(f"Report saved: {out_report}")

    return class_stats


# ── Phase 3: Side-by-side comparison ─────────────────────────────────────────

def generate_comparison(stp_stats: dict[str, dict],
                        chem_results_csv: Path = _CHEM_RESULTS,
                        out_report: Path = _CMP_REPORT) -> None:
    """Produce a side-by-side chem_target vs STP comparison report."""
    if not chem_results_csv.exists():
        print(f"WARNING: {chem_results_csv} not found — skipping comparison.")
        return

    # Load chem_target results
    chem_stats: dict[str, dict] = {}
    with open(chem_results_csv, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            # Column name may be true_target_class or true_class
            cls = (row.get("true_target_class") or
                   row.get("true_class") or "unknown")
            if cls not in chem_stats:
                chem_stats[cls] = {"n": 0, "top1": 0, "top3": 0, "mrr_sum": 0.0}
            s = chem_stats[cls]
            s["n"] += 1
            # top1_hit / top3_hit may be "True"/"False" or "1"/"0"
            t1 = row.get("top1_hit", "0")
            t3 = row.get("top3_hit", "0")
            s["top1"] += 1 if t1 in ("True", "1") else 0
            s["top3"] += 1 if t3 in ("True", "1") else 0
            try:
                s["mrr_sum"] += float(row.get("mrr", 0))
            except ValueError:
                pass

    # Compute totals
    def _totals(stats: dict) -> tuple[int, int, int, float]:
        n = sum(s["n"] for s in stats.values())
        t1 = sum(s["top1"] for s in stats.values())
        t3 = sum(s["top3"] for s in stats.values())
        mrr = sum(s["mrr_sum"] for s in stats.values()) / max(n, 1)
        return n, t1, t3, mrr

    sn, st1, st3, smrr = _totals(stp_stats)
    cn, ct1, ct3, cmrr = _totals(chem_stats)

    all_classes = sorted(set(stp_stats) | set(chem_stats))

    lines: list[str] = [
        "=" * 80,
        "  chem_target  vs  SwissTargetPrediction — Side-by-Side Comparison",
        "=" * 80,
        "",
        "Both benchmarks use the same ChEMBL curated compound set (≥1 µM activity).",
        "Evaluation: class-level Top-1/Top-3/MRR using ChEMBL ID + class label mapping.",
        "",
        "OVERALL SUMMARY",
        "-" * 40,
        f"  {'Metric':<26} {'chem_target':>12} {'STP':>12}",
        "  " + "-" * 52,
        f"  {'Total compounds':<26} {cn:>12} {sn:>12}",
        f"  {'Top-1 accuracy':<26} {ct1/cn*100:>11.1f}% {st1/sn*100:>11.1f}%",
        f"  {'Top-3 accuracy':<26} {ct3/cn*100:>11.1f}% {st3/sn*100:>11.1f}%",
        f"  {'Mean Reciprocal Rank':<26} {cmrr:>12.3f} {smrr:>12.3f}",
        "",
        "PER-CLASS BREAKDOWN",
        "-" * 40,
        f"  {'Target class':<22}  {'N':>4}  "
        f"{'cT-1':>5}  {'cT-3':>5}  {'cMRR':>5}  |  "
        f"{'STP-1':>5}  {'STP-3':>5}  {'SMRR':>5}",
        "  " + "-" * 76,
    ]

    for cls in all_classes:
        c = chem_stats.get(cls, {"n": 0, "top1": 0, "top3": 0, "mrr_sum": 0.0})
        s = stp_stats.get(cls, {"n": 0, "top1": 0, "top3": 0, "mrr_sum": 0.0})
        n = c["n"] or s["n"]
        if n == 0:
            continue
        cn2 = c["n"] or 1
        sn2 = s["n"] or 1
        lines.append(
            f"  {cls:<22}  {n:>4}  "
            f"{c['top1']/cn2*100:>4.0f}%  {c['top3']/cn2*100:>4.0f}%  "
            f"{c['mrr_sum']/cn2:>5.3f}  |  "
            f"{s['top1']/sn2*100:>4.0f}%  {s['top3']/sn2*100:>4.0f}%  "
            f"{s['mrr_sum']/sn2:>5.3f}"
        )

    lines += [
        "",
        "INTERPRETATION NOTES",
        "-" * 40,
        "  chem_target bias: BioLiP/PDB structural circularity (overestimates compounds",
        "    with deposited crystal structures).",
        "  STP bias:         ChEMBL fingerprint similarity (compounds from ChEMBL",
        "    may be near-duplicates of STP training data, overestimating performance).",
        "  Both tools are best used together: chem_target provides mechanism-level",
        "    insights (which FGs drive binding); STP provides rapid ligand-based screening.",
        "",
        "FILES",
        "-" * 40,
        f"  STP results:         {_STP_RESULTS}",
        f"  chem_target results: {_CHEM_RESULTS}",
        "",
        "=" * 80,
        "",
    ]

    report_text = "\n".join(lines)
    out_report.write_text(report_text, encoding="utf-8")
    print(report_text)
    print(f"Comparison report saved: {out_report}")


# ── Load helpers ─────────────────────────────────────────────────────────────

def _load_curated_compounds(limit: int | None = None) -> list[dict]:
    """Load compounds.csv from the curated benchmark download."""
    if not _CURATED_CSV.exists():
        raise FileNotFoundError(
            f"Curated compounds not found: {_CURATED_CSV}\n"
            f"Run first: python run_benchmark.py download --mode curated"
        )
    compounds: list[dict] = []
    with open(_CURATED_CSV, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("smiles"):
                compounds.append(row)
    if limit:
        compounds = compounds[:limit]
    return compounds


def _load_stp_results() -> list[dict]:
    """Load previously saved per-compound STP results."""
    if not _STP_RESULTS.exists():
        raise FileNotFoundError(
            f"STP results not found: {_STP_RESULTS}\n"
            f"Run without --report-only first."
        )
    results: list[dict] = []
    with open(_STP_RESULTS, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            # Convert bool strings
            row["stp_top1"] = row.get("stp_top1", "False") == "True"
            row["stp_top3"] = row.get("stp_top3", "False") == "True"
            try:
                row["stp_mrr"] = float(row.get("stp_mrr", 0))
            except ValueError:
                row["stp_mrr"] = 0.0
            try:
                row["stp_n_preds"] = int(row.get("stp_n_preds", 0))
            except ValueError:
                row["stp_n_preds"] = 0
            results.append(row)
    return results


def _rebuild_stats_from_results(results: list[dict]) -> dict[str, dict]:
    """Recompute per-class stats from saved per-compound results."""
    stats: dict[str, dict] = {}
    for r in results:
        cls = r.get("true_class", "unknown")
        if cls not in stats:
            stats[cls] = {"n": 0, "top1": 0, "top3": 0, "mrr_sum": 0.0}
        s = stats[cls]
        s["n"]       += 1
        s["top1"]    += int(r.get("stp_top1", False))
        s["top3"]    += int(r.get("stp_top3", False))
        s["mrr_sum"] += float(r.get("stp_mrr", 0.0))
    return stats


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    """Parse CLI args and run the requested phases."""
    parser = argparse.ArgumentParser(
        description="SwissTargetPrediction comparison benchmark"
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Limit to first N compounds (for quick testing)",
    )
    parser.add_argument(
        "--report-only", action="store_true",
        help="Skip STP queries; regenerate report from cached stp_results.csv",
    )
    parser.add_argument(
        "--compare", action="store_true",
        help="Generate side-by-side comparison (requires both benchmarks already run)",
    )
    parser.add_argument(
        "--delay", type=float, default=_REQUEST_DELAY,
        help=f"Seconds between STP requests (default {_REQUEST_DELAY})",
    )
    args = parser.parse_args()

    if args.report_only:
        print("=== Loading cached STP results ===")
        results = _load_stp_results()
        stp_stats = _rebuild_stats_from_results(results)
        generate_report(results)
    elif args.compare:
        print("=== Generating comparison report ===")
        results = _load_stp_results()
        stp_stats = _rebuild_stats_from_results(results)
        generate_comparison(stp_stats)
    else:
        # Full run
        print("=== PHASE: LOAD ===")
        compounds = _load_curated_compounds(limit=args.limit)
        print(f"Loaded {len(compounds)} compounds")

        print("\n=== PHASE: STP QUERIES ===")
        print(f"Delay between requests: {args.delay}s")
        results = run_stp_queries(compounds, delay=args.delay)

        print("\n=== PHASE: REPORT ===")
        stp_stats = generate_report(results)

        if _CHEM_RESULTS.exists():
            print("\n=== PHASE: COMPARISON ===")
            generate_comparison(stp_stats)

    print("\nAll done.")


if __name__ == "__main__":
    main()
