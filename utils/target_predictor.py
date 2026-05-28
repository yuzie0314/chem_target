"""Target class and residue interaction prediction from a SMILES string.

Workflow
--------
1. Detect functional groups (SMARTS-based) in the input SMILES.
2. Load the pre-built FG × residue co-occurrence table (db/fg_residue_table.csv).
   rows  = 20 standard amino acids (1-letter code)
   cols  = FG names from FG_SMARTS
   value = binding-event count from BioLiP 2.0
3. Score each residue by **z-score-normalising** each FG column then summing
   (prevents high-frequency generic FGs from drowning specific pharmacophores).
4. Vote target classes using known_target_classes from db/fg_database.json,
   weighted by **IDF** (specific targets ranked above generic labels).
5. Return a structured result dict + a formatted text report.

Scoring details
---------------
Residue scoring (z-score normalisation, default on):
    For each detected FG, normalise its column to mean=0 / std=1 across all 20 AAs
    before summing.  This ensures that a highly represented FG like Hydroxyl
    (>11 000 GLY events) does not completely overshadow Steroid (~1 400 PHE events).

Target class scoring (IDF weighting, default on):
    weight(tc) = log( N_all_FGs / N_FGs_that_annotate_this_tc )
    High weight → target class is specific (few FGs annotate it, e.g. VKORC1, tubulin).
    Low weight  → target class is generic (many FGs annotate it, e.g. kinase, GPCR).
    final_score = vote_count × weight

Prerequisites
-------------
Build the FG × residue table once (uses local BioLiP_nr.txt.gz):

    python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz

Outputs: db/fg_residue_table.csv  (tracked in git after first build)

Quick test with the first 1000 entries:

    python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz --top 1000

Standalone usage
----------------
    python utils/target_predictor.py "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin
    python utils/target_predictor.py "CC(=O)Oc1ccccc1C(=O)O" --top 5
"""

import json
import sys
from collections import Counter, defaultdict
from math import log
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from utils.fg_detector import detect_smarts  # noqa: E402

# ── Paths ──────────────────────────────────────────────────────────────────────
TABLE_PATH = _ROOT / "db" / "fg_residue_table.csv"
FG_DB_PATH = _ROOT / "db" / "fg_database.json"

# ── Amino acid code lookup ─────────────────────────────────────────────────────
AA_1TO3: dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


# ── Data loaders ───────────────────────────────────────────────────────────────

def load_table(path: Path = TABLE_PATH) -> pd.DataFrame:
    """Load the pre-built FG × residue co-occurrence table.

    The table has:
        index  = residue (1-letter AA code)
        cols   = "residue_name" + FG names from FG_SMARTS
        values = integer co-occurrence counts from BioLiP 2.0

    Raises:
        FileNotFoundError: with instructions to build the table if missing.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"\nFG-residue table not found: {path}\n\n"
            "Build it first (requires db/BioLiP_nr.txt.gz):\n"
            "  python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz\n\n"
            "Quick test (first 1 000 entries):\n"
            "  python utils/interaction_analyzer.py --local db/BioLiP_nr.txt.gz --top 1000"
        )
    df = pd.read_csv(path, index_col="residue")
    return df


def load_fg_db(path: Path = FG_DB_PATH) -> dict:
    """Load fg_database.json and return the functional_groups sub-dict.

    Returns an empty dict if the file is missing (graceful fallback).
    """
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return data.get("functional_groups", {})


# ── Scoring helpers ────────────────────────────────────────────────────────────

def _compute_target_idf(fg_db: dict) -> dict[str, float]:
    """Compute inverse-document-frequency weights for target classes.

    IDF(tc) = log( N_all_FGs / N_FGs_that_annotate_tc )

    Interpretation:
        High IDF → target is specific (few FGs list it, e.g. VKORC1, tubulin).
        Low IDF  → target is generic (many FGs list it, e.g. kinase, GPCR).

    Args:
        fg_db: functional_groups dict loaded from fg_database.json.

    Returns:
        Dict mapping target class name → IDF weight (float).
    """
    n_fgs = len(fg_db)
    if n_fgs == 0:
        return {}
    tc_count: Counter = Counter()
    for entry in fg_db.values():
        for tc in entry.get("known_target_classes", []):
            tc_count[tc] += 1
    return {tc: log(n_fgs / count) for tc, count in tc_count.items()}


def predict_residues(
    fgs_detected: list[str],
    table: pd.DataFrame,
    top_n: int = 20,
    normalize: bool = True,
) -> pd.DataFrame:
    """Score residues by summing BioLiP co-occurrence counts for detected FGs.

    When *normalize* is True (default), each FG column is z-score-normalised
    (column mean → 0, column std → 1) before summing.  This prevents dominant
    FGs (e.g. Hydroxyl: ~11 000 GLY hits) from overshadowing structurally
    specific FGs (e.g. Steroid: ~1 400 PHE hits) when a compound contains both.

    Args:
        fgs_detected: FG names matching FG_SMARTS keys.
        table:        Numeric FG × residue table (residue_name column excluded).
        top_n:        Return only the top-N residues by score (0 = return all).
        normalize:    Apply per-FG z-score normalisation before summing.

    Returns:
        DataFrame with columns: residue_name | score | contributing_fgs
        Index: residue (1-letter AA code).  Rows with score ≤ 0 are dropped.
    """
    valid_fgs = [fg for fg in fgs_detected if fg in table.columns]

    if not valid_fgs:
        return pd.DataFrame(
            columns=["residue_name", "score", "contributing_fgs"]
        )

    sub = table[valid_fgs].astype(float)

    if normalize:
        # Per-column z-score.  Avoid division by zero for all-zero FG columns.
        col_std = sub.std()
        col_std[col_std == 0] = 1.0
        sub = (sub - sub.mean()) / col_std

    score_series = sub.sum(axis=1)

    result = pd.DataFrame({
        "residue_name": score_series.index.map(lambda aa: AA_1TO3.get(aa, aa)),
        "score": score_series,
        "contributing_fgs": ", ".join(valid_fgs),
    })
    result.index.name = "residue"
    result = result[result["score"] > 0].sort_values("score", ascending=False)

    if top_n and top_n > 0:
        result = result.head(top_n)

    return result


def predict_target_classes(
    fgs_detected: list[str],
    fg_db: dict,
    use_idf: bool = True,
) -> pd.DataFrame:
    """Vote target classes using known_target_classes from fg_database.json.

    Each detected FG casts 1 vote to every target class it annotates.
    When *use_idf* is True (default), votes are multiplied by IDF weight:

        score(tc) = vote_count × log( N_all_FGs / N_FGs_annotating_tc )

    This gives higher final scores to specific targets (e.g. VKORC1, tubulin,
    antimalarial target) than to generic labels (e.g. kinase, GPCR) even when
    the generic label accumulates more raw votes.

    Args:
        fgs_detected: FG names matching FG_SMARTS keys.
        fg_db:        Loaded functional_groups dict from fg_database.json.
        use_idf:      Apply IDF weighting (default True).

    Returns:
        DataFrame with columns: target_class | score | votes | evidence_fgs
        Sorted by score descending.  Empty DataFrame if no annotations found.
    """
    votes: Counter = Counter()
    evidence: dict[str, list[str]] = defaultdict(list)

    idf: dict[str, float] = _compute_target_idf(fg_db) if use_idf else {}

    for fg in fgs_detected:
        entry = fg_db.get(fg, {})
        for tc in entry.get("known_target_classes", []):
            votes[tc] += 1
            evidence[tc].append(fg)

    if not votes:
        return pd.DataFrame(
            columns=["target_class", "score", "votes", "evidence_fgs"]
        )

    rows = []
    for tc, count in votes.most_common():
        weight = idf.get(tc, 1.0) if use_idf else 1.0
        rows.append({
            "target_class": tc,
            "score": round(count * weight, 3),
            "votes": count,
            "evidence_fgs": ", ".join(evidence[tc]),
        })

    df = (
        pd.DataFrame(rows)
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    return df


# ── Main prediction pipeline ───────────────────────────────────────────────────

def predict(
    smiles: str,
    top_residues: int = 10,
    normalize: bool = True,
    use_idf: bool = True,
    table_path: Path = TABLE_PATH,
    db_path: Path = FG_DB_PATH,
) -> dict:
    """Full target prediction pipeline for one SMILES string.

    Args:
        smiles:       Input molecule as SMILES.
        top_residues: How many top-scoring residues to return.
        normalize:    Z-score-normalise residue scores (default True).
        use_idf:      IDF-weight target class votes (default True).
        table_path:   Path to db/fg_residue_table.csv.
        db_path:      Path to db/fg_database.json.

    Returns:
        dict with keys:
          smiles             (str)
          fgs_detected       (list[str])
          residue_scores     (pd.DataFrame)  — residue | residue_name | score | contributing_fgs
          target_class_votes (pd.DataFrame)  — target_class | score | votes | evidence_fgs
          warning            (str | None)    — set on SMILES parse failure or missing table
    """
    result: dict = {
        "smiles": smiles,
        "fgs_detected": [],
        "residue_scores": pd.DataFrame(),
        "target_class_votes": pd.DataFrame(),
        "warning": None,
    }

    # 1. Detect FGs (SMARTS-based, consistent with the BioLiP table)
    fgs = detect_smarts(smiles)
    if not fgs:
        result["warning"] = (
            f"No SMARTS-matched functional groups detected in: {smiles}\n"
            "  Check that the SMILES is valid and that FG_SMARTS covers the scaffold."
        )
        return result
    result["fgs_detected"] = fgs

    # 2. Load FG × residue table
    try:
        table = load_table(table_path)
    except FileNotFoundError as exc:
        result["warning"] = str(exc)
        return result

    # 3. Score residues (drop residue_name label column before scoring)
    numeric_table = table.drop(columns=["residue_name"], errors="ignore")
    result["residue_scores"] = predict_residues(
        fgs, numeric_table, top_n=top_residues, normalize=normalize
    )

    # 4. Vote + weight target classes
    fg_db = load_fg_db(db_path)
    result["target_class_votes"] = predict_target_classes(
        fgs, fg_db, use_idf=use_idf
    )

    return result


# ── Report formatter ───────────────────────────────────────────────────────────

def format_report(pred: dict, compound_name: str = "") -> str:
    """Format a prediction result dict as a human-readable text report.

    Args:
        pred:          Result dict returned by predict().
        compound_name: Optional display name for the compound.

    Returns:
        Formatted multi-line string.
    """
    label = compound_name or pred["smiles"]
    title = f"Target Prediction — {label}"
    bar   = "=" * len(title)
    lines = [bar, title, bar]

    if pred["warning"]:
        lines.append(f"\n⚠  {pred['warning']}")
        return "\n".join(lines)

    lines.append(f"\n  SMILES : {pred['smiles']}")

    # FGs detected
    fgs = pred["fgs_detected"]
    lines.append(f"\n  Functional groups detected ({len(fgs)}):")
    for fg in fgs:
        lines.append(f"    • {fg}")

    # Residue scores
    lines.append("\n  Top binding residues (z-score-normalised BioLiP FG×residue):")
    rs = pred["residue_scores"]
    if rs.empty:
        lines.append("    (no residue matches in table)")
    else:
        lines.append(f"    {'Res':>4}  {'Name':>5}  {'Score':>8}")
        lines.append(f"    {'---':>4}  {'----':>5}  {'-----':>8}")
        for aa, row in rs.iterrows():
            lines.append(
                f"    {aa:>4}  {row['residue_name']:>5}  {row['score']:>8.3f}"
            )

    # Target class votes (IDF-weighted)
    lines.append("\n  Predicted target classes (IDF-weighted FG votes):")
    tc = pred["target_class_votes"]
    if tc.empty:
        lines.append("    (no annotations in fg_database.json for detected FGs)")
    else:
        lines.append(
            f"    {'Target class':<30}  {'Score':>6}  {'Votes':>5}  Evidence FGs"
        )
        lines.append(
            f"    {'-'*30}  {'------':>6}  {'-----':>5}  {'-'*30}"
        )
        for _, row in tc.iterrows():
            lines.append(
                f"    {row['target_class']:<30}  {row['score']:>6.2f}"
                f"  {int(row['votes']):>5}  {row['evidence_fgs']}"
            )

    lines.append("")
    return "\n".join(lines)


# ── Standalone entry point ─────────────────────────────────────────────────────

def main() -> None:
    """CLI wrapper: predict targets for a single SMILES string."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Predict target residues and protein classes for a SMILES string"
    )
    parser.add_argument("smiles", help="Input SMILES string")
    parser.add_argument("--name",  default="",  help="Compound display name")
    parser.add_argument("--top",   type=int, default=10,
                        help="Number of top residues to show (default: 10)")
    parser.add_argument("--no-normalize", action="store_true",
                        help="Disable z-score normalisation of residue scores")
    parser.add_argument("--no-idf", action="store_true",
                        help="Disable IDF weighting of target class votes")
    parser.add_argument("--table", default=str(TABLE_PATH),
                        help="Path to fg_residue_table.csv")
    parser.add_argument("--db",    default=str(FG_DB_PATH),
                        help="Path to fg_database.json")
    args = parser.parse_args()

    pred = predict(
        args.smiles,
        top_residues=args.top,
        normalize=not args.no_normalize,
        use_idf=not args.no_idf,
        table_path=Path(args.table),
        db_path=Path(args.db),
    )
    print(format_report(pred, compound_name=args.name))


if __name__ == "__main__":
    main()
