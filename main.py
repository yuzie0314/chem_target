"""chem_target — entry point.

Subcommands
-----------
fg       Functional group analysis (abundance table + SVG images)
predict  Target prediction (residue scores + target class votes)

Usage examples
--------------
  # FG analysis
  python main.py fg --input data/compounds.csv
  python main.py fg --input data/compounds.csv --output fg_table.csv --images output/images

  # Target prediction — single SMILES
  python main.py predict --smiles "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin

  # Target prediction — batch CSV, save summary + individual reports
  python main.py predict --input data/compounds.csv
  python main.py predict --input data/compounds.csv --output output/predictions.csv
  python main.py predict --input data/compounds.csv --output output/predictions.csv \\
      --report-dir output/reports --top 10
"""

import argparse
import csv
import sys
from pathlib import Path


# ── fg subcommand ──────────────────────────────────────────────────────────────

def run_fg(args: argparse.Namespace) -> None:
    """Functional group analysis pipeline."""
    from utils.io_handler import read_file
    from utils.fg_detector import detect_smarts_table
    from utils.visualizer import draw_compounds

    input_file = args.input
    if not Path(input_file).exists():
        print(f"Error: input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading compounds from {input_file} ...")
    compounds = read_file(input_file, name_property=args.name_property)
    print(f"  {len(compounds)} compound(s) loaded.")

    print("Detecting functional groups ...")
    df = detect_smarts_table(compounds)

    df.to_csv(args.output, index=True)
    print(f"Abundance table saved: {args.output}")
    print(df.to_string())

    print("\nRendering SVG images ...")
    draw_compounds(compounds, output_dir=args.images)

    print("\nDone.")


# ── predict subcommand ─────────────────────────────────────────────────────────

_CSV_TARGETS   = 3   # top-N target classes in summary CSV
_CSV_RESIDUES  = 3   # top-N residues in summary CSV


def _prediction_to_row(name: str, smiles: str, pred: dict) -> dict:
    """Convert a single prediction result dict to a flat CSV row.

    Columns:
        compound        — display name
        smiles          — input SMILES
        n_fgs           — number of FGs detected
        fgs             — FG names, comma-separated
        target_1..3     — top-3 predicted target classes
        score_1..3      — IDF-weighted scores for those targets
        votes_1..3      — raw vote counts
        residue_1..3    — top-3 binding residues (3-letter AA code)
        res_score_1..3  — z-score-normalised residue scores
        warning         — any error / missing-FG message

    Args:
        name:   Compound display name.
        smiles: Input SMILES string.
        pred:   Result dict from predict().

    Returns:
        Flat dict suitable for csv.DictWriter.
    """
    row: dict = {
        "compound":  name,
        "smiles":    smiles,
        "n_fgs":     len(pred["fgs_detected"]),
        "fgs":       ", ".join(pred["fgs_detected"]),
        "warning":   pred.get("warning") or "",
    }

    # Target class columns
    tc = pred["target_class_votes"]
    for i in range(1, _CSV_TARGETS + 1):
        if not tc.empty and len(tc) >= i:
            r = tc.iloc[i - 1]
            row[f"target_{i}"]  = r["target_class"]
            row[f"score_{i}"]   = round(r["score"], 3)
            row[f"votes_{i}"]   = int(r["votes"])
        else:
            row[f"target_{i}"] = ""
            row[f"score_{i}"]  = ""
            row[f"votes_{i}"]  = ""

    # Residue columns
    rs = pred["residue_scores"]
    for i in range(1, _CSV_RESIDUES + 1):
        if not rs.empty and len(rs) >= i:
            r = rs.iloc[i - 1]
            row[f"residue_{i}"]    = r["residue_name"]
            row[f"res_score_{i}"]  = round(r["score"], 3)
        else:
            row[f"residue_{i}"]   = ""
            row[f"res_score_{i}"] = ""

    return row


def run_predict(args: argparse.Namespace) -> None:
    """Target prediction pipeline.

    For each compound, runs the full prediction and:
      - Prints a formatted text report to stdout.
      - (If --output)     appends a one-row CSV summary entry.
      - (If --report-dir) saves a <name>.txt individual report.
    """
    from utils.io_handler import read_file
    from utils.target_predictor import predict, format_report

    # ── Build compound dict ────────────────────────────────────────────────────
    if args.smiles:
        compounds = {args.name or args.smiles: args.smiles}
    else:
        if not args.input:
            print(
                "Error: supply either --input <file> or --smiles <SMILES>",
                file=sys.stderr,
            )
            sys.exit(1)
        if not Path(args.input).exists():
            print(f"Error: input file not found: {args.input}", file=sys.stderr)
            sys.exit(1)
        print(f"Reading compounds from {args.input} ...")
        compounds = read_file(args.input, name_property=args.name_property)
        print(f"  {len(compounds)} compound(s) loaded.\n")

    # ── Optional output paths ──────────────────────────────────────────────────
    csv_path    = Path(args.output)     if args.output     else None
    report_dir  = Path(args.report_dir) if args.report_dir else None

    if report_dir:
        report_dir.mkdir(parents=True, exist_ok=True)

    # ── CSV header (written once) ──────────────────────────────────────────────
    csv_fieldnames = (
        ["compound", "smiles", "n_fgs", "fgs"]
        + [f"target_{i}"  for i in range(1, _CSV_TARGETS  + 1)]
        + [f"score_{i}"   for i in range(1, _CSV_TARGETS  + 1)]
        + [f"votes_{i}"   for i in range(1, _CSV_TARGETS  + 1)]
        + [f"residue_{i}" for i in range(1, _CSV_RESIDUES + 1)]
        + [f"res_score_{i}" for i in range(1, _CSV_RESIDUES + 1)]
        + ["warning"]
    )

    csv_fh     = None
    csv_writer = None
    if csv_path:
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        csv_fh = open(csv_path, "w", newline="", encoding="utf-8")
        csv_writer = csv.DictWriter(csv_fh, fieldnames=csv_fieldnames)
        csv_writer.writeheader()

    # ── Predict each compound ──────────────────────────────────────────────────
    try:
        for name, smiles in compounds.items():
            pred    = predict(smiles, top_residues=args.top)
            report  = format_report(pred, compound_name=name)

            # stdout
            print(report)

            # individual text report
            if report_dir:
                safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
                (report_dir / f"{safe_name}.txt").write_text(report, encoding="utf-8")

            # CSV summary row
            if csv_writer:
                csv_writer.writerow(_prediction_to_row(name, smiles, pred))
    finally:
        if csv_fh:
            csv_fh.close()

    if csv_path:
        print(f"\nSummary CSV saved: {csv_path}")
    if report_dir:
        print(f"Individual reports saved: {report_dir}/")


# ── Argument parser ────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="chem_target",
        description="Compound-to-target prediction tool",
    )
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    subparsers.required = True

    # ── fg ────────────────────────────────────────────────────────────────────
    fg_p = subparsers.add_parser(
        "fg",
        help="Functional group analysis — abundance table + SVG images",
    )
    fg_p.add_argument(
        "--input", default="data/compounds.csv",
        help="Input CSV or SDF file (default: data/compounds.csv)",
    )
    fg_p.add_argument(
        "--output", default="fg_abundance_table.csv",
        help="Output CSV path (default: fg_abundance_table.csv)",
    )
    fg_p.add_argument(
        "--images", default="output/images",
        help="Output directory for SVG images (default: output/images)",
    )
    fg_p.add_argument(
        "--name-property", default=None, dest="name_property",
        help="SDF property tag to use as compound name (e.g. 'ChEMBL_ID'). "
             "Ignored for CSV. Falls back to molecule title line if not set.",
    )

    # ── predict ───────────────────────────────────────────────────────────────
    pred_p = subparsers.add_parser(
        "predict",
        help="Target prediction — residue scores + target class votes",
    )
    pred_p.add_argument(
        "--input", default=None,
        help="Input CSV or SDF file",
    )
    pred_p.add_argument(
        "--name-property", default=None, dest="name_property",
        help="SDF property tag to use as compound name (e.g. 'ChEMBL_ID'). "
             "Ignored for CSV. Falls back to molecule title line if not set.",
    )
    pred_p.add_argument(
        "--smiles", default=None,
        help="Single SMILES string (alternative to --input)",
    )
    pred_p.add_argument(
        "--name", default="",
        help="Compound display name (used with --smiles)",
    )
    pred_p.add_argument(
        "--top", type=int, default=10,
        help="Top N residues to show per compound (default: 10)",
    )
    pred_p.add_argument(
        "--output", default=None,
        help="Save a one-row-per-compound summary CSV to this path "
             "(e.g. output/predictions.csv). Directory is created if needed.",
    )
    pred_p.add_argument(
        "--report-dir", default=None, dest="report_dir",
        help="Save a <compound_name>.txt full report for each compound "
             "to this directory (e.g. output/reports).",
    )

    return parser


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "fg":
        run_fg(args)
    elif args.command == "predict":
        run_predict(args)
