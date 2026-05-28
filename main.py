"""chem_target — entry point.

Subcommands
-----------
fg       Functional group analysis (abundance table + SVG images)
predict  Target prediction (residue scores + target class votes)

Usage examples
--------------
  python main.py fg --input data/compounds.csv
  python main.py fg --input data/compounds.csv --output my_table.csv --images output/images

  python main.py predict --input data/compounds.csv
  python main.py predict --input data/compounds.csv --top 5

  # Single SMILES (predict only):
  python main.py predict --smiles "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin
"""

import argparse
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

def run_predict(args: argparse.Namespace) -> None:
    """Target prediction pipeline."""
    from utils.io_handler import read_file
    from utils.target_predictor import predict, format_report

    # Build compound dict from file or single --smiles flag
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

    for name, smiles in compounds.items():
        pred = predict(smiles, top_residues=args.top)
        print(format_report(pred, compound_name=name))


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
        help="Top N residues to show (default: 10)",
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
