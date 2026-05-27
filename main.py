"""chem_target — entry point.

Usage:
    python main.py
    python main.py --input data/compounds.csv --output fg_abundance_table.csv
"""

import argparse
import sys
from pathlib import Path

from utils.io_handler import read_file
from utils.fg_detector import detect
from utils.visualizer import draw_compounds


def main(input_file: str, output_csv: str, image_dir: str) -> None:
    """Run the full functional group analysis pipeline."""
    if not Path(input_file).exists():
        print(f"Error: input file not found: {input_file}")
        sys.exit(1)

    print(f"Reading compounds from {input_file} ...")
    compounds = read_file(input_file)
    print(f"  {len(compounds)} compound(s) loaded.")

    print("Detecting functional groups ...")
    df = detect(compounds)

    df.to_csv(output_csv, index=True)
    print(f"Abundance table saved: {output_csv}")
    print(df.to_string())

    print("\nRendering SVG images ...")
    draw_compounds(compounds, output_dir=image_dir)

    print("\nDone.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Functional group analysis pipeline")
    parser.add_argument(
        "--input", default="data/compounds.csv", help="Input CSV file"
    )
    parser.add_argument(
        "--output", default="fg_abundance_table.csv", help="Output CSV file"
    )
    parser.add_argument(
        "--images", default="output/images", help="Output directory for SVG images"
    )
    args = parser.parse_args()

    main(args.input, args.output, args.images)
