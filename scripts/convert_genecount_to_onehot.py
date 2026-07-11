#!/usr/bin/env python3
"""Convert OrthoFinder gene count TSV to one-hot (presence/absence) format.

Reads Orthogroups.GeneCount.tsv where each cell contains an integer gene count,
and outputs a TSV where each cell is 1 (gene present, count > 0) or 0 (absent).
The 'Total' column is dropped from the output.

Usage:
    python convert_genecount_to_onehot.py [INPUT_TSV] [OUTPUT_TSV]

    If no arguments are provided, uses default project paths.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def convert_genecount_to_onehot(input_path: Path, output_path: Path) -> None:
    """Convert gene count matrix to one-hot presence/absence matrix.

    Args:
        input_path: Path to Orthogroups.GeneCount.tsv
        output_path: Path to write one-hot output TSV
    """
    df = pd.read_csv(input_path, sep="\t", index_col=0)

    if "Total" in df.columns:
        df = df.drop(columns=["Total"])

    one_hot = (df > 0).astype(int)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    one_hot.to_csv(output_path, sep="\t")
    print(f"Converted {len(one_hot)} orthogroups × {len(one_hot.columns)} species")
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert OrthoFinder gene count TSV to one-hot format."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="Input Orthogroups.GeneCount.tsv path",
    )
    parser.add_argument(
        "output",
        nargs="?",
        default=None,
        help="Output one-hot TSV path",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent

    input_path = (
        Path(args.input)
        if args.input
        else project_root
        / "Outputs"
        / "OrthoFinder"
        / "Results_Mar03"
        / "Orthogroups"
        / "Orthogroups.GeneCount.tsv"
    )
    output_path = (
        Path(args.output)
        if args.output
        else project_root
        / "Notebooks"
        / "Tables_set2"
        / "one-hot-orthogroups-gene-count.tsv"
    )

    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    convert_genecount_to_onehot(input_path, output_path)


if __name__ == "__main__":
    main()
