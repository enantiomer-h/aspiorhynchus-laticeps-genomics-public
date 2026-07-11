#!/usr/bin/env python3
"""
Find common enrichment terms across different evolutionary analyses.

This script identifies and summarizes GO and KEGG KO terms that are
significantly enriched across multiple CATEGORIES of analysis types
(CAFE, KaKs, unique). Terms must span at least two different categories
to be considered "common" - this filters out same-category overlaps
(e.g., multiple KaKs subcategories) in favor of biologically meaningful
cross-category overlaps.

Analysis Categories:
  - CAFE: CAFE_contracted, CAFE_expanded
  - KaKs: KaKs_positive_tentative, KaKs_relaxed_purifying, KaKs_strong_purifying
  - unique: unique_nonzero

Usage:
    python find_common_enrichment_terms.py \\
        --input-dir /path/to/enrichment/results \\
        --output-dir /path/to/output \\
        --p-threshold 0.05 \\
        --min-categories 2
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


# Analysis type names (extracted from filenames)
ANALYSIS_TYPES = [
    "CAFE_contracted",
    "CAFE_expanded",
    "KaKs_positive_tentative",
    "KaKs_relaxed_purifying",
    "KaKs_strong_purifying",
    "unique_nonzero",
]

# Map analysis types to their parent categories
# Cross-category overlaps are more biologically meaningful than within-category
ANALYSIS_CATEGORIES = {
    "CAFE_contracted": "CAFE",
    "CAFE_expanded": "CAFE",
    "KaKs_positive_tentative": "KaKs",
    "KaKs_relaxed_purifying": "KaKs",
    "KaKs_strong_purifying": "KaKs",
    "unique_nonzero": "unique",
}

# GO categories
GO_CATEGORIES = ["BP", "CC", "MF"]

# Required columns for enrichment results
REQUIRED_COLUMNS = ["ID", "Description", "p.adjust", "Count", "geneID"]


def get_unique_categories(analyses: List[str]) -> set:
    """Return the set of unique categories for a list of analysis names.

    This is used to filter for cross-category overlaps, which are more
    biologically meaningful than within-category overlaps.
    """
    return {ANALYSIS_CATEGORIES[a] for a in analyses}


def load_ko_enrichment(input_dir: Path, analysis_type: str) -> pd.DataFrame:
    """Load a KO enrichment TSV file."""
    filename = f"KO_enrichment_{analysis_type}_results.tsv"
    filepath = input_dir / filename

    if not filepath.exists():
        print(f"  Warning: {filename} not found, skipping")
        return pd.DataFrame()

    df = pd.read_csv(filepath, sep="\t")

    # Verify required columns exist
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        print(f"  Warning: {filename} missing columns {missing_cols}, skipping")
        return pd.DataFrame()

    return df


def load_go_enrichment(input_dir: Path, analysis_type: str, category: str) -> pd.DataFrame:
    """Load a GO enrichment CSV file."""
    filename = f"GO_enrichment_{analysis_type}_{category}_results.csv"
    filepath = input_dir / filename

    if not filepath.exists():
        print(f"  Warning: {filename} not found, skipping")
        return pd.DataFrame()

    df = pd.read_csv(filepath)

    # Verify required columns exist
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        print(f"  Warning: {filename} missing columns {missing_cols}, skipping")
        return pd.DataFrame()

    return df


def find_common_terms(
    data: Dict[str, pd.DataFrame],
    min_categories: int = 2,
    p_threshold: float = 0.05
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Find enrichment terms common across multiple analysis categories.

    Terms must appear in at least `min_categories` different analysis categories
    (CAFE, KaKs, unique) to be considered "common". This filters out same-category
    overlaps (e.g., KaKs_positive + KaKs_relaxed) which are less biologically
    meaningful than cross-category overlaps (e.g., CAFE_expanded + KaKs_positive).

    Args:
        data: Dictionary mapping analysis_name -> enrichment DataFrame
        min_categories: Minimum number of different categories a term must span
        p_threshold: Maximum adjusted p-value to consider significant

    Returns:
        Tuple of (summary_df, detailed_df, presence_matrix_df)
    """
    # Filter each dataframe by p-value threshold and collect term info
    term_occurrences: Dict[str, List[dict]] = {}

    for analysis_name, df in data.items():
        if df.empty:
            continue

        # Filter by p-value threshold
        sig_df = df[df["p.adjust"] <= p_threshold].copy()

        for _, row in sig_df.iterrows():
            term_id = row["ID"]
            if term_id not in term_occurrences:
                term_occurrences[term_id] = []

            term_occurrences[term_id].append({
                "analysis": analysis_name,
                "description": row["Description"],
                "p_adjust": row["p.adjust"],
                "count": row["Count"],
                "gene_ids": row["geneID"],
            })

    # Filter to terms spanning >= min_categories different analysis categories
    # This ensures cross-category overlaps (biologically meaningful) rather than
    # same-category overlaps (e.g., multiple KaKs sub-categories)
    common_terms = {
        term_id: occurrences
        for term_id, occurrences in term_occurrences.items()
        if len(get_unique_categories([occ["analysis"] for occ in occurrences])) >= min_categories
    }

    if not common_terms:
        # Return empty dataframes with correct structure
        summary_cols = ["ID", "Description", "num_categories", "categories", "num_analyses", "analyses", "min_p_adjust", "max_p_adjust", "total_count"]
        detailed_cols = ["ID", "Description", "analysis", "p_adjust", "count", "gene_ids"]
        presence_cols = ["ID", "Description"] + list(data.keys())
        return (
            pd.DataFrame(columns=summary_cols),
            pd.DataFrame(columns=detailed_cols),
            pd.DataFrame(columns=presence_cols),
        )

    # Build summary table (wide format)
    summary_rows = []
    for term_id, occurrences in common_terms.items():
        # Get description (use first occurrence, should be consistent)
        description = occurrences[0]["description"]
        analyses = sorted([occ["analysis"] for occ in occurrences])
        categories = get_unique_categories(analyses)
        p_values = [occ["p_adjust"] for occ in occurrences]
        counts = [occ["count"] for occ in occurrences]

        summary_rows.append({
            "ID": term_id,
            "Description": description,
            "num_categories": len(categories),
            "categories": "; ".join(sorted(categories)),
            "num_analyses": len(occurrences),
            "analyses": "; ".join(analyses),
            "min_p_adjust": min(p_values),
            "max_p_adjust": max(p_values),
            "total_count": sum(counts),
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values(
        by=["num_categories", "num_analyses", "min_p_adjust"],
        ascending=[False, False, True]
    ).reset_index(drop=True)

    # Build detailed table (long format)
    detailed_rows = []
    for term_id, occurrences in common_terms.items():
        for occ in occurrences:
            detailed_rows.append({
                "ID": term_id,
                "Description": occ["description"],
                "analysis": occ["analysis"],
                "p_adjust": occ["p_adjust"],
                "count": occ["count"],
                "gene_ids": occ["gene_ids"],
            })

    detailed_df = pd.DataFrame(detailed_rows)
    detailed_df = detailed_df.sort_values(
        by=["ID", "analysis"]
    ).reset_index(drop=True)

    # Build presence matrix
    all_analyses = sorted(data.keys())
    presence_rows = []
    for term_id, occurrences in common_terms.items():
        description = occurrences[0]["description"]
        present_analyses = set(occ["analysis"] for occ in occurrences)

        row = {"ID": term_id, "Description": description}
        for analysis in all_analyses:
            row[analysis] = 1 if analysis in present_analyses else 0
        presence_rows.append(row)

    presence_df = pd.DataFrame(presence_rows)
    # Sort by total presence count (descending)
    presence_df["_total"] = presence_df[all_analyses].sum(axis=1)
    presence_df = presence_df.sort_values(by="_total", ascending=False).drop(columns=["_total"]).reset_index(drop=True)

    return summary_df, detailed_df, presence_df


def print_summary_report(
    ko_summary: pd.DataFrame,
    go_summaries: Dict[str, pd.DataFrame],
    min_categories: int,
    p_threshold: float,
) -> None:
    """Print a summary report to console."""
    print("\n" + "=" * 70)
    print("COMMON ENRICHMENT TERMS - SUMMARY REPORT")
    print("=" * 70)
    print(f"\nParameters: min_categories={min_categories}, p_threshold={p_threshold}")
    print("(Only terms spanning multiple analysis categories are shown)")

    # KO summary
    print(f"\n{'─' * 70}")
    print("KEGG KO PATHWAYS")
    print(f"{'─' * 70}")
    if ko_summary.empty:
        print("  No common KO terms found")
    else:
        print(f"  Total common terms: {len(ko_summary)}")
        print(f"\n  Top shared pathways (by number of analyses):")
        for _, row in ko_summary.head(10).iterrows():
            print(f"    {row['ID']}: {row['Description'][:50]}")
            print(f"      Present in {row['num_analyses']} analyses: {row['analyses']}")

    # GO summaries by category
    for category in GO_CATEGORIES:
        category_name = {
            "BP": "Biological Process",
            "CC": "Cellular Component",
            "MF": "Molecular Function",
        }[category]

        print(f"\n{'─' * 70}")
        print(f"GO {category} - {category_name}")
        print(f"{'─' * 70}")

        summary = go_summaries.get(category, pd.DataFrame())
        if summary.empty:
            print("  No common terms found")
        else:
            print(f"  Total common terms: {len(summary)}")
            print(f"\n  Top shared terms (by number of analyses):")
            for _, row in summary.head(5).iterrows():
                desc_truncated = row['Description'][:50] + "..." if len(row['Description']) > 50 else row['Description']
                print(f"    {row['ID']}: {desc_truncated}")
                print(f"      Present in {row['num_analyses']} analyses: {row['analyses']}")

    # Highlight biologically interesting overlaps
    print(f"\n{'─' * 70}")
    print("BIOLOGICALLY INTERESTING OVERLAPS")
    print(f"{'─' * 70}")

    # Look for terms shared between CAFE_expanded and KaKs_positive_tentative
    interesting_combo = ["CAFE_expanded", "KaKs_positive_tentative"]

    if not ko_summary.empty:
        interesting_ko = ko_summary[
            ko_summary["analyses"].apply(
                lambda x: all(a in x for a in interesting_combo)
            )
        ]
        if not interesting_ko.empty:
            print(f"\n  KO pathways in both {' and '.join(interesting_combo)}:")
            for _, row in interesting_ko.iterrows():
                print(f"    • {row['ID']}: {row['Description']}")

    for category in GO_CATEGORIES:
        summary = go_summaries.get(category, pd.DataFrame())
        if not summary.empty:
            interesting_go = summary[
                summary["analyses"].apply(
                    lambda x: all(a in x for a in interesting_combo)
                )
            ]
            if not interesting_go.empty:
                print(f"\n  GO {category} terms in both {' and '.join(interesting_combo)}:")
                for _, row in interesting_go.head(5).iterrows():
                    desc_truncated = row['Description'][:60] + "..." if len(row['Description']) > 60 else row['Description']
                    print(f"    • {row['ID']}: {desc_truncated}")

    print("\n" + "=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Find common enrichment terms across evolutionary analyses"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing enrichment result files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write output files",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.05,
        help="Maximum adjusted p-value threshold (default: 0.05)",
    )
    parser.add_argument(
        "--min-categories",
        type=int,
        default=2,
        help="Minimum number of different analysis categories (CAFE, KaKs, unique) a term must span (default: 2)",
    )

    args = parser.parse_args()

    # Validate input directory
    if not args.input_dir.exists():
        print(f"Error: Input directory does not exist: {args.input_dir}")
        sys.exit(1)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"P-value threshold: {args.p_threshold}")
    print(f"Minimum categories: {args.min_categories}")

    # Process KO enrichment files
    print("\n--- Loading KO enrichment files ---")
    ko_data: Dict[str, pd.DataFrame] = {}
    for analysis_type in ANALYSIS_TYPES:
        print(f"Loading {analysis_type}...")
        df = load_ko_enrichment(args.input_dir, analysis_type)
        if not df.empty:
            ko_data[analysis_type] = df
            print(f"  Loaded {len(df)} terms")

    print("\n--- Finding common KO terms ---")
    ko_summary, ko_detailed, ko_presence = find_common_terms(
        ko_data, args.min_categories, args.p_threshold
    )
    print(f"  Found {len(ko_summary)} common terms")

    # Save KO results
    ko_summary.to_csv(args.output_dir / "KO_common_terms_summary.csv", index=False)
    ko_detailed.to_csv(args.output_dir / "KO_common_terms_detailed.csv", index=False)
    ko_presence.to_csv(args.output_dir / "KO_presence_matrix.csv", index=False)

    # Process GO enrichment files for each category
    go_summaries: Dict[str, pd.DataFrame] = {}

    for category in GO_CATEGORIES:
        print(f"\n--- Loading GO {category} enrichment files ---")
        go_data: Dict[str, pd.DataFrame] = {}
        for analysis_type in ANALYSIS_TYPES:
            print(f"Loading {analysis_type}...")
            df = load_go_enrichment(args.input_dir, analysis_type, category)
            if not df.empty:
                go_data[analysis_type] = df
                print(f"  Loaded {len(df)} terms")

        print(f"\n--- Finding common GO {category} terms ---")
        go_summary, go_detailed, go_presence = find_common_terms(
            go_data, args.min_categories, args.p_threshold
        )
        print(f"  Found {len(go_summary)} common terms")

        go_summaries[category] = go_summary

        # Save GO results
        go_summary.to_csv(args.output_dir / f"GO_{category}_common_terms_summary.csv", index=False)
        go_detailed.to_csv(args.output_dir / f"GO_{category}_common_terms_detailed.csv", index=False)
        go_presence.to_csv(args.output_dir / f"GO_{category}_presence_matrix.csv", index=False)

    # Print summary report
    print_summary_report(ko_summary, go_summaries, args.min_categories, args.p_threshold)

    print(f"\nOutput files written to: {args.output_dir}")
    print("Done!")


if __name__ == "__main__":
    main()
