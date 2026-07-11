#!/usr/bin/env python3
"""CAFE-Synteny Integration Visualization.

Cross-references CAFE5 gene family evolution results with MCScanX synteny
blocks and tandem duplications. Generates publication-quality visualizations
and summary statistics.

Extracted from the CAFE-Synteny Integration section of 6-ultimate-MCScanX.qmd.

Pipeline:
    1. Load CAFE family p-values and per-species change counts
    2. Identify significantly expanded/contracted families for focal species
    3. Map family IDs to individual genes via OrthoFinder Orthogroups.tsv
    4. Convert gene IDs to MCScanX format (strip .t#, add species prefix)
    5. Parse MCScanX .collinearity files to extract synteny gene pairs
    6. Parse MCScanX .tandem files to extract tandem duplicated genes
    7. Calculate overlaps (CAFE genes in synteny / in tandem)
    8. Generate matplotlib bar chart (publication-quality)
    9. Export integration tables as TSV

Usage:
    # Auto-detect all paths from project config:
    python scripts/run_mcscanx_pairwise-visualization.py

    # Explicit paths:
    python scripts/run_mcscanx_pairwise-visualization.py \\
        --mcscanx-dirs ./Outputs_set3/MCScanX_Results/dm_dr \\
                       ./Outputs_set3/MCScanX_Results/dm_al \\
        --focal-species Diptychus_maculatus \\
        --cafe-family-results ./Outputs_set3/CAFE/cafe_results/single_lambda/Base_family_results.txt \\
        --cafe-change-tab ./Outputs_set3/CAFE/cafe_results/single_lambda/Base_change.tab \\
        --orthogroups-tsv ./Outputs/OrthoFinder/Results_Mar03/Orthogroups/Orthogroups.tsv \\
        --output-dir ./Outputs_set3/MCScanX_Results/integration \\
        --figures-dir ./Notebooks/Figures_set3

Author: Extracted from 6-ultimate-MCScanX.qmd (CAFE-Synteny Integration section)
"""

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend for headless environments
import matplotlib.pyplot as plt

# Embed TrueType fonts (type 42) instead of the matplotlib default Type 3.
# Type 3 fonts are rendered as outlines that Adobe Illustrator cannot edit;
# type 42 keeps text as live, selectable, replaceable text objects. Set at
# module import so the plotting helpers below inherit it even when imported
# by a standalone wrapper.
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

import numpy as np
import pandas as pd

PYCIRCLIZE_AVAILABLE = False
try:
    from pycirclize import Circos

    PYCIRCLIZE_AVAILABLE = True
except ImportError:
    pass

PYGENOMEVIZ_AVAILABLE = False
try:
    from pygenomeviz import GenomeViz

    PYGENOMEVIZ_AVAILABLE = True
except ImportError:
    pass


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def get_species_prefix(species_name: str) -> str:
    """Derive a 2-letter MCScanX prefix from a species name.

    Takes the first letter of genus + first letter of species epithet,
    lowercased.

    Args:
        species_name: e.g. 'Diptychus_maculatus'

    Returns:
        Prefix string, e.g. 'dm'
    """
    parts = species_name.split("_")
    if len(parts) < 2:
        print(
            f"Error: Species name must be 'Genus_species' format: {species_name}",
            file=sys.stderr,
        )
        sys.exit(1)
    return (parts[0][0] + parts[1][0]).lower()


def get_display_name(species_name: str, abbreviated: bool = False) -> str:
    """Convert species name to display format.

    Args:
        species_name: e.g. 'Diptychus_maculatus'
        abbreviated: If True, abbreviate genus (e.g. 'D. maculatus')

    Returns:
        Display name string
    """
    parts = species_name.split("_")
    genus, epithet = parts[0], parts[1]
    if abbreviated:
        return f"{genus[0]}. {epithet}"
    return f"{genus} {epithet}"


def get_abbreviation(species_name: str) -> str:
    """Get 3-letter abbreviation from species name.

    Args:
        species_name: e.g. 'Diptychus_maculatus'

    Returns:
        3-letter abbreviation, e.g. 'Dip'
    """
    return species_name.split("_")[0][:3]


def convert_to_mcscanx_format(gene_id: str, species_prefix: str) -> str:
    """Convert OrthoFinder gene ID to MCScanX-compatible format.

    OrthoFinder format: g123.t1 (gene.transcript)
    MCScanX format: dm_g123 (species_gene, no transcript suffix)

    Args:
        gene_id: OrthoFinder gene ID (e.g. 'g123.t1')
        species_prefix: Species prefix (e.g. 'dm')

    Returns:
        MCScanX-compatible gene ID (e.g. 'dm_g123')
    """
    gene_base = re.sub(r"\.t\d+$", "", gene_id)
    return f"{species_prefix}_{gene_base}"


def natural_sort_key(s: str) -> List[object]:
    """Generate natural-sort key for chromosome strings.

    Example:
        dm1, dm2, ..., dm10 (instead of dm1, dm10, dm2)

    Args:
        s: Chromosome name.

    Returns:
        Mixed list of text and numeric components for sorting.
    """
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"([0-9]+)", s)]


def species_to_latex(species_name: Optional[str], fallback_tag: str) -> str:
    """Format species name as italicized LaTeX label.

    Args:
        species_name: Genus_species string when available.
        fallback_tag: Prefix or short tag used when name is unavailable.

    Returns:
        LaTeX italic label string.
    """
    if species_name and "_" in species_name:
        parts = species_name.split("_")
        if len(parts) >= 2:
            return f"$\\it{{{parts[0][0]}.\\ {parts[1]}}}$"
    return f"$\\it{{{fallback_tag.upper()}}}$"


def load_chromosome_sizes_from_bed(bed_path: Path) -> Dict[str, int]:
    """Load chromosome sizes from MCScanX BED file.

    MCScanX BED format: chr, start, end, gene_id (tab-separated).
    Chromosome size is inferred as max(end) per chromosome.

    Args:
        bed_path: BED path.

    Returns:
        Mapping of chromosome name to chromosome size (bp).
    """
    if not bed_path.exists():
        print(f"  Warning: BED file not found: {bed_path}", file=sys.stderr)
        return {}

    bed_df = pd.read_csv(
        bed_path,
        sep="\t",
        names=["chr", "start", "end", "gene"],
        usecols=[0, 1, 2, 3],
    )
    if bed_df.empty:
        return {}
    return bed_df.groupby("chr")["end"].max().astype(int).to_dict()


def calculate_gene_density(
    bed_path: Path, window_size: int = 1_000_000
) -> Dict[str, List[Tuple[int, int]]]:
    """Calculate gene density per chromosome from BED file.

    Args:
        bed_path: MCScanX BED file path.
        window_size: Window size in bp (default: 1 Mb).

    Returns:
        Dict mapping chromosome to list of (window_start, gene_count).
    """
    if not bed_path.exists():
        print(f"  Warning: BED file not found: {bed_path}", file=sys.stderr)
        return {}

    bed_df = pd.read_csv(
        bed_path,
        sep="\t",
        names=["chr", "start", "end", "gene"],
        usecols=[0, 1, 2, 3],
    )
    if bed_df.empty:
        return {}

    density: Dict[str, List[Tuple[int, int]]] = {}
    for chrom in bed_df["chr"].unique():
        chrom_df = bed_df[bed_df["chr"] == chrom].copy()
        max_pos = int(chrom_df["end"].max())
        midpoints = (chrom_df["start"] + chrom_df["end"]) / 2

        windows: List[Tuple[int, int]] = []
        for window_start in range(0, max_pos + window_size, window_size):
            window_end = window_start + window_size
            count = int(((midpoints >= window_start) & (midpoints < window_end)).sum())
            windows.append((window_start, count))
        density[chrom] = windows

    return density


def parse_collinearity_for_dotplot(
    collinearity_path: Path,
    bed1_path: Path,
    bed2_path: Path,
) -> pd.DataFrame:
    """Parse MCScanX collinearity file into Oxford plot coordinates.

    Args:
        collinearity_path: Path to .collinearity file.
        bed1_path: BED for species 1.
        bed2_path: BED for species 2.

    Returns:
        DataFrame with gene1, chr1, pos1, gene2, chr2, pos2, orientation.
    """
    columns = ["gene1", "chr1", "pos1", "gene2", "chr2", "pos2", "orientation"]
    if (
        not collinearity_path.exists()
        or not bed1_path.exists()
        or not bed2_path.exists()
    ):
        print(
            "  Warning: Missing input for dotplot parsing "
            f"({collinearity_path}, {bed1_path}, {bed2_path})",
            file=sys.stderr,
        )
        return pd.DataFrame(columns=columns)

    bed1_df = pd.read_csv(
        bed1_path,
        sep="\t",
        names=["chr", "start", "end", "gene"],
        usecols=[0, 1, 2, 3],
    )
    bed2_df = pd.read_csv(
        bed2_path,
        sep="\t",
        names=["chr", "start", "end", "gene"],
        usecols=[0, 1, 2, 3],
    )

    pos1 = {
        row["gene"]: {
            "chr": row["chr"],
            "start": int(row["start"]),
            "end": int(row["end"]),
        }
        for _, row in bed1_df.iterrows()
    }
    pos2 = {
        row["gene"]: {
            "chr": row["chr"],
            "start": int(row["start"]),
            "end": int(row["end"]),
        }
        for _, row in bed2_df.iterrows()
    }

    pairs = []
    current_chr_pair = None
    orientation = "plus"

    with open(collinearity_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("## Alignment"):
                parts = line.split()
                current_chr_pair = next((p for p in parts if "&" in p), None)
                orientation = "plus" if "plus" in line else "minus"
            elif line and not line.startswith("#") and current_chr_pair:
                parts = line.split("\t")
                if len(parts) >= 3:
                    gene1 = parts[1].strip()
                    gene2 = parts[2].strip()
                    if gene1 in pos1 and gene2 in pos2:
                        pairs.append(
                            {
                                "gene1": gene1,
                                "chr1": pos1[gene1]["chr"],
                                "pos1": (pos1[gene1]["start"] + pos1[gene1]["end"])
                                // 2,
                                "gene2": gene2,
                                "chr2": pos2[gene2]["chr"],
                                "pos2": (pos2[gene2]["start"] + pos2[gene2]["end"])
                                // 2,
                                "orientation": orientation,
                            }
                        )

    if not pairs:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(pairs)


def parse_synteny_blocks_for_circos(
    collinearity_path: Path,
    bed1_path: Path,
    bed2_path: Path,
    chr1_prefix: str,
    chr2_prefix: str,
    min_genes: int = 10,
) -> List[Dict[str, object]]:
    """Parse MCScanX collinearity into synteny blocks for circos ribbons.

    Args:
        collinearity_path: Path to .collinearity file.
        bed1_path: BED for species 1.
        bed2_path: BED for species 2.
        chr1_prefix: Prefix for species 1 chromosomes.
        chr2_prefix: Prefix for species 2 chromosomes.
        min_genes: Minimum genes per block to keep.

    Returns:
        List of synteny block dicts with coordinates and orientation.
    """
    if (
        not collinearity_path.exists()
        or not bed1_path.exists()
        or not bed2_path.exists()
    ):
        print(
            "  Warning: Missing input for circos synteny parsing "
            f"({collinearity_path}, {bed1_path}, {bed2_path})",
            file=sys.stderr,
        )
        return []

    bed1_df = pd.read_csv(
        bed1_path,
        sep="\t",
        names=["chr", "start", "end", "gene"],
        usecols=[0, 1, 2, 3],
    )
    bed2_df = pd.read_csv(
        bed2_path,
        sep="\t",
        names=["chr", "start", "end", "gene"],
        usecols=[0, 1, 2, 3],
    )

    pos1 = {
        row["gene"]: {
            "chr": row["chr"],
            "start": int(row["start"]),
            "end": int(row["end"]),
        }
        for _, row in bed1_df.iterrows()
    }
    pos2 = {
        row["gene"]: {
            "chr": row["chr"],
            "start": int(row["start"]),
            "end": int(row["end"]),
        }
        for _, row in bed2_df.iterrows()
    }

    blocks: List[Dict[str, object]] = []
    current_block: Optional[Dict[str, object]] = None

    with open(collinearity_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("## Alignment"):
                if current_block and len(current_block["genes"]) >= min_genes:  # type: ignore[index]
                    genes = current_block["genes"]  # type: ignore[index]
                    starts1 = [g[0] for g in genes if g[0] is not None]
                    ends1 = [g[1] for g in genes if g[1] is not None]
                    starts2 = [g[2] for g in genes if g[2] is not None]
                    ends2 = [g[3] for g in genes if g[3] is not None]
                    if starts1 and ends1 and starts2 and ends2:
                        blocks.append(
                            {
                                "chr1": current_block["chr1"],
                                "start1": int(min(starts1)),
                                "end1": int(max(ends1)),
                                "chr2": current_block["chr2"],
                                "start2": int(min(starts2)),
                                "end2": int(max(ends2)),
                                "n_genes": len(genes),
                                "orientation": current_block["orientation"],
                            }
                        )

                chr_pair_parts = next(
                    (
                        part.split("&")
                        for part in line.split()
                        if "&" in part and len(part.split("&")) == 2
                    ),
                    None,
                )
                orientation = "minus" if "minus" in line else "plus"

                if chr_pair_parts:
                    chr1, chr2 = chr_pair_parts
                    if chr1.startswith(chr1_prefix) and chr2.startswith(chr2_prefix):
                        current_block = {
                            "chr1": chr1,
                            "chr2": chr2,
                            "orientation": orientation,
                            "genes": [],
                        }
                    else:
                        current_block = None
                else:
                    current_block = None
            elif line and not line.startswith("#") and current_block:
                parts = line.split("\t")
                if len(parts) >= 3:
                    gene1 = parts[1].strip()
                    gene2 = parts[2].strip()
                    start1 = pos1.get(gene1, {}).get("start")
                    end1 = pos1.get(gene1, {}).get("end")
                    start2 = pos2.get(gene2, {}).get("start")
                    end2 = pos2.get(gene2, {}).get("end")
                    current_block["genes"].append((start1, end1, start2, end2))  # type: ignore[index]

    if current_block and len(current_block["genes"]) >= min_genes:  # type: ignore[index]
        genes = current_block["genes"]  # type: ignore[index]
        starts1 = [g[0] for g in genes if g[0] is not None]
        ends1 = [g[1] for g in genes if g[1] is not None]
        starts2 = [g[2] for g in genes if g[2] is not None]
        ends2 = [g[3] for g in genes if g[3] is not None]
        if starts1 and ends1 and starts2 and ends2:
            blocks.append(
                {
                    "chr1": current_block["chr1"],
                    "start1": int(min(starts1)),
                    "end1": int(max(ends1)),
                    "chr2": current_block["chr2"],
                    "start2": int(min(starts2)),
                    "end2": int(max(ends2)),
                    "n_genes": len(genes),
                    "orientation": current_block["orientation"],
                }
            )

    return blocks


def create_genome_circos(
    chrom_sizes: Dict[str, int],
    gene_density: Dict[str, List[Tuple[int, int]]],
    species_label: str,
    output_path: Path,
) -> bool:
    """Create a single-genome circos plot with chromosome and gene-density tracks."""
    if not PYCIRCLIZE_AVAILABLE:
        print(
            "  Warning: pyCirclize not installed; skipping genome circos.",
            file=sys.stderr,
        )
        return False
    if not chrom_sizes:
        print(
            "  Warning: Empty chromosome sizes; skipping genome circos.",
            file=sys.stderr,
        )
        return False

    sorted_chroms = dict(
        sorted(chrom_sizes.items(), key=lambda x: natural_sort_key(x[0]))
    )
    circos = Circos(sorted_chroms, space=2)

    n_chroms = len(sorted_chroms)
    colors = plt.cm.tab20(np.linspace(0, 1, min(n_chroms, 20)))
    if n_chroms > 20:
        colors = plt.cm.nipy_spectral(np.linspace(0, 1, n_chroms))

    all_densities: List[int] = []
    for chrom in sorted_chroms.keys():
        if chrom in gene_density:
            all_densities.extend([d[1] for d in gene_density[chrom]])
    max_density = max(all_densities) if all_densities else 1

    for idx, sector in enumerate(circos.sectors):
        color = colors[idx % len(colors)]
        outer = sector.add_track((95, 100))
        outer.axis(fc=color, alpha=0.8)
        if idx % 2 == 0:
            outer.text(sector.name, fontsize=6, r=107)

        density_track = sector.add_track((80, 93))
        density_track.axis(fc="white", ec="gray", lw=0.5)
        if sector.name in gene_density:
            density_data = gene_density[sector.name]
            x_vals = np.array(
                [min(d[0] + 500000, sector.size - 1) for d in density_data]
            )
            y_vals = np.array([d[1] for d in density_data])
            y_norm = (y_vals / max_density) * 10
            density_track.fill_between(x_vals, y_norm, fc="steelblue", alpha=0.7)

    circos.text(
        f"{species_label}\\nGenome Overview ({len(sorted_chroms)} chromosomes)",
        r=40,
        size=10,
    )

    fig = circos.plotfig()
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="steelblue", alpha=0.7, label="Gene density (genes per Mb)"),
        Patch(
            facecolor="lightgray",
            alpha=0.8,
            edgecolor="black",
            linewidth=0.5,
            label="Chromosome",
        ),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower right",
        bbox_to_anchor=(0.95, 0.05),
        fontsize=8,
        frameon=True,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return True


def create_synteny_circos(
    focal_chroms: Dict[str, int],
    comp_chroms: Dict[str, int],
    synteny_blocks: List[Dict[str, object]],
    focal_tag: str,
    comp_tag: str,
    focal_label: str,
    comp_label: str,
    output_path: Path,
) -> bool:
    """Create dual-genome circos with synteny ribbons between species."""
    if not PYCIRCLIZE_AVAILABLE:
        print(
            "  Warning: pyCirclize not installed; skipping synteny circos.",
            file=sys.stderr,
        )
        return False
    if not focal_chroms or not comp_chroms:
        print(
            "  Warning: Missing chromosome sizes; skipping synteny circos.",
            file=sys.stderr,
        )
        return False

    focal_filtered = {
        k: v
        for k, v in sorted(focal_chroms.items(), key=lambda x: natural_sort_key(x[0]))
        if v > 5_000_000
    }
    comp_filtered = {
        k: v
        for k, v in sorted(comp_chroms.items(), key=lambda x: natural_sort_key(x[0]))
        if v > 5_000_000
    }
    if not focal_filtered or not comp_filtered:
        print(
            "  Warning: Chromosomes filtered out (<5 Mb); skipping synteny circos.",
            file=sys.stderr,
        )
        return False

    dual_sectors = {f"{focal_tag}_{k}": v for k, v in focal_filtered.items()}
    dual_sectors.update({f"{comp_tag}_{k}": v for k, v in comp_filtered.items()})
    circos = Circos(dual_sectors, space=1)

    focal_color = "#3498db"
    comp_color = "#e74c3c"
    sector_order = list(dual_sectors.keys())

    for sector in circos.sectors:
        is_focal = sector.name.startswith(focal_tag)
        color = focal_color if is_focal else comp_color
        outer = sector.add_track((90, 100))
        outer.axis(fc=color, alpha=0.7)

        idx = sector_order.index(sector.name)
        display_name = sector.name.replace(f"{focal_tag}_", "").replace(
            f"{comp_tag}_", ""
        )
        if idx % 3 == 0:
            outer.text(display_name, fontsize=5, r=105)

        inner = sector.add_track((75, 88))
        inner.axis(fc="white", ec="gray", lw=0.5)

    for block in synteny_blocks:
        focal_sector = f"{focal_tag}_{block['chr1']}"
        comp_sector = f"{comp_tag}_{block['chr2']}"
        if focal_sector in dual_sectors and comp_sector in dual_sectors:
            ribbon_color = "#3498db" if block["orientation"] == "plus" else "#e74c3c"
            circos.link(
                (focal_sector, int(block["start1"]), int(block["end1"])),
                (comp_sector, int(block["start2"]), int(block["end2"])),
                color=ribbon_color,
                alpha=0.3,
            )

    circos.text(f"Synteny: {focal_label} vs {comp_label}\n(MCScanX)", r=35, size=9)

    fig = circos.plotfig()
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#3498db", alpha=0.7, label="Plus orientation"),
        Patch(facecolor="#e74c3c", alpha=0.7, label="Minus orientation"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower right",
        bbox_to_anchor=(0.95, 0.05),
        fontsize=8,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return True


def create_oxford_plot(
    pairs_df: pd.DataFrame,
    chr1_sizes: Dict[str, int],
    chr2_sizes: Dict[str, int],
    species1_name: str,
    species2_name: str,
    output_path: Path,
    chr1_prefix: str,
    chr2_prefix: str,
) -> bool:
    """Create Oxford (dot) plot from syntenic pair coordinates."""
    if pairs_df.empty or not chr1_sizes or not chr2_sizes:
        print("  Warning: Missing Oxford plot inputs; skipping.", file=sys.stderr)
        return False

    sorted_chr1 = sorted(
        [c for c in chr1_sizes.keys() if c.startswith(chr1_prefix)],
        key=natural_sort_key,
    )
    sorted_chr2 = sorted(
        [c for c in chr2_sizes.keys() if c.startswith(chr2_prefix)],
        key=natural_sort_key,
    )
    if not sorted_chr1 or not sorted_chr2:
        print(
            "  Warning: No matching chromosome prefixes for Oxford plot.",
            file=sys.stderr,
        )
        return False

    chr1_cumsum: Dict[str, int] = {}
    cumsum = 0
    for chrom in sorted_chr1:
        chr1_cumsum[chrom] = cumsum
        cumsum += int(chr1_sizes[chrom])
    total_size1 = cumsum

    chr2_cumsum: Dict[str, int] = {}
    cumsum = 0
    for chrom in sorted_chr2:
        chr2_cumsum[chrom] = cumsum
        cumsum += int(chr2_sizes[chrom])
    total_size2 = cumsum

    plot_df = pairs_df.copy()
    plot_df["genome_pos1"] = plot_df.apply(
        lambda row: chr1_cumsum.get(row["chr1"], 0) + row["pos1"]
        if row["chr1"] in chr1_cumsum
        else None,
        axis=1,
    )
    plot_df["genome_pos2"] = plot_df.apply(
        lambda row: chr2_cumsum.get(row["chr2"], 0) + row["pos2"]
        if row["chr2"] in chr2_cumsum
        else None,
        axis=1,
    )
    plot_df = plot_df.dropna(subset=["genome_pos1", "genome_pos2"])

    if plot_df.empty:
        print("  Warning: No valid coordinates after filtering; skipping Oxford plot.")
        return False

    fig, ax = plt.subplots(figsize=(12, 12))
    colors = (
        plot_df["orientation"]
        .map({"plus": "#3498db", "minus": "#e74c3c"})
        .fillna("gray")
    )
    ax.scatter(
        plot_df["genome_pos1"] / 1e6,
        plot_df["genome_pos2"] / 1e6,
        c=colors,
        s=1,
        alpha=0.5,
        rasterized=True,
    )

    for i, chrom in enumerate(sorted_chr1):
        x = chr1_cumsum[chrom] / 1e6
        ax.axvline(x, color="gray", linewidth=0.3, alpha=0.5)
        if i % 5 == 0:
            mid_x = (chr1_cumsum[chrom] + chr1_sizes[chrom] / 2) / 1e6
            ax.text(
                mid_x,
                -total_size2 / 1e6 * 0.02,
                chrom.replace(chr1_prefix, ""),
                ha="center",
                va="top",
                fontsize=8,
                rotation=45,
            )

    for i, chrom in enumerate(sorted_chr2):
        y = chr2_cumsum[chrom] / 1e6
        ax.axhline(y, color="gray", linewidth=0.3, alpha=0.5)
        if i % 5 == 0:
            mid_y = (chr2_cumsum[chrom] + chr2_sizes[chrom] / 2) / 1e6
            ax.text(
                -total_size1 / 1e6 * 0.02,
                mid_y,
                chrom.replace(chr2_prefix, ""),
                ha="right",
                va="center",
                fontsize=8,
            )

    ax.set_xlabel(f"{species1_name} (Mb)", fontsize=14, fontweight="bold")
    ax.set_ylabel(f"{species2_name} (Mb)", fontsize=14, fontweight="bold")
    ax.set_title(
        f"Oxford Plot: {species1_name} vs {species2_name}\\nSyntenic Gene Pairs",
        fontsize=16,
        fontweight="bold",
    )
    ax.set_xlim(0, total_size1 / 1e6)
    ax.set_ylim(0, total_size2 / 1e6)

    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#3498db",
            markersize=8,
            label="Plus strand",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#e74c3c",
            markersize=8,
            label="Minus strand (inverted)",
        ),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=10)
    ax.tick_params(axis="both", labelsize=12)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return True


def create_linear_ideogram(
    chrom_sizes: Dict[str, int],
    gene_density: Dict[str, List[Tuple[int, int]]],
    species_label: str,
    output_path: Path,
) -> bool:
    """Create linear chromosome ideogram with gene-density segments."""
    if not chrom_sizes or not gene_density:
        print("  Warning: Missing ideogram inputs; skipping.", file=sys.stderr)
        return False

    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap

    chromosomes = dict(
        sorted(chrom_sizes.items(), key=lambda x: natural_sort_key(x[0]))
    )
    gene_density_ideogram: Dict[str, np.ndarray] = {}
    for chrom in chromosomes.keys():
        if chrom in gene_density:
            density_data = gene_density[chrom]
            gene_density_ideogram[chrom] = np.array([d[1] for d in density_data])
        else:
            n_windows = max(1, chromosomes[chrom] // 1_000_000)
            gene_density_ideogram[chrom] = np.zeros(n_windows)

    fig_height = max(10, len(chromosomes) * 0.72)
    fig, axes = plt.subplots(
        len(chromosomes),
        1,
        figsize=(14, fig_height),
        gridspec_kw={"hspace": 0.2},
    )

    density_cmap = LinearSegmentedColormap.from_list(
        "density", ["white", "steelblue", "darkblue"]
    )
    max_size = max(chromosomes.values())
    all_density_values: List[float] = []
    for chrom in chromosomes.keys():
        all_density_values.extend(gene_density_ideogram[chrom])
    max_density_val = max(all_density_values) if all_density_values else 1

    for idx, (chrom, size) in enumerate(chromosomes.items()):
        ax = axes[idx] if len(chromosomes) > 1 else axes
        scale = size / max_size

        chrom_height = 0.4
        chrom_rect = mpatches.FancyBboxPatch(
            (0, 0.3),
            scale,
            chrom_height,
            boxstyle="round,pad=0.01,rounding_size=0.02",
            facecolor="lightgray",
            edgecolor="black",
            linewidth=1,
        )
        ax.add_patch(chrom_rect)

        density = gene_density_ideogram[chrom]
        n_windows = len(density)
        if n_windows > 0:
            window_width = scale / n_windows
            for i, d in enumerate(density):
                x_start = i * window_width
                color = density_cmap(d / max_density_val)
                rect = mpatches.Rectangle(
                    (x_start, 0.3),
                    window_width,
                    chrom_height * 0.5,
                    facecolor=color,
                    edgecolor="none",
                )
                ax.add_patch(rect)

        ax.text(
            -0.02, 0.5, chrom, ha="right", va="center", fontsize=9, fontweight="bold"
        )
        ax.text(
            scale + 0.01,
            0.5,
            f"{size / 1e6:.0f} Mb",
            ha="left",
            va="center",
            fontsize=8,
        )
        ax.set_xlim(-0.08, 1.1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    fig.suptitle(
        f"{species_label}: Linear Chromosome Ideogram\\nwith Gene Density",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )
    legend_elements = [
        mpatches.Patch(facecolor="lightgray", edgecolor="black", label="Chromosome"),
        mpatches.Patch(
            facecolor="steelblue", edgecolor="none", label="Gene density (high)"
        ),
        mpatches.Patch(facecolor="white", edgecolor="none", label="Gene density (low)"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=3,
        frameon=True,
        fontsize=9,
        bbox_to_anchor=(0.5, 0.01),
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.96])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_cafe_results(
    family_results_path: Path,
    change_tab_path: Path,
    focal_species: str,
    pvalue_threshold: float = 0.05,
) -> Tuple[List[str], List[str]]:
    """Load CAFE5 results and identify significant families for focal species.

    Args:
        family_results_path: Path to Base_family_results.txt
        change_tab_path: Path to Base_change.tab
        focal_species: Focal species name (e.g. 'Diptychus_maculatus')
        pvalue_threshold: Significance threshold (default 0.05)

    Returns:
        Tuple of (expanded_family_ids, contracted_family_ids)
    """
    print("Loading CAFE results...")

    # Load family-level p-values
    cafe_family = pd.read_csv(family_results_path, sep="\t")
    # First column header may start with '#'; standardize
    cafe_family.columns = ["FamilyID", "pvalue", "Significant"]

    # Load per-family per-species changes
    cafe_changes = pd.read_csv(change_tab_path, sep="\t")

    # Find the column for focal species
    focal_col = None
    for col in cafe_changes.columns:
        if focal_species in col:
            focal_col = col
            break

    if focal_col is None:
        print(
            f"  Error: Focal species '{focal_species}' not found in change tab columns:",
            file=sys.stderr,
        )
        print(f"  Available: {list(cafe_changes.columns)}", file=sys.stderr)
        sys.exit(1)

    abbreviation = get_abbreviation(focal_species)
    change_col = f"{abbreviation}_Change"

    # Merge
    cafe_merged = cafe_family.copy()
    cafe_merged[change_col] = cafe_changes[focal_col].values

    # Filter significant families
    expanded = cafe_merged[
        (cafe_merged["pvalue"] < pvalue_threshold) & (cafe_merged[change_col] > 0)
    ]["FamilyID"].tolist()

    contracted = cafe_merged[
        (cafe_merged["pvalue"] < pvalue_threshold) & (cafe_merged[change_col] < 0)
    ]["FamilyID"].tolist()

    display = get_display_name(focal_species, abbreviated=True)
    print(f"  CAFE-expanded families in {display}: {len(expanded)}")
    print(f"  CAFE-contracted families in {display}: {len(contracted)}")

    return expanded, contracted


def load_orthogroup_genes(
    orthogroups_path: Path,
    family_ids: List[str],
    focal_species: str,
) -> List[Tuple[str, str]]:
    """Extract focal species genes for a list of orthogroup families.

    Args:
        orthogroups_path: Path to Orthogroups.tsv
        family_ids: List of orthogroup family IDs (e.g. ['OG0000001', ...])
        focal_species: Focal species column name in Orthogroups.tsv

    Returns:
        List of (gene_id, family_id) tuples
    """
    orthogroups_df = pd.read_csv(orthogroups_path, sep="\t")

    if focal_species not in orthogroups_df.columns:
        print(
            f"  Warning: '{focal_species}' column not found in Orthogroups.tsv",
            file=sys.stderr,
        )
        print(
            f"  Available columns: {list(orthogroups_df.columns)[:10]}...",
            file=sys.stderr,
        )
        return []

    genes = []
    for family in family_ids:
        family_row = orthogroups_df[orthogroups_df["Orthogroup"] == family]
        if not family_row.empty and pd.notna(family_row[focal_species].values[0]):
            gene_list = str(family_row[focal_species].values[0]).split(", ")
            genes.extend([(gene, family) for gene in gene_list])

    return genes


def parse_collinearity_file(collinearity_path: Path) -> pd.DataFrame:
    """Parse MCScanX collinearity file to extract synteny blocks.

    Format:
        ## Alignment 0: score=17535.0 e_value=0 N=355 dm1&dm2 plus
          0-  0:	dm_g6	dm_g483	  8e-97

    Args:
        collinearity_path: Path to .collinearity file

    Returns:
        DataFrame with columns: Block, Gene1, Gene2
    """
    blocks = []
    current_block = None

    if not collinearity_path.exists():
        print(f"  Collinearity file not found: {collinearity_path}", file=sys.stderr)
        return pd.DataFrame()

    with open(collinearity_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("## Alignment"):
                parts = line.split()
                if len(parts) >= 4:
                    current_block = parts[2]  # Block number (e.g. "0:")
            elif line and not line.startswith("#") and current_block:
                parts = line.split("\t")
                if len(parts) >= 3:
                    gene1 = parts[1].split(":")[0] if ":" in parts[1] else parts[1]
                    gene2 = parts[2].split(":")[0] if ":" in parts[2] else parts[2]
                    blocks.append(
                        {
                            "Block": current_block,
                            "Gene1": gene1.strip(),
                            "Gene2": gene2.strip(),
                        }
                    )

    return pd.DataFrame(blocks)


def parse_tandem_file(tandem_path: Path) -> Set[str]:
    """Parse MCScanX tandem duplication file.

    Each line contains comma-separated gene pairs that are tandem duplicates.

    Args:
        tandem_path: Path to .tandem file

    Returns:
        Set of all tandem-duplicated gene IDs
    """
    genes = set()
    if not tandem_path.exists():
        return genes

    with open(tandem_path, "r") as f:
        for line in f:
            parts = line.strip().split(",")
            genes.update(g.strip() for g in parts if g.strip())

    return genes


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def run_cafe_synteny_integration(
    mcscanx_dirs: List[Path],
    focal_species: str,
    cafe_family_results: Path,
    cafe_change_tab: Path,
    orthogroups_tsv: Path,
    output_dir: Path,
    figures_dir: Path,
    pvalue_threshold: float = 0.05,
    plots: str = "all",
    prefix_species_map: Optional[Dict[str, str]] = None,
) -> None:
    """Run the complete CAFE-Synteny integration analysis.

    Args:
        mcscanx_dirs: List of MCScanX result directories (each containing
                      {name}.collinearity, {name}.tandem, etc.)
        focal_species: Focal species name (e.g. 'Diptychus_maculatus')
        cafe_family_results: Path to Base_family_results.txt
        cafe_change_tab: Path to Base_change.tab
        orthogroups_tsv: Path to Orthogroups.tsv
        output_dir: Directory for integration TSV outputs
        figures_dir: Directory for figure outputs
        pvalue_threshold: CAFE significance threshold
        plots: Plot mode selector (all, cafe, oxford, circos, ideogram)
        prefix_species_map: Optional mapping of prefix -> Genus_species
    """
    focal_prefix = get_species_prefix(focal_species)
    focal_display = get_display_name(focal_species, abbreviated=True)

    print(f"\n{'=' * 60}")
    print(f"  CAFE-Synteny Integration Analysis")
    print(f"{'=' * 60}")
    print(
        f"  Focal species: {get_display_name(focal_species)} ({focal_prefix.upper()})"
    )
    print(f"  MCScanX dirs: {len(mcscanx_dirs)}")
    print(f"  P-value threshold: {pvalue_threshold}")
    print(f"{'=' * 60}\n")

    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1: Load CAFE results
    # ------------------------------------------------------------------
    print("--- Step 1/6: Load CAFE Results ---\n")
    cafe_expanded, cafe_contracted = load_cafe_results(
        cafe_family_results, cafe_change_tab, focal_species, pvalue_threshold
    )

    if not cafe_expanded and not cafe_contracted:
        print("  No significant CAFE families found. Aborting.", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 2: Map families to genes via OrthoFinder
    # ------------------------------------------------------------------
    print("\n--- Step 2/6: Map Families to Genes ---\n")
    print(f"Loading OrthoFinder orthogroups from: {orthogroups_tsv}")

    expanded_genes = load_orthogroup_genes(
        orthogroups_tsv, cafe_expanded, focal_species
    )
    contracted_genes = load_orthogroup_genes(
        orthogroups_tsv, cafe_contracted, focal_species
    )

    print(f"  Genes in CAFE-expanded families: {len(expanded_genes)}")
    print(f"  Genes in CAFE-contracted families: {len(contracted_genes)}")

    # ------------------------------------------------------------------
    # Step 3: Convert gene IDs to MCScanX format
    # ------------------------------------------------------------------
    print("\n--- Step 3/6: Convert Gene IDs ---\n")

    expanded_genes_mcscanx = [
        (convert_to_mcscanx_format(gene, focal_prefix), family)
        for gene, family in expanded_genes
    ]
    contracted_genes_mcscanx = [
        (convert_to_mcscanx_format(gene, focal_prefix), family)
        for gene, family in contracted_genes
    ]

    # Save MCScanX-compatible gene lists
    expanded_df = pd.DataFrame(
        expanded_genes_mcscanx, columns=["Gene_ID_MCScanX", "Orthogroup"]
    )
    contracted_df = pd.DataFrame(
        contracted_genes_mcscanx, columns=["Gene_ID_MCScanX", "Orthogroup"]
    )
    expanded_df["Gene_ID_Original"] = [g[0] for g in expanded_genes]
    contracted_df["Gene_ID_Original"] = [g[0] for g in contracted_genes]

    expanded_tsv = output_dir / "cafe_expanded_genes_mcscanx.tsv"
    contracted_tsv = output_dir / "cafe_contracted_genes_mcscanx.tsv"
    expanded_df.to_csv(expanded_tsv, sep="\t", index=False)
    contracted_df.to_csv(contracted_tsv, sep="\t", index=False)

    print(f"  Saved {expanded_tsv.name}: {len(expanded_df)} genes")
    print(f"  Saved {contracted_tsv.name}: {len(contracted_df)} genes")

    if expanded_genes:
        orig = expanded_genes[0][0]
        conv = expanded_genes_mcscanx[0][0]
        print(f"  Format example: {orig} -> {conv}")

    # ------------------------------------------------------------------
    # Step 4: Parse synteny blocks from all comparisons
    # ------------------------------------------------------------------
    print("\n--- Step 4/6: Parse Synteny Blocks ---\n")

    all_synteny_dfs = []
    for mcscanx_dir in mcscanx_dirs:
        # Auto-detect comparison name from directory
        comp_name = mcscanx_dir.name
        col_file = mcscanx_dir / f"{comp_name}.collinearity"

        if not col_file.exists():
            # Try to find any .collinearity file in the directory
            col_files = list(mcscanx_dir.glob("*.collinearity"))
            if col_files:
                col_file = col_files[0]
                comp_name = col_file.stem
            else:
                print(
                    f"  No .collinearity file found in {mcscanx_dir}", file=sys.stderr
                )
                continue

        df = parse_collinearity_file(col_file)
        if not df.empty:
            df["Comparison"] = comp_name
            all_synteny_dfs.append(df)
            print(f"  Loaded {len(df)} gene pairs from {col_file.name}")

    if all_synteny_dfs:
        synteny_df = pd.concat(all_synteny_dfs, ignore_index=True)
        print(f"\n  Total synteny gene pairs: {len(synteny_df):,}")
    else:
        print("  No synteny data found. Cannot proceed.", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 5: Calculate CAFE-Synteny overlap
    # ------------------------------------------------------------------
    print("\n--- Step 5/6: CAFE-Synteny Overlap Analysis ---\n")

    expanded_ids = [g[0] for g in expanded_genes_mcscanx]
    contracted_ids = [g[0] for g in contracted_genes_mcscanx]

    # All genes in synteny blocks
    synteny_genes = set(synteny_df["Gene1"].tolist() + synteny_df["Gene2"].tolist())

    # Filter to focal species genes only (starting with focal prefix)
    focal_synteny_genes = {g for g in synteny_genes if g.startswith(f"{focal_prefix}_")}

    expanded_in_synteny = len(set(expanded_ids) & focal_synteny_genes)
    contracted_in_synteny = len(set(contracted_ids) & focal_synteny_genes)

    n_expanded = len(expanded_ids)
    n_contracted = len(contracted_ids)

    print(
        f"  Total {focal_display} genes in synteny blocks: {len(focal_synteny_genes):,}"
    )
    if n_expanded > 0:
        print(
            f"  CAFE-expanded genes in synteny: "
            f"{expanded_in_synteny}/{n_expanded} "
            f"({expanded_in_synteny / n_expanded * 100:.1f}%)"
        )
    if n_contracted > 0:
        print(
            f"  CAFE-contracted genes in synteny: "
            f"{contracted_in_synteny}/{n_contracted} "
            f"({contracted_in_synteny / n_contracted * 100:.1f}%)"
        )

    # ------------------------------------------------------------------
    # Step 5b: Tandem duplication overlap
    # ------------------------------------------------------------------
    print("\n  --- Tandem Duplication Analysis ---")

    all_tandem_genes: Set[str] = set()
    for mcscanx_dir in mcscanx_dirs:
        comp_name = mcscanx_dir.name
        tandem_file = mcscanx_dir / f"{comp_name}.tandem"

        if not tandem_file.exists():
            tandem_files = list(mcscanx_dir.glob("*.tandem"))
            if tandem_files:
                tandem_file = tandem_files[0]

        if tandem_file.exists():
            genes = parse_tandem_file(tandem_file)
            all_tandem_genes.update(genes)
            print(f"  Loaded tandem genes from {tandem_file.name}")

    focal_tandem_genes = {
        g for g in all_tandem_genes if g.startswith(f"{focal_prefix}_")
    }
    print(f"  Total {focal_display} tandem genes: {len(focal_tandem_genes)}")

    expanded_tandem = 0
    contracted_tandem = 0
    if focal_tandem_genes:
        expanded_tandem = len(set(expanded_ids) & focal_tandem_genes)
        contracted_tandem = len(set(contracted_ids) & focal_tandem_genes)

        if n_expanded > 0:
            print(
                f"  CAFE-expanded genes in tandem duplications: "
                f"{expanded_tandem}/{n_expanded} "
                f"({expanded_tandem / n_expanded * 100:.1f}%)"
            )
        if n_contracted > 0:
            print(
                f"  CAFE-contracted genes in tandem duplications: "
                f"{contracted_tandem}/{n_contracted} "
                f"({contracted_tandem / n_contracted * 100:.1f}%)"
            )

        if expanded_tandem > contracted_tandem:
            print(
                f"\n  ★ Insight: CAFE-expanded families show more tandem duplications,"
            )
            print(
                f"    suggesting recent gene amplification events in "
                f"{focal_display} lineage."
            )

    # ------------------------------------------------------------------
    # Step 6: Generate visualizations
    # ------------------------------------------------------------------
    print("\n--- Step 6/6: Generate Visualizations ---\n")

    # Publication-quality matplotlib settings
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "font.family": "serif",
            "font.size": 14,
            "axes.titlesize": 16,
            "axes.titleweight": "bold",
            "axes.labelsize": 14,
            "axes.labelweight": "bold",
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "pdf.fonttype": 42,  # editable TrueType text in Illustrator (not Type 3)
            "ps.fonttype": 42,
        }
    )

    # --- Figure 1: CAFE-Synteny Overlap Bar Chart ---
    fig, ax = plt.subplots(figsize=(10, 6))

    categories = ["CAFE-Expanded", "CAFE-Contracted"]
    in_synteny = [expanded_in_synteny, contracted_in_synteny]
    not_in_synteny = [
        n_expanded - expanded_in_synteny,
        n_contracted - contracted_in_synteny,
    ]

    x = np.arange(len(categories))
    width = 0.35

    bars1 = ax.bar(
        x - width / 2,
        in_synteny,
        width,
        label="In Synteny Blocks",
        color="#2E8B57",
    )
    bars2 = ax.bar(
        x + width / 2,
        not_in_synteny,
        width,
        label="Not in Synteny",
        color="#DC143C",
        alpha=0.6,
    )

    ax.set_ylabel("Number of Genes", fontsize=14)
    ax.set_title(
        "CAFE-Significant Genes in Synteny Blocks", fontsize=16, fontweight="bold"
    )
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=12)
    ax.legend(fontsize=12)

    # Add percentage labels
    for bar, total in zip(bars1, [n_expanded, n_contracted]):
        if total > 0:
            height = bar.get_height()
            ax.annotate(
                f"{height / total * 100:.1f}%",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=11,
            )

    plt.tight_layout()

    synteny_png = figures_dir / "cafe_synteny_overlap.png"
    synteny_pdf = figures_dir / "cafe_synteny_overlap.pdf"
    fig.savefig(synteny_png, dpi=300, bbox_inches="tight")
    fig.savefig(synteny_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved {synteny_png.name} ({synteny_png.stat().st_size / 1024:.1f} KB)")
    print(f"  ✓ Saved {synteny_pdf.name}")

    # --- Figure 2: Tandem Duplication Bar Chart (if data available) ---
    if focal_tandem_genes and (n_expanded > 0 or n_contracted > 0):
        fig2, ax2 = plt.subplots(figsize=(10, 6))

        in_tandem = [expanded_tandem, contracted_tandem]
        not_in_tandem = [
            n_expanded - expanded_tandem,
            n_contracted - contracted_tandem,
        ]

        bars_t1 = ax2.bar(
            x - width / 2,
            in_tandem,
            width,
            label="Tandem Duplicated",
            color="#4169E1",
        )
        bars_t2 = ax2.bar(
            x + width / 2,
            not_in_tandem,
            width,
            label="Not Tandem",
            color="#FF8C00",
            alpha=0.6,
        )

        ax2.set_ylabel("Number of Genes", fontsize=14)
        ax2.set_title(
            "CAFE-Significant Genes in Tandem Duplications",
            fontsize=16,
            fontweight="bold",
        )
        ax2.set_xticks(x)
        ax2.set_xticklabels(categories, fontsize=12)
        ax2.legend(fontsize=12)

        for bar, total in zip(bars_t1, [n_expanded, n_contracted]):
            if total > 0:
                height = bar.get_height()
                ax2.annotate(
                    f"{height / total * 100:.1f}%",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=11,
                )

        plt.tight_layout()

        tandem_png = figures_dir / "cafe_tandem_overlap.png"
        tandem_pdf = figures_dir / "cafe_tandem_overlap.pdf"
        fig2.savefig(tandem_png, dpi=300, bbox_inches="tight")
        fig2.savefig(tandem_pdf, bbox_inches="tight")
        plt.close(fig2)
        print(
            f"  ✓ Saved {tandem_png.name} ({tandem_png.stat().st_size / 1024:.1f} KB)"
        )
        print(f"  ✓ Saved {tandem_pdf.name}")

    # --- Figure 3: Combined summary bar chart ---
    fig3, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: Synteny
    ax_syn = axes[0]
    bar_data_syn = {
        "Expanded": [expanded_in_synteny, n_expanded - expanded_in_synteny],
        "Contracted": [contracted_in_synteny, n_contracted - contracted_in_synteny],
    }
    bottom_syn = np.zeros(2)
    colors_syn = ["#2E8B57", "#DC143C"]
    labels_syn = ["In Synteny", "Not in Synteny"]

    for i, label in enumerate(labels_syn):
        vals = [bar_data_syn["Expanded"][i], bar_data_syn["Contracted"][i]]
        ax_syn.bar(
            categories,
            vals,
            bottom=bottom_syn,
            label=label,
            color=colors_syn[i],
            alpha=0.8 if i == 0 else 0.5,
        )
        bottom_syn += vals

    ax_syn.set_ylabel("Number of Genes")
    ax_syn.set_title("Synteny Block Overlap")
    ax_syn.legend(fontsize=10)

    # Right panel: Tandem
    ax_tan = axes[1]
    colors_tan = ["#4169E1", "#FF8C00"]
    labels_tan = ["Tandem Duplicated", "Not Tandem"]
    bottom_tan = np.zeros(2)

    tan_data = {
        "Expanded": [expanded_tandem, n_expanded - expanded_tandem],
        "Contracted": [contracted_tandem, n_contracted - contracted_tandem],
    }

    for i, label in enumerate(labels_tan):
        vals = [tan_data["Expanded"][i], tan_data["Contracted"][i]]
        ax_tan.bar(
            categories,
            vals,
            bottom=bottom_tan,
            label=label,
            color=colors_tan[i],
            alpha=0.8 if i == 0 else 0.5,
        )
        bottom_tan += vals

    ax_tan.set_ylabel("Number of Genes")
    ax_tan.set_title("Tandem Duplication Overlap")
    ax_tan.legend(fontsize=10)

    fig3.suptitle(
        f"CAFE-Synteny Integration: {focal_display}",
        fontsize=18,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()

    combined_png = figures_dir / "cafe_synteny_combined.png"
    combined_pdf = figures_dir / "cafe_synteny_combined.pdf"
    fig3.savefig(combined_png, dpi=300, bbox_inches="tight")
    fig3.savefig(combined_pdf, bbox_inches="tight")
    plt.close(fig3)
    print(
        f"  ✓ Saved {combined_png.name} ({combined_png.stat().st_size / 1024:.1f} KB)"
    )
    print(f"  ✓ Saved {combined_pdf.name}")

    # ------------------------------------------------------------------
    # Step 7: State-of-the-Art Visualizations
    # ------------------------------------------------------------------
    print("\n--- Step 7: State-of-the-Art Visualizations ---\n")

    stoa_dir = figures_dir / "StateOfTheArt"
    stoa_dir.mkdir(parents=True, exist_ok=True)

    for mcscanx_dir in mcscanx_dirs:
        comp_name = mcscanx_dir.name
        prefix1 = comp_name.split("_")[0]
        prefix2 = comp_name.split("_")[1] if "_" in comp_name else prefix1

        bed1_path = mcscanx_dir / f"{prefix1}-filtered-modified.bed"
        bed2_path = mcscanx_dir / f"{prefix2}-filtered-modified.bed"
        col_path = mcscanx_dir / f"{comp_name}.collinearity"

        if not col_path.exists():
            col_files = list(mcscanx_dir.glob("*.collinearity"))
            if col_files:
                col_path = col_files[0]

        print(f"  Generating plots for {comp_name}...")

        chr1_sizes = (
            load_chromosome_sizes_from_bed(bed1_path) if bed1_path.exists() else {}
        )
        chr2_sizes = (
            load_chromosome_sizes_from_bed(bed2_path) if bed2_path.exists() else {}
        )

        sp1_name = species_to_latex(None, prefix1.upper())
        sp2_name = species_to_latex(None, prefix2.upper())

        # --- Oxford / Dot Plot (matplotlib, no extra deps) ---
        if col_path.exists() and bed1_path.exists() and bed2_path.exists():
            pairs_df = parse_collinearity_for_dotplot(col_path, bed1_path, bed2_path)
            if not pairs_df.empty:
                inter_pairs = pairs_df[
                    pairs_df["chr1"].str.startswith(prefix1)
                    & pairs_df["chr2"].str.startswith(prefix2)
                ]
                if not inter_pairs.empty:
                    oxford_path = stoa_dir / f"oxford_plot_{comp_name}.png"
                    create_oxford_plot(
                        inter_pairs,
                        chr1_sizes,
                        chr2_sizes,
                        sp1_name,
                        sp2_name,
                        oxford_path,
                        chr1_prefix=prefix1,
                        chr2_prefix=prefix2,
                    )
                else:
                    print(f"    No inter-species pairs for Oxford plot ({comp_name})")
            else:
                print(f"    No syntenic pairs found for Oxford plot ({comp_name})")

        # --- Linear Ideogram (matplotlib, no extra deps) ---
        if chr1_sizes:
            density1 = calculate_gene_density(bed1_path) if bed1_path.exists() else {}
            if density1:
                ideo_path = stoa_dir / f"linear_ideogram_{prefix1}.png"
                create_linear_ideogram(chr1_sizes, density1, sp1_name, ideo_path)

        if chr2_sizes and prefix2 != prefix1:
            density2 = calculate_gene_density(bed2_path) if bed2_path.exists() else {}
            if density2:
                ideo_path = stoa_dir / f"linear_ideogram_{prefix2}.png"
                create_linear_ideogram(chr2_sizes, density2, sp2_name, ideo_path)

        # --- Genome Circos (requires pyCirclize) ---
        if PYCIRCLIZE_AVAILABLE:
            if chr1_sizes:
                density1 = (
                    calculate_gene_density(bed1_path) if bed1_path.exists() else {}
                )
                circos1_path = stoa_dir / f"genome_circos_{prefix1}.png"
                create_genome_circos(chr1_sizes, density1, sp1_name, circos1_path)

            if chr2_sizes and prefix2 != prefix1:
                density2 = (
                    calculate_gene_density(bed2_path) if bed2_path.exists() else {}
                )
                circos2_path = stoa_dir / f"genome_circos_{prefix2}.png"
                create_genome_circos(chr2_sizes, density2, sp2_name, circos2_path)

            # --- Dual-Genome Synteny Circos ---
            if (
                col_path.exists()
                and chr1_sizes
                and chr2_sizes
                and bed1_path.exists()
                and bed2_path.exists()
            ):
                syn_blocks = parse_synteny_blocks_for_circos(
                    col_path,
                    bed1_path,
                    bed2_path,
                    chr1_prefix=prefix1,
                    chr2_prefix=prefix2,
                    min_genes=10,
                )
                if syn_blocks:
                    syn_circos_path = stoa_dir / f"synteny_circos_{comp_name}.png"
                    create_synteny_circos(
                        chr1_sizes,
                        chr2_sizes,
                        syn_blocks,
                        focal_tag=prefix1.upper(),
                        comp_tag=prefix2.upper(),
                        focal_label=sp1_name,
                        comp_label=sp2_name,
                        output_path=syn_circos_path,
                    )
                else:
                    print(f"    No synteny blocks >= 10 genes for circos ({comp_name})")
        else:
            print("    pyCirclize not installed — skipping circos plots")
            print("    Install with: pip install pycirclize")

    # ------------------------------------------------------------------
    # Save summary table
    # ------------------------------------------------------------------
    summary_data = {
        "Category": ["CAFE-Expanded", "CAFE-Contracted"],
        "Total_Genes": [n_expanded, n_contracted],
        "In_Synteny": [expanded_in_synteny, contracted_in_synteny],
        "Not_In_Synteny": [
            n_expanded - expanded_in_synteny,
            n_contracted - contracted_in_synteny,
        ],
        "Synteny_Pct": [
            expanded_in_synteny / n_expanded * 100 if n_expanded > 0 else 0,
            contracted_in_synteny / n_contracted * 100 if n_contracted > 0 else 0,
        ],
        "Tandem_Duplicated": [expanded_tandem, contracted_tandem],
        "Tandem_Pct": [
            expanded_tandem / n_expanded * 100 if n_expanded > 0 else 0,
            contracted_tandem / n_contracted * 100 if n_contracted > 0 else 0,
        ],
    }
    summary_df = pd.DataFrame(summary_data)
    summary_tsv = output_dir / "cafe_synteny_integration_summary.tsv"
    summary_df.to_csv(summary_tsv, sep="\t", index=False)
    print(f"\n  ✓ Saved {summary_tsv.name}")

    # ------------------------------------------------------------------
    # Print final summary
    # ------------------------------------------------------------------
    print(f"\n{'=' * 60}")
    print(f"    CAFE-SYNTENY INTEGRATION SUMMARY")
    print(f"{'=' * 60}")

    n_unique_blocks = len(synteny_df["Block"].unique()) if not synteny_df.empty else 0
    print(f"\nSynteny Analysis:")
    print(f"  • Total synteny blocks analyzed: {n_unique_blocks:,}")
    print(f"  • {focal_display} genes in synteny: {len(focal_synteny_genes):,}")
    if n_expanded > 0:
        print(
            f"  • CAFE-expanded genes in synteny: "
            f"{expanded_in_synteny} ({expanded_in_synteny / n_expanded * 100:.1f}%)"
        )
    if n_contracted > 0:
        print(
            f"  • CAFE-contracted genes in synteny: "
            f"{contracted_in_synteny} ({contracted_in_synteny / n_contracted * 100:.1f}%)"
        )

    print(f"\nCAFE Families:")
    print(f"  • Significantly expanded families: {len(cafe_expanded)}")
    print(f"  • Significantly contracted families: {len(cafe_contracted)}")

    print(f"\nBiological Interpretation:")
    print(f"  • Genes in synteny blocks represent conserved genomic regions")
    print(f"  • CAFE-expanded genes NOT in synteny may represent novel duplications")
    print(f"  • CAFE-contracted genes in synteny may indicate ancestral gene loss")

    print(f"\nOutput Files:")
    print(f"  • {expanded_tsv}")
    print(f"  • {contracted_tsv}")
    print(f"  • {summary_tsv}")
    print(f"  • {synteny_png}")
    print(f"  • {synteny_pdf}")
    if focal_tandem_genes:
        tandem_png_path = figures_dir / "cafe_tandem_overlap.png"
        if tandem_png_path.exists():
            print(f"  • {tandem_png_path}")
    print(f"  • {combined_png}")
    print(f"  • {combined_pdf}")
    print(f"\n{'=' * 60}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CAFE-Synteny integration visualization for MCScanX results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Auto-detect from project config:
  %(prog)s

  # Explicit paths:
  %(prog)s \\
      --mcscanx-dirs ./Outputs_set3/MCScanX_Results/dm_dr \\
      --focal-species Diptychus_maculatus \\
      --cafe-family-results ./Outputs_set3/CAFE/cafe_results/single_lambda/Base_family_results.txt \\
      --cafe-change-tab ./Outputs_set3/CAFE/cafe_results/single_lambda/Base_change.tab \\
      --orthogroups-tsv ./Outputs/OrthoFinder/Results_Mar03/Orthogroups/Orthogroups.tsv \\
      --output-dir ./output/integration \\
      --figures-dir ./output/figures

  # Multiple MCScanX comparison directories:
  %(prog)s \\
      --mcscanx-dirs ./results/dm_dr ./results/dm_al ./results/dm_tp \\
      --focal-species Diptychus_maculatus \\
      ...
""",
    )

    parser.add_argument(
        "--mcscanx-dirs",
        type=Path,
        nargs="+",
        default=None,
        help=(
            "One or more MCScanX result directories containing "
            "{name}.collinearity and {name}.tandem files. "
            "Default: auto-detect from config."
        ),
    )
    parser.add_argument(
        "--focal-species",
        default=None,
        help=(
            "Focal species name (Genus_species format). "
            "Default: auto-detect from config."
        ),
    )
    parser.add_argument(
        "--cafe-family-results",
        type=Path,
        default=None,
        help="Path to CAFE Base_family_results.txt. Default: from config.",
    )
    parser.add_argument(
        "--cafe-change-tab",
        type=Path,
        default=None,
        help="Path to CAFE Base_change.tab. Default: from config.",
    )
    parser.add_argument(
        "--orthogroups-tsv",
        type=Path,
        default=None,
        help="Path to OrthoFinder Orthogroups.tsv. Default: from config.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for integration TSV outputs. Default: {mcscanx_results}/integration",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=None,
        help="Directory for figure outputs. Default: from config or ./Figures",
    )
    parser.add_argument(
        "--pvalue",
        type=float,
        default=0.05,
        help="CAFE significance p-value threshold (default: %(default)s)",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Resolve paths: try config system, then defaults
    # ------------------------------------------------------------------
    # Type-safe config function holders
    _cfg_get_path: Optional[object] = None
    _cfg_get_focal: Optional[object] = None
    _cfg_get_comps: Optional[object] = None

    try:
        sys.path.insert(0, "./Notebooks/config")
        from load_config import (
            get_path as _get_path,
            get_focal_species as _get_focal_species,
            get_mcscanx_comparisons as _get_mcscanx_comparisons,
        )

        _cfg_get_path = _get_path
        _cfg_get_focal = _get_focal_species
        _cfg_get_comps = _get_mcscanx_comparisons
        print("Config system loaded.")
    except (ImportError, Exception) as e:
        print(f"Config system not available ({e}); using explicit arguments.")

    # Focal species
    focal_species = args.focal_species
    if focal_species is None and _cfg_get_focal is not None:
        focal_species = _cfg_get_focal()  # type: ignore[operator]
    if focal_species is None:
        print(
            "Error: --focal-species is required (config not available).",
            file=sys.stderr,
        )
        sys.exit(1)

    # CAFE family results
    cafe_family_results = args.cafe_family_results
    if cafe_family_results is None and _cfg_get_path is not None:
        path_str = _cfg_get_path("cafe.single_lambda.family_results")  # type: ignore[operator]
        if path_str:
            cafe_family_results = Path(path_str)
    if cafe_family_results is None or not cafe_family_results.exists():
        print(
            f"Error: --cafe-family-results is required. Got: {cafe_family_results}",
            file=sys.stderr,
        )
        sys.exit(1)

    # CAFE change tab
    cafe_change_tab = args.cafe_change_tab
    if cafe_change_tab is None and _cfg_get_path is not None:
        path_str = _cfg_get_path("cafe.single_lambda.change_tab")  # type: ignore[operator]
        if path_str:
            cafe_change_tab = Path(path_str)
    if cafe_change_tab is None or not cafe_change_tab.exists():
        print(
            f"Error: --cafe-change-tab is required. Got: {cafe_change_tab}",
            file=sys.stderr,
        )
        sys.exit(1)

    # OrthoFinder orthogroups
    orthogroups_tsv = args.orthogroups_tsv
    if orthogroups_tsv is None and _cfg_get_path is not None:
        path_str = _cfg_get_path("orthofinder.orthogroups_tsv")  # type: ignore[operator]
        if path_str:
            orthogroups_tsv = Path(path_str)
    if orthogroups_tsv is None or not orthogroups_tsv.exists():
        print(
            f"Error: --orthogroups-tsv is required. Got: {orthogroups_tsv}",
            file=sys.stderr,
        )
        sys.exit(1)

    # MCScanX directories
    mcscanx_dirs = args.mcscanx_dirs
    if (
        mcscanx_dirs is None
        and _cfg_get_path is not None
        and _cfg_get_comps is not None
    ):
        mcscanx_results_str = _cfg_get_path("outputs.mcscanx_results")  # type: ignore[operator]
        mcscanx_db_str = _cfg_get_path("database.mcscanx_db")  # type: ignore[operator]
        comparisons = _cfg_get_comps()  # type: ignore[operator]

        # Try both possible locations for MCScanX results
        dirs: List[Path] = []
        for comp in comparisons:
            for base in [mcscanx_results_str, mcscanx_db_str]:
                if base:
                    candidate = Path(base) / comp["name"]
                    if candidate.exists():
                        dirs.append(candidate)
                        break

        if dirs:
            mcscanx_dirs = dirs
        else:
            # Fallback: search for any collinearity files under results dir
            if mcscanx_results_str:
                results_path = Path(mcscanx_results_str)
                if results_path.exists():
                    mcscanx_dirs = [
                        d
                        for d in results_path.iterdir()
                        if d.is_dir() and list(d.glob("*.collinearity"))
                    ]

    if not mcscanx_dirs:
        print(
            "Error: --mcscanx-dirs is required (no directories found).", file=sys.stderr
        )
        sys.exit(1)

    # Validate MCScanX directories
    valid_dirs: List[Path] = []
    for d in mcscanx_dirs:
        if d.exists() and d.is_dir():
            valid_dirs.append(d)
        else:
            print(f"  Warning: MCScanX dir not found: {d}", file=sys.stderr)
    mcscanx_dirs = valid_dirs

    if not mcscanx_dirs:
        print("Error: No valid MCScanX result directories found.", file=sys.stderr)
        sys.exit(1)

    # Output directory
    output_dir = args.output_dir
    if output_dir is None:
        # Place integration outputs alongside MCScanX results
        output_dir = mcscanx_dirs[0].parent / "integration"

    # Figures directory
    figures_dir = args.figures_dir
    if figures_dir is None and _cfg_get_path is not None:
        path_str = _cfg_get_path("outputs.figures")  # type: ignore[operator]
        if path_str:
            figures_dir = Path(path_str)
    if figures_dir is None:
        figures_dir = output_dir / "figures"

    # ------------------------------------------------------------------
    # Run the analysis
    # ------------------------------------------------------------------
    run_cafe_synteny_integration(
        mcscanx_dirs=mcscanx_dirs,
        focal_species=focal_species,
        cafe_family_results=cafe_family_results,
        cafe_change_tab=cafe_change_tab,
        orthogroups_tsv=orthogroups_tsv,
        output_dir=output_dir,
        figures_dir=figures_dir,
        pvalue_threshold=args.pvalue,
    )


if __name__ == "__main__":
    main()
