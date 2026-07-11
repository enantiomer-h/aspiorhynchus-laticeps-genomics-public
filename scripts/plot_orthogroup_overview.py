#!/usr/bin/env python3
"""Figure 2 — orthogroup overview (PI revisions p19 + p20).

Two fixes requested by the PI:

  * p19 ("normalized result has errors → verify the normalization"):
    the previous Orthogroups_SpeciesOverlaps_Normalized applied a *global*
    MinMaxScaler across the whole pairwise matrix. That is misleading — the
    diagonal (each species' own orthogroup total) dominates the global min/max,
    so off-diagonal cells are squashed and no longer comparable. We replace it
    with the **Jaccard similarity** J(A,B)=|A∩B|/|A∪B| over the orthogroup
    presence/absence sets, which is the principled, symmetric 0–1 metric for
    "how similar are two species' orthogroup repertoires".

  * p20 ("too much information, magnitudes differ greatly → extract the key
    signal: number of species-specific orthogroups"): we add a focused
    horizontal bar of the number of species-specific orthogroups per species
    (orthogroups present in exactly one species).

Outputs (set-agnostic — both sets share the Mar03 OrthoFinder run; written to
both Figures_set2/ and Figures_set3/):
  Orthogroups_SpeciesOverlaps.{png,pdf,svg}             raw shared-OG counts
  Orthogroups_SpeciesOverlaps_Normalized.{png,pdf,svg}  Jaccard similarity
  orthogroup_species_specific_counts.{png,pdf,svg}      key signal (p20)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import (  # noqa: E402
    OF_RESULTS, SETS, italic_label, apply_pub_rcparams, save_fig_both_sets,
)

GENECOUNT_TSV = OF_RESULTS / "Orthogroups/Orthogroups.GeneCount.tsv"


def load_presence() -> pd.DataFrame:
    gc = pd.read_csv(GENECOUNT_TSV, sep="\t", index_col=0)
    if "Total" in gc.columns:
        gc = gc.drop(columns=["Total"])
    return gc > 0  # orthogroup x species boolean presence matrix


def overlap_counts(presence: pd.DataFrame) -> pd.DataFrame:
    """Symmetric matrix: cell[i,j] = # orthogroups containing both species i,j."""
    mat = presence.astype(int)
    counts = mat.T.dot(mat)
    return counts


def jaccard_matrix(presence: pd.DataFrame) -> pd.DataFrame:
    """Symmetric Jaccard similarity over orthogroup presence sets."""
    mat = presence.astype(int)
    inter = mat.T.dot(mat)                       # |A ∩ B|
    sizes = np.diag(inter.values)                # |A|
    union = sizes[:, None] + sizes[None, :] - inter.values  # |A|+|B|-|A∩B|
    with np.errstate(divide="ignore", invalid="ignore"):
        jac = np.where(union > 0, inter.values / union, 0.0)
    return pd.DataFrame(jac, index=inter.index, columns=inter.columns)


def heatmap(matrix: pd.DataFrame, title: str, basename: str,
            fmt: str, cbar_label: str, cmap: str) -> None:
    labels = [italic_label(s) for s in matrix.index]
    fig, ax = plt.subplots(figsize=(13, 11))
    sns.heatmap(
        matrix, ax=ax, cmap=cmap, square=True, linewidths=0.5, linecolor="white",
        annot=True, fmt=fmt, annot_kws={"size": 10},
        xticklabels=labels, yticklabels=labels,
        cbar_kws={"label": cbar_label, "shrink": 0.8},
    )
    ax.set_xticklabels(labels, rotation=45, ha="right", fontstyle="italic")
    ax.set_yticklabels(labels, rotation=0, fontstyle="italic")
    ax.set_title(title)
    fig.tight_layout()
    written = save_fig_both_sets(fig, basename, close=True)
    print(f"wrote {basename}.{{png,pdf,svg}} -> {len(written)} files")


def species_specific_bar(presence: pd.DataFrame) -> None:
    """Key signal (p20): # orthogroups present in exactly one species."""
    species_specific_mask = presence.sum(axis=1) == 1
    counts = presence.loc[species_specific_mask].sum(axis=0).sort_values()

    fig, ax = plt.subplots(figsize=(12, 8))
    focal = {SETS["set2"]["focal"], SETS["set3"]["focal"]}
    colors = ["#d7191c" if s in focal else "#4575b4" for s in counts.index]
    bars = ax.barh(range(len(counts)), counts.values, color=colors,
                   edgecolor="black", linewidth=0.4)
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels([italic_label(s) for s in counts.index])
    ax.set_xlabel("number of species-specific orthogroups")
    ax.spines[["top", "right"]].set_visible(False)
    xmax = counts.max()
    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_width() + xmax * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{int(val):,}", va="center", ha="left", fontsize=11)
    ax.set_xlim(0, xmax * 1.12)
    # Legend explaining the two focal-species highlights.
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#d7191c", edgecolor="black",
              label="focal species (" + italic_label("Aspiorhynchus_laticeps", abbrev=True)
                    + ", " + italic_label("Diptychus_maculatus", abbrev=True) + ")"),
        Patch(facecolor="#4575b4", edgecolor="black", label="other species"),
    ], loc="lower right", frameon=False)
    fig.tight_layout()
    written = save_fig_both_sets(fig, "orthogroup_species_specific_counts", close=True)
    print(f"wrote orthogroup_species_specific_counts.{{png,pdf,svg}} -> {len(written)} files")


def main() -> None:
    apply_pub_rcparams()
    presence = load_presence()
    print(f"Loaded presence matrix: {presence.shape[0]} orthogroups x {presence.shape[1]} species")

    heatmap(overlap_counts(presence),
            "Shared orthogroups between species (counts)",
            "Orthogroups_SpeciesOverlaps", fmt="d",
            cbar_label="shared orthogroups", cmap="YlGnBu")

    heatmap(jaccard_matrix(presence),
            "Orthogroup repertoire similarity (Jaccard)",
            "Orthogroups_SpeciesOverlaps_Normalized", fmt=".2f",
            cbar_label="Jaccard similarity", cmap="viridis")

    species_specific_bar(presence)


if __name__ == "__main__":
    main()
