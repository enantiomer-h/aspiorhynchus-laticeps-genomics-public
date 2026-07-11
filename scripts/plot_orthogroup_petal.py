#!/usr/bin/env python3
"""Summary supplementary — orthogroup petal / flower diagram (results-summary p4).

PI request: "Build a petal/flower diagram: center = shared gene families, petals
= species-specific gene families; mind the layout and enlarge fonts."  No prior
code existed; this follows the standard pan-genome "flower plot".

  * center hub  = number of CORE orthogroups (present in all 12 species)
  * each petal  = one species, labelled with its species-specific orthogroup
                  count (orthogroups present in that species only)

Data: Outputs/OrthoFinder/Results_Mar03/Orthogroups.GeneCount.tsv
Set-agnostic (same Mar03 run) → written to both Figures_set2/ and Figures_set3/.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import (  # noqa: E402
    OF_RESULTS, SETS, SPECIES, italic_label, apply_pub_rcparams, save_fig_both_sets,
)

GENECOUNT_TSV = OF_RESULTS / "Orthogroups/Orthogroups.GeneCount.tsv"


def main() -> None:
    apply_pub_rcparams()
    gc = pd.read_csv(GENECOUNT_TSV, sep="\t", index_col=0)
    if "Total" in gc.columns:
        gc = gc.drop(columns=["Total"])
    presence = gc[SPECIES] > 0

    core = int((presence.sum(axis=1) == len(SPECIES)).sum())
    specific = {sp: int(((presence.sum(axis=1) == 1) & presence[sp]).sum()) for sp in SPECIES}

    n = len(SPECIES)
    focal = {SETS["set2"]["focal"], SETS["set3"]["focal"]}
    cmap = plt.cm.tab20(np.linspace(0, 1, n))

    fig, ax = plt.subplots(figsize=(13, 13))
    petal_r = 1.55       # distance of petal centre from origin
    petal_len = 2.7      # petal major axis (radial)
    petal_wid = 0.95     # petal minor axis (tangential)

    for i, sp in enumerate(SPECIES):
        theta = np.pi / 2 - 2 * np.pi * i / n          # start at top, go clockwise
        cx, cy = petal_r * np.cos(theta), petal_r * np.sin(theta)
        ell = Ellipse((cx, cy), width=petal_len, height=petal_wid,
                      angle=np.degrees(theta), facecolor=cmap[i],
                      edgecolor=("#1A237E" if sp in focal else "black"),
                      lw=(3.0 if sp in focal else 1.0), alpha=0.78, zorder=2)
        ax.add_patch(ell)

        # species-specific count near the outer half of the petal
        lr = petal_r + petal_len * 0.28
        ax.text(lr * np.cos(theta), lr * np.sin(theta), f"{specific[sp]:,}",
                ha="center", va="center", fontsize=13, fontweight="bold", zorder=3)

        # species name outside the petal, rotated to read outward
        nr = petal_r + petal_len * 0.62
        deg = np.degrees(theta)
        rot = deg if -90 <= deg <= 90 else deg + 180
        ha = "left" if np.cos(theta) >= 0 else "right"
        ax.text(nr * np.cos(theta), nr * np.sin(theta), italic_label(sp),
                ha=ha, va="center", fontsize=14, rotation=rot,
                rotation_mode="anchor", zorder=3)

    # central hub: core orthogroups
    ax.add_patch(Circle((0, 0), 0.95, facecolor="white", edgecolor="black", lw=1.5, zorder=4))
    ax.text(0, 0.16, f"{core:,}", ha="center", va="center", fontsize=20,
            fontweight="bold", zorder=5)
    ax.text(0, -0.30, "core\northogroups", ha="center", va="center", fontsize=12, zorder=5)

    ax.set_xlim(-4.7, 4.7)
    ax.set_ylim(-4.7, 4.7)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Orthogroup distribution: core families (center) and "
                 "species-specific families (petals)", fontsize=16, fontweight="bold")
    fig.tight_layout()
    written = save_fig_both_sets(fig, "orthogroup_petal", close=True)
    print(f"core orthogroups = {core:,}")
    for sp in SPECIES:
        print(f"  {sp:32s} species-specific = {specific[sp]:,}")
    print(f"wrote orthogroup_petal.{{png,pdf,svg}} -> {len(written)} files")


if __name__ == "__main__":
    main()
