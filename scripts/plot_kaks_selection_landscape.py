#!/usr/bin/env python3
"""Figure 6 (left panel) — genome-wide selection landscape (Ka/Ks = omega).

Replaces the previous `kaks_volcano_plot`, which was degenerate: it plotted the
codeml *branch-model* results, of which only ~9 orthogroups were analysable and
none were significant (the focal-foreground branch model found ~0 families with
codon alignments). That plot conveyed no selection signal.

Instead we summarise the GENOME-WIDE selection screen
(GenomewideSelection/summary/genome_wide_selection_master.tsv, thousands of
orthogroups with a per-OG omega) as a distribution over selection regimes,
highlighting the positive-selection tail (omega > 1). This is the real,
informative "distribution of selected genes" the PI's Fig 6 left panel calls for.

Output (per set): kaks_selection_landscape.{png,pdf,svg}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SETS, italic_label, apply_pub_rcparams, save_fig  # noqa: E402

# selection regimes by mean omega (matches the enrichment category cut-offs)
BANDS = [
    (0.0, 0.1, "#2166AC", "strong purifying\n(ω ≤ 0.1)"),
    (0.1, 0.5, "#67A9CF", "moderate purifying\n(0.1 < ω ≤ 0.5)"),
    (0.5, 1.0, "#FDB863", "relaxed\n(0.5 < ω ≤ 1)"),
    (1.0, 3.0, "#D6604D", "positive\n(ω > 1)"),
]
XCAP = 3.0


def main() -> None:
    apply_pub_rcparams()
    for s, paths in SETS.items():
        master = paths["outputs_dir"] / "GenomewideSelection/summary/genome_wide_selection_master.tsv"
        if not master.exists():
            print(f"[{s}] master not found: {master}; skipping")
            continue
        df = pd.read_csv(master, sep="\t", usecols=["Orthogroup", "mean_omega", "max_omega"])
        omega = pd.to_numeric(df["mean_omega"], errors="coerce").dropna()
        n_pos = int((omega > 1).sum())   # positive by MEAN omega (matches the shaded bands)
        focal = paths["focal"]

        fig, ax = plt.subplots(figsize=(12, 8))
        clipped = omega.clip(upper=XCAP)
        ax.hist(clipped, bins=60, color="#999999", edgecolor="white", linewidth=0.3, zorder=2)
        # shade selection regimes
        ymax = ax.get_ylim()[1]
        for lo, hi, color, label in BANDS:
            ax.axvspan(lo, min(hi, XCAP), color=color, alpha=0.18, zorder=1)
            n = int(((omega > lo) & (omega <= hi)).sum())
            xmid = (lo + min(hi, XCAP)) / 2
            ax.text(xmid, ymax * 0.96, f"n={n:,}", ha="center", va="top",
                    fontsize=12, fontweight="bold", color=color)
        ax.axvline(1.0, color="#B2182B", linestyle="--", linewidth=2, zorder=3)
        ax.set_xlim(0, XCAP)
        ax.set_xlabel(r"selection pressure  $\omega$ = Ka/Ks  (mean per orthogroup)")
        ax.set_ylabel("number of orthogroups")
        ax.set_title(f"$\\it{{{focal.split('_')[0][0]}.\\ {focal.split('_')[1]}}}$: "
                     f"genome-wide selection landscape  "
                     f"({len(omega):,} orthogroups; mean ω per orthogroup)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(handles=[Patch(facecolor=c, alpha=0.18, edgecolor="none",
                                 label=lab.replace("\n", " ")) for _, _, c, lab in BANDS],
                  loc="upper right", frameon=False)
        fig.tight_layout()
        written = save_fig(fig, paths["figures_dir"], "kaks_selection_landscape", close=True)
        print(f"[{s}] {focal}: {len(omega):,} OGs, {n_pos:,} positive (ω>1) -> {len(written)} files")


if __name__ == "__main__":
    main()
