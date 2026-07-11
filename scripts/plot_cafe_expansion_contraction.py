#!/usr/bin/env python3
"""Figure 3 — CAFE gene-family expansion / contraction (PI revisions).

PI feedback (set3 p43 + set2 p38): "merge the two figures into one; highlight
A. laticeps and D. maculatus; put contraction counts in red on the left and
expansion counts after the bars on the right; enlarge text."  And (p45/p40):
"the two small panels are duplicates → keep only one; enlarge fonts and tidy."

Because set2 and set3 share the same CAFE run (identical species tree), the
per-species expansion/contraction counts are identical across sets — so a
single chart that highlights BOTH focal species is the correct "merge".

Data: Outputs_set{2,3}/CAFE/cafe_results/single_lambda/Base_clade_results.txt
      (columns: #Taxon_ID  Increase  Decrease, one row per tree node).

Outputs (per set):
  cafe_expansions_contractions.{png,pdf,svg}        diverging butterfly chart
  cafe_species_expansion_contraction.{png,pdf,svg}  (same canonical chart)
  cafe_species_changes_barplot.{png,pdf,svg}        single tidy NET-change panel
"""
import re
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import (  # noqa: E402
    SETS, SPECIES, italic_label, apply_pub_rcparams, save_fig,
)

CLADE_REL = "CAFE/cafe_results/single_lambda/Base_clade_results.txt"
GREEN, RED, GREY = "#27AE60", "#E74C3C", "#95A5A6"
SPECIES_SET = set(SPECIES)


def load_clade_results(outputs_dir: Path) -> pd.DataFrame:
    """Return per-extant-species DataFrame: Species, Expansions, Contractions."""
    path = outputs_dir / CLADE_REL
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    # Strip the CAFE node-id suffix "<N>"; keep only extant (named) species.
    df["Species"] = df["Taxon_ID"].apply(lambda t: re.sub(r"<\d+>$", "", t))
    df = df[df["Species"].isin(SPECIES_SET)].copy()
    df = df.rename(columns={"Increase": "Expansions", "Decrease": "Contractions"})
    df["Net"] = df["Expansions"] - df["Contractions"]
    return df[["Species", "Expansions", "Contractions", "Net"]].reset_index(drop=True)


def _ylabel(species: str, focal: set[str]) -> str:
    lab = italic_label(species)
    return r"$\bf{\rightarrow}$ " + lab if species in focal else lab


def diverging_chart(df: pd.DataFrame, focal: set[str], figures_dir: Path,
                    basename: str) -> None:
    """Butterfly chart: contractions (red, left) vs expansions (green, right)."""
    d = df.sort_values("Expansions", ascending=True).reset_index(drop=True)
    y = range(len(d))
    fig, ax = plt.subplots(figsize=(13, 9))

    ax.barh(y, d["Expansions"], color=GREEN, edgecolor="black", linewidth=0.5,
            label="Expansions", zorder=3)
    ax.barh(y, -d["Contractions"], color=RED, edgecolor="black", linewidth=0.5,
            label="Contractions", zorder=3)
    ax.axvline(0, color="black", linewidth=1.0, zorder=4)

    xmax_r = d["Expansions"].max()
    xmax_l = d["Contractions"].max()
    pad = 0.02 * (xmax_r + xmax_l)
    for i, row in d.iterrows():
        # Expansion count to the RIGHT of its bar.
        ax.text(row["Expansions"] + pad, i, f"{int(row['Expansions']):,}",
                va="center", ha="left", color=GREEN, fontsize=12, fontweight="bold")
        # Contraction count in RED to the LEFT of its bar.
        ax.text(-row["Contractions"] - pad, i, f"{int(row['Contractions']):,}",
                va="center", ha="right", color=RED, fontsize=12, fontweight="bold")

    ax.set_yticks(list(y))
    ax.set_yticklabels([_ylabel(s, focal) for s in d["Species"]])
    # Highlight the two focal species' tick labels.
    for tick, s in zip(ax.get_yticklabels(), d["Species"]):
        if s in focal:
            tick.set_fontweight("bold")
            tick.set_color("#222222")
    ax.set_xlim(-xmax_l * 1.30, xmax_r * 1.30)
    ax.set_xlabel("← contracted families        gene families        expanded families →")
    # Show |value| on the x-axis (bars are signed only for the diverging layout).
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels([f"{abs(int(t)):,}" for t in ax.get_xticks()])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(loc="lower right", frameon=True, framealpha=0.9)
    fig.tight_layout()
    save_fig(fig, figures_dir, basename, close=True)


def net_change_panel(df: pd.DataFrame, focal: set[str], figures_dir: Path,
                     basename: str) -> None:
    """Single tidy panel: NET change (expansions - contractions) per species."""
    d = df.sort_values("Net", ascending=True).reset_index(drop=True)
    y = range(len(d))
    colors = [RED if n < 0 else GREEN for n in d["Net"]]
    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.barh(y, d["Net"], color=colors, edgecolor="black", linewidth=0.5, zorder=3)
    # Outline the two focal species so they stand out.
    for bar, s in zip(bars, d["Species"]):
        if s in focal:
            bar.set_edgecolor("#1A237E")
            bar.set_linewidth(2.2)
    ax.axvline(0, color="black", linewidth=1.0)
    ax.set_yticks(list(y))
    ax.set_yticklabels([_ylabel(s, focal) for s in d["Species"]])
    span = max(abs(d["Net"].min()), abs(d["Net"].max()))
    for i, n in enumerate(d["Net"]):
        ax.text(n + (0.015 * span if n >= 0 else -0.015 * span), i,
                f"{int(n):+,}", va="center", ha="left" if n >= 0 else "right",
                fontsize=12, fontweight="bold")
    ax.set_xlim(-span * 1.25, span * 1.25)
    ax.set_xlabel("net change in gene families (expansions − contractions)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_fig(fig, figures_dir, basename, close=True)


def main() -> None:
    apply_pub_rcparams()
    focal = {SETS["set2"]["focal"], SETS["set3"]["focal"]}
    for s, paths in SETS.items():
        df = load_clade_results(paths["outputs_dir"])
        figdir = paths["figures_dir"]
        diverging_chart(df, focal, figdir, "cafe_expansions_contractions")
        diverging_chart(df, focal, figdir, "cafe_species_expansion_contraction")
        net_change_panel(df, focal, figdir, "cafe_species_changes_barplot")
        # tidy table alongside
        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        df.to_csv(paths["tables_dir"] / "cafe_species_expansion_contraction.csv", index=False)
        focal_rows = df[df["Species"].isin(focal)]
        print(f"[{s}] {len(df)} species; focal:\n{focal_rows.to_string(index=False)}")


if __name__ == "__main__":
    main()
