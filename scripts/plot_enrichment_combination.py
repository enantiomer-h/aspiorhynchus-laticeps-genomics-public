#!/usr/bin/env python3
"""Figure 7 — CAFE × KaKs joint-enrichment combination barcharts.

Faithful extraction of the notebook-8 `create_combination_barchart` (so the
PI-reviewed content + drafted caption are preserved) routed through the shared
editable/large-font helpers (adds SVG + editable PDF text — the only thing the
PI's "confirm final + fonts" request needs).

Reads the shared combination data in `Outputs/EnrichmentComparison/` (produced by
`scripts/find_common_enrichment_terms.py`) and writes the four-colour barcharts to
both `Figures_set{2,3}/EnrichmentComparison/`.
"""
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import PROJECT_ROOT, SETS, apply_pub_rcparams, save_fig  # noqa: E402

ENRICH_DIR = PROJECT_ROOT / "Outputs/EnrichmentComparison"
COMBINATION_COLORS = {
    "Expanded + Positive": "#E41A1C",     # adaptive expansion
    "Expanded + Purifying": "#377EB8",    # essential expansion
    "Contracted + Positive": "#FF7F00",   # adaptive contraction
    "Contracted + Purifying": "#4DAF4A",  # essential contraction
}
ONTOLOGIES = {
    "KO": "KEGG Pathways: CAFE × KaKs Enrichment",
    "GO_BP": "GO Biological Process: CAFE × KaKs Enrichment",
    "GO_MF": "GO Molecular Function: CAFE × KaKs Enrichment",
    "GO_CC": "GO Cellular Component: CAFE × KaKs Enrichment",
}
MAX_TERMS = 15


def get_combination_label(analyses_str):
    if pd.isna(analyses_str):
        return None
    a = [x.strip().lower() for x in str(analyses_str).split(";")]
    exp, con = any("expanded" in x for x in a), any("contracted" in x for x in a)
    pos, pur = any("positive" in x for x in a), any("purifying" in x for x in a)
    if exp and pos:
        return "Expanded + Positive"
    if exp and pur:
        return "Expanded + Purifying"
    if con and pos:
        return "Contracted + Positive"
    if con and pur:
        return "Contracted + Purifying"
    return None


def preprocess(summary_df):
    if summary_df.empty or "analyses" not in summary_df.columns:
        return pd.DataFrame()
    df = summary_df.copy()
    df["combination"] = df["analyses"].apply(get_combination_label)
    df = df[df["combination"].notna()].copy()
    df["neg_log10_pval"] = -np.log10(df["min_p_adjust"].clip(lower=1e-50))
    return df


def combination_barchart(df, title):
    if df.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, f"No CAFE × KaKs combinations for {title}", ha="center", va="center")
        ax.set_axis_off()
        return fig
    d = df.nsmallest(MAX_TERMS, "min_p_adjust").copy()
    d["wrapped"] = d["Description"].apply(lambda x: textwrap.fill(str(x), width=45))
    fig, ax = plt.subplots(figsize=(12, 10))
    for y, (_, row) in zip(range(len(d)), d.iterrows()):
        ax.barh(y, row["neg_log10_pval"], height=0.7,
                color=COMBINATION_COLORS.get(row["combination"], "gray"),
                alpha=0.85, edgecolor="black", linewidth=0.5)
        ax.text(row["neg_log10_pval"] + 0.15, y, f"n={int(row['total_count'])}",
                va="center", ha="left", fontsize=13, fontweight="bold")
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d["wrapped"], fontsize=13, linespacing=0.95)
    ax.set_xlabel("$-\\log_{10}$(adjusted p-value)", fontsize=14, fontweight="bold")
    ax.set_title(title, fontsize=18, fontweight="bold", pad=15)
    ax.invert_yaxis()
    sig = -np.log10(0.05)
    ax.axvline(sig, color="red", linestyle="--", linewidth=1.5, alpha=0.7)
    ax.text(sig + 0.1, -0.5, "p = 0.05", fontsize=12, color="red", va="bottom")
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, linestyle="--", alpha=0.3)
    ax.set_xlim(0, d["neg_log10_pval"].max() * 1.15)
    ax.legend(handles=[mpatches.Patch(facecolor=c, edgecolor="black", linewidth=0.5, label=l)
                       for l, c in COMBINATION_COLORS.items()],
              title="CAFE × KaKs Combination", loc="lower right", frameon=True,
              fontsize=12, title_fontsize=13)
    fig.tight_layout()
    return fig


def main() -> None:
    apply_pub_rcparams()
    for ont, title in ONTOLOGIES.items():
        f = ENRICH_DIR / f"{ont}_common_terms_summary.csv"
        if not f.exists():
            print(f"skip {ont}: {f} missing")
            continue
        proc = preprocess(pd.read_csv(f))
        fig = combination_barchart(proc, title)
        for s, paths in SETS.items():
            outdir = paths["figures_dir"] / "EnrichmentComparison"
            save_fig(fig, outdir, f"{ont}_combination_barchart")
        plt.close(fig)
        ncombo = 0 if proc.empty else len(proc)
        print(f"{ont}: {ncombo} CAFE×KaKs terms -> {ont}_combination_barchart.{{png,pdf,svg}} (both sets)")


if __name__ == "__main__":
    main()
