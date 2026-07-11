#!/usr/bin/env python3
"""Build per-species stacked bar of gene-family categories from OrthoFinder outputs.

Categories (matching reference paper Fig 3A, Zhao et al. Sci Data 2026):
  Single-copy orthologs / Multi-copy orthologs / Other orthologs / Unique paralogs / Unclustered

Outputs are written to BOTH set2 and set3 directories (set2 and set3 share
the same 12-species roster from the Mar03 OrthoFinder run).
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import OF_RESULTS, SETS, italic_label, apply_pub_rcparams, save_fig  # noqa: E402

GENECOUNT_TSV = OF_RESULTS / "Orthogroups/Orthogroups.GeneCount.tsv"
SINGLECOPY_TXT = OF_RESULTS / "Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
PERSPECIES_TSV = OF_RESULTS / "Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv"


def main() -> None:
    gc = pd.read_csv(GENECOUNT_TSV, sep="\t", index_col=0)
    if "Total" in gc.columns:
        gc = gc.drop(columns=["Total"])
    species = list(gc.columns)

    single_copy_set = set(SINGLECOPY_TXT.read_text().split())
    species_specific_mask = (gc > 0).sum(axis=1) == 1
    single_copy_mask = gc.index.isin(single_copy_set)
    shared_mask = (~species_specific_mask) & (~single_copy_mask)

    perspecies_full = pd.read_csv(PERSPECIES_TSV, sep="\t", header=None, dtype=str)
    header_row = perspecies_full.iloc[0, 1:].tolist()
    rows_by_label = {row[0]: row[1:] for row in perspecies_full.itertuples(index=False)}

    def num_row(label: str) -> pd.Series:
        return pd.to_numeric(pd.Series(rows_by_label[label], index=header_row)).reindex(species).fillna(0)

    total_genes = num_row("Number of genes")
    unassigned = num_row("Number of unassigned genes")

    sc_genes = gc.loc[single_copy_mask, species].sum(axis=0)
    sp_specific_genes = gc.loc[species_specific_mask, species].sum(axis=0)
    shared_subset = gc.loc[shared_mask, species]
    multi_copy_genes = shared_subset.where(shared_subset >= 2, 0).sum(axis=0)
    other_ortholog_genes = shared_subset.where(shared_subset == 1, 0).sum(axis=0)

    cats = pd.DataFrame({
        "Single-copy orthologs": sc_genes,
        "Multi-copy orthologs": multi_copy_genes,
        "Other orthologs": other_ortholog_genes,
        "Unique paralogs": sp_specific_genes,
        "Unclustered genes": unassigned,
    })

    delta = total_genes - cats.sum(axis=1)
    print("Per-species accounting (total - sum_of_categories, expect ~0):")
    print(delta.to_string())

    # Publication style: large fonts + editable PDF/SVG text (shared helper).
    apply_pub_rcparams()

    colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    fig, ax = plt.subplots(figsize=(12, 8))
    cats.plot(kind="barh", stacked=True, ax=ax, color=colors, edgecolor="black", linewidth=0.4)
    # PI (results-summary p3): x-axis label = "numbers"; no top title.
    ax.set_xlabel("numbers")
    ax.set_ylabel("")
    ax.set_yticklabels([italic_label(s) for s in cats.index])
    ax.invert_yaxis()
    ax.legend(title="Gene family category", loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    cats_with_total = cats.assign(Total=total_genes)
    for s, paths in SETS.items():
        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        written = save_fig(fig, paths["figures_dir"], "gene_family_categories_bar")
        csv = paths["tables_dir"] / "gene_family_categories.csv"
        cats_with_total.to_csv(csv)
        print(f"[{s}] wrote {', '.join(p.name for p in written)}, {csv.name}")
    plt.close(fig)


if __name__ == "__main__":
    main()
