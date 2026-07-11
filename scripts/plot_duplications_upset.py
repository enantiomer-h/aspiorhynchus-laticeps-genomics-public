#!/usr/bin/env python3
"""Summary supplementary — multi-species gene-duplication UpSet (set3 p132).

PI request: "Genome-duplication-event analysis → add more species and build a
multi-species UpSet plot; must include zebrafish, A. laticeps, and D. maculatus."
No prior code existed; this adapts the orthogroup UpSet pattern
(scripts/plot_orthogroup_venn.py) to gene-duplication events.

Method: OrthoFinder's Gene_Duplication_Events/Duplications.tsv records every
inferred duplication and the species-tree node it maps to. Rows whose
"Species Tree Node" is an extant species are species-specific (terminal)
duplications. We build an orthogroup x species indicator of "had >=1 terminal
duplication" and UpSet the intersections across a roster that includes the three
required species plus the schizothoracine ingroup.

Data: Outputs/OrthoFinder/Results_Mar03/Gene_Duplication_Events/Duplications.tsv
Set-agnostic (same Mar03 run) → written to both Figures_set2/ and Figures_set3/.
"""
import sys
import warnings
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import (  # noqa: E402
    OF_RESULTS, SETS, italic_label, apply_pub_rcparams, save_fig,
)

DUPLICATIONS_TSV = OF_RESULTS / "Gene_Duplication_Events/Duplications.tsv"

# Roster for the UpSet — MUST include zebrafish, A. laticeps, D. maculatus (PI),
# plus the schizothoracine ingroup for context.
ROSTER = [
    "Danio_rerio",
    "Aspiorhynchus_laticeps",
    "Diptychus_maculatus",
    "Gymnocypris_eckloni",
    "Oxygymnocypris_stewartii",
    "Schizopygopsis_younghusbandi",
]
MIN_SUBSET = 500  # hide tiny intersections so the plot stays legible


def main() -> None:
    apply_pub_rcparams()
    roster_set = set(ROSTER)

    # Read only the two needed columns from the 89 MB file.
    dup = pd.read_csv(DUPLICATIONS_TSV, sep="\t",
                      usecols=["Orthogroup", "Species Tree Node"])
    term = dup[dup["Species Tree Node"].isin(roster_set)]
    print(f"Total duplication rows={len(dup):,}; terminal (roster species)={len(term):,}")

    # orthogroup x species boolean: did this species have >=1 terminal dup in this OG?
    pairs = term[["Orthogroup", "Species Tree Node"]].drop_duplicates()
    ind = (pairs.assign(v=True)
                .pivot_table(index="Orthogroup", columns="Species Tree Node",
                             values="v", fill_value=False))
    for sp in ROSTER:                       # ensure all roster cols present
        if sp not in ind.columns:
            ind[sp] = False
    ind = ind[ROSTER].astype(bool)
    print(f"Orthogroups with >=1 terminal duplication in any roster species: {len(ind):,}")
    per_species = ind.sum(axis=0)
    for sp in ROSTER:
        print(f"  {sp:32s} OGs with terminal dup = {int(per_species[sp]):,}")

    ind.columns = [italic_label(s, abbrev=True) for s in ROSTER]
    upset_input = from_indicators(ind.columns.tolist(), data=ind)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        fig = plt.figure(figsize=(16, 9))
        UpSet(upset_input, subset_size="count", sort_by="cardinality",
              show_counts=True, min_subset_size=MIN_SUBSET, element_size=None).plot(fig=fig)
        fig.suptitle("Gene-duplication sharing across species "
                     f"(orthogroups with species-specific duplications; subsets ≥{MIN_SUBSET})",
                     fontsize=15, fontweight="bold")
        # Enlarge the intersection-size count labels so they don't read as cramped.
        for ax_ in fig.get_axes():
            for txt in ax_.texts:
                txt.set_fontsize(11)

    for s, paths in SETS.items():
        written = save_fig(fig, paths["figures_dir"], "duplications_upset")
        print(f"[{s}] wrote {', '.join(p.name for p in written)}")
    plt.close(fig)


if __name__ == "__main__":
    main()
