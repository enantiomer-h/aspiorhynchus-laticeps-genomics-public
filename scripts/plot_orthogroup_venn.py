#!/usr/bin/env python3
"""Venn + UpSet plot of orthogroup sharing across schizothoracines.

Produces TWO Venn diagrams (one per set, focused on the set's focal species)
plus a shared UpSet plot covering all 5 schizothoracines.

  - set2: focal = A. laticeps; Venn = AL + Diptychus + Gymnocypris
  - set3: focal = D. maculatus; Venn = DM + AL + Gymnocypris
  - shared UpSet = AL + DM + Gymnocypris + Oxygymnocypris + Schizopygopsis
"""
import sys
import warnings
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from upsetplot import UpSet, from_indicators

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import OF_RESULTS, SETS, italic_label, apply_pub_rcparams, save_fig  # noqa: E402

GENECOUNT_TSV = OF_RESULTS / "Orthogroups/Orthogroups.GeneCount.tsv"

VENN_BY_SET = {
    "set2": ["Aspiorhynchus_laticeps", "Diptychus_maculatus", "Gymnocypris_eckloni"],
    "set3": ["Diptychus_maculatus", "Aspiorhynchus_laticeps", "Gymnocypris_eckloni"],
}
UPSET_SPECIES = [
    "Aspiorhynchus_laticeps",
    "Diptychus_maculatus",
    "Gymnocypris_eckloni",
    "Oxygymnocypris_stewartii",
    "Schizopygopsis_younghusbandi",
]


def make_venn(presence: pd.DataFrame, species_list: list[str], focal: str,
              title: str, figures_dir: Path, basename: str) -> None:
    sets = [set(presence.index[presence[sp]]) for sp in species_list]
    # Larger canvas + larger fonts so the PI's "overlapping numbers" complaint
    # is resolved (small subset regions get more room).
    fig, ax = plt.subplots(figsize=(12, 11))
    v = venn3(sets, set_labels=[italic_label(s) for s in species_list],
              set_colors=("#e41a1c", "#377eb8", "#4daf4a"), alpha=0.55, ax=ax)
    venn3_circles(sets, linewidth=1.4, ax=ax)
    for txt in v.set_labels or []:
        if txt is not None:
            txt.set_fontsize(20)
            txt.set_fontweight("bold")
    for txt in v.subset_labels or []:
        if txt is not None:
            txt.set_fontsize(16)
    # Nudge the three tiny pairwise-only labels outward so they stop colliding
    # with the central triple-overlap count.
    for sid, dy in (("100", 0.0), ("010", 0.0), ("001", 0.0),
                    ("110", 0.06), ("101", -0.06), ("011", 0.06)):
        lbl = v.get_label_by_id(sid)
        if lbl is not None:
            x, y = lbl.get_position()
            lbl.set_position((x, y + dy))
    ax.set_title(title)
    fig.tight_layout()
    save_fig(fig, figures_dir, basename, close=True)


def main() -> None:
    apply_pub_rcparams()

    gc = pd.read_csv(GENECOUNT_TSV, sep="\t", index_col=0)
    if "Total" in gc.columns:
        gc = gc.drop(columns=["Total"])
    presence = gc > 0

    for s, paths in SETS.items():
        paths["figures_dir"].mkdir(parents=True, exist_ok=True)
        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        species_list = VENN_BY_SET[s]
        focal = paths["focal"]
        focal_short = italic_label(focal, abbrev=True).replace("$", "")
        title = f"Shared orthogroups: {italic_label(focal)} and 2 related schizothoracines"
        make_venn(presence, species_list, focal, title,
                  paths["figures_dir"], f"orthogroup_venn3_{s}")
        print(f"[{s}] wrote orthogroup_venn3_{s}.{{png,pdf,svg}} (focal {focal})")

        sets_data = {sp: set(presence.index[presence[sp]]) for sp in species_list}
        rows = []
        for sp in species_list:
            others = [o for o in species_list if o != sp]
            only = sets_data[sp] - set().union(*(sets_data[o] for o in others))
            rows.append({"Species": sp, "Total OGs containing": len(sets_data[sp]),
                         f"OGs in {sp} only (vs other 2)": len(only)})
        all_three = sets_data[species_list[0]] & sets_data[species_list[1]] & sets_data[species_list[2]]
        rows.append({"Species": "ALL THREE shared", "Total OGs containing": len(all_three), f"OGs in {sp} only (vs other 2)": ""})
        pd.DataFrame(rows).to_csv(paths["tables_dir"] / f"orthogroup_venn3_summary_{s}.csv", index=False)

    bool_df = presence[UPSET_SPECIES].copy()
    bool_df.columns = [italic_label(s, abbrev=True) for s in UPSET_SPECIES]
    upset_input = from_indicators(bool_df.columns.tolist(), data=bool_df)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        fig2 = plt.figure(figsize=(16, 9))
        UpSet(upset_input, subset_size="count", sort_by="cardinality",
              show_counts=True, min_subset_size=200, element_size=None).plot(fig=fig2)
        fig2.suptitle("Orthogroup sharing among the 5 schizothoracines (subsets ≥200)",
                      fontsize=17, fontweight="bold")
    for s, paths in SETS.items():
        save_fig(fig2, paths["figures_dir"], "orthogroup_upset_schizothoracines")
        print(f"[{s}] wrote orthogroup_upset_schizothoracines.{{png,pdf,svg}}")
    plt.close(fig2)

    rows = []
    for sp in UPSET_SPECIES:
        rows.append({"Species": sp, "OGs containing species": int(presence[sp].sum())})
    all_5 = presence[UPSET_SPECIES].all(axis=1).sum()
    for s, paths in SETS.items():
        focal = paths["focal"]
        if focal in UPSET_SPECIES:
            others = [o for o in UPSET_SPECIES if o != focal]
            focal_only = (presence[focal] & ~presence[others].any(axis=1)).sum()
        else:
            focal_only = "n/a"
        keys = pd.DataFrame([
            {"Subset": "OGs in ALL 5 schizothoracines", "Count": int(all_5)},
            {"Subset": f"OGs in focal {focal} ONLY (vs other 4)", "Count": focal_only},
        ])
        keys.to_csv(paths["tables_dir"] / f"orthogroup_upset_key_intersections_{s}.csv", index=False)
        print(f"[{s}] all-5-shared={all_5}, focal-{focal}-only={focal_only}")


if __name__ == "__main__":
    main()
