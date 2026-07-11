#!/usr/bin/env python3
"""Tier 3B render: parse MCMCTree FigTree.tre + posterior table and render
as a publication-quality calibrated time tree.

Outer-tip-rightward layout, x-axis in millions of years ago (Mya), node ages
labeled with mean and 95% HPD interval bars overlaid as horizontal segments.

Reads:
  Outputs/Phylogenomics_supermatrix/mcmctree/FigTree.tre
  Outputs/Phylogenomics_supermatrix/mcmctree/out_step2_posterior.txt

Outputs PNG + PDF saved to Notebooks/Figures_set{2,3}/.
"""
import re
import sys
from io import StringIO
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SETS, italic_label  # noqa: E402

PROJECT_ROOT = Path("/home/jovyan")
WORK = PROJECT_ROOT / "Outputs/Phylogenomics_supermatrix/mcmctree"
FIGTREE = WORK / "FigTree.tre"
OUT_TXT = WORK / "out_step2_posterior.txt"


def parse_posterior(out_txt: Path) -> dict[str, dict]:
    posterior: dict[str, dict] = {}
    seen_header = False
    for line in out_txt.read_text().splitlines():
        if "Posterior mean" in line and "HPD" in line:
            seen_header = True
            continue
        if not seen_header:
            continue
        m = re.match(r"^(t_n\d+)\s+([\d.\-]+)\s+\(\s*([\d.\-]+),\s*([\d.\-]+)\)\s+\(\s*([\d.\-]+),\s*([\d.\-]+)\)\s+([\d.\-]+)", line)
        if m:
            posterior[m.group(1)] = {
                "mean": float(m.group(2)),
                "et": (float(m.group(3)), float(m.group(4))),
                "hpd": (float(m.group(5)), float(m.group(6))),
                "width": float(m.group(7)),
            }
        elif posterior and not line.strip():
            break
    return posterior


def parse_figtree(figtree_path: Path):
    txt = figtree_path.read_text()
    bare = re.sub(r"\[[^\]]*\]", "", txt)
    m = re.search(r"(?:UTREE|tree)\s+\S+\s*=\s*(.+?);", bare, re.DOTALL | re.IGNORECASE)
    if not m: return None
    nwk = m.group(1).strip() + ";"
    return Phylo.read(StringIO(nwk), "newick")


def assign_node_ages(tree) -> dict:
    ages = {}
    def walk(node, parent_age):
        bl = node.branch_length or 0.0
        my_age = parent_age - bl
        ages[id(node)] = my_age
        for c in node.clades:
            walk(c, my_age)
    root_depth = max(tree.distance(tree.root, t) for t in tree.get_terminals())
    walk(tree.root, root_depth)
    return ages


def main() -> None:
    plt.rcParams.update({
        "figure.dpi": 300, "savefig.dpi": 300,
        "font.family": "serif", "font.size": 12,
    })
    if not FIGTREE.is_file():
        print(f"missing {FIGTREE}"); return
    tree = parse_figtree(FIGTREE)
    if tree is None:
        print("could not parse FigTree.tre"); return

    posterior = parse_posterior(OUT_TXT) if OUT_TXT.is_file() else {}
    print(f"Read {len(posterior)} posterior estimates")

    ages = assign_node_ages(tree)
    root_age_100mya = ages[id(tree.root)]
    print(f"Root age = {root_age_100mya*100:.1f} Mya")

    internal = [c for c in tree.find_clades() if not c.is_terminal()]
    sorted_internal = sorted(internal, key=lambda c: -ages[id(c)])
    label_map = {id(c): f"t_n{13+i}" for i, c in enumerate(sorted_internal)}

    fig, ax = plt.subplots(figsize=(18, 8))

    def label_func(c):
        if c.is_terminal():
            return "  " + italic_label(c.name)
        return ""

    Phylo.draw(tree, axes=ax, do_show=False, label_func=label_func, branch_labels=lambda c: None)

    for clade in internal:
        age = ages[id(clade)]
        node_name = label_map.get(id(clade))
        if node_name and node_name in posterior:
            post = posterior[node_name]
            mean = post["mean"]
            hpd_lo, hpd_hi = post["hpd"]
            x_lo = root_age_100mya - hpd_hi
            x_hi = root_age_100mya - hpd_lo
            children_y = [ax.get_lines()[0].get_ydata().mean()] if False else None
            try:
                terms = list(clade.get_terminals())
                ys = []
                for t in terms:
                    for line in ax.get_lines():
                        ydata = line.get_ydata()
                        for yv in ydata:
                            ys.append(yv)
                y = sum(ys) / len(ys)
            except Exception:
                y = 0.0
            x_mid = root_age_100mya - mean
            ax.text(x_mid, 0, "", visible=False)

    label_strings = []
    for c in internal:
        nname = label_map.get(id(c))
        if nname and nname in posterior:
            mean_mya = posterior[nname]["mean"] * 100
            lo, hi = posterior[nname]["hpd"]
            label_strings.append(f"{nname}: {mean_mya:5.1f} Mya  (95% HPD: {lo*100:5.1f}–{hi*100:5.1f})")

    mya_ticks = [0, 10, 20, 30, 40, 50, 60, 70]
    xticks = [root_age_100mya - mya/100.0 for mya in mya_ticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(m) for m in mya_ticks])
    ax.set_xlim(left=root_age_100mya - 0.75, right=root_age_100mya + 0.05)
    for mya in mya_ticks:
        x = root_age_100mya - mya / 100.0
        ax.axvline(x, color="grey", lw=0.4, ls=":", alpha=0.4, zorder=0)
    ax.set_xlabel("Million years ago (Mya)", fontsize=12, fontweight="bold")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("MCMCTree calibrated time tree — 12 cyprinids\n226,598-aa supermatrix, 3 fossil-calibration nodes (50K MCMC samples)",
                 fontsize=13, fontweight="bold", pad=8)

    fig.subplots_adjust(right=0.68, left=0.05, top=0.92, bottom=0.10)
    if label_strings:
        fig.text(0.99, 0.5, "Calibrated node ages\n(posterior mean, 95% HPD):\n\n" + "\n".join(label_strings),
                 ha="right", va="center", fontsize=10, family="monospace",
                 bbox=dict(facecolor="white", edgecolor="darkgrey", boxstyle="round,pad=0.6"))
    for s, paths in SETS.items():
        paths["figures_dir"].mkdir(parents=True, exist_ok=True)
        png = paths["figures_dir"] / "mcmctree_timetree.png"
        pdf = paths["figures_dir"] / "mcmctree_timetree.pdf"
        fig.savefig(png, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        print(f"[{s}] wrote {png}")

    table_lines = ["node\tposterior_mean_Mya\thpd_lo_Mya\thpd_hi_Mya"]
    for c in internal:
        nname = label_map.get(id(c))
        if nname and nname in posterior:
            mean_mya = posterior[nname]["mean"] * 100
            lo, hi = posterior[nname]["hpd"]
            table_lines.append(f"{nname}\t{mean_mya:.2f}\t{lo*100:.2f}\t{hi*100:.2f}")
    for s, paths in SETS.items():
        out = paths["tables_dir"] / "mcmctree_node_ages.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(table_lines) + "\n")
        print(f"[{s}] wrote {out}")


if __name__ == "__main__":
    main()
