#!/usr/bin/env python3
"""Render the IQ-TREE 3A supermatrix tree as a publication figure.

Re-roots on Danio rerio (basal cyprinid in this dataset) and labels internal
nodes with the triple-support string SH-aLRT / UFBoot / sCFL (%) emitted by
IQ-TREE 3 when run with `-bnni -alrt 1000 -B 1000` followed by `--scfl 100`.

Produces PNG + PDF saved to Notebooks/Figures_set{2,3}/.
"""
import sys
from pathlib import Path
from io import StringIO
import matplotlib.pyplot as plt
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SETS, italic_label  # noqa: E402

PROJECT_ROOT = Path("/home/jovyan")
TREEFILE = PROJECT_ROOT / "Outputs/Phylogenomics_supermatrix/iqtree_run.scfl.cf.tree"
ROOT_TAXON = "Danio_rerio"


def label_func(clade):
    if clade.is_terminal():
        name = clade.name or ""
        return italic_label(name)
    return ""


def support_label(clade):
    if clade.is_terminal():
        return None
    raw = clade.name
    if not raw or "/" not in raw:
        return None
    parts = raw.split("/")
    if len(parts) == 3:
        shalrt, ufboot, scfl = parts
        try:
            scfl_f = float(scfl)
        except ValueError:
            return raw
        return f"{shalrt}/{ufboot}/{scfl_f:.0f}"
    return raw


def main() -> None:
    plt.rcParams.update({
        "figure.dpi": 300, "savefig.dpi": 300,
        "font.family": "serif", "font.size": 12,
    })

    nwk_text = TREEFILE.read_text()
    tree = Phylo.read(StringIO(nwk_text), "newick")
    print(f"Read tree with {len(tree.get_terminals())} terminals")

    target = next((t for t in tree.get_terminals() if t.name == ROOT_TAXON), None)
    if target is not None:
        tree.root_with_outgroup(target)
        print(f"Rooted on {ROOT_TAXON}")

    fig, ax = plt.subplots(figsize=(13, 8.5))
    Phylo.draw(
        tree, axes=ax, do_show=False,
        label_func=label_func,
        branch_labels=support_label,
    )
    ax.set_title(
        "IQ-TREE 3 partitioned supermatrix phylogeny\n"
        "(500 OGs, 226,598 aa columns; SH-aLRT / UFBoot / sCFL on internal nodes)",
        fontsize=13, fontweight="bold", pad=12,
    )
    ax.set_xlabel("Substitutions per site", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_yticks([])
    fig.tight_layout()

    for s, paths in SETS.items():
        paths["figures_dir"].mkdir(parents=True, exist_ok=True)
        png = paths["figures_dir"] / "supermatrix_iqtree_rendered.png"
        pdf = paths["figures_dir"] / "supermatrix_iqtree_rendered.pdf"
        fig.savefig(png, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        print(f"[{s}] wrote {png}")
        print(f"[{s}] wrote {pdf}")


if __name__ == "__main__":
    main()
