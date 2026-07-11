#!/usr/bin/env python3
"""Regenerate all synteny circos figures into Submission/SyntenyCircos.

This is a thin, CAFE-free wrapper around the plotting helpers in
``run_mcscanx_pairwise-visualization.py``. It exists because the main script's
``main()`` requires CAFE inputs (family results, change tab, orthogroups) that
are irrelevant to the circos plots, and because we need one extra figure
(``al_dm``) that is not in any config.

What it produces (PDF + PNG each, editable in Illustrator via pdf.fonttype=42):
    synteny_circos_al_al   self-synteny, Aspiorhynchus laticeps
    synteny_circos_dm_dm   self-synteny, Diptychus maculatus
    synteny_circos_al_dr   A. laticeps vs Danio rerio
    synteny_circos_dm_dr   D. maculatus vs D. rerio
    synteny_circos_al_sm   A. laticeps vs Schizothorax macropogon (section 7.2.7)
    synteny_circos_al_dm   A. laticeps vs D. maculatus  (NEW)

``al_dm`` is the symmetric counterpart of the already-computed ``dm_al``
comparison: synteny ribbons are undirected, so we parse the existing
``dm_al.collinearity`` (column order dm, al) and swap each block's endpoints so
that A. laticeps becomes the focal (blue inner) ring, matching the orientation
of the other al-centric figures. No MCScanX/BLASTP recomputation is needed.

Run inside the project container:
    docker exec gpu-jupyter conda run -n env_biopython \\
        python /home/jovyan/scripts/generate_synteny_circos_submission.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Publication defaults (project convention) + editable TrueType fonts. Must be set before
# any figure is drawn. pdf.fonttype=42 keeps text live/selectable in Illustrator
# instead of matplotlib's default Type 3 outlines.
plt.rcParams.update(
    {
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "font.family": "serif",
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.titleweight": "bold",
        "axes.labelsize": 14,
        "axes.labelweight": "bold",
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "legend.title_fontsize": 13,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

# ---------------------------------------------------------------------------
# Locate the project root and import helpers from the hyphenated module.
# This file lives in <root>/scripts/. When run in the container the cwd is
# /home/jovyan, but we resolve relative to this file so it works either way.
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
MODULE_PATH = SCRIPT_DIR / "run_mcscanx_pairwise-visualization.py"

_spec = importlib.util.spec_from_file_location("mcscanx_viz", MODULE_PATH)
if _spec is None or _spec.loader is None:
    sys.exit(f"Error: cannot load module from {MODULE_PATH}")
mcscanx_viz = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mcscanx_viz)

load_chromosome_sizes_from_bed = mcscanx_viz.load_chromosome_sizes_from_bed
parse_synteny_blocks_for_circos = mcscanx_viz.parse_synteny_blocks_for_circos
create_synteny_circos = mcscanx_viz.create_synteny_circos

if not getattr(mcscanx_viz, "PYCIRCLIZE_AVAILABLE", False):
    sys.exit("Error: pyCirclize is not available in this environment.")

OUTPUT_DIR = PROJECT_ROOT / "Submission" / "SyntenyCircos"

# Italic species labels (mathtext, per project convention).
AL = "$\\it{A.\\ laticeps}$"
DM = "$\\it{D.\\ maculatus}$"
DR = "$\\it{D.\\ rerio}$"
SM = "$\\it{S.\\ macropogon}$"


def swap_blocks(blocks: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """Swap focal/comp endpoints so species 2 becomes the focal (chr1) side.

    Synteny ribbons are undirected; swapping the two ends re-orients the figure
    without recomputing anything. Orientation (plus/minus) is preserved.
    """
    swapped: List[Dict[str, object]] = []
    for b in blocks:
        swapped.append(
            {
                "chr1": b["chr2"],
                "start1": b["start2"],
                "end1": b["end2"],
                "chr2": b["chr1"],
                "start2": b["start1"],
                "end2": b["end1"],
                "n_genes": b["n_genes"],
                "orientation": b["orientation"],
            }
        )
    return swapped


def render_pair(
    out_name: str,
    mcscanx_dir: Path,
    col_name: str,
    bed_focal: str,
    prefix_focal: str,
    label_focal: str,
    bed_comp: str,
    prefix_comp: str,
    label_comp: str,
    swap: bool = False,
) -> bool:
    """Parse one comparison and render its synteny circos into OUTPUT_DIR.

    When ``swap`` is True the collinearity columns are (comp, focal) on disk, so
    we parse in disk order then swap the blocks so ``focal`` is the chr1 side.
    """
    col_path = mcscanx_dir / col_name
    bed_focal_path = mcscanx_dir / bed_focal
    bed_comp_path = mcscanx_dir / bed_comp

    missing = [
        str(p)
        for p in (col_path, bed_focal_path, bed_comp_path)
        if not p.exists()
    ]
    if missing:
        print(f"  [SKIP] {out_name}: missing input(s): {', '.join(missing)}")
        return False

    focal_chroms = load_chromosome_sizes_from_bed(bed_focal_path)
    comp_chroms = load_chromosome_sizes_from_bed(bed_comp_path)

    if swap:
        # On disk: chr1 = comp prefix, chr2 = focal prefix.
        blocks = parse_synteny_blocks_for_circos(
            col_path,
            bed_comp_path,
            bed_focal_path,
            chr1_prefix=prefix_comp,
            chr2_prefix=prefix_focal,
            min_genes=10,
        )
        blocks = swap_blocks(blocks)
    else:
        blocks = parse_synteny_blocks_for_circos(
            col_path,
            bed_focal_path,
            bed_comp_path,
            chr1_prefix=prefix_focal,
            chr2_prefix=prefix_comp,
            min_genes=10,
        )

    if not blocks:
        print(f"  [SKIP] {out_name}: no synteny blocks >= 10 genes")
        return False

    out_path = OUTPUT_DIR / f"synteny_circos_{out_name}.png"
    ok = create_synteny_circos(
        focal_chroms=focal_chroms,
        comp_chroms=comp_chroms,
        synteny_blocks=blocks,
        focal_tag=prefix_focal.upper(),
        comp_tag=prefix_comp.upper(),
        focal_label=label_focal,
        comp_label=label_comp,
        output_path=out_path,
    )
    if ok:
        print(
            f"  [DONE] {out_name}: {len(blocks)} blocks -> "
            f"{out_path.name} (+ .pdf)"
        )
    else:
        print(f"  [FAIL] {out_name}: create_synteny_circos returned False")
    return ok


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")

    set3 = PROJECT_ROOT / "Outputs_set3" / "MCScanX_Results"
    db = PROJECT_ROOT / "DB" / "MCScanX"

    # (out_name, dir, collinearity, bed_focal, prefix_focal, label_focal,
    #  bed_comp, prefix_comp, label_comp, swap)
    jobs = [
        (
            "al_al", set3 / "al_al", "al_al.collinearity",
            "al-filtered-modified.bed", "al", AL,
            "al-filtered-modified.bed", "al", AL, False,
        ),
        (
            "dm_dm", set3 / "dm_dm", "dm_dm.collinearity",
            "dm-filtered-modified.bed", "dm", DM,
            "dm-filtered-modified.bed", "dm", DM, False,
        ),
        (
            "al_dr", set3 / "al_dr", "al_dr.collinearity",
            "al-filtered-modified.bed", "al", AL,
            "dr-filtered-modified.bed", "dr", DR, False,
        ),
        (
            "dm_dr", set3 / "dm_dr", "dm_dr.collinearity",
            "dm-filtered-modified.bed", "dm", DM,
            "dr-filtered-modified.bed", "dr", DR, False,
        ),
        (
            "al_sm", db / "al_sm", "al_sm.collinearity",
            "al-filtered-modified.bed", "al", AL,
            "sm-filtered-modified.bed", "sm", SM, False,
        ),
        # NEW: al_dm from the already-computed dm_al data. Despite the directory
        # name, MCScanX wrote the cross-species blocks in al-first column order
        # (al&dm), so al is already the focal (chr1) side -- no swap needed.
        (
            "al_dm", db / "dm_al", "dm_al.collinearity",
            "al-filtered-modified.bed", "al", AL,
            "dm-filtered-modified.bed", "dm", DM, False,
        ),
    ]

    results: Dict[str, bool] = {}
    for job in jobs:
        out_name = job[0]
        print(f"\n--- {out_name} ---")
        results[out_name] = render_pair(
            out_name=job[0],
            mcscanx_dir=job[1],
            col_name=job[2],
            bed_focal=job[3],
            prefix_focal=job[4],
            label_focal=job[5],
            bed_comp=job[6],
            prefix_comp=job[7],
            label_comp=job[8],
            swap=job[9],
        )

    print("\n=== Summary ===")
    done = [k for k, v in results.items() if v]
    skipped = [k for k, v in results.items() if not v]
    print(f"  Generated ({len(done)}): {', '.join(done) if done else '-'}")
    print(f"  Skipped   ({len(skipped)}): {', '.join(skipped) if skipped else '-'}")
    if skipped:
        sys.exit(1)


if __name__ == "__main__":
    main()
