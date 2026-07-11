"""Shared helper: canonical species roster + per-set output paths.

Reflects Notebooks/config/paths-set2.yaml and paths-set3.yaml.
Both sets share the same 12-species roster (Mar03 OrthoFinder run); they
differ only in focal species and output directories.
"""
from pathlib import Path

PROJECT_ROOT = Path("/home/jovyan")

SPECIES = [
    "Aspiorhynchus_laticeps",
    "Carassius_auratus",
    "Cyprinus_carpio",
    "Danio_rerio",
    "Diptychus_maculatus",
    "Gymnocypris_eckloni",
    "Oxygymnocypris_stewartii",
    "Schizopygopsis_younghusbandi",
    "Sinocyclocheilus_grahami",
    "Triplophysa_pappenheimi",
    "Triplophysa_tibetana",
    "Triplophysa_yaopeizhii",
]

OF_RESULTS = PROJECT_ROOT / "Outputs/OrthoFinder/Results_Mar03"
BRAKER_MASKED = PROJECT_ROOT / "Outputs/Preprocessing/BRAKER_MASKED"

SETS = {
    "set2": {"focal": "Aspiorhynchus_laticeps",
             "figures_dir": PROJECT_ROOT / "Notebooks/Figures_set2",
             "tables_dir": PROJECT_ROOT / "Notebooks/Tables_set2",
             "outputs_dir": PROJECT_ROOT / "Outputs_set2"},
    "set3": {"focal": "Diptychus_maculatus",
             "figures_dir": PROJECT_ROOT / "Notebooks/Figures_set3",
             "tables_dir": PROJECT_ROOT / "Notebooks/Tables_set3",
             "outputs_dir": PROJECT_ROOT / "Outputs_set3"},
}


def braker_gff(species: str) -> Path:
    return BRAKER_MASKED / species / "braker.gff3"


def braker_cds(species: str) -> Path:
    return BRAKER_MASKED / species / "braker.codingseq"


def italic_label(species: str, abbrev: bool = False) -> str:
    parts = species.split("_")
    if abbrev and len(parts) == 2:
        parts = [parts[0][0] + ".", parts[1]]
    return r"$\it{" + r"\ ".join(parts) + r"}$"


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def save_to_both_sets(save_fn, basename: str, kind: str = "figure") -> list[Path]:
    """Call save_fn(target_path) for the same artifact in both set2 and set3 dirs."""
    written = []
    for s in ("set2", "set3"):
        d = SETS[s][f"{kind}s_dir"] if kind != "figure" else SETS[s]["figures_dir"]
        d.mkdir(parents=True, exist_ok=True)
        target = d / basename
        save_fn(target)
        written.append(target)
    return written


# =============================================================================
# Publication style helpers (shared across all plot_*.py scripts)
# =============================================================================
# Centralises the manuscript figure standards so every script is consistent:
#   * large fonts (publication minimums) that survive journal size reduction
#   * EDITABLE vector text in PDF/SVG (pdf.fonttype=42 -> embedded TrueType,
#     svg.fonttype='none' -> text stays as <text>, not outlined paths). This is
#     the fix for the PI's recurring "PDF text/legend not editable" complaint.

PUB_RCPARAMS = {
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
    "font.family": "serif",
    "font.size": 14,            # base
    "axes.titlesize": 16,       # title
    "axes.titleweight": "bold",
    "axes.labelsize": 14,       # axis labels
    "axes.labelweight": "bold",
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "legend.title_fontsize": 13,
    # --- editable vector text ---
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def apply_pub_rcparams(extra: dict | None = None) -> None:
    """Apply the project publication rcParams (large fonts + editable PDF/SVG).

    Call once near the top of a plotting script. Pass ``extra`` to override or
    add keys (e.g. ``apply_pub_rcparams({'font.family': 'sans-serif'})``).
    """
    import matplotlib.pyplot as plt  # local import: keeps this module light

    params = dict(PUB_RCPARAMS)
    if extra:
        params.update(extra)
    plt.rcParams.update(params)


def save_fig(fig, figures_dir: Path, basename: str,
             formats: tuple[str, ...] = ("png", "pdf", "svg"),
             close: bool = False, **savefig_kwargs) -> list[Path]:
    """Save ``fig`` as ``basename.{png,pdf,svg}`` into ``figures_dir``.

    PNG is raster (300 dpi); PDF/SVG carry editable vector text thanks to the
    rcParams set by :func:`apply_pub_rcparams`. Returns the written paths.
    """
    import matplotlib.pyplot as plt

    figures_dir.mkdir(parents=True, exist_ok=True)
    kwargs = {"bbox_inches": "tight"}
    kwargs.update(savefig_kwargs)
    written = []
    for ext in formats:
        target = figures_dir / f"{basename}.{ext}"
        fig.savefig(target, **kwargs)
        written.append(target)
    if close:
        plt.close(fig)
    return written


def save_fig_both_sets(fig, basename: str,
                       formats: tuple[str, ...] = ("png", "pdf", "svg"),
                       close: bool = False, **savefig_kwargs) -> list[Path]:
    """Save the same set-agnostic figure into BOTH set2 and set3 figure dirs."""
    written = []
    for s in ("set2", "set3"):
        written += save_fig(fig, SETS[s]["figures_dir"], basename,
                            formats=formats, close=False, **savefig_kwargs)
    if close:
        import matplotlib.pyplot as plt
        plt.close(fig)
    return written
