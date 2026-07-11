#!/usr/bin/env python3
"""Tier 2C: Multi-track circos for the focal species (AL = set2, DM = set3).

Tracks (outer → inner):
  a. Chromosome ideogram (top-25 sequences by length)
  b. Gene density (genes per Mb in 1-Mb bins)
  c. Repeat density (% repeat per Mb)
  d. GC ratio (% per Mb)

ncRNA track will be added in a follow-up once Tier 2A finishes.
"""
import sys
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
from Bio import SeqIO
from pycirclize import Circos

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SETS, italic_label  # noqa: E402

PROJECT_ROOT = Path("/home/jovyan")
WINDOW = 1_000_000

FOCAL = {
    "Aspiorhynchus_laticeps": {
        "fna": PROJECT_ROOT / "DB/GENOME_Comparison/Aspiorhynchus_laticeps/ncbi_dataset/data/GCA_023376895.1/GCA_023376895.1_ASM2337689v1_genomic.fna",
        "gene_gff": PROJECT_ROOT / "DB/GENOME_Comparison/Aspiorhynchus_laticeps/AL.chromosome.gff3",
        "repeat_gff": PROJECT_ROOT / "Outputs/RepeatMasker_chromosomal/Aspiorhynchus_laticeps/GCA_023376895.1_ASM2337689v1_genomic.fna.out.gff",
        "set": "set2",
        "remap_gff_to_fna": "by_size_order",
        "repeat_remap": "case_insensitive",
    },
    "Diptychus_maculatus": {
        "fna": PROJECT_ROOT / "DB/GENOME_Comparison/Diptychus_maculatus/DM.fasta",
        "gene_gff": PROJECT_ROOT / "DB/GENOME_Comparison/Diptychus_maculatus/DM_final.evm.gff3",
        "repeat_gff": PROJECT_ROOT / "DB/GENOME_Comparison/Diptychus_maculatus/DM.fa.out.gff",
        "set": "set3",
        "remap_gff_to_fna": "case_insensitive",
    },
}


def top_n_fna_seqs(fna: Path, n: int = 25, min_mb: float = 5.0):
    items = [(r.id, len(r.seq)) for r in SeqIO.parse(str(fna), "fasta")]
    items.sort(key=lambda x: -x[1])
    sel = [(rid, L) for rid, L in items[:n] if L >= min_mb * 1e6]
    return sel


def build_remap(top_fna: list[tuple[str, int]], gff_chrs: set[str], how: str) -> dict[str, str]:
    """Return dict: gff_chr_id → canonical FNA id (string), for chrs we keep."""
    fna_ids = [r for r, _ in top_fna]

    if how == "case_insensitive":
        fna_lc = {r.lower(): r for r in fna_ids}
        return {gc: fna_lc[gc.lower()] for gc in gff_chrs if gc.lower() in fna_lc}

    if how == "by_size_order":
        chr_pat = sorted([c for c in gff_chrs if c.lower().startswith("chr") and c[3:].isdigit()],
                         key=lambda s: int(s[3:]))
        return {gc: fna_id for gc, fna_id in zip(chr_pat, fna_ids)}

    return {gc: fna_id for gc, fna_id in zip(sorted(gff_chrs), fna_ids)}


def per_chr_bin_counts(gff: Path, accept: dict[str, str], window: int) -> dict[str, dict[int, int]]:
    counts = defaultdict(lambda: defaultdict(int))
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[2] != "gene": continue
            ch = f[0]
            if ch not in accept: continue
            try: start = int(f[3])
            except ValueError: continue
            counts[accept[ch]][start // window] += 1
    return counts


def per_chr_bin_repeat_bp(gff: Path, accept: dict[str, str], window: int, max_bp_per_chr: dict[str, int]) -> dict[str, dict[int, int]]:
    bp = defaultdict(lambda: defaultdict(int))
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5: continue
            ch = f[0]
            if ch not in accept: continue
            try:
                s = int(f[3]); e = int(f[4])
            except ValueError:
                continue
            fid = accept[ch]
            chr_max = max_bp_per_chr.get(fid, e)
            s = max(1, min(s, chr_max)); e = max(1, min(e, chr_max))
            for w in range(s // window, e // window + 1):
                lo = max(s, w * window); hi = min(e, (w + 1) * window)
                bp[fid][w] += max(0, hi - lo + 1)
    return bp


def per_chr_gc(fna: Path, fna_ids: list[str], window: int) -> dict[str, list[float]]:
    keep = set(fna_ids)
    out: dict[str, list[float]] = {fid: [] for fid in fna_ids}
    for rec in SeqIO.parse(str(fna), "fasta"):
        if rec.id not in keep: continue
        seq = str(rec.seq).upper()
        nbins = (len(seq) + window - 1) // window
        bins = []
        for i in range(nbins):
            chunk = seq[i*window:(i+1)*window]
            gc = chunk.count("G") + chunk.count("C")
            ats = chunk.count("A") + chunk.count("T") + gc
            bins.append(100.0 * gc / ats if ats > 0 else 0.0)
        out[rec.id] = bins
    return out


def safe_mid(w: int, chr_len: int) -> float:
    midpoint = w * WINDOW + WINDOW / 2
    return min(midpoint, chr_len - 1)


def build_circos(focal: str, info: dict) -> None:
    print(f"=== {focal} ===")
    top = top_n_fna_seqs(info["fna"], n=25, min_mb=5.0)
    print(f"  picked {len(top)} chromosomes from FNA: max={top[0][1]/1e6:.1f} Mb, total={sum(L for _,L in top)/1e6:.1f} Mb")
    fna_ids = [r for r, _ in top]
    fna_lens = dict(top)

    gff_chrs = set()
    with info["gene_gff"].open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f = line.split("\t")
            if len(f) > 2 and f[2] == "gene":
                gff_chrs.add(f[0])

    accept = build_remap(top, gff_chrs, info["remap_gff_to_fna"])
    print(f"  gene gff→fna mapping covers {len(accept)} of {len(gff_chrs)} GFF chrs")
    if not accept:
        print(f"  [error] {focal}: no GFF chrs map to FNA top-N. GFF examples: {list(gff_chrs)[:5]}; FNA examples: {fna_ids[:3]}")
        return

    repeat_chrs = set()
    with info["repeat_gff"].open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f = line.split("\t")
            if len(f) > 0:
                repeat_chrs.add(f[0])
    repeat_remap_mode = info.get("repeat_remap", info["remap_gff_to_fna"])
    repeat_accept = build_remap(top, repeat_chrs, repeat_remap_mode)
    print(f"  repeat gff→fna mapping covers {len(repeat_accept)} of {len(repeat_chrs)} repeat-GFF chrs")

    print("  computing gene density ...")
    gd = per_chr_bin_counts(info["gene_gff"], accept, WINDOW)
    print("  computing repeat density ...")
    rd = per_chr_bin_repeat_bp(info["repeat_gff"], repeat_accept, WINDOW, fna_lens)
    print("  computing GC ratio ...")
    gc = per_chr_gc(info["fna"], fna_ids, WINDOW)

    sectors = {fid: fna_lens[fid] for fid in fna_ids}
    circos = Circos(sectors=sectors, space=2)
    circos.text(italic_label(focal), size=18, r=10, deg=270)

    palette = {"gene": "#2c7bb6", "rep": "#d7191c", "gc": "#1b9e77"}

    for sector in circos.sectors:
        fid = sector.name
        L = fna_lens[fid]
        sector.axis(fc="lightgrey", lw=0.6)
        label = fid.split(".")[0] if fid.startswith("CM") else fid
        sector.text(label, r=110, size=8, color="black")

        nbins = max(1, (L + WINDOW - 1) // WINDOW)
        x = [safe_mid(w, L) for w in range(nbins)]

        gene_y = [gd.get(fid, {}).get(w, 0) for w in range(nbins)]
        track_g = sector.add_track((78, 95), r_pad_ratio=0.05)
        track_g.axis(fc="white", ec="grey", lw=0.3)
        if any(gene_y):
            track_g.bar(x, gene_y, width=WINDOW * 0.9, color=palette["gene"], align="center")

        rep_y = [100.0 * rd.get(fid, {}).get(w, 0) / WINDOW for w in range(nbins)]
        track_r = sector.add_track((58, 75), r_pad_ratio=0.05)
        track_r.axis(fc="white", ec="grey", lw=0.3)
        if any(rep_y):
            track_r.bar(x, rep_y, width=WINDOW * 0.9, color=palette["rep"], align="center")

        gc_y = gc.get(fid, [])
        track_c = sector.add_track((38, 55), r_pad_ratio=0.05)
        track_c.axis(fc="white", ec="grey", lw=0.3)
        if gc_y:
            n = min(len(gc_y), len(x))
            track_c.line(x[:n], gc_y[:n], color=palette["gc"], lw=0.8)

    fig = circos.plotfig(figsize=(12, 12), dpi=300)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=palette["gene"], label="Gene density (genes / Mb)"),
        plt.Rectangle((0, 0), 1, 1, color=palette["rep"], label="Repeat density (% / Mb)"),
        plt.Line2D([0], [0], color=palette["gc"], lw=2, label="GC ratio (%)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=11, frameon=False, bbox_to_anchor=(0.5, 0.02))
    fig.suptitle(f"Genome circos: {italic_label(focal)}", y=0.97, fontsize=15, fontweight="bold")

    s = info["set"]
    out_dir = SETS[s]["figures_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"focal_circos_{s}.png"
    pdf = out_dir / f"focal_circos_{s}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")
    print(f"  wrote {pdf}")


def main() -> None:
    plt.rcParams.update({"font.family": "serif", "font.size": 12})
    for focal, info in FOCAL.items():
        for k in ("fna", "gene_gff", "repeat_gff"):
            if not info[k].is_file():
                print(f"[skip] {focal}: missing {k}={info[k]}"); break
        else:
            try:
                build_circos(focal, info)
            except Exception as exc:
                print(f"[error] {focal}: {exc}")
                import traceback; traceback.print_exc()


if __name__ == "__main__":
    main()
