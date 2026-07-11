#!/usr/bin/env python3
"""Per-species gene-structure statistics from BRAKER_MASKED braker.gff3.

For each species in set2/set3 (same 12-species roster), compute:
  - Number of genes (gene records)
  - Number of mRNAs
  - Mean gene length (bp)
  - Mean CDS length per gene (sum of CDS records, per longest mRNA)
  - Mean exon length (bp)
  - Mean intron length (bp)
  - Mean exons per gene (longest mRNA)

Output: one CSV written to BOTH Tables_set2 and Tables_set3.
"""
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SPECIES, SETS, braker_gff  # noqa: E402


def parse_gff_stats(gff_path: Path) -> dict[str, float]:
    gene_lengths: list[int] = []
    mrna_to_gene: dict[str, str] = {}
    mrna_cds_lengths: dict[str, list[int]] = defaultdict(list)
    mrna_exon_lengths: dict[str, list[int]] = defaultdict(list)
    mrna_intron_lengths: dict[str, list[int]] = defaultdict(list)
    n_mrna = 0

    with gff_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            try:
                start = int(f[3]); end = int(f[4])
            except ValueError:
                continue
            length = end - start + 1
            attrs = dict(p.split("=", 1) for p in f[8].rstrip(";").split(";") if "=" in p)

            if ftype == "gene":
                gene_lengths.append(length)
            elif ftype == "mRNA":
                n_mrna += 1
                mid = attrs.get("ID")
                pid = attrs.get("Parent")
                if mid and pid:
                    mrna_to_gene[mid] = pid
            elif ftype == "CDS":
                pid = attrs.get("Parent")
                if pid:
                    mrna_cds_lengths[pid].append(length)
            elif ftype == "exon":
                pid = attrs.get("Parent")
                if pid:
                    mrna_exon_lengths[pid].append(length)
            elif ftype == "intron":
                pid = attrs.get("Parent")
                if pid:
                    mrna_intron_lengths[pid].append(length)

    gene_to_longest_mrna: dict[str, str] = {}
    mrna_total_cds = {m: sum(L) for m, L in mrna_cds_lengths.items()}
    for mid, gid in mrna_to_gene.items():
        cur = gene_to_longest_mrna.get(gid)
        if cur is None or mrna_total_cds.get(mid, 0) > mrna_total_cds.get(cur, 0):
            gene_to_longest_mrna[gid] = mid

    longest_mrnas = list(gene_to_longest_mrna.values())
    cds_lengths_per_gene = [mrna_total_cds.get(m, 0) for m in longest_mrnas]
    exons_per_gene = [len(mrna_exon_lengths.get(m, [])) for m in longest_mrnas]

    all_exon_lengths = [L for m in longest_mrnas for L in mrna_exon_lengths.get(m, [])]
    all_intron_lengths = [L for m in longest_mrnas for L in mrna_intron_lengths.get(m, [])]

    n_genes = len(gene_lengths)

    def mean(xs):
        return float(sum(xs) / len(xs)) if xs else 0.0

    return {
        "Number of genes": n_genes,
        "Number of mRNAs": n_mrna,
        "Mean gene length (bp)": round(mean(gene_lengths), 1),
        "Mean CDS length per gene (bp)": round(mean(cds_lengths_per_gene), 1),
        "Mean exon length (bp)": round(mean(all_exon_lengths), 1),
        "Mean intron length (bp)": round(mean(all_intron_lengths), 1),
        "Mean exons per gene": round(mean(exons_per_gene), 2),
    }


def main() -> None:
    rows = []
    for sp in SPECIES:
        gff = braker_gff(sp)
        if not gff.is_file():
            print(f"[skip] {sp}: {gff} missing")
            continue
        print(f"[parse] {sp}: {gff}")
        stats = parse_gff_stats(gff)
        stats["Species"] = sp
        rows.append(stats)

    df = pd.DataFrame(rows).set_index("Species")
    cols = [
        "Number of genes", "Number of mRNAs",
        "Mean gene length (bp)", "Mean CDS length per gene (bp)",
        "Mean exon length (bp)", "Mean intron length (bp)",
        "Mean exons per gene",
    ]
    df = df[cols]
    print()
    print(df.to_string())

    for s, paths in SETS.items():
        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        out = paths["tables_dir"] / "gene_structure_stats.csv"
        df.to_csv(out)
        print(f"[{s}] wrote {out}")


if __name__ == "__main__":
    main()
