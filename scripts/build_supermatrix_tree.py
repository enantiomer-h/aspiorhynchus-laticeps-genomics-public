#!/usr/bin/env python3
"""Build a partitioned protein supermatrix from near-single-copy OGs and run IQ-TREE 2.

Strategy:
  1. From Orthogroups.GeneCount.tsv, select OGs that are "near-single-copy":
       count == 1 in >=10 of the 12 species (allows 2 missing/extra-copy species).
  2. For each surviving OG, read its MSA from MultipleSequenceAlignments/.
  3. Pick exactly one sequence per species (longest non-gap content if multiple).
  4. Gap-pad species missing from an OG to the OG's alignment length.
  5. Concatenate across OGs into a per-species supermatrix.
  6. Run IQ-TREE 2 with `-m MFP -B 1000` for ML tree + ultrafast bootstrap.

Output goes to BOTH set2 and set3 dirs since the input species roster is identical.
"""
import sys
import subprocess
from pathlib import Path
import pandas as pd
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import OF_RESULTS, SETS, SPECIES, PROJECT_ROOT  # noqa: E402

GENECOUNT_TSV = OF_RESULTS / "Orthogroups/Orthogroups.GeneCount.tsv"
MSA_DIR = OF_RESULTS / "MultipleSequenceAlignments"

WORK_DIR = PROJECT_ROOT / "Outputs/Phylogenomics_supermatrix"
SUPERMATRIX_FA = WORK_DIR / "supermatrix.faa"
IQTREE_PREFIX = WORK_DIR / "iqtree_run"
MIN_SPECIES_PRESENT = 10
MAX_OGS = 500  # cap for runtime


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    name, seq_chunks = None, []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_chunks)))
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            records.append((name, "".join(seq_chunks)))
    return records


def species_of(header: str) -> str | None:
    for sp in SPECIES:
        if header.startswith(sp + "_"):
            return sp
    return None


def main() -> None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Project species ({len(SPECIES)}):", SPECIES)

    print("Filtering OGs from GeneCount.tsv ...")
    gc = pd.read_csv(GENECOUNT_TSV, sep="\t", index_col=0)
    if "Total" in gc.columns:
        gc = gc.drop(columns=["Total"])
    gc = gc[SPECIES]
    eq1 = (gc == 1).sum(axis=1)
    candidate_ogs = gc.index[eq1 >= MIN_SPECIES_PRESENT].tolist()
    print(f"OGs with >={MIN_SPECIES_PRESENT}/{len(SPECIES)} species at exactly 1 copy: {len(candidate_ogs)}")

    if len(candidate_ogs) > MAX_OGS:
        scores = eq1[candidate_ogs].sort_values(ascending=False)
        candidate_ogs = scores.head(MAX_OGS).index.tolist()
        print(f"Capped to top {MAX_OGS} OGs by species coverage")

    per_species_chunks: dict[str, list[str]] = {sp: [] for sp in SPECIES}
    partitions: list[tuple[str, int, int]] = []
    cursor = 0
    n_used = 0
    n_skipped_no_msa = 0
    for og in candidate_ogs:
        msa_path = MSA_DIR / f"{og}.fa"
        if not msa_path.is_file():
            n_skipped_no_msa += 1
            continue
        records = read_fasta(msa_path)
        if not records:
            continue
        aln_len = len(records[0][1])

        per_sp_seqs: dict[str, list[tuple[int, str]]] = defaultdict(list)
        for hdr, seq in records:
            sp = species_of(hdr)
            if sp is None:
                continue
            ungapped = sum(1 for c in seq if c not in ("-", "?", "X", "*"))
            per_sp_seqs[sp].append((ungapped, seq))

        chosen: dict[str, str] = {}
        for sp in SPECIES:
            if per_sp_seqs[sp]:
                per_sp_seqs[sp].sort(reverse=True)
                chosen[sp] = per_sp_seqs[sp][0][1].replace("*", "-")
            else:
                chosen[sp] = "-" * aln_len

        present_species = sum(1 for sp in SPECIES if sum(1 for c in chosen[sp] if c != "-") > 0)
        if present_species < MIN_SPECIES_PRESENT:
            continue

        for sp in SPECIES:
            per_species_chunks[sp].append(chosen[sp])
        partitions.append((og, cursor + 1, cursor + aln_len))
        cursor += aln_len
        n_used += 1

    print(f"Used {n_used} OGs in supermatrix; skipped {n_skipped_no_msa} (no MSA on disk)")
    if n_used == 0:
        sys.exit("ERROR: zero OGs survived filters")

    print(f"Supermatrix length: {cursor} aa columns")
    with SUPERMATRIX_FA.open("w") as out:
        for sp in SPECIES:
            seq = "".join(per_species_chunks[sp])
            out.write(f">{sp}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")
    print(f"Wrote supermatrix to {SUPERMATRIX_FA}")

    part_path = WORK_DIR / "partitions.nex"
    with part_path.open("w") as out:
        out.write("#nexus\nbegin sets;\n")
        for i, (og, s, e) in enumerate(partitions, 1):
            out.write(f"  charset part{i} = {s}-{e};\n")
        out.write("end;\n")
    print(f"Wrote partitions to {part_path}")

    pd.DataFrame(partitions, columns=["OG", "start", "end"]).to_csv(
        WORK_DIR / "partitions_table.csv", index=False)

    cmd = [
        "iqtree3",
        "-s", str(SUPERMATRIX_FA),
        "-p", str(part_path),
        "-m", "MFP",
        "-B", "1000",
        "-bnni",
        "-alrt", "1000",
        "-T", "AUTO",
        "--prefix", str(IQTREE_PREFIX),
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    treefile = Path(str(IQTREE_PREFIX) + ".treefile")
    if not treefile.is_file():
        sys.exit(f"ERROR: IQ-TREE did not produce {treefile}")
    print(f"IQ-TREE done. Tree at {treefile}")

    for s, paths in SETS.items():
        paths["figures_dir"].mkdir(parents=True, exist_ok=True)
        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        dest_tree = paths["figures_dir"] / "supermatrix_iqtree.treefile"
        dest_tree.write_text(treefile.read_text())
        print(f"[{s}] copied tree to {dest_tree}")
        log_src = Path(str(IQTREE_PREFIX) + ".iqtree")
        if log_src.is_file():
            (paths["tables_dir"] / "supermatrix_iqtree.log.txt").write_text(log_src.read_text())


if __name__ == "__main__":
    main()
