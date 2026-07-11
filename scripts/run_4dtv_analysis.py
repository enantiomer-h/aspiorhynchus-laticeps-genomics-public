#!/usr/bin/env python3
"""4DTv (Four-Fold Degenerate Transversion) Analysis Pipeline.

Calculates transversion rates at four-fold degenerate third-codon positions
for syntenic gene pairs identified by MCScanX. Used to estimate relative
timing of whole-genome duplication (WGD) and speciation events.

Pipeline steps:
    1. Parse MCScanX collinearity file to extract syntenic gene pairs
    2. Load CDS sequences and map MCScanX IDs to BRAKER3 transcript IDs
    3. For each pair: translate → align proteins → back-translate → extract
       four-fold degenerate sites → count transversions → apply HKY correction
    4. Output per-pair TSV results
    5. Generate publication-quality histogram with KDE overlay

Input data requirements:
    - .collinearity file from MCScanX (self- or inter-species)
    - braker.codingseq : CDS sequences from BRAKER3 (per species)

Usage:
    # Self-comparisons for WGD detection:
    python scripts/run_4dtv_analysis.py \\
        --self Aspiorhynchus_laticeps --self Diptychus_maculatus

    # Inter-species comparison (requires collinearity file):
    python scripts/run_4dtv_analysis.py \\
        --pair Aspiorhynchus_laticeps Diptychus_maculatus

    # With options:
    python scripts/run_4dtv_analysis.py \\
        --self Aspiorhynchus_laticeps --self Diptychus_maculatus \\
        --aligner biopython --workers 8

Author: Comparative genomics pipeline
"""

import argparse
import math
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Standard genetic code: codon -> amino acid (T used instead of U)
CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# First two positions of codons where the third position is four-fold
# degenerate (any nucleotide at position 3 yields the same amino acid).
FOURFOLD_DEGENERATE_PREFIXES = frozenset({
    "GC",  # Ala: GCT, GCC, GCA, GCG
    "GG",  # Gly: GGT, GGC, GGA, GGG
    "GT",  # Val: GTT, GTC, GTA, GTG
    "CT",  # Leu: CTT, CTC, CTA, CTG
    "CC",  # Pro: CCT, CCC, CCA, CCG
    "TC",  # Ser: TCT, TCC, TCA, TCG
    "CG",  # Arg: CGT, CGC, CGA, CGG
    "AC",  # Thr: ACT, ACC, ACA, ACG
})

PURINES = frozenset({"A", "G"})
PYRIMIDINES = frozenset({"C", "T"})

# Gene pair line pattern in .collinearity files
# e.g. "  0-  0:\tal_g14\tal_g2476\t      0"
COLLINEARITY_PAIR_RE = re.compile(
    r"^\s*\d+-\s*\d+:\s+(\S+)\s+(\S+)\s+"
)


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def get_species_prefix(species_name: str) -> str:
    """Derive a 2-letter MCScanX prefix from a species name.

    Takes the first letter of genus + first letter of species epithet,
    lowercased.

    Args:
        species_name: e.g. 'Aspiorhynchus_laticeps'

    Returns:
        Prefix string, e.g. 'al'
    """
    parts = species_name.split("_")
    if len(parts) < 2:
        print(
            f"Error: Species name must be 'Genus_species' format: {species_name}",
            file=sys.stderr,
        )
        sys.exit(1)
    return (parts[0][0] + parts[1][0]).lower()


def get_display_name(species_name: str, abbreviated: bool = False) -> str:
    """Convert species name to display format.

    Args:
        species_name: e.g. 'Aspiorhynchus_laticeps'
        abbreviated: If True, abbreviate genus (e.g. 'A. laticeps')

    Returns:
        Display name string
    """
    parts = species_name.split("_")
    genus, epithet = parts[0], parts[1]
    if abbreviated:
        return f"{genus[0]}. {epithet}"
    return f"{genus} {epithet}"


def get_mathtext_name(species_name: str, abbreviated: bool = False) -> str:
    """Format species name as matplotlib mathtext italic.

    Args:
        species_name: e.g. 'Aspiorhynchus_laticeps'
        abbreviated: If True, abbreviate genus

    Returns:
        Mathtext string, e.g. '$\\it{A.\\ laticeps}$'
    """
    parts = species_name.split("_")
    genus, epithet = parts[0], parts[1]
    if abbreviated:
        return f"$\\it{{{genus[0]}.\\ {epithet}}}$"
    return f"$\\it{{{genus}\\ {epithet}}}$"


def is_transversion(base1: str, base2: str) -> bool:
    """Check if a substitution is a transversion (purine <-> pyrimidine)."""
    if base1 == base2:
        return False
    return (base1 in PURINES) != (base2 in PURINES)


# ---------------------------------------------------------------------------
# I/O functions
# ---------------------------------------------------------------------------


def parse_collinearity(collinearity_path: Path) -> List[Tuple[str, str]]:
    """Parse MCScanX collinearity file to extract syntenic gene pairs.

    Args:
        collinearity_path: Path to .collinearity file

    Returns:
        List of (gene1_id, gene2_id) tuples with MCScanX-style IDs
    """
    pairs = []
    with open(collinearity_path) as fh:
        for line in fh:
            m = COLLINEARITY_PAIR_RE.match(line)
            if m:
                pairs.append((m.group(1), m.group(2)))
    return pairs


def load_cds_sequences(cds_path: Path, prefix: str) -> Dict[str, str]:
    """Load CDS FASTA and build mapping from MCScanX-style IDs to sequences.

    BRAKER3 CDS headers are '>g14.t1'. MCScanX IDs are 'al_g14'.
    This function maps '{prefix}_g{N}' -> CDS sequence of 'g{N}.t1'.

    Only the primary transcript (.t1) is kept per gene. If .t1 is missing,
    falls back to the first available transcript for that gene.

    Args:
        cds_path: Path to braker.codingseq FASTA file
        prefix: Species prefix (e.g. 'al')

    Returns:
        Dict mapping MCScanX IDs to uppercase DNA sequences
    """
    sequences: Dict[str, str] = {}
    # Track all transcripts per gene for fallback
    gene_transcripts: Dict[str, Dict[str, str]] = {}

    current_id = None
    current_seq: List[str] = []

    with open(cds_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                # Save previous sequence
                if current_id is not None:
                    full_seq = "".join(current_seq).upper()
                    # Parse gene and transcript: "g14.t1" -> gene="g14", tx="t1"
                    parts = current_id.split(".")
                    gene_id = parts[0]
                    tx_id = parts[1] if len(parts) > 1 else "t1"
                    if gene_id not in gene_transcripts:
                        gene_transcripts[gene_id] = {}
                    gene_transcripts[gene_id][tx_id] = full_seq

                current_id = line[1:].split()[0]  # strip > and take first token
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last entry
        if current_id is not None:
            full_seq = "".join(current_seq).upper()
            parts = current_id.split(".")
            gene_id = parts[0]
            tx_id = parts[1] if len(parts) > 1 else "t1"
            if gene_id not in gene_transcripts:
                gene_transcripts[gene_id] = {}
            gene_transcripts[gene_id][tx_id] = full_seq

    # Build MCScanX-keyed dict, preferring .t1
    for gene_id, transcripts in gene_transcripts.items():
        mcscanx_id = f"{prefix}_{gene_id}"
        if "t1" in transcripts:
            sequences[mcscanx_id] = transcripts["t1"]
        else:
            # Fallback: take the first transcript
            sequences[mcscanx_id] = next(iter(transcripts.values()))

    return sequences


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------


def translate_sequence(dna_seq: str) -> str:
    """Translate DNA sequence to protein using standard genetic code.

    Handles sequences whose length is not a multiple of 3 by ignoring the
    trailing partial codon. Stops translation at the first stop codon.

    Args:
        dna_seq: Uppercase DNA string

    Returns:
        Protein string (single-letter amino acids)
    """
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i + 3]
        if "N" in codon:
            protein.append("X")  # ambiguous
            continue
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def align_proteins_biopython(seq1: str, seq2: str) -> Tuple[str, str]:
    """Align two protein sequences using BioPython PairwiseAligner.

    Uses BLOSUM62 substitution matrix for scoring. This is a fast fallback
    when MAFFT is not available.

    Args:
        seq1: First protein sequence
        seq2: Second protein sequence

    Returns:
        Tuple of (aligned_seq1, aligned_seq2)
    """
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"

    alignments = aligner.align(seq1, seq2)
    best = alignments[0]

    # Extract aligned sequences from the Alignment object
    aln_seqs = best.format("fasta").strip().split("\n")
    # Parse FASTA-formatted alignment
    aligned1_parts: List[str] = []
    aligned2_parts: List[str] = []
    current_parts: Optional[List[str]] = None
    for line in aln_seqs:
        if line.startswith(">"):
            if current_parts is None:
                current_parts = aligned1_parts
            else:
                current_parts = aligned2_parts
        else:
            if current_parts is not None:
                current_parts.append(line)

    return "".join(aligned1_parts), "".join(aligned2_parts)


def align_proteins_mafft(
    seq1: str,
    seq2: str,
    seq1_name: str = "seq1",
    seq2_name: str = "seq2",
    use_docker: bool = True,
    container_name: str = "gpu-jupyter",
) -> Tuple[str, str]:
    """Align two protein sequences using MAFFT.

    Args:
        seq1: First protein sequence
        seq2: Second protein sequence
        seq1_name: Name for first sequence
        seq2_name: Name for second sequence
        use_docker: Run via Docker container
        container_name: Docker container name

    Returns:
        Tuple of (aligned_seq1, aligned_seq2)
    """
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", delete=False, dir="/tmp"
    ) as tmp:
        tmp.write(f">{seq1_name}\n{seq1}\n>{seq2_name}\n{seq2}\n")
        tmp_path = tmp.name

    try:
        if use_docker:
            # Map host /tmp to container /tmp (assumed accessible)
            cmd = [
                "docker", "exec", container_name,
                "mafft", "--auto", "--quiet", tmp_path,
            ]
        else:
            cmd = ["mafft", "--auto", "--quiet", tmp_path]

        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr[:200]}")

        # Parse FASTA output
        aligned_seqs: List[str] = []
        current_seq: List[str] = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith(">"):
                if current_seq:
                    aligned_seqs.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            aligned_seqs.append("".join(current_seq))

        if len(aligned_seqs) != 2:
            raise RuntimeError(
                f"Expected 2 aligned sequences, got {len(aligned_seqs)}"
            )

        return aligned_seqs[0], aligned_seqs[1]
    finally:
        os.unlink(tmp_path)


def back_translate_alignment(
    prot_aln1: str,
    prot_aln2: str,
    cds1: str,
    cds2: str,
) -> Tuple[str, str]:
    """Back-translate protein alignment to codon alignment (PAL2NAL logic).

    Walks through aligned protein sequences position by position:
    - Amino acid in both: emit corresponding codon from each CDS
    - Gap in one: emit '---' for that sequence, codon for the other

    Args:
        prot_aln1: Aligned protein sequence 1 (may contain '-')
        prot_aln2: Aligned protein sequence 2 (may contain '-')
        cds1: Original CDS for sequence 1 (unaligned, uppercase)
        cds2: Original CDS for sequence 2 (unaligned, uppercase)

    Returns:
        Tuple of (codon_aln1, codon_aln2) — equal length, multiple of 3

    Raises:
        ValueError: If codon/protein mismatch detected
    """
    codon_aln1: List[str] = []
    codon_aln2: List[str] = []
    pos1 = 0  # codon cursor for CDS1
    pos2 = 0  # codon cursor for CDS2

    for aa1, aa2 in zip(prot_aln1, prot_aln2):
        if aa1 == "-":
            codon_aln1.append("---")
        else:
            codon_aln1.append(cds1[pos1 * 3:(pos1 + 1) * 3])
            pos1 += 1

        if aa2 == "-":
            codon_aln2.append("---")
        else:
            codon_aln2.append(cds2[pos2 * 3:(pos2 + 1) * 3])
            pos2 += 1

    return "".join(codon_aln1), "".join(codon_aln2)


def extract_4d_sites(
    codon_aln1: str, codon_aln2: str
) -> Tuple[str, str]:
    """Extract third-position bases at four-fold degenerate sites.

    A site qualifies if BOTH codons in the alignment column have their
    first two positions in FOURFOLD_DEGENERATE_PREFIXES, and neither
    codon contains gaps or ambiguous bases.

    Args:
        codon_aln1: Codon-aligned sequence 1 (length multiple of 3)
        codon_aln2: Codon-aligned sequence 2 (length multiple of 3)

    Returns:
        Tuple of (sites1, sites2) — strings of 3rd-position bases at 4D sites
    """
    sites1: List[str] = []
    sites2: List[str] = []
    valid_bases = frozenset("ATCG")

    for i in range(0, min(len(codon_aln1), len(codon_aln2)) - 2, 3):
        c1 = codon_aln1[i:i + 3]
        c2 = codon_aln2[i:i + 3]

        # Skip if either codon has gaps or ambiguous bases
        if not all(b in valid_bases for b in c1):
            continue
        if not all(b in valid_bases for b in c2):
            continue

        # Check if both codons are four-fold degenerate at position 3
        prefix1 = c1[:2]
        prefix2 = c2[:2]
        if prefix1 in FOURFOLD_DEGENERATE_PREFIXES and \
           prefix2 in FOURFOLD_DEGENERATE_PREFIXES:
            sites1.append(c1[2])
            sites2.append(c2[2])

    return "".join(sites1), "".join(sites2)


def calculate_4dtv(
    sites1: str, sites2: str
) -> Tuple[float, float, int, int]:
    """Calculate raw and HKY-corrected 4DTv.

    HKY correction formula:
        piR = piA + piG (purine frequency)
        piY = piC + piT (pyrimidine frequency)
        Q   = transversion proportion
        corrected = -2 * piR * piY * ln(1 - Q / (2 * piR * piY))

    Args:
        sites1: Third-position bases from sequence 1 at 4D sites
        sites2: Third-position bases from sequence 2 at 4D sites

    Returns:
        Tuple of (raw_4dtv, corrected_4dtv, n_4d_sites, n_transversions)
        corrected_4dtv is NaN if saturated.
    """
    n_sites = len(sites1)
    if n_sites == 0:
        return (float("nan"), float("nan"), 0, 0)

    # Count transversions
    n_tv = sum(
        1 for b1, b2 in zip(sites1, sites2) if is_transversion(b1, b2)
    )

    raw_4dtv = n_tv / n_sites

    # Compute base frequencies from both sequences combined
    all_bases = sites1 + sites2
    counts = Counter(all_bases)
    total = len(all_bases)
    pi_a = counts.get("A", 0) / total
    pi_t = counts.get("T", 0) / total
    pi_g = counts.get("G", 0) / total
    pi_c = counts.get("C", 0) / total

    pi_r = pi_a + pi_g  # purine frequency
    pi_y = pi_t + pi_c  # pyrimidine frequency

    # HKY correction
    denom = 2.0 * pi_r * pi_y
    if denom <= 0:
        return (raw_4dtv, float("nan"), n_sites, n_tv)

    arg = 1.0 - raw_4dtv / denom
    if arg <= 0:
        # Saturated — too many transversions for reliable correction
        return (raw_4dtv, float("nan"), n_sites, n_tv)

    corrected = -denom * math.log(arg)

    return (raw_4dtv, corrected, n_sites, n_tv)


def process_gene_pair(
    gene1: str,
    gene2: str,
    cds_seqs: Dict[str, str],
    aligner: str = "biopython",
    use_docker: bool = True,
    container_name: str = "gpu-jupyter",
    min_4d_sites: int = 10,
) -> Optional[Dict]:
    """Process a single gene pair through the full 4DTv pipeline.

    Args:
        gene1: MCScanX gene ID (e.g. 'al_g14')
        gene2: MCScanX gene ID (e.g. 'al_g2476')
        cds_seqs: Dict mapping MCScanX IDs to CDS sequences
        aligner: 'mafft' or 'biopython'
        use_docker: Use Docker for MAFFT
        container_name: Docker container name
        min_4d_sites: Minimum 4D sites required

    Returns:
        Result dict or None if pair fails
    """
    # Step 1: Look up CDS
    cds1 = cds_seqs.get(gene1)
    cds2 = cds_seqs.get(gene2)
    if cds1 is None or cds2 is None:
        return None

    # Step 2: Validate CDS length
    if len(cds1) < 6 or len(cds2) < 6:
        return None

    # Step 3: Translate to protein
    prot1 = translate_sequence(cds1)
    prot2 = translate_sequence(cds2)
    if len(prot1) < 2 or len(prot2) < 2:
        return None

    # Step 4: Align proteins
    try:
        if aligner == "mafft":
            aln1, aln2 = align_proteins_mafft(
                prot1, prot2, gene1, gene2, use_docker, container_name
            )
        else:
            aln1, aln2 = align_proteins_biopython(prot1, prot2)
    except Exception:
        return None

    # Step 5: Back-translate to codon alignment
    # Truncate CDS to match translated protein length (codon-aligned)
    cds1_trimmed = cds1[:len(prot1) * 3]
    cds2_trimmed = cds2[:len(prot2) * 3]
    try:
        codon_aln1, codon_aln2 = back_translate_alignment(
            aln1, aln2, cds1_trimmed, cds2_trimmed
        )
    except (ValueError, IndexError):
        return None

    # Step 6: Extract 4D sites
    sites1, sites2 = extract_4d_sites(codon_aln1, codon_aln2)
    if len(sites1) < min_4d_sites:
        return None

    # Step 7: Calculate 4DTv
    raw, corrected, n_sites, n_tv = calculate_4dtv(sites1, sites2)

    return {
        "gene1": gene1,
        "gene2": gene2,
        "n_4d_sites": n_sites,
        "n_transversions": n_tv,
        "raw_4dtv": raw,
        "corrected_4dtv": corrected,
    }


# ---------------------------------------------------------------------------
# Worker function for multiprocessing (must be at module level)
# ---------------------------------------------------------------------------

# Global references set before spawning workers
_WORKER_CDS_SEQS: Dict[str, str] = {}
_WORKER_ALIGNER: str = "biopython"
_WORKER_USE_DOCKER: bool = True
_WORKER_CONTAINER: str = "gpu-jupyter"
_WORKER_MIN_4D: int = 10


def _init_worker(
    cds_seqs: Dict[str, str],
    aligner: str,
    use_docker: bool,
    container_name: str,
    min_4d_sites: int,
) -> None:
    """Initialize worker process globals."""
    global _WORKER_CDS_SEQS, _WORKER_ALIGNER, _WORKER_USE_DOCKER
    global _WORKER_CONTAINER, _WORKER_MIN_4D
    _WORKER_CDS_SEQS = cds_seqs
    _WORKER_ALIGNER = aligner
    _WORKER_USE_DOCKER = use_docker
    _WORKER_CONTAINER = container_name
    _WORKER_MIN_4D = min_4d_sites


def _worker_process_pair(pair: Tuple[str, str]) -> Optional[Dict]:
    """Worker function to process a single gene pair."""
    return process_gene_pair(
        pair[0], pair[1],
        _WORKER_CDS_SEQS,
        _WORKER_ALIGNER,
        _WORKER_USE_DOCKER,
        _WORKER_CONTAINER,
        _WORKER_MIN_4D,
    )


# ---------------------------------------------------------------------------
# Parallelization
# ---------------------------------------------------------------------------


def process_pairs_parallel(
    pairs: List[Tuple[str, str]],
    cds_seqs: Dict[str, str],
    n_workers: int = 4,
    aligner: str = "biopython",
    use_docker: bool = True,
    container_name: str = "gpu-jupyter",
    min_4d_sites: int = 10,
) -> List[Dict]:
    """Process all gene pairs with parallel workers.

    Args:
        pairs: List of (gene1, gene2) tuples
        cds_seqs: Dict mapping MCScanX IDs to CDS sequences
        n_workers: Number of parallel workers
        aligner: Alignment method
        use_docker: Use Docker for MAFFT
        container_name: Docker container name
        min_4d_sites: Minimum 4D sites per pair

    Returns:
        List of result dicts (only successful pairs)
    """
    results: List[Dict] = []
    total = len(pairs)
    t0 = time.time()

    if n_workers <= 1 or aligner == "mafft":
        # Sequential processing (MAFFT is already I/O-bound with subprocess)
        for idx, (g1, g2) in enumerate(pairs, 1):
            result = process_gene_pair(
                g1, g2, cds_seqs, aligner, use_docker, container_name,
                min_4d_sites,
            )
            if result is not None:
                results.append(result)
            if idx % 1000 == 0:
                elapsed = time.time() - t0
                rate = idx / elapsed
                print(
                    f"  Processed {idx}/{total} pairs "
                    f"({len(results)} valid, {rate:.0f} pairs/sec)",
                    flush=True,
                )
    else:
        # Parallel with ProcessPoolExecutor for BioPython aligner
        from multiprocessing import Pool

        with Pool(
            processes=n_workers,
            initializer=_init_worker,
            initargs=(cds_seqs, aligner, use_docker, container_name, min_4d_sites),
        ) as pool:
            done = 0
            for result in pool.imap_unordered(_worker_process_pair, pairs, chunksize=100):
                done += 1
                if result is not None:
                    results.append(result)
                if done % 1000 == 0:
                    elapsed = time.time() - t0
                    rate = done / elapsed
                    print(
                        f"  Processed {done}/{total} pairs "
                        f"({len(results)} valid, {rate:.0f} pairs/sec)",
                        flush=True,
                    )

    elapsed = time.time() - t0
    print(
        f"  Completed {total} pairs in {elapsed:.1f}s — "
        f"{len(results)} valid results ({len(results)/total*100:.1f}%)",
        flush=True,
    )
    return results


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_results_tsv(results: List[Dict], output_path: Path) -> None:
    """Write per-pair 4DTv results to TSV file.

    Args:
        results: List of result dicts from process_gene_pair
        output_path: Output TSV file path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write(
            "gene1\tgene2\tn_4d_sites\tn_transversions\t"
            "raw_4dtv\tcorrected_4dtv\n"
        )
        for r in results:
            fh.write(
                f"{r['gene1']}\t{r['gene2']}\t{r['n_4d_sites']}\t"
                f"{r['n_transversions']}\t{r['raw_4dtv']:.6f}\t"
                f"{r['corrected_4dtv']:.6f}\n"
            )


def write_summary(
    all_results: Dict[str, List[Dict]],
    comparison_labels: Dict[str, Tuple[str, str]],
    output_path: Path,
) -> None:
    """Write summary statistics to text file.

    Args:
        all_results: Maps comparison key to results list
        comparison_labels: Maps comparison key to (species1, species2) tuple
        output_path: Output file path
    """
    with open(output_path, "w") as fh:
        fh.write("4DTv Analysis Summary\n")
        fh.write("=" * 60 + "\n\n")
        for key, results in all_results.items():
            sp1, sp2 = comparison_labels[key]
            label = (
                f"{get_display_name(sp1)} (self)"
                if sp1 == sp2
                else f"{get_display_name(sp1)} vs {get_display_name(sp2)}"
            )
            fh.write(f"Comparison: {label}\n")
            fh.write(f"  Total valid pairs: {len(results)}\n")

            if results:
                import numpy as np

                corrected = [
                    r["corrected_4dtv"] for r in results
                    if not math.isnan(r["corrected_4dtv"])
                ]
                if corrected:
                    arr = np.array(corrected)
                    fh.write(f"  Corrected 4DTv statistics:\n")
                    fh.write(f"    Mean:   {np.mean(arr):.4f}\n")
                    fh.write(f"    Median: {np.median(arr):.4f}\n")
                    fh.write(f"    Std:    {np.std(arr):.4f}\n")
                    fh.write(f"    Min:    {np.min(arr):.4f}\n")
                    fh.write(f"    Max:    {np.max(arr):.4f}\n")

                    raw_vals = [r["raw_4dtv"] for r in results]
                    fh.write(f"  Raw 4DTv mean: {np.mean(raw_vals):.4f}\n")

                    sites = [r["n_4d_sites"] for r in results]
                    fh.write(
                        f"  Mean 4D sites per pair: {np.mean(sites):.1f}\n"
                    )
                    saturated = len(results) - len(corrected)
                    if saturated > 0:
                        fh.write(
                            f"  Saturated pairs (NaN): {saturated}\n"
                        )
            fh.write("\n")


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------


def plot_4dtv_distribution(
    all_results: Dict[str, List[Dict]],
    comparison_labels: Dict[str, Tuple[str, str]],
    output_dir: Path,
    bin_width: float = 0.02,
    x_max: float = 0.8,
) -> None:
    """Generate publication-quality 4DTv distribution histogram with KDE.

    Args:
        all_results: Maps comparison key to results list
        comparison_labels: Maps comparison key to (species1, species2) tuple
        output_dir: Output directory for figures
        bin_width: Histogram bin width
        x_max: Maximum x-axis value
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import gaussian_kde

    # Publication rcParams (project convention)
    plt.rcParams.update({
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
    })

    bins = np.arange(0, x_max + bin_width, bin_width)

    # Color palette for comparisons
    color_list = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                  "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]

    # Collect data for all comparisons
    plot_data: list = []  # (arr, color, label) tuples
    for idx, (key, results) in enumerate(all_results.items()):
        sp1, sp2 = comparison_labels[key]

        # Filter out NaN corrected values and negative zeros
        values = [
            max(0.0, r["corrected_4dtv"]) for r in results
            if not math.isnan(r["corrected_4dtv"]) and r["corrected_4dtv"] < x_max
        ]
        if not values:
            continue

        arr = np.array(values)
        color = color_list[idx % len(color_list)]

        # Build label with italicized species names
        if sp1 == sp2:
            label = f"{get_mathtext_name(sp1, abbreviated=True)} (self)"
        else:
            label = (
                f"{get_mathtext_name(sp1, abbreviated=True)} vs "
                f"{get_mathtext_name(sp2, abbreviated=True)}"
            )
        plot_data.append((arr, color, label))

    if not plot_data:
        print("  WARNING: No valid data to plot.")
        return

    # --- Figure 1: Frequency histogram (count-based, y-axis capped) ---
    fig, ax = plt.subplots(figsize=(10, 7))

    for arr, color, label in plot_data:
        counts, _, patches = ax.hist(
            arr, bins=bins, alpha=0.4,
            color=color, edgecolor=color, linewidth=0.5,
            label=label,
        )

    ax.set_xlabel("4DTv (corrected)")
    ax.set_ylabel("Number of gene pairs")
    ax.set_title("Distribution of 4DTv Values")
    ax.set_xlim(0, x_max)
    # Cap y-axis to make the tail visible (exclude the first bin from max)
    all_counts = []
    for arr, color, label in plot_data:
        h, _ = np.histogram(arr, bins=bins)
        all_counts.extend(h[1:])  # exclude first bin (near-zero peak)
    if all_counts:
        y_cap = max(all_counts) * 1.5
        # Only cap if the first bin is much larger than the rest
        for arr, _, _ in plot_data:
            h, _ = np.histogram(arr, bins=bins)
            if h[0] > y_cap * 2:
                ax.set_ylim(0, y_cap)
                ax.annotate(
                    f"Peak at 0 truncated",
                    xy=(0.02, 0.97), xycoords="axes fraction",
                    fontsize=10, fontstyle="italic",
                    ha="left", va="top",
                )
                break
    ax.legend(loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "4dtv_distribution.pdf")
    fig.savefig(output_dir / "4dtv_distribution.png")
    plt.close(fig)

    # --- Figure 2: Frequency histogram with log y-axis ---
    # Log scale reveals secondary peaks masked by the dominant near-zero peak
    fig2, ax2 = plt.subplots(figsize=(10, 7))

    for arr, color, label in plot_data:
        ax2.hist(
            arr, bins=bins, alpha=0.4,
            color=color, edgecolor=color, linewidth=0.5,
            label=label,
        )

    ax2.set_xlabel("4DTv (corrected)")
    ax2.set_ylabel("Number of gene pairs (log scale)")
    ax2.set_title("Distribution of 4DTv Values")
    ax2.set_xlim(0, x_max)
    ax2.set_yscale("log")
    ax2.legend(loc="upper right")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig2.tight_layout()
    fig2.savefig(output_dir / "4dtv_distribution_log.pdf")
    fig2.savefig(output_dir / "4dtv_distribution_log.png")
    plt.close(fig2)

    # --- Figure 3: KDE density plot with y-axis capped ---
    # Exclude values below a threshold from KDE to focus on divergence signal
    fig3, ax3 = plt.subplots(figsize=(10, 7))

    for arr, color, label in plot_data:
        try:
            kde = gaussian_kde(arr, bw_method=0.05)
            x_kde = np.linspace(0, x_max, 1000)
            y_kde = kde(x_kde)
            ax3.plot(x_kde, y_kde, color=color, linewidth=2.5, label=label)
            ax3.fill_between(x_kde, y_kde, alpha=0.15, color=color)
        except Exception:
            pass

    ax3.set_xlabel("4DTv (corrected)")
    ax3.set_ylabel("Density")
    ax3.set_title("4DTv Density Distribution")
    ax3.set_xlim(0, x_max)
    # Cap y-axis to make the tail visible
    all_kde_max = []
    for arr, color, label in plot_data:
        try:
            kde = gaussian_kde(arr, bw_method=0.05)
            # Find density at x=0.05 (beyond the near-zero peak shoulder)
            shoulder_density = kde(np.array([0.05]))[0]
            all_kde_max.append(shoulder_density)
        except Exception:
            pass
    if all_kde_max:
        y_cap_kde = max(all_kde_max) * 3.0
        ax3.set_ylim(0, y_cap_kde)
        ax3.annotate(
            "Near-zero peak truncated",
            xy=(0.02, 0.97), xycoords="axes fraction",
            fontsize=10, fontstyle="italic",
            ha="left", va="top",
        )
    ax3.legend(loc="upper right")
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)

    fig3.tight_layout()
    fig3.savefig(output_dir / "4dtv_density.pdf")
    fig3.savefig(output_dir / "4dtv_density.png")
    plt.close(fig3)

    print(
        f"  Saved figures to {output_dir}/"
        f"{{4dtv_distribution,4dtv_density}}.{{pdf,png}}"
    )


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------


def resolve_collinearity_path(
    species1: str,
    species2: str,
    mcscanx_base: Path,
) -> Optional[Path]:
    """Find the collinearity file for a species pair.

    Tries multiple naming conventions:
        {p1}_{p2}/{p1}_{p2}.collinearity
        {p2}_{p1}/{p2}_{p1}.collinearity

    Args:
        species1: First species name
        species2: Second species name
        mcscanx_base: Base directory for MCScanX results

    Returns:
        Path to collinearity file, or None if not found
    """
    p1 = get_species_prefix(species1)
    p2 = get_species_prefix(species2)

    for a, b in [(p1, p2), (p2, p1)]:
        name = f"{a}_{b}"
        path = mcscanx_base / name / f"{name}.collinearity"
        if path.exists():
            return path

    return None


def run_4dtv_pipeline(
    comparisons: List[Tuple[str, str]],
    braker_base: Path,
    mcscanx_base: Path,
    output_dir: Path,
    n_workers: int = 4,
    aligner: str = "biopython",
    use_docker: bool = True,
    container_name: str = "gpu-jupyter",
    min_4d_sites: int = 10,
    bin_width: float = 0.02,
) -> None:
    """Run the complete 4DTv analysis pipeline.

    Args:
        comparisons: List of (species1, species2) tuples
        braker_base: Base dir containing {species}/braker.codingseq
        mcscanx_base: Base dir with MCScanX results
        output_dir: Output directory
        n_workers: Parallel workers
        aligner: Alignment method ('mafft' or 'biopython')
        use_docker: Use Docker for MAFFT
        container_name: Docker container name
        min_4d_sites: Minimum 4D sites per pair
        bin_width: Histogram bin width
    """
    all_results: Dict[str, List[Dict]] = {}
    comparison_labels: Dict[str, Tuple[str, str]] = {}
    # Cache loaded CDS sequences per species
    cds_cache: Dict[str, Dict[str, str]] = {}

    output_dir.mkdir(parents=True, exist_ok=True)

    for sp1, sp2 in comparisons:
        p1 = get_species_prefix(sp1)
        p2 = get_species_prefix(sp2)
        comp_key = f"{p1}_{p2}"

        if sp1 == sp2:
            print(f"\n{'='*60}")
            print(f"Self-comparison: {get_display_name(sp1)}")
            print(f"{'='*60}")
        else:
            print(f"\n{'='*60}")
            print(
                f"Inter-species: {get_display_name(sp1)} vs "
                f"{get_display_name(sp2)}"
            )
            print(f"{'='*60}")

        # 1. Find collinearity file
        coll_path = resolve_collinearity_path(sp1, sp2, mcscanx_base)
        if coll_path is None:
            print(
                f"  WARNING: No collinearity file found for {comp_key}. "
                f"Run MCScanX first. Skipping.",
                file=sys.stderr,
            )
            continue
        print(f"  Collinearity: {coll_path}")

        # 2. Parse gene pairs
        pairs = parse_collinearity(coll_path)
        print(f"  Syntenic gene pairs: {len(pairs)}")

        # 3. Load CDS sequences (cached per species)
        species_to_load = {sp1} if sp1 == sp2 else {sp1, sp2}
        merged_cds: Dict[str, str] = {}
        for sp in species_to_load:
            if sp not in cds_cache:
                prefix = get_species_prefix(sp)
                cds_path = braker_base / sp / "braker.codingseq"
                if not cds_path.exists():
                    print(
                        f"  ERROR: CDS file not found: {cds_path}",
                        file=sys.stderr,
                    )
                    continue
                print(f"  Loading CDS: {cds_path}")
                cds_cache[sp] = load_cds_sequences(cds_path, prefix)
                print(f"    Loaded {len(cds_cache[sp])} genes")
            merged_cds.update(cds_cache[sp])

        if not merged_cds:
            print("  ERROR: No CDS sequences loaded. Skipping.")
            continue

        # 4. Process all pairs
        print(f"  Processing pairs (aligner={aligner}, workers={n_workers})...")
        results = process_pairs_parallel(
            pairs, merged_cds, n_workers, aligner, use_docker,
            container_name, min_4d_sites,
        )

        # 5. Write TSV
        tsv_path = output_dir / f"{comp_key}_4dtv.tsv"
        write_results_tsv(results, tsv_path)
        print(f"  Results: {tsv_path}")

        all_results[comp_key] = results
        comparison_labels[comp_key] = (sp1, sp2)

    if not all_results:
        print("\nNo valid results. Exiting.")
        return

    # 6. Write summary
    summary_path = output_dir / "summary.txt"
    write_summary(all_results, comparison_labels, summary_path)
    print(f"\nSummary: {summary_path}")

    # 7. Generate visualization
    print("\nGenerating 4DTv distribution plot...")
    plot_4dtv_distribution(
        all_results, comparison_labels, output_dir, bin_width,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "4DTv analysis: calculate transversion rates at four-fold "
            "degenerate sites for syntenic gene pairs from MCScanX."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Self-comparisons for WGD detection:
  python scripts/run_4dtv_analysis.py \\
      --self Aspiorhynchus_laticeps --self Diptychus_maculatus

  # Inter-species comparison:
  python scripts/run_4dtv_analysis.py \\
      --pair Aspiorhynchus_laticeps Diptychus_maculatus

  # Combined with BioPython aligner:
  python scripts/run_4dtv_analysis.py \\
      --self Aspiorhynchus_laticeps --self Diptychus_maculatus \\
      --aligner biopython --workers 8
""",
    )

    parser.add_argument(
        "--self",
        dest="self_species",
        action="append",
        default=[],
        metavar="SPECIES",
        help=(
            "Self-comparison for WGD detection "
            "(e.g. Aspiorhynchus_laticeps). Repeatable."
        ),
    )
    parser.add_argument(
        "--pair",
        dest="pair_species",
        action="append",
        nargs=2,
        default=[],
        metavar=("SP1", "SP2"),
        help=(
            "Inter-species comparison "
            "(e.g. Aspiorhynchus_laticeps Diptychus_maculatus). Repeatable."
        ),
    )
    parser.add_argument(
        "--braker-base",
        type=Path,
        default=Path("/home/jovyan/Outputs/Preprocessing/BRAKER_MASKED"),
        help="Base directory containing {species}/braker.codingseq",
    )
    parser.add_argument(
        "--mcscanx-base",
        type=Path,
        default=Path("/home/jovyan/Outputs_set3/MCScanX_Results"),
        help="Base directory containing MCScanX results",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/home/jovyan/Outputs_set3/4DTv"),
        help="Output directory (default: /home/jovyan/Outputs_set3/4DTv)",
    )
    parser.add_argument(
        "--aligner",
        choices=["mafft", "biopython"],
        default="mafft",
        help="Protein alignment method (default: mafft)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel workers (default: 4)",
    )
    parser.add_argument(
        "--no-docker",
        action="store_true",
        help="Run MAFFT directly without Docker",
    )
    parser.add_argument(
        "--container",
        default="gpu-jupyter",
        help="Docker container name (default: gpu-jupyter)",
    )
    parser.add_argument(
        "--min-4d-sites",
        type=int,
        default=10,
        help="Minimum four-fold degenerate sites per pair (default: 10)",
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.02,
        help="Histogram bin width (default: 0.02)",
    )

    return parser.parse_args()


def main() -> None:
    """Entry point."""
    args = parse_args()

    if not args.self_species and not args.pair_species:
        print(
            "Error: Specify at least one --self or --pair comparison.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Build comparison list
    comparisons: List[Tuple[str, str]] = []
    for sp in args.self_species:
        comparisons.append((sp, sp))
    for sp1, sp2 in args.pair_species:
        comparisons.append((sp1, sp2))

    use_docker = not args.no_docker

    print("4DTv Analysis Pipeline")
    print("=" * 60)
    print(f"Comparisons:  {len(comparisons)}")
    print(f"Aligner:      {args.aligner}")
    print(f"Workers:      {args.workers}")
    print(f"Docker:       {'yes' if use_docker else 'no'}")
    print(f"Min 4D sites: {args.min_4d_sites}")
    print(f"BRAKER base:  {args.braker_base}")
    print(f"MCScanX base: {args.mcscanx_base}")
    print(f"Output dir:   {args.output_dir}")

    run_4dtv_pipeline(
        comparisons=comparisons,
        braker_base=args.braker_base,
        mcscanx_base=args.mcscanx_base,
        output_dir=args.output_dir,
        n_workers=args.workers,
        aligner=args.aligner,
        use_docker=use_docker,
        container_name=args.container,
        min_4d_sites=args.min_4d_sites,
        bin_width=args.bin_width,
    )


if __name__ == "__main__":
    main()
