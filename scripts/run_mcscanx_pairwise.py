#!/usr/bin/env python3
"""MCScanX Pairwise Synteny Analysis Pipeline.

Complete pipeline for pairwise synteny analysis between any two species using
MCScanX. Extracted from 6-ultimate-MCScanX.qmd notebook into a standalone,
reusable script.

Pipeline steps:
    1. Create combined protein file with species-prefixed headers
    2. Build BLAST database
    3. Run all-vs-all BLASTP search
    4. Convert GFF3 annotations to BED format
    5. Process BED files (standardize chromosome names, filter scaffolds)
    6. Create MCScanX input files (.gff and .blast)
    7. Run MCScanX synteny detection
    8. Generate visualizations (dot plot, dual synteny, circle, bar)
    9. Report results summary

Input data requirements (per species):
    - braker.aa   : Protein sequences from BRAKER3
    - braker.gff3 : Gene annotations from BRAKER3

Usage:
    # Basic usage (auto-discovers files under --braker-base):
    python scripts/run_mcscanx_pairwise.py \\
        Aspiorhynchus_laticeps Schizothorax_macropogon

    # With explicit paths:
    python scripts/run_mcscanx_pairwise.py \\
        Aspiorhynchus_laticeps Schizothorax_macropogon \\
        --braker-base /home/jovyan/Outputs/Preprocessing/BRAKER_UNMASKED \\
        --output-dir /home/jovyan/Outputs/MCScanX_Pairwise/al_sm

    # Skip BLAST (reuse existing results) and skip visualizations:
    python scripts/run_mcscanx_pairwise.py \\
        Aspiorhynchus_laticeps Schizothorax_macropogon \\
        --skip-blast --skip-viz

Author: Extracted from 6-ultimate-MCScanX.qmd
"""

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


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


def extract_chromosome_number(chrom_name: str) -> str:
    """Extract chromosome/scaffold number from NCBI long names.

    Handles various naming conventions:
        - chromosome_1, chromosome 1
        - chr1
        - scaffold_1
        - CM041793.1_...chromosome_1...
        - Last number in string (fallback)

    Args:
        chrom_name: Raw chromosome name from GFF3/BED file

    Returns:
        Chromosome number as string
    """
    patterns = [
        r"chromosome[_\s]+(\d+)",
        r"chr(\d+)",
        r"scaffold[_\s]+(\d+)",
        r"CM\d+\.\d+.*chromosome[_\s]+(\d+)",
        r"(\d+)$",
    ]
    for pattern in patterns:
        match = re.search(pattern, chrom_name, re.IGNORECASE)
        if match:
            return match.group(1)
    return str(abs(hash(chrom_name)) % 1000)


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------


def create_pairwise_dataset(
    species1_aa: Path,
    species2_aa: Path,
    prefix1: str,
    prefix2: str,
    output_faa: Path,
) -> None:
    """Create combined protein dataset with species-prefixed headers.

    BRAKER3 outputs use sequential numbering (g1, g2, ...) per species.
    Without prefixes, both species would have g1, g2, etc., causing MCScanX
    to fail with 0 matches imported.

    Header format: >g1.t1 -> >al_g1.t1

    Args:
        species1_aa: Path to first species braker.aa
        species2_aa: Path to second species braker.aa
        prefix1: Species prefix for first species
        prefix2: Species prefix for second species
        output_faa: Path to write combined FASTA
    """
    print("Creating combined protein dataset...")
    output_faa.parent.mkdir(parents=True, exist_ok=True)

    file_list = [
        (species1_aa, prefix1),
        (species2_aa, prefix2),
    ]

    with open(output_faa, "w") as outfile:
        for input_file, species_prefix in file_list:
            if not input_file.exists():
                print(f"  WARNING: {input_file} not found", file=sys.stderr)
                continue

            seq_count = 0
            with open(input_file, "r") as infile:
                for line in infile:
                    if line.startswith(">"):
                        header = line[1:].strip()
                        outfile.write(f">{species_prefix}_{header}\n")
                        seq_count += 1
                    else:
                        outfile.write(line)

            size_mb = input_file.stat().st_size / (1024 * 1024)
            print(
                f"  Added {species_prefix.upper()}: {seq_count:,} sequences ({size_mb:.1f} MB)"
            )
            print(f"    Header format: >g1.t1 -> >{species_prefix}_g1.t1")

    final_size = output_faa.stat().st_size / (1024 * 1024)
    print(f"  Combined output: {output_faa} ({final_size:.1f} MB)")


def create_blast_database(
    faa_file: Path,
    output_prefix: Path,
    makeblastdb_cmd: str,
) -> bool:
    """Build BLAST protein database.

    Note: -parse_seqids is NOT used because BRAKER3 outputs use simple
    headers (>g1.t1) that don't follow NCBI's SeqID format requirements.

    Args:
        faa_file: Input combined FASTA file
        output_prefix: Output database prefix
        makeblastdb_cmd: Path to makeblastdb binary

    Returns:
        True if database creation succeeded
    """
    print(f"Building BLAST database from {faa_file.name}...")
    log_file = output_prefix.with_suffix(".log")

    cmd = (
        f"{makeblastdb_cmd} -in {faa_file} -dbtype prot "
        f"-out {output_prefix} -logfile {log_file} "
        f"-title {output_prefix.stem}"
    )
    print(f"  Command: {cmd}")

    exit_code = os.system(cmd)

    if exit_code == 0:
        db_files = [output_prefix.with_suffix(ext) for ext in [".phr", ".pin", ".psq"]]
        existing = [f for f in db_files if f.exists()]
        print(f"  ✓ Database created ({len(existing)}/3 files)")
        return True
    else:
        print(f"  ✗ Database creation failed (exit code {exit_code})", file=sys.stderr)
        return False


def run_blast_search(
    query_faa: Path,
    db_prefix: Path,
    output_file: Path,
    blastp_cmd: str,
    threads: int = 8,
    evalue: float = 1e-10,
) -> bool:
    """Run all-vs-all BLASTP search.

    Args:
        query_faa: Query FASTA file
        db_prefix: BLAST database prefix
        output_file: Output file path (tabular format 6)
        blastp_cmd: Path to blastp binary
        threads: Number of threads
        evalue: E-value threshold

    Returns:
        True if BLAST search succeeded
    """
    print("Running BLASTP search...")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"{blastp_cmd} -query {query_faa} -db {db_prefix} "
        f"-out {output_file} -evalue {evalue} -num_threads {threads} "
        f"-outfmt 6 -num_alignments 5"
    )
    print(f"  Parameters: E-value {evalue}, {threads} threads, format 6")
    print(f"  Output: {output_file}")

    start = time.time()
    exit_code = os.system(cmd)
    elapsed = time.time() - start

    if exit_code == 0 and output_file.exists():
        size_mb = output_file.stat().st_size / (1024 * 1024)
        with open(output_file, "r") as f:
            hit_count = sum(1 for _ in f)
        print(
            f"  ✓ BLAST completed in {elapsed:.1f}s ({hit_count:,} hits, {size_mb:.1f} MB)"
        )
        return True
    else:
        print(f"  ✗ BLAST search failed (exit code {exit_code})", file=sys.stderr)
        return False


def convert_gff3_to_bed(gff3_path: Path, bed_path: Path) -> bool:
    """Convert GFF3 annotation to BED format for MCScanX.

    Extracts gene features, converts 1-based GFF3 coordinates to 0-based BED.

    Args:
        gff3_path: Input GFF3 file
        bed_path: Output BED file

    Returns:
        True if conversion succeeded
    """
    print(f"Converting {gff3_path.name} to BED format...")

    if not gff3_path.exists():
        print(f"  ERROR: Input file not found: {gff3_path}", file=sys.stderr)
        return False

    gene_count = 0
    skipped = 0

    try:
        with open(gff3_path, "r") as infile, open(bed_path, "w") as outfile:
            for line in infile:
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    skipped += 1
                    continue

                chrom = parts[0]
                start = str(int(parts[3]) - 1)  # 1-based -> 0-based
                end = parts[4]
                attributes = parts[8]

                gene_id = None
                for field in attributes.split(";"):
                    if field.startswith("ID="):
                        gene_id = field[3:]
                        break

                if gene_id:
                    outfile.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")
                    gene_count += 1
                else:
                    skipped += 1

        print(f"  ✓ Converted {gene_count} genes (skipped {skipped} lines)")
        return True

    except Exception as e:
        print(f"  ERROR during conversion: {e}", file=sys.stderr)
        return False


def process_bed_file(
    input_bed: Path,
    output_bed: Path,
    species_prefix: str,
) -> int:
    """Standardize BED file for MCScanX.

    Processes BRAKER3-derived BED files:
    1. Extracts chromosome numbers from long NCBI names
    2. Adds species prefix to chromosome names (e.g. al1, sm2)
    3. Adds species prefix to gene IDs (e.g. al_g1)
    4. Filters out scaffolds/contigs (keeps chromosomes only)

    Args:
        input_bed: Input BED file from GFF3 conversion
        output_bed: Output processed BED file
        species_prefix: Species prefix (e.g. 'al', 'sm')

    Returns:
        Number of genes in the processed file
    """
    print(f"  Processing {species_prefix.upper()} BED file...")

    if not input_bed.exists():
        print(f"    ERROR: Input file not found: {input_bed}", file=sys.stderr)
        return 0

    with open(input_bed, "r") as infile:
        lines = infile.readlines()

    processed_lines = []
    chrom_mapping: Dict[str, str] = {}
    scaffold_count = 0

    for line in lines:
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue

        orig_chrom = parts[0]
        start = parts[1]
        end = parts[2]
        gene_id = parts[3]

        # Filter scaffolds/contigs — keep only chromosomes
        chrom_lower = orig_chrom.lower()
        if (
            "scaffold" in chrom_lower
            or "contig" in chrom_lower
            or "_ctg" in chrom_lower
            or "_ptg" in chrom_lower
        ):
            scaffold_count += 1
            continue

        chrom_num = extract_chromosome_number(orig_chrom)
        new_chrom = f"{species_prefix}{chrom_num}"
        new_gene_id = f"{species_prefix}_{gene_id}"

        if orig_chrom not in chrom_mapping:
            chrom_mapping[orig_chrom] = new_chrom

        processed_lines.append(f"{new_chrom}\t{start}\t{end}\t{new_gene_id}\n")

    with open(output_bed, "w") as outfile:
        outfile.writelines(processed_lines)

    print(f"    Processed {len(processed_lines)} genes")
    print(f"    Filtered {scaffold_count} scaffolds/contigs")
    print(f"    Unique chromosomes: {len(chrom_mapping)}")
    if chrom_mapping:
        for i, (orig, new) in enumerate(list(chrom_mapping.items())[:3]):
            orig_short = orig[:40] + "..." if len(orig) > 40 else orig
            print(f"      {orig_short} -> {new}")
        if len(chrom_mapping) > 3:
            print(f"      ... and {len(chrom_mapping) - 3} more")

    return len(processed_lines)


def create_mcscanx_inputs(
    comparison_dir: Path,
    comparison_name: str,
    prefix1: str,
    prefix2: str,
    blast_source: Path,
) -> bool:
    """Create final MCScanX input files (.gff and .blast).

    Combines processed BED files into MCScanX GFF format and strips
    transcript suffixes from BLAST hits to match GFF gene IDs.

    MCScanX GFF format: chromosome<TAB>gene_id<TAB>start<TAB>end

    Args:
        comparison_dir: Directory containing processed BED files
        comparison_name: Comparison identifier (e.g. 'al_sm')
        prefix1: First species prefix
        prefix2: Second species prefix
        blast_source: Path to raw BLAST output file

    Returns:
        True if all input files created successfully
    """
    print(f"Creating MCScanX input files for {comparison_name}...")

    combined_bed = comparison_dir / f"{comparison_name}.bed"
    combined_gff = comparison_dir / f"{comparison_name}.gff"
    target_blast = comparison_dir / f"{comparison_name}.blast"

    species_files = [
        comparison_dir / f"{prefix1}-filtered-modified.bed",
        comparison_dir / f"{prefix2}-filtered-modified.bed",
    ]

    # Combine BED files
    total_genes = 0
    with open(combined_bed, "w") as outfile:
        for species_file in species_files:
            if species_file.exists():
                with open(species_file, "r") as infile:
                    lines = infile.readlines()
                    outfile.writelines(lines)
                    total_genes += len(lines)
                    sp_name = species_file.stem.split("-")[0]
                    print(f"  Added {sp_name}: {len(lines)} genes")
            else:
                print(f"  WARNING: {species_file} not found", file=sys.stderr)

    # Reformat BED -> MCScanX GFF (chromosome, gene_id, start, end)
    with open(combined_bed, "r") as infile, open(combined_gff, "w") as outfile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                reordered = [parts[0], parts[3], parts[1], parts[2]]
                outfile.write("\t".join(reordered) + "\n")

    # Process BLAST: strip .t# transcript suffix to match GFF gene IDs
    if blast_source.exists():
        processed_count = 0
        with open(blast_source, "r") as infile, open(target_blast, "w") as outfile:
            for line in infile:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    parts[0] = re.sub(r"\.t\d+$", "", parts[0])
                    parts[1] = re.sub(r"\.t\d+$", "", parts[1])
                    outfile.write("\t".join(parts) + "\n")
                    processed_count += 1
        blast_size = target_blast.stat().st_size / (1024 * 1024)
        print(f"  Processed BLAST file: {processed_count:,} hits, {blast_size:.1f} MB")
        print(f"    (Stripped .t# transcript suffixes to match GFF gene IDs)")
    else:
        print(f"  WARNING: BLAST file not found: {blast_source}", file=sys.stderr)

    print(f"  Combined {total_genes} genes in {combined_gff}")

    # Verify
    required = [combined_gff, target_blast]
    all_present = all(f.exists() for f in required)

    if all_present:
        print(f"  ✓ Input files ready for MCScanX")
    else:
        print(f"  ✗ Missing required files", file=sys.stderr)

    return all_present


def run_mcscanx(input_prefix: Path, mcscanx_cmd: str) -> bool:
    """Execute MCScanX synteny detection.

    Args:
        input_prefix: Path prefix for .gff and .blast files
        mcscanx_cmd: Path to MCScanX binary

    Returns:
        True if MCScanX analysis succeeded
    """
    gff_file = Path(f"{input_prefix}.gff")
    blast_file = Path(f"{input_prefix}.blast")

    if not gff_file.exists() or not blast_file.exists():
        print(f"✗ Missing input files:", file=sys.stderr)
        print(f"  GFF:   {gff_file} (exists={gff_file.exists()})", file=sys.stderr)
        print(f"  BLAST: {blast_file} (exists={blast_file.exists()})", file=sys.stderr)
        return False

    with open(gff_file, "r") as f:
        gene_count = sum(1 for _ in f)
    blast_size = blast_file.stat().st_size / (1024 * 1024)

    print(f"Running MCScanX...")
    print(f"  Input genes: {gene_count}")
    print(f"  BLAST file: {blast_size:.1f} MB")

    cmd = f"{mcscanx_cmd} {input_prefix}"
    print(f"  Command: {cmd}")

    start = time.time()
    exit_code = os.system(cmd)
    elapsed = time.time() - start

    if exit_code != 0:
        print(f"  ✗ MCScanX failed (exit code {exit_code})", file=sys.stderr)
        return False

    print(f"  ✓ MCScanX completed in {elapsed:.1f}s")

    # Report results
    collinearity_file = Path(f"{input_prefix}.collinearity")
    tandem_file = Path(f"{input_prefix}.tandem")

    if collinearity_file.exists():
        with open(collinearity_file, "r") as f:
            lines = f.readlines()
        blocks = len([l for l in lines if l.startswith("##")])
        pairs = len([l for l in lines if not l.startswith("#")])
        print(f"  Collinearity blocks: {blocks}")
        print(f"  Syntenic gene pairs: {pairs}")

    if tandem_file.exists():
        with open(tandem_file, "r") as f:
            tandem_count = sum(1 for _ in f)
        print(f"  Tandem gene pairs: {tandem_count}")

    return True


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------


def compile_visualization_tools(tools_dir: Path) -> bool:
    """Compile Java-based MCScanX visualization tools.

    Args:
        tools_dir: Path to MCScanX downstream_analyses directory

    Returns:
        True if all tools compiled successfully
    """
    print("Compiling Java visualization tools...")

    java_files = [
        "dual_synteny_plotter.java",
        "dot_plotter.java",
        "circle_plotter.java",
        "bar_plotter.java",
        "family_tree_plotter.java",
        "Cubic.java",
    ]

    success_count = 0
    for java_file in java_files:
        file_path = tools_dir / java_file
        if file_path.exists():
            cmd = f"javac {java_file}"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, cwd=str(tools_dir)
            )
            if result.returncode == 0:
                print(f"  ✓ {java_file}")
                success_count += 1
            else:
                print(f"  ✗ {java_file} compilation failed")
        else:
            print(f"  ✗ {java_file} not found")

    print(f"  Compiled {success_count}/{len(java_files)} tools")
    return success_count == len(java_files)


def extract_chromosome_info(
    gff_file: Path,
    prefix1: str,
    prefix2: str,
) -> Dict[str, List[str]]:
    """Extract chromosome info from MCScanX GFF file.

    Args:
        gff_file: Path to combined .gff file
        prefix1: First species prefix
        prefix2: Second species prefix

    Returns:
        Dict with 'species1_chroms' and 'species2_chroms' lists
    """
    if not gff_file.exists():
        print(f"  GFF file not found: {gff_file}", file=sys.stderr)
        return {}

    chromosomes = set()
    with open(gff_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                chromosomes.add(parts[0])

    sp1_chroms = sorted(
        [c for c in chromosomes if c.startswith(prefix1)],
        key=lambda x: int(re.sub(r"[^0-9]", "", x) or "0"),
    )
    sp2_chroms = sorted(
        [c for c in chromosomes if c.startswith(prefix2)],
        key=lambda x: int(re.sub(r"[^0-9]", "", x) or "0"),
    )

    print(f"  {prefix1.upper()} chromosomes: {len(sp1_chroms)}")
    print(f"  {prefix2.upper()} chromosomes: {len(sp2_chroms)}")

    return {
        "species1_chroms": sp1_chroms,
        "species2_chroms": sp2_chroms,
    }


def create_control_files(
    comparison_dir: Path,
    chrom_info: Dict[str, List[str]],
    label1: str,
    label2: str,
) -> None:
    """Create MCScanX visualization control files.

    Creates dot.ctl, circle.ctl, dual_synteny.ctl, and bar.ctl.

    Args:
        comparison_dir: Output directory
        chrom_info: Chromosome info from extract_chromosome_info()
        label1: Display label for species 1
        label2: Display label for species 2
    """
    sp1_chroms = chrom_info["species1_chroms"]
    sp2_chroms = chrom_info["species2_chroms"]
    sp1_str = ",".join(sp1_chroms)
    sp2_str = ",".join(sp2_chroms)

    # dot.ctl
    with open(comparison_dir / "dot.ctl", "w") as f:
        f.write("1200\t//dimension (in pixels) of x axis\n")
        f.write("1200\t//dimension (in pixels) of y axis\n")
        f.write(f"{sp1_str}\t//chromosomes in x axis ({label1})\n")
        f.write(f"{sp2_str}\t//chromosomes in y axis ({label2})\n")
    print(f"  Created dot.ctl")

    # circle.ctl — interleave chromosomes for better visualization
    interleaved = []
    for i in range(max(len(sp1_chroms), len(sp2_chroms))):
        if i < len(sp1_chroms):
            interleaved.append(sp1_chroms[i])
        if i < len(sp2_chroms):
            interleaved.append(sp2_chroms[i])
    with open(comparison_dir / "circle.ctl", "w") as f:
        f.write("1000\t//plot width and height (in pixels)\n")
        f.write(f"{','.join(interleaved)}\t//chromosomes in the circle\n")
    print(f"  Created circle.ctl")

    # dual_synteny.ctl
    with open(comparison_dir / "dual_synteny.ctl", "w") as f:
        f.write("800\t//plot width (in pixels)\n")
        f.write("1200\t//plot height (in pixels)\n")
        f.write(f"{sp1_str}\t//chromosomes in the left column ({label1})\n")
        f.write(f"{sp2_str}\t//chromosomes in the right column ({label2})\n")
    print(f"  Created dual_synteny.ctl")

    # bar.ctl
    with open(comparison_dir / "bar.ctl", "w") as f:
        f.write("1200\t//dimension (in pixels) of x axis\n")
        f.write("800\t//dimension (in pixels) of y axis\n")
        f.write(f"{sp1_str}\t//reference chromosomes ({label1})\n")
        f.write(f"{sp2_str}\t//target chromosomes ({label2})\n")
    print(f"  Created bar.ctl")


def generate_visualizations(
    comparison_dir: Path,
    comparison_name: str,
    tools_dir: Path,
) -> None:
    """Generate MCScanX Java visualizations.

    Creates dot plot, dual synteny plot, circle plot, and bar plot.

    Args:
        comparison_dir: Directory with MCScanX results and control files
        comparison_name: Comparison identifier (e.g. 'al_sm')
        tools_dir: Path to MCScanX downstream_analyses directory
    """
    print("Generating visualizations...")

    gff_file = comparison_dir / f"{comparison_name}.gff"
    col_file = comparison_dir / f"{comparison_name}.collinearity"

    if not gff_file.exists() or not col_file.exists():
        print("  Missing GFF or collinearity file, skipping.", file=sys.stderr)
        return

    viz_types = [
        {
            "name": "Dot Plot",
            "class": "dot_plotter",
            "control": "dot.ctl",
            "output": "dot.PNG",
        },
        {
            "name": "Dual Synteny Plot",
            "class": "dual_synteny_plotter",
            "control": "dual_synteny.ctl",
            "output": "dual_synteny.PNG",
        },
        {
            "name": "Circle Plot",
            "class": "circle_plotter",
            "control": "circle.ctl",
            "output": "circle.PNG",
        },
        {
            "name": "Bar Plot",
            "class": "bar_plotter",
            "control": "bar.ctl",
            "output": "bar.PNG",
        },
    ]

    for viz in viz_types:
        ctl_file = comparison_dir / viz["control"]
        out_file = comparison_dir / viz["output"]

        if not ctl_file.exists():
            print(f"  ✗ {viz['name']}: control file {viz['control']} not found")
            continue

        cmd = (
            f"java -cp {tools_dir} {viz['class']} "
            f"-g {gff_file} -s {col_file} "
            f"-c {ctl_file} -o {out_file}"
        )
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0 and out_file.exists():
            size_kb = out_file.stat().st_size / 1024
            print(f"  ✓ {viz['name']} ({size_kb:.1f} KB)")
        else:
            print(f"  ✗ {viz['name']} failed")
            if result.stderr:
                first_line = result.stderr.strip().split("\n")[0][:80]
                print(f"    Error: {first_line}")


# ---------------------------------------------------------------------------
# Results reporting
# ---------------------------------------------------------------------------


def analyze_results(comparison_dir: Path, comparison_name: str) -> None:
    """Print comprehensive results analysis.

    Args:
        comparison_dir: Directory with MCScanX results
        comparison_name: Comparison identifier
    """
    print(f"\n{'=' * 60}")
    print(f"  RESULTS SUMMARY: {comparison_name.upper()}")
    print(f"{'=' * 60}\n")

    col_file = comparison_dir / f"{comparison_name}.collinearity"
    tandem_file = comparison_dir / f"{comparison_name}.tandem"

    if col_file.exists():
        with open(col_file, "r") as f:
            lines = f.readlines()

        blocks = [l for l in lines if l.startswith("##")]
        pairs = [l for l in lines if not l.startswith("#")]

        print(f"Synteny Analysis:")
        print(f"  Synteny blocks: {len(blocks)}")
        print(f"  Syntenic gene pairs: {len(pairs)}")

        if blocks:
            block_sizes = []
            for block in blocks:
                if "N=" in block:
                    try:
                        size = int(block.split("N=")[1].split()[0])
                        block_sizes.append(size)
                    except (ValueError, IndexError):
                        pass

            if block_sizes:
                avg = sum(block_sizes) / len(block_sizes)
                print(f"  Average block size: {avg:.1f} genes")
                print(f"  Largest block: {max(block_sizes)} genes")
                print(f"  Smallest block: {min(block_sizes)} genes")
    else:
        print("  No collinearity file found.", file=sys.stderr)

    if tandem_file.exists():
        with open(tandem_file, "r") as f:
            tandem_count = sum(1 for _ in f)
        print(f"\nTandem Duplication:")
        print(f"  Tandem gene pairs: {tandem_count}")

    # Visualizations
    viz_files = ["dot.PNG", "dual_synteny.PNG", "circle.PNG", "bar.PNG"]
    available = []
    for vf in viz_files:
        vp = comparison_dir / vf
        if vp.exists():
            size_kb = vp.stat().st_size / 1024
            available.append(f"{vf} ({size_kb:.1f} KB)")

    print(f"\nVisualizations:")
    if available:
        for v in available:
            print(f"  ✓ {v}")
    else:
        print("  None generated")

    # HTML files
    html_files = [f for f in comparison_dir.iterdir() if f.suffix == ".html"]
    if html_files:
        print(f"\nHTML Reports: {len(html_files)} files")

    print(f"\nOutput directory: {comparison_dir}")
    print(f"{'=' * 60}\n")


# ---------------------------------------------------------------------------
# Main pipeline orchestration
# ---------------------------------------------------------------------------


def run_pipeline(
    species1: str,
    species2: str,
    braker_base: Path,
    output_dir: Path,
    blastp_cmd: str,
    makeblastdb_cmd: str,
    mcscanx_cmd: str,
    mcscanx_tools_dir: Optional[Path],
    threads: int = 8,
    evalue: float = 1e-10,
    skip_blast: bool = False,
    skip_viz: bool = False,
) -> None:
    """Run the complete MCScanX pairwise synteny pipeline.

    Args:
        species1: First species name (Genus_species)
        species2: Second species name (Genus_species)
        braker_base: Base directory with {species}/braker.aa and braker.gff3
        output_dir: Output directory for all results
        blastp_cmd: Path to blastp binary
        makeblastdb_cmd: Path to makeblastdb binary
        mcscanx_cmd: Path to MCScanX binary
        mcscanx_tools_dir: Path to MCScanX downstream_analyses dir (for viz)
        threads: Number of BLAST threads
        evalue: BLAST E-value threshold
        skip_blast: Skip BLAST if output exists
        skip_viz: Skip visualization generation
    """
    prefix1 = get_species_prefix(species1)
    prefix2 = get_species_prefix(species2)
    comparison_name = f"{prefix1}_{prefix2}"

    # Resolve input file paths
    sp1_aa = braker_base / species1 / "braker.aa"
    sp2_aa = braker_base / species2 / "braker.aa"
    sp1_gff3 = braker_base / species1 / "braker.gff3"
    sp2_gff3 = braker_base / species2 / "braker.gff3"

    print(f"\n{'=' * 60}")
    print(f"  MCScanX Pairwise Synteny Analysis")
    print(f"{'=' * 60}")
    print(f"  Species 1: {get_display_name(species1)} ({prefix1.upper()})")
    print(f"  Species 2: {get_display_name(species2)} ({prefix2.upper()})")
    print(f"  Comparison: {comparison_name}")
    print(f"  Output: {output_dir}")
    print(f"{'=' * 60}\n")

    # Verify input files
    print("Verifying input files...")
    missing = []
    for fp in [sp1_aa, sp2_aa, sp1_gff3, sp2_gff3]:
        if fp.exists():
            size_mb = fp.stat().st_size / (1024 * 1024)
            print(f"  ✓ {fp} ({size_mb:.1f} MB)")
        else:
            print(f"  ✗ {fp} (NOT FOUND)", file=sys.stderr)
            missing.append(fp)

    if missing:
        print(
            f"\nError: {len(missing)} required input files missing. Aborting.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Create output directories
    blast_db_dir = output_dir / "blast_db"
    blastp_dir = output_dir / "blastp"
    comparison_dir = output_dir / comparison_name
    for d in [blast_db_dir, blastp_dir, comparison_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Step 1: Create combined protein dataset
    # -----------------------------------------------------------------------
    print(f"\n--- Step 1/7: Combined Protein Dataset ---\n")
    combined_faa = blast_db_dir / f"{comparison_name}.faa"
    create_pairwise_dataset(sp1_aa, sp2_aa, prefix1, prefix2, combined_faa)

    # -----------------------------------------------------------------------
    # Step 2: Build BLAST database
    # -----------------------------------------------------------------------
    print(f"\n--- Step 2/7: BLAST Database ---\n")
    db_prefix = blast_db_dir / comparison_name
    blast_output = blastp_dir / f"{comparison_name}.blast"

    if skip_blast and blast_output.exists():
        print(f"  Skipping database creation (--skip-blast, output exists)")
    else:
        if not create_blast_database(combined_faa, db_prefix, makeblastdb_cmd):
            print("Error: BLAST database creation failed. Aborting.", file=sys.stderr)
            sys.exit(1)

    # -----------------------------------------------------------------------
    # Step 3: BLASTP search
    # -----------------------------------------------------------------------
    print(f"\n--- Step 3/7: BLASTP Search ---\n")
    if skip_blast and blast_output.exists():
        size_mb = blast_output.stat().st_size / (1024 * 1024)
        print(f"  Skipping BLAST (--skip-blast, output exists: {size_mb:.1f} MB)")
    else:
        if not run_blast_search(
            combined_faa, db_prefix, blast_output, blastp_cmd, threads, evalue
        ):
            print("Error: BLAST search failed. Aborting.", file=sys.stderr)
            sys.exit(1)

    # -----------------------------------------------------------------------
    # Step 4: GFF3 -> BED conversion
    # -----------------------------------------------------------------------
    print(f"\n--- Step 4/7: GFF3 to BED Conversion ---\n")
    sp1_bed = comparison_dir / f"{prefix1}.bed"
    sp2_bed = comparison_dir / f"{prefix2}.bed"
    convert_gff3_to_bed(sp1_gff3, sp1_bed)
    convert_gff3_to_bed(sp2_gff3, sp2_bed)

    # -----------------------------------------------------------------------
    # Step 5: Process BED files (standardize chromosome names)
    # -----------------------------------------------------------------------
    print(f"\n--- Step 5/7: BED Processing ---\n")
    sp1_processed = comparison_dir / f"{prefix1}-filtered-modified.bed"
    sp2_processed = comparison_dir / f"{prefix2}-filtered-modified.bed"

    sp1_genes = process_bed_file(sp1_bed, sp1_processed, prefix1)
    sp2_genes = process_bed_file(sp2_bed, sp2_processed, prefix2)
    print(f"\n  Total: {sp1_genes} + {sp2_genes} = {sp1_genes + sp2_genes} genes")

    # -----------------------------------------------------------------------
    # Step 6: Create MCScanX input files
    # -----------------------------------------------------------------------
    print(f"\n--- Step 6/7: MCScanX Inputs ---\n")
    if not create_mcscanx_inputs(
        comparison_dir, comparison_name, prefix1, prefix2, blast_output
    ):
        print("Error: MCScanX input creation failed. Aborting.", file=sys.stderr)
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Step 7: Run MCScanX
    # -----------------------------------------------------------------------
    print(f"\n--- Step 7/7: MCScanX Analysis ---\n")
    input_prefix = comparison_dir / comparison_name
    if not run_mcscanx(input_prefix, mcscanx_cmd):
        print("Error: MCScanX analysis failed. Aborting.", file=sys.stderr)
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Optional: Visualizations
    # -----------------------------------------------------------------------
    if not skip_viz and mcscanx_tools_dir:
        print(f"\n--- Visualization ---\n")

        # Compile Java tools
        compile_success = compile_visualization_tools(mcscanx_tools_dir)

        if compile_success:
            gff_file = comparison_dir / f"{comparison_name}.gff"
            chrom_info = extract_chromosome_info(gff_file, prefix1, prefix2)

            if chrom_info:
                label1 = get_display_name(species1, abbreviated=True)
                label2 = get_display_name(species2, abbreviated=True)
                create_control_files(comparison_dir, chrom_info, label1, label2)
                generate_visualizations(
                    comparison_dir, comparison_name, mcscanx_tools_dir
                )
        else:
            print("  Java compilation failed; skipping visualizations.")
    elif skip_viz:
        print("\n  Skipping visualizations (--skip-viz)")

    # -----------------------------------------------------------------------
    # Results summary
    # -----------------------------------------------------------------------
    analyze_results(comparison_dir, comparison_name)
    print("Pipeline completed successfully.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MCScanX pairwise synteny analysis between two species.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Basic usage (uses project default paths):
  %(prog)s Aspiorhynchus_laticeps Schizothorax_macropogon

  # With explicit braker base directory:
  %(prog)s Aspiorhynchus_laticeps Schizothorax_macropogon \\
      --braker-base /home/jovyan/Outputs/Preprocessing/BRAKER_UNMASKED

  # Reuse existing BLAST results and skip visualizations:
  %(prog)s Aspiorhynchus_laticeps Schizothorax_macropogon \\
      --skip-blast --skip-viz

  # Custom output directory and thread count:
  %(prog)s Aspiorhynchus_laticeps Triplophysa_bombifrons \\
      --output-dir ./my_synteny_results --threads 16
""",
    )

    parser.add_argument(
        "species1",
        help="First species (Genus_species format, e.g. Aspiorhynchus_laticeps)",
    )
    parser.add_argument(
        "species2",
        help="Second species (Genus_species format, e.g. Schizothorax_macropogon)",
    )

    # Input paths
    parser.add_argument(
        "--braker-base",
        type=Path,
        default=None,
        help=(
            "Base directory containing {species}/braker.aa and braker.gff3. "
            "Default: auto-detect from config or ./Outputs/Preprocessing/BRAKER_UNMASKED"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory for all results. "
            "Default: ./Outputs/MCScanX_Pairwise/{prefix1}_{prefix2}"
        ),
    )

    # Tool paths
    parser.add_argument(
        "--blastp",
        default="./Apps/Blast/ncbi-blast-2.16.0+/bin/blastp",
        help="Path to blastp binary (default: %(default)s)",
    )
    parser.add_argument(
        "--makeblastdb",
        default="./Apps/Blast/ncbi-blast-2.16.0+/bin/makeblastdb",
        help="Path to makeblastdb binary (default: %(default)s)",
    )
    parser.add_argument(
        "--mcscanx",
        default="./Apps/MCScanX/V1/MCScanX-1.0.0/MCScanX",
        help="Path to MCScanX binary (default: %(default)s)",
    )
    parser.add_argument(
        "--mcscanx-tools",
        type=Path,
        default=Path("./Apps/MCScanX/V1/MCScanX-1.0.0/downstream_analyses"),
        help="Path to MCScanX downstream_analyses dir (default: %(default)s)",
    )

    # Pipeline options
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads for BLAST (default: %(default)s)",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-10,
        help="BLAST E-value threshold (default: %(default)s)",
    )
    parser.add_argument(
        "--skip-blast",
        action="store_true",
        help="Skip BLAST steps if output already exists",
    )
    parser.add_argument(
        "--skip-viz",
        action="store_true",
        help="Skip Java visualization generation",
    )

    args = parser.parse_args()

    # Resolve braker base directory
    braker_base = args.braker_base
    if braker_base is None:
        # Try config system first
        try:
            sys.path.insert(0, "./Notebooks/config")
            from load_config import get_path

            braker_base_str = get_path("cds.masked_dir") or get_path("cds.base_dir")
            if braker_base_str:
                braker_base = Path(braker_base_str)
        except (ImportError, Exception):
            pass

        # Fallback to default project structure
        if braker_base is None or not braker_base.exists():
            script_dir = Path(__file__).resolve().parent
            project_root = script_dir.parent
            braker_base = project_root / "Outputs" / "Preprocessing" / "BRAKER_UNMASKED"

    if not braker_base.exists():
        print(f"Error: BRAKER base directory not found: {braker_base}", file=sys.stderr)
        print("Specify --braker-base explicitly.", file=sys.stderr)
        sys.exit(1)

    # Resolve output directory
    output_dir = args.output_dir
    if output_dir is None:
        prefix1 = get_species_prefix(args.species1)
        prefix2 = get_species_prefix(args.species2)
        script_dir = Path(__file__).resolve().parent
        project_root = script_dir.parent
        output_dir = (
            project_root / "Outputs" / "MCScanX_Pairwise" / f"{prefix1}_{prefix2}"
        )

    # Run the pipeline
    run_pipeline(
        species1=args.species1,
        species2=args.species2,
        braker_base=braker_base,
        output_dir=output_dir,
        blastp_cmd=args.blastp,
        makeblastdb_cmd=args.makeblastdb,
        mcscanx_cmd=args.mcscanx,
        mcscanx_tools_dir=args.mcscanx_tools if not args.skip_viz else None,
        threads=args.threads,
        evalue=args.evalue,
        skip_blast=args.skip_blast,
        skip_viz=args.skip_viz,
    )


if __name__ == "__main__":
    main()
