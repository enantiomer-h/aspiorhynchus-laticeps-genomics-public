#!/bin/bash
# Tier 2A complement: cmscan against an Rfam subset (rRNA models only) for
# eukaryote 5S/5.8S/SSU/LSU rRNA detection. Faster than full Rfam scan (~10
# min per genome). Complements barrnap (which lacks a euk-default HMM).
set -uo pipefail

ENV="env_biopython"
RFAM="/home/jovyan/DB/Rfam/Rfam_rRNA.cm"
OUT="/home/jovyan/Outputs/ncRNA"
mkdir -p "$OUT"

declare -A GENOMES
GENOMES[Aspiorhynchus_laticeps]="/home/jovyan/DB/GENOME_Comparison/Aspiorhynchus_laticeps/ncbi_dataset/data/GCA_023376895.1/GCA_023376895.1_ASM2337689v1_genomic.fna"
GENOMES[Carassius_auratus]="/home/jovyan/DB/GENOME_Comparison/Carassius_auratus/ncbi_dataset/data/GCA_003368295.1/GCA_003368295.1_ASM336829v1_genomic.fna"
GENOMES[Cyprinus_carpio]="/home/jovyan/DB/GENOME_Comparison/Cyprinus_carpio/ncbi_dataset/data/GCA_018340385.1/GCA_018340385.1_ASM1834038v1_genomic.fna"
GENOMES[Danio_rerio]="/home/jovyan/DB/GENOME_Comparison/Danio_rerio/ncbi_dataset/data/GCF_000002035.6/GCF_000002035.6_GRCz11_genomic.fna"
GENOMES[Diptychus_maculatus]="/home/jovyan/DB/GENOME_Comparison/Diptychus_maculatus/DM.fasta"
GENOMES[Gymnocypris_eckloni]="/home/jovyan/DB/GENOME_Comparison/Gymnocypris_eckloni/ncbi_dataset/data/GCA_027564155.1/GCA_027564155.1_ASM2756415v1_genomic.fna"
GENOMES[Oxygymnocypris_stewartii]="/home/jovyan/DB/GENOME_Comparison/Oxygymnocypris_stewartii/ncbi_dataset/data/GCA_003573665.1/GCA_003573665.1_Novo_Ost_1.0_genomic.fna"
GENOMES[Sinocyclocheilus_grahami]="/home/jovyan/DB/GENOME_Comparison/Sinocyclocheilus_grahami/ncbi_dataset/data/GCA_001515645.1/GCA_001515645.1_SAMN03320097.WGS_v1.1_genomic.fna"
GENOMES[Triplophysa_tibetana]="/home/jovyan/DB/GENOME_Comparison/Triplophysa_tibetana/ncbi_dataset/data/GCA_008369825.1/GCA_008369825.1_ASM836982v1_genomic.fna"

for sp in "${!GENOMES[@]}"; do
  fna="${GENOMES[$sp]}"
  out="$OUT/$sp/cmscan_rRNA.tbl"
  if [[ -s "$out" ]]; then echo "[skip] $sp"; continue; fi
  if [[ ! -s "$fna" ]]; then echo "[skip missing] $sp $fna"; continue; fi
  mkdir -p "$OUT/$sp"
  echo "[rRNA cmscan] $sp $(date +%H:%M:%S)"
  conda run -n "$ENV" --no-capture-output cmscan \
    --cpu 4 --tblout "$out" --noali \
    "$RFAM" "$fna" \
    > "$OUT/$sp/cmscan_rRNA.out" 2> "$OUT/$sp/cmscan_rRNA.log" || echo "[err] $sp"
  echo "[rRNA cmscan] $sp DONE $(date +%H:%M:%S)"
done

echo "ALL rRNA scans DONE $(date)"
