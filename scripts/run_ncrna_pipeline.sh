#!/bin/bash
# Tier 2A: Comparative ncRNA inventory — tRNAscan-SE + barrnap on 9 species with genome FNAs.
# Skips Infernal+Rfam scan (too slow for this session); Rfam is cmpressed for follow-up.
set -uo pipefail

CONDA_ENV="env_biopython"
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

run_barrnap () {
  local sp="$1" fna="$2"
  local out="$OUT/$sp/barrnap.gff3"
  if [[ -s "$out" ]]; then echo "[skip barrnap] $sp"; return; fi
  mkdir -p "$OUT/$sp"
  echo "[barrnap] $sp $(date +%H:%M:%S)"
  conda run -n "$CONDA_ENV" --no-capture-output barrnap --kingdom euk --threads 4 "$fna" > "$out" 2> "$OUT/$sp/barrnap.log"
  echo "[barrnap] $sp DONE $(date +%H:%M:%S)"
}

run_trnascan () {
  local sp="$1" fna="$2"
  local out="$OUT/$sp/trnascan.tsv"
  if [[ -s "$out" ]]; then echo "[skip trnascan] $sp"; return; fi
  mkdir -p "$OUT/$sp"
  echo "[trnascan] $sp $(date +%H:%M:%S)"
  conda run -n "$CONDA_ENV" --no-capture-output tRNAscan-SE -E --thread 4 -o "$out" -f "$OUT/$sp/trnascan.struct" -m "$OUT/$sp/trnascan.stats" "$fna" 2> "$OUT/$sp/trnascan.log" || echo "[trnascan ERR] $sp"
  echo "[trnascan] $sp DONE $(date +%H:%M:%S)"
}

for sp in "${!GENOMES[@]}"; do
  fna="${GENOMES[$sp]}"
  if [[ ! -s "$fna" ]]; then echo "[skip] $sp: missing $fna"; continue; fi
  run_barrnap "$sp" "$fna"
  run_trnascan "$sp" "$fna"
done

echo "ALL DONE $(date)"
