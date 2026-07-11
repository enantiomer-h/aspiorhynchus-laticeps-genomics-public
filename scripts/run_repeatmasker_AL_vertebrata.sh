#!/bin/bash
# Tier 2B follow-up: re-run AL RepeatMasker with the broader -species vertebrata
# library (Dfam-actinopterygii was too sparse — only 4.38% masked). Vertebrata
# spans ~50K HMMs vs actinopterygii's ~few hundred.
set -uo pipefail

ENV="genome_annotation_env"
SP="Aspiorhynchus_laticeps"
FNA="/home/jovyan/DB/GENOME_Comparison/Aspiorhynchus_laticeps/ncbi_dataset/data/GCA_023376895.1/GCA_023376895.1_ASM2337689v1_genomic.fna"
OUT="/home/jovyan/Outputs/RepeatMasker_chromosomal_vertebrata/${SP}"
mkdir -p "$OUT"

if [[ -s "$OUT/$(basename $FNA).tbl" ]]; then
  echo "[skip] $SP: existing $OUT/$(basename $FNA).tbl"
  exit 0
fi

echo "[RM-vert] $SP $(date)"
conda run -n "$ENV" --no-capture-output \
  RepeatMasker \
    -species "vertebrata" \
    -pa 2 \
    -gff \
    -dir "$OUT" \
    "$FNA" 2>&1 | tail -250
echo "[RM-vert] $SP DONE $(date)"
ls -la "$OUT/"
