#!/bin/bash
# Tier 2B/2C support: Run RepeatMasker on AL chromosomal assembly so that
# (a) we get a chromosome-coordinate repeatmasker output for the AL circos
#     repeat track, and (b) the *.tbl summary feeds the comparative TE
#     classification table.
#
# Memory: ~2 GB; CPU: -pa 2 (uses ~8 cores) — leaves headroom for emapper +
# ncRNA pipelines. Runtime: ~12-24 h on a 1.2 Gb fish genome.
set -uo pipefail

ENV="genome_annotation_env"
SP="Aspiorhynchus_laticeps"
FNA="/home/jovyan/DB/GENOME_Comparison/Aspiorhynchus_laticeps/ncbi_dataset/data/GCA_023376895.1/GCA_023376895.1_ASM2337689v1_genomic.fna"
OUT="/home/jovyan/Outputs/RepeatMasker_chromosomal/${SP}"
mkdir -p "$OUT"

if [[ -s "$OUT/$(basename $FNA).tbl" ]]; then
  echo "[skip] $SP: already has $OUT/$(basename $FNA).tbl"
  exit 0
fi

echo "[RM] $SP $(date)"
conda run -n "$ENV" --no-capture-output \
  RepeatMasker \
    -species "actinopterygii" \
    -pa 2 \
    -gff \
    -dir "$OUT" \
    "$FNA" 2>&1 | tail -200
echo "[RM] $SP DONE $(date)"
echo "Outputs:"
ls -la "$OUT/"
