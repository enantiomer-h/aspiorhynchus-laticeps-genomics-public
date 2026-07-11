#!/bin/bash
# Run eggnog-mapper on the 10 non-focal species for Tier 3C.
# Focal species (AL, DM) already have emapper runs — we only do the rest.
set -euo pipefail

SPECIES=(
  "Carassius_auratus"
  "Cyprinus_carpio"
  "Danio_rerio"
  "Gymnocypris_eckloni"
  "Oxygymnocypris_stewartii"
  "Schizopygopsis_younghusbandi"
  "Sinocyclocheilus_grahami"
  "Triplophysa_pappenheimi"
  "Triplophysa_tibetana"
  "Triplophysa_yaopeizhii"
)

OUT_BASE="/home/jovyan/Outputs_set2/Eggnog/per_species"
mkdir -p "$OUT_BASE"

for sp in "${SPECIES[@]}"; do
  CDS="/home/jovyan/Outputs/Preprocessing/BRAKER_MASKED/${sp}/braker.codingseq"
  OUT="$OUT_BASE/${sp}"
  if [[ -f "${OUT}.emapper.annotations" ]]; then
    echo "[skip] ${sp}: existing annotations at ${OUT}.emapper.annotations"
    continue
  fi
  if [[ ! -f "$CDS" ]]; then
    echo "[skip] ${sp}: CDS missing at $CDS"
    continue
  fi
  echo "[run] ${sp}: $(date)"
  conda run -n eggnog_env --no-capture-output emapper.py \
    --data_dir /home/jovyan/DB/Eggnog \
    --itype CDS \
    -i "$CDS" \
    --output "$OUT" \
    -m diamond \
    --cpu 12 \
    --override
  echo "[done] ${sp}: $(date)"
done

echo "ALL DONE $(date)"
