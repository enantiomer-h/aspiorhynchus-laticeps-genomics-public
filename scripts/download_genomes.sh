#!/usr/bin/env bash
#
# download_genomes.sh — fetch the 12 comparative-set genome assemblies from NCBI.
#
# Reads data/accessions.tsv and downloads each assembly with the NCBI Datasets
# CLI into DB/GENOME_Comparison/<Species>/, matching the layout the notebooks
# and scripts expect (…/ncbi_dataset/data/<accession>/*.fna). BRAKER3 then
# re-annotates these genomes (see Notebooks/0-ultimate-preprocessing.qmd).
#
# Requirements: the NCBI 'datasets' CLI and 'unzip' on PATH.
#   conda install -c conda-forge ncbi-datasets-cli
#
# Usage:
#   bash scripts/download_genomes.sh [OUTDIR]
# OUTDIR defaults to DB/GENOME_Comparison relative to the repository root.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ACCESSIONS="${REPO_ROOT}/data/accessions.tsv"
OUTDIR="${1:-${REPO_ROOT}/DB/GENOME_Comparison}"

command -v datasets >/dev/null 2>&1 || {
  echo "ERROR: NCBI 'datasets' CLI not found. Install with:" >&2
  echo "       conda install -c conda-forge ncbi-datasets-cli" >&2
  exit 1
}
[[ -f "$ACCESSIONS" ]] || { echo "ERROR: $ACCESSIONS not found." >&2; exit 1; }

mkdir -p "$OUTDIR"
skipped=()

# Read tab-separated rows, skipping comments and the header line.
while IFS=$'\t' read -r species accession _rest; do
  [[ -z "${species:-}" || "${species:0:1}" == "#" || "$species" == "species" ]] && continue

  case "$accession" in
    GCA_*|GCF_*) : ;;
    *)  # placeholder (e.g. NEEDS_TABLE_S1) — cannot download yet
        skipped+=("$species ($accession)")
        continue ;;
  esac

  dest="${OUTDIR}/${species}"
  if [[ -d "${dest}/ncbi_dataset/data/${accession}" ]]; then
    echo "[skip] ${species}: ${accession} already present"
    continue
  fi

  echo "[get ] ${species}: ${accession}"
  tmpzip="$(mktemp --suffix=.zip)"
  datasets download genome accession "$accession" --include genome --filename "$tmpzip"
  mkdir -p "$dest"
  unzip -o -q "$tmpzip" -d "$dest"
  rm -f "$tmpzip"
done < "$ACCESSIONS"

if ((${#skipped[@]})); then
  echo
  echo "WARNING: the following species have no downloadable accession yet and were skipped:" >&2
  printf '  - %s\n' "${skipped[@]}" >&2
  echo "Fill in their accessions in data/accessions.tsv (from the manuscript Table S1)," >&2
  echo "then re-run this script." >&2
fi

echo "Done. Genomes written under: ${OUTDIR}"
