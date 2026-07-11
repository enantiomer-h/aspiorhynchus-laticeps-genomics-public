#!/usr/bin/env python3
"""Aggregate ncRNA outputs (tRNAscan-SE + cmscan rRNA) into a per-species table.

Output: Notebooks/Tables_set{2,3}/ncrna_inventory.csv
Format mirrors reference paper Table 7.
"""
import sys
from pathlib import Path
import re
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SETS  # noqa: E402

PROJECT_ROOT = Path("/home/jovyan")
NCRNA_DIR = PROJECT_ROOT / "Outputs/ncRNA"

SPECIES_AVAILABLE = [
    "Aspiorhynchus_laticeps", "Carassius_auratus", "Cyprinus_carpio",
    "Danio_rerio", "Diptychus_maculatus", "Gymnocypris_eckloni",
    "Oxygymnocypris_stewartii", "Sinocyclocheilus_grahami", "Triplophysa_tibetana",
]


def parse_trnascan(tsv: Path) -> dict[str, int]:
    """Count tRNAs by isotype from a tRNAscan-SE -o output (header: Sequence Name, ...)."""
    if not tsv.is_file() or tsv.stat().st_size == 0:
        return {"tRNA total": 0}
    n_total = 0
    with tsv.open() as fh:
        for line in fh:
            if line.startswith(("#", "Sequence", "Name")) or not line.strip(): continue
            parts = line.split()
            if len(parts) < 5: continue
            n_total += 1
    return {"tRNA total": n_total}


def parse_cmscan_rRNA(tbl: Path) -> dict[str, int]:
    """Count rRNAs by class from a cmscan --tblout file."""
    counts = {"5S": 0, "5.8S": 0, "SSU/18S": 0, "LSU/28S": 0, "SSU_bact": 0, "LSU_bact": 0}
    if not tbl.is_file() or tbl.stat().st_size == 0:
        return {"rRNA " + k: 0 for k in counts}
    with tbl.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            parts = re.split(r"\s+", line.strip(), maxsplit=17)
            if len(parts) < 17: continue
            target = parts[0]
            if "5S_rRNA" in target: counts["5S"] += 1
            elif "5_8S_rRNA" in target: counts["5.8S"] += 1
            elif "SSU_rRNA_eukarya" in target: counts["SSU/18S"] += 1
            elif "LSU_rRNA_eukarya" in target: counts["LSU/28S"] += 1
            elif "SSU_rRNA_bacteria" in target: counts["SSU_bact"] += 1
            elif "LSU_rRNA_bacteria" in target: counts["LSU_bact"] += 1
    return {"rRNA " + k: v for k, v in counts.items()}


def main() -> None:
    rows = []
    for sp in SPECIES_AVAILABLE:
        sp_dir = NCRNA_DIR / sp
        row = {"Species": sp}
        row.update(parse_trnascan(sp_dir / "trnascan.tsv"))
        row.update(parse_cmscan_rRNA(sp_dir / "cmscan_rRNA.tbl"))
        rows.append(row)

    df = pd.DataFrame(rows).set_index("Species")
    print(df.to_string())

    for s, paths in SETS.items():
        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        out = paths["tables_dir"] / "ncrna_inventory.csv"
        df.to_csv(out)
        print(f"[{s}] wrote {out}")


if __name__ == "__main__":
    main()
