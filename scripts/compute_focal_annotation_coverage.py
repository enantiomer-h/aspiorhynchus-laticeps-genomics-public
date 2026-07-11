#!/usr/bin/env python3
"""Per-focal-species functional annotation coverage from EggNOG-mapper outputs.

Each set's eggnog run was performed on the focal species' BRAKER_MASKED CDS:
  - set2 → Aspiorhynchus_laticeps
  - set3 → Diptychus_maculatus

Produces a Table 10-style coverage summary (analog of Zhao et al. 2026, Table 10)
for each focal species, written to Tables_set{2,3}/annotation_coverage_focal.csv.

The other 10 species in each set's roster do NOT have EggNOG annotations on
disk (would require running emapper on each proteome). That's documented in
the report as deferred work.
"""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SETS, PROJECT_ROOT, braker_cds  # noqa: E402

EGGNOG = {
    "set2": PROJECT_ROOT / "Outputs_set2/Eggnog/Eggnog.emapper.annotations",
    "set3": PROJECT_ROOT / "Outputs_set3/Eggnog/Eggnog.emapper.annotations",
}

COLUMNS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]


def count_proteins_in_cds(cds_fa: Path) -> int:
    if not cds_fa.is_file():
        return -1
    n = 0
    with cds_fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def has_value(s: str) -> bool:
    return s not in ("-", "", "nan") and not (isinstance(s, float) and pd.isna(s))


def coverage_row(focal: str, eggnog_path: Path) -> dict:
    cds_count = count_proteins_in_cds(braker_cds(focal))
    df = pd.read_csv(eggnog_path, sep="\t", comment="#", header=None, names=COLUMNS, dtype=str, keep_default_na=False)
    n = len(df)

    def pct_with(col: str) -> tuple[int, float]:
        m = df[col].apply(has_value).sum()
        return int(m), round(100.0 * m / cds_count, 1) if cds_count > 0 else 0.0

    eggnog_n, eggnog_pct = (n, round(100.0 * n / cds_count, 1) if cds_count > 0 else 0.0)
    cog_n, cog_pct = pct_with("COG_category")
    desc_n, desc_pct = pct_with("Description")
    name_n, name_pct = pct_with("Preferred_name")
    go_n, go_pct = pct_with("GOs")
    ko_n, ko_pct = pct_with("KEGG_ko")
    path_n, path_pct = pct_with("KEGG_Pathway")
    pfam_n, pfam_pct = pct_with("PFAMs")
    cazy_n, cazy_pct = pct_with("CAZy")
    ec_n, ec_pct = pct_with("EC")

    return {
        "Focal species": focal,
        "Total CDS (from BRAKER_MASKED)": cds_count,
        "EggNOG annotated": f"{eggnog_n} ({eggnog_pct}%)",
        "COG category": f"{cog_n} ({cog_pct}%)",
        "Description": f"{desc_n} ({desc_pct}%)",
        "Preferred gene name": f"{name_n} ({name_pct}%)",
        "GO terms": f"{go_n} ({go_pct}%)",
        "KEGG KO": f"{ko_n} ({ko_pct}%)",
        "KEGG Pathway": f"{path_n} ({path_pct}%)",
        "PFAMs": f"{pfam_n} ({pfam_pct}%)",
        "CAZy": f"{cazy_n} ({cazy_pct}%)",
        "EC numbers": f"{ec_n} ({ec_pct}%)",
    }


def main() -> None:
    for s, paths in SETS.items():
        focal = paths["focal"]
        eggnog_path = EGGNOG[s]
        if not eggnog_path.is_file():
            print(f"[{s}] missing {eggnog_path}; skip")
            continue
        print(f"[{s}] focal={focal}, eggnog={eggnog_path}")
        row = coverage_row(focal, eggnog_path)
        df = pd.DataFrame([row]).set_index("Focal species")
        print(df.T.to_string())

        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        out = paths["tables_dir"] / "annotation_coverage_focal.csv"
        df.to_csv(out)
        print(f"[{s}] wrote {out}")
        print()


if __name__ == "__main__":
    main()
