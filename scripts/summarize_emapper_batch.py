#!/usr/bin/env python3
"""Aggregate per-species emapper outputs (set2 + set3 focals + 10 batch species)
into one Tier 3C coverage table covering all 12 species.

Reads from Outputs_set2/Eggnog/per_species/{species}.emapper.annotations
plus the original focal-only outputs at Outputs_set{2,3}/Eggnog/Eggnog.emapper.annotations.
"""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SPECIES, SETS, braker_cds  # noqa: E402

PROJECT_ROOT = Path("/home/jovyan")
PER_SPECIES_DIR = PROJECT_ROOT / "Outputs_set2/Eggnog/per_species"
FOCAL_FILES = {
    "Aspiorhynchus_laticeps": PROJECT_ROOT / "Outputs_set2/Eggnog/Eggnog.emapper.annotations",
    "Diptychus_maculatus": PROJECT_ROOT / "Outputs_set3/Eggnog/Eggnog.emapper.annotations",
}
COLUMNS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]


def cds_count(species: str) -> int:
    cds = braker_cds(species)
    if not cds.is_file(): return 0
    return sum(1 for line in cds.open() if line.startswith(">"))


def has(s: str) -> bool: return s not in ("", "-", "nan")


def coverage_row(species: str, ann_path: Path) -> dict:
    n_cds = cds_count(species)
    if n_cds == 0:
        return {"Species": species, "Total CDS": 0}
    df = pd.read_csv(ann_path, sep="\t", comment="#", header=None, names=COLUMNS, dtype=str, keep_default_na=False)
    n = len(df)

    def pf(col: str) -> str:
        m = df[col].apply(has).sum()
        return f"{int(m)} ({100.0*m/n_cds:.1f}%)"

    return {
        "Species": species,
        "Total CDS": n_cds,
        "EggNOG annotated": f"{n} ({100.0*n/n_cds:.1f}%)",
        "COG cat": pf("COG_category"),
        "Description": pf("Description"),
        "Preferred name": pf("Preferred_name"),
        "GO": pf("GOs"),
        "KEGG_ko": pf("KEGG_ko"),
        "KEGG_Pathway": pf("KEGG_Pathway"),
        "PFAMs": pf("PFAMs"),
        "CAZy": pf("CAZy"),
        "EC": pf("EC"),
    }


def main() -> None:
    rows = []
    for sp in SPECIES:
        if sp in FOCAL_FILES and FOCAL_FILES[sp].is_file():
            ann = FOCAL_FILES[sp]
        else:
            ann = PER_SPECIES_DIR / f"{sp}.emapper.annotations"
        if not ann.is_file():
            print(f"[skip] {sp}: missing {ann}")
            continue
        print(f"[parse] {sp}: {ann}")
        rows.append(coverage_row(sp, ann))

    df = pd.DataFrame(rows).set_index("Species")
    print()
    print(df.to_string())

    for s, paths in SETS.items():
        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        out = paths["tables_dir"] / "annotation_coverage_all_species.csv"
        df.to_csv(out)
        print(f"[{s}] wrote {out}")


if __name__ == "__main__":
    main()
