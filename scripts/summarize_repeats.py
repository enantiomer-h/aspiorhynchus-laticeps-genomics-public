#!/usr/bin/env python3
"""Tier 2B: Aggregate per-species RepeatMasker .tbl summaries into a
comparative TE-classification table (analog of reference paper Tables 5–6).

Reads from:
  - Outputs/RepeatMasker_chromosomal/<species>/<*.fna.tbl>  (this-session reruns)
  - Outputs/Preprocessing/RepeatMasker/<species>/*.fna.tbl  (legacy)

Per species, extracts the % of genome occupied by each TE class
(DNA / LINE / SINE / LTR / RC/Helitron / Other / Total).

Output: Notebooks/Tables_set{2,3}/te_classification.csv
"""
import re
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _set_paths import SPECIES, SETS  # noqa: E402

PROJECT_ROOT = Path("/home/jovyan")
SOURCES = [
    PROJECT_ROOT / "Outputs/RepeatMasker_chromosomal_vertebrata",
    PROJECT_ROOT / "Outputs/RepeatMasker_chromosomal",
    PROJECT_ROOT / "Outputs/Preprocessing/RepeatMasker",
]


PCT_RE = re.compile(r"([\d.]+)\s*%")


def parse_tbl(tbl_path: Path) -> dict[str, float]:
    """Extract {category: % of genome} from a RepeatMasker .tbl."""
    cats = {
        "Total length (bp)": None,
        "GC level": None,
        "Bases masked %": None,
        "Retroelements %": None,
        "  SINEs %": None,
        "  Penelope %": None,
        "  LINEs %": None,
        "  LTR elements %": None,
        "DNA elements %": None,
        "Rolling-circles %": None,
        "Unclassified %": None,
        "Total interspersed repeats %": None,
        "Small RNA %": None,
        "Satellites %": None,
        "Simple repeats %": None,
        "Low complexity %": None,
    }

    text = tbl_path.read_text(errors="ignore").splitlines()
    for line in text:
        ls = line.strip()
        if ls.startswith("total length:"):
            m = re.search(r"total length:\s+([\d,]+) bp", line)
            if m: cats["Total length (bp)"] = int(m.group(1).replace(",", ""))
        elif ls.startswith("GC level:"):
            m = re.search(r"GC level:\s+([\d.]+) %", line)
            if m: cats["GC level"] = float(m.group(1))
        elif ls.startswith("bases masked:"):
            m = PCT_RE.search(line)
            if m: cats["Bases masked %"] = float(m.group(1))
        elif "Retroelements" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["Retroelements %"] = float(m.group(1))
        elif "SINEs:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["  SINEs %"] = float(m.group(1))
        elif "Penelope" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["  Penelope %"] = float(m.group(1))
        elif "LINEs:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["  LINEs %"] = float(m.group(1))
        elif "LTR elements:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["  LTR elements %"] = float(m.group(1))
        elif ls.startswith("DNA transposons") or "DNA transposons:" in line:
            m = PCT_RE.search(line)
            if m: cats["DNA elements %"] = float(m.group(1))
        elif "Rolling-circles" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["Rolling-circles %"] = float(m.group(1))
        elif "Unclassified:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["Unclassified %"] = float(m.group(1))
        elif "Total interspersed repeats:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["Total interspersed repeats %"] = float(m.group(1))
        elif "Small RNA:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["Small RNA %"] = float(m.group(1))
        elif "Satellites:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["Satellites %"] = float(m.group(1))
        elif "Simple repeats:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["Simple repeats %"] = float(m.group(1))
        elif "Low complexity:" in line and "%" in line:
            m = PCT_RE.search(line)
            if m: cats["Low complexity %"] = float(m.group(1))
    return cats


def find_tbl(species: str) -> Path | None:
    for src in SOURCES:
        d = src / species
        if not d.is_dir(): continue
        tbls = list(d.glob("**/*.tbl"))
        tbls = [t for t in tbls if t.stat().st_size > 100]
        if tbls:
            return tbls[0]
    return None


def main() -> None:
    rows = []
    for sp in SPECIES:
        tbl = find_tbl(sp)
        if tbl is None:
            print(f"[skip] {sp}: no .tbl found")
            continue
        print(f"[parse] {sp}: {tbl}")
        d = parse_tbl(tbl)
        d["Species"] = sp
        rows.append(d)

    if not rows:
        print("No RepeatMasker .tbl files found yet.")
        return

    df = pd.DataFrame(rows).set_index("Species")
    print()
    print(df.to_string())

    for s, paths in SETS.items():
        paths["tables_dir"].mkdir(parents=True, exist_ok=True)
        out = paths["tables_dir"] / "te_classification.csv"
        df.to_csv(out)
        print(f"[{s}] wrote {out}")


if __name__ == "__main__":
    main()
