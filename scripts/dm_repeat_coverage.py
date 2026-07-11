"""Quick: DM non-overlapping repeat coverage from .fa.out.gff."""
from collections import defaultdict
from pathlib import Path

gff = Path("/home/jovyan/DB/GENOME_Comparison/Diptychus_maculatus/DM.fa.out.gff")
ivs = defaultdict(list)
with gff.open() as fh:
    for line in fh:
        if line.startswith("#") or not line.strip(): continue
        f = line.split("\t")
        if len(f) < 5: continue
        try:
            ivs[f[0]].append((int(f[3]), int(f[4])))
        except ValueError:
            continue

total = 0
for ch, lst in ivs.items():
    lst.sort()
    cur_s, cur_e = lst[0]
    for s, e in lst[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    total += cur_e - cur_s + 1

genome = 1506815489
print(f"DM non-overlapping repeat bp = {total:,}")
print(f"DM genome size              = {genome:,} bp")
print(f"DM masked %                 = {100.0*total/genome:.2f}%")
