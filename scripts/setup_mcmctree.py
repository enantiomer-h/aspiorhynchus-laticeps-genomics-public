#!/usr/bin/env python3
"""Tier 3B: Set up MCMCTree divergence dating from the 3A IQ-TREE supermatrix.

Two-step PAML MCMCTree workflow with approximate likelihood (the standard for
large AA matrices):

  Step 1: mcmctree --usedata=3 to compute branch-length Hessian (`in.BV`).
  Step 2: mcmctree --usedata=2 to run the MCMC posterior using `in.BV`.

Both control files + the calibrated rooted tree + the PHYLIP-format alignment
are written under Outputs/Phylogenomics_supermatrix/mcmctree/.

Calibrations (in units of 100 Mya; soft bounds):
  - Root (Cyprinidae crown, separating Danio-clade from rest):  49-150 Mya
    -> 'B(0.49, 1.5, 0.025, 0.025)'
  - Cyprininae (Carassius+Cyprinus+Sinocyclocheilus):           28-50 Mya
    -> 'B(0.28, 0.50, 0.025, 0.025)'
  - Schizothoracinae (Aspiorhynchus + 4 schizothoracines):      5-25 Mya
    -> 'B(0.05, 0.25, 0.025, 0.025)'

Topology rooted on (Danio + Triplophysa-clade) | (everything else).
"""
import sys
from pathlib import Path
from collections import OrderedDict

PROJECT_ROOT = Path("/home/jovyan")
SUP_FAA = PROJECT_ROOT / "Outputs/Phylogenomics_supermatrix/supermatrix.faa"
WORK = PROJECT_ROOT / "Outputs/Phylogenomics_supermatrix/mcmctree"
WORK.mkdir(parents=True, exist_ok=True)

PHYLIP_OUT = WORK / "supermatrix.phy"
TREE_OUT = WORK / "calibrated.tre"
CTL_STEP1 = WORK / "mcmctree_step1.ctl"
CTL_STEP2 = WORK / "mcmctree_step2.ctl"


def faa_to_phylip(faa: Path, phy: Path) -> tuple[int, int]:
    seqs = OrderedDict()
    cur_name, cur = None, []
    with faa.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_name is not None:
                    seqs[cur_name] = "".join(cur)
                cur_name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line)
        if cur_name is not None:
            seqs[cur_name] = "".join(cur)
    n = len(seqs)
    L = len(next(iter(seqs.values())))
    assert all(len(s) == L for s in seqs.values()), "unequal seq lengths"
    with phy.open("w") as out:
        out.write(f"  {n}  {L}\n")
        for name, seq in seqs.items():
            short = name[:30]
            out.write(f"{short.ljust(31)}{seq}\n")
    return n, L


CALIBRATED_TREE = """\
12 1
((Danio_rerio,(Triplophysa_pappenheimi,(Triplophysa_tibetana,Triplophysa_yaopeizhii))),(((Carassius_auratus,Cyprinus_carpio),Sinocyclocheilus_grahami)'B(0.28,0.50,0.025,0.025)',(Aspiorhynchus_laticeps,(Diptychus_maculatus,((Gymnocypris_eckloni,Schizopygopsis_younghusbandi),Oxygymnocypris_stewartii)))'B(0.05,0.25,0.025,0.025)'))'B(0.49,1.50,0.025,0.025)';
"""


def write_tree() -> None:
    TREE_OUT.write_text(CALIBRATED_TREE)


def write_ctl(usedata: int, ctl_path: Path, alignment: Path, tree: Path, mcmc_label: str) -> None:
    body = f"""          seed = -1
       seqfile = {alignment}
      treefile = {tree}
       mcmcfile = {WORK}/mcmc_{mcmc_label}.txt
       outfile  = {WORK}/out_{mcmc_label}.txt

         ndata = 1
       seqtype = 2     * 2 = aa
       usedata = {usedata}     * 0 prior; 1 exact lik; 2 approx; 3 outBV
        clock = 2     * 1 strict, 2 ind. rates, 3 corr. rates
       RootAge = 'B(0.49,1.50,0.025,0.025)'

       model = 2     * JTT+Gamma
       alpha = 0.5
        ncatG = 5
     cleandata = 0

   BDparas = 1 1 0.1 C
   kappa_gamma = 6 2
   alpha_gamma = 1 1
   rgene_gamma = 2 20 1
   sigma2_gamma = 1 10 1

   finetune = 1: .1 .1 .1 .1 .1 .1
       print = 1
      burnin = 100000
    sampfreq = 100
     nsample = 50000
"""
    ctl_path.write_text(body)


def main() -> None:
    n, L = faa_to_phylip(SUP_FAA, PHYLIP_OUT)
    print(f"Wrote PHYLIP: {n} species × {L:,} aa columns to {PHYLIP_OUT}")
    write_tree()
    print(f"Wrote calibrated rooted tree to {TREE_OUT}")
    write_ctl(3, CTL_STEP1, PHYLIP_OUT, TREE_OUT, "step1_hessian")
    write_ctl(2, CTL_STEP2, PHYLIP_OUT, TREE_OUT, "step2_posterior")
    print(f"Wrote {CTL_STEP1}")
    print(f"Wrote {CTL_STEP2}")
    print()
    print("Run step 1 (Hessian, ~30-60 min for 226K aa cols):")
    print(f"  cd {WORK} && conda run -n env_paml mcmctree {CTL_STEP1.name}")
    print("After step 1 produces 'rst2', rename it to 'in.BV' for step 2:")
    print(f"  mv rst2 in.BV")
    print("Then run step 2 (MCMC posterior, 1-3 h):")
    print(f"  cd {WORK} && conda run -n env_paml mcmctree {CTL_STEP2.name}")


if __name__ == "__main__":
    main()
