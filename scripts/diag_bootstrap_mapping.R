#!/usr/bin/env Rscript
# Diagnostic: can the supermatrix ML bootstrap values be mapped onto the
# chronogram (OrthoFinder species tree) nodes? Only clades present in BOTH trees
# can receive a support value. Reports per-chronogram-node match + RF distance.
suppressMessages(library(ape))

OF      <- "/home/jovyan/Outputs/OrthoFinder/Results_Mar03"
chrono  <- read.tree(file.path(OF, "Species_Tree/SpeciesTree_chronogram_MYA.tre"))
super   <- read.tree("/home/jovyan/Notebooks/Figures_set2/supermatrix_iqtree.treefile")

cat("same tip set:", setequal(chrono$tip.label, super$tip.label), "\n")
cat("RF distance (unrooted):", dist.topo(unroot(chrono), unroot(super)), "\n")
cat("(RF = 0 means identical unrooted topology)\n\n")

clade_sig <- function(tree, node) paste(sort(extract.clade(tree, node)$tip.label), collapse = "|")

# supermatrix: clade signature -> support label
su_int <- (Ntip(super) + 1):(Ntip(super) + super$Nnode)
su_map <- setNames(super$node.label, sapply(su_int, clade_sig, tree = super))

ch_int <- (Ntip(chrono) + 1):(Ntip(chrono) + chrono$Nnode)
cat(sprintf("%-4s %-9s %s\n", "node", "support", "clade (descendant tips)"))
for (nd in ch_int) {
  sig <- clade_sig(chrono, nd)
  sup <- if (sig %in% names(su_map)) su_map[[sig]] else "—— CONFLICT"
  tips <- strsplit(sig, "\\|")[[1]]
  short <- if (length(tips) <= 4) paste(tips, collapse = ", ")
           else paste0(length(tips), " spp: ", paste(substr(tips,1,12), collapse=","))
  cat(sprintf("%-4d %-9s %s\n", nd, sup, short))
}
matched <- sum(sapply(ch_int, function(nd) clade_sig(chrono, nd) %in% names(su_map)))
cat(sprintf("\nmatched %d / %d internal nodes\n", matched, length(ch_int)))
