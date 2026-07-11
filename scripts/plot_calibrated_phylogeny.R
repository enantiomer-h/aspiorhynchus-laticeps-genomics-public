#!/usr/bin/env Rscript
# Figure 1 — calibrated phylogeny composite (PI revisions, set3 P14-16).
#
# Keep (per PI): species relationships, RED node ages, bottom geological epochs,
# key internal-node DUPLICATION counts ("Dup: 83 / Dup: 1459").
# Fixes:
#  (1) node-age labels: the notebook computed labels as `root_age - x` AFTER
#      revts() flips x to negative, inflating every label by ~root_age and making
#      the numbers inconsistent with the geological epoch band. Here ages are
#      precomputed (root_age - node_depth) and joined by node id — no revts math.
#  (2) per-node duplication counts: parsed from the OrthoFinder duplications tree
#      node labels (N0_83, N1_1459, ...) and mapped onto the chronogram nodes by
#      matching descendant-tip sets, then drawn as "Dup: N".
#
# Inputs:
#   Outputs/.../Species_Tree/SpeciesTree_chronogram_MYA.tre          (ultrametric)
#   Outputs/.../Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support.txt
# Output (both set dirs): SpeciesTree_chronogram_geological_ggtree.{png,pdf,svg}
suppressMessages({
  library(ape); library(ggtree); library(ggplot2); library(deeptime); library(treeio)
})

OF        <- "/home/jovyan/Outputs/OrthoFinder/Results_Mar03"
CHRONO    <- file.path(OF, "Species_Tree/SpeciesTree_chronogram_MYA.tre")
DUPTREE   <- file.path(OF, "Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support.txt")
OUT_DIRS  <- c("/home/jovyan/Notebooks/Figures_set2", "/home/jovyan/Notebooks/Figures_set3")

chrono <- read.tree(CHRONO)
dup    <- read.tree(DUPTREE)
# supermatrix ML tree carries ultrafast-bootstrap (UFBoot, %) support on its nodes
super  <- read.tree("/home/jovyan/Notebooks/Figures_set2/supermatrix_iqtree.treefile")

# clade signature = sorted descendant tip labels, joined
clade_sig <- function(tree, node, strip = FALSE) {
  tips <- extract.clade(tree, node)$tip.label
  if (strip) tips <- sub("_\\d+$", "", tips)   # drop trailing _<count> on dup-tree tips
  paste(sort(tips), collapse = "|")
}

# --- per-node duplication counts from the dup tree's internal-node labels ------
ntip_d <- Ntip(dup)
dup_node_ids <- (ntip_d + 1):(ntip_d + dup$Nnode)
dup_counts   <- as.integer(sub("^N\\d+_", "", dup$node.label))   # "N1_1459" -> 1459
dup_by_clade <- setNames(dup_counts, sapply(dup_node_ids, clade_sig, tree = dup, strip = TRUE))

# --- UFBoot support from the supermatrix ML tree, matched by BIPARTITION -------
# The supermatrix tree's unrooted topology is identical to the chronogram's
# (Robinson-Foulds distance = 0); only the rooting differs. Support is a property
# of a branch (bipartition), so each chronogram branch maps to a supermatrix
# branch by matching its descendant set OR that set's complement. Values are %.
all_tips    <- sort(chrono$tip.label)
su_node_ids <- (Ntip(super) + 1):(Ntip(super) + super$Nnode)
su_by_clade <- setNames(super$node.label, sapply(su_node_ids, clade_sig, tree = super))
support_for <- function(nd) {
  tips <- extract.clade(chrono, nd)$tip.label
  v <- su_by_clade[paste(sort(tips), collapse = "|")]
  if (is.na(v)) v <- su_by_clade[paste(sort(setdiff(all_tips, tips)), collapse = "|")]
  if (is.na(v) || v == "") NA_character_ else v
}

# --- node ages + duplication counts mapped onto the chronogram -----------------
depths   <- node.depth.edgelength(chrono)
root_age <- max(depths)
ntip_c   <- Ntip(chrono)
nnode_c  <- chrono$Nnode
internal <- (ntip_c + 1):(ntip_c + nnode_c)

node_df <- data.frame(node = seq_len(ntip_c + nnode_c),
                      age  = round(root_age - depths, 1),
                      dup  = NA_integer_,
                      bs   = NA_character_, stringsAsFactors = FALSE)
for (nd in internal) {
  sig <- clade_sig(chrono, nd, strip = FALSE)
  if (!is.na(dup_by_clade[sig])) node_df$dup[node_df$node == nd] <- dup_by_clade[sig]
  node_df$bs[node_df$node == nd] <- support_for(nd)
}
node_df$age_lab <- ifelse(node_df$node %in% internal, sprintf("%.1f", node_df$age), NA)
node_df$dup_lab <- ifelse(!is.na(node_df$dup), paste0("Dup: ", node_df$dup), NA)
node_df$bs_lab  <- ifelse(node_df$node %in% internal & !is.na(node_df$bs), node_df$bs, NA)
cat("Root age:", round(root_age, 1), "Ma | dup counts at",
    sum(!is.na(node_df$dup)), "/", nnode_c, "nodes | UFBoot support at",
    sum(!is.na(node_df$bs)), "/", nnode_c, "nodes\n")
print(node_df[node_df$node %in% internal, c("node", "age", "dup", "bs")], row.names = FALSE)

# --- build the composite -------------------------------------------------------
root_age_rounded <- ceiling(root_age / 10) * 10
axis_breaks <- seq(0, root_age_rounded, by = 20)

p <- ggtree(chrono, size = 1.6, color = "#2E4057") %<+% node_df +
  geom_tiplab(aes(label = gsub("_", " ", label)), fontface = "italic",
              size = 7, offset = root_age * 0.015, color = "#1A1A1A") +
  geom_nodepoint(color = "#4A90D9", size = 4, alpha = 0.85) +
  # red node ages ABOVE each internal node (correct values, joined by node id)
  geom_text2(aes(subset = !isTip, label = age_lab),
             size = 5.2, color = "#CC4444", fontface = "bold", nudge_y = 0.34) +
  # blue duplication counts BELOW each internal node
  geom_text2(aes(subset = !isTip, label = dup_lab),
             size = 5.0, color = "#1B5E9B", fontface = "bold", nudge_y = -0.34) +
  # green ultrafast-bootstrap support on the branch entering each internal node
  geom_text2(aes(subset = !isTip, label = bs_lab),
             size = 4.3, color = "#2E7D32", fontface = "bold",
             nudge_x = -root_age * 0.022, nudge_y = 0.16, hjust = 1) +
  theme_tree2()

p_geo <- revts(p) +
  coord_geo(xlim = c(-(root_age_rounded + 6), root_age * 0.34),
            ylim = c(-2, ntip_c + 1.5), neg = TRUE,
            abbrv = list(TRUE, FALSE), pos = list("bottom", "bottom"),
            dat = list("epochs", "periods"), size = list(3.0, 4.5),
            center_end_labels = TRUE) +
  scale_x_continuous(name = "Divergence time (Million years ago)",
                     breaks = -rev(axis_breaks), labels = rev(axis_breaks)) +
  theme(axis.title.x = element_text(size = 20, face = "bold", family = "serif"),
        axis.text.x  = element_text(size = 18, color = "#333333"),
        plot.title   = element_text(size = 24, face = "bold", hjust = 0.5, family = "serif"),
        plot.subtitle = element_text(size = 15, hjust = 0.5, color = "#666666", family = "serif"),
        plot.margin = margin(10, 40, 10, 10)) +
  labs(title = "Time-calibrated phylogeny with gene-duplication events",
       subtitle = paste0("Red = node age (Ma); green = supermatrix ML ultrafast-bootstrap support (%); ",
                         "blue = inferred gene duplications; root age ", round(root_age, 1), " Ma"))

for (d in OUT_DIRS) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(d, "SpeciesTree_chronogram_geological_ggtree")
  ggsave(paste0(base, ".png"), p_geo, width = 20, height = 14, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), p_geo, width = 20, height = 14, bg = "white", device = "pdf")
  ggsave(paste0(base, ".svg"), p_geo, width = 20, height = 14, bg = "white", device = svglite::svglite)
  cat("wrote", base, ".{png,pdf,svg}\n")
}
cat("DONE\n")
