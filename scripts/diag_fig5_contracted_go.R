#!/usr/bin/env Rscript
# Diagnostic for Figure 5: why is the dm (set3) CAFE-contracted GO enrichment empty?
# Replicates the notebook's contracted_OG computation, measures gene->GID mapping
# coverage, and runs enrichGO at q=0.05 and q=0.10 for (a) all-species genes (as
# the notebook does) and (b) dm-only genes. Prints counts only — writes nothing.
suppressMessages({
  library(dplyr); library(clusterProfiler); library(org.Dmaculatus.eg.db)
})
setwd("/home/jovyan")
source("/home/jovyan/Notebooks/config/load_config.R")
load_paths(config_path = "/home/jovyan/Notebooks/config/paths-set2.yaml", force_reload = TRUE)

ORTHOGROUPS_TSV   <- get_path("orthofinder.orthogroups_tsv")
CAFE_FAMILY_RES   <- get_path("cafe.single_lambda.family_results")
CAFE_CHANGE_TAB   <- get_path("cafe.single_lambda.change_tab")
focal             <- get_focal_species()
orgdb             <- org.Dmaculatus.eg.db

fam <- read.table(CAFE_FAMILY_RES, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
colnames(fam)[1] <- "FamilyID"
chg <- read.table(CAFE_CHANGE_TAB, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
asp_col <- grep(focal, colnames(chg), value = TRUE)[1]
cat("focal change column:", asp_col, "\n")

merged <- merge(fam, chg[, c("FamilyID", asp_col)], by = "FamilyID")
contracted_OG <- merged$FamilyID[merged$pvalue < 0.05 & merged[[asp_col]] < 0]
cat("significant contracted OGs (p<0.05 & change<0):", length(contracted_OG), "\n")

og <- read.delim(ORTHOGROUPS_TSV, header = TRUE, check.names = FALSE)
rows <- og[og$Orthogroup %in% contracted_OG, ]

collect <- function(df, cols) {
  g <- c()
  for (cn in cols) g <- c(g, unlist(strsplit(as.character(df[[cn]]), ", ")))
  unique(g[g != "" & !is.na(g)])
}
all_cols <- setdiff(colnames(rows), "Orthogroup")
genes_all <- collect(rows, all_cols)
genes_dm  <- collect(rows, focal)

gid_keys <- keys(orgdb, keytype = "GID")
cat(sprintf("genes (all species) = %d ; map to GID = %d (%.1f%%)\n",
            length(genes_all), sum(genes_all %in% gid_keys),
            100 * mean(genes_all %in% gid_keys)))
cat(sprintf("genes (dm only)      = %d ; map to GID = %d (%.1f%%)\n",
            length(genes_dm), sum(genes_dm %in% gid_keys),
            100 * mean(genes_dm %in% gid_keys)))

run <- function(genes, q, label) {
  res <- tryCatch(enrichGO(gene = genes, OrgDb = orgdb, keyType = "GID",
                           ont = "BP", pAdjustMethod = "BH",
                           pvalueCutoff = 0.05, qvalueCutoff = q),
                  error = function(e) {cat("  ERROR:", conditionMessage(e), "\n"); NULL})
  n <- if (is.null(res)) 0 else nrow(res)
  cat(sprintf("  [%s] enrichGO BP q<=%.2f -> %d terms\n", label, q, n))
}
cat("--- all-species genes ---\n"); run(genes_all, 0.05, "all"); run(genes_all, 0.10, "all")
cat("--- dm-only genes ---\n");     run(genes_dm, 0.05, "dm");   run(genes_dm, 0.10, "dm")
cat("DONE\n")
