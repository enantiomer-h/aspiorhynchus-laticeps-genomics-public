#!/usr/bin/env Rscript
# CAFE-based KEGG (KO) enrichment for the focal species, with the same
# focal-only gene fix as the GO script (run_cafe_go_enrichment.R). Produces
# editable PNG+PDF+SVG barplot & dotplot with large fonts.
#
# Usage:
#   Rscript run_cafe_ko_enrichment.R [config.yaml] [Category1,Category2,...]
#     Categories: CAFE_Expanded (default), CAFE_Contracted, unique_nonzero
#   NOTE: do NOT regenerate set3 CAFE_Contracted KO — that figure is PI-approved.
suppressMessages({
  library(dplyr); library(clusterProfiler); library(ggplot2); library(enrichplot)
})

args        <- commandArgs(trailingOnly = TRUE)
default_cfg <- Sys.getenv("GENOMICS_CONFIG_PATH", "/home/jovyan/Notebooks/config/paths-set2.yaml")
config_path <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_cfg
categories  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else c("CAFE_Expanded")
P_CUT <- 0.05; Q_CUT <- 0.05; Q_RELAX <- 0.10; SHOW <- 15

setwd("/home/jovyan")
source("/home/jovyan/Notebooks/config/load_config.R")
suppressMessages(load_paths(config_path = config_path, force_reload = TRUE))

focal      <- get_focal_species()
genus      <- get_focal_genus(); epithet <- get_focal_species_epithet()
focal_disp <- paste0(substr(genus, 1, 1), ". ", epithet)
FIG        <- get_path("outputs.figures")
TBL        <- get_path("outputs.tables")
KO_TBL     <- file.path(TBL, "KO_Enrichment"); dir.create(KO_TBL, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("config=%s | focal=%s\n", config_path, focal))

# ---- KEGG term2gene / term2name ----------------------------------------------
g2p <- read.table(get_path("eggnog.gene2pathway"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
term2gene <- g2p[, c("Pathway", "GID")]
load(get_path("eggnog.kegg_info"))   # provides pathway2name
pathway2name$Name <- gsub("\\s*\\[BR:ko[0-9]+\\]", "", pathway2name$Name)

# ---- OG category sets (same as GO script) ------------------------------------
fam <- read.table(get_path("cafe.single_lambda.family_results"), header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
colnames(fam)[1] <- "FamilyID"
chg <- read.table(get_path("cafe.single_lambda.change_tab"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
asp_col <- grep(focal, colnames(chg), value = TRUE)[1]
merged  <- merge(fam, chg[, c("FamilyID", asp_col)], by = "FamilyID")
og_tsv  <- read.delim(get_path("orthofinder.orthogroups_tsv"), header = TRUE, check.names = FALSE)
gc_count <- read.table(get_path("orthofinder.orthogroups_count"), header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
sp_cols <- setdiff(colnames(gc_count), c("Orthogroup", "Total")); others <- setdiff(sp_cols, focal)
og_sets <- list(
  CAFE_Contracted = merged$FamilyID[merged$pvalue < 0.05 & merged[[asp_col]] < 0],
  CAFE_Expanded   = merged$FamilyID[merged$pvalue < 0.05 & merged[[asp_col]] > 0],
  unique_nonzero  = gc_count$Orthogroup[gc_count[[focal]] > 0 & rowSums(gc_count[, others] != 0) == 0]
)
focal_genes_of <- function(ogs) {
  g <- unlist(strsplit(as.character(og_tsv[og_tsv$Orthogroup %in% ogs, focal]), ", "))
  unique(g[g != "" & !is.na(g)])
}

pub_theme <- theme_minimal(base_size = 16) + theme(
  text = element_text(size = 16), plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 14, color = "black"),
  legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 13),
  panel.grid.minor = element_blank())
save_editable <- function(plt, outdir, basename, w = 12, h = 9) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, paste0(basename, ".png")), plt, width = w, height = h, dpi = 300)
  ggsave(file.path(outdir, paste0(basename, ".pdf")), plt, width = w, height = h, device = "pdf")
  ggsave(file.path(outdir, paste0(basename, ".svg")), plt, width = w, height = h, device = svglite::svglite)
}

enrich_one <- function(category) {
  genes <- focal_genes_of(og_sets[[category]])
  cat(sprintf("[%s/KO] OGs=%d focal-genes=%d\n", category, length(og_sets[[category]]), length(genes)))
  res <- tryCatch(enricher(gene = genes, TERM2GENE = term2gene, TERM2NAME = pathway2name,
                           pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, pAdjustMethod = "BH", minGSSize = 1),
                  error = function(e) {cat("  enricher error:", conditionMessage(e), "\n"); NULL})
  qused <- Q_CUT
  if (!is.null(res) && nrow(res) == 0) {
    res <- enricher(gene = genes, TERM2GENE = term2gene, TERM2NAME = pathway2name,
                    pvalueCutoff = P_CUT, qvalueCutoff = Q_RELAX, pAdjustMethod = "BH", minGSSize = 1)
    qused <- Q_RELAX
  }
  write.table(as.data.frame(res),
              file.path(KO_TBL, sprintf("KO_enrichment_%s_results.tsv", category)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  if (is.null(res) || nrow(res) == 0) { cat("  -> no significant KO pathways; no plot\n"); return(invisible(0)) }
  cat(sprintf("  -> %d pathways (FDR<=%.2f)\n", nrow(res), qused))
  cat_label <- gsub("_", "-", category)
  outdir <- file.path(FIG, "KO_Enrichment", sprintf("%s_KO_Plots", category))
  bp <- barplot(res, showCategory = SHOW, label_format = 45) + pub_theme +
    ggtitle(bquote(italic(.(focal_disp))~.(sprintf(": %s families — KEGG (bar)", cat_label))))
  dp <- dotplot(res, showCategory = SHOW, label_format = 45) + pub_theme +
    ggtitle(bquote(italic(.(focal_disp))~.(sprintf(": %s families — KEGG (dot)", cat_label))))
  save_editable(bp, outdir, "KO_enrichment_barplot")
  save_editable(dp, outdir, "KO_enrichment_dotplot")
  invisible(nrow(res))
}

for (category in categories) enrich_one(category)
cat("DONE\n")
