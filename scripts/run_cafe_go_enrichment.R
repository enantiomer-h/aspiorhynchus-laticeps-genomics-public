#!/usr/bin/env Rscript
# CAFE-based GO enrichment for the focal species, with the Figure-5 bug fix.
#
# BUG (diagnosed): the notebook's enrich_orthogroups() pooled genes from ALL 12
# species and mapped them to the focal species' OrgDb. Because BRAKER gene IDs
# (g<n>.t<n>) are reused across species, only ~24% of the 88k pooled genes hit
# the dm OrgDb (mostly coincidental collisions) -> enrichGO returned 0 terms for
# dm contracted families. FIX: extract ONLY the focal species' genes from the
# target orthogroups (90% map) -> 354 BP terms at FDR<0.05, no cutoff change.
#
# Usage:
#   Rscript run_cafe_go_enrichment.R [config.yaml] [Category1,Category2,...]
#     config.yaml : default /home/jovyan/Notebooks/config/paths-set2.yaml
#     Categories  : CAFE_Contracted (default), CAFE_Expanded, unique_nonzero
#
# Outputs editable PNG+PDF+SVG barplot & dotplot into
#   <figures_dir>/GO_Enrichment/<Category>_<ONT>_Plots/
# and the result table into <tables_dir>/GO_Enrichment/.
suppressMessages({
  library(dplyr); library(clusterProfiler); library(ggplot2); library(enrichplot)
})

args        <- commandArgs(trailingOnly = TRUE)
default_cfg <- Sys.getenv("GENOMICS_CONFIG_PATH", "/home/jovyan/Notebooks/config/paths-set2.yaml")
config_path <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_cfg
categories  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else c("CAFE_Contracted")
ONTS        <- c("BP", "MF", "CC")
P_CUT <- 0.05; Q_CUT <- 0.05; Q_RELAX <- 0.10; SHOW <- 15

setwd("/home/jovyan")
source("/home/jovyan/Notebooks/config/load_config.R")
suppressMessages(load_paths(config_path = config_path, force_reload = TRUE))

focal       <- get_focal_species()
genus       <- get_focal_genus(); epithet <- get_focal_species_epithet()
focal_disp  <- paste0(substr(genus, 1, 1), ". ", epithet)   # e.g. "D. maculatus"
orgdb_name  <- get_orgdb_package_name()
suppressMessages(library(orgdb_name, character.only = TRUE))
orgdb       <- get(orgdb_name)

FIG     <- get_path("outputs.figures")
TBL     <- get_path("outputs.tables")
GO_TBL  <- file.path(TBL, "GO_Enrichment"); dir.create(GO_TBL, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("config=%s | focal=%s | orgdb=%s\n", config_path, focal, orgdb_name))

# ---- Orthogroup category sets -------------------------------------------------
fam <- read.table(get_path("cafe.single_lambda.family_results"), header = TRUE,
                  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
colnames(fam)[1] <- "FamilyID"
chg <- read.table(get_path("cafe.single_lambda.change_tab"), header = TRUE,
                  sep = "\t", stringsAsFactors = FALSE)
asp_col <- grep(focal, colnames(chg), value = TRUE)[1]
merged  <- merge(fam, chg[, c("FamilyID", asp_col)], by = "FamilyID")

og_tsv   <- read.delim(get_path("orthofinder.orthogroups_tsv"), header = TRUE, check.names = FALSE)
gc_count <- read.table(get_path("orthofinder.orthogroups_count"), header = TRUE,
                       sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

og_sets <- list(
  CAFE_Contracted = merged$FamilyID[merged$pvalue < 0.05 & merged[[asp_col]] < 0],
  CAFE_Expanded   = merged$FamilyID[merged$pvalue < 0.05 & merged[[asp_col]] > 0]
)
# species-specific: focal > 0 and all other species == 0
sp_cols <- setdiff(colnames(gc_count), c("Orthogroup", "Total"))
others  <- setdiff(sp_cols, focal)
ss_mask <- gc_count[[focal]] > 0 & rowSums(gc_count[, others] != 0) == 0
og_sets$unique_nonzero <- gc_count$Orthogroup[ss_mask]

focal_genes_of <- function(ogs) {
  rows <- og_tsv[og_tsv$Orthogroup %in% ogs, focal]
  g <- unlist(strsplit(as.character(rows), ", "))
  unique(g[g != "" & !is.na(g)])
}

# ---- editable, large-font plotting -------------------------------------------
pub_theme <- theme_minimal(base_size = 16) + theme(
  text         = element_text(size = 16),
  plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
  axis.title   = element_text(size = 16, face = "bold"),
  axis.text    = element_text(size = 14, color = "black"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text  = element_text(size = 13),
  panel.grid.minor = element_blank()
)

save_editable <- function(plt, outdir, basename, w = 12, h = 9) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, paste0(basename, ".png")), plt, width = w, height = h, dpi = 300)
  ggsave(file.path(outdir, paste0(basename, ".pdf")), plt, width = w, height = h, device = "pdf")  # editable text
  ggsave(file.path(outdir, paste0(basename, ".svg")), plt, width = w, height = h, device = svglite::svglite)
}

enrich_one <- function(category, ont) {
  ogs   <- og_sets[[category]]
  genes <- focal_genes_of(ogs)
  cat(sprintf("[%s/%s] OGs=%d focal-genes=%d\n", category, ont, length(ogs), length(genes)))
  ego <- tryCatch(enrichGO(gene = genes, OrgDb = orgdb, keyType = "GID", ont = ont,
                           pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT),
                  error = function(e) {cat("  enrichGO error:", conditionMessage(e), "\n"); NULL})
  qused <- Q_CUT
  if (!is.null(ego) && nrow(ego) == 0) {     # relax once, per project decision
    ego <- enrichGO(gene = genes, OrgDb = orgdb, keyType = "GID", ont = ont,
                    pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_RELAX)
    qused <- Q_RELAX
  }
  write.csv(as.data.frame(ego),
            file.path(GO_TBL, sprintf("GO_enrichment_%s_%s_results.csv", category, ont)),
            row.names = FALSE)
  if (is.null(ego) || nrow(ego) == 0) {
    cat(sprintf("  -> no significant terms (even at q<=%.2f); no plot\n", Q_RELAX)); return(invisible(0))
  }
  cat(sprintf("  -> %d terms (FDR<=%.2f)\n", nrow(ego), qused))
  cat_label <- gsub("_", "-", category)
  ttl <- function(kind) bquote(italic(.(focal_disp))~.(sprintf(": %s families — GO %s (%s)", cat_label, ont, kind)))
  outdir <- file.path(FIG, "GO_Enrichment", sprintf("%s_%s_Plots", category, ont))
  bp <- barplot(ego, showCategory = SHOW, label_format = 45) + pub_theme + ggtitle(ttl("bar"))
  dp <- dotplot(ego, showCategory = SHOW, label_format = 45) + pub_theme + ggtitle(ttl("dot"))
  save_editable(bp, outdir, "GO_enrichment_barplot")
  save_editable(dp, outdir, "GO_enrichment_dotplot")
  invisible(nrow(ego))
}

for (category in categories) for (ont in ONTS) enrich_one(category, ont)
cat("DONE\n")
