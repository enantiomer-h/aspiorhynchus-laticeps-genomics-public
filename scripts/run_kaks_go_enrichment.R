#!/usr/bin/env Rscript
# Figure 6 (right panel) — GO enrichment of genes under selection (KaKs categories).
# Same focal-only gene fix + editable/large-font output as run_cafe_go_enrichment.R.
#
# KaKs categories from the genome-wide selection screen
# (Outputs_set{2,3}/GenomewideSelection/summary/genome_wide_selection_master.tsv):
#   positive_tentative : max_omega > 1            (PI: al positively selected)
#   strong_purifying   : mean_omega <= 0.1        (PI: dm strong purifying)
#
# Usage: Rscript run_kaks_go_enrichment.R [config.yaml] [cat1,cat2]
#   Categories: positive_tentative (default), strong_purifying
suppressMessages({
  library(dplyr); library(data.table); library(clusterProfiler); library(ggplot2); library(enrichplot)
})

args        <- commandArgs(trailingOnly = TRUE)
default_cfg <- Sys.getenv("GENOMICS_CONFIG_PATH", "/home/jovyan/Notebooks/config/paths-set2.yaml")
config_path <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_cfg
categories  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else c("positive_tentative", "strong_purifying")
ONTS <- c("BP", "MF", "CC"); P_CUT <- 0.05; Q_CUT <- 0.05; Q_RELAX <- 0.10; SHOW <- 15

setwd("/home/jovyan")
source("/home/jovyan/Notebooks/config/load_config.R")
suppressMessages(load_paths(config_path = config_path, force_reload = TRUE))

focal      <- get_focal_species()
genus      <- get_focal_genus(); epithet <- get_focal_species_epithet()
focal_disp <- paste0(substr(genus, 1, 1), ". ", epithet)
orgdb_name <- get_orgdb_package_name()
suppressMessages(library(orgdb_name, character.only = TRUE)); orgdb <- get(orgdb_name)
FIG <- get_path("outputs.figures"); TBL <- get_path("outputs.tables")
GO_TBL <- file.path(TBL, "GO_Enrichment"); dir.create(GO_TBL, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("config=%s | focal=%s | orgdb=%s\n", config_path, focal, orgdb_name))

# --- KaKs OG categories from the genome-wide selection master -----------------
master_path <- file.path(get_path("outputs.genomewide_selection.summary"), "genome_wide_selection_master.tsv")
master <- fread(master_path, select = c("Orthogroup", "mean_omega", "max_omega"))
og_sets <- list(
  positive_tentative = master$Orthogroup[!is.na(master$max_omega)  & master$max_omega > 1],
  strong_purifying   = master$Orthogroup[!is.na(master$mean_omega) & master$mean_omega <= 0.1]
)
# label used for the output directory name (matches existing KaKs_*_Plots convention)
DIR_LABEL <- c(positive_tentative = "KaKs_Positive_Tentative",
               strong_purifying   = "KaKs_Strong_Purifying")

og_tsv <- read.delim(get_path("orthofinder.orthogroups_tsv"), header = TRUE, check.names = FALSE)
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

enrich_one <- function(category, ont) {
  genes <- focal_genes_of(og_sets[[category]])
  cat(sprintf("[%s/%s] OGs=%d focal-genes=%d\n", category, ont, length(og_sets[[category]]), length(genes)))
  ego <- tryCatch(enrichGO(gene = genes, OrgDb = orgdb, keyType = "GID", ont = ont,
                           pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT),
                  error = function(e) {cat("  enrichGO error:", conditionMessage(e), "\n"); NULL})
  qused <- Q_CUT
  if (!is.null(ego) && nrow(ego) == 0) {
    ego <- enrichGO(gene = genes, OrgDb = orgdb, keyType = "GID", ont = ont,
                    pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_RELAX); qused <- Q_RELAX
  }
  write.csv(as.data.frame(ego),
            file.path(GO_TBL, sprintf("GO_enrichment_%s_%s_results.csv", category, ont)), row.names = FALSE)
  if (is.null(ego) || nrow(ego) == 0) { cat("  -> no significant terms; no plot\n"); return(invisible(0)) }
  cat(sprintf("  -> %d terms (FDR<=%.2f); top: %s\n", nrow(ego), qused,
              paste(head(as.data.frame(ego)$Description, 4), collapse = " | ")))
  outdir <- file.path(FIG, "GO_Enrichment", sprintf("%s_%s_Plots", DIR_LABEL[[category]], ont))
  sel_lab <- ifelse(category == "positive_tentative", "positive selection", "strong purifying")
  ttl <- function(kind) bquote(italic(.(focal_disp))~.(sprintf(": %s — GO %s (%s)", sel_lab, ont, kind)))
  bp <- barplot(ego, showCategory = SHOW, label_format = 45) + pub_theme + ggtitle(ttl("bar"))
  dp <- dotplot(ego, showCategory = SHOW, label_format = 45) + pub_theme + ggtitle(ttl("dot"))
  save_editable(bp, outdir, "GO_enrichment_barplot")
  save_editable(dp, outdir, "GO_enrichment_dotplot")
  invisible(nrow(ego))
}

for (category in categories) for (ont in ONTS) enrich_one(category, ont)
cat("DONE\n")
