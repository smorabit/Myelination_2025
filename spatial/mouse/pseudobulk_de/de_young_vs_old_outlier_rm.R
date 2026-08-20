#!/usr/bin/env Rscript
# DE analysis: Young vs Old per oligo cell type, with the OUTLIER sample
# (Old 2, the canonical outlier) excluded.
#
# Why: Old 2 has the smallest pseudobulk in every cell type under P5
# (OPC=14, Inter=5, Mature=10 lesion cells; total raw counts 393–964
# across 50 genes), and it sits visibly outlying on PC1×PC2 of the
# combined-VST PCA. This script rebuilds the DE estimates with that
# sample excluded as a sensitivity check.
#
# Outputs are written next to the standard outputs but suffixed
# `_outlier_rm` (file_basename + "_outlier_rm" + extension), e.g.:
#   data/de_results/<strategy>/OPC_de_outlier_rm.csv
#   data/de_results/<strategy>/_combined_de_outlier_rm.csv
#
# Usage: Rscript scripts/de_young_vs_old_outlier_rm.R <strategy>
#   strategy: stringent_p5 (default) | stringent | stringent_p4 | ...
#
# Design after exclusion: 5 Young + 4 Old (n = 9 samples per cell type).

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
strategy <- if (length(args) >= 1) args[1] else "stringent_p5"

OUTLIER_TAG <- "outlier_rm"

# Repo/data root: $SPATIAL_REPO_ROOT if set, else the working dir (driver cds here).
repo_root <- Sys.getenv("SPATIAL_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- getwd()
repo_root <- normalizePath(repo_root)

# Outlier sample (Old 2) from the samplesheet; no real IDs are stored in this repo.
ss_path <- Sys.getenv("SPATIAL_SAMPLESHEET", unset = file.path(repo_root, "samplesheet.csv"))
if (!file.exists(ss_path)) {
  stop(sprintf("samplesheet not found at %s; set SPATIAL_SAMPLESHEET (see spatial/shared/samplesheet.template.csv)", ss_path))
}
ss <- read.csv(ss_path, stringsAsFactors = FALSE)
OUTLIER_SAMPLE <- ss$sample_id[ss$species == "mouse" & ss$status == "outlier"][1]

pb_dir   <- file.path(repo_root, "data", "pseudobulk", strategy)
out_dir  <- file.path(repo_root, "data", "de_results", strategy)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_types <- c("OPC", "Intermediate_Oligo", "Mature_Oligo")

# Same gene_classes mapping as the standard DE script (50 panel genes).
gene_classes <- tibble::tribble(
  ~gene,       ~gene_class,
  "Olig2",     "marker",       "Sox10",     "marker",
  "Pdgfra",    "marker",       "Ptprz1",    "marker",
  "Pcdh15",    "marker",       "Enpp6",     "marker",
  "Mbp",       "marker",       "Opalin",    "marker",
  "Mog",       "marker",
  "Bach2",     "primary_TF",   "Elf2",      "primary_TF",
  "Foxk2",     "primary_TF",   "Bhlhe41",   "primary_TF",
  "Nr6a1",     "primary_TF",   "Sox8",      "primary_TF",
  "Stat3",     "primary_TF",   "Sox5",      "primary_TF",
  "Klk6",      "primary_TF",
  "Arhgef12",  "Bach2_target", "Nr3c1",     "Bach2_target",
  "Fermt2",    "Bach2_target", "Cenpb",     "Bach2_target",
  "Slc39a3",   "Elf2_target",  "Gpatch4",   "Elf2_target",
  "Samd8",     "Elf2_target",
  "Aacs",      "Foxk2_target", "Abl2",      "Foxk2_target",
  "Uggt1",     "Foxk2_target",
  "Slc38a6",   "Bhlhe41_target", "Tada1",   "Bhlhe41_target",
  "Hspbap1",   "Bhlhe41_target",
  "Anks3",     "Nr6a1_target", "Bbs2",      "Nr6a1_target",
  "Nfib",      "Nr6a1_target", "Naa40",     "Nr6a1_target",
  "Selenoh",   "Sox8_target",  "Hsp90b1",   "Sox8_target",
  "Eif1b",     "Sox8_target",
  "Nup214",    "Stat3_target", "Trrap",     "Stat3_target",
  "Wdsub1",    "Stat3_target",
  "Pknox2",    "Sox5_inducer", "Rora",      "Sox5_inducer",
  "Zeb1",      "Sox5_inducer",
  "Plag1",     "Klk6_inducer", "Nkx2-9",    "Klk6_inducer",
  "Mitf",      "Klk6_inducer",
  "Cpa1",      "pancreas_control", "Spink1", "pancreas_control",
  "Nupr1",     "pancreas_or_stress"
)

run_de_for_celltype <- function(cell_type, strategy, pb_dir, out_dir) {
  counts_path  <- file.path(pb_dir, paste0(cell_type, "_counts.csv"))
  coldata_path <- file.path(pb_dir, paste0(cell_type, "_coldata.csv"))
  if (!file.exists(counts_path) || !file.exists(coldata_path)) {
    cat(sprintf("  [skip] %s: counts or coldata missing\n", cell_type))
    return(NULL)
  }

  counts  <- read_csv(counts_path, show_col_types = FALSE) %>% column_to_rownames("gene")
  coldata <- read_csv(coldata_path, show_col_types = FALSE) %>%
    mutate(age_group = factor(age_group, levels = c("Young", "Old"))) %>%
    column_to_rownames("sample_id")

  # Drop the outlier sample (if present)
  excluded_n_cells <- NA
  if (OUTLIER_SAMPLE %in% rownames(coldata)) {
    excluded_n_cells <- coldata[OUTLIER_SAMPLE, "n_cells"]
    coldata <- coldata[rownames(coldata) != OUTLIER_SAMPLE, , drop = FALSE]
    counts  <- counts[, colnames(counts) != OUTLIER_SAMPLE, drop = FALSE]
    cat(sprintf("  [exclude] %s: dropped %s (n_cells=%s)\n",
                cell_type, OUTLIER_SAMPLE, excluded_n_cells))
  } else {
    cat(sprintf("  [warn] %s: %s not found in coldata; nothing dropped\n",
                cell_type, OUTLIER_SAMPLE))
  }

  counts <- counts[, rownames(coldata), drop = FALSE]
  stopifnot(all(colnames(counts) == rownames(coldata)))

  cat(sprintf("  %s: %d samples (%d Young, %d Old), %d genes\n",
              cell_type, ncol(counts),
              sum(coldata$age_group == "Young"),
              sum(coldata$age_group == "Old"),
              nrow(counts)))

  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(counts),
    colData   = coldata,
    design    = ~ age_group
  )

  dds <- tryCatch(
    DESeq(dds, fitType = "parametric"),
    error = function(e) {
      cat(sprintf("    parametric fit failed (%s); retrying fitType='local'\n",
                  conditionMessage(e)))
      DESeq(dds, fitType = "local")
    }
  )

  res_raw <- results(dds, contrast = c("age_group", "Old", "Young"))
  res_shrunk <- tryCatch(
    lfcShrink(dds, contrast = c("age_group", "Old", "Young"),
              type = "normal", quiet = TRUE),
    error = function(e) {
      cat(sprintf("    lfcShrink failed: %s\n", conditionMessage(e)))
      res_raw
    }
  )

  de_df <- as.data.frame(res_raw) %>%
    rownames_to_column("gene") %>%
    mutate(
      log2FC_shrunk = res_shrunk$log2FoldChange[match(gene, rownames(res_shrunk))],
      cell_type = cell_type,
      strategy = strategy,
      excluded_sample = OUTLIER_SAMPLE
    ) %>%
    left_join(gene_classes, by = "gene") %>%
    mutate(gene_class = ifelse(is.na(gene_class), "other", gene_class)) %>%
    select(gene, gene_class, cell_type, strategy, excluded_sample,
           baseMean, log2FoldChange, log2FC_shrunk, lfcSE, stat, pvalue, padj)

  out_path <- file.path(out_dir, paste0(cell_type, "_de_", OUTLIER_TAG, ".csv"))
  write_csv(de_df, out_path)
  cat(sprintf("    wrote %s\n", out_path))

  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  vst_mat <- assay(vsd)
  vst_df <- as.data.frame(vst_mat) %>% rownames_to_column("gene")
  vst_path <- file.path(out_dir, paste0(cell_type, "_vst_", OUTLIER_TAG, ".csv"))
  write_csv(vst_df, vst_path)
  cat(sprintf("    wrote %s\n", vst_path))

  list(de = de_df, vst = vst_df, dds = dds, vsd = vsd)
}

cat(sprintf("Strategy: %s    OUTLIER excluded: %s\n", strategy, OUTLIER_SAMPLE))
cat(sprintf("Pseudobulk dir: %s\n", pb_dir))
cat(sprintf("Output dir: %s    suffix: _%s\n", out_dir, OUTLIER_TAG))
cat("\n")

results_list <- list()
for (ct in cell_types) {
  results_list[[ct]] <- run_de_for_celltype(ct, strategy, pb_dir, out_dir)
}

# Combined long-form DE table
combined <- bind_rows(lapply(results_list, function(x) if (is.null(x)) NULL else x$de))
combined_path <- file.path(out_dir, paste0("_combined_de_", OUTLIER_TAG, ".csv"))
write_csv(combined, combined_path)
cat(sprintf("\nCombined DE table: %s (%d rows)\n", combined_path, nrow(combined)))

cat(sprintf("\n=== Top hypothesis-gene hits per cell type (padj < 0.1; %s excluded) ===\n",
            OUTLIER_SAMPLE))
hypothesis_classes <- c("primary_TF", "Bach2_target", "Elf2_target", "Foxk2_target",
                        "Bhlhe41_target", "Nr6a1_target", "Sox8_target", "Stat3_target",
                        "Sox5_inducer", "Klk6_inducer")
for (ct in cell_types) {
  cat(sprintf("\n--- %s ---\n", ct))
  hits <- combined %>%
    filter(cell_type == ct, gene_class %in% hypothesis_classes,
           !is.na(padj), padj < 0.1) %>%
    arrange(padj)
  if (nrow(hits) == 0) {
    cat("  (no hypothesis genes with padj < 0.1)\n")
  } else {
    print(as.data.frame(hits %>% select(gene, gene_class, log2FoldChange, log2FC_shrunk, padj)))
  }
}

cat(sprintf("\n=== Marker findings (padj < 0.1; %s excluded) ===\n", OUTLIER_SAMPLE))
for (ct in cell_types) {
  cat(sprintf("\n--- %s ---\n", ct))
  hits <- combined %>%
    filter(cell_type == ct, gene_class == "marker",
           !is.na(padj), padj < 0.1) %>%
    arrange(padj)
  if (nrow(hits) == 0) {
    cat("  (no markers with padj < 0.1)\n")
  } else {
    print(as.data.frame(hits %>% select(gene, log2FoldChange, log2FC_shrunk, pvalue, padj)))
  }
}

cat("\nDone.\n")
