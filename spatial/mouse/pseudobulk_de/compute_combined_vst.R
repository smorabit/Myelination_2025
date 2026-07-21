#!/usr/bin/env Rscript
# Combined-celltype VST: builds a single DESeq2 dataset with all 30 (sample x
# cell_type) pseudobulk profiles, fits dispersions across the combined data,
# and writes the variance-stabilising-transformed expression matrix.
#
# This avoids the per-cell-type-fit failure mode that affects strategies with
# all-zero genes in one cell type (e.g. P5 OPC has 3 always-zero genes:
# Enpp6, Mog, Opalin -- by the precedence's "Mog/Opalin -> Mature" routing).
# Combined fit lands all profiles on a single normalisation scale and avoids
# DESeq2's parametric dispersion fit clamping to a constant for OPC.
#
# Usage: Rscript scripts/compute_combined_vst.R <strategy>
#   strategy: stringent | stringent_p4 | stringent_p5 | brain-ref | strict | legacy
#
# Inputs : data/pseudobulk/{strategy}/{cell_type}_counts.csv
#          data/pseudobulk/{strategy}/{cell_type}_coldata.csv
# Output : data/de_results/{strategy}/_combined_vst.csv
#          (rows = genes, columns = "{sample_id}__{cell_type}", values = VST)

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
strategy <- if (length(args) >= 1) args[1] else "stringent"
# Optional 2nd arg: a sample_id to exclude from the combined VST. When given,
# the output file is suffixed with "_outlier_rm" so it doesn't overwrite the
# full-cohort VST. Used to mirror de_young_vs_old_outlier_rm.R.
exclude_sample <- if (length(args) >= 2) args[2] else NA_character_
# Optional 3rd arg: a per-(sample x cell_type) floor on n_cells. When given,
# combos with n_cells < min_cells are dropped before the combined fit. Output
# suffix becomes "_floor<min_cells>". Mirrors de_young_vs_old_floor10.R.
min_cells_arg <- if (length(args) >= 3) suppressWarnings(as.integer(args[3])) else NA_integer_

# Repo/data root: $SPATIAL_REPO_ROOT if set, else the working dir (driver cds here).
repo_root <- Sys.getenv("SPATIAL_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- getwd()
repo_root <- normalizePath(repo_root)

pb_dir <- file.path(repo_root, "data", "pseudobulk", strategy)
out_dir <- file.path(repo_root, "data", "de_results", strategy)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cell_types <- c("OPC", "Intermediate_Oligo", "Mature_Oligo")
count_list <- list()
coldata_list <- list()

for (ct in cell_types) {
  counts_path <- file.path(pb_dir, paste0(ct, "_counts.csv"))
  coldata_path <- file.path(pb_dir, paste0(ct, "_coldata.csv"))
  if (!file.exists(counts_path)) stop(paste("Missing", counts_path))
  if (!file.exists(coldata_path)) stop(paste("Missing", coldata_path))
  counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  coldata <- read.csv(coldata_path)
  # Re-key columns as "{sample_id}__{cell_type}" to disambiguate
  new_cols <- paste0(colnames(counts), "__", ct)
  colnames(counts) <- new_cols
  coldata$profile_id <- paste0(coldata$sample_id, "__", ct)
  coldata$cell_type <- ct
  rownames(coldata) <- coldata$profile_id
  count_list[[ct]] <- counts
  coldata_list[[ct]] <- coldata
}

# Build combined matrix (assume same gene order across all 3 — they use the
# same panel)
genes <- rownames(count_list[[1]])
for (ct in cell_types[-1]) {
  if (!identical(rownames(count_list[[ct]]), genes)) {
    stop(paste("Gene order mismatch in", ct))
  }
}
combined_counts <- do.call(cbind, count_list)
combined_coldata <- do.call(rbind, coldata_list)
combined_coldata <- combined_coldata[colnames(combined_counts), ]

# Optional outlier removal: drop all profiles whose sample_id == exclude_sample
out_suffix <- ""
if (!is.na(exclude_sample) && nzchar(exclude_sample)) {
  keep_mask <- combined_coldata$sample_id != exclude_sample
  n_dropped <- sum(!keep_mask)
  combined_counts <- combined_counts[, keep_mask, drop = FALSE]
  combined_coldata <- combined_coldata[keep_mask, , drop = FALSE]
  out_suffix <- "_outlier_rm"
  cat(sprintf("Dropped %d profile(s) for sample_id=%s\n", n_dropped, exclude_sample))
}

# Optional per-combo floor: drop combos with n_cells < min_cells_arg
if (!is.na(min_cells_arg) && min_cells_arg > 0) {
  keep_mask <- combined_coldata$n_cells >= min_cells_arg
  dropped_rows <- combined_coldata[!keep_mask, c("sample_id", "cell_type", "n_cells")]
  if (nrow(dropped_rows) > 0) {
    cat(sprintf("Per-combo floor n_cells>=%d: dropping %d profile(s):\n",
                min_cells_arg, nrow(dropped_rows)))
    print(dropped_rows)
  } else {
    cat(sprintf("Per-combo floor n_cells>=%d: nothing to drop\n", min_cells_arg))
  }
  combined_counts <- combined_counts[, keep_mask, drop = FALSE]
  combined_coldata <- combined_coldata[keep_mask, , drop = FALSE]
  out_suffix <- paste0("_floor", min_cells_arg)
}

cat(sprintf("strategy: %s\n", strategy))
cat(sprintf("combined matrix: %d genes x %d profiles\n",
            nrow(combined_counts), ncol(combined_counts)))
cat(sprintf("cell types: %s\n", paste(unique(combined_coldata$cell_type), collapse = ", ")))

# Build DESeqDataSet. Design includes age_group + cell_type to model both;
# VST itself is design-aware via dispersion estimation.
combined_coldata$age_group <- factor(combined_coldata$age_group, levels = c("Young", "Old"))
combined_coldata$cell_type <- factor(combined_coldata$cell_type, levels = cell_types)

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(combined_counts),
  colData = combined_coldata,
  design = ~ age_group + cell_type
)

# Estimate size factors and dispersions on the combined data.
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = "parametric")

# VST. Use blind = FALSE so the design-aware dispersion estimates feed VST.
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_mat <- assay(vsd)

# Write CSV: rows = genes, columns = profile_id ("{sample_id}__{cell_type}")
out_path <- file.path(out_dir, paste0("_combined_vst", out_suffix, ".csv"))
out_df <- data.frame(gene = rownames(vst_mat), vst_mat, check.names = FALSE)
write.csv(out_df, out_path, row.names = FALSE)

cat(sprintf("Wrote %s (%d rows x %d cols)\n", out_path, nrow(out_df), ncol(out_df)))
