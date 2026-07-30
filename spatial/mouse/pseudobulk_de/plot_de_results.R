#!/usr/bin/env Rscript
# DE results visualization. Args: <strategy> <suffix_tag> <barplot_method>
# (canonical: stringent_p5 outlier_rm counts, then a second vst run for volcanoes/heatmap).
#   - Per-gene barplot: 3 cell types on x-axis, Young vs Old fill, individual sample dots
#   - Volcano plot per cell type
#   - Heatmap of log2FC across hypothesis genes x cell types
#
# Inputs (suffixed by the chosen strategy/tag):
#   data/de_results/{strategy}/_combined_de_{tag}.csv
#   data/de_results/{strategy}/_combined_vst_{tag}.csv
#   data/pseudobulk/{strategy}/{cell_type}_coldata.csv
#
# Outputs (to docs/images/):
#   de_barplot_{gene}_{strategy}_{tag}[_counts].png   one per hypothesis gene
#   de_volcano_{cell_type}_{strategy}_{tag}.png       one per cell type
#   de_heatmap_log2fc_{strategy}_{tag}.png

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
strategy <- if (length(args) >= 1) args[1] else "stringent"
# Optional 2nd arg = suffix tag (e.g. "outlier_rm"). When given, the script
# reads `_combined_de_<tag>.csv`, `_combined_vst_<tag>.csv` (or per-celltype
# `<ct>_vst_<tag>.csv`) and writes images suffixed `_<strategy>_<tag>.png`.
suffix_tag <- if (length(args) >= 2) args[2] else ""
suffix_in <- if (nzchar(suffix_tag)) paste0("_", suffix_tag) else ""
suffix_out <- if (nzchar(suffix_tag)) paste0("_", suffix_tag) else ""

# Optional 3rd arg = barplot value source:
#   "vst"     (default; uses combined VST from compute_combined_vst.R)
#   "logcpm"  (per-profile log2(CPM+1) on raw pseudobulk counts)
#   "cpm"     (per-profile raw CPM, no log transform)
#   "counts"  (per-profile mean counts per cell — most biologically intuitive;
#              matches the unit used on the TF_expression.md heatmap)
# When non-vst, per-gene barplot images get a matching suffix on the filename,
# and y-axis label / subtitle reflect the value source. Volcanoes and the
# log2FC heatmap are unaffected.
barplot_method <- if (length(args) >= 3) args[3] else "vst"
if (!barplot_method %in% c("vst", "logcpm", "cpm", "counts")) {
  stop(sprintf("3rd arg must be 'vst', 'logcpm', 'cpm', or 'counts'; got '%s'",
               barplot_method))
}
barplot_suffix <- switch(barplot_method,
                         "logcpm" = "_logcpm",
                         "cpm"    = "_cpm",
                         "counts" = "_counts",
                         "")

# Repo/data root: $SPATIAL_REPO_ROOT if set, else the working dir (driver cds here).
repo_root <- Sys.getenv("SPATIAL_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- getwd()
repo_root <- normalizePath(repo_root)

de_dir   <- file.path(repo_root, "data", "de_results", strategy)
pb_dir   <- file.path(repo_root, "data", "pseudobulk", strategy)
img_dir  <- file.path(repo_root, "docs", "images")
dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)

cell_types <- c("OPC", "Intermediate_Oligo", "Mature_Oligo")
cell_type_labels <- c(
  OPC = "OPC",
  Intermediate_Oligo = "Intermediate",
  Mature_Oligo = "Mature"
)

age_colors <- c(Young = "#56B4E9", Old = "#D55E00")

# ---------- Load data ----------
de <- read_csv(file.path(de_dir, paste0("_combined_de", suffix_in, ".csv")), show_col_types = FALSE)

# ---------- TF-group classification + shared colour palette ----------
# Each primary TF sits at the head of its own group (TF + its targets/inducers).
# Group colours are shared between the heatmap row-annotation strip and the
# volcano dot colours so visual identification is consistent across plots.
# Markers (cell-type-defining, not expected to shift Young vs Old) get black;
# non-hypothesis controls get neutral grey.
tf_group_order <- c("Bach2", "Elf2", "Foxk2", "Bhlhe41", "Nr6a1",
                    "Sox8", "Stat3", "Sox5", "Klk6")
# Colourblind-safe TF-group palette. The 9 TFs use Paul Tol's "muted" qualitative
# scheme (distinguishable under deuteranopia/protanopia/tritanopia); mirrors
# plot_style.TF_GROUP_COLORS on the Python side so every figure shares one map.
tf_group_colors <- c(
  Bach2             = "#CC6677",   # rose
  Elf2              = "#332288",   # indigo
  Foxk2             = "#DDCC77",   # sand
  Bhlhe41           = "#AA4499",   # purple
  Nr6a1             = "#88CCEE",   # cyan
  Sox8              = "#882255",   # wine (headline)
  Stat3             = "#117733",   # green
  Sox5              = "#999933",   # olive
  Klk6              = "#44AA99",   # teal (headline)
  marker            = "#000000",   # black, cell-type markers, expect no Y vs O shift
  `pancreas/stress` = "#888888",   # mid-grey
  other             = "#cccccc"    # light grey
)

de <- de %>%
  mutate(
    tf_group = case_when(
      gene_class == "primary_TF"      ~ gene,           # primary TFs map to their own name
      gene_class == "Bach2_target"    ~ "Bach2",
      gene_class == "Elf2_target"     ~ "Elf2",
      gene_class == "Foxk2_target"    ~ "Foxk2",
      gene_class == "Bhlhe41_target"  ~ "Bhlhe41",
      gene_class == "Nr6a1_target"    ~ "Nr6a1",
      gene_class == "Sox8_target"     ~ "Sox8",
      gene_class == "Stat3_target"    ~ "Stat3",
      gene_class == "Sox5_inducer"    ~ "Sox5",
      gene_class == "Klk6_inducer"    ~ "Klk6",
      gene_class == "marker"          ~ "marker",
      gene_class %in% c("pancreas_control", "pancreas_or_stress") ~ "pancreas/stress",
      TRUE                            ~ "other"
    )
  )

# Load VST in long-form. Prefer the combined-celltype VST (single DESeq2 fit
# across all 30 sample × cell_type profiles, written by
# scripts/compute_combined_vst.R) when available — this avoids the per-cell-
# type-fit failure mode that clamps OPC VST values to a constant when the
# class has multiple all-zero genes (e.g. P5 OPC has Enpp6, Mog, Opalin all
# zero by construction). Fall back to per-cell-type VST if the combined
# file isn't present.
combined_vst_path <- file.path(de_dir, paste0("_combined_vst", suffix_in, ".csv"))
if (file.exists(combined_vst_path)) {
  cat(sprintf("Using combined VST: %s\n", combined_vst_path))
  cvst <- read_csv(combined_vst_path, show_col_types = FALSE)
  vst_long <- cvst %>%
    pivot_longer(-gene, names_to = "profile_col", values_to = "vst")
  # Column names are "{ct}.{sample_id}__{ct}" (R cbind on a named list
  # prepends the list name with a "." separator, then we appended "__{ct}"
  # in compute_combined_vst.R for explicitness). Recover (sample_id, ct).
  parse_profile <- function(s) {
    # Strip trailing "__<ct>" suffix
    parts <- strsplit(s, "__", fixed = TRUE)[[1]]
    ct <- tail(parts, 1)
    rest <- paste(head(parts, -1), collapse = "__")
    # Strip leading "<ct>." prefix if present (R cbind list-name artefact)
    pfx <- paste0(ct, ".")
    if (startsWith(rest, pfx)) {
      rest <- substring(rest, nchar(pfx) + 1)
    }
    list(sample_id = rest, cell_type = ct)
  }
  parsed <- lapply(unique(vst_long$profile_col), parse_profile)
  parsed_df <- tibble(
    profile_col = unique(vst_long$profile_col),
    sample_id = sapply(parsed, function(x) x$sample_id),
    cell_type = sapply(parsed, function(x) x$cell_type)
  )
  vst_long <- vst_long %>%
    inner_join(parsed_df, by = "profile_col") %>%
    select(-profile_col)
  # Add age_group from any per-celltype coldata file
  age_lookup <- bind_rows(lapply(cell_types, function(ct) {
    p <- file.path(pb_dir, paste0(ct, "_coldata.csv"))
    if (!file.exists(p)) return(NULL)
    read_csv(p, show_col_types = FALSE) %>% select(sample_id, age_group)
  })) %>% distinct(sample_id, age_group) %>%
    mutate(age_group = factor(age_group, levels = c("Young", "Old")))
  vst_long <- vst_long %>% inner_join(age_lookup, by = "sample_id")
} else {
  cat(sprintf("Combined VST not found (%s); falling back to per-celltype VST.\n",
              combined_vst_path))
  vst_long <- bind_rows(lapply(cell_types, function(ct) {
    vst_path     <- file.path(de_dir, paste0(ct, "_vst", suffix_in, ".csv"))
    coldata_path <- file.path(pb_dir, paste0(ct, "_coldata.csv"))
    if (!file.exists(vst_path) || !file.exists(coldata_path)) {
      return(NULL)
    }
    vst <- read_csv(vst_path, show_col_types = FALSE)
    coldata <- read_csv(coldata_path, show_col_types = FALSE) %>%
      mutate(age_group = factor(age_group, levels = c("Young", "Old")))
    vst %>%
      pivot_longer(-gene, names_to = "sample_id", values_to = "vst") %>%
      inner_join(coldata, by = "sample_id") %>%
      mutate(cell_type = ct)
  }))
}

# Optional: replace vst_long with log2(CPM+1) values for barplots only.
# (Heatmap and volcanoes use the `de` table directly — log2FC and padj from
# DESeq2 raw counts — and are unaffected.) Library size = column sum across
# all 50 panel genes per pseudobulk profile.
if (barplot_method %in% c("logcpm", "cpm", "counts")) {
  cat(sprintf(
    "Barplot value source: %s per profile (raw counts).\n",
    switch(barplot_method,
           "logcpm" = "log2(CPM+1)",
           "cpm"    = "CPM",
           "counts" = "mean counts per cell")
  ))
  cpm_long <- bind_rows(lapply(cell_types, function(ct) {
    counts_path  <- file.path(pb_dir, paste0(ct, "_counts.csv"))
    coldata_path <- file.path(pb_dir, paste0(ct, "_coldata.csv"))
    if (!file.exists(counts_path) || !file.exists(coldata_path)) return(NULL)
    counts  <- read_csv(counts_path, show_col_types = FALSE) %>% column_to_rownames("gene")
    coldata <- read_csv(coldata_path, show_col_types = FALSE) %>%
      mutate(age_group = factor(age_group, levels = c("Young", "Old")))
    # Optional outlier removal applies here too (consume suffix_tag)
    if (suffix_tag == "outlier_rm") {
      # Outlier sample from the samplesheet; no real IDs are stored in this repo.
      ss_path <- Sys.getenv("SPATIAL_SAMPLESHEET", unset = file.path(repo_root, "samplesheet.csv"))
      if (!file.exists(ss_path)) {
        stop(sprintf("samplesheet not found at %s; set SPATIAL_SAMPLESHEET (see spatial/shared/samplesheet.template.csv)", ss_path))
      }
      ss <- read.csv(ss_path, stringsAsFactors = FALSE)
      excl <- ss$sample_id[ss$species == "mouse" & ss$status == "outlier"][1]
      keep <- coldata$sample_id != excl
      coldata <- coldata[keep, ]
      counts  <- counts[, coldata$sample_id, drop = FALSE]
    }
    values <- switch(barplot_method,
      "logcpm" = {
        lib <- colSums(counts)
        log2(sweep(counts, 2, lib, "/") * 1e6 + 1)
      },
      "cpm" = {
        lib <- colSums(counts)
        sweep(counts, 2, lib, "/") * 1e6
      },
      "counts" = {
        # mean counts per cell: divide each profile's gene counts by its n_cells
        nc <- setNames(coldata$n_cells, coldata$sample_id)[colnames(counts)]
        sweep(counts, 2, nc, "/")
      }
    )
    val_df <- as.data.frame(values)
    val_df$gene <- rownames(val_df)
    val_df %>%
      pivot_longer(-gene, names_to = "sample_id", values_to = "vst") %>%
      inner_join(coldata, by = "sample_id") %>%
      mutate(cell_type = ct)
  }))
  vst_long <- cpm_long
} else {
  cat(sprintf("Barplot value source: VST (combined-celltype DESeq2 if available).\n"))
}

# Hypothesis-gene set ordering (drives the gene tab order in the docs page)
hypothesis_classes <- c("primary_TF", "Bach2_target", "Elf2_target", "Foxk2_target",
                        "Bhlhe41_target", "Nr6a1_target", "Sox8_target", "Stat3_target",
                        "Sox5_inducer", "Klk6_inducer")
all_hypothesis_genes <- de %>%
  filter(gene_class %in% hypothesis_classes) %>%
  distinct(gene, gene_class) %>%
  arrange(factor(gene_class, levels = hypothesis_classes), gene) %>%
  pull(gene)

# ---------- Per-gene barplots ----------
plot_gene_barplot <- function(g) {
  d <- vst_long %>%
    filter(gene == g) %>%
    mutate(cell_type = factor(cell_type, levels = cell_types,
                              labels = cell_type_labels[cell_types]))
  if (nrow(d) == 0) return(invisible(NULL))

  # Per-cell-type stats from DE table
  ann <- de %>%
    filter(gene == g, cell_type %in% cell_types) %>%
    mutate(
      cell_type = factor(cell_type, levels = cell_types,
                         labels = cell_type_labels[cell_types]),
      sig_label = case_when(
        is.na(pvalue) ~ "",
        pvalue < 0.001 ~ sprintf("p=%.1e", pvalue),
        TRUE ~ sprintf("p=%.3f", pvalue)
      ),
      lfc_label = sprintf("log2FC=%.2f", log2FoldChange)
    )

  summary_df <- d %>%
    group_by(cell_type, age_group) %>%
    summarise(mean = mean(vst), se = sd(vst) / sqrt(n()), n = n(), .groups = "drop")

  y_max <- max(d$vst, na.rm = TRUE)
  y_top <- y_max + (y_max - min(d$vst, na.rm = TRUE)) * 0.15

  # Get gene class for title — rewrite for publication (replace underscores
  # with spaces so "primary_TF" → "primary TF", "Bach2_target" → "Bach2 target").
  gc_raw <- (de %>% filter(gene == g) %>% pull(gene_class))[1]
  gc <- gsub("_", " ", gc_raw)

  p <- ggplot() +
    geom_col(
      data = summary_df,
      aes(x = cell_type, y = mean, fill = age_group),
      position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85,
      color = "black", linewidth = 0.3
    ) +
    geom_errorbar(
      data = summary_df,
      aes(x = cell_type, ymin = mean - se, ymax = mean + se, group = age_group),
      position = position_dodge(width = 0.8), width = 0.25
    ) +
    geom_jitter(
      data = d,
      aes(x = cell_type, y = vst, fill = age_group),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      shape = 21, size = 2, alpha = 0.9, stroke = 0.4
    ) +
    geom_text(
      data = ann,
      aes(x = cell_type, y = y_top, label = paste0(lfc_label, "\n", sig_label)),
      size = 3.4, lineheight = 0.95
    ) +
    scale_fill_manual(values = age_colors, name = "Age") +
    labs(
      title = bquote(italic(.(g)) ~ "(" * .(gc) * ")"),
      x = NULL,
      y = switch(barplot_method,
        "logcpm" = "log₂(CPM + 1)",
        "cpm"    = "CPM",
        "counts" = "mean counts / cell",
        "VST-normalised expression"
      )
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 13),
      axis.title.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9)
    )

  # Publication-quality output: PNG (raster, for docs site) + PDF (vector).
  out_base <- file.path(img_dir, sprintf("de_barplot_%s_%s%s%s",
                                         gsub("-", "_", g), strategy, suffix_out, barplot_suffix))
  ggsave(paste0(out_base, ".png"), p, width = 4.5, height = 4.8, dpi = 200)
  ggsave(paste0(out_base, ".pdf"), p, width = 4.5, height = 4.8, device = cairo_pdf)
  invisible(out_base)
}

cat(sprintf("Plotting per-gene barplots (%d genes)...\n", length(all_hypothesis_genes)))
for (g in all_hypothesis_genes) {
  plot_gene_barplot(g)
}
cat(sprintf("  wrote %d barplots to %s\n", length(all_hypothesis_genes),
            file.path("docs", "images")))

# ---------- Volcano plots, one per cell type ----------
plot_volcano <- function(ct) {
  d <- de %>% filter(cell_type == ct, !is.na(pvalue))
  if (nrow(d) == 0) return(invisible(NULL))

  d <- d %>%
    mutate(
      neglog10p = -log10(pvalue),
      tf_group = factor(tf_group,
                        levels = c(tf_group_order, "marker", "pancreas/stress", "other")),
      # Label: every primary TF (their own tf_group is their name → in tf_group_order),
      # every marker (always labelled — user wants to see they don't move), and any
      # other gene with p < 0.05.
      label = ifelse(
        pvalue < 0.05
          | gene_class == "primary_TF"
          | gene_class == "marker",
        gene, ""
      )
    )

  p <- ggplot(d, aes(x = log2FoldChange, y = neglog10p, color = tf_group)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey80") +
    geom_point(size = 2.5, alpha = 0.85) +
    ggrepel::geom_text_repel(
      aes(label = label), size = 3, max.overlaps = 30, box.padding = 0.3,
      segment.alpha = 0.5, fontface = "italic"
    ) +
    scale_color_manual(values = tf_group_colors, name = "TF group",
                       drop = FALSE) +
    labs(
      title = sprintf("Old vs Young, %s", cell_type_labels[ct]),
      x = "log2 fold change (Old / Young)",
      y = expression(-log[10](p))
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 13),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 9)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

  out_path <- file.path(img_dir, sprintf("de_volcano_%s_%s%s.png", ct, strategy, suffix_out))
  ggsave(out_path, p, width = 8, height = 6.5, dpi = 200)
  invisible(out_path)
}

if (barplot_method == "vst") {
  cat("\nPlotting volcanoes...\n")
  for (ct in cell_types) {
    out <- plot_volcano(ct)
    if (!is.null(out)) cat(sprintf("  %s -> %s\n", ct, basename(out)))
  }
} else {
  cat("\n[skip volcanoes — barplot_method=logcpm; volcanoes use log2FC only and are unchanged.]\n")
}

# ---------- log2FC heatmap (hypothesis genes x cell types) ----------
# Genes are grouped by TF group (primary TF + its targets/inducers in one
# block). The primary TF sits at the top of each block; targets/inducers
# are hierarchically clustered (Euclidean / complete linkage) so the most-
# similar log2FC profiles sit adjacent within the block. The row-annotation
# strip uses the same colour palette as the volcano dots, so TF identity is
# consistent across plots. Uses un-shrunk DESeq2 `log2FoldChange` (matches
# the prose / TL;DR / barplot annotations / volcano captions).
plot_heatmap <- function() {
  d <- de %>%
    filter(gene_class %in% hypothesis_classes, cell_type %in% cell_types) %>%
    select(gene, gene_class, tf_group, cell_type, log2FoldChange)

  mat_wide <- d %>%
    pivot_wider(names_from = cell_type, values_from = log2FoldChange)

  # Order within each tf_group: primary TF first, then hierarchically cluster
  # the targets/inducers. Groups emitted in tf_group_order so the heatmap
  # reads top-to-bottom as Bach2 / Elf2 / Foxk2 / ... / Klk6.
  ordered_rows <- list()
  for (grp in tf_group_order) {
    sub <- mat_wide %>% filter(tf_group == grp)
    if (nrow(sub) == 0) next
    primary_row <- sub %>% filter(gene_class == "primary_TF")
    target_rows <- sub %>% filter(gene_class != "primary_TF")
    if (nrow(target_rows) > 1) {
      tmat <- as.matrix(target_rows[, cell_types])
      rownames(tmat) <- target_rows$gene
      tmat[is.na(tmat)] <- 0
      hc <- hclust(dist(tmat), method = "complete")
      target_rows <- target_rows[hc$order, ]
    }
    ordered_rows[[grp]] <- bind_rows(primary_row, target_rows)
  }
  mat_ordered <- bind_rows(ordered_rows)

  rownames_use <- mat_ordered$gene
  m <- as.matrix(mat_ordered[, cell_types])
  rownames(m) <- rownames_use
  colnames(m) <- cell_type_labels[colnames(m)]

  ann_row <- data.frame(
    `TF group` = factor(mat_ordered$tf_group, levels = tf_group_order),
    row.names = rownames_use,
    check.names = FALSE
  )
  ann_colors <- list(`TF group` = tf_group_colors[tf_group_order])

  # Gap positions: cumulative group sizes (row indices where each TF block ends).
  group_sizes <- sapply(tf_group_order, function(g) sum(mat_ordered$tf_group == g))
  group_sizes <- group_sizes[group_sizes > 0]
  gaps_row <- head(cumsum(group_sizes), -1)

  out_base <- file.path(img_dir, sprintf("de_heatmap_log2fc_%s%s", strategy, suffix_out))
  # 7 x 11 in canvas: wide enough that the title fits without wrapping and the
  # OPC/Intermediate/Mature x-axis labels are not truncated. Title is a single
  # short line; the "grouped by TF" detail is implicit in the row-strip colours
  # + annotation legend (and described in the figure caption / methods).
  hm_args <- list(
    m,
    cluster_rows = FALSE,  # rows ordered manually within tf_group
    cluster_cols = FALSE,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    annotation_names_row = FALSE,
    gaps_row = gaps_row,
    color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101),
    breaks = seq(-2.5, 2.5, length.out = 102),
    main = "log2FC (Old / Young)",
    fontsize_row = 10, fontsize_col = 12,
    fontsize = 13,
    display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8,
    silent = FALSE,  # silent=TRUE blanks the active device (PNG/PDF would be empty)
    treeheight_row = 0, treeheight_col = 0
  )
  # PNG (raster, for docs site)
  png(paste0(out_base, ".png"), width = 1400, height = 2200, res = 200)
  do.call(pheatmap::pheatmap, hm_args)
  dev.off()
  # PDF (vector, for publication)
  cairo_pdf(paste0(out_base, ".pdf"), width = 7, height = 11)
  do.call(pheatmap::pheatmap, hm_args)
  dev.off()
  invisible(paste0(out_base, ".png"))
}

if (barplot_method == "vst") {
  cat("\nPlotting heatmap...\n")
  hm_path <- plot_heatmap()
  cat(sprintf("  -> %s\n", basename(hm_path)))
} else {
  cat("\n[skip heatmap — barplot_method=logcpm; heatmap uses log2FC only and is unchanged.]\n")
}

cat("\nAll plots written to docs/images/\n")
