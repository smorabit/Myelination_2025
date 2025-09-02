###==================================
# Date: 01/29/2025 - 02/14/2025
# Robin JM Franklin
# Project: CNS_remyelination
# Ze
# conda activate SeuratV4
###==================================

###============###============###============###============
# NOTE Date 02/12/2025
###============###============###============###============

# set.seed(123456)
# # library(Seurat)
# # library(Signac)
# library(tidyverse)
# library(cowplot)
# library(ggrepel)
# library(GenomicRanges)
# library(RColorBrewer)
# # loading library
# library(viridis)
# library(dplyr)
# # install.packages("igraph")
# library(igraph)
# library(qgraph)
#
# library(ComplexHeatmap)
#

# PATH
# data_EX_in
data_in <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/data/20240830OverexpressionbulkRNAseq/'
tfPATH <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/data/fromSam/stuff_ze/'
data_out <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/Analysis/bulkRNA/deg/'

TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs


#================ NOTE

#================ NOTE Plotting

library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggrastr)

TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs

# Define colors for the volcano plot
color1 = '#FEBEB1'
color2 = '#B3E2E2'
color3 = "seagreen"
color4 = "darkorchid3"
color5 = "#5F6897"

# Parameters for the volcano plot
nlabel_volcano = 10
volcano_width = 5
volcano_height = 5
make_volcano = TRUE
raster = TRUE
dpi_num = 800

# Loop through each TF in the TFnames list
for (cur_tf in TFnames) {
  cat("Processing TF:", cur_tf, "\n")

  # Load bulk RNA data
  bulkRNA <- read.delim(paste0(data_in, 'GFP_vs_', cur_tf, '.txt'))
  bulkRNA$gene <- bulkRNA$Feature
  bulkRNA$Log2FC <- bulkRNA[[paste0("Log2.Fold.Change..GFP_vs_", cur_tf, ".")]]
  bulkRNA$p_val_adj <- bulkRNA[[paste0("FDR..GFP_vs_", cur_tf, ".")]]
  bulkRNA$p_val <- bulkRNA[[paste0("P.value..GFP_vs_", cur_tf, ".")]]

  # when have the duplicated gene names, only keep the row with the largest absolute log2FC
  bulkRNA <- bulkRNA %>%
    group_by(gene) %>%
    filter(abs(Log2FC) == max(abs(Log2FC))) %>%
    ungroup() %>% as.data.frame()
  #

  # Filter significant DEGs
  # signif_DEG <- subset(bulkRNA, p_val_adj < 0.1)
  signif_DEG <- subset(bulkRNA, p_val < 0.1)
  deg_bulkRNA <- signif_DEG$gene

  # Load ODC targets data for Old and Young
  ODC_Old_targets <- read.csv(paste0(tfPATH, "ODC_Old_targets_", cur_tf, ".csv"))
  ODC_Young_targets <- read.csv(paste0(tfPATH, "ODC_Young_targets_", cur_tf, ".csv"))

  # Filter out missing targets
  ODC_Old_targets <- ODC_Old_targets %>% filter(target %in% bulkRNA$Feature)
  ODC_Young_targets <- ODC_Young_targets %>% filter(target %in% bulkRNA$Feature)

  # Extract primary and secondary genes for Old and Young
  primary_genes_Old <- ODC_Old_targets %>% subset(depth == "1") %>% .$target
  secondary_genes_Old <- subset(ODC_Old_targets, depth == "2") %>% .$target
  primary_genes_Young <- ODC_Young_targets %>% subset(depth == "1") %>% .$target
  secondary_genes_Young <- subset(ODC_Young_targets, depth == "2") %>% .$target

  # Remove the current TF from the primary and secondary genes
  primary_genes_Old <- setdiff(primary_genes_Old, cur_tf)
  secondary_genes_Old <- setdiff(secondary_genes_Old, cur_tf)
  primary_genes_Young <- setdiff(primary_genes_Young, cur_tf)
  secondary_genes_Young <- setdiff(secondary_genes_Young, cur_tf)

  # Check for duplicated genes
  duplicated_genes_Old <- intersect(primary_genes_Old, secondary_genes_Old)
  duplicated_genes_Young <- intersect(primary_genes_Young, secondary_genes_Young)

  if (length(duplicated_genes_Old) > 0) {
    cat("Duplicated genes found for", cur_tf, "in Old targets:\n")
    print(length(duplicated_genes_Old))
  }

  if (length(duplicated_genes_Young) > 0) {
    cat("Duplicated genes found for", cur_tf, "in Young targets:\n")
    print(length(duplicated_genes_Young) )
  }

  secondary_genes_Old <- setdiff(secondary_genes_Old, primary_genes_Old)
  secondary_genes_Young <- setdiff(secondary_genes_Young, primary_genes_Young)

  # Find intersecting genes
  intersect_gene_Old <- intersect(c(primary_genes_Old, secondary_genes_Old, cur_tf), deg_bulkRNA)
  intersect_gene_Young <- intersect(c(primary_genes_Young, secondary_genes_Young, cur_tf), deg_bulkRNA)

  shared_gene <- intersect(intersect_gene_Young, intersect_gene_Old)
  intersect_gene_Young <- setdiff(intersect_gene_Young, shared_gene)
  intersect_gene_Old <- setdiff(intersect_gene_Old, shared_gene)

  # Prepare data for volcano plot
  # Prepare anno for volcano plot
  filtered_cur_degs <- bulkRNA[bulkRNA$gene %in% c(shared_gene, intersect_gene_Young, intersect_gene_Old), ]
  top_thresh_genes <- filtered_cur_degs %>% top_n(nlabel_volcano, wt = Log2FC) %>% .$gene
  bottom_thresh_genes <- filtered_cur_degs %>% top_n(-1 * nlabel_volcano, wt = Log2FC) %>% .$gene

  bulkRNA$anno <- ifelse(bulkRNA$gene %in% c(top_thresh_genes, bottom_thresh_genes), bulkRNA$gene, NA)
  bulkRNA$anno <- ifelse(bulkRNA$gene %in% cur_tf, cur_tf, bulkRNA$anno)
  bulkRNA$anno <- ifelse(bulkRNA$p_val >= 0.1, NA, bulkRNA$anno)

  # Prepare color for volcano plot
  bulkRNA$color <- ifelse(bulkRNA$p_val >= 0.1, 'gray', ifelse(bulkRNA$Log2FC > 0, color1, color2)) # p_val_adj
  bulkRNA$color <- ifelse(bulkRNA$gene %in% shared_gene, color5, bulkRNA$color)
  bulkRNA$color <- ifelse(bulkRNA$gene %in% intersect_gene_Young, color4, bulkRNA$color)
  bulkRNA$color <- ifelse(bulkRNA$gene %in% intersect_gene_Old, color3, bulkRNA$color)
  bulkRNA$color <- ifelse(bulkRNA$gene %in% cur_tf, 'red', bulkRNA$color)

  # Generate volcano plot
  plot_degs <- bulkRNA
  p <- plot_degs %>%
    ggplot(aes(x = Log2FC, y = -log10(p_val))) + # p_val_adj
    geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
    ggrastr::rasterise(geom_point(alpha = 0.5, color = plot_degs$color), dpi = dpi_num) +
    geom_point(inherit.aes = FALSE, data = subset(plot_degs, !is.na(anno)), aes(Log2FC, -log10(p_val)), # p_val_adj
               fill = subset(plot_degs, !is.na(anno)) %>% .$color, shape = 21, size = 2, color = 'black') +
    geom_text_repel(aes(label = anno), color = 'black', fontface = 'italic', min.segment.length = 0, max.overlaps = 100) +
    xlim(-1 * max(abs(plot_degs$Log2FC)) - 0.1, max(abs(plot_degs$Log2FC)) + 0.1) +
    ggtitle(paste0("Bulk RNAseq DEG analysis: GFP vs ", cur_tf)) +
    theme(
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = 'bottom'
    )

  # Save the volcano plot
  pdf(paste0(data_out, '/volcano_Bulk_RNA_', cur_tf, '.pdf'), width = volcano_width, height = volcano_height, useDingbats = FALSE)
  print(p)
  dev.off()

  genenames <- unique(c(cur_tf, top_thresh_genes, bottom_thresh_genes))
  sink(paste0(data_out, "/volcano_Bulk_RNA_", cur_tf, "_plotted_genenames.txt"))
  print(genenames)
  sink()

  cat("Processing complete for TF:", cur_tf, "\n")

}
