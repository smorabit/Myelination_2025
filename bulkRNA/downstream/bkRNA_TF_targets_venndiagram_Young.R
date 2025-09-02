###==================================
# Date: 01/29/2025 - 02/14/2025
# Robin JM Franklin
# Project: CNS_remyelination
# Ze
###==================================

###============###============###============###============
# NOTE Date 02/12/2025
###============###============###============###============

# NOTE loading library
set.seed(5)

library(ComplexHeatmap)
library(viridis)
library(dplyr)
library(MetBrewer)
library(NatParksPalettes)

library(cowplot)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(ggplot2)

theme_set(theme_cowplot())

library(dplyr)


###============
# PATH
data_in <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/data/'
outPATH <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/Analysis/'
dataPATH <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/data/processed_data/'
tfPATH <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/data/fromSam/stuff_ze/'

###============
# Load data

TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs

###============
# Load necessary libraries
library(dplyr)
library(eulerr)
library(GeneOverlap)

###============
# Loop over each TF
###============

for (cur_tf in TFnames) {

  cat("\n===========================================\n")
  cat("Processing TF:", cur_tf, "\n")
  cat("===========================================\n")

  ###============
  # Load data from Bulk RNA-seq
  cat("Step 1: Loading bulk RNA-seq data for", cur_tf, "...\n")
  bulkRNA <- read.delim(paste0(data_in, "20240830OverexpressionbulkRNAseq/", "GFP_vs_", cur_tf, ".txt"))

  ###============
  # Subset and process data
  cat("Step 2: Subsetting and processing data for", cur_tf, "...\n")
  bulkRNA_sub <- bulkRNA

  # bulkRNA_sub$log2FC_GFPvsTF <- bulkRNA_sub$Log2.Fold.Change..GFP_vs_.cur_tf.
  bulkRNA_sub$log2FC_GFPvsTF <- bulkRNA_sub[[paste0("Log2.Fold.Change..GFP_vs_", cur_tf, ".")]]
  # bulkRNA_sub$FDR <- bulkRNA_sub$FDR..GFP_vs_.cur_tf.
  bulkRNA_sub$FDR <- bulkRNA_sub[[paste0("FDR..GFP_vs_", cur_tf, ".")]]
  # bulkRNA_sub$Pvalue <- bulkRNA_sub$P.value..GFP_vs_.cur_tf.
  bulkRNA_sub$Pvalue <- bulkRNA_sub[[paste0("P.value..GFP_vs_", cur_tf, ".")]]

  bulkRNA_sub <- bulkRNA_sub %>% select(Feature, log2FC_GFPvsTF, Pvalue, FDR) #  Pvalue,
  # bulkRNA_sub <- subset(bulkRNA_sub, FDR < 0.1)
  bulkRNA_sub <- subset(bulkRNA_sub, Pvalue < 0.1)

  bulkRNA_sub$gene <- bulkRNA_sub$Feature

  # deg_bulkRNA <- bulkRNA_sub$gene[bulkRNA_sub$log2FC_GFPvsTF > 0.5 | bulkRNA_sub$log2FC_GFPvsTF < -0.5]
  # deg_bulkRNA <- bulkRNA_sub$gene[bulkRNA_sub$log2FC_GFPvsTF > 0 | bulkRNA_sub$log2FC_GFPvsTF < -0]
  deg_bulkRNA <- bulkRNA_sub$gene

  ###============
  # Load data from Sam
  cat("Step 3: Loading ODC targets data for", cur_tf, "...\n")
  ODC_Young_targets <- read.csv(paste0(tfPATH, "ODC_Young_targets_", cur_tf, ".csv"))

  missing_targets <- ODC_Young_targets$target[!(ODC_Young_targets$target %in% bulkRNA$Feature)]
  ODC_Young_targets <- ODC_Young_targets %>% filter(!(target %in% missing_targets))

  ###============
  # Extract primary and secondary genes
  cat("Step 4: Extracting primary and secondary genes for", cur_tf, "...\n")
  primary_genes <- ODC_Young_targets %>% subset(depth == "1") %>% .$target
  secondary_genes <- subset(ODC_Young_targets, depth == "2") %>% .$target

  primary_genes <- setdiff(primary_genes, cur_tf)
  secondary_genes <- setdiff(secondary_genes, cur_tf)

  ###============
  # Check for overlapping genes
  cat("Step 5: Checking for duplicated genes for", cur_tf, "...\n")
  duplicated_genes <- intersect(primary_genes, secondary_genes)

  if (length(duplicated_genes) > 0) {
    cat("Duplicated genes found for", cur_tf, ":\n")
    print(duplicated_genes)
  } else {
    cat("No duplicated genes found for", cur_tf, "\n")
  }

  secondary_genes <- setdiff(secondary_genes, primary_genes)

  ###============
  # Find intersecting genes
  cat("Step 6: Finding intersecting genes for", cur_tf, "...\n")
  intersect_gene <- intersect(c(primary_genes, secondary_genes, cur_tf), deg_bulkRNA)

  ###============
  # Create Euler diagram
  cat("Step 7: Creating Euler diagram for", cur_tf, "...\n")
  p <- euler(c(
    primary_secondary_target = length(c(primary_genes, secondary_genes, cur_tf)),
    bulkRNA = length(deg_bulkRNA),
    "primary_secondary_target&bulkRNA" = length(intersect_gene)
  ))

  plot_p <- plot(
    p,
    fills = c("seagreen", "#72C2D6", '#50A696'),  #
    labels = list(font = 1, cex = 1),  # Label font and size
    main = paste("Overlap for", cur_tf)  # Title
  )

  ###============
  # Save Euler diagram
  cat("Step 8: Saving Euler diagram for", cur_tf, "...\n")
  ggsave(paste0(outPATH, "bulkRNA/ODC_Young/", "GFPvs", cur_tf, "_GeneOverlap_w", cur_tf,"_targets_venn.pdf"), plot = plot_p, width = 8, height = 6)

  ###============
  # Gene Overlap Analysis
  cat("Step 9: Performing gene overlap analysis for", cur_tf, "...\n")
  go.obj <- newGeneOverlap(c(primary_genes, secondary_genes, cur_tf),
                           deg_bulkRNA,
                           genome.size = length(unique(bulkRNA$Feature))
  )  # Fixed missing parenthesis


  go.obj <- testGeneOverlap(go.obj)

  ###============
  # Save gene overlap results
  cat("Step 10: Saving gene overlap results for", cur_tf, "...\n")
  sink(paste0(outPATH, "bulkRNA/ODC_Young/", "GFPvs", cur_tf, "_GeneOverlap_w", cur_tf,"_targets_calculation.txt"))
  print(go.obj)
  sink()

  cat("Processing complete for TF:", cur_tf, "\n")
  cat("===========================================\n\n")
}
