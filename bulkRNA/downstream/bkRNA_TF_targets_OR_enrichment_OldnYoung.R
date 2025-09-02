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
###============
# Load necessary libraries
library(dplyr)
library(eulerr)
library(GeneOverlap)

###============
# NOTE Loading Targets from bulk RNA
# Define transcription factors
TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs

# Initialize a list to store intersected genes for each TF
bulkRNA_genes_list <- list()

# Loop through each TF
for (cur_tf in TFnames) {

  ###============
  # Load data from Bulk RNA-seq
  bulkRNA <- read.delim(paste0(data_in, "20240830OverexpressionbulkRNAseq/", "GFP_vs_", cur_tf, ".txt"))

  ###============
  # Subset and process data
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
  deg_bulkRNA <- bulkRNA_sub$gene

  bulkRNA_genes_list[[cur_tf]] <- deg_bulkRNA

}

# Print the final list of intersected genes for each TF
print(bulkRNA_genes_list)

names(bulkRNA_genes_list)


###============
# NOTE Loading Targets from Young
# Define transcription factors
TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs

# Initialize a list to store intersected genes for each TF
ODC_Young_genes_list <- list()

# Loop through each TF
for (cur_tf in TFnames) {

  # Load data
  ODC_Young_targets <- read.csv(paste0(tfPATH, "ODC_Young_targets_", cur_tf, ".csv"))
  bulkRNA <- read.delim(paste0(data_in, "20240830OverexpressionbulkRNAseq/", "GFP_vs_", cur_tf, ".txt"))

  # Remove missing targets not present in bulkRNA
  missing_targets <- ODC_Young_targets$target[!(ODC_Young_targets$target %in% bulkRNA$Feature)]
  ODC_Young_targets <- ODC_Young_targets %>% filter(!(target %in% missing_targets))

  # Extract primary and secondary genes
  primary_genes <- ODC_Young_targets %>% subset(depth == "1") %>% .$target
  secondary_genes <- subset(ODC_Young_targets, depth == "2") %>% .$target

  # Remove the current TF from the gene lists
  primary_genes <- setdiff(primary_genes, cur_tf)
  secondary_genes <- setdiff(secondary_genes, cur_tf)

  # Identify and remove duplicated genes
  duplicated_genes <- intersect(primary_genes, secondary_genes)

  if (length(duplicated_genes) > 0) {
    cat("Duplicated genes found for", cur_tf, ":\n")
    print(duplicated_genes)
  } else {
    cat("No duplicated genes found for", cur_tf, "\n")
  }

  # Ensure secondary genes don't include primary genes
  secondary_genes <- setdiff(secondary_genes, primary_genes)

  # Find intersected genes with deg_bulkRNA
  # collected_genes <- intersect(c(primary_genes, secondary_genes, cur_tf), deg_bulkRNA)
  collected_genes <- unique(c(primary_genes, secondary_genes, cur_tf))

  # Store in list with TF as key
  ODC_Young_genes_list[[cur_tf]] <- collected_genes
}

# Print the final list of intersected genes for each TF
print(ODC_Young_genes_list)



###============
# NOTE Loading Targets from Old
# Define transcription factors
TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs

# Initialize a list to store intersected genes for each TF
ODC_Old_genes_list <- list()

# Loop through each TF
for (cur_tf in TFnames) {

  # Load data
  ODC_Old_targets <- read.csv(paste0(tfPATH, "ODC_Old_targets_", cur_tf, ".csv"))
  bulkRNA <- read.delim(paste0(data_in, "20240830OverexpressionbulkRNAseq/", "GFP_vs_", cur_tf, ".txt"))

  # Remove missing targets not present in bulkRNA
  missing_targets <- ODC_Old_targets$target[!(ODC_Old_targets$target %in% bulkRNA$Feature)]
  ODC_Old_targets <- ODC_Old_targets %>% filter(!(target %in% missing_targets))

  # Extract primary and secondary genes
  primary_genes <- ODC_Old_targets %>% subset(depth == "1") %>% .$target
  secondary_genes <- subset(ODC_Old_targets, depth == "2") %>% .$target

  # Remove the current TF from the gene lists
  primary_genes <- setdiff(primary_genes, cur_tf)
  secondary_genes <- setdiff(secondary_genes, cur_tf)

  # Identify and remove duplicated genes
  duplicated_genes <- intersect(primary_genes, secondary_genes)

  if (length(duplicated_genes) > 0) {
    cat("Duplicated genes found for", cur_tf, ":\n")
    print(duplicated_genes)
  } else {
    cat("No duplicated genes found for", cur_tf, "\n")
  }

  # Ensure secondary genes don't include primary genes
  secondary_genes <- setdiff(secondary_genes, primary_genes)

  # Find intersected genes with deg_bulkRNA
  collected_genes <- unique(c(primary_genes, secondary_genes, cur_tf))

  # Store in list with TF as key
  ODC_Old_genes_list[[cur_tf]] <- collected_genes
}

# Print the final list of intersected genes for each TF
print(ODC_Old_genes_list)




# gene <- bulkRNA_genes_list[["Bach2"]]
# length(gene)



#
#===============================
# NOTE overlapping
# TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs


# Create a list
List1 <- list(
  Bach2_YoungTargets = ODC_Young_genes_list[["Bach2"]],
  Bhlhe41_YoungTargets = ODC_Young_genes_list[["Bhlhe41"]],
  Elf2_YoungTargets = ODC_Young_genes_list[["Elf2"]],
  Foxk2_YoungTargets = ODC_Young_genes_list[["Foxk2"]],
  Nr6a1_YoungTargets = ODC_Young_genes_list[["Nr6a1"]],
  Sox8_YoungTargets = ODC_Young_genes_list[["Sox8"]],
  Stat3_YoungTargets = ODC_Young_genes_list[["Stat3"]]
)


# Create a list
List2 <- list(
  Bach2_OldTargets = ODC_Old_genes_list[["Bach2"]],
  Bhlhe41_OldTargets = ODC_Old_genes_list[["Bhlhe41"]],
  Elf2_OldTargets = ODC_Old_genes_list[["Elf2"]],
  Foxk2_OldTargets = ODC_Old_genes_list[["Foxk2"]],
  Nr6a1_OldTargets = ODC_Old_genes_list[["Nr6a1"]],
  Sox8_OldTargets = ODC_Old_genes_list[["Sox8"]],
  Stat3_OldTargets = ODC_Old_genes_list[["Stat3"]]
)


List3 <- list(
  GFPvsBach2 = bulkRNA_genes_list[["Bach2"]],
  GFPvsBhlhe41 = bulkRNA_genes_list[["Bhlhe41"]],
  GFPvsElf2 = bulkRNA_genes_list[["Elf2"]],
  GFPvsFoxk2 = bulkRNA_genes_list[["Foxk2"]],
  GFPvsNr6a1 = bulkRNA_genes_list[["Nr6a1"]],
  GFPvsSox8 = bulkRNA_genes_list[["Sox8"]],
  GFPvsStat3 = bulkRNA_genes_list[["Stat3"]]
)

#,
# PiD_peakset = unique(Peakset$nearestGene))

# names(List2)
bulkRNA <- read.delim(paste0(data_in, "20240830OverexpressionbulkRNAseq/", "GFP_vs_Bach2.txt"))
nrow(bulkRNA)

List1.background <- unique(bulkRNA$Feature)

List2.background <- unique(bulkRNA$Feature)

List3.background <- unique(bulkRNA$Feature)



# #===================================
# # NOTE
# #===================================
##### DEX Enrichment
library(WGCNA)

setwd("/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/Analysis/bulkRNA/")

source('/dfs7/swaruplab/zechuas/Projects/PiD_2021/snRNA_data/Analysis/Figures/DEG_DAR/Code/ORA.R') ##in the codes folder

#==========================
# # NOTE List 1 and List 3

ORA.OR = matrix(NA,nrow=length(List1),ncol=length(List3));
      colnames(ORA.OR) = names(List3);
      rownames(ORA.OR) = names(List1);
      ORA.P = matrix(NA,nrow=length(List1),ncol=length(List3));
      colnames(ORA.P) = names(List3);
      rownames(ORA.P) = names(List1);

      for (i in 1:length(List3)) {
        for (j in 1:length(List1)) {
          result = ORA(List1[[j]],List3[[i]],List1.background,List3.background);
          ORA.OR[j,i] = result[1];
          ORA.P[j,i] = result[2];
  }
}

ORA.OR<-apply(ORA.OR,2,as.numeric)
dim(ORA.OR)<-dim(ORA.P)
head(ORA.OR)



FDRmat.Array <- matrix(p.adjust(ORA.P,method="fdr"),nrow=nrow(ORA.P),ncol=ncol(ORA.P))
rownames(FDRmat.Array)=rownames(ORA.P)
colnames(FDRmat.Array)=colnames(ORA.P)

ORA.P=matrix(as.numeric(ORA.P),nrow=nrow(ORA.P),ncol=ncol(ORA.P))
ORA.OR=matrix(as.numeric(ORA.OR),nrow=nrow(ORA.OR),ncol=ncol(ORA.OR))
rownames(ORA.P) <- rownames(ORA.OR) <- rownames(FDRmat.Array)
colnames(ORA.P) <- colnames(ORA.OR) <- colnames(FDRmat.Array)

dispMat <- ORA.OR ## You can change this to be just log2(Bmat) if you want the color to reflect the odds ratios
## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  ORA.OR
txtMat[FDRmat.Array >0.05] <- ""
txtMat[FDRmat.Array <0.05&FDRmat.Array >0.01] <- "*"
txtMat[FDRmat.Array <0.01&FDRmat.Array >0.005] <- "**"
txtMat[FDRmat.Array <0.005] <- "***"
# txtMat[FDRmat.Array >0.01] <- ""
# txtMat[FDRmat.Array <0.01&FDRmat.Array >0.005] <- "*"
# txtMat[FDRmat.Array <0.005&FDRmat.Array >0.001] <- "**"
# txtMat[FDRmat.Array <0.001] <- "***"

txtMat1 <- signif(ORA.OR,2)
txtMat1[txtMat1<1] <- "" # txtMat1<1.5
# txtMat1[txtMat1<2] <- ""

textMatrix1 = paste( txtMat1, '\n', txtMat , sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( ORA.P),nrow=nrow( ORA.P))


save(List1, List3, List1.background, List3.background, file= "ODC_Young_n_bulkRNA_enrichment_data.rda")

pdf("ODC_Young_n_bulkRNA.pdf", width=7,height=12)
labeledHeatmap(Matrix=dispMat,
		yLabels=rownames(dispMat),
		yColorLabels=TRUE,
		xLabels= colnames(dispMat),
		colors=viridis(50), # blueWhiteRed(40),
		textMatrix = textMatrix1,
		cex.lab.x=1.0,
		zlim=c(-0.1,3), # changed from zlim=c(-0.1,3),
		main="ODC Young and bulk RNA Enrichment Heatmap")
dev.off()




#==========================
# # NOTE List 2 and List 3

ORA.OR = matrix(NA,nrow=length(List2),ncol=length(List3));
      colnames(ORA.OR) = names(List3);
      rownames(ORA.OR) = names(List2);
      ORA.P = matrix(NA,nrow=length(List2),ncol=length(List3));
      colnames(ORA.P) = names(List3);
      rownames(ORA.P) = names(List2);

      for (i in 1:length(List3)) {
        for (j in 1:length(List2)) {
          result = ORA(List2[[j]],List3[[i]],List2.background,List3.background);
          ORA.OR[j,i] = result[1];
          ORA.P[j,i] = result[2];
  }
}

ORA.OR<-apply(ORA.OR,2,as.numeric)
dim(ORA.OR)<-dim(ORA.P)
head(ORA.OR)



FDRmat.Array <- matrix(p.adjust(ORA.P,method="fdr"),nrow=nrow(ORA.P),ncol=ncol(ORA.P))
rownames(FDRmat.Array)=rownames(ORA.P)
colnames(FDRmat.Array)=colnames(ORA.P)

ORA.P=matrix(as.numeric(ORA.P),nrow=nrow(ORA.P),ncol=ncol(ORA.P))
ORA.OR=matrix(as.numeric(ORA.OR),nrow=nrow(ORA.OR),ncol=ncol(ORA.OR))
rownames(ORA.P) <- rownames(ORA.OR) <- rownames(FDRmat.Array)
colnames(ORA.P) <- colnames(ORA.OR) <- colnames(FDRmat.Array)

dispMat <- ORA.OR ## You can change this to be just log2(Bmat) if you want the color to reflect the odds ratios
## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  ORA.OR
txtMat[FDRmat.Array >0.05] <- ""
txtMat[FDRmat.Array <0.05&FDRmat.Array >0.01] <- "*"
txtMat[FDRmat.Array <0.01&FDRmat.Array >0.005] <- "**"
txtMat[FDRmat.Array <0.005] <- "***"
# txtMat[FDRmat.Array >0.01] <- ""
# txtMat[FDRmat.Array <0.01&FDRmat.Array >0.005] <- "*"
# txtMat[FDRmat.Array <0.005&FDRmat.Array >0.001] <- "**"
# txtMat[FDRmat.Array <0.001] <- "***"

txtMat1 <- signif(ORA.OR,2)
txtMat1[txtMat1<1] <- "" # txtMat1<1.5
# txtMat1[txtMat1<2] <- ""

textMatrix1 = paste( txtMat1, '\n', txtMat , sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( ORA.P),nrow=nrow( ORA.P))


save(List2, List3, List2.background, List3.background, file= "ODC_Old_n_bulkRNA_enrichment_data.rda")

pdf("ODC_Old_n_bulkRNA.pdf", width=7,height=12)
labeledHeatmap(Matrix=dispMat,
		yLabels=rownames(dispMat),
		yColorLabels=TRUE,
		xLabels= colnames(dispMat),
		colors=viridis(50), # blueWhiteRed(40),
		textMatrix = textMatrix1,
		cex.lab.x=1.0,
		zlim=c(-0.1,3), # changed from zlim=c(-0.1,3),
		main="ODC Old and bulk RNA Enrichment Heatmap")
dev.off()
