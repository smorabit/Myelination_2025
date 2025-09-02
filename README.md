
# Dysregulation of transcriptional networks regulating oligodendrogenesis in age-related decline in CNS remyelination

Dimas, Morabito, & Rawji et al. 2025 (TODO: Link to the manuscript here)

## Abstract

*In demyelinating diseases like multiple sclerosis, efficient remyelination is critical for functional recovery. Remyelination loses efficiency with age, and is linked to progressive disability. The gene regulatory network underlying remyelination, and how it is altered with aging, remains unclear. Here we present a comparative single-nucleus RNA and ATAC sequencing analysis of remyelination in young and aged mice. We identified gene modules dynamically expressed throughout oligodendrocyte maturation, revealing age-dependent changes in key processes related to myelination. Multi-omic analysis allowed us to map the regulatory network driving efficient remyelination within oligodendrocyte lineage cells in young mice. We highlight key transcription factors in the network dysregulated with age, and we describe similar TF dysregulations in human MS lesions. Modifying the expression of these TFs in primary oligodendrocyte progenitor cells impacts proliferation and differentiation. These findings provide a foundational understanding of this regenerative process in the context of ageing and in chronic demyelinating diseases.*

## About this repository

This repository contains the code used for the data analysis associated with **Dimas, Morabito, Rawji *et al.* bioRxiv (2025)**. The sections of the `README` below follow the same flow as the mannuscript, linking to the relevant script(s) for each of the data analysis steps. 

These scripts can be used to reproduce the results starting from the raw sequencing reads downloaded from GEO, or from the processed data objects. Please note that the file paths included in these scripts are relative to the UCI HPC3 or the CNAG cluster, and these paths must be updated in order to run these scripts on a different machine. 

## 🖥️ Data generated in this study 

Raw and processed snATAC-seq and snRNA-seq datasets have been deposited on the NCBI Gene Expression Omnibus (GEO). At this time, the data remains private and will be released upon publication of this study. A data access code has been provided to the peer reviewers of this study.

## 🧬 Pre-processing raw sequencing data

snATAC-seq and snRNA-seq were performed **separately** in mouse spinal cord sections using 10X Genomics Chromium. Pre-processing (i.e. mapping, quantification, etc) of the sequencing data from these assays was performed using `cellranger-atac count` for snATAC-seq and `cellranger count` for snRNA-seq. 

* [Pre-processing snATAC-seq data with `cellranger-atac count`](snATAC/preprocessing/cellranger-atac_count.sub)
* [Pre-processing snRNA-seq data with `cellranger count`](snRNA/preprocessing/cellranger_count.sub)

## 🗺️ snATAC-seq + snRNA-seq clustering and integration (Fig. 1)

### 🌳 Clustering

We first performed clustering analysis for ATAC and RNA separately. This clustering analysis includes a batch correction step to correct for the effect of different sequencing batches, using **[harmony](https://github.com/immunogenomics/harmony)** for snRNA-seq and **[Liger](https://github.com/welch-lab/liger)** for snATAC-seq.

* [snATAC-seq clustering analysis with ArchR + Harmony](snATAC/clustering/snATAC_clustering.Rmd)
* [snRNA-seq clustering analysis with Seurat + Liger](snRNA/clustering/snRNA_clustering.Rmd)
* [snRNA-seq cluster marker genes](snRNA/clustering/cluster_markers.sub)

### ⚓ Integration

We next integrated snATAC-seq and snRNA-seq into a joint cellular embedding using Seurat. 

* [snATAC-seq + snRNA-seq integration with Seurat](integration/integration.Rmd)
* [Differential cell abundance analysis with MiloR](integration/diff_cell_abundance_miloR.Rmd)


## 🧐 Multi-omic oligodendrocyte trajectory analysis (Fig. 2)

Here we perform a multi-omic pseudotime trajectory analysis using monocle3 to study the transcriptional and epigenomic dynamics of the oligodendrocyte lineage in remyelination.

* [snATAC + snRNA oligodendrocyte trajectory analysis with monocle3](integration/trajectory/trajectory_monocle3.Rmd)
* [Identifying trajectory DEGs with generalized linear models](integration/trajectory/trajectory_DEGs.Rmd)
* [Modeling trajectory expression dynamics using recurrent variational autoencoders (RVAEs)](integration/RVAE_modeling.ipynb)

We used **[hdWGCNA](https://smorabit.github.io/hdWGCNA/)** to construct a gene co-expression network in the oligodendrocyte lineage, identifying 10 gene modules capturing key gene expression patterns across the trajectory. We also approximated epigenomic activity of these modules using the snATAC-seq dataset.

* [Co-expression network analysis with hdWGCNA](integration/networks/hdWGCNA.Rmd)

We next used **[cicero](https://cole-trapnell-lab.github.io/cicero-release/)** to identify linkages between promoter regions and non-coding cis-regulatory elements (putative enhancers), and we compared the strength of these linkages across ageing.

* [Enhancer-promoter links via co-accessibility analysis](integration/networks/coaccessibility_cicero.Rmd)

## 🕸️ Transcription factor networks and optimal pathfinding (Fig. 3)

We leverage epigenomic and transcriptomic measurements to construct the gene-regulatory network (GRN) underlying the oligodendrocyte lineage throughout remyelination in young and aged mice. We leveraged TF motif presence in promoter elements and linked enhancer elements from co-accessibility analysis in order to define directed links between TFs and downstream target genes. These links are assembled into a full GRN, and XGBoost regression was used to prioritize the most likely TF-gene links and to define TF Regulons.

* [Construct TF regulatory networks](integration/networks/TF_network_construct.Rmd)
* [Define the top target genes (Regulons) for each TF](integration/networks/TF_network_regulons.Rmd)

We sought to identify the processive regulatory events, or *cascades* of TF-TF activation that underly the cellular transitions from progenitors to mature oligodendrocytes. We propose a unique computational strategy to accomplish this goal. In the young and aged GRNs, we performed *optimal pathfinding analysis*, using the strongest directed TF-TF links to find the most likely path from "start" TFs (expressed in progenitors) to "end" TFs (expressed in mature).

* [Identify TF regulatory cascades via optimal pathfinding analysis](integration/networks/TF_pathfinding.Rmd)

## ⭐ Dynamically altered TFs in ageing and disease (Fig. 4)

We next sought to prioritize TFs that are altered in ageing mice, and altered in human multiple sclerosis (MS). First, we collected and re-processed two snRNA-seq datasets of human MS.

* [Preprocessing the Schirmer *et al.* 2019 dataset](human/Schirmer/preprocessing/)
* [Preprocessing the Absinta *et al.* 2021 dataset](human/Absinta/preprocessing/)
* [Clustering the Schirmer *et al.* 2019 dataset](human/Schirmer/process_Schirmer.ipynb)
* [Clustering the Absinta *et al.* 2021 dataset](human/Absinta/preprocessing/)

Next we performed differential expression analysis across ageing in the mouse dataset, and between MS tissue vs. healthy in the human datasets.

* [Differential expression analysis with ageing in mice](snRNA/DEGs/ageing_DEGs.sub)
* [Differential expression analysis in the human datasets](human/DEGs/)

Finally, we synthesize this information to prioritize TFs for further study with downstream experiments. 

* [Prioritize TFs for downstream experiments, additional plotting](integration/prioritize_TFs.Rmd)

## ⚡ Bulk RNA-seq in cell models (Fig. 5)

We performed bulk RNA-seq experiments in cell models where we over-expressed TFs of interest (Bach2, Bhlhe41, Elf2, Foxk2, Nr6a1, Sox8, and Stat3).

* [TODO: Bulk RNA-seq data processing](bulkRNA/processing/)
* [TODO: Bulk RNA-seq differential expression analysis](bulkRNA/DEGs/)
* [Bulk RNA-seq downstream plotting and comparison with single-nucleus data](bulkRNA/downstream/)

## Miscellaneous plotting and helper scripts

Here are scripts that were used for additional plotting, and scripts containing helper functions.

* [Miscellaneous plotting](misc/plotting_supplementary.Rmd)
* [Helper scripts](misc/helpers/)
