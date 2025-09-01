
# Dysregulation of transcriptional networks regulating oligodendrogenesis in age-related decline in CNS remyelination

Dimas, Morabito, & Rawji et al. 2025 (TODO: Link to the manuscript here)

## Abstract

*In demyelinating diseases such as multiple sclerosis, the body's ability to replace damaged myelin sheaths—a process termed remyelination—is crucial for functional recovery. While remyelination is robust in young individuals, it becomes less efficient with age, and is linked to progressive disability. Despite its importance, the precise gene regulatory network that controls this complex process and how it fails with aging remains largely unknown. To address this gap, we created a comprehensive single-nucleus RNA and ATAC sequencing dataset of remyelination in young and aged mice using the lysolecithin model. Our integrated analysis allowed us to map the GRN that drives efficient remyelination by oligodendrocyte lineage cells in the CNS of young mice. We identified 10 dynamic gene modules with distinct biological functions and temporal expression patterns throughout the oligodendrocyte maturation trajectory. Importantly, our analysis revealed significant age-dependent changes in the activity of these modules, highlighting specific gene networks that are dysregulated with age.  We then pinpointed key transcription factors at critical nodes of this network that are dysregulated with age. We validated these findings by showing that similar TF dysregulations occur in human MS lesions by integrating our data with publicly available datasets. Functional validation in primary oligodendrocyte progenitor cells confirmed that altering the expression of these TFs significantly impacts OPC function and their ability to proliferate and differentiate. This study is the first to comprehensively map the GRN of remyelination and pinpoint TFs that act as central nodes and are key contributors to age-related remyelination failure. Our findings provide a foundational understanding of this critical regenerative process and age-related changes relevant to older patients and those with chronic MS who currently have limited therapeutic options.*

## About this repository

TODO: Briefly describe the layout of the repo.

## 🖥️ Data generated in this study 

Raw and processed snATAC-seq and snRNA-seq datasets have been deposited on the NCBI Gene Expression Omnibus (GEO). At this time, the data remains private and will be released upon publication of this study. A data access code has been provided to the peer reviewers of this study.

## 🧬 Processing raw sequencing data

* [Mapping and quantifying snATAC-seq data with `cellranger-atac count`](snATAC/preprocessing/cellranger-atac_count.sub)
* [Mapping and quantifying snRNA-seq data with `cellranger count`](snRNA/preprocessing/cellranger_count.sub)

## snATAC-seq + snRNA-seq clustering and integration (Fig. 1)

### 🌳 Clustering

* [snATAC-seq clustering analysis with ArchR](snATAC/clustering/snATAC_clustering.Rmd)
* [snRNA-seq clustering analysis with Seurat + Liger](snRNA/clustering/snRNA_clustering.Rmd)
* [snRNA-seq cluster marker genes](snRNA/clustering/cluster_markers.sub)

### ⚓ Integration

* [snATAC-seq + snRNA-seq integration with Seurat](integration/integration.Rmd)
* [Differential cell abundance analysis with MiloR](integration/diff_cell_abundance_miloR.Rmd)

### Miscellaneous plotting 

* TODO

## 🧐 Multi-omic oligodendrocyte trajectory analysis (Fig. 2)

* [snATAC + snRNA oligodendrocyte trajectory analysis with monocle3](integration/trajectory/trajectory_monocle3.Rmd)
* [Identifying trajectory DEGs with generalized linear models](integration/trajectory/trajectory_DEGs.Rmd)
* [Enhancer-promoter links via co-accessibility analysis](integration/networks/coaccessibility_cicero.Rmd)
* [Modeling trajectory expression dynamics using recurrent variational autoencoders (RVAEs)](integration/RVAE_modeling.ipynb)
* [Co-expression network analysis with hdWGCNA](integration/networks/hdWGCNA.Rmd)

## Transcription factor networks and optimal pathfinding (Fig. 3)

* [Construct TF regulatory networks](integration/networks/TF_network_construct.Rmd)
* [Define the top target genes (Regulons) for each TF](integration/networks/TF_network_regulons.Rmd)
* [Identify TF regulatory cascades via optimal pathfinding analysis](integration/networks/TF_pathfinding.Rmd)

## Dynamically altered TFs in ageing and disease (Fig. 4)

* Re-processing human snRNA-seq data
    * Pre-processing with Kallisto (TODO: get these files from the cross-disease repo).
    * Clustering analysis (TODO: Need to get the .ipynb files from the UCI HPC. Possibly already have them in the cross disease repo?).
* [Mouse snRNA-seq DEGs](snRNA/DEGs/) 
* [Human snRNA-seq DEGs](human/DEGs/)
* Prioritize TFs

## Bulk RNA-seq (Fig. 5)