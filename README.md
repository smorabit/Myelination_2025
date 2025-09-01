
# Dysregulation of transcriptional networks regulating oligodendrogenesis in age-related decline in CNS remyelination

Dimas, Morabito, & Rawji et al. 2025 (TODO: Link to the manuscript here)


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
* [Perform optimal pathfinding analysis](integration/networks/TF_pathfinding.Rmd)

## Dynamically altered TFs in ageing and disease (Fig. 4)

* Mouse snRNA-seq DEGs 
* Human snRNA-seq DEGs 
* Prioritize TFs
* Plot TF Regulatory networks 

## Bulk RNA-seq (Fig. 5)