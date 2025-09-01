
# Dysregulation of transcriptional networks regulating oligodendrogenesis in age-related decline in CNS remyelination

Dimas, Morabito, & Rawji et al. 2025 (TODO: Link to the manuscript here)


## Data generated in this study 

Raw and processed snATAC-seq and snRNA-seq datasets have been deposited on the NCBI Gene Expression Omnibus (GEO). At this time, the data remains private and will be released upon publication of this study. A data access code has been provided to the peer reviewers of this study.

## Processing raw sequencing data

* [Mapping and quantifying snATAC-seq data with `cellranger-atac count`](snATAC/preprocessing/cellranger-atac_count.sub)
* [Mapping and quantifying snRNA-seq data with `cellranger count`](snRNA/preprocessing/cellranger_count.sub)

## snATAC-seq + snRNA-seq clustering and integration (Fig. 1)

### Clustering

* [snATAC-seq clustering analysis with ArchR](snATAC/clustering/snATAC_clustering.Rmd)
* [snRNA-seq clustering analysis with Seurat + Liger](snRNA/clustering/snRNA_clustering.Rmd)

### Integration

* [snATAC-seq + snRNA-seq integration with Seurat](integration/integration.Rmd)
* [Differential cell abundance analysis with MiloR](integration/diff_cell_abundance_miloR.Rmd)


## Multi-omic oligodendrocyte trajectory analysis (Fig. 2)

* Defining the trajectory 
* Trajectory DEGs 
* Co-accessibility 
* Co-expression network analysis 
