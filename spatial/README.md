# Spatial transcriptomics (10x Xenium)

Code for the spinal-cord Xenium spatial analysis in **[Dimas, Morabito & Rawji *et al.* bioRxiv (2025)](https://www.biorxiv.org/content/10.1101/2025.11.14.688494v1)**. This is a curated subset of the working analysis repo: it keeps only the code behind figures that appear in the paper (main panels plus the mouse and human spatial supplementary figures).

The scripts are a Python + R mix, unlike the rest of this repository. As with the other modalities, file paths are relative to a data root that you set per machine (see below), and the data will be released on GEO on publication.

## Layout

```
spatial/
  run_analysis.sh          master run order (compute -> figures)
  shared/                  modules imported everywhere (plot style, geometry, human panel, cell typing)
  mouse/
    staging/               per-sample cell QC (from staged AnnData)
    composition/           oligo-lineage proportions Young-vs-Old (differentiation delay)
    pseudobulk_de/         pseudobulk build, DESeq2 young-vs-old, PCA + DE figures
    expression/            expression evidence + signal-to-noise, aggregated dot plot
    manuscript_panel.py    composite manuscript figure panels (A-D + suppl human maps)
    extract_panel_composition.py   gene-panel design tables
  human/                   human TF evidence, panel ranking, cross-PCF, spatial maps
```

Shared modules live in `shared/`; every script adds that directory to its path at
runtime, so scripts run from any depth without a `PYTHONPATH` set.

## Running

Set a data root and the two transcript directories, then run the master driver:

```bash
micromamba activate xenium-processing          # Python: scanpy, anndata, shapely, ...
Rscript -e 'renv::restore(prompt = FALSE)'     # R: DESeq2, mgcv, ggplot2, ...

export SPATIAL_REPO_ROOT=/path/to/data-root    # directory that holds data/ (+ docs/, Manuscript/)
export HUMAN_TX_DIR=/path/to/human/transcripts
export XENIUM_RAW_DIR=/path/to/mouse/xenium/bundle
bash spatial/run_analysis.sh
```

`SPATIAL_REPO_ROOT` is where inputs are read (`data/...`) and figures are written
(`docs/images/`, `Manuscript/manuscript_figures/`); it defaults to this repo root
if unset. Override the interpreters with `SPATIAL_PY` / `SPATIAL_RSCRIPT`. Each
figure driver can also run on its own once the compute tables are staged:
`bash spatial/mouse/regenerate_figures.sh`, `.../human/regenerate_human.sh`,
`.../mouse/regenerate_panel.sh`.

## Figures and the scripts behind them

**Manuscript panels (main).** [`mouse/manuscript_panel.py`](mouse/manuscript_panel.py) generates panels A (DAPI morphology), B (single-cell TF view), C (aggregated dot plot, via [`mouse/expression/plot_expression_aggregated_dotplot.py`](mouse/expression/plot_expression_aggregated_dotplot.py)), D (Sox8 DE barplot) and the supplementary human TF maps (via [`human/plot_human_tf_spatial_maps.py`](human/plot_human_tf_spatial_maps.py)). Gene-panel design tables come from [`mouse/extract_panel_composition.py`](mouse/extract_panel_composition.py).

**Mouse spatial (supplementary).**

- Cell-type annotation: the markers defining the three oligo-lineage cell types (OPC / Intermediate / Mature) are in [`shared/annotate_celltypes.py`](shared/annotate_celltypes.py) (a marker table, not a generated figure).
- Age-related delay in oligodendrocyte differentiation reproduced (3 boxplots, one per cell type, Young vs Old within-lineage proportions): [`mouse/composition/plot_lesion_lineage_proportions.py`](mouse/composition/plot_lesion_lineage_proportions.py) on [`compute_lesion_celltype_composition.py`](mouse/composition/compute_lesion_celltype_composition.py).
- 50-gene panel + secondary targets/inducers: [`mouse/extract_panel_composition.py`](mouse/extract_panel_composition.py).
- Genuine TF expression in day-14 lesions, signal-to-noise vs platform noise: [`mouse/expression/plot_signal_to_noise.py`](mouse/expression/plot_signal_to_noise.py) on [`compute_expression_evidence.py`](mouse/expression/compute_expression_evidence.py).
- Per-TF differential expression (Young vs Old) with downstream targets/inducers: barplots + per-cell-type volcanoes + cross-cell-type log2FC heatmap [`mouse/pseudobulk_de/plot_de_results.R`](mouse/pseudobulk_de/plot_de_results.R) and PCA [`mouse/pseudobulk_de/plot_pca_celltypes_groups.py`](mouse/pseudobulk_de/plot_pca_celltypes_groups.py), on [`de_young_vs_old_outlier_rm.R`](mouse/pseudobulk_de/de_young_vs_old_outlier_rm.R) + [`compute_combined_vst.R`](mouse/pseudobulk_de/compute_combined_vst.R) over [`build_pseudobulk.py`](mouse/pseudobulk_de/build_pseudobulk.py).

**Human spatial (supplementary).**

- Expressed above background (boxplot): [`human/plot_human_tf_signal_to_noise.py`](human/plot_human_tf_signal_to_noise.py) on [`compute_human_tf_evidence.py`](human/compute_human_tf_evidence.py).
- Where the factors rank across the panel: [`human/plot_human_panel_ranking.py`](human/plot_human_panel_ranking.py).
- Which cell types express each factor (Sox8 & Klk6 localize to oligo-lineage cells by proximity to cell-type markers): [`human/compute_cross_pcf.py`](human/compute_cross_pcf.py) (cross-PCF; computes and renders the cell-type localization heatmaps).

## Notes on this curation

- Scope: only code behind figures in the manuscript supplementary is included. Sample identifiers are not stored here; the cohort is read from an untracked `samplesheet.csv` (copy `shared/samplesheet.template.csv` and fill in real IDs), via [`shared/cohort.py`](shared/cohort.py).
- Paths in the individual script headers describe the original working-repo layout; the drivers above are the canonical entry points and set the environment for you.
