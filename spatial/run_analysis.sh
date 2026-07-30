#!/usr/bin/env bash
#
# run_analysis.sh -- documented run order for the spinal-cord Xenium spatial
# analysis, from staged per-sample AnnData through the paper figures.
# -----------------------------------------------------------------------------
# Compute stages write intermediate tables under <data root>/data/; the figure
# drivers (Stage 6) render from those. Python steps need the analysis env
# (scanpy/anndata/shapely/...); R steps need the DESeq2/mgcv stack.
#
# Config (override via environment):
#   SPATIAL_REPO_ROOT  data root that contains data/ (+ docs/, Manuscript/)  (default: repo root)
#   HUMAN_TX_DIR       human Xenium transcripts directory      (needed Stage 4 + panel)
#   XENIUM_RAW_DIR     raw mouse Xenium bundle                 (needed for panel A/B)
#   SPATIAL_PY         python interpreter                      (default: python)
#   SPATIAL_RSCRIPT    Rscript interpreter                     (default: Rscript)
#
# Reproduce from a clean checkout after staging data/ and restoring the envs:
#   micromamba activate xenium-processing
#   Rscript -e 'renv::restore(prompt = FALSE)'
#   export SPATIAL_REPO_ROOT=/path/to/data-root HUMAN_TX_DIR=... XENIUM_RAW_DIR=...
#   bash spatial/run_analysis.sh
# -----------------------------------------------------------------------------
set -euo pipefail

SPATIAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"      # spatial/
DATA_ROOT="${SPATIAL_REPO_ROOT:-$(cd "${SPATIAL_DIR}/.." && pwd)}"
export SPATIAL_REPO_ROOT="${DATA_ROOT}"
cd "${DATA_ROOT}"

PY="${SPATIAL_PY:-python}"
RSCRIPT="${SPATIAL_RSCRIPT:-Rscript}"
export PYTHONPATH="${SPATIAL_DIR}/shared:${PYTHONPATH:-}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mplcache_xenium}"; mkdir -p "${MPLCONFIGDIR}"
export KMP_DUPLICATE_LIB_OK=TRUE
STRAT=stringent_p5
M="${SPATIAL_DIR}/mouse"
H="${SPATIAL_DIR}/human"

echo "== Stage 0: per-sample cell QC =="
# Inputs: staged per-sample AnnData under data/spatial_anndata/ (download the
# processed objects from GEO into that directory; see spatial/README.md).
"${PY}" "${M}/staging/qc_spatial_cells.py"                     # QC summary + figures

echo "== Stage 1: cell typing + pseudobulk + lesion composition (Python) =="
"${PY}" "${SPATIAL_DIR}/shared/annotate_celltypes.py" --strategy "${STRAT}"  # -> data/celltype_refined/{sample}_celltype_refined_<STRAT>.csv
"${PY}" "${M}/pseudobulk_de/build_pseudobulk.py" --strategies "${STRAT}" --per_cell_min_transcripts 10
"${PY}" "${M}/composition/compute_lesion_celltype_composition.py" --strategy "${STRAT}"  # oligo-lineage proportions

echo "== Stage 2: differential expression (R) =="
# Outlier sample id comes from the samplesheet (cohort.py), the single source of truth.
OUTLIER="$("${PY}" -c 'import cohort; print(cohort.outlier_sample())')"
"${RSCRIPT}" "${M}/pseudobulk_de/de_young_vs_old_outlier_rm.R" "${STRAT}"              # canonical DE (outlier removed)
"${RSCRIPT}" "${M}/pseudobulk_de/compute_combined_vst.R" "${STRAT}"            # full-cohort VST (PCA)
"${RSCRIPT}" "${M}/pseudobulk_de/compute_combined_vst.R" "${STRAT}" "${OUTLIER}"  # outlier-removed VST (volcanoes + log2FC heatmap)

echo "== Stage 3: expression evidence (Python) =="
"${PY}" "${M}/expression/compute_expression_evidence.py"

echo "== Stage 4: human TF evidence + cross-PCF (Python) =="
"${PY}" "${H}/compute_human_tf_evidence.py"
"${PY}" "${H}/compute_cross_pcf.py"                            # which cell types express each factor

echo "== Stage 5: figures =="
bash "${M}/regenerate_figures.sh"
bash "${H}/regenerate_human.sh"
bash "${M}/regenerate_panel.sh"

echo
echo "Done. Figures under ${DATA_ROOT}/docs/images/ and ${DATA_ROOT}/Manuscript/."
