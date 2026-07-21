#!/usr/bin/env bash
#
# regenerate_figures.sh -- regenerate the MOUSE spatial figures (supplementary).
# -----------------------------------------------------------------------------
# All plotting scripts read pre-computed tables from <data root>/data/ (staged
# from GEO, or produced by run_analysis.sh). Human figures are handled by
# ../human/regenerate_human.sh; the manuscript panel by ./regenerate_panel.sh.
#
# Config (override via environment):
#   SPATIAL_REPO_ROOT  data root that contains data/ + docs/  (default: repo root)
#   SPATIAL_PY         python interpreter                     (default: python)
#   SPATIAL_RSCRIPT    Rscript interpreter                    (default: Rscript)
#
# Usage:
#   micromamba activate xenium-processing
#   export SPATIAL_REPO_ROOT=/path/to/data-root
#   bash spatial/mouse/regenerate_figures.sh
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"       # spatial/mouse
SPATIAL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"                    # spatial/
DATA_ROOT="${SPATIAL_REPO_ROOT:-$(cd "${SPATIAL_DIR}/.." && pwd)}"
export SPATIAL_REPO_ROOT="${DATA_ROOT}"
cd "${DATA_ROOT}"

PY="${SPATIAL_PY:-python}"
RSCRIPT="${SPATIAL_RSCRIPT:-Rscript}"
export PYTHONPATH="${SPATIAL_DIR}/shared:${PYTHONPATH:-}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mplcache_xenium}"; mkdir -p "${MPLCONFIGDIR}"
export KMP_DUPLICATE_LIB_OK=TRUE

echo "Data root : ${DATA_ROOT}"
echo "Python    : ${PY}"
echo

# --- Cell-type composition: oligo-lineage differentiation (Young vs Old) -----
# 3 boxplots (OPC / Intermediate / Mature), within-lineage proportions Young vs Old.
"${PY}" "${SCRIPT_DIR}/composition/plot_lesion_lineage_proportions.py" --strategies stringent_p5
cp "data/lesion_celltype_proportions/figures/lineage_young_vs_old_stringent_p5.png" \
   "docs/images/lineage_young_vs_old_stringent_p5.png" 2>/dev/null || true

# --- Expression evidence ----------------------------------------------------
"${PY}" "${SCRIPT_DIR}/expression/plot_signal_to_noise.py"                # signal-to-noise boxplot
"${PY}" "${SCRIPT_DIR}/expression/plot_expression_aggregated_dotplot.py"  # aggregated dot plot (also panel C)

# --- Differential expression ------------------------------------------------
"${PY}" "${SCRIPT_DIR}/pseudobulk_de/plot_pca_celltypes_groups.py"        # PCA of lesion pseudobulks
"${RSCRIPT}" "${SCRIPT_DIR}/pseudobulk_de/plot_de_results.R" stringent_p5 outlier_rm counts  # per-TF barplots
"${RSCRIPT}" "${SCRIPT_DIR}/pseudobulk_de/plot_de_results.R" stringent_p5 outlier_rm vst     # volcanoes + log2FC heatmap

echo
echo "Mouse spatial figures regenerated under ${DATA_ROOT}/docs/images/."
