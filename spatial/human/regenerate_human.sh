#!/usr/bin/env bash
#
# regenerate_human.sh -- regenerate the HUMAN spatial figures (supplementary).
# -----------------------------------------------------------------------------
# Panel ranking + signal-to-noise read data/xenium_human/.../human_tf_evidence.csv
# (produced by compute_human_tf_evidence.py, run_analysis.sh Stage 4). The cross-PCF
# and spatial-map figures read the human transcripts under HUMAN_TX_DIR directly.
#
# Config (override via environment):
#   SPATIAL_REPO_ROOT  data root that contains data/ + docs/   (default: repo root)
#   HUMAN_TX_DIR       human Xenium transcripts directory       (REQUIRED)
#   SPATIAL_PY         python interpreter                       (default: python)
#
# Usage:
#   micromamba activate xenium-processing
#   export SPATIAL_REPO_ROOT=/path/to/data-root
#   export HUMAN_TX_DIR=/path/to/human/transcripts
#   bash spatial/human/regenerate_human.sh
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"       # spatial/human
SPATIAL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"                    # spatial/
DATA_ROOT="${SPATIAL_REPO_ROOT:-$(cd "${SPATIAL_DIR}/.." && pwd)}"
export SPATIAL_REPO_ROOT="${DATA_ROOT}"
cd "${DATA_ROOT}"

PY="${SPATIAL_PY:-python}"
export PYTHONPATH="${SPATIAL_DIR}/shared:${PYTHONPATH:-}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mplcache_xenium}"; mkdir -p "${MPLCONFIGDIR}"
export KMP_DUPLICATE_LIB_OK=TRUE
: "${HUMAN_TX_DIR:?set HUMAN_TX_DIR to the human Xenium transcripts directory}"

echo "Data root    : ${DATA_ROOT}"
echo "HUMAN_TX_DIR : ${HUMAN_TX_DIR}"
echo

"${PY}" "${SCRIPT_DIR}/plot_human_panel_ranking.py"          # where the factors rank across the panel
"${PY}" "${SCRIPT_DIR}/plot_human_tf_signal_to_noise.py"     # expressed above background (boxplot)
"${PY}" "${SCRIPT_DIR}/plot_human_tf_spatial_maps.py"        # TF spatial maps (per section)
# The cross-PCF cell-type localization heatmaps are computed+rendered by
# compute_cross_pcf.py in run_analysis.sh (Stage 4), not here.

echo
echo "Human spatial figures regenerated under ${DATA_ROOT}/docs/images/."
