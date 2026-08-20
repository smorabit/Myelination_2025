#!/usr/bin/env bash
#
# regenerate_panel.sh -- regenerate the manuscript figure panels (A-D + suppl).
# -----------------------------------------------------------------------------
# Runs manuscript_panel.py, which emits each component figure as an editable
# vector PDF (+ PNG) to <data root>/Manuscript/manuscript_figures/. Panels A-D
# are mouse; the supplementary human TF maps delegate to ../human/. The panel
# script assembles the mouse dot plot (mouse/expression) and the human spatial
# maps (human/) as imported modules.
#
# Config (override via environment):
#   SPATIAL_REPO_ROOT  data root with data/ + Manuscript/   (default: repo root)
#   XENIUM_RAW_DIR     raw mouse Xenium bundle (morphology) (REQUIRED for panel A/B)
#   HUMAN_TX_DIR       human transcripts (suppl maps)       (REQUIRED for suppl)
#   MANUSCRIPT_ROOT    output subdir under data root        (default: Manuscript)
#   SPATIAL_PY         python interpreter                   (default: python)
#
# Usage:
#   micromamba activate xenium-processing
#   export SPATIAL_REPO_ROOT=/path/to/data-root
#   export XENIUM_RAW_DIR=/path/to/mouse/xenium/bundle
#   export HUMAN_TX_DIR=/path/to/human/transcripts
#   bash spatial/mouse/regenerate_panel.sh
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"       # spatial/mouse
SPATIAL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"                    # spatial/
DATA_ROOT="${SPATIAL_REPO_ROOT:-$(cd "${SPATIAL_DIR}/.." && pwd)}"
export SPATIAL_REPO_ROOT="${DATA_ROOT}"
cd "${DATA_ROOT}"

PY="${SPATIAL_PY:-python}"
export PYTHONPATH="${SPATIAL_DIR}/shared:${PYTHONPATH:-}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mplcache_xenium}"; mkdir -p "${MPLCONFIGDIR}"
export KMP_DUPLICATE_LIB_OK=TRUE
: "${XENIUM_RAW_DIR:?set XENIUM_RAW_DIR to the raw mouse Xenium bundle}"
: "${HUMAN_TX_DIR:?set HUMAN_TX_DIR to the human Xenium transcripts directory}"

echo "Data root      : ${DATA_ROOT}"
echo "XENIUM_RAW_DIR : ${XENIUM_RAW_DIR}"
echo "HUMAN_TX_DIR   : ${HUMAN_TX_DIR}"
echo

"${PY}" "${SCRIPT_DIR}/manuscript_panel.py" all --tf-set top9

echo
echo "Done. Figures in ${DATA_ROOT}/${MANUSCRIPT_ROOT:-Manuscript}/manuscript_figures/"
