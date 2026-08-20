#!/usr/bin/env python3
"""Shared config + transcript loader for the human Xenium TF-expression analysis.

This is the human analog of ``spatial_common.py`` for the *mouse* study, but it is
deliberately **segmentation-free**: it reads the Xenium ``transcripts.parquet``
directly and never touches cell boundaries. Cell segmentation on these three human
sections is unreliable (weak DAPI -> median 1-3 transcripts/cell, 32-92% empty
cells, 60-73% of transcripts unassigned), so any cell-level clustering would fit
noise. Working at the transcript level keeps 100% of the molecules and their
locations.

Per-sample input: ``{TX_DIR}/{sample}_transcripts.parquet`` (staged repo layout)
or ``{TX_DIR}/{sample}/transcripts.parquet`` (raw Xenium output-dir layout).
Columns used downstream: ``feature_name, qv, x_location, y_location,
overlaps_nucleus, nucleus_distance, is_gene, codeword_category``.

The transcript source directory is configurable via ``HUMAN_TX_DIR`` so the same
scripts run wherever the human transcripts are staged:

    HUMAN_TX_DIR=/path/to/human/transcripts python compute_human_tf_evidence.py

Panel note: these samples use the human ``hBrain`` panel (316 detected genes). All
9 primary TFs, 28/29 mouse targets (NKX2-9 has no clean human ortholog and is
absent), and all 9 oligo-lineage markers are present (verified against the
transcript feature names).
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

import cohort
from spatial_paths import repo_root

REPO_ROOT = repo_root()
# Transcript/metrics source. Defaults to the repo staging dir; override with
# HUMAN_TX_DIR to point at the human transcripts directory.
TX_DIR = Path(
    os.environ.get("HUMAN_TX_DIR", REPO_ROOT / "data" / "xenium_human" / "transcripts")
)
OUTPUT_DIR = REPO_ROOT / "data" / "xenium_human" / "expression_evidence"

QV_MIN = 20  # standard Xenium high-quality (Phred Q20) threshold

# Xenium codeword_category values that are NOT biological genes.
NEG_PROBE = (
    "negative_control_probe"  # non-targeting probes -> non-specific binding floor
)
NEG_CODEWORD = "negative_control_codeword"  # unused codewords -> optical/decoding floor

# --- samples ---------------------------------------------------------------
# Sourced from the samplesheet (shared/cohort.py, gitignored samplesheet.csv); no real
# section identifiers are stored in this repository. Each entry has `dir` (the Xenium
# output directory name), `label` (short region tag), `slide`, and `flagged` (dim capture
# failure: image QC failed, ~34k transcripts, ~92% empty cells, ~11.5% negative-control).
SAMPLES = cohort.human_samples()

# --- gene groups (human symbols, verified present on the hBrain panel) ------
# 9 primary transcription factors (the hypothesis genes).
PRIMARY_TFS = [
    "BACH2",
    "ELF2",
    "FOXK2",
    "BHLHE41",
    "NR6A1",
    "SOX8",
    "STAT3",
    "SOX5",
    "KLK6",
]

# 28 downstream targets / inducers (NKX2-9 dropped: no human ortholog on panel).
TARGETS = [
    "ARHGEF12",
    "NR3C1",
    "FERMT2",
    "CENPB",
    "SLC39A3",
    "GPATCH4",
    "SAMD8",
    "AACS",
    "ABL2",
    "UGGT1",
    "SLC38A6",
    "TADA1",
    "HSPBAP1",
    "ANKS3",
    "BBS2",
    "NFIB",
    "NAA40",
    "SELENOH",
    "HSP90B1",
    "EIF1B",
    "NUP214",
    "TRRAP",
    "WDSUB1",
    "PKNOX2",
    "RORA",
    "ZEB1",
    "PLAG1",
    "MITF",
]

# 9 oligodendrocyte-lineage markers.
OLIGO_MARKERS = [
    "OLIG2",
    "SOX10",
    "PDGFRA",
    "PTPRZ1",
    "PCDH15",
    "ENPP6",
    "MBP",
    "OPALIN",
    "MOG",
]

# Positive-control reference: abundant brain/myelin + astrocyte genes that MUST be
# expressed if the assay worked. Xenium ships no positive-control *probe*, so these
# on-panel high-abundance genes are the positive anchor. Several overlap OLIGO_MARKERS
# (MBP, MOG) -- gene_class() keeps their marker class; is_positive_control() flags them.
POSITIVE_CONTROLS = [
    "MBP",
    "PLP1",
    "MOG",
    "MOBP",
    "MAG",
    "CLDN11",
    "MAL",
    "ERMN",
    "AQP4",
    "GJA1",
]

# Proliferation / cell-cycle genes: a biological "expected-low" reference, since brain
# is largely post-mitotic. The tight mitotic markers rank in the bottom ~13% of the
# panel. PCNA is the intended exception (it also does DNA repair in non-dividing cells,
# so it stays high) and is kept in the set to show the reference behaves like real biology.
PROLIFERATION_MARKERS = ["MKI67", "TOP2A", "CENPF", "CDK1", "CCNB2", "CCNA1", "PCNA"]


# Oligodendrocyte-lineage typing rules, mirroring the mouse stringent_p5 strategy:
# a cell must co-express the lineage gate (OLIG2 AND SOX10), then a stage marker.
# p5 drops MBP from the mature set (its mRNA diffuses far in Xenium) and applies
# precedence Intermediate > Mature > OPC.
LINEAGE_GATE = ["OLIG2", "SOX10"]
OPC_MARKERS = ["PDGFRA", "PTPRZ1", "PCDH15"]
INTERMEDIATE_MARKERS = ["ENPP6"]
MATURE_MARKERS = ["MOG", "OPALIN"]


def gene_class(gene: str) -> str:
    """Single-label class for a gene, by priority (TF > target > marker > other)."""
    if gene in PRIMARY_TFS:
        return "primary_tf"
    if gene in TARGETS:
        return "target"
    if gene in OLIGO_MARKERS:
        return "oligo_marker"
    return "other"


def is_positive_control(gene: str) -> bool:
    return gene in POSITIVE_CONTROLS


def _resolve(sample_dir: str, staged_name: str, nested_name: str) -> Path:
    """Resolve a per-sample file across the staged and raw output-dir layouts."""
    staged = TX_DIR / f"{sample_dir}_{staged_name}"
    if staged.exists():
        return staged
    nested = TX_DIR / sample_dir / nested_name
    if nested.exists():
        return nested
    raise FileNotFoundError(
        f"{nested_name} for {sample_dir} not found under {TX_DIR} "
        f"(looked for {staged.name} and {sample_dir}/{nested_name}). "
        f"Set HUMAN_TX_DIR or stage the data."
    )


def tx_path(sample_dir: str) -> Path:
    return _resolve(sample_dir, "transcripts.parquet", "transcripts.parquet")


def metrics_path(sample_dir: str) -> Path:
    return _resolve(sample_dir, "metrics_summary.csv", "metrics_summary.csv")


def cell_matrix_path(sample_dir: str) -> Path:
    """Path to the Xenium cell x feature matrix (cell_feature_matrix.h5)."""
    return _resolve(sample_dir, "cell_feature_matrix.h5", "cell_feature_matrix.h5")


def load_transcripts(sample_dir: str, columns: list[str] | None = None) -> pd.DataFrame:
    """Read one sample's transcripts.parquet (optionally a column subset)."""
    return pd.read_parquet(tx_path(sample_dir), columns=columns)


if __name__ == "__main__":
    # Smoke test: confirm each sample's parquet resolves and report row counts.
    print(f"TX_DIR: {TX_DIR}\n")
    print(f"{'label':<26} {'flagged':>7} {'parquet':>14}")
    for s in SAMPLES:
        try:
            p = tx_path(s["dir"])
            n = pd.read_parquet(p, columns=["transcript_id"]).shape[0]
            print(f"{s['label']:<26} {str(s['flagged']):>7} {n:>14,}")
        except FileNotFoundError as exc:
            print(f"{s['label']:<26} {str(s['flagged']):>7}   MISSING: {exc}")
