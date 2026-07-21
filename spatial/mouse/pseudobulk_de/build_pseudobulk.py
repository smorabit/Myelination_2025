#!/usr/bin/env python3
"""Build pseudobulk count matrices per (cell_type, strategy) for downstream DE.

For each strategy and each oligo-lineage cell type:
  1. For each sample, find cells that are (in the lesion mask) AND (of this
     cell type under the strategy's refined labels) AND with transcript_counts
     >= --per_cell_min_transcripts (default 10, per-cell QC floor).
  2. Subset the per-sample spatial AnnData to those cells, take adata.raw
     (all 50 panel genes), sum counts across cells -> one row per sample.
  3. Drop (sample, cell_type) combos with < FLOOR cells (after pcmt filter).

Outputs (per strategy, under data/pseudobulk/{strategy}/):
  {cell_type}_counts.csv     genes (rows) x samples (cols), integer counts
  {cell_type}_coldata.csv    sample, age_group, chip, n_cells, n_before_pcmt,
                             per_cell_min_transcripts
  _drop_log.csv              dropped (sample, cell_type) combos with reason

Canonical pipeline uses --per_cell_min_transcripts 10 per the 2026-05-11 spike
(plans/2026-05-11_SPIKE_tf-zero-expression-filtering_report.md). The pre-spike
behaviour (no per-cell QC) corresponds to --per_cell_min_transcripts 5, which
matches the cell-calling pipeline's min_mols_per_cell floor.
"""

from __future__ import annotations

# --- locate deposited module dirs (works from any depth) ------------------
import sys as _sys
from pathlib import Path as _Path

for _p in _Path(__file__).resolve().parents:
    if (_p / "shared" / "spatial_common.py").exists():
        for _sub in ("shared",):
            if (_p / _sub).is_dir() and str(_p / _sub) not in _sys.path:
                _sys.path.insert(0, str(_p / _sub))
        break
import cohort
from spatial_paths import repo_root


import argparse
import re
import sys

import anndata as ad
import numpy as np
import pandas as pd

REPO_ROOT = repo_root()
SPATIAL_DIR = REPO_ROOT / "data" / "spatial_anndata"
LESION_DIR = REPO_ROOT / "data" / "cells_in_lesions"
REFINED_DIR = REPO_ROOT / "data" / "celltype_refined"
OUTPUT_BASE = REPO_ROOT / "data" / "pseudobulk"

OLIGO_TYPES = ["OPC", "Intermediate_Oligo", "Mature_Oligo"]

# Cohort from the samplesheet (shared/cohort.py); no real sample IDs stored here.
YOUNG_SAMPLES = cohort.young_samples()
OLD_SAMPLES = cohort.old_samples()

LESION_FILE_RE = re.compile(
    r"output-XETG\w+__(\d+__R\d+_\d+)__\d+__\d+_lesion(?:_check)?\.csv$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--strategies",
        nargs="+",
        default=["stringent"],
        help="Strategies to process (default: stringent).",
    )
    parser.add_argument(
        "--floor",
        type=int,
        default=5,
        help="Drop (sample, cell_type) combos below this cell-count floor (default: 5).",
    )
    parser.add_argument(
        "--per_cell_min_transcripts",
        "--pcmt",
        type=int,
        default=10,
        help=(
            "Per-cell QC floor on transcript_counts (panel-gene total per cell). "
            "Cells below this are excluded before pseudobulk aggregation. "
            "Default 10 is the canonical value per the 2026-05-11 spike "
            "(plans/2026-05-11_SPIKE_tf-zero-expression-filtering_report.md); "
            "set to 5 to reproduce pre-migration pseudobulks."
        ),
    )
    return parser.parse_args()


def load_lesion_cell_ids() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for path in sorted(LESION_DIR.glob("output-XETG*__*_lesion*.csv")):
        m = LESION_FILE_RE.search(path.name)
        if not m:
            continue
        sample_id = m.group(1)
        df = pd.read_csv(path, comment="#")
        if "Cell ID" not in df.columns:
            continue
        ids = set(df["Cell ID"].astype(str))
        if not ids:
            continue
        out[sample_id] = ids
    return out


def age_for(sample_id: str) -> str:
    return cohort.age_for(sample_id)


def chip_for(sample_id: str) -> str:
    return sample_id.split("__")[0]


def aggregate(
    adata_raw: ad.AnnData,
    refined: pd.DataFrame,
    lesion_ids: set[str],
    cell_type: str,
    per_cell_min_transcripts: int,
) -> tuple[pd.Series, int, int]:
    refined_match = refined.loc[
        refined["celltype_refined"] == cell_type, "cell_id"
    ].astype(str)
    obs_names = adata_raw.obs_names.astype(str)
    in_lesion = obs_names.isin(lesion_ids)
    in_ct = obs_names.isin(set(refined_match))
    base_mask = in_lesion & in_ct
    n_before_pcmt = int(base_mask.sum())
    if "transcript_counts" not in adata_raw.obs.columns:
        raise SystemExit(
            "adata.obs is missing 'transcript_counts' column needed for per-cell QC filter"
        )
    transcript_counts = adata_raw.obs["transcript_counts"].to_numpy()
    pass_pcmt = transcript_counts >= per_cell_min_transcripts
    mask = base_mask & pass_pcmt
    n_cells = int(mask.sum())
    if n_cells == 0:
        return pd.Series(dtype=int), 0, n_before_pcmt
    sub = adata_raw[mask, :]
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    counts = pd.Series(X.sum(axis=0).astype(int), index=list(sub.var_names))
    return counts, n_cells, n_before_pcmt


def process_strategy(
    strategy: str,
    floor: int,
    per_cell_min_transcripts: int,
    lesion_ids_by_sample: dict[str, set[str]],
) -> None:
    out_dir = OUTPUT_BASE / strategy
    out_dir.mkdir(parents=True, exist_ok=True)

    per_celltype_counts: dict[str, dict[str, pd.Series]] = {
        ct: {} for ct in OLIGO_TYPES
    }
    per_celltype_coldata: dict[str, list[dict]] = {ct: [] for ct in OLIGO_TYPES}
    drops: list[dict] = []

    samples = sorted(YOUNG_SAMPLES | OLD_SAMPLES)
    print(f"  samples to process: {len(samples)}")

    for sample_id in samples:
        spatial_path = SPATIAL_DIR / f"{sample_id}_spatial_with_annotations.h5ad"
        refined_path = REFINED_DIR / f"{sample_id}_celltype_refined_{strategy}.csv"
        if not spatial_path.exists():
            print(f"    [skip] {sample_id}: spatial AnnData missing")
            continue
        if not refined_path.exists():
            print(f"    [skip] {sample_id}: refined CSV missing for {strategy}")
            continue
        if sample_id not in lesion_ids_by_sample:
            print(f"    [skip] {sample_id}: no lesion mask")
            continue

        adata = ad.read_h5ad(spatial_path)
        adata_raw = adata.raw.to_adata() if adata.raw is not None else adata
        refined = pd.read_csv(refined_path)
        lesion_ids = lesion_ids_by_sample[sample_id]

        for ct in OLIGO_TYPES:
            counts, n_cells, n_before_pcmt = aggregate(
                adata_raw, refined, lesion_ids, ct, per_cell_min_transcripts
            )
            cells_dropped_by_pcmt = n_before_pcmt - n_cells
            if cells_dropped_by_pcmt > 0:
                print(
                    f"    [pcmt] {sample_id} {ct}: kept {n_cells}/{n_before_pcmt} "
                    f"(dropped {cells_dropped_by_pcmt} with transcript_counts < {per_cell_min_transcripts})"
                )
            if n_cells < floor:
                drops.append(
                    {
                        "strategy": strategy,
                        "sample_id": sample_id,
                        "cell_type": ct,
                        "n_cells": n_cells,
                        "n_before_pcmt": n_before_pcmt,
                        "floor": floor,
                        "per_cell_min_transcripts": per_cell_min_transcripts,
                        "reason": f"below {floor}-cell floor",
                    }
                )
                continue
            per_celltype_counts[ct][sample_id] = counts
            per_celltype_coldata[ct].append(
                {
                    "sample_id": sample_id,
                    "age_group": age_for(sample_id),
                    "chip": chip_for(sample_id),
                    "n_cells": n_cells,
                    "n_before_pcmt": n_before_pcmt,
                    "per_cell_min_transcripts": per_cell_min_transcripts,
                }
            )

    print()
    for ct in OLIGO_TYPES:
        if not per_celltype_counts[ct]:
            print(f"  {ct}: NO usable samples — all dropped")
            continue
        counts_df = pd.DataFrame(per_celltype_counts[ct])
        # ensure consistent gene ordering
        counts_df = counts_df.sort_index()
        coldata_df = pd.DataFrame(per_celltype_coldata[ct]).set_index("sample_id")
        # align column order
        counts_df = counts_df[coldata_df.index]
        counts_path = out_dir / f"{ct}_counts.csv"
        coldata_path = out_dir / f"{ct}_coldata.csv"
        counts_df.to_csv(counts_path, index_label="gene")
        coldata_df.to_csv(coldata_path)
        n_y = (coldata_df["age_group"] == "Young").sum()
        n_o = (coldata_df["age_group"] == "Old").sum()
        print(
            f"  {ct}: {counts_df.shape[1]} samples ({n_y} Young + {n_o} Old), "
            f"{counts_df.shape[0]} genes -> {counts_path.relative_to(REPO_ROOT)}"
        )

    if drops:
        drops_df = pd.DataFrame(drops)
        drops_path = out_dir / "_drop_log.csv"
        drops_df.to_csv(drops_path, index=False)
        print(f"\n  Dropped {len(drops)} (sample, cell_type) combo(s):")
        print(drops_df.to_string(index=False))


def main() -> int:
    args = parse_args()
    lesion_ids_by_sample = load_lesion_cell_ids()
    if not lesion_ids_by_sample:
        raise SystemExit(f"No usable lesion CSVs in {LESION_DIR}")

    print(f"Lesion masks loaded for {len(lesion_ids_by_sample)} samples")
    print(f"Floor: ≥ {args.floor} cells per (sample, cell_type)")
    print(f"Per-cell QC: transcript_counts ≥ {args.per_cell_min_transcripts}")
    print(f"Strategies: {args.strategies}")

    for strategy in args.strategies:
        print(f"\n=== {strategy} ===")
        process_strategy(
            strategy, args.floor, args.per_cell_min_transcripts, lesion_ids_by_sample
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
