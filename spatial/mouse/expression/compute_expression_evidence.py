#!/usr/bin/env python3
"""Compute per-(sample × cell_type × gene) expression evidence for all 50 panel genes.

Independent of the DE analysis: produces a descriptive long-form table of
expression magnitude (pseudobulk CPM, log2(CPM+1)) AND per-cell detection
prevalence (% cells with count ≥ 1) for every panel gene in every (sample ×
cell_type) pseudobulk that entered the canonical DE.

Cell population: `in_lesion ∩ celltype_refined ∈ {OPC, Intermediate_Oligo,
Mature_Oligo} ∩ transcript_counts ≥ 10` (the canonical pseudobulk population
under `outlier_rm + pcmt=10`). All 10 samples included — Old 2's sparser
pseudobulks appear transparently in lower detection fractions.

Outputs (under data/expression_evidence/):
  expression_per_sample_celltype_<strategy>.csv
    sample_id, sample_label, age_group, sample_is_old2, cell_type, gene,
    gene_class, n_cells, n_cells_detected, mean_cpm, detection_fraction,
    log2_cpm_plus1
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
from scipy.sparse import issparse

REPO_ROOT = repo_root()
SPATIAL_DIR = REPO_ROOT / "data" / "spatial_anndata"
REFINED_DIR = REPO_ROOT / "data" / "celltype_refined"
LESION_DIR = REPO_ROOT / "data" / "cells_in_lesions"
OUT_DIR = REPO_ROOT / "data" / "expression_evidence"

OLIGO_TYPES = ["OPC", "Intermediate_Oligo", "Mature_Oligo"]
DEFAULT_STRATEGY = "stringent_p5"
DEFAULT_PCMT = 10
# Cohort labels from the samplesheet (shared/cohort.py); no real sample IDs stored here.
OLD2_ID = cohort.outlier_sample()
SAMPLE_LABELS = cohort.sample_labels()

GENE_CLASS = {
    "Olig2": "marker",
    "Sox10": "marker",
    "Pdgfra": "marker",
    "Ptprz1": "marker",
    "Pcdh15": "marker",
    "Enpp6": "marker",
    "Mbp": "marker",
    "Opalin": "marker",
    "Mog": "marker",
    "Bach2": "primary_TF",
    "Elf2": "primary_TF",
    "Foxk2": "primary_TF",
    "Bhlhe41": "primary_TF",
    "Nr6a1": "primary_TF",
    "Sox8": "primary_TF",
    "Stat3": "primary_TF",
    "Sox5": "primary_TF",
    "Klk6": "primary_TF",
    "Arhgef12": "Bach2_target",
    "Nr3c1": "Bach2_target",
    "Fermt2": "Bach2_target",
    "Cenpb": "Bach2_target",
    "Slc39a3": "Elf2_target",
    "Gpatch4": "Elf2_target",
    "Samd8": "Elf2_target",
    "Aacs": "Foxk2_target",
    "Abl2": "Foxk2_target",
    "Uggt1": "Foxk2_target",
    "Slc38a6": "Bhlhe41_target",
    "Tada1": "Bhlhe41_target",
    "Hspbap1": "Bhlhe41_target",
    "Anks3": "Nr6a1_target",
    "Bbs2": "Nr6a1_target",
    "Nfib": "Nr6a1_target",
    "Naa40": "Nr6a1_target",
    "Selenoh": "Sox8_target",
    "Hsp90b1": "Sox8_target",
    "Eif1b": "Sox8_target",
    "Nup214": "Stat3_target",
    "Trrap": "Stat3_target",
    "Wdsub1": "Stat3_target",
    "Pknox2": "Sox5_inducer",
    "Rora": "Sox5_inducer",
    "Zeb1": "Sox5_inducer",
    "Plag1": "Klk6_inducer",
    "Nkx2-9": "Klk6_inducer",
    "Mitf": "Klk6_inducer",
    "Cpa1": "pancreas_control",
    "Spink1": "pancreas_control",
    "Nupr1": "pancreas_or_stress",
}

LESION_FILE_RE = re.compile(
    r"output-XETG\w+__(\d+__R\d+_\d+)__\d+__\d+_lesion(?:_check)?\.csv$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--strategy",
        default=DEFAULT_STRATEGY,
        help=f"Cell-type strategy (default: {DEFAULT_STRATEGY}).",
    )
    parser.add_argument(
        "--per_cell_min_transcripts",
        "--pcmt",
        type=int,
        default=DEFAULT_PCMT,
        help=f"Per-cell transcript_counts floor (default: {DEFAULT_PCMT}).",
    )
    return parser.parse_args()


def load_lesion_cells(sample_id: str) -> set[str]:
    files = list(LESION_DIR.glob(f"output-XETG*__{sample_id}__*_lesion*.csv"))
    cells: set[str] = set()
    for f in files:
        df = pd.read_csv(f, comment="#")
        cells |= set(df["Cell ID"].astype(str))
    return cells


def load_refined(sample_id: str, strategy: str) -> pd.DataFrame:
    p = REFINED_DIR / f"{sample_id}_celltype_refined_{strategy}.csv"
    return pd.read_csv(p)


def age_for(sample_id: str) -> str:
    return cohort.age_for(sample_id)


def process_sample(sample_id: str, strategy: str, pcmt: int) -> list[dict]:
    """Return one row per (cell_type, gene) for this sample."""
    adata = ad.read_h5ad(SPATIAL_DIR / f"{sample_id}_spatial_with_annotations.h5ad")
    obs = adata.obs.copy().reset_index()
    if obs.columns[0] != "cell_id":
        obs = obs.rename(columns={obs.columns[0]: "cell_id"})
    obs["cell_id"] = obs["cell_id"].astype(str)
    refined = load_refined(sample_id, strategy)
    obs = obs.merge(refined[["cell_id", "celltype_refined"]], on="cell_id", how="left")
    obs["in_lesion"] = obs["cell_id"].isin(load_lesion_cells(sample_id))

    raw = adata.raw
    raw_genes = list(raw.var.index)
    X = raw.X.toarray() if issparse(raw.X) else np.asarray(raw.X)

    rows: list[dict] = []
    for ct in OLIGO_TYPES:
        mask = (
            (obs["celltype_refined"] == ct)
            & obs["in_lesion"]
            & (obs["transcript_counts"] >= pcmt)
        ).to_numpy()
        n_cells = int(mask.sum())
        if n_cells == 0:
            for gene in raw_genes:
                rows.append(
                    {
                        "sample_id": sample_id,
                        "sample_label": SAMPLE_LABELS.get(sample_id, sample_id),
                        "age_group": age_for(sample_id),
                        "sample_is_old2": sample_id == OLD2_ID,
                        "cell_type": ct,
                        "gene": gene,
                        "gene_class": GENE_CLASS.get(gene, "other"),
                        "n_cells": 0,
                        "n_cells_detected": 0,
                        "mean_cpm": float("nan"),
                        "mean_counts_per_cell": float("nan"),
                        "detection_fraction": float("nan"),
                        "log2_cpm_plus1": float("nan"),
                        "log2_counts_plus1": float("nan"),
                    }
                )
            continue

        sub = X[mask]  # (n_cells, n_genes)
        sum_per_gene = sub.sum(axis=0)  # vector length n_genes
        library_size = float(sum_per_gene.sum())
        if library_size == 0:
            mean_cpm_per_gene = np.zeros(len(raw_genes))
        else:
            mean_cpm_per_gene = sum_per_gene / library_size * 1e6

        n_detected_per_gene = (sub >= 1).sum(axis=0).astype(int)

        for j, gene in enumerate(raw_genes):
            cpm = float(mean_cpm_per_gene[j])
            sum_gene = float(sum_per_gene[j])
            counts_per_cell = sum_gene / n_cells
            n_det = int(n_detected_per_gene[j])
            rows.append(
                {
                    "sample_id": sample_id,
                    "sample_label": SAMPLE_LABELS.get(sample_id, sample_id),
                    "age_group": age_for(sample_id),
                    "sample_is_old2": sample_id == OLD2_ID,
                    "cell_type": ct,
                    "gene": gene,
                    "gene_class": GENE_CLASS.get(gene, "other"),
                    "n_cells": n_cells,
                    "n_cells_detected": n_det,
                    "mean_cpm": cpm,
                    "mean_counts_per_cell": counts_per_cell,
                    "detection_fraction": n_det / n_cells,
                    "log2_cpm_plus1": float(np.log2(cpm + 1)),
                    "log2_counts_plus1": float(np.log2(counts_per_cell + 1)),
                }
            )

        # --- Negative-control pseudo-rows ---
        # Aggregate per-cell totals (sum across all probes in each Xenium
        # negative-control category, computed over the same in-lesion +
        # cell-type-restricted cells as the real genes). Reported as
        # mean_counts_per_cell = category_total_per_cell / n_cells.
        # mean_cpm and log2_cpm_plus1 are NaN — CPM is not meaningful for
        # an aggregate-of-many-probes pseudo-gene.
        for nc_col, nc_label in (
            ("control_probe_counts", "control_probe"),
            ("genomic_control_counts", "genomic_control"),
            ("unassigned_codeword_counts", "unassigned_codeword"),
        ):
            if nc_col not in obs.columns:
                continue
            nc_vals = obs.loc[mask, nc_col].to_numpy()
            nc_sum = float(nc_vals.sum())
            nc_mean_per_cell = nc_sum / n_cells
            nc_n_det = int((nc_vals >= 1).sum())
            rows.append(
                {
                    "sample_id": sample_id,
                    "sample_label": SAMPLE_LABELS.get(sample_id, sample_id),
                    "age_group": age_for(sample_id),
                    "sample_is_old2": sample_id == OLD2_ID,
                    "cell_type": ct,
                    "gene": nc_label,
                    "gene_class": "negative_control",
                    "n_cells": n_cells,
                    "n_cells_detected": nc_n_det,
                    "mean_cpm": float("nan"),
                    "mean_counts_per_cell": nc_mean_per_cell,
                    "detection_fraction": nc_n_det / n_cells,
                    "log2_cpm_plus1": float("nan"),
                    "log2_counts_plus1": float(np.log2(nc_mean_per_cell + 1)),
                }
            )
    return rows


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_samples = sorted(SAMPLE_LABELS.keys())
    print(f"Strategy : {args.strategy}")
    print(f"pcmt     : {args.per_cell_min_transcripts}")
    print(f"Samples  : {len(all_samples)}")
    print()

    all_rows: list[dict] = []
    for sid in all_samples:
        print(f"  processing {sid} ({SAMPLE_LABELS.get(sid, sid)})")
        all_rows.extend(
            process_sample(sid, args.strategy, args.per_cell_min_transcripts)
        )

    df = pd.DataFrame(all_rows)
    out_path = OUT_DIR / f"expression_per_sample_celltype_{args.strategy}.csv"
    df.to_csv(out_path, index=False)
    print()
    print(f"Wrote {len(df)} rows -> {out_path.relative_to(REPO_ROOT)}")

    # Quick sanity printout
    print(
        "\n=== n_cells per (sample × cell_type) — should match canonical pseudobulks ==="
    )
    n_table = (
        df.drop_duplicates(["sample_id", "cell_type"])
        .pivot(index="sample_id", columns="cell_type", values="n_cells")
        .reindex(columns=OLIGO_TYPES)
        .reindex(index=all_samples)
    )
    n_table.insert(0, "label", [SAMPLE_LABELS.get(s, s) for s in n_table.index])
    print(n_table.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
