#!/usr/bin/env python3
"""Transcript-level expression evidence for the human Xenium TF panel.

Segmentation-free: counts high-quality (QV>=20) transcripts per gene straight from
transcripts.parquet, then compares each gene to the Xenium negative-control-probe
noise floor. Because cell segmentation on these sections is unreliable (weak DAPI),
we never assign transcripts to cells for the headline metric; we only use the
nucleus mask for an orthogonal "is this TF nuclear?" sanity check.

Outputs (under data/xenium_human/expression_evidence/):
  - human_tf_evidence.csv   one row per (sample x gene) for the full panel + two
                            negative-control pseudo-rows per sample. Columns:
                            sample_dir, sample_label, flagged, gene, gene_class,
                            is_positive_control, count, panel_rank, panel_pct,
                            neg_probe_mean, fold_over_neg, z_over_neg,
                            nuclear_fraction, median_nucleus_distance,
                            n_gene_transcripts
  - human_qc_summary.csv    per-sample QC metrics pulled from metrics_summary.csv,
                            for the report's segmentation section.

Run:
  micromamba activate xenium-processing
  HUMAN_TX_DIR=/path/to/human/transcripts python scripts/compute_human_tf_evidence.py
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


import numpy as np
import pandas as pd

import human_tf_common as C

# Metrics we surface for the segmentation/QC section of the report.
QC_FIELDS = [
    "total_high_quality_decoded_transcripts",
    "fraction_transcripts_decoded_q20",
    "num_cells_detected",
    "fraction_empty_cells",
    "fraction_transcripts_assigned",
    "median_transcripts_per_cell",
    "median_genes_per_cell",
    "nuclear_transcripts_per_100um2",
    "negative_control_probe_counts_per_control_per_cell",
    "region_area",
]


def qc_summary() -> pd.DataFrame:
    """One row per sample of the headline QC metrics from metrics_summary.csv."""
    rows = []
    for s in C.SAMPLES:
        m = pd.read_csv(C.metrics_path(s["dir"])).iloc[0]
        row = {"sample_label": s["label"], "slide": s["slide"], "flagged": s["flagged"]}
        for f in QC_FIELDS:
            row[f] = m.get(f, np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def evidence_for_sample(s: dict) -> pd.DataFrame:
    """Full-panel transcript-count evidence for one sample."""
    tx = C.load_transcripts(
        s["dir"],
        columns=[
            "feature_name",
            "qv",
            "is_gene",
            "codeword_category",
            "overlaps_nucleus",
            "nucleus_distance",
            "cell_id",
        ],
    )
    tx = tx[tx["qv"] >= C.QV_MIN]  # keep only high-quality molecules
    tx["assigned"] = tx["cell_id"] != "UNASSIGNED"  # landed inside a segmented cell

    # --- negative-control-probe noise floor (per-feature counts) ---
    ncp = tx[tx["codeword_category"] == C.NEG_PROBE]
    ncp_per_feature = ncp["feature_name"].value_counts()
    n_ncp_feat = int(ncp_per_feature.size)
    # mean/std across the individual probes; guard the dim sample where probes may be sparse
    neg_mean = float(ncp_per_feature.mean()) if n_ncp_feat else np.nan
    neg_std = float(ncp_per_feature.std(ddof=0)) if n_ncp_feat > 1 else np.nan

    ncc = tx[tx["codeword_category"] == C.NEG_CODEWORD]
    ncc_per_feature = ncc["feature_name"].value_counts()
    negcw_mean = float(ncc_per_feature.mean()) if ncc_per_feature.size else np.nan

    # --- per-gene counts over the real panel ---
    genes = tx[tx["is_gene"]]
    counts = genes["feature_name"].value_counts()  # descending
    # nuclear stats per gene (mask-based, segmentation-light)
    nuc = genes.groupby("feature_name", observed=True).agg(
        nuclear_fraction=("overlaps_nucleus", "mean"),
        median_nucleus_distance=("nucleus_distance", "median"),
        assigned_fraction=("assigned", "mean"),
    )

    n_genes = int(counts.size)
    ranks = counts.rank(ascending=False, method="min").astype(int)

    recs = []
    for gene, cnt in counts.items():
        cnt = int(cnt)
        fold = cnt / neg_mean if neg_mean and neg_mean > 0 else np.nan
        z = (cnt - neg_mean) / neg_std if neg_std and neg_std > 0 else np.nan
        recs.append(
            {
                "sample_dir": s["dir"],
                "sample_label": s["label"],
                "flagged": s["flagged"],
                "gene": gene,
                "gene_class": C.gene_class(gene),
                "is_positive_control": C.is_positive_control(gene),
                "count": cnt,
                "panel_rank": int(ranks[gene]),
                "panel_pct": round(100.0 * (n_genes - ranks[gene] + 1) / n_genes, 1),
                "neg_probe_mean": round(neg_mean, 3)
                if not np.isnan(neg_mean)
                else np.nan,
                "fold_over_neg": round(fold, 2) if not np.isnan(fold) else np.nan,
                "z_over_neg": round(z, 1) if not np.isnan(z) else np.nan,
                "nuclear_fraction": round(float(nuc.loc[gene, "nuclear_fraction"]), 3),
                "assigned_fraction": round(
                    float(nuc.loc[gene, "assigned_fraction"]), 3
                ),
                "median_nucleus_distance": round(
                    float(nuc.loc[gene, "median_nucleus_distance"]), 2
                ),
                "n_gene_transcripts": int(len(genes)),
            }
        )

    # negative-control pseudo-rows (per-feature mean count for the category)
    for cat, mean_val, nfeat in [
        (C.NEG_PROBE, neg_mean, n_ncp_feat),
        (C.NEG_CODEWORD, negcw_mean, int(ncc_per_feature.size)),
    ]:
        recs.append(
            {
                "sample_dir": s["dir"],
                "sample_label": s["label"],
                "flagged": s["flagged"],
                "gene": cat,
                "gene_class": "negative_control",
                "is_positive_control": False,
                "count": round(mean_val, 2) if not np.isnan(mean_val) else np.nan,
                "panel_rank": np.nan,
                "panel_pct": np.nan,
                "neg_probe_mean": round(neg_mean, 3)
                if not np.isnan(neg_mean)
                else np.nan,
                "fold_over_neg": 1.0 if cat == C.NEG_PROBE else np.nan,
                "z_over_neg": 0.0 if cat == C.NEG_PROBE else np.nan,
                "nuclear_fraction": np.nan,
                "median_nucleus_distance": np.nan,
                "n_gene_transcripts": int(len(genes)),
                "n_control_features": nfeat,
            }
        )

    return pd.DataFrame(recs)


def main() -> None:
    C.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    qc = qc_summary()
    qc_out = C.OUTPUT_DIR / "human_qc_summary.csv"
    qc.to_csv(qc_out, index=False)
    print(f"wrote {qc_out}  ({len(qc)} samples)")

    frames = [evidence_for_sample(s) for s in C.SAMPLES]
    ev = pd.concat(frames, ignore_index=True)
    ev_out = C.OUTPUT_DIR / "human_tf_evidence.csv"
    ev.to_csv(ev_out, index=False)
    print(f"wrote {ev_out}  ({len(ev)} rows)")

    # console sanity summary: primary TFs per sample, fold over neg floor
    tf = ev[ev["gene_class"] == "primary_tf"]
    print("\nPrimary TF fold-over-negative-probe (per sample):")
    pivot = tf.pivot_table(
        index="gene", columns="sample_label", values="fold_over_neg", observed=True
    )
    print(pivot.reindex(C.PRIMARY_TFS).round(1).to_string())


if __name__ == "__main__":
    main()
