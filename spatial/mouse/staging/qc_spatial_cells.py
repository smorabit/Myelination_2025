#!/usr/bin/env python3
"""Phase 0.1-0.2 per-sample cell QC for the spatial TF lesion analysis.

Loads each spinal-cord sample via the shared three-way join in
`scripts/spatial_common.py` (spatial AnnData x refined cell types x lesion-mask
membership, use_raw=True -> 50-gene panel) and produces per-sample QC tables and
distribution figures. The default cohort is COHORT_10 (all 10 samples) so that
Old 2 (the canonical outlier) is shown and its exclusion can
be RE-CONFIRMED here rather than silently re-derived.

Inputs (per sample, read by spatial_common.load_sample_cells):
  data/spatial_anndata/{sample}_spatial_with_annotations.h5ad
  data/celltype_refined/{sample}_celltype_refined_{strategy}.csv
  data/cells_in_lesions/output-XETG<serial>__{sample}__*_lesion[_check].csv

Outputs:
  data/spatial/spatial_qc_summary.csv        one row per sample (cell + lesion QC)
  data/spatial/spatial_qc_tf_detection.csv   per-TF detection fractions
  docs/images/spatial_qc_transcript_counts.{png,pdf}
  docs/images/spatial_qc_n_genes.{png,pdf}
  docs/images/spatial_qc_cell_area.{png,pdf}
  docs/images/spatial_qc_lesion_oligo_counts.{png,pdf}

Console: an "Old 2 re-confirmation" verdict comparing Old 2's in-lesion oligo
counts and total cell counts against the cohort (epistemically honest: states
ambiguity rather than overstating).

Run (project env; the micromamba run wrapper hits a lockfile so use the binary):
  KMP_DUPLICATE_LIB_OK=TRUE PYTHONPATH=scripts \
    python \
    scripts/qc_spatial_cells.py
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


import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

import cohort
from spatial_common import (
    COHORT_9,
    COHORT_10,
    OLIGO_TYPES,
    OUTLIER_SAMPLE,
    OUTPUT_DIR,
    PRIMARY_TFS,
    REPO_ROOT,
    age_for,
    load_sample_cells,
)

IMAGES_DIR = REPO_ROOT / "docs" / "images"
DEFAULT_STRATEGY = "stringent_p5"

# Shared palettes.
AGE_COLORS = {"Young": "#6a5acd", "Old": "#2e8b57"}  # slateblue / seagreen
CELL_TYPE_COLORS = {
    "OPC": "#b23aee",  # darkorchid2
    "Intermediate_Oligo": "#4169e1",  # royalblue
    "Mature_Oligo": "#00ff7f",  # springgreen
}
CELL_TYPE_LABELS = {
    "OPC": "OPC",
    "Intermediate_Oligo": "Intermediate",
    "Mature_Oligo": "Mature",
}

# Human-readable sample labels from the samplesheet (shared/cohort.py); no real IDs here.
SAMPLE_LABELS = cohort.sample_labels()

# QC metrics rendered as per-sample distribution figures.
DIST_METRICS = {
    "transcript_counts": ("Transcript counts per cell", True),  # (axis label, log y)
    "n_genes": ("Genes detected per cell", False),
    "cell_area": ("Cell area (um^2)", False),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--strategy",
        default=DEFAULT_STRATEGY,
        help=f"Refined-celltype strategy suffix (default: {DEFAULT_STRATEGY}).",
    )
    parser.add_argument(
        "--cohort",
        choices=["9", "10"],
        default="10",
        help=(
            "Cohort to QC: 9 (canonical, Old 2 dropped) or 10 (all; default) "
            "so Old 2 is shown and its exclusion can be re-confirmed."
        ),
    )
    return parser.parse_args()


def sample_label(sample_id: str) -> str:
    return SAMPLE_LABELS.get(sample_id, sample_id)


def order_cohort(cohort: list[str]) -> list[str]:
    """Order samples Young-then-Old, by their numeric label index."""

    def key(sid: str) -> tuple[int, int]:
        lab = sample_label(sid)
        age_rank = 0 if lab.startswith("Young") else 1
        try:
            idx = int(lab.split()[-1])
        except (ValueError, IndexError):
            idx = 99
        return (age_rank, idx)

    return sorted(cohort, key=key)


def detection_fraction(matrix: np.ndarray) -> np.ndarray:
    """Fraction of cells (rows) with count > 0 per gene (column)."""
    if matrix.shape[0] == 0:
        return np.full(matrix.shape[1], np.nan)
    return np.asarray((matrix > 0).sum(axis=0)).ravel() / matrix.shape[0]


def to_dense(x) -> np.ndarray:
    return x.toarray() if hasattr(x, "toarray") else np.asarray(x)


def build_summary_row(sample_id: str, adata) -> dict:
    obs = adata.obs
    required = [
        "transcript_counts",
        "n_genes",
        "cell_area",
        "nucleus_area",
        "in_lesion",
        "celltype_refined",
    ]
    missing = [c for c in required if c not in obs.columns]
    if missing:
        raise ValueError(f"{sample_id}: obs missing required columns {missing}")

    in_lesion = obs["in_lesion"].to_numpy(dtype=bool)
    les = obs.loc[in_lesion]
    les_ct = les["celltype_refined"].value_counts()

    return {
        "sample_id": sample_id,
        "sample_label": sample_label(sample_id),
        "age_group": age_for(sample_id),
        "is_outlier": sample_id == OUTLIER_SAMPLE,
        "n_cells_total": int(adata.n_obs),
        "n_in_lesion": int(in_lesion.sum()),
        "n_OPC_lesion": int(les_ct.get("OPC", 0)),
        "n_Intermediate_lesion": int(les_ct.get("Intermediate_Oligo", 0)),
        "n_Mature_lesion": int(les_ct.get("Mature_Oligo", 0)),
        "median_transcript_counts": float(obs["transcript_counts"].median()),
        "mean_transcript_counts": float(obs["transcript_counts"].mean()),
        "median_n_genes": float(obs["n_genes"].median()),
        "median_cell_area": float(obs["cell_area"].median()),
        "median_nucleus_area": float(obs["nucleus_area"].median()),
        "lesion_area_um2": adata.uns.get("lesion_area_um2"),
        "lesion_export_flagged": bool(adata.uns.get("lesion_export_flagged", False)),
        "n_lesion_unmatched": int(adata.uns.get("n_lesion_unmatched", 0)),
        "n_cells_unannotated": int(adata.uns.get("n_cells_unannotated", 0)),
    }


def build_tf_detection(sample_id: str, adata, tf_idx: dict[str, int]) -> list[dict]:
    """Per-TF detection fraction, section-wide and within in-lesion oligo cells."""
    X = to_dense(adata.X)
    in_lesion = adata.obs["in_lesion"].to_numpy(dtype=bool)
    is_oligo = adata.obs["celltype_refined"].isin(OLIGO_TYPES).to_numpy(dtype=bool)
    les_oligo = in_lesion & is_oligo

    rows = []
    for tf in PRIMARY_TFS:
        col = tf_idx[tf]
        all_counts = X[:, col]
        det_all = float((all_counts > 0).mean()) if adata.n_obs else np.nan
        if les_oligo.sum() > 0:
            lo_counts = X[les_oligo, col]
            det_lo = float((lo_counts > 0).mean())
            mean_lo = float(lo_counts.mean())
        else:
            det_lo = np.nan
            mean_lo = np.nan
        rows.append(
            {
                "sample_id": sample_id,
                "sample_label": sample_label(sample_id),
                "age_group": age_for(sample_id),
                "is_outlier": sample_id == OUTLIER_SAMPLE,
                "tf": tf,
                "n_cells_total": int(adata.n_obs),
                "n_lesion_oligo": int(les_oligo.sum()),
                "detection_frac_all": det_all,
                "mean_counts_all": float(all_counts.mean()) if adata.n_obs else np.nan,
                "detection_frac_lesion_oligo": det_lo,
                "mean_counts_lesion_oligo": mean_lo,
            }
        )
    return rows


def _label_with_flag(sid: str) -> str:
    lab = sample_label(sid)
    return f"{lab} *" if sid == OUTLIER_SAMPLE else lab


def plot_distribution(
    per_sample: dict[str, np.ndarray],
    order: list[str],
    metric: str,
    axis_label: str,
    log_y: bool,
) -> None:
    """One box per sample (Young-then-Old), coloured by age; Old 2 flagged."""
    data = [per_sample[sid] for sid in order]
    fig, ax = plt.subplots(figsize=(max(7, 0.85 * len(order)), 4.5))
    bp = ax.boxplot(
        data,
        widths=0.6,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": "black", "linewidth": 1.2},
    )
    for sid, box in zip(order, bp["boxes"]):
        box.set_facecolor(AGE_COLORS[age_for(sid)])
        box.set_alpha(0.75)
        if sid == OUTLIER_SAMPLE:
            box.set_edgecolor("red")
            box.set_linewidth(2.2)
            box.set_hatch("///")
        else:
            box.set_edgecolor("black")
            box.set_linewidth(0.8)

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(
        [_label_with_flag(sid) for sid in order], rotation=45, ha="right", fontsize=9
    )
    for tick, sid in zip(ax.get_xticklabels(), order):
        if sid == OUTLIER_SAMPLE:
            tick.set_color("red")
            tick.set_fontweight("bold")
    if log_y:
        ax.set_yscale("log")
    ax.set_ylabel(axis_label, fontsize=10)
    ax.set_title(f"Per-sample {axis_label} (* = Old 2, canonical outlier)", fontsize=11)
    legend_handles = [
        Patch(facecolor=AGE_COLORS["Young"], alpha=0.75, label="Young"),
        Patch(facecolor=AGE_COLORS["Old"], alpha=0.75, label="Old"),
        Patch(facecolor="white", edgecolor="red", hatch="///", label="Old 2 (outlier)"),
    ]
    ax.legend(handles=legend_handles, fontsize=8, loc="best")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()

    base = IMAGES_DIR / f"spatial_qc_{metric}"
    fig.savefig(f"{base}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {base.relative_to(REPO_ROOT)}.png + .pdf")


def plot_lesion_oligo_counts(summary: pd.DataFrame, order: list[str]) -> None:
    """Grouped bars of in-lesion oligo counts per sample per type (counts are small)."""
    summary = summary.set_index("sample_id").reindex(order)
    count_cols = {
        "OPC": "n_OPC_lesion",
        "Intermediate_Oligo": "n_Intermediate_lesion",
        "Mature_Oligo": "n_Mature_lesion",
    }
    n = len(order)
    x = np.arange(n)
    width = 0.26

    fig, ax = plt.subplots(figsize=(max(7, 0.95 * n), 4.5))
    for k, (ct, col) in enumerate(count_cols.items()):
        vals = summary[col].to_numpy()
        ax.bar(
            x + (k - 1) * width,
            vals,
            width=width,
            color=CELL_TYPE_COLORS[ct],
            edgecolor="black",
            linewidth=0.5,
            label=CELL_TYPE_LABELS[ct],
        )
        for xi, v in zip(x + (k - 1) * width, vals):
            if not np.isnan(v):
                ax.text(xi, v, f"{int(v)}", ha="center", va="bottom", fontsize=6)

    # Mark Old 2 with a red span behind its group.
    for xi, sid in zip(x, order):
        if sid == OUTLIER_SAMPLE:
            ax.axvspan(xi - 0.5, xi + 0.5, color="red", alpha=0.08, zorder=0)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [_label_with_flag(sid) for sid in order], rotation=45, ha="right", fontsize=9
    )
    for tick, sid in zip(ax.get_xticklabels(), order):
        if sid == OUTLIER_SAMPLE:
            tick.set_color("red")
            tick.set_fontweight("bold")
    ax.set_ylabel("In-lesion oligo-lineage cells", fontsize=10)
    ax.set_title(
        "In-lesion oligo-lineage cell counts per sample (* = Old 2, outlier)",
        fontsize=11,
    )
    ax.legend(fontsize=8, title="Cell type")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()

    base = IMAGES_DIR / "spatial_qc_lesion_oligo_counts"
    fig.savefig(f"{base}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {base.relative_to(REPO_ROOT)}.png + .pdf")


def old2_verdict(summary: pd.DataFrame) -> None:
    """Print an epistemically honest Old 2 re-confirmation verdict."""
    print("\n=== Old 2 re-confirmation ===")
    if OUTLIER_SAMPLE not in set(summary["sample_id"]):
        print(
            "  Old 2 is NOT in this cohort (run with --cohort 10 to re-confirm). "
            "Skipping verdict."
        )
        return

    old2 = summary.set_index("sample_id").loc[OUTLIER_SAMPLE]
    others = summary[summary["sample_id"] != OUTLIER_SAMPLE]
    old_others = others[others["age_group"] == "Old"]

    old2_oligo = (
        old2["n_OPC_lesion"] + old2["n_Intermediate_lesion"] + old2["n_Mature_lesion"]
    )

    def span(s: pd.Series) -> str:
        return f"min={s.min():.0f}, median={s.median():.0f}, max={s.max():.0f}"

    other_oligo = (
        others["n_OPC_lesion"]
        + others["n_Intermediate_lesion"]
        + others["n_Mature_lesion"]
    )
    old_other_oligo = (
        old_others["n_OPC_lesion"]
        + old_others["n_Intermediate_lesion"]
        + old_others["n_Mature_lesion"]
    )

    print(f"  Old 2 total cells          : {int(old2['n_cells_total'])}")
    print(f"    cohort (other 9)         : {span(others['n_cells_total'])}")
    print(f"  Old 2 in-lesion cells      : {int(old2['n_in_lesion'])}")
    print(f"    cohort (other 9)         : {span(others['n_in_lesion'])}")
    print(
        f"  Old 2 in-lesion oligo cells: {int(old2_oligo)} "
        f"(OPC={int(old2['n_OPC_lesion'])}, "
        f"Inter={int(old2['n_Intermediate_lesion'])}, "
        f"Mature={int(old2['n_Mature_lesion'])})"
    )
    print(f"    cohort (other 9)         : {span(other_oligo)}")
    print(f"    other Old (n=4)          : {span(old_other_oligo)}")
    print(f"  Old 2 lesion_export_flagged: {bool(old2['lesion_export_flagged'])}")
    print(f"  Old 2 n_lesion_unmatched   : {int(old2['n_lesion_unmatched'])}")

    # Honest verdict: is Old 2 an extreme low outlier on the oligo axis?
    is_min_oligo = old2_oligo <= other_oligo.min()
    below_old_min = old2_oligo < old_other_oligo.min()
    print("\n  Verdict:")
    if is_min_oligo and below_old_min:
        print(
            "    Old 2 has the FEWEST in-lesion oligo cells of all 10 samples and "
            "falls below the range of the other Old samples. On the cell-count axis "
            "the data SUPPORTS the canonical exclusion: Old 2 contributes too few "
            "in-lesion oligo cells for stable per-celltype spatial estimates."
        )
    elif is_min_oligo:
        print(
            "    Old 2 has the fewest in-lesion oligo cells overall but is within the "
            "spread of the other Old samples. The exclusion is DEFENSIBLE but the "
            "cell-count gap to the next-lowest sample is the deciding factor; inspect "
            "the figures before treating it as clear-cut."
        )
    else:
        print(
            "    On in-lesion oligo cell counts Old 2 is NOT the lowest sample. The "
            "cell-count axis alone does NOT cleanly justify the exclusion; the "
            "canonical drop likely rests on other evidence (e.g. expression/QC in the "
            "existing pseudobulk pipeline). Reported here, NOT re-derived. Treat as "
            "ambiguous on this axis."
        )
    print(
        "    Note: this QC re-confirms an existing canonical decision; it does not "
        "re-derive it. n=10 (with Old 2), single lesion per sample."
    )


def main() -> int:
    args = parse_args()
    cohort = COHORT_10 if args.cohort == "10" else COHORT_9
    order = order_cohort(cohort)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    IMAGES_DIR.mkdir(parents=True, exist_ok=True)

    print(
        f"Strategy: {args.strategy}   cohort: {len(order)} samples "
        f"({'all 10, Old 2 shown' if args.cohort == '10' else 'canonical 9'})\n"
    )

    summary_rows: list[dict] = []
    tf_rows: list[dict] = []
    dist_data: dict[str, dict[str, np.ndarray]] = {m: {} for m in DIST_METRICS}

    tf_idx: dict[str, int] | None = None
    for sid in order:
        adata = load_sample_cells(sid, strategy=args.strategy)
        if adata.n_obs == 0:
            raise ValueError(f"{sid}: zero annotated cells after join — cannot QC")

        # TF presence check (fail explicitly if any primary TF is absent).
        if tf_idx is None:
            var_names = list(adata.var_names)
            missing_tfs = [tf for tf in PRIMARY_TFS if tf not in var_names]
            if missing_tfs:
                raise ValueError(
                    f"{sid}: primary TFs missing from var_names: {missing_tfs}. "
                    f"Panel has {len(var_names)} genes; expected all 9 PRIMARY_TFS."
                )
            tf_idx = {tf: var_names.index(tf) for tf in PRIMARY_TFS}
        else:
            missing_tfs = [tf for tf in PRIMARY_TFS if tf not in set(adata.var_names)]
            if missing_tfs:
                raise ValueError(f"{sid}: primary TFs missing: {missing_tfs}")

        summary_rows.append(build_summary_row(sid, adata))
        tf_rows.extend(build_tf_detection(sid, adata, tf_idx))
        for metric in DIST_METRICS:
            dist_data[metric][sid] = adata.obs[metric].to_numpy(dtype=float)

        flag = " [OUTLIER]" if sid == OUTLIER_SAMPLE else ""
        print(
            f"  {sample_label(sid):<8} {sid:<22} n={adata.n_obs:>6} "
            f"in_lesion={int(adata.obs['in_lesion'].sum()):>5}{flag}"
        )

    summary = pd.DataFrame(summary_rows)
    tf_det = pd.DataFrame(tf_rows)

    summary_path = OUTPUT_DIR / "spatial_qc_summary.csv"
    tf_path = OUTPUT_DIR / "spatial_qc_tf_detection.csv"
    summary.to_csv(summary_path, index=False)
    tf_det.to_csv(tf_path, index=False)
    print(f"\nWrote {summary_path.relative_to(REPO_ROOT)} ({len(summary)} rows)")
    print(f"Wrote {tf_path.relative_to(REPO_ROOT)} ({len(tf_det)} rows)")

    # --- Figures ---
    print("\nFigures:")
    for metric, (axis_label, log_y) in DIST_METRICS.items():
        plot_distribution(dist_data[metric], order, metric, axis_label, log_y)
    plot_lesion_oligo_counts(summary, order)

    # --- TF detection summary (console) ---
    print("\n=== Primary-TF detection (mean across cohort) ===")
    tf_means = (
        tf_det.groupby("tf")[["detection_frac_all", "detection_frac_lesion_oligo"]]
        .mean()
        .reindex(PRIMARY_TFS)
    )
    print(f"  {'TF':<10} {'det_all':>9} {'det_lesion_oligo':>18}")
    for tf, row in tf_means.iterrows():
        lo = row["detection_frac_lesion_oligo"]
        lo_str = "  n/a" if np.isnan(lo) else f"{lo:>17.3f}"
        print(f"  {tf:<10} {row['detection_frac_all']:>9.3f} {lo_str}")

    old2_verdict(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
