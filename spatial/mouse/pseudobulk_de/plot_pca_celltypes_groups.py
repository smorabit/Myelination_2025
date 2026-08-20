#!/usr/bin/env python3
"""PCA over the 50-gene panel × 30 (sample × cell-type) lesion pseudobulk profiles.

Defaults to the combined-celltype DESeq2 VST (`data/de_results/<strategy>/_combined_vst.csv`,
written by `scripts/compute_combined_vst.R`) so all 30 profiles share a single
normalisation scale (avoids the per-celltype-fit clamping artefact that affects
P5 OPC under the per-celltype VST). Falls back to log2(CPM+1) on the raw
pseudobulk counts if the combined VST file isn't present.

The figure has three panels:
  - Top-left: PC1 vs PC2 scatter, point color = cell type, shape = age group,
    samples labelled.
  - Top-right: PC1 vs PC3 scatter (same encodings).
  - Bottom: scree (proportion of variance explained per PC, first 10 PCs).

Inputs : data/de_results/<strategy>/_combined_vst.csv  (preferred)
         data/pseudobulk/<strategy>/<cell_type>_counts.csv  (fallback)
         data/pseudobulk/<strategy>/<cell_type>_coldata.csv
Output : docs/images/pca_celltypes_groups_<strategy>.png
         data/diagnostics/pca_celltypes_groups_<strategy>.csv  (PC scores + loadings)
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
from spatial_paths import repo_root


import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import plot_style

REPO_ROOT = repo_root()
DE_DIR = REPO_ROOT / "data" / "de_results"
PB_DIR = REPO_ROOT / "data" / "pseudobulk"
DIAG_DIR = REPO_ROOT / "data" / "diagnostics"
FIG_DIR = REPO_ROOT / "docs" / "images"

CELL_TYPES = ("OPC", "Intermediate_Oligo", "Mature_Oligo")
CELLTYPE_COLORS = {
    "OPC": plot_style.CELLTYPE_COLORS["OPC"],
    "Intermediate_Oligo": plot_style.CELLTYPE_COLORS["Intermediate"],
    "Mature_Oligo": plot_style.CELLTYPE_COLORS["Mature"],
}
GROUP_MARKERS = {"Young": "o", "Old": "^"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--strategy",
        default="stringent_p5",
        help="Strategy directory (default: stringent_p5; the primary analysis).",
    )
    p.add_argument(
        "--method",
        choices=["combined_vst", "log_cpm"],
        default="combined_vst",
        help=(
            "Input expression values for PCA. 'combined_vst' (default): read "
            "_combined_vst.csv from data/de_results/<strategy>/. 'log_cpm': "
            "compute log2(CPM+1) per profile from raw pseudobulk counts."
        ),
    )
    return p.parse_args()


def load_combined_vst(strategy: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (expr: 50 × 30, profile_meta: 30 × {sample, cell_type, age_group})."""
    vst_path = DE_DIR / strategy / "_combined_vst.csv"
    if not vst_path.exists():
        raise SystemExit(
            f"Combined VST CSV not found: {vst_path}.\n"
            f"Run `Rscript scripts/compute_combined_vst.R {strategy}` first."
        )
    vst = pd.read_csv(vst_path).set_index("gene")
    # Build sample → age_group lookup from any per-celltype coldata
    age_lookup = {}
    for ct in CELL_TYPES:
        coldata = pd.read_csv(PB_DIR / strategy / f"{ct}_coldata.csv")
        for _, row in coldata.iterrows():
            age_lookup[row["sample_id"]] = row["age_group"]

    profile_meta = []
    new_cols = []
    for col in vst.columns:
        if "__" not in col:
            continue
        sample_id, ct = col.rsplit("__", 1)
        if ct not in CELL_TYPES:
            continue
        prefix = f"{ct}."
        if sample_id.startswith(prefix):
            sample_id = sample_id[len(prefix) :]
        new_cols.append((col, sample_id, ct))
        profile_meta.append(
            {
                "profile_col": col,
                "sample_id": sample_id,
                "cell_type": ct,
                "age_group": age_lookup.get(sample_id, "Unknown"),
            }
        )
    meta_df = pd.DataFrame(profile_meta).set_index("profile_col")
    return vst, meta_df


def load_log_cpm(strategy: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute log2(CPM+1) per profile on raw pseudobulk counts."""
    pieces = []
    profile_meta = []
    for ct in CELL_TYPES:
        counts_path = PB_DIR / strategy / f"{ct}_counts.csv"
        coldata_path = PB_DIR / strategy / f"{ct}_coldata.csv"
        if not counts_path.exists() or not coldata_path.exists():
            raise SystemExit(f"Missing inputs for {ct} under {strategy}")
        counts = pd.read_csv(counts_path).set_index("gene")
        coldata = pd.read_csv(coldata_path).set_index("sample_id")
        lib = counts.sum(axis=0).replace(0, np.nan)
        cpm = counts.div(lib, axis=1) * 1e6
        log_cpm = np.log2(cpm + 1)
        # Disambiguate column names with cell type suffix
        new_cols = [f"{c}__{ct}" for c in log_cpm.columns]
        log_cpm.columns = new_cols
        pieces.append(log_cpm)
        for sid in coldata.index:
            profile_meta.append(
                {
                    "profile_col": f"{sid}__{ct}",
                    "sample_id": sid,
                    "cell_type": ct,
                    "age_group": coldata.at[sid, "age_group"],
                }
            )
    expr = pd.concat(pieces, axis=1)
    meta_df = pd.DataFrame(profile_meta).set_index("profile_col")
    return expr, meta_df


def main() -> int:
    args = parse_args()
    plot_style.apply()
    strategy = args.strategy
    method = args.method
    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    if method == "combined_vst":
        expr, meta = load_combined_vst(strategy)
        method_label = "DESeq2 combined VST"
    else:
        expr, meta = load_log_cpm(strategy)
        method_label = "log₂(CPM + 1)"

    # Align meta to expr columns
    meta = meta.reindex(expr.columns)
    if meta["sample_id"].isna().any():
        raise SystemExit("Some profile columns lack metadata.")

    n_genes = expr.shape[0]
    n_profiles = expr.shape[1]
    print(f"Strategy: {strategy}    method: {method_label}")
    print(f"Input matrix: {n_genes} genes × {n_profiles} (sample × cell-type) profiles")
    print()

    # PCA — samples are rows, genes are features. Z-score per gene first
    # (StandardScaler on transposed matrix). PCA on the scaled matrix.
    X = expr.T.to_numpy()  # 30 × 50
    scaler = StandardScaler(with_mean=True, with_std=True)
    X_scaled = scaler.fit_transform(X)
    n_components = min(10, n_profiles, n_genes)
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(X_scaled)
    var_explained = pca.explained_variance_ratio_

    # Save PC scores + meta. Build score_df from meta (already aligned to
    # expr.columns) so profile_col / sample_id / cell_type / age_group are
    # always present.
    score_df = meta.copy()
    score_df.insert(0, "profile_col", score_df.index)
    for i in range(n_components):
        score_df[f"PC{i + 1}"] = scores[:, i]
    score_df = score_df.reset_index(drop=True)
    score_df.to_csv(DIAG_DIR / f"pca_celltypes_groups_{strategy}.csv", index=False)

    # Loadings (gene × PC) — useful for inspection
    loadings = pd.DataFrame(
        pca.components_.T,
        index=expr.index,
        columns=[f"PC{i + 1}" for i in range(n_components)],
    )
    loadings_path = DIAG_DIR / f"pca_loadings_{strategy}.csv"
    loadings.to_csv(loadings_path)

    # Plot
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.2], hspace=0.35, wspace=0.28)
    ax_pc12 = fig.add_subplot(gs[0, 0])
    ax_pc13 = fig.add_subplot(gs[0, 1])
    ax_scree = fig.add_subplot(gs[1, :])

    def _plot_scatter(ax, x_idx, y_idx):
        for ct in CELL_TYPES:
            for grp in ("Young", "Old"):
                mask = (meta["cell_type"] == ct) & (meta["age_group"] == grp)
                if not mask.any():
                    continue
                idx = meta.index[mask]
                xs = score_df.loc[score_df["profile_col"].isin(idx), f"PC{x_idx + 1}"]
                ys = score_df.loc[score_df["profile_col"].isin(idx), f"PC{y_idx + 1}"]
                ax.scatter(
                    xs,
                    ys,
                    s=140,
                    c=CELLTYPE_COLORS[ct],
                    marker=GROUP_MARKERS[grp],
                    edgecolor="white",
                    linewidth=1.0,
                    alpha=0.92,
                    zorder=3,
                )
        # Sample labels
        for _, row in score_df.iterrows():
            sid_short = (
                row["sample_id"].split("__")[1]
                if "__" in row["sample_id"]
                else row["sample_id"]
            )
            ax.annotate(
                sid_short,
                xy=(row[f"PC{x_idx + 1}"], row[f"PC{y_idx + 1}"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=6,
                alpha=0.65,
                zorder=4,
            )
        ax.axhline(0, color="grey", lw=0.5, alpha=0.4, zorder=1)
        ax.axvline(0, color="grey", lw=0.5, alpha=0.4, zorder=1)
        ax.set_xlabel(
            f"PC{x_idx + 1}  ({var_explained[x_idx] * 100:.1f}%)", fontsize=10
        )
        ax.set_ylabel(
            f"PC{y_idx + 1}  ({var_explained[y_idx] * 100:.1f}%)", fontsize=10
        )
        ax.grid(True, linestyle=":", alpha=0.3, zorder=0)

    _plot_scatter(ax_pc12, 0, 1)
    ax_pc12.set_title("PC1 vs PC2", fontsize=11)
    _plot_scatter(ax_pc13, 0, 2)
    ax_pc13.set_title("PC1 vs PC3", fontsize=11)

    # Scree (first 10 PCs or fewer)
    ax_scree.bar(
        range(1, n_components + 1),
        var_explained * 100,
        color=plot_style.PALETTE[4],
        edgecolor="white",
    )
    ax_scree.set_xlabel("Principal component", fontsize=10)
    ax_scree.set_ylabel("Variance explained (%)", fontsize=10)
    ax_scree.set_title("Scree: variance explained per PC", fontsize=11)
    for i, v in enumerate(var_explained * 100):
        ax_scree.text(i + 1, v + 0.5, f"{v:.1f}%", ha="center", va="bottom", fontsize=8)
    ax_scree.set_xticks(range(1, n_components + 1))
    ax_scree.grid(axis="y", linestyle=":", alpha=0.3)

    # Combined legend (cell type + age group)
    legend_handles = []
    for ct in CELL_TYPES:
        legend_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="s",
                linestyle="",
                markersize=11,
                markerfacecolor=CELLTYPE_COLORS[ct],
                markeredgecolor="white",
                label=ct.replace("_", " "),
            )
        )
    for grp, marker in GROUP_MARKERS.items():
        legend_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker=marker,
                linestyle="",
                markersize=11,
                markerfacecolor="grey",
                markeredgecolor="white",
                label=grp,
            )
        )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=5,
        bbox_to_anchor=(0.5, 0.99),
        frameon=False,
        fontsize=9,
    )

    fig.suptitle(
        f"PCA: {n_profiles} lesion pseudobulks, {n_genes}-gene panel ({strategy})",
        y=0.94,
    )

    fig_path = FIG_DIR / f"pca_celltypes_groups_{strategy}.png"
    plot_style.save(fig, f"pca_celltypes_groups_{strategy}")

    print("Variance explained (first 10 PCs):")
    for i, v in enumerate(var_explained * 100):
        print(f"  PC{i + 1:>2}: {v:5.2f}%")
    print()
    print(f"Wrote {fig_path.relative_to(REPO_ROOT)}")
    print(
        f"Wrote {(DIAG_DIR / f'pca_celltypes_groups_{strategy}.csv').relative_to(REPO_ROOT)}"
    )
    print(f"Wrote {loadings_path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
