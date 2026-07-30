#!/usr/bin/env python3
"""Signal-to-noise companion figure for the TF expression-evidence argument.

Reads `data/expression_evidence/expression_per_sample_celltype_<strategy>.csv`
(produced by `scripts/compute_expression_evidence.py`, which must already have
been extended to include negative-control pseudo-rows) and renders a single
box+strip figure comparing the distribution of mean counts per cell across 4
biological/technical categories:

  Negative controls (Xenium platform: control_probe, genomic_control,
                     unassigned_codeword — pooled)
  Targets / inducers (29 secondary hypothesis genes)
  Primary TFs (9 hypothesis TFs)
  Markers (9 cell-type-defining genes)

Each category gets one boxplot summarising mean_counts_per_cell across the
30 (sample × cell_type) pseudobulks × N_genes-in-category individual entries.
Individual entries are overlaid as jittered points so a reader can see the
spread within each category.

Y-axis is log-scaled with a small positive epsilon offset so zero values
remain plottable.

Output: docs/images/signal_to_noise_distribution_<strategy>.png and .pdf
(~7 × 5 in canvas).
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

import plot_style

REPO_ROOT = repo_root()
IN_DIR = REPO_ROOT / "data" / "expression_evidence"
OUT_DIR = REPO_ROOT / "docs" / "images"
DEFAULT_STRATEGY = "stringent_p5"

# Category buckets — readable label → set of gene_class values from the CSV.
# Pancreas controls are excluded: Cpa1 and Nupr1 turn out to express in
# spinal-cord oligo lineage, so they aren't clean biological negatives —
# the Xenium platform negative-control categories are the only true noise
# floor in this dataset.
CATEGORY_BUCKETS = [
    ("Negative\ncontrols", {"negative_control"}),
    (
        "Targets /\ninducers",
        {
            "Bach2_target",
            "Elf2_target",
            "Foxk2_target",
            "Bhlhe41_target",
            "Nr6a1_target",
            "Sox8_target",
            "Stat3_target",
            "Sox5_inducer",
            "Klk6_inducer",
        },
    ),
    ("Primary\nTFs", {"primary_TF"}),
    ("Markers", {"marker"}),
]

# Box colours: grey for noise; red for hypothesis genes (primary TFs +
# their targets/inducers, same colour because they're conceptually one
# category — TF + downstream); blue for markers as a biological reference.
CATEGORY_COLORS = {
    "Negative\ncontrols": plot_style.PALETTE[7],  # grey
    "Targets /\ninducers": plot_style.PALETTE[5],  # vermillion (same as primary TFs)
    "Primary\nTFs": plot_style.PALETTE[5],  # vermillion
    "Markers": plot_style.PALETTE[4],  # blue
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--strategy",
        default=DEFAULT_STRATEGY,
        help=f"Strategy suffix (default: {DEFAULT_STRATEGY})",
    )
    parser.add_argument(
        "--include-old2",
        action="store_true",
        help="Include Old 2 pseudobulks in the distributions (default: exclude).",
    )
    return parser.parse_args()


def bucket_for(gene_class: str) -> str | None:
    for label, classes in CATEGORY_BUCKETS:
        if gene_class in classes:
            return label
    return None


def main() -> int:
    args = parse_args()
    plot_style.apply()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(IN_DIR / f"expression_per_sample_celltype_{args.strategy}.csv")
    if not args.include_old2:
        df = df[~df["sample_is_old2"]]

    df["category"] = df["gene_class"].map(bucket_for)
    df = df.dropna(subset=["category"])

    # Order categories left → right
    cat_order = [label for label, _ in CATEGORY_BUCKETS]
    df["category"] = pd.Categorical(df["category"], categories=cat_order, ordered=True)

    # Epsilon for log-scale (zero → epsilon so boxplot doesn't lose data)
    eps = 1e-3
    df["plot_value"] = df["mean_counts_per_cell"].fillna(0).clip(lower=eps)

    # Aggregate stats per category for annotation
    stats = df.groupby("category", observed=True)["mean_counts_per_cell"].agg(
        ["median", "min", "max", "count"]
    )

    # --- Figure ---
    fig, ax = plt.subplots(figsize=(6.0, 5.0))

    box_data = [df.loc[df["category"] == cat, "plot_value"].values for cat in cat_order]
    bp = ax.boxplot(
        box_data,
        positions=range(len(cat_order)),
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.6),
        boxprops=dict(linewidth=1),
        whiskerprops=dict(linewidth=1),
        capprops=dict(linewidth=1),
    )
    for patch, cat in zip(bp["boxes"], cat_order):
        patch.set_facecolor(CATEGORY_COLORS[cat])
        patch.set_alpha(0.35)

    # Jittered points
    rng = np.random.default_rng(42)
    for i, cat in enumerate(cat_order):
        vals = df.loc[df["category"] == cat, "plot_value"].values
        jitter = rng.uniform(-0.18, 0.18, len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            s=8,
            alpha=0.45,
            color=CATEGORY_COLORS[cat],
            edgecolor="none",
        )

    # Noise-floor reference line. The all-NC median collapses to 0 because
    # 2 of the 3 NC categories (genomic_control, unassigned_codeword) are
    # essentially all-zero in this dataset. Use control_probe specifically
    # as the meaningful noise floor — it's the category with non-zero
    # background and is the standard Xenium signal-to-noise reference.
    cp_median = float(
        df.loc[df["gene"] == "control_probe", "mean_counts_per_cell"]
        .replace(0, np.nan)
        .dropna()
        .median()
    )
    if not np.isnan(cp_median) and cp_median > 0:
        ax.axhline(
            y=cp_median,
            color="grey",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
            zorder=0,
        )
        # Label OUTSIDE the right edge of the axes so it never overlaps data.
        # clip_on=False allows the text to render in the figure margin.
        ax.text(
            len(cat_order) - 0.4,
            cp_median,
            f"  control_probe\n  noise floor\n  ({cp_median:.2f})",
            fontsize=8,
            color="grey",
            ha="left",
            va="center",
            clip_on=False,
        )

    ax.set_yscale("log")
    ax.set_ylim(eps * 0.5, 100)
    ax.set_xlim(-0.6, len(cat_order) - 0.4)
    ax.set_xticks(range(len(cat_order)))
    ax.set_xticklabels(cat_order, fontsize=11)
    ax.set_ylabel("Mean counts per cell (log scale)", fontsize=11)
    ax.grid(True, axis="y", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)

    # n annotation inside each box, placed at the top of the plot area so it
    # never collides with x-tick labels or with the boxes themselves.
    for i, cat in enumerate(cat_order):
        n = int(stats.loc[cat, "count"])
        ax.text(
            i,
            60,
            f"n = {n}",
            ha="center",
            va="center",
            fontsize=9,
            color="grey",
        )

    fig.suptitle("Expression magnitude across gene categories")

    suffix = "" if not args.include_old2 else "_with_old2"
    plot_style.save(fig, f"signal_to_noise_distribution_{args.strategy}{suffix}")

    # Console summary
    print("\n=== Per-category mean counts/cell summary ===")
    print(stats.round(3).to_string())
    print(
        "\nNote: the all-NC median collapses to 0 because 2 of 3 NC categories "
        "(genomic_control, unassigned_codeword) are essentially all-zero in this "
        "dataset. The meaningful noise floor is the control_probe median."
    )
    print(f"control_probe non-zero median: {cp_median:.3f} transcripts/cell")
    markers_median = float(stats.loc["Markers", "median"])
    primary_tfs_label = "Primary\nTFs"
    primary_tfs_median = float(stats.loc[primary_tfs_label, "median"])
    if cp_median > 0:
        print(
            f"Markers median:    {markers_median:.3f}  ({markers_median / cp_median:.0f}× above control_probe noise)"
        )
        print(
            f"Primary TFs median: {primary_tfs_median:.3f}  ({primary_tfs_median / cp_median:.0f}× above control_probe noise)"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
