#!/usr/bin/env python3
"""Headline figure: are the TFs expressed above background in human tissue?

Reads data/xenium_human/expression_evidence/human_tf_evidence.csv and, for each
sample, plots the distribution of high-quality transcript counts per gene, split
into biological categories, against the Xenium negative-control-probe noise floor
(dashed line). This is transcript-level and segmentation-free.

Story the figure tells, left to right: all panel genes sit above the noise floor,
the downstream targets and the 9 primary TFs sit solidly in the expressed range,
and the positive-control brain genes (MBP, PLP1, AQP4, ...) are highest of all.

Output: docs/images/human_tf_signal_to_noise.png and .pdf
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


import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_style

REPO_ROOT = repo_root()
IN_CSV = (
    REPO_ROOT
    / "data"
    / "xenium_human"
    / "expression_evidence"
    / "human_tf_evidence.csv"
)
OUT_DIR = REPO_ROOT / "docs" / "images"

# Left → right by increasing expected expression. (label, colour)
CATEGORIES = [
    ("All other\npanel genes", "#999999"),
    ("Downstream\ntargets", "#E69F00"),
    ("Primary\nTFs", "#D55E00"),
    ("Positive controls\n(brain genes)", "#0072B2"),
]


def bucket_for(row: pd.Series) -> str | None:
    if row["gene_class"] == "negative_control":
        return None
    if row["is_positive_control"]:
        return "Positive controls\n(brain genes)"
    if row["gene_class"] == "primary_tf":
        return "Primary\nTFs"
    if row["gene_class"] == "target":
        return "Downstream\ntargets"
    return "All other\npanel genes"


def main() -> int:
    plot_style.apply()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(IN_CSV)
    df = df[df["gene_class"] != "negative_control"].copy()
    df["category"] = df.apply(bucket_for, axis=1)

    samples = df["sample_label"].drop_duplicates().tolist()
    cat_labels = [c for c, _ in CATEGORIES]
    cat_color = dict(CATEGORIES)
    eps = 0.5  # log floor so zero-count genes stay plottable

    # neg-probe floor per sample (stored on every row as neg_probe_mean)
    floor = df.groupby("sample_label")["neg_probe_mean"].first()

    fig, axes = plt.subplots(
        1, len(samples), figsize=(4.2 * len(samples), 5.2), sharey=True
    )
    if len(samples) == 1:
        axes = [axes]

    rng = np.random.default_rng(42)
    for ax, sample in zip(axes, samples):
        sdf = df[df["sample_label"] == sample]
        box_data = [
            sdf.loc[sdf["category"] == cat, "count"].clip(lower=eps).values
            for cat in cat_labels
        ]
        bp = ax.boxplot(
            box_data,
            positions=range(len(cat_labels)),
            widths=0.6,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color="black", linewidth=1.5),
        )
        for patch, cat in zip(bp["boxes"], cat_labels):
            patch.set_facecolor(cat_color[cat])
            patch.set_alpha(0.4)
        for i, cat in enumerate(cat_labels):
            vals = sdf.loc[sdf["category"] == cat, "count"].clip(lower=eps).values
            jitter = rng.uniform(-0.18, 0.18, len(vals))
            ax.scatter(
                np.full(len(vals), i) + jitter,
                vals,
                s=7,
                alpha=0.4,
                color=cat_color[cat],
                edgecolor="none",
            )

        f = float(floor.get(sample, np.nan))
        if not np.isnan(f) and f > 0:
            ax.axhline(f, color="grey", linestyle="--", linewidth=1, zorder=0)
            ax.text(
                -0.4,
                f,
                "noise floor",
                fontsize=7.5,
                color="grey",
                ha="left",
                va="bottom",
            )

        ax.set_yscale("log")
        ax.set_xlim(-0.6, len(cat_labels) - 0.4)
        ax.set_xticks(range(len(cat_labels)))
        ax.set_xticklabels(cat_labels, fontsize=8.5)
        flag = " (dim, unreliable)" if sdf["flagged"].iloc[0] else ""
        ax.set_title(f"{sample}{flag}", fontsize=10)
        ax.grid(True, axis="y", linestyle=":", alpha=0.4)
        ax.set_axisbelow(True)

    axes[0].set_ylabel("Transcripts per gene (log scale)", fontsize=11)
    fig.suptitle("Transcription factor expression above background")
    fig.subplots_adjust(wspace=0.08)

    plot_style.save(fig, "human_tf_signal_to_noise")
    return 0


if __name__ == "__main__":
    sys.exit(main())
