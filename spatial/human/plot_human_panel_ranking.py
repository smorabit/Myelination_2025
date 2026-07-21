#!/usr/bin/env python3
"""Full-panel expression ranking for one sample.

Reads data/xenium_human/expression_evidence/human_tf_evidence.csv and draws all
316 panel genes ranked by high-quality transcript count (log scale), so a reader
can see which genes are expressed in this human tissue and where the 9 primary TFs
and the positive-control brain genes fall relative to everything else.

Genes are drawn as a thin grey ranked curve; primary TFs (red) and positive
controls (blue) are highlighted and labelled. The negative-control-probe noise
floor is a dashed line.

Default sample is the first non-flagged section in the samplesheet.

Output: docs/images/human_panel_ranking.png and .pdf
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

import cohort

REPO_ROOT = repo_root()
IN_CSV = (
    REPO_ROOT
    / "data"
    / "xenium_human"
    / "expression_evidence"
    / "human_tf_evidence.csv"
)
OUT_DIR = REPO_ROOT / "docs" / "images"
_HUMAN = cohort.human_samples()
DEFAULT_SAMPLE = next(
    (s["label"] for s in _HUMAN if not s["flagged"]),
    _HUMAN[0]["label"] if _HUMAN else "",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--sample", default=DEFAULT_SAMPLE, help="sample_label to plot")
    p.add_argument(
        "--proliferation",
        action="store_true",
        help="Variant: also mark proliferation/cell-cycle genes (expected-low reference).",
    )
    return p.parse_args()


def main() -> int:
    import human_tf_common as C
    import plot_style

    args = parse_args()
    plot_style.apply()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(IN_CSV)
    df = df[
        (df["sample_label"] == args.sample) & (df["gene_class"] != "negative_control")
    ]
    if df.empty:
        print(f"no rows for sample {args.sample!r}", file=sys.stderr)
        return 1
    df = df.sort_values("count", ascending=False).reset_index(drop=True)
    df["x"] = np.arange(1, len(df) + 1)
    floor = float(df["neg_probe_mean"].iloc[0])

    fig, ax = plt.subplots(figsize=(11, 5.5))
    fig.subplots_adjust(top=0.88, bottom=0.12, left=0.07, right=0.98)

    # all genes as a thin grey stem/curve
    ax.plot(
        df["x"], df["count"].clip(lower=0.5), color="#999999", linewidth=1.2, zorder=1
    )

    tf = df[df["gene_class"] == "primary_tf"]
    pc = df[df["is_positive_control"]]
    ax.scatter(
        tf["x"], tf["count"], color="#D55E00", s=45, zorder=3, label="Primary TFs"
    )
    ax.scatter(
        pc["x"],
        pc["count"],
        color="#0072B2",
        s=45,
        zorder=3,
        label="Positive controls (brain genes)",
    )

    pr = (
        df[df["gene"].isin(C.PROLIFERATION_MARKERS)]
        if args.proliferation
        else df.iloc[0:0]
    )
    if args.proliferation:
        ax.scatter(
            pr["x"],
            pr["count"],
            color="#009E73",
            s=45,
            zorder=3,
            label="Proliferation genes (expected low)",
        )

    # label TFs, positive controls, and (in the variant) proliferation genes
    label_color = {
        **{g: "#0072B2" for g in pc["gene"]},
        **{g: "#D55E00" for g in tf["gene"]},
    }
    if args.proliferation:
        label_color.update({g: "#009E73" for g in pr["gene"]})
    # stagger label heights by rank so crowded tails (TFs among proliferation genes) don't overlap
    labelled = (
        pd.concat([tf, pc, pr])
        .drop_duplicates("gene")
        .sort_values("x")
        .reset_index(drop=True)
    )
    for i, r in labelled.iterrows():
        ax.annotate(
            r["gene"],
            (r["x"], r["count"]),
            xytext=(0, 8 + 12 * (i % 3)),
            textcoords="offset points",
            fontsize=7.5,
            color=label_color.get(r["gene"], "#333333"),
            ha="center",
            rotation=45,
        )

    if floor > 0:
        ax.axhline(floor, color="grey", linestyle="--", linewidth=1, zorder=0)
        ax.text(
            len(df),
            floor,
            " noise floor",
            fontsize=8,
            color="grey",
            ha="right",
            va="bottom",
        )

    ax.set_yscale("log")
    ax.set_xlabel("Panel genes ranked by expression (highest → lowest)", fontsize=11)
    ax.set_ylabel("Transcripts per gene (log scale)", fontsize=11)
    ax.set_xlim(0, len(df) + 1)
    ax.grid(True, axis="y", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", fontsize=9)

    flag = " (dim, unreliable)" if bool(df["flagged"].iloc[0]) else ""
    fig.suptitle(
        f"Where TFs rank among all {len(df)} panel genes ({args.sample}{flag})"
    )
    name = (
        "human_panel_ranking_proliferation"
        if args.proliferation
        else "human_panel_ranking"
    )
    plot_style.save(fig, name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
