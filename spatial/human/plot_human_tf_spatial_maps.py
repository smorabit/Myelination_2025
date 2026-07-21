#!/usr/bin/env python3
"""Spatial maps of TF transcripts across the tissue (segmentation-free).

For one sample, plots the x/y location of every high-quality (QV>=20) transcript
of each gene as a small-multiples grid. Default genes = the 9 primary TFs plus a
myelin positive-control (MBP). This shows the TFs are expressed across the tissue
section, not confined to an artefact, without relying on cell segmentation.

Reads transcripts.parquet directly via human_tf_common (honours HUMAN_TX_DIR).
Default sample is the first non-flagged section in the samplesheet.

Two layout wrinkles this script handles:
  - Some sections carry two physically separated tissue blocks on one capture
    area (stacked along the long axis). Forcing
    both blocks into one shared bounding box under equal aspect leaves large empty
    margins ("gaps") in every panel. `--block a|b` crops to a single block, split at
    the density valley between the two.
  - With `match_aspect` (default in the builder), each panel's axes box is set to
    the data aspect ratio, so equal-aspect data fills the box with no dead margin.

Output: docs/images/human_tf_spatial_maps[_<sample>].png and .pdf
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

import human_tf_common as C
import plot_style

REPO_ROOT = repo_root()
OUT_DIR = REPO_ROOT / "docs" / "images"
DEFAULT_SAMPLE_LABEL = next(
    (s["label"] for s in C.SAMPLES if not s["flagged"]),
    C.SAMPLES[0]["label"] if C.SAMPLES else "",
)
DEFAULT_GENES = C.PRIMARY_TFS + [
    "MBP"
]  # MBP = abundant myelin marker (tissue architecture)
REF_GENES = {"MBP"}  # MBP drawn blue (myelin marker); TFs drawn orange
MAX_POINTS = 40_000  # per-gene render cap (downsample for speed/file size)

REF_COLOR = "#0072B2"
TF_COLOR = "#D55E00"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sample", default=DEFAULT_SAMPLE_LABEL, help="sample_label to plot"
    )
    p.add_argument(
        "--genes",
        default=None,
        help="comma-separated gene list (default: 9 primary TFs + MBP)",
    )
    p.add_argument(
        "--block",
        default="all",
        choices=["all", "a", "b"],
        help="crop to one tissue block split at the density valley along the long "
        "axis: 'a' = lower-coordinate block, 'b' = higher-coordinate block "
        "(default: all = no crop, legacy behaviour).",
    )
    p.add_argument("--ncol", type=int, default=5, help="panels per row")
    return p.parse_args()


def resolve_genes(arg: str | None) -> list[str]:
    if not arg:
        return list(DEFAULT_GENES)
    return [g.strip() for g in arg.split(",") if g.strip()]


def _long_axis_gap(
    coord: np.ndarray,
    nbins: int = 80,
    edge_frac: float = 0.1,
    window: int = 6,
    rise: float = 1.5,
):
    """Return the coordinate of the density valley separating two tissue blocks.

    Two cores on one capture area show as two density modes along the long axis with
    a dip between them. The second core can itself be low-density (a sparse plateau),
    so a global threshold sees "gap continuing to the edge" and a smoothed global
    minimum sits on the tail. Instead this finds the deepest *interior local minimum*
    whose density rises by >= `rise`x within `window` bins on BOTH sides (a real dip
    flanked by tissue, not a monotone tail). Returns that split coordinate, else None.
    """
    lo, hi = float(coord.min()), float(coord.max())
    if hi <= lo:
        return None
    h, edges = np.histogram(coord, bins=nbins, range=(lo, hi))
    i0, i1 = int(nbins * edge_frac), int(nbins * (1 - edge_frac))
    best_j, best_v = None, np.inf
    for j in range(max(i0, 1), min(i1, nbins - 1)):
        if not (h[j] <= h[j - 1] and h[j] <= h[j + 1]):
            continue  # not a local minimum
        left_hi = h[max(0, j - window) : j].max(initial=0)
        right_hi = h[j + 1 : min(nbins, j + 1 + window)].max(initial=0)
        floor = max(h[j], 1)
        if left_hi >= rise * floor and right_hi >= rise * floor and h[j] < best_v:
            best_j, best_v = j, h[j]
    if best_j is None:
        return None
    return 0.5 * (edges[best_j] + edges[best_j + 1])


def crop_to_block(tx, block: str):
    """Crop transcripts to one tissue block; return (subset, label, split_info)."""
    if block == "all":
        return tx, "all", None
    xs = tx["x_location"].to_numpy()
    ys = tx["y_location"].to_numpy()
    long_is_y = (ys.max() - ys.min()) >= (xs.max() - xs.min())
    coord = ys if long_is_y else xs
    axis = "y" if long_is_y else "x"
    thr = _long_axis_gap(coord)
    if thr is None:
        return tx, "all", None  # no clean split -> keep everything
    lower = tx[coord < thr]
    upper = tx[coord >= thr]
    sel = lower if block == "a" else upper
    info = {
        "axis": axis,
        "threshold": float(thr),
        "n_lower": int(len(lower)),
        "n_upper": int(len(upper)),
    }
    return sel, block, info


def build_spatial_maps(
    sample_label: str = DEFAULT_SAMPLE_LABEL,
    genes: list[str] | None = None,
    block: str = "all",
    ncol: int = 5,
    figsize: tuple[float, float] | None = None,
    point_size: float = 0.6,
    alpha: float = 0.3,
    max_points: int = MAX_POINTS,
    match_aspect: bool = True,
    title: str | None = None,
):
    """Build the small-multiples transcript-location figure; return (fig, info).

    `match_aspect=True` sets each axes box to the (cropped) data aspect so equal-
    aspect points fill the panel with no dead margin (the gap fix). The CLI's legacy
    path calls with match_aspect=False to reproduce the original report figure.
    """
    genes = genes or list(DEFAULT_GENES)
    sample = next((s for s in C.SAMPLES if s["label"] == sample_label), None)
    if sample is None:
        raise SystemExit(f"unknown sample {sample_label!r}")

    tx = C.load_transcripts(
        sample["dir"],
        columns=["feature_name", "qv", "is_gene", "x_location", "y_location"],
    )
    tx = tx[(tx["qv"] >= C.QV_MIN) & (tx["is_gene"])]
    tx, block_used, split_info = crop_to_block(tx, block)

    nrow = int(np.ceil(len(genes) / ncol))
    if figsize is None:
        figsize = (2.6 * ncol, 2.6 * nrow)
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize)
    axes = np.atleast_1d(axes).ravel()

    rng = np.random.default_rng(0)
    xmin, xmax = tx["x_location"].min(), tx["x_location"].max()
    ymin, ymax = tx["y_location"].min(), tx["y_location"].max()
    box_aspect = (ymax - ymin) / (xmax - xmin) if xmax > xmin else 1.0

    for ax, gene in zip(axes, genes):
        g = tx[tx["feature_name"] == gene]
        n = len(g)
        if n > max_points:
            g = g.iloc[rng.choice(n, max_points, replace=False)]
        is_ref = gene in REF_GENES
        ax.scatter(
            g["x_location"],
            g["y_location"],
            s=point_size,
            alpha=alpha,
            color=REF_COLOR if is_ref else TF_COLOR,
            edgecolor="none",
        )
        ax.set_title(f"{gene}  (n={n:,})")  # inherit rcParams (uniform across panels)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal")
        if match_aspect:
            ax.set_box_aspect(box_aspect)
        ax.set_xticks([])
        ax.set_yticks([])

    for ax in axes[len(genes) :]:
        ax.axis("off")

    if title is None:
        flag = " (dim, unreliable)" if sample["flagged"] else ""
        blk = f", block {block_used}" if block_used != "all" else ""
        title = f"TF transcript locations across the tissue ({sample_label}{blk}{flag})"
    # Anchor the title a CONSTANT ~0.2 in from the top (and start the axes ~0.5 in
    # down), in inches rather than figure fraction. Otherwise a tall figure (e.g. the
    # single-column variant) turns the same fractional gap into a large empty band
    # between the title and the first panel.
    fh = float(fig.get_size_inches()[1])
    # Leave ~0.85 in above the first panel so the suptitle clears the per-panel gene
    # titles; place the suptitle ~0.3 in from the top. Constant inches, so a tall
    # single-column figure doesn't open a big fractional gap under the title. A single
    # column (ncol==1) needs extra row gap so each gene title clears the panel above.
    hspace = 0.55 if ncol == 1 else 0.2
    fig.subplots_adjust(top=1 - 0.85 / fh, hspace=hspace, wspace=0.05)
    fig.suptitle(title, y=1 - 0.3 / fh, va="top")
    return fig, {"block_used": block_used, "split_info": split_info, "n_tx": len(tx)}


def main() -> int:
    args = parse_args()
    plot_style.apply()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Legacy CLI path reproduces the original report figure when run with defaults
    # (all genes, no crop, square panels on a shared bbox).
    fig, info = build_spatial_maps(
        sample_label=args.sample,
        genes=resolve_genes(args.genes),
        block=args.block,
        ncol=args.ncol,
        match_aspect=(args.block != "all"),
    )
    if info["split_info"]:
        print(f"  block split: {info['split_info']}")

    suffix = (
        ""
        if args.sample == DEFAULT_SAMPLE_LABEL
        else "_" + args.sample.replace("/", "_")
    )
    if args.block != "all":
        suffix += f"_block{args.block}"
    plot_style.save(fig, f"human_tf_spatial_maps{suffix}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
