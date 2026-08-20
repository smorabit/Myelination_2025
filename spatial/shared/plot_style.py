#!/usr/bin/env python3
"""Shared publication figure style for every Python plot in this project.

Single source of truth so all figures share fonts, sizes, colours and export
settings. Import and call `apply()` once at the top of a plotting script, use
`PALETTE` / `CELLTYPE_COLORS` for consistent colours, and `save(fig, name)` to
write a PNG (200 dpi) and a vector PDF with editable text.

    import plot_style
    plot_style.apply()
    ...
    plot_style.save(fig, "human_tf_cross_pcf")

Design choices for publication:
  - Colourblind-safe categorical palette (Okabe-Ito).
  - Top/right spines off, no default grid (add a light y-grid per plot if needed).
  - pdf/ps fonttype 42 so text stays editable in Illustrator/Inkscape.
  - No on-figure caption/footer text: methods notes belong in the paper caption.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

from spatial_paths import repo_root

REPO_ROOT = repo_root()
FIG_DIR = REPO_ROOT / "docs" / "images"

# Okabe-Ito colourblind-safe categorical palette (orange, sky blue, bluish green,
# yellow, blue, vermillion, reddish purple, grey). Black is kept out of the cycle
# so lines never collide with axis/reference black.
PALETTE = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#999999",
]

# Diverging colormap for enrichment / log2FC heatmaps (blue = low, red = high).
DIVERGING_CMAP = "RdBu_r"

# Consistent colours for the oligodendrocyte lineage across composition, lineage
# and spatial figures. Blue -> green -> orange follows OPC -> intermediate -> mature.
CELLTYPE_COLORS = {
    "OPC": "#0072B2",
    "Intermediate": "#009E73",
    "Intermediate Oligo": "#009E73",
    "Mature": "#D55E00",
    "Mature Oligo": "#D55E00",
    "Oligo": "#D55E00",
}

# Fixed group colours for young vs old comparisons.
GROUP_COLORS = {"Young": "#56B4E9", "Old": "#D55E00"}

# TF-group colours for the DE volcano / log2FC heatmap / expression heatmap+dotplot.
# The 9 primary TFs use Paul Tol's "muted" qualitative palette, which is designed to
# stay mutually distinguishable under deuteranopia, protanopia and tritanopia (Tol,
# "Colour Schemes", SRON/EPS-TN, https://personal.sron.nl/~pault/). The old Set1-style
# map had a red/green (Bach2/Foxk2) pair that fails for red-green colourblindness. The
# two headline factors Sox8 (wine) and Klk6 (teal) are given well-separated hues. The
# neutral categories are greys/black separated by lightness (they never co-occur in one
# figure, so pancreas/stress and negative ctrl can share the mid-grey band).
TF_GROUP_COLORS = {
    "Bach2": "#CC6677",  # rose
    "Elf2": "#332288",  # indigo
    "Foxk2": "#DDCC77",  # sand
    "Bhlhe41": "#AA4499",  # purple
    "Nr6a1": "#88CCEE",  # cyan
    "Sox8": "#882255",  # wine (headline)
    "Stat3": "#117733",  # green
    "Sox5": "#999933",  # olive
    "Klk6": "#44AA99",  # teal (headline)
    "marker": "#000000",  # black
    "pancreas/stress": "#888888",  # mid grey
    "negative ctrl": "#888888",  # mid grey (never shares a figure with pancreas/stress)
    "other": "#CCCCCC",  # light grey
}


def apply() -> None:
    """Install the house style into matplotlib rcParams (idempotent)."""
    mpl.rcParams.update(
        {
            # fonts
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.titleweight": "normal",
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 13,
            # frame: drop the top/right spines, no default grid
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "axes.linewidth": 0.8,
            # legend: no box, so it doesn't fence off part of the plot
            "legend.frameon": False,
            # colour cycle
            "axes.prop_cycle": mpl.cycler(color=PALETTE),
            "image.cmap": DIVERGING_CMAP,
            # export: crisp raster, editable vector text
            "figure.dpi": 150,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save(fig, name: str, fig_dir: Path | None = None) -> None:
    """Save `fig` as <fig_dir>/<name>.png (200 dpi) and .pdf (editable text)."""
    out_dir = fig_dir or FIG_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / name
    fig.savefig(f"{out}.png")
    fig.savefig(f"{out}.pdf")
    plt.close(fig)
    try:
        rel = out_dir.relative_to(REPO_ROOT)
    except ValueError:
        rel = out_dir
    print(f"wrote {rel}/{name}.png + .pdf")
