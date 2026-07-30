#!/usr/bin/env python3
"""Nature Comms revision — standalone manuscript figure-panel generator.

Generates the individual component figures (A-D) for one composite panel, each as an
editable vector PDF (+ PNG preview) in the project house style. They are handed to an
illustrator to assemble; this script exists so any figure can be regenerated / resized
by editing the CONFIG block below and re-running, with no other code changes.

Figures:
  A  DAPI morphology, zoomed on the lesion (per sample). True Xenium morphology image.
  B  Zoomed single-cell view: cell outlines coloured by oligo lineage + grey celltype-
     marker transcripts + 3 colour-coded TF transcripts per row. TFs ranked by
     expression; rows = TF ranks 1-3 / 4-6 / 7-9 (top9) or 1-3 / 4-6 (top6). Rendered
     per sample (single column) and as a composed Young|Old pair.
  C  Aggregated TF-expression dot plot, rescaled (delegates to
     plot_expression_aggregated_dotplot.build_dotplot).
  D  Sox8 differential-expression barplot: mean counts/cell, Young vs Old, across the
     3 oligo-lineage cell types (SEM bars, sample dots, padj stars; outlier-removed set).
  suppl  Human TF spatial maps (former Panel D; MBP + TFs), cropped to one tissue block
     (delegates to plot_human_tf_spatial_maps.build_spatial_maps).

Environment / data sources (read-only):
  - Mouse raw Xenium bundle (morphology, boundaries, transcripts):
      XENIUM_RAW_DIR  (default = the raw Xenium bundle path below)
  - Human transcripts (Fig D):
      HUMAN_TX_DIR    (e.g. /path/to/human/transcripts)
  - Processed per-sample AnnData + refined celltypes + lesion masks: repo data/ (via
    spatial_common).

Run (under the xenium-processing env):
    micromamba activate xenium-processing
    export XENIUM_RAW_DIR=/path/to/mouse/xenium/bundle
    export HUMAN_TX_DIR=/path/to/human/transcripts
    python manuscript_panel.py all
    python manuscript_panel.py B --tf-set top6 --samples <MOUSE_SAMPLE_ID>

Outputs: Manuscript/manuscript_figures/*.{pdf,png}
"""

from __future__ import annotations

# --- locate deposited module dirs (works from any depth) ------------------
import sys as _sys
from pathlib import Path as _Path

for _p in _Path(__file__).resolve().parents:
    if (_p / "shared" / "spatial_common.py").exists():
        for _sub in (
            "shared",
            "mouse/expression",
            "human",
        ):
            if (_p / _sub).is_dir() and str(_p / _sub) not in _sys.path:
                _sys.path.insert(0, str(_p / _sub))
        break
from spatial_paths import repo_root


import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import zarr
from matplotlib.collections import PolyCollection

import cohort  # noqa: E402
import plot_style  # noqa: E402
import spatial_common as sc  # noqa: E402
from annotate_celltypes import (  # noqa: E402
    CORE_MARKERS,
    INTER_MARKERS,
    MATURE_EXTRAS,
    OPC_EXTRAS,
)
import plot_expression_aggregated_dotplot as dotplot  # noqa: E402
import plot_human_tf_spatial_maps as humanmaps  # noqa: E402

# ============================ CONFIG (edit here) ============================
REPO_ROOT = repo_root()
# Top-level output root is env-overridable (MANUSCRIPT_ROOT, default `Manuscript`);
# figures land in Manuscript/manuscript_figures/.
OUT_DIR = (
    REPO_ROOT / os.environ.get("MANUSCRIPT_ROOT", "Manuscript") / "manuscript_figures"
)

# Raw mouse Xenium bundle (per-sample output-XETG<serial>__<sid>__* subdirs).
XENIUM_RAW_DIR = Path(
    os.environ.get(
        "XENIUM_RAW_DIR",
        "/path/to/mouse/xenium/bundle",
    )
)
PIXEL_SIZE_FALLBACK = 0.2125  # µm/px; per-sample value read from experiment.xenium

# sample_id -> short filename-safe label (e.g. "Young1"), from the samplesheet
# (shared/cohort.py). No real sample IDs are stored in this repository.
SAMPLE_LABELS = {
    sid: lbl.replace(" ", "") for sid, lbl in cohort.sample_labels().items()
}
# Render A/B across all cohort-9 samples by default (user picks the winner later).
AB_YOUNG = [s for s in sc.COHORT_9 if sc.age_for(s) == "Young"]
AB_OLD = [s for s in sc.COHORT_9 if sc.age_for(s) == "Old"]
# Composed Young|Old pairs for Fig B. Default = every Young x Old combination so the
# best pairing can be chosen later; trim this list to reduce output.
B_PAIRS = [(y, o) for y in AB_YOUNG for o in AB_OLD]

# Per-figure canvas sizes in inches, sized for panel placement. Edit to resize.
# --- Multipanel layout ------------------------------------------------------
# Composite target: Nature Comms double-column, 180 mm = 7.087 in wide. Each panel is
# rendered at the PHYSICAL FOOTPRINT it occupies in the composite (from the layout
# mock-up), and every panel shares one font size (FONT_PT). Place each at 100 % in
# Illustrator (no rescaling) and the text matches across panels. Edit these to resize.
FONT_PT = 6.0  # shared body font (pt) at final print size — Nature-typical
FONT_TITLE_PT = 7.0  # panel titles / suptitles, one step up
FIG_SIZES = {
    # A = one morphology crop, cropped square so all samples share output height; two
    # (Young, Old) sit side by side. Slightly bigger than the first draft.
    "A": (1.2, 1.2),
    # B = one sample column of square ROI tiles; the Young|Old pair is 2x this width.
    # Height is the 3-row height (top6 auto-scales to 2/3).
    "B_col": (0.98, 3.1),
    "C": (2.8, 6.03),  # centre dot plot; height 10% shorter than the 6.7 draft
    # D = Sox8 DE barplot: small block under the dot plot (C). Compact footprint.
    "D": (1.4, 1.26),
    # Supplementary (former Panel D): human TF spatial maps — row + tall composite variants.
    "SUPPL_HUMAN": (7.0, 1.9),
    "SUPPL_HUMAN_vertical": (1.9, 7.2),
}
DOTPLOT_DOT_RANGE = 48  # dot-area scale for the ~2.8 in-wide manuscript dot plot
DOTPLOT_CBAR_FRAC = 0.5  # colorbar occupies this fraction of the panel height (centred)

# Fig B options.
ROI_SIZE_UM = 240.0  # side length of the zoom window
ROI_MIN_OLIGO = 8  # min annotated oligo-lineage cells for a window to qualify
ROI_TARGET_RNORM = (
    0.55  # target lesion depth (0=centroid, 1=boundary) for comparability
)
ROI_RNORM_SIGMA = 0.22  # tolerance around the target depth
ROI_MANUAL: dict[str, tuple[float, float]] = {}  # sid -> (x0, y0) lower-left override
LESION_PAD_UM = 60.0  # padding around lesion bbox for ROI search

# Biological groupings of the 9 TFs into rows (replaces expression-rank ordering, which
# crowded all the abundant TFs into the top row). Each grouping is an ordered list of
# (row_label, [3 TFs]); top9 shows all 3 rows, top6 the first 2. Grounded in the docs
# (de_results.md, spatial_tf_analysis.md, TF_expression.md). Rendered for each grouping
# so the framing can be chosen.
B_GROUPINGS = {
    # A — molecular family / oligo-lineage. Best within-row colour separation; keeps the
    # two oligo-specific, human-conserved factors (Sox8, Klk6) together.
    # Short row labels (rotated on the left edge; must fit one row's height).
    "A": [
        ("Oligo lineage", ["Sox8", "Sox5", "Klk6"]),
        ("bHLH/bZIP/NR", ["Bhlhe41", "Bach2", "Nr6a1"]),
        ("STAT/ETS/FOX", ["Stat3", "Elf2", "Foxk2"]),
    ],
    # C — direction of change with age (the study's narrative; Sox8 leads "down in Old").
    "C": [
        ("Down in Old", ["Sox8", "Sox5", "Bach2"]),
        ("Up in Old", ["Bhlhe41", "Stat3", "Klk6"]),
        ("Flat / low", ["Elf2", "Foxk2", "Nr6a1"]),
    ],
}
B_GROUPINGS_TO_RENDER = ["A", "C"]


# TF dots use the canonical per-TF group colours (Paul Tol muted) shared with the
# dot plot / DE figures, so a TF is the same colour everywhere in the manuscript.
def tf_color(tf: str) -> str:
    return plot_style.TF_GROUP_COLORS.get(tf, "#333333")


TF_DOT_SIZE = 3.5  # scatter point² for TF transcripts (small so cells show through)
MARKER_DOT_SIZE = 1.2  # scatter point² for celltype-marker transcripts (stage-coloured)
LEGEND_COL_IN = 0.62  # width (in) of the per-row TF legend column, right of each row
LEGEND_FACECOLOR = "white"
LEGEND_ALPHA = 0.75  # semi-opaque so legends stay readable over dense transcripts
LINEAGE_COLOR = {
    "OPC": plot_style.CELLTYPE_COLORS["OPC"],
    "Intermediate_Oligo": plot_style.CELLTYPE_COLORS["Intermediate"],
    "Mature_Oligo": plot_style.CELLTYPE_COLORS["Mature"],
}
OTHER_CELL_COLOR = "#dddddd"
# Celltype-marker transcript dots, coloured by the lineage stage they mark. Muted /
# neutral so they read as evidence under the vivid TF dots; each stage keeps the hue
# family of its outline (dot = desaturated version of the ring).
MARKER_STAGE_COLORS = {
    "pan": "#9a9a9a",  # grey — pan-lineage (Olig2, Sox10)
    "OPC": "#6d8fb5",  # muted blue
    "Intermediate": "#6fa06f",  # muted green
    "Mature": "#9c6b4a",  # brown
}
MARKER_STAGE = {
    **{g: "pan" for g in CORE_MARKERS},
    **{g: "OPC" for g in OPC_EXTRAS},
    **{g: "Intermediate" for g in INTER_MARKERS},
    **{g: "Mature" for g in MATURE_EXTRAS},
}
# Mbp excluded by default: its mRNA diffuses 5-13x further than Mog/Opalin and is
# dropped from stringent_p5 mature typing, so showing it would overstate mature ID.
MARKER_GENES = [g for g in MARKER_STAGE if g != "Mbp"]

# Fig D (Sox8 DE barplot). Mean counts/cell, Young vs Old, across the 3 oligo-lineage
# cell types. Bars use the canonical age colours (plot_style.GROUP_COLORS); the DE-excluded
# outlier sample is dropped so the bar means and the padj come from the same sample set.
BARPLOT_GENE = "Sox8"
KEY_CELLTYPES = ["OPC", "Intermediate_Oligo", "Mature_Oligo"]
CELLTYPE_LABELS = {
    "OPC": "OPC",
    "Intermediate_Oligo": "Intermediate",
    "Mature_Oligo": "Mature",
}
# Short x-tick labels for the small Panel-D footprint (full names collide at ~1.4 in wide).
CELLTYPE_LABELS_SHORT = {
    "OPC": "OPC",
    "Intermediate_Oligo": "Int",
    "Mature_Oligo": "Mat",
}
EXPRESSION_CSV = (
    REPO_ROOT
    / "data"
    / "expression_evidence"
    / "expression_per_sample_celltype_stringent_p5.csv"
)
DE_CSV = (
    REPO_ROOT / "data" / "de_results" / "stringent_p5" / "_combined_de_outlier_rm.csv"
)

# Supplementary (former Panel D): human TF spatial maps. Default section = first
# non-flagged human sample from the samplesheet (shared/cohort.py); override via CLI.
_HUMAN = cohort.human_samples()
HUMAN_SAMPLE = next(
    (s["label"] for s in _HUMAN if not s["flagged"]),
    _HUMAN[0]["label"] if _HUMAN else "",
)
HUMAN_BLOCK = "a"  # 'a' = dense block (y<valley); flip to 'b' if that is bottom_144
HUMAN_GENES = ["MBP", "STAT3", "SOX8", "SOX5", "KLK6", "BHLHE41"]
HUMAN_POINT_SIZE = 0.8  # small: the vertical composite panels are only ~1.9 in wide
HUMAN_POINT_ALPHA = 0.4
# ===========================================================================


# ------------------------- shared helpers -------------------------
def apply_manuscript_style() -> None:
    """House style + a single manuscript font size across every panel.

    Overrides the larger plot_style defaults so all panels share FONT_PT. Inline
    per-call font sizes are avoided elsewhere so text stays uniform; render each panel
    at its FIG_SIZES footprint and place at 100 % for matching fonts in the composite.
    """
    plot_style.apply()
    mpl.rcParams.update(
        {
            "font.size": FONT_PT,
            "axes.titlesize": FONT_TITLE_PT,
            "axes.labelsize": FONT_PT,
            "xtick.labelsize": FONT_PT,
            "ytick.labelsize": FONT_PT,
            "legend.fontsize": FONT_PT,
            "legend.title_fontsize": FONT_PT,
            "figure.titlesize": FONT_TITLE_PT,
            "axes.linewidth": 0.5,
        }
    )


def sample_raw_dir(sid: str) -> Path:
    hits = sorted(XENIUM_RAW_DIR.glob(f"output-XETG*__{sid}__*"))
    if not hits:
        raise FileNotFoundError(
            f"no raw bundle for {sid} under {XENIUM_RAW_DIR} (set XENIUM_RAW_DIR)"
        )
    return hits[0]


def read_pixel_size(raw_dir: Path) -> float:
    try:
        d = json.loads((raw_dir / "experiment.xenium").read_text())
        return float(d.get("pixel_size", PIXEL_SIZE_FALLBACK))
    except Exception:
        return PIXEL_SIZE_FALLBACK


def lesion_bbox(adata, pad: float) -> tuple[float, float, float, float]:
    xy = adata.obs.loc[adata.obs["in_lesion"], ["x_centroid", "y_centroid"]].to_numpy(
        dtype=float
    )
    if len(xy) == 0:
        raise ValueError("no in-lesion cells")
    x0, y0 = xy.min(0) - pad
    x1, y1 = xy.max(0) + pad
    return float(x0), float(y0), float(x1), float(y1)


def add_scalebar(ax, length_um: float, label: str, color: str = "black") -> None:
    """Draw a scale bar bottom-left. Position is in AXES FRACTION and the bar length is
    the correct data-fraction, so it never spills outside the panel (a fixed data-length
    bar overflows on small crops / letterboxed images)."""
    dx = abs(ax.get_xlim()[1] - ax.get_xlim()[0])
    frac = min(length_um / dx, 0.85) if dx else 0.3  # bar length as axes fraction
    x_start, y = 0.06, 0.07  # near the lower-left of the axes box
    ax.plot(
        [x_start, x_start + frac],
        [y, y],
        transform=ax.transAxes,
        color=color,
        lw=1.3,
        solid_capstyle="butt",
        clip_on=False,
    )
    ax.text(
        x_start + frac / 2,
        y + 0.025,
        label,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=FONT_PT,
        color=color,
    )


def read_dapi_crop(raw_dir: Path, bbox_um, px: float, target_px: int = 2200):
    """Return (crop_array, extent) for the lesion bbox from the DAPI OME-TIFF.

    Picks a pyramid level so the crop's larger side is ~target_px (keeps memory and
    file size sane); reads only that window. extent is in microns for imshow with
    origin='upper' (y increases downward, matching Xenium image/centroid coords).
    """
    x0, y0, x1, y1 = bbox_um
    store = tifffile.imread(
        raw_dir / "morphology_focus" / "ch0000_dapi.ome.tif", aszarr=True
    )
    z = zarr.open(store, mode="r")
    levels = sorted((k for k in z.keys()), key=int) if hasattr(z, "keys") else ["0"]
    # full-res pixel window
    c0, c1 = int(x0 / px), int(x1 / px)
    r0, r1 = int(y0 / px), int(y1 / px)
    maxdim = max(c1 - c0, r1 - r0)
    lvl = 0
    while lvl + 1 < len(levels) and maxdim / (2 ** (lvl + 1)) > target_px:
        lvl += 1
    arr = z[levels[lvl]] if hasattr(z, "keys") else z
    f = 2**lvl
    crop = np.asarray(arr[max(0, r0 // f) : r1 // f, max(0, c0 // f) : c1 // f])
    extent = [c0 * px, c1 * px, r1 * px, r0 * px]  # left, right, bottom, top
    return crop, extent


def load_hq_transcripts(raw_dir: Path) -> pd.DataFrame:
    tx = pd.read_parquet(
        raw_dir / "transcripts.parquet",
        columns=[
            "feature_name",
            "x_location",
            "y_location",
            "qv",
            "is_gene",
            "cell_id",
        ],
    )
    return tx[(tx["qv"] >= 20) & (tx["is_gene"])]


def rank_tfs(strategy: str = sc.DEFAULT_STRATEGY) -> list[str]:
    """Primary TFs ordered by mean transcripts/cell in oligo lineage (Old-2 excluded)."""
    df = pd.read_csv(
        REPO_ROOT
        / "data"
        / "expression_evidence"
        / f"expression_per_sample_celltype_{strategy}.csv"
    )
    sub = df[
        df["gene"].isin(sc.PRIMARY_TFS)
        & df["cell_type"].isin(sc.OLIGO_TYPES)
        & (~df["sample_is_old2"])
    ]
    order = (
        sub.groupby("gene")["mean_counts_per_cell"].mean().sort_values(ascending=False)
    )
    return order.index.tolist()


# ------------------------- Figure A -------------------------
def _square_bbox(bbox):
    """Expand a bbox to a centred square so every A crop has the same (1:1) aspect and
    therefore the same output height when placed at a fixed square figure size."""
    x0, y0, x1, y1 = bbox
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    s = max(x1 - x0, y1 - y0) / 2
    return cx - s, cy - s, cx + s, cy + s


def figure_A(sid: str, show_lesion: bool = False):
    adata = sc.load_sample_cells(sid, oligo_only=False)
    bbox = _square_bbox(lesion_bbox(adata, pad=LESION_PAD_UM))
    raw = sample_raw_dir(sid)
    px = read_pixel_size(raw)
    crop, extent = read_dapi_crop(raw, bbox, px)

    fig, ax = plt.subplots(figsize=FIG_SIZES["A"])
    # DAPI is low-signal fluorescence: most pixels are background, so a plain
    # percentile stretch leaves tissue midtones dark. Normalise to [2, 99.7] pct
    # then apply gamma 0.7 to lift midtones so morphology reads clearly.
    if crop.size:
        lo, hi = np.percentile(crop, [2, 99.7])
        norm = np.clip((crop.astype(float) - lo) / max(hi - lo, 1.0), 0, 1) ** 0.7
    else:
        norm = crop.astype(float)
    ax.imshow(norm, cmap="gray", vmin=0, vmax=1, origin="upper", extent=extent)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])  # inverted (y down)
    ax.set_aspect("equal")
    ax.axis("off")
    if show_lesion:
        les = adata.obs.loc[
            adata.obs["in_lesion"], ["x_centroid", "y_centroid"]
        ].to_numpy()
        ax.scatter(les[:, 0], les[:, 1], s=1, c="#D55E00", alpha=0.3, linewidths=0)
    add_scalebar(ax, 200.0, "200 µm", color="white")
    return fig


# ------------------------- Figure B -------------------------
def _radial_norm(px: float, py: float, centroid, hull) -> float:
    """Normalised lesion depth of a point: 0 at the hull centroid, 1 at its boundary.

    Casts a ray from the centroid through the point to the hull boundary; returns
    |p - centroid| / |boundary - centroid| along that ray (clipped to [0, 1.3]). This
    lets Fig B sample a comparable core-to-edge position across samples of different
    lesion size/shape, so a Young-vs-Old difference is not confounded by depth.
    """
    from shapely.geometry import LineString

    cx, cy = centroid
    d = float(np.hypot(px - cx, py - cy))
    if d < 1e-6 or hull is None or hull.is_empty:
        return 0.0
    ux, uy = (px - cx) / d, (py - cy) / d
    ray = LineString([(cx, cy), (cx + ux * 1e5, cy + uy * 1e5)])
    inter = ray.intersection(hull.boundary)
    if inter.is_empty:
        return 1.3
    coords = np.asarray(
        [list(g.coords)[0] for g in getattr(inter, "geoms", [inter])], dtype=float
    ).reshape(-1, 2)
    r_edge = float(np.hypot(coords[:, 0] - cx, coords[:, 1] - cy).max())
    return min(d / r_edge, 1.3) if r_edge > 1e-6 else 1.3


def find_roi(
    tf_xy,
    tf_feat,
    oligo_xy,
    hull,
    centroid,
    bbox,
    size_um,
    min_oligo,
    target_rnorm,
    rnorm_sigma,
):
    """Choose the zoom window by three criteria (the reason Fig B is comparable):

      1. TF signal      — n TF transcripts x (n distinct TFs + 1) in the window,
      2. celltype cover — n annotated oligo-lineage cells in the window (gate + factor),
      3. lesion depth   — window centre near `target_rnorm` core-to-edge position,
                          weighted by a Gaussian of width `rnorm_sigma`.

    Returns (x0, y0) lower-left in microns. Falls back to the best signal window if
    none clear the oligo-cell gate.
    """
    x0, y0, x1, y1 = bbox
    step = size_um / 3.0
    xs = np.arange(x0, max(x0 + step, x1 - size_um) + step, step)
    ys = np.arange(y0, max(y0 + step, y1 - size_um) + step, step)
    ox, oy = oligo_xy[:, 0], oligo_xy[:, 1]
    tx_, ty_ = tf_xy[:, 0], tf_xy[:, 1]
    best, best_score = None, -1.0
    best_nogate, best_score_nogate = None, -1.0
    for wx in xs:
        for wy in ys:
            omask = (ox >= wx) & (ox < wx + size_um) & (oy >= wy) & (oy < wy + size_um)
            n_oligo = int(omask.sum())
            tmask = (
                (tx_ >= wx) & (tx_ < wx + size_um) & (ty_ >= wy) & (ty_ < wy + size_um)
            )
            ntf = int(tmask.sum())
            ndist = len(np.unique(tf_feat[tmask])) if ntf else 0
            rnorm = _radial_norm(wx + size_um / 2, wy + size_um / 2, centroid, hull)
            depth_w = float(np.exp(-(((rnorm - target_rnorm) / rnorm_sigma) ** 2)))
            signal = ntf * (ndist + 1) * (1 + n_oligo)
            score = signal * depth_w
            if score > best_score_nogate:
                best_score_nogate, best_nogate = score, (float(wx), float(wy))
            if n_oligo >= min_oligo and score > best_score:
                best_score, best = score, (float(wx), float(wy))
    return best if best is not None else best_nogate


def _cell_polys_in_roi(cb: pd.DataFrame, celltype: pd.Series, cell_ids):
    """Split ROI cells into (background non-oligo polys) and (oligo polys, colours).

    Oligo-lineage cells carry the celltype-annotation signal the panel is about, so
    they are drawn thick/coloured on top; non-oligo cells are faint background.
    """
    sub = cb[cb["cell_id"].isin(set(cell_ids))]
    bg_polys, oligo_polys, oligo_colors = [], [], []
    for cid, grp in sub.groupby("cell_id"):
        verts = grp[["vertex_x", "vertex_y"]].to_numpy()
        lineage = str(celltype.get(cid, "other"))
        if lineage in LINEAGE_COLOR:
            oligo_polys.append(verts)
            oligo_colors.append(LINEAGE_COLOR[lineage])
        else:
            bg_polys.append(verts)
    return bg_polys, oligo_polys, oligo_colors


def _draw_roi_panel(
    ax, roi, cb, celltype, adata, txb, row_genes, size_um, row_label=None
):
    """One zoomed panel: lineage cell outlines + stage-coloured marker dots + TF dots."""
    rx, ry = roi
    obs = adata.obs
    in_roi = (
        (obs["x_centroid"] >= rx)
        & (obs["x_centroid"] < rx + size_um)
        & (obs["y_centroid"] >= ry)
        & (obs["y_centroid"] < ry + size_um)
    )
    bg_polys, oligo_polys, oligo_colors = _cell_polys_in_roi(
        cb, celltype, obs.index[in_roi]
    )
    if bg_polys:  # non-oligo cells: faint background context
        ax.add_collection(
            PolyCollection(
                bg_polys,
                facecolors="none",
                edgecolors="#e2e2e2",
                linewidths=0.2,
                zorder=1,
            )
        )
    if oligo_polys:  # oligo lineage: coloured outline (the annotation cue)
        ax.add_collection(
            PolyCollection(
                oligo_polys,
                facecolors="none",
                edgecolors=oligo_colors,
                linewidths=0.6,
                zorder=2,
            )
        )
    win = txb[
        (txb["x_location"] >= rx)
        & (txb["x_location"] < rx + size_um)
        & (txb["y_location"] >= ry)
        & (txb["y_location"] < ry + size_um)
    ]
    mk = win[win["feature_name"].isin(MARKER_GENES)]
    mk_colors = [MARKER_STAGE_COLORS[MARKER_STAGE[g]] for g in mk["feature_name"]]
    ax.scatter(
        mk["x_location"],
        mk["y_location"],
        s=MARKER_DOT_SIZE,
        c=mk_colors,
        alpha=0.6,
        linewidths=0,
        zorder=1.5,
    )
    for (
        gene
    ) in row_genes:  # TF transcripts (canonical group colour), no in-panel legend
        g = win[win["feature_name"] == gene]
        ax.scatter(
            g["x_location"],
            g["y_location"],
            s=TF_DOT_SIZE,
            c=tf_color(gene),
            edgecolors="white",
            linewidths=0.2,
            zorder=3,
        )
    ax.set_xlim(rx, rx + size_um)
    ax.set_ylim(ry + size_um, ry)  # y down to match image convention
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    if row_label:  # biological group label as a rotated row descriptor (left edge)
        ax.set_ylabel(row_label, fontsize=FONT_PT, rotation=90, labelpad=3)


def _row_tf_legend(ax, row_genes):
    """Draw one row's 3-TF colour key into a dedicated (legend-only) axes on the right."""
    ax.axis("off")
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            ls="",
            mfc=tf_color(g),
            mec="white",
            mew=0.3,
            ms=5,
            label=g,
        )
        for g in row_genes
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        frameon=False,
        handletextpad=0.3,
        labelspacing=0.4,
        borderpad=0.0,
    )


_PREP_CACHE: dict[tuple, dict] = {}  # (sid, genes) -> prepped dict; reused across pairs


def _prep_sample_B(sid: str, genes: list[str]):
    """Load everything Fig B needs for one sample and choose its ROI (cached).

    Cached per (sid, genes) so the many Young x Old pair figures reuse the per-sample
    load + ROI instead of re-reading AnnData / transcripts / boundaries each time.
    """
    key = (sid, tuple(sorted(genes)))  # ROI depends on the gene set, not its order
    if key in _PREP_CACHE:
        return _PREP_CACHE[key]
    adata = sc.load_sample_cells(sid, oligo_only=False)
    celltype = pd.Series(
        adata.obs["celltype_refined"].to_numpy(), index=adata.obs_names.astype(str)
    )
    raw = sample_raw_dir(sid)
    tx = load_hq_transcripts(raw)
    bx = lesion_bbox(adata, pad=LESION_PAD_UM)
    txb = tx[
        (tx["x_location"] >= bx[0])
        & (tx["x_location"] <= bx[2])
        & (tx["y_location"] >= bx[1])
        & (tx["y_location"] <= bx[3])
    ]
    tf_tx = txb[txb["feature_name"].isin(genes)]
    # lesion geometry for the radial-depth term + oligo-cell coordinates for coverage
    les_xy = adata.obs.loc[
        adata.obs["in_lesion"], ["x_centroid", "y_centroid"]
    ].to_numpy(dtype=float)
    hull, _ = sc.alpha_shape(les_xy)
    centroid = (float(hull.centroid.x), float(hull.centroid.y))
    oligo_mask = adata.obs["celltype_refined"].isin(sc.OLIGO_TYPES).to_numpy()
    oligo_xy = adata.obs.loc[oligo_mask, ["x_centroid", "y_centroid"]].to_numpy(
        dtype=float
    )
    roi = ROI_MANUAL.get(sid) or find_roi(
        tf_tx[["x_location", "y_location"]].to_numpy(),
        tf_tx["feature_name"].to_numpy(),
        oligo_xy,
        hull,
        centroid,
        bx,
        ROI_SIZE_UM,
        ROI_MIN_OLIGO,
        ROI_TARGET_RNORM,
        ROI_RNORM_SIGMA,
    )
    cb = pd.read_parquet(
        raw / "cell_boundaries.parquet", columns=["cell_id", "vertex_x", "vertex_y"]
    )
    rnorm = _radial_norm(
        roi[0] + ROI_SIZE_UM / 2, roi[1] + ROI_SIZE_UM / 2, centroid, hull
    )
    d = {
        "adata": adata,
        "celltype": celltype,
        "txb": txb,
        "cb": cb,
        "roi": roi,
        "rnorm": rnorm,
    }
    _PREP_CACHE[key] = d
    return d


def _rows_for(grouping_rows, n_tf):
    """First n_tf//3 (row_label, [genes]) rows of a grouping (3 for top9, 2 for top6)."""
    return grouping_rows[: n_tf // 3]


def _genes_of(rows):
    return [g for _lbl, gs in rows for g in gs]


def figure_B_single(sid: str, grouping_rows, n_tf: int):
    rows = _rows_for(grouping_rows, n_tf)
    d = _prep_sample_B(sid, _genes_of(rows))
    w, h_full = FIG_SIZES["B_col"]
    h = h_full * (len(rows) / 3.0)
    # grid: [plot | TF legend]; legend lives outside the data panel.
    lw_ratio = LEGEND_COL_IN / w
    fig = plt.figure(figsize=(w + LEGEND_COL_IN, h))
    gs = fig.add_gridspec(
        len(rows), 2, width_ratios=[1, lw_ratio], hspace=0.06, wspace=0.02
    )
    for i, (_label, row_genes) in enumerate(rows):
        ax = fig.add_subplot(gs[i, 0])
        _draw_roi_panel(
            ax,
            d["roi"],
            d["cb"],
            d["celltype"],
            d["adata"],
            d["txb"],
            row_genes,
            ROI_SIZE_UM,
        )
        _row_tf_legend(fig.add_subplot(gs[i, 1]), row_genes)
        if i == len(rows) - 1:
            add_scalebar(ax, 50.0, "50 µm")
    fig.suptitle(SAMPLE_LABELS.get(sid, sid), fontsize=FONT_TITLE_PT)
    fig.subplots_adjust(top=0.95, bottom=0.02, left=0.02, right=0.98)
    print(f"    {SAMPLE_LABELS.get(sid, sid)}: ROI depth r_norm={d['rnorm']:.2f}")
    return fig, d["roi"]


def figure_B_pair(young_sid: str, old_sid: str, grouping_rows, n_tf: int):
    rows = _rows_for(grouping_rows, n_tf)
    genes = _genes_of(rows)
    dy = _prep_sample_B(young_sid, genes)
    do = _prep_sample_B(old_sid, genes)
    w, h_full = FIG_SIZES["B_col"]
    h = h_full * (len(rows) / 3.0)
    # grid per row: [Young | Old | TF legend]; legend to the RIGHT of the two plots.
    lw_ratio = LEGEND_COL_IN / w
    fig = plt.figure(figsize=(2 * w + LEGEND_COL_IN, h))
    gs = fig.add_gridspec(
        len(rows), 3, width_ratios=[1, 1, lw_ratio], hspace=0.06, wspace=0.03
    )
    first_ax = {}
    for i, (_label, row_genes) in enumerate(rows):
        for j, d in enumerate((dy, do)):
            ax = fig.add_subplot(gs[i, j])
            _draw_roi_panel(
                ax,
                d["roi"],
                d["cb"],
                d["celltype"],
                d["adata"],
                d["txb"],
                row_genes,
                ROI_SIZE_UM,
            )
            first_ax[(i, j)] = ax
        _row_tf_legend(fig.add_subplot(gs[i, 2]), row_genes)
    # Short column headers (specific sample IDs are in the filename).
    first_ax[(0, 0)].set_title("Young", fontsize=FONT_TITLE_PT)
    first_ax[(0, 1)].set_title("Old", fontsize=FONT_TITLE_PT)
    add_scalebar(first_ax[(len(rows) - 1, 0)], 50.0, "50 µm")
    fig.subplots_adjust(top=0.95, bottom=0.02, left=0.02, right=0.98)
    return fig


def figure_B_legend():
    """Standalone Fig B key, sized to sit UNDER Panel B (not the whole figure).

    Stacked so it stays ~Panel-B wide: top row = cell-border rings (the lineage call),
    below = lineage-marker transcript dots, coloured by stage, with the gene names.
    Emitted once for the illustrator to place; keeps the data panels legend-free.
    """
    genes_for = lambda st: [g for g in MARKER_GENES if MARKER_STAGE[g] == st]  # noqa: E731
    # ~Panel B pair width (2 cols + legend col).
    b_width = 2 * FIG_SIZES["B_col"][0] + LEGEND_COL_IN
    fig, (axT, axB) = plt.subplots(
        2, 1, figsize=(b_width, 1.5), gridspec_kw={"height_ratios": [1.0, 2.2]}
    )
    for a in (axT, axB):
        a.axis("off")

    # TOP: cell-border rings, one horizontal row.
    ring_handles = [
        plt.Line2D(
            [0], [0], marker="o", ls="", mfc="none", mec=LINEAGE_COLOR[k], mew=1.4, ms=8
        )
        for k in ("OPC", "Intermediate_Oligo", "Mature_Oligo")
    ]
    axT.legend(
        ring_handles,
        ["OPC", "Intermediate", "Mature"],
        loc="center",
        ncol=3,
        title="Cell border (lineage call)",
        frameon=False,
        fontsize=FONT_PT,
        title_fontsize=FONT_PT,
        handletextpad=0.3,
        columnspacing=0.9,
        borderpad=0.1,
    )

    # BELOW: lineage-marker dots + gene names, single column (gene lists are long).
    dot_handles, dot_labels = [], []
    for lbl, stage in (
        ("pan-lineage", "pan"),
        ("OPC", "OPC"),
        ("Intermediate", "Intermediate"),
        ("Mature", "Mature"),
    ):
        dot_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="o",
                ls="",
                mfc=MARKER_STAGE_COLORS[stage],
                mec="none",
                ms=6,
            )
        )
        dot_labels.append(f"{lbl}: {', '.join(genes_for(stage))}")
    axB.legend(
        dot_handles,
        dot_labels,
        loc="center",
        ncol=1,
        title="Lineage-marker transcripts (dots)",
        frameon=False,
        fontsize=FONT_PT,
        title_fontsize=FONT_PT,
        handletextpad=0.3,
        labelspacing=0.4,
        borderpad=0.1,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.03, hspace=0.1)
    return fig


# ------------------------- Figure C -------------------------
def figure_C():
    fig, _agg, _meta = dotplot.build_dotplot(
        figsize=FIG_SIZES["C"],
        dot_size_range=DOTPLOT_DOT_RANGE,
        cbar_height_frac=DOTPLOT_CBAR_FRAC,
        title_y=0.95,
        top=0.925,
        font_pt=FONT_PT,
    )
    return fig


# ------------------------- Figure D (Sox8 DE barplot) -------------------------
def _padj_star(padj) -> str:
    """Significance marker from a BH-adjusted p-value."""
    if padj is None or pd.isna(padj):
        return "n.s."
    if padj < 0.001:
        return "***"
    if padj < 0.01:
        return "**"
    if padj < 0.05:
        return "*"
    return "n.s."


def figure_D():
    """Panel D: Sox8 differential-expression barplot.

    Mean transcript counts per cell for Sox8 in Young vs Old, across the three
    oligo-lineage cell types. Bars = cohort mean with SEM error bars; individual
    sample means overlaid as dots; padj significance (from the canonical
    outlier-removed DE) marked per cell type. The DE-excluded outlier sample is
    dropped so the bar means and the padj share one sample set.
    """
    expr = pd.read_csv(EXPRESSION_CSV)
    de = pd.read_csv(DE_CSV)
    de_g = de[de["gene"] == BARPLOT_GENE]
    excluded = set(de_g["excluded_sample"].dropna().astype(str).unique())
    padj_by_ct = de_g.set_index("cell_type")["padj"].to_dict()

    e = expr[
        (expr["gene"] == BARPLOT_GENE)
        & (expr["cell_type"].isin(KEY_CELLTYPES))
        & (~expr["sample_id"].astype(str).isin(excluded))
    ].copy()

    ages = ["Young", "Old"]
    x = np.arange(len(KEY_CELLTYPES))
    width = 0.38
    fig, ax = plt.subplots(figsize=FIG_SIZES["D"])
    bar_tops = np.zeros(
        len(KEY_CELLTYPES)
    )  # tallest element per group, for star placement

    for i, age in enumerate(ages):
        means, sems = [], []
        for ct in KEY_CELLTYPES:
            vals = e[(e["cell_type"] == ct) & (e["age_group"] == age)][
                "mean_counts_per_cell"
            ].to_numpy()
            means.append(float(np.mean(vals)) if len(vals) else 0.0)
            sems.append(
                float(np.std(vals, ddof=1) / np.sqrt(len(vals)))
                if len(vals) > 1
                else 0.0
            )
        offset = (i - 0.5) * width  # Young left, Old right
        ax.bar(
            x + offset,
            means,
            width,
            yerr=sems,
            color=plot_style.GROUP_COLORS[age],
            label=age,
            edgecolor="black",
            linewidth=0.4,
            error_kw={
                "elinewidth": 0.6,
                "capsize": 1.6,
                "capthick": 0.6,
                "ecolor": "#333333",
            },
            zorder=2,
        )
        for j, ct in enumerate(KEY_CELLTYPES):
            vals = e[(e["cell_type"] == ct) & (e["age_group"] == age)][
                "mean_counts_per_cell"
            ].to_numpy()
            if not len(vals):
                continue
            jit = (
                np.linspace(-1, 1, len(vals)) * 0.25 * width
                if len(vals) > 1
                else np.zeros(1)
            )
            ax.scatter(
                np.full(len(vals), x[j] + offset) + jit,
                vals,
                s=4.0,
                facecolor="white",
                edgecolor="#333333",
                linewidth=0.4,
                zorder=3,
            )
            bar_tops[j] = max(bar_tops[j], means[j] + sems[j], float(vals.max()))

    ymax = float(bar_tops.max())
    for j, ct in enumerate(KEY_CELLTYPES):
        star = _padj_star(padj_by_ct.get(ct))
        y = bar_tops[j] + 0.04 * ymax
        xl, xr = x[j] - 0.5 * width, x[j] + 0.5 * width
        ax.plot(
            [xl, xl, xr, xr],
            [y, y + 0.02 * ymax, y + 0.02 * ymax, y],
            lw=0.6,
            color="#333333",
        )
        ax.text(
            x[j],
            y + 0.03 * ymax,
            star,
            ha="center",
            va="bottom",
            fontsize=FONT_PT,
            color="#333333",
        )

    ax.set_xticks(x)
    ax.set_xticklabels([CELLTYPE_LABELS_SHORT[ct] for ct in KEY_CELLTYPES])
    ax.set_ylabel(f"{BARPLOT_GENE} counts / cell")
    ax.set_ylim(0, ymax * 1.30)
    ax.set_yticks([0, 2, 4])
    ax.margins(x=0.10)
    # Legend as a single row ABOVE the axes so it never overlaps bars in the small panel;
    # no in-panel title (the panel letter/caption is added in the composite).
    ax.legend(
        frameon=False,
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        handlelength=0.9,
        handletextpad=0.4,
        columnspacing=0.9,
        borderpad=0.0,
        fontsize=FONT_PT,
    )
    fig.tight_layout()
    return fig


# --- Supplementary (former Panel D): human TF spatial maps ---
def figure_suppl_human(ncol: int, figsize):
    fig, info = humanmaps.build_spatial_maps(
        sample_label=HUMAN_SAMPLE,
        genes=HUMAN_GENES,
        block=HUMAN_BLOCK,
        ncol=ncol,
        figsize=figsize,
        point_size=HUMAN_POINT_SIZE,
        alpha=HUMAN_POINT_ALPHA,
        match_aspect=True,
        title="Human TF spatial distribution",
    )
    print(f"  Suppl human block: {info['block_used']} | split: {info['split_info']}")
    return fig


# ------------------------- orchestration -------------------------
def run_A(samples):
    for sid in samples:
        fig = figure_A(sid)
        plot_style.save(
            fig, f"panelA_morphology_{SAMPLE_LABELS.get(sid, sid)}", fig_dir=OUT_DIR
        )


def run_B(samples, tf_set, grouping_name):
    rows = B_GROUPINGS[grouping_name]
    n_tf = 6 if tf_set == "top6" else 9
    tag = f"group{grouping_name}_{tf_set}"
    print(f"  {tag}: rows = {[lbl for lbl, _ in rows[: n_tf // 3]]}")
    for sid in samples:
        fig, _ = figure_B_single(sid, rows, n_tf)
        plot_style.save(
            fig, f"panelB_{tag}_{SAMPLE_LABELS.get(sid, sid)}", fig_dir=OUT_DIR
        )
    # every Young x Old pair whose members are both in this run
    for y, o in B_PAIRS:
        if y not in samples or o not in samples:
            continue
        fig = figure_B_pair(y, o, rows, n_tf)
        plot_style.save(
            fig,
            f"panelB_{tag}_pair_{SAMPLE_LABELS.get(y, y)}_{SAMPLE_LABELS.get(o, o)}",
            fig_dir=OUT_DIR,
        )


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("figure", choices=["A", "B", "C", "D", "suppl", "all"])
    p.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="sample IDs for A/B (default: cohort-9)",
    )
    p.add_argument("--tf-set", choices=["top9", "top6", "both"], default="top9")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    apply_manuscript_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ab_samples = args.samples if args.samples else (AB_YOUNG + AB_OLD)
    tf_sets = ["top9", "top6"] if args.tf_set == "both" else [args.tf_set]

    if args.figure in ("A", "all"):
        print("Figure A (morphology)…")
        run_A(ab_samples)
    if args.figure in ("B", "all"):
        print("Figure B (single-cell expression)…")
        for grp in B_GROUPINGS_TO_RENDER:
            for ts in tf_sets:
                run_B(ab_samples, ts, grp)
        plot_style.save(figure_B_legend(), "panelB_legend_lineage", fig_dir=OUT_DIR)
    if args.figure in ("C", "all"):
        print("Figure C (dot plot)…")
        plot_style.save(figure_C(), "panelC_dotplot", fig_dir=OUT_DIR)
    if args.figure in ("D", "all"):
        print("Figure D (Sox8 DE barplot)…")
        plot_style.save(figure_D(), "panelD_sox8_de_barplot", fig_dir=OUT_DIR)
    if args.figure in ("suppl", "all"):
        print("Supplementary (human TF spatial maps; former Panel D)…")
        plot_style.save(
            figure_suppl_human(
                len(HUMAN_GENES), FIG_SIZES["SUPPL_HUMAN"]
            ),  # single row
            f"supplD_human_tf_maps_block{HUMAN_BLOCK}",
            fig_dir=OUT_DIR,
        )
        plot_style.save(
            figure_suppl_human(1, FIG_SIZES["SUPPL_HUMAN_vertical"]),
            f"supplD_human_tf_maps_block{HUMAN_BLOCK}_vertical",
            fig_dir=OUT_DIR,
        )
    print(f"done -> {OUT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
