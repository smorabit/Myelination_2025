#!/usr/bin/env python3
"""Aggregated expression-evidence dot plot — publication summary companion.

Reads `data/expression_evidence/expression_per_sample_celltype_<strategy>.csv`
(produced by `scripts/compute_expression_evidence.py`) and renders one
publication-ready dot plot:
  - dot size  = mean detection fraction (aggregated across samples within
                each age × cell_type group)
  - dot color = mean transcripts per cell (mean of per-sample means, same
                `outlier_rm + pcmt=10` cohort as the canonical DE)
  - rows      = 47 panel genes (9 markers + 9 primary TFs + 29 targets/inducers)
                + 3 negative controls. Order: TF blocks (Bach2 → Klk6) with
                primary TF at the head of each block → markers → controls.
  - columns   = 6 aggregated pseudobulks ordered by cell type with Young/Old
                adjacent within each cell-type block:
                [OPC-Y, OPC-O, Int-Y, Int-O, Mat-Y, Mat-O].
                Old 2 excluded by default (5 Young + 4 Old samples averaged).

This is the publication-ready companion to the per-sample heatmap produced by
`plot_expression_heatmap.py`. The per-sample heatmap remains the QC view; this
script provides a compact ~6.5 × 13 in main-text panel.

Output:
  docs/images/expression_evidence_dotplot_<strategy>_counts_no_old2.png + .pdf
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
from matplotlib.colors import LogNorm

import plot_style

REPO_ROOT = repo_root()
IN_DIR = REPO_ROOT / "data" / "expression_evidence"
OUT_DIR = REPO_ROOT / "docs" / "images"
DEFAULT_STRATEGY = "stringent_p5"

OLIGO_TYPES = ["OPC", "Intermediate_Oligo", "Mature_Oligo"]
CELL_TYPE_LABELS = {
    "OPC": "OPC",
    "Intermediate_Oligo": "Intermediate",
    "Mature_Oligo": "Mature",
}
AGE_ORDER = ["Young", "Old"]

TF_GROUP_ORDER = [
    "Bach2",
    "Elf2",
    "Foxk2",
    "Bhlhe41",
    "Nr6a1",
    "Sox8",
    "Stat3",
    "Sox5",
    "Klk6",
]

# Markers ordered by oligo-lineage cell-type assignment so visual scanning
# top-to-bottom matches the cell-type column order (OPC → Int → Mature).
# Pan-lineage TFs first (Olig2, Sox10), then OPC-specific (Pdgfra, Ptprz1,
# Pcdh15), then the single committed-progenitor / new-oligo marker (Enpp6),
# then mature myelinating markers (Mbp, Mog, Opalin grouped together).
MARKER_ORDER = [
    "Olig2",
    "Sox10",
    "Pdgfra",
    "Ptprz1",
    "Pcdh15",
    "Enpp6",
    "Mbp",
    "Mog",
    "Opalin",
]

GENE_CLASS_TO_TF_GROUP = {
    "Bach2_target": "Bach2",
    "Elf2_target": "Elf2",
    "Foxk2_target": "Foxk2",
    "Bhlhe41_target": "Bhlhe41",
    "Nr6a1_target": "Nr6a1",
    "Sox8_target": "Sox8",
    "Stat3_target": "Stat3",
    "Sox5_inducer": "Sox5",
    "Klk6_inducer": "Klk6",
    "marker": "marker",
    "negative_control": "negative ctrl",
}

# Colourblind-safe TF-group palette, shared across all figures (Paul Tol 'muted').
TF_GROUP_COLORS = plot_style.TF_GROUP_COLORS

CELL_TYPE_COLORS = {
    "OPC": plot_style.CELLTYPE_COLORS["OPC"],
    "Intermediate_Oligo": plot_style.CELLTYPE_COLORS["Intermediate"],
    "Mature_Oligo": plot_style.CELLTYPE_COLORS["Mature"],
}
CELL_TYPE_TEXT_COLOR = {
    "OPC": "white",
    "Intermediate_Oligo": "white",
    "Mature_Oligo": "white",
}
AGE_COLORS = plot_style.GROUP_COLORS

DROPPED_CLASSES = {"pancreas_control", "pancreas_or_stress"}

# Dot-size scaling: s is points² in matplotlib scatter (diameter = sqrt(s)).
# Maximise the range across 0–100 % so "no expression" vs "fully expressed"
# is visually unambiguous: 0 % → ~1.4 pt pinpoint; 100 % → ~13 pt diameter
# (fills the ~15 pt tall tile in an 11.5 in tall canvas without overflow).
DOT_SIZE_MIN = 2
DOT_SIZE_RANGE = 170

# Colour-scale bounds for the log normaliser. vmin=0.05 puts the
# control_probe noise floor (0.27 transcripts/cell median) at ~13 % of the
# colormap and a meaningful "moderate expression" value of 2 at ~70 %.
# Values < vmin clip to the dark end (essentially indistinguishable from 0).
COLOR_VMIN = 0.05
COLOR_VMAX = 10.0
# Reference noise floor (median of non-zero control_probe means from the
# signal-to-noise analysis). Marked on the colorbar as the "platform noise"
# boundary so the colorbar itself anchors the above-baseline argument.
NOISE_FLOOR = 0.27


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--strategy",
        default=DEFAULT_STRATEGY,
        help=f"Strategy suffix (default: {DEFAULT_STRATEGY}).",
    )
    parser.add_argument(
        "--include-old2",
        action="store_true",
        help="Include Old 2 in the aggregation (default: exclude, matching canonical cohort).",
    )
    return parser.parse_args()


def assign_tf_group(row: pd.Series) -> str:
    if row["gene_class"] == "primary_TF":
        return row["gene"]
    return GENE_CLASS_TO_TF_GROUP.get(row["gene_class"], "other")


def gene_block_order(genes_df: pd.DataFrame) -> list[str]:
    """TF blocks (Bach2 → Klk6) → markers → negative controls."""
    out: list[str] = []
    for grp in TF_GROUP_ORDER:
        block = genes_df[genes_df["tf_group"] == grp]
        primary = block[block["gene_class"] == "primary_TF"]["gene"].tolist()
        non_primary = sorted(
            block[block["gene_class"] != "primary_TF"]["gene"].tolist()
        )
        out.extend(primary + non_primary)
    # Markers in developmental / cell-type order (not alphabetical) so the
    # rows scan vertically in the same direction as the OPC → Int → Mature
    # column layout.
    markers_present = set(genes_df[genes_df["gene_class"] == "marker"]["gene"])
    markers_ordered = [g for g in MARKER_ORDER if g in markers_present]
    # Defensive: append any markers in the data but missing from MARKER_ORDER
    # (alphabetical, so the failure mode is "extra marker shows up at the end"
    # rather than "silently dropped").
    leftover = sorted(markers_present - set(markers_ordered))
    out.extend(markers_ordered + leftover)
    nc_order = ["control_probe", "genomic_control", "unassigned_codeword"]
    nc_present = genes_df[genes_df["gene_class"] == "negative_control"]["gene"].tolist()
    out.extend([g for g in nc_order if g in nc_present])
    return out


def hex_to_rgb01(c: str) -> tuple[float, float, float]:
    return (int(c[1:3], 16) / 255, int(c[3:5], 16) / 255, int(c[5:7], 16) / 255)


def build_dotplot(
    strategy: str = DEFAULT_STRATEGY,
    include_old2: bool = False,
    figsize: tuple[float, float] = (7.5, 11.5),
    dot_size_min: float = DOT_SIZE_MIN,
    dot_size_range: float = DOT_SIZE_RANGE,
    cbar_height_frac: float = 1.0,
    title_y: float = 0.985,
    top: float | None = None,
    font_pt: float | None = None,
):
    """Build the aggregated expression dot plot and return (fig, agg, meta).

    Parametrized so the manuscript-panel module can request a smaller / narrower
    canvas (and a proportionally smaller dot-size range so tiles don't overlap)
    without duplicating the layout code. `cbar_height_frac` shrinks the colorbar to
    a centred fraction of the panel height (1.0 = full height, as in the report);
    `title_y` / `top` tighten the gap between the suptitle and the dot grid. `font_pt`
    sets one uniform font size (used by the manuscript panel so text matches the other
    panels); left None it keeps the original report sizes. The CLI path uses the
    defaults, so its output is byte-for-byte the same. The caller is responsible for
    `plot_style.apply()` and for saving.
    """
    # Resolve the handful of text sizes: report defaults, or all tied to font_pt.
    if font_pt is None:
        fs_ct, fs_age, fs_row = 11, 8, 7
        fs_leg_title, fs_leg = 9, 8
        fs_cbar_label, fs_cbar_tick, fs_noise = 9, 8, 7
    else:
        # Cell-type strip label ("Intermediate") is the longest text over the narrowest
        # block, so drop it a step to fit; everything else stays at font_pt.
        fs_ct = max(font_pt - 1, 4)
        fs_age = fs_row = font_pt
        fs_leg_title = fs_leg = font_pt
        fs_cbar_label = fs_cbar_tick = font_pt
        fs_noise = max(font_pt - 1, 4)
    df = pd.read_csv(IN_DIR / f"expression_per_sample_celltype_{strategy}.csv")
    df["tf_group"] = df.apply(assign_tf_group, axis=1)
    df = df[~df["gene_class"].isin(DROPPED_CLASSES)].reset_index(drop=True)

    if not include_old2:
        df = df[~df["sample_is_old2"]].reset_index(drop=True)
        cohort_note = "Old 2 excluded"
    else:
        cohort_note = "all 10 samples"

    # --- Aggregate: simple mean of per-sample means ---
    agg = (
        df.groupby(["gene", "age_group", "cell_type"], observed=True)
        .agg(
            mean_counts_per_cell=("mean_counts_per_cell", "mean"),
            detection_fraction=("detection_fraction", "mean"),
            gene_class=("gene_class", "first"),
            tf_group=("tf_group", "first"),
        )
        .reset_index()
    )

    # --- Row + column orders ---
    gene_meta = agg.drop_duplicates("gene")[
        ["gene", "gene_class", "tf_group"]
    ].reset_index(drop=True)
    row_order = gene_block_order(gene_meta)
    gene_to_tf_group = dict(zip(gene_meta["gene"], gene_meta["tf_group"]))

    col_order: list[tuple[str, str]] = []
    for ct in OLIGO_TYPES:
        for age in AGE_ORDER:
            col_order.append((ct, age))
    col_keys = [f"{ct}__{age}" for ct, age in col_order]

    agg["col_key"] = agg["cell_type"] + "__" + agg["age_group"]
    pivot_color = agg.pivot(
        index="gene", columns="col_key", values="mean_counts_per_cell"
    ).reindex(index=row_order, columns=col_keys)
    pivot_size = agg.pivot(
        index="gene", columns="col_key", values="detection_fraction"
    ).reindex(index=row_order, columns=col_keys)

    n_rows = len(row_order)
    n_cols = len(col_keys)

    # --- Figure layout ---
    # 7.5 × 11.5 in canvas: narrower than the per-sample heatmap because we
    # collapsed 27 → 6 columns, but wide enough that the "Intermediate"
    # cell-type label has comfortable margin. Height reduced from 13 → 11.5
    # so each tile is ~15 pt tall (still ample for 7 pt italic row labels)
    # — slightly more compact for journal-panel use.
    fig = plt.figure(figsize=figsize)
    gs_kw = {} if top is None else {"top": top}
    # ncols=6: [row-strip, dots, gap, size-legend, gap, colorbar]. The extra gap column
    # between the size legend and the colorbar keeps them from touching on the narrow
    # manuscript canvas.
    gs = fig.add_gridspec(
        nrows=3,
        ncols=6,
        height_ratios=[1.8, 0.7, n_rows],
        width_ratios=[0.6, n_cols, 0.4, 1.1, 0.55, 0.4],
        hspace=0.05,
        wspace=0.08,
        **gs_kw,
    )
    ax_ct = fig.add_subplot(gs[0, 1])
    ax_age = fig.add_subplot(gs[1, 1])
    ax_row = fig.add_subplot(gs[2, 0])
    ax = fig.add_subplot(gs[2, 1])
    ax_size = fig.add_subplot(gs[2, 3])
    if cbar_height_frac >= 1.0:
        ax_cbar = fig.add_subplot(gs[2, 5])
    else:
        # Centre a shorter colorbar in the right column: an over-long bar reads as
        # visual clutter, especially on the narrower manuscript canvas.
        pad = (1.0 - cbar_height_frac) / 2.0
        sub = gs[2, 5].subgridspec(
            3, 1, height_ratios=[pad, cbar_height_frac, pad], hspace=0.0
        )
        ax_cbar = fig.add_subplot(sub[1, 0])

    # --- Main panel: scatter ---
    cmap = plt.get_cmap("viridis")
    # Log scale so the biologically meaningful 0–2 transcripts/cell range
    # gets ample colormap range instead of being compressed into the darkest
    # 20 % of viridis. Values below vmin clip to dark (effectively "0").
    norm = LogNorm(vmin=COLOR_VMIN, vmax=COLOR_VMAX, clip=True)

    xs: list[int] = []
    ys: list[int] = []
    ss: list[float] = []
    cs: list[tuple[float, float, float, float]] = []
    n_missing = 0
    for i, gene in enumerate(row_order):
        for j in range(n_cols):
            v = pivot_color.iloc[i, j]
            d = pivot_size.iloc[i, j]
            if pd.isna(v) or pd.isna(d):
                n_missing += 1
                continue
            xs.append(j)
            ys.append(i)
            ss.append(dot_size_min + dot_size_range * float(d))
            # Clip to vmin floor for log-norm; values of 0 become the
            # darkest cmap colour instead of NaN.
            cs.append(cmap(norm(max(float(v), COLOR_VMIN))))

    ax.scatter(xs, ys, s=ss, c=cs, edgecolor="black", linewidth=0.3)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)

    # Light vertical dividers between cell-type blocks (between cols 1|2 and 3|4)
    for x in [1.5, 3.5]:
        ax.axvline(x=x, color="lightgrey", linewidth=0.8, zorder=0)
    # Light horizontal dividers between TF-group blocks
    row_groups = [gene_to_tf_group[g] for g in row_order]
    for k in range(1, n_rows):
        if row_groups[k] != row_groups[k - 1]:
            ax.axhline(y=k - 0.5, color="lightgrey", linewidth=0.6, zorder=0)

    # --- Row strip + gene labels on its left ---
    row_strip = np.zeros((n_rows, 1, 3))
    for i, g in enumerate(row_order):
        row_strip[i, 0] = hex_to_rgb01(
            TF_GROUP_COLORS.get(gene_to_tf_group[g], "#cccccc")
        )
    ax_row.imshow(row_strip, aspect="auto", interpolation="nearest", origin="upper")
    ax_row.set_yticks(range(n_rows))
    ax_row.set_yticklabels(row_order, fontsize=fs_row, fontstyle="italic")
    ax_row.yaxis.tick_left()
    ax_row.set_xticks([])
    ax_row.tick_params(left=False, right=False, length=0, pad=2)
    for sp in ax_row.spines.values():
        sp.set_visible(False)
    for k in range(1, n_rows):
        if row_groups[k] != row_groups[k - 1]:
            ax_row.axhline(y=k - 0.5, color="white", linewidth=1.2)

    # --- Cell-type strip (top) ---
    ct_strip = np.zeros((1, n_cols, 3))
    for j, (ct, _) in enumerate(col_order):
        ct_strip[0, j] = hex_to_rgb01(CELL_TYPE_COLORS[ct])
    ax_ct.imshow(ct_strip, aspect="auto", interpolation="nearest", origin="upper")
    for k, ct in enumerate(OLIGO_TYPES):
        center = 2 * k + 0.5
        ax_ct.text(
            center,
            0,
            CELL_TYPE_LABELS[ct],
            ha="center",
            va="center",
            fontsize=fs_ct,
            color=CELL_TYPE_TEXT_COLOR[ct],
            fontweight="bold",
        )
    ax_ct.set_xticks([])
    ax_ct.set_yticks([])
    for sp in ax_ct.spines.values():
        sp.set_visible(False)
    for x in [1.5, 3.5]:
        ax_ct.axvline(x=x, color="white", linewidth=2.0)

    # --- Age strip ---
    age_strip = np.zeros((1, n_cols, 3))
    for j, (_, age) in enumerate(col_order):
        age_strip[0, j] = hex_to_rgb01(AGE_COLORS[age])
    ax_age.imshow(age_strip, aspect="auto", interpolation="nearest", origin="upper")
    for j, (_, age) in enumerate(col_order):
        ax_age.text(
            j,
            0,
            age,
            ha="center",
            va="center",
            fontsize=fs_age,
            color="black",
            fontweight="bold",
        )
    ax_age.set_xticks([])
    ax_age.set_yticks([])
    for sp in ax_age.spines.values():
        sp.set_visible(False)
    for x in [1.5, 3.5]:
        ax_age.axvline(x=x, color="white", linewidth=2.0)

    # --- Size legend ---
    ax_size.set_axis_off()
    ax_size.text(
        0.5,
        0.97,
        "Detection\n(% cells)",
        ha="center",
        va="top",
        fontsize=fs_leg_title,
        fontweight="bold",
        transform=ax_size.transAxes,
    )
    ref_pcts = [0.10, 0.25, 0.50, 0.75, 1.00]
    n_ref = len(ref_pcts)
    y_positions = np.linspace(0.85, 0.55, n_ref)
    for pct, y in zip(ref_pcts, y_positions):
        ax_size.scatter(
            0.25,
            y,
            s=dot_size_min + dot_size_range * pct,
            c="lightgrey",
            edgecolor="black",
            linewidth=0.3,
            transform=ax_size.transAxes,
        )
        ax_size.text(
            0.55,
            y,
            f"{int(pct * 100)}%",
            ha="left",
            va="center",
            fontsize=fs_leg,
            transform=ax_size.transAxes,
        )

    # --- Colorbar (log scale, with noise-floor reference tick) ---
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = plt.colorbar(sm, cax=ax_cbar)
    cb.set_label("Mean counts per cell\n(log scale)", fontsize=fs_cbar_label)
    cb.set_ticks([0.05, 0.1, 0.3, 1.0, 3.0, 10.0])
    cb.set_ticklabels(["≤0.05", "0.1", "0.3", "1", "3", "≥10"])
    cb.ax.tick_params(labelsize=fs_cbar_tick)
    # Mark the platform noise floor (control_probe non-zero median) with a
    # thin horizontal line + side label so a reader can visually anchor
    # "above-baseline" without reading text. NB: matplotlib's vertical
    # colorbars map data range to y-coordinates of cb.ax, so axhline works.
    cb.ax.axhline(y=NOISE_FLOOR, color="red", linewidth=0.9, linestyle="--")
    # Label to the LEFT of the bar so it doesn't collide with the right-side tick labels.
    cb.ax.text(
        -0.15,
        NOISE_FLOOR,
        "noise\nfloor ",
        transform=cb.ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=fs_noise,
        color="red",
    )

    # --- Title ---
    fig.suptitle("Mean transcripts per cell: Young vs Old", y=title_y)

    return fig, agg, {"cohort_note": cohort_note, "n_missing": n_missing}


def main() -> int:
    args = parse_args()
    plot_style.apply()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, agg, meta = build_dotplot(
        strategy=args.strategy, include_old2=args.include_old2
    )

    suffix = "" if args.include_old2 else "_no_old2"
    out_base = OUT_DIR / f"expression_evidence_dotplot_{args.strategy}_counts{suffix}"
    plot_style.save(fig, out_base.name)
    print(f"  cohort: {meta['cohort_note']}")
    if meta["n_missing"]:
        print(
            f"  warning: {meta['n_missing']} (gene, age, cell_type) tiles had no data"
        )

    # --- Spot-check ---
    sox8 = agg[
        (agg["gene"] == "Sox8")
        & (agg["cell_type"] == "Intermediate_Oligo")
        & (agg["age_group"] == "Young")
    ]
    if not sox8.empty:
        v = sox8.iloc[0]["mean_counts_per_cell"]
        d = sox8.iloc[0]["detection_fraction"]
        print(f"  spot-check Sox8/Intermediate/Young: mean_counts={v:.3f}, det={d:.3f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
