#!/usr/bin/env python3
"""Main analysis: all-pairs cross-PCF co-localisation at cell scale, with statistics.

Segmentation-free spatial co-localisation of transcription factors with cell-type
markers and with each other, at a single cell-scale radius, tested for significance
by label permutation. Run on both good human sections.

Effect size (density-aware): among a gene's transcript neighbours within R um, is a
second gene over-represented relative to its overall abundance? Both parts scale with
local density, so it is not confounded by the fact that abundant myelin genes make
white matter dense (which a global null wrongly reads as TFs being "depleted").

    enrich_ij   = (fraction of gene i's neighbours that are gene j) / (gene j global fraction)
    log2_enrich = log2(enrich_ij)
    p           empirical two-sided p from a label-permutation null (N_PERM=1000; enrich
                centres on 1 under the null), resolving to 1/(1+N_PERM). A normal-approx
                p (p_znorm) is reported alongside for reference only.
                q = Benjamini-Hochberg over TF x gene pairs, on the empirical p.

log2_enrich > 0 = gene j enriched in gene i's neighbourhood; < 0 = excluded.

Outputs:
  data/xenium_human/cross_pcf/cross_pcf_pairs_<sample>.csv   (TF x all genes, per sample)
  data/xenium_human/cross_pcf/cross_pcf_summary.csv          (significant in BOTH sections)
  docs/images/human_tf_cross_pcf_markers.png / .pdf          (TF x marker heatmap)
  docs/images/human_tf_cross_pcf_tftf.png / .pdf             (TF-TF heatmap)

Run:
  HUMAN_TX_DIR=/path/to/human/transcripts python scripts/compute_cross_pcf.py
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


import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import norm

import human_tf_common as C
import plot_style
from spatial_common import alpha_shape

TABLE_DIR = C.REPO_ROOT / "data" / "xenium_human" / "cross_pcf"
FIG_DIR = C.REPO_ROOT / "docs" / "images"
R = 10.0  # um, cell scale
N_PERM = 1000  # label shuffles; empirical two-sided p resolves to 1/(1+N_PERM)
HULL_SUBSAMPLE = 40_000
SEED = 0

# Cell-type marker groups shown on the TF x marker heatmap. Restricted to markers
# actually on the hBrain panel (verified against human_tf_evidence.csv); the panel is
# TF-focused, so classic markers like RBFOX3, GFAP, CLDN5, C1QA are absent and cannot
# be added. SLC17A6 (neuron), ITGAM (microglia), MAL/ERMN (mature myelin) are on-panel
# and included here to broaden the cell-type read-out.
MARKER_GROUPS = [
    ("OPC", ["PDGFRA", "PTPRZ1", "PCDH15", "CSPG4"]),
    ("Oligo (pan)", ["OLIG1", "OLIG2", "SOX10"]),
    ("Intermediate", ["ENPP6"]),
    (
        "Mature",
        ["MBP", "MOG", "OPALIN", "MAG", "PLP1", "CLDN11", "MOBP", "MAL", "ERMN"],
    ),
    ("Astrocyte", ["AQP4", "GJA1", "SOX9"]),
    ("Microglia", ["AIF1", "CX3CR1", "P2RY12", "ITGAM"]),
    ("Neuron", ["SLC17A7", "SLC17A6", "GAD1", "GAD2"]),
    ("Endothelial", ["PECAM1", "FLT1"]),
]
MARKERS = [g for _, gs in MARKER_GROUPS for g in gs]


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR q-values for a 1D array of p-values (NaNs preserved)."""
    q = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    qv = np.empty(n)
    qv[order] = ranked
    q[ok] = qv
    return q


def compute_sample(sample: dict, rng: np.random.Generator) -> dict:
    tx = C.load_transcripts(
        sample["dir"],
        columns=["feature_name", "qv", "is_gene", "x_location", "y_location"],
    )
    tx = tx[(tx["qv"] >= C.QV_MIN) & (tx["is_gene"])]
    xy = tx[["x_location", "y_location"]].to_numpy()

    sub = xy[rng.choice(len(xy), min(HULL_SUBSAMPLE, len(xy)), replace=False)]
    hull, _ = alpha_shape(sub)
    eroded = hull.buffer(-R)

    import shapely

    keep = shapely.contains_xy(eroded, xy[:, 0], xy[:, 1])
    xy = xy[keep]
    cats = pd.Categorical(tx["feature_name"].to_numpy()[keep])
    codes = cats.codes.astype(np.int64)
    names = list(cats.categories)
    n = len(names)

    pairs = cKDTree(xy).query_pairs(r=R, output_type="ndarray")
    pa, pb = pairs[:, 0], pairs[:, 1]

    counts = np.bincount(codes, minlength=n).astype(float)
    N = len(codes)
    global_frac = counts / N

    def enrich(cd: np.ndarray) -> np.ndarray:
        # symmetric co-occurrence counts, then density-aware compositional enrichment:
        # (fraction of gene i's neighbours that are gene j) / (gene j's global fraction).
        flat = (
            np.bincount(cd[pa] * n + cd[pb], minlength=n * n)
            .reshape(n, n)
            .astype(float)
        )
        m = flat + flat.T
        with np.errstate(divide="ignore", invalid="ignore"):
            return (m / m.sum(1, keepdims=True)) / global_frac

    obs = enrich(codes)
    with np.errstate(divide="ignore", invalid="ignore"):
        log2e = np.clip(
            np.log2(obs), -3.0, 3.0
        )  # full grid; zero co-occ -> -inf, floored

    # Label-permutation null: shuffle labels on the fixed pair list. The enrichment is
    # density-normalised, so under the null it centres on 1 (log2 0). Only the TF rows
    # are ever read downstream (tf_gene_table / cell), so collect the null for those
    # rows only -- (N_PERM, n_tf, n) float32 is ~11 MB at N_PERM=1000, vs ~400 MB for a
    # full n x n null. The canonical significance is an empirical two-sided p from this
    # null; the z (normal-approx) is kept alongside for reference only.
    idx = {nm: k for k, nm in enumerate(names)}
    tfs = [t for t in C.PRIMARY_TFS if t in idx]
    tf_rows = [idx[t] for t in tfs]

    null = np.empty((N_PERM, len(tf_rows), n), dtype=np.float32)
    perm = codes.copy()
    for k in range(N_PERM):
        rng.shuffle(perm)
        null[k] = enrich(perm)[tf_rows]

    obs_tf = obs[tf_rows, :]
    mean = null.mean(0)
    std = np.sqrt(np.maximum(null.var(0), 1e-12))
    with np.errstate(divide="ignore", invalid="ignore"):
        z_tf = (obs_tf - mean) / std
        # empirical two-sided p: how often a permuted enrichment deviates from the null
        # centre (log2 0) at least as far as the observed does. Resolves to 1/(1+N_PERM);
        # NaN comparisons are False, so this is well defined for zero co-occurrence too.
        obs_dev = np.abs(np.log2(obs_tf))
        null_dev = np.abs(np.log2(null))
    p_emp_tf = (1 + (null_dev >= obs_dev[None]).sum(0)) / (1 + N_PERM)
    p_z_tf = 2.0 * np.minimum(norm.sf(np.abs(z_tf)), 0.5)

    # scatter TF-row stats into full n x n grids (NaN elsewhere; only TF rows are read)
    z = np.full((n, n), np.nan)
    p = np.full((n, n), np.nan)
    p_znorm = np.full((n, n), np.nan)
    z[tf_rows, :] = z_tf
    p[tf_rows, :] = p_emp_tf
    p_znorm[tf_rows, :] = p_z_tf

    return {
        "label": sample["label"],
        "names": names,
        "idx": idx,
        "counts": counts.astype(int),
        "log2e": log2e,
        "z": z,
        "p": p,
        "p_znorm": p_znorm,
    }


def tf_gene_table(res: dict) -> pd.DataFrame:
    idx, names = res["idx"], res["names"]
    tfs = [t for t in C.PRIMARY_TFS if t in idx]
    rows = []
    for tf in tfs:
        i = idx[tf]
        for k, nm in enumerate(names):
            if k == i:
                continue
            rows.append(
                {
                    "sample_label": res["label"],
                    "tf": tf,
                    "gene": nm,
                    "gene_class": C.gene_class(nm),
                    "is_positive_control": C.is_positive_control(nm),
                    "n_gene": int(res["counts"][k]),
                    "log2_enrich": round(float(res["log2e"][i, k]), 3),
                    "z": round(float(res["z"][i, k]), 2),
                    "p": float(res["p"][i, k]),  # canonical: empirical two-sided
                    "p_znorm": float(res["p_znorm"][i, k]),  # reference: normal-approx
                }
            )
    df = pd.DataFrame(rows)
    df["q"] = bh_qvalues(df["p"].to_numpy())
    return df


def heatmap(mat, sig, rows_lbl, cols_lbl, title, fname, figsize, groups=None):
    finite = mat[np.isfinite(mat)]
    vmax = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 1.0
    vmax = vmax or 1.0
    mat = np.clip(mat, -vmax, vmax)  # keeps NaN (diagonal) as NaN
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cols_lbl)))
    ax.set_xticklabels(cols_lbl, rotation=90, fontsize=8)
    ax.set_yticks(range(len(rows_lbl)))
    ax.set_yticklabels(rows_lbl, fontsize=9)
    if groups:
        x = 0
        for name, gs in groups:
            ax.axvline(x - 0.5, color="white", linewidth=2)
            ax.text(
                x + len(gs) / 2 - 0.5, -0.55, name, ha="center", va="bottom", fontsize=7
            )
            x += len(gs)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if sig[i, j]:
                ax.text(j, i, "*", ha="center", va="center", fontsize=9)
    # pad the title so it sits clear above the group labels
    ax.set_title(title, fontsize=10, pad=24 if groups else 10)
    # colorbar last, and no subplots_adjust afterwards (that would stretch the
    # heatmap back over the colorbar); bbox_inches="tight" handles the margins.
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04).set_label(
        "log2 neighbourhood enrichment", fontsize=8
    )
    plot_style.save(fig, fname)


def main() -> int:
    global R
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--radius", type=float, default=10.0, help="neighbourhood radius in um"
    )
    R = ap.parse_args().radius
    suffix = "" if R == 10.0 else f"_{int(R)}um"

    plot_style.apply()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    good = [s for s in C.SAMPLES if not s["flagged"]]

    results, tables = {}, {}
    for s in good:
        res = compute_sample(s, rng)
        results[s["label"]] = res
        t = tf_gene_table(res)
        t.to_csv(
            TABLE_DIR / f"cross_pcf_pairs_{s['label'].replace('/', '_')}{suffix}.csv",
            index=False,
        )
        tables[s["label"]] = t
        print(
            f"{s['label']}: {len(t)} TF-gene pairs, {int((t['q'] < 0.05).sum())} sig (q<0.05)"
        )

    # summary: significant, same direction, in BOTH sections
    a, b = [tables[s["label"]] for s in good]
    m = a.merge(
        b, on=["tf", "gene", "gene_class", "is_positive_control"], suffixes=("_A", "_B")
    )
    both = m[
        (m["q_A"] < 0.05)
        & (m["q_B"] < 0.05)
        & (np.sign(m["log2_enrich_A"]) == np.sign(m["log2_enrich_B"]))
    ].copy()
    both["log2_enrich_mean"] = (both["log2_enrich_A"] + both["log2_enrich_B"]) / 2
    both = both.sort_values(["tf", "log2_enrich_mean"], ascending=[True, False])
    cols = [
        "tf",
        "gene",
        "gene_class",
        "is_positive_control",
        "log2_enrich_A",
        "log2_enrich_B",
        "log2_enrich_mean",
        "q_A",
        "q_B",
    ]
    both[cols].to_csv(TABLE_DIR / f"cross_pcf_summary{suffix}.csv", index=False)
    print(f"summary: {len(both)} TF-gene co-localisations significant in BOTH sections")

    # heatmaps use mean log2_enrich across samples; star = sig in both, same direction
    ra, rb = results[good[0]["label"]], results[good[1]["label"]]
    tfs = [t for t in C.PRIMARY_TFS if t in ra["idx"] and t in rb["idx"]]

    def cell(res, r_gene, c_gene):
        return res["log2e"][res["idx"][r_gene], res["idx"][c_gene]], res["p"][
            res["idx"][r_gene], res["idx"][c_gene]
        ]

    # BH per sample already in tables; approximate sig-in-both via q<0.05 from tables
    def sig_both(r_gene, c_gene):
        qa = a[(a.tf == r_gene) & (a.gene == c_gene)]["q"]
        qb = b[(b.tf == r_gene) & (b.gene == c_gene)]["q"]
        if len(qa) and len(qb):
            la, _ = cell(ra, r_gene, c_gene)
            lb, _ = cell(rb, r_gene, c_gene)
            return bool(
                qa.iloc[0] < 0.05 and qb.iloc[0] < 0.05 and np.sign(la) == np.sign(lb)
            )
        return False

    present_mk = [mk for mk in MARKERS if mk in ra["idx"] and mk in rb["idx"]]
    groups_present = [
        (nm, [g for g in gs if g in present_mk]) for nm, gs in MARKER_GROUPS
    ]
    groups_present = [(nm, gs) for nm, gs in groups_present if gs]
    mat_mk = np.array(
        [
            [np.mean([cell(ra, tf, mk)[0], cell(rb, tf, mk)[0]]) for mk in present_mk]
            for tf in tfs
        ]
    )
    sig_mk = np.array([[sig_both(tf, mk) for mk in present_mk] for tf in tfs])
    heatmap(
        mat_mk,
        sig_mk,
        tfs,
        present_mk,
        f"TF co-localisation with cell-type markers ({R:g} µm)",
        f"human_tf_cross_pcf_markers{suffix}",
        (0.42 * len(present_mk) + 2.5, 5.5),
        groups=groups_present,
    )

    mat_tf = np.array(
        [
            [
                np.mean([cell(ra, x, y)[0], cell(rb, x, y)[0]]) if x != y else np.nan
                for y in tfs
            ]
            for x in tfs
        ]
    )
    sig_tf = np.array([[sig_both(x, y) if x != y else False for y in tfs] for x in tfs])
    # sharing a cell is mutual, so show a symmetric matrix: average the two directions
    # of the (directional) enrichment, and star a pair if either direction is significant.
    mat_tf = (mat_tf + mat_tf.T) / 2
    sig_tf = sig_tf | sig_tf.T
    heatmap(
        mat_tf,
        sig_tf,
        tfs,
        tfs,
        f"TF-TF co-localisation ({R:g} µm)",
        f"human_tf_cross_pcf_tftf{suffix}",
        (6.5, 5.5),
    )

    # console: top significant colocalisers per TF (both sections)
    print("\nTop co-localisers per TF, significant in both sections (log2 enrich):")
    for tf in tfs:
        sub = both[both.tf == tf].head(6)
        s = "  ".join(f"{r.gene}({r.log2_enrich_mean:+.1f})" for r in sub.itertuples())
        print(f"  {tf:8s} {s if s else '(none)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
