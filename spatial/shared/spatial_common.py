#!/usr/bin/env python3
"""Shared loaders + cohort constants for the spatial TF analysis (Phase 0+).

Single source of truth for the three-way join every spatial script depends on:

    spatial AnnData (coords + counts)  ×  refined cell types  ×  lesion-mask membership

The per-sample inputs are:
  - data/spatial_anndata/{sample}_spatial_with_annotations.h5ad
        obs: x_centroid, y_centroid, transcript_counts, cell_area, nucleus_area, ...
        adata.raw -> all 50 panel genes (adata.X is the 46-gene processed matrix)
        NOTE: obsm has only X_pca / X_umap — there is NO obsm['spatial'];
              load_sample_cells() constructs it from x_centroid/y_centroid so that
              squidpy (Phase 4) works.
  - data/celltype_refined/{sample}_celltype_refined_{strategy}.csv  (cell_id, celltype_refined, ...)
  - data/cells_in_lesions/output-XETG<serial>__{sample}__*_lesion[_check].csv  (Cell ID, #Area)

Cohort constants and lesion-mask parsing mirror build_pseudobulk.py and
compute_lesion_celltype_composition.py so spatial cells are read identically to the
existing pseudobulk pipeline. (Those two scripts each redefine these small constants
independently; this module is the single definition for the *new* spatial scripts.)
"""

from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from shapely.geometry import MultiPoint, Polygon
from shapely.ops import unary_union

import cohort
from spatial_paths import repo_root

REPO_ROOT = repo_root()
SPATIAL_DIR = REPO_ROOT / "data" / "spatial_anndata"
LESION_DIR = REPO_ROOT / "data" / "cells_in_lesions"
REFINED_DIR = REPO_ROOT / "data" / "celltype_refined"
OUTPUT_DIR = REPO_ROOT / "data" / "spatial"

DEFAULT_STRATEGY = "stringent_p5"

OLIGO_TYPES = ["OPC", "Intermediate_Oligo", "Mature_Oligo"]

# 9 primary transcription factors (the spatial hypothesis genes).
PRIMARY_TFS = [
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

# Cohort identifiers come from the samplesheet (shared/cohort.py, gitignored
# samplesheet.csv); no real sample IDs are stored in this repository. Young 2 was the
# `_lesion_check`-flagged export; Old 2 is the canonical outlier dropped from COHORT_9.
YOUNG_SAMPLES = cohort.young_samples()
OLD_SAMPLES = cohort.old_samples()
OUTLIER_SAMPLE = (
    cohort.outlier_sample()
)  # Old 2 — dropped in the canonical 9-sample cohort
COHORT_9 = cohort.cohort9()  # 5 Young + 4 Old (Old 2 excluded)
COHORT_10 = cohort.cohort10()

# output-XETG<serial>__{sample}__{date}__{time}_lesion[_check].csv
LESION_FILE_RE = re.compile(
    r"output-XETG\w+__(\d+__R\d+_\d+)__\d+__\d+_lesion(?:_check)?\.csv$"
)


def age_for(sample_id: str) -> str:
    return cohort.age_for(sample_id)


def chip_for(sample_id: str) -> str:
    return sample_id.split("__")[0]


def alpha_shape(
    points: np.ndarray,
    alpha_radius: float | None = None,
    percentile: float = 99.0,
):
    """Concave hull (alpha shape) of 2D points via scipy Delaunay + shapely.

    Dependency-light alternative to the `alphashape` package (scipy + shapely only).

    Algorithm:
      1. Delaunay-triangulate the points.
      2. For each triangle compute its circumradius R = (a*b*c) / (4*area), where
         a/b/c are the edge lengths. Degenerate (zero-area / collinear) triangles
         get R = +inf and are dropped.
      3. Keep triangles whose circumradius R <= `alpha_radius`; shapely
         `unary_union` of those triangle polygons is the concave hull.

    `alpha_radius` semantics: larger radius -> more triangles kept -> hull closer to
    the convex hull; smaller radius -> tighter, more concave hull. When
    `alpha_radius is None` it is chosen data-adaptively as the `percentile`-th (default
    99th) percentile of the finite triangle circumradii. The 99th percentile was chosen
    empirically over this cohort: lower percentiles (90/95) fragment a single contiguous
    lesion annotation into many disconnected components (a Delaunay sliver-triangle
    artifact, not real fragmentation) and undercount area (ratio vs export #Area ~0.4-0.6),
    whereas the 99th percentile yields a single connected component for every sample with
    area ratios ~0.6-0.8 of the hand-drawn export #Area. (A ratio < 1 is expected: the
    cell-derived hull cannot extend past the outermost cells and excludes cell-sparse
    interior gaps that the hand-drawn polygon includes.)

    Edge cases (never raises):
      - fewer than 4 points, or all points collinear (Delaunay fails / zero finite
        circumradii) -> fall back to the convex hull of `MultiPoint(points)`.

    Returns:
      (geometry, info) where geometry is a shapely (Multi)Polygon (or a lower-dim
      geometry only in pathological <3-point inputs) and info is a dict:
        {"alpha_radius": float|None, "fallback_convex": bool, "n_triangles_kept": int,
         "n_triangles_total": int, "reason": str}.
    """
    pts = np.asarray(points, dtype=float)
    info: dict = {
        "alpha_radius": alpha_radius,
        "fallback_convex": False,
        "n_triangles_kept": 0,
        "n_triangles_total": 0,
        "reason": "",
    }

    def _convex_fallback(reason: str):
        info["fallback_convex"] = True
        info["reason"] = reason
        return MultiPoint([tuple(p) for p in pts]).convex_hull, info

    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"alpha_shape expects (N, 2) points, got shape {pts.shape}")
    if pts.shape[0] < 4:
        return _convex_fallback("fewer than 4 points")

    try:
        tri = Delaunay(pts)
    except Exception as exc:  # collinear / degenerate input
        return _convex_fallback(f"Delaunay failed ({exc})")

    simplices = tri.simplices
    info["n_triangles_total"] = int(len(simplices))
    if len(simplices) == 0:
        return _convex_fallback("no triangles (degenerate)")

    ia, ib, ic = simplices[:, 0], simplices[:, 1], simplices[:, 2]
    pa, pb, pc = pts[ia], pts[ib], pts[ic]
    a = np.linalg.norm(pb - pc, axis=1)
    b = np.linalg.norm(pa - pc, axis=1)
    c = np.linalg.norm(pa - pb, axis=1)
    # signed area via cross product; abs for triangle area
    area = 0.5 * np.abs(
        (pb[:, 0] - pa[:, 0]) * (pc[:, 1] - pa[:, 1])
        - (pc[:, 0] - pa[:, 0]) * (pb[:, 1] - pa[:, 1])
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        circum_r = (a * b * c) / (4.0 * area)
    circum_r[~np.isfinite(circum_r)] = np.inf  # guard zero-area triangles

    finite = circum_r[np.isfinite(circum_r)]
    if finite.size == 0:
        return _convex_fallback("all triangles degenerate (collinear points)")

    radius = (
        float(np.percentile(finite, percentile))
        if alpha_radius is None
        else float(alpha_radius)
    )
    info["alpha_radius"] = radius

    keep = circum_r <= radius
    info["n_triangles_kept"] = int(keep.sum())
    if not keep.any():
        return _convex_fallback("no triangles below alpha_radius")

    triangles = [Polygon(pts[simplices[k]]) for k in np.nonzero(keep)[0]]
    hull = unary_union(triangles)
    if hull.is_empty or hull.area == 0:
        return _convex_fallback("union empty/zero-area")
    return hull, info


def find_lesion_csv(sample_id: str) -> Path | None:
    """Return the lesion-mask CSV for a sample (there is exactly one per sample)."""
    matches = [
        p
        for p in sorted(LESION_DIR.glob("output-XETG*__*_lesion*.csv"))
        if (m := LESION_FILE_RE.search(p.name)) and m.group(1) == sample_id
    ]
    if not matches:
        return None
    return matches[0]


def parse_lesion_csv(path: Path) -> tuple[set[str], float | None, bool]:
    """Return (lesion_cell_ids, lesion_area_um2, is_check_flagged).

    Mirrors compute_lesion_celltype_composition.parse_lesion_csv, which (unlike
    build_pseudobulk.load_lesion_cell_ids) also reads the `#Area` header line — needed
    for the Phase 0.3 morphometrics cross-check.
    """
    is_check = "_lesion_check.csv" in path.name
    area: float | None = None
    with path.open() as fh:
        for line in fh:
            if line.startswith("#Area"):
                area = float(line.split(":", 1)[1].strip())
            elif not line.startswith("#"):
                break
    df = pd.read_csv(path, comment="#")
    if "Cell ID" not in df.columns:
        raise ValueError(f"`Cell ID` column not found in {path.name}")
    cell_ids = set(df["Cell ID"].astype(str))
    return cell_ids, area, is_check


def load_sample_cells(
    sample_id: str,
    strategy: str = DEFAULT_STRATEGY,
    oligo_only: bool = False,
    use_raw: bool = True,
) -> ad.AnnData:
    """Load one sample's cells with the three-way join applied.

    Returns an AnnData whose obs carries the joined annotations and whose obsm['spatial']
    is set from the Xenium centroids. By default keeps ALL annotated cells (the lesion
    mask spans every cell type, so geometry must use all lesion cells, not just oligos);
    pass oligo_only=True to restrict to OPC/Intermediate_Oligo/Mature_Oligo.

    Added obs columns:
      celltype_refined, in_lesion (bool), sample_id, age_group, chip
    Added uns:
      lesion_area_um2, lesion_export_flagged, n_lesion_unmatched, n_lesion_total

    use_raw=True returns the 50-gene panel matrix (adata.raw promoted to X); use_raw=False
    keeps the 46-gene processed X. Either way obs/obsm annotations are identical.
    """
    spatial_path = SPATIAL_DIR / f"{sample_id}_spatial_with_annotations.h5ad"
    refined_path = REFINED_DIR / f"{sample_id}_celltype_refined_{strategy}.csv"
    if not spatial_path.exists():
        raise FileNotFoundError(f"spatial AnnData missing: {spatial_path}")
    if not refined_path.exists():
        raise FileNotFoundError(f"refined CSV missing: {refined_path}")

    adata = ad.read_h5ad(spatial_path)
    if use_raw:
        if adata.raw is None:
            raise ValueError(
                f"{sample_id}: adata.raw is None — cannot get 50-gene panel"
            )
        adata = adata.raw.to_adata()  # obs is carried over unchanged

    obs_names = adata.obs_names.astype(str)

    # --- join refined cell types (inner: only annotated cells kept) ---
    refined = pd.read_csv(refined_path)
    refined.index = refined["cell_id"].astype(str)
    ct = refined["celltype_refined"].reindex(obs_names)
    n_unannotated = int(ct.isna().sum())
    adata.obs["celltype_refined"] = ct.to_numpy()

    # --- join lesion membership ---
    lesion_path = find_lesion_csv(sample_id)
    if lesion_path is None:
        raise FileNotFoundError(f"no lesion mask CSV for {sample_id}")
    lesion_ids, area, is_check = parse_lesion_csv(lesion_path)
    in_lesion = np.asarray(obs_names.isin(lesion_ids))
    adata.obs["in_lesion"] = in_lesion
    n_unmatched = len(lesion_ids) - int(
        in_lesion.sum()
    )  # mask IDs absent from pipeline

    # --- cohort metadata ---
    adata.obs["sample_id"] = sample_id
    adata.obs["age_group"] = age_for(sample_id)
    adata.obs["chip"] = chip_for(sample_id)

    # --- coordinates -> obsm['spatial'] (does not exist in the h5ad) ---
    for col in ("x_centroid", "y_centroid"):
        if col not in adata.obs.columns:
            raise ValueError(f"{sample_id}: obs missing '{col}'")
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy(
        dtype=float
    )

    adata.uns["lesion_area_um2"] = area
    adata.uns["lesion_export_flagged"] = is_check
    adata.uns["n_lesion_total"] = len(lesion_ids)
    adata.uns["n_lesion_unmatched"] = n_unmatched
    adata.uns["n_cells_unannotated"] = n_unannotated

    # drop cells with no refined label (not in refined CSV), then optionally oligo-only
    adata = adata[~pd.isna(adata.obs["celltype_refined"])].copy()
    if oligo_only:
        adata = adata[adata.obs["celltype_refined"].isin(OLIGO_TYPES)].copy()

    return adata


def iter_cohort(
    cohort: list[str] | None = None,
    strategy: str = DEFAULT_STRATEGY,
    oligo_only: bool = False,
    use_raw: bool = True,
):
    """Yield (sample_id, AnnData) for each sample in the cohort (default COHORT_9)."""
    for sample_id in cohort if cohort is not None else COHORT_9:
        yield (
            sample_id,
            load_sample_cells(
                sample_id, strategy=strategy, oligo_only=oligo_only, use_raw=use_raw
            ),
        )


if __name__ == "__main__":
    # Smoke test: load each cohort sample and report the join health.
    import sys

    cohort = COHORT_10 if "--all" in sys.argv else COHORT_9
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Strategy: {DEFAULT_STRATEGY}   cohort: {len(cohort)} samples\n")
    print(
        f"{'sample':<22} {'age':<6} {'cells':>7} {'in_les':>7} "
        f"{'OPC':>5} {'Int':>5} {'Mat':>5} {'unmatch':>8} {'flag':>5} {'area_um2':>12}"
    )
    for sid in cohort:
        a = load_sample_cells(sid)
        les = a[a.obs["in_lesion"]]
        cts = les.obs["celltype_refined"].value_counts()
        print(
            f"{sid:<22} {age_for(sid):<6} {a.n_obs:>7} {int(a.obs['in_lesion'].sum()):>7} "
            f"{cts.get('OPC', 0):>5} {cts.get('Intermediate_Oligo', 0):>5} "
            f"{cts.get('Mature_Oligo', 0):>5} {a.uns['n_lesion_unmatched']:>8} "
            f"{'chk' if a.uns['lesion_export_flagged'] else '-':>5} "
            f"{a.uns['lesion_area_um2']:>12.0f}"
        )
