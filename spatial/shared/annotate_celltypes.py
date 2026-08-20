#!/usr/bin/env python3
"""Refine oligodendrocyte-lineage cell-type labels using Xenium marker rules.

Applies the marker rules from the project description to each per-sample
spatial AnnData and produces a refined per-cell category in
{OPC, Intermediate_Oligo, Mature_Oligo, other}.

Marker rules (from docs/README.md):

    OPC                = Olig2+ Sox10+ AND (Pdgfra+ OR Ptprz1+ OR Pcdh15+)
    Intermediate_Oligo = Olig2+ Sox10+ AND Enpp6+
    Mature_Oligo       = Olig2+ Sox10+ AND (Mbp+ OR Opalin+ OR Mog+)
    other              = none of the above

Four strategies are supported via `--strategy`. Default is **`stringent`**.
Primary analysis runs **both** stringent and brain-ref; the other two are
sensitivity-check / reference variants.

    stringent (default, primary) — the reference cell-typing criteria taken literally:
        OPC                = Olig2+ Sox10+ AND (Pdgfra+ OR Ptprz1+ OR Pcdh15+)
        Intermediate_Oligo = Olig2+ Sox10+ AND Enpp6+
        Mature_Oligo       = Olig2+ Sox10+ AND (Mbp+ OR Opalin+ OR Mog+)
        other              = anything else
        Resolution: Intermediate > OPC > Mature (Enpp6 gates the transitional
        call; among Enpp6-negative cells, OPC-marker positivity beats
        mature-marker positivity because Mbp probe transcripts diffuse beyond
        their cell of origin in Xenium and inflate the mature signal).

        This is the docx specification verbatim — "at least one of" the
        cell-type-specific markers. The recommended primary because the call
        is defensible from the docx criteria alone, with no reliance on the
        snRNAseq reference. The Intermediate>OPC>Mature precedence is a
        choice we make on top of the docx (which doesn't specify resolution),
        designed to be Mbp-diffusion-tolerant.

    brain-ref (recall complement, run alongside primary):
        OPC                = brain-ref celltype_l1 == 'OPC'
        Intermediate_Oligo = brain-ref OLG AND Enpp6+
        Mature_Oligo       = brain-ref OLG AND (Mbp+ OR Opalin+ OR Mog+) AND Enpp6-
        other              = anything else

        Trusts the brain-reference probabilistic classifier (trained on
        snRNAseq with many genes) for the high-level OPC vs OLG split, then
        applies the docx markers within OLG. Higher recall than stringent
        (~3,229 OPC / 4,256 Inter / 14,715 Mature cohort-wide). Robust to
        the Mbp diffusion problem (Mbp is detected in 93 % of brain-ref
        OPCs in sample R4 -- a Xenium issue with myelin transcripts spreading
        beyond their cell of origin), but a TF-DE signal that surfaces only
        here and not in stringent should be checked for that artefact.

    strict (sensitivity check):
        OPC                = Olig2+ Sox10+ AND OPC-marker+ AND mature-marker-
        Mature_Oligo       = Olig2+ Sox10+ AND mature-marker+ AND OPC-marker-
        Intermediate_Oligo = Olig2+ Sox10+ AND Enpp6+ (overrides OPC/Mature if both)
        other              = anything else

        Pure marker-driven classification with mutually-exclusive OPC and
        Mature pools (no requirement that ALL markers be detected, just that
        opposing-fate markers are absent). Yields tiny OPC pools (~47 cells
        across 10 samples) because of the Mbp diffusion problem.

    legacy (sensitivity check):
        Same per-rule definitions as 'strict' but Mature > Intermediate >
        OPC precedence (a cell with both OPC and mature markers is called
        Mature). Even more permissive; useful only for back-compat / range.

Per-rule boolean flags (`is_opc_by_markers`, `is_intermediate_by_markers`,
`is_mature_by_markers`) are kept in the output so anyone can re-resolve
with a different rule without re-running.

Counts come from `adata.raw` so the 4 panel genes filtered out of `.X`
(Cpa1, Plag1, Selenoh, Spink1) are still accessible -- although none of
the 9 cell-typing markers happen to be among the dropped 4, this keeps
the script's data source consistent with the downstream DE pipeline.

Inputs:  data/spatial_anndata/{sample_id}_spatial_with_annotations.h5ad
Outputs: data/celltype_refined/{sample_id}_celltype_refined_{strategy}.csv
         data/celltype_refined/_summary_{strategy}.csv
         data/celltype_refined/_overlap_with_brain_reference_{strategy}.csv
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from spatial_paths import repo_root

REPO_ROOT = repo_root()
DEFAULT_INPUT_DIR = REPO_ROOT / "data" / "spatial_anndata"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "celltype_refined"

CORE_MARKERS = ("Olig2", "Sox10")
OPC_EXTRAS = ("Pdgfra", "Ptprz1", "Pcdh15")
INTER_MARKERS = ("Enpp6",)
MATURE_EXTRAS = ("Mbp", "Opalin", "Mog")

ALL_MARKERS = sorted(set(CORE_MARKERS + OPC_EXTRAS + INTER_MARKERS + MATURE_EXTRAS))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directory of input .h5ad files (default: {DEFAULT_INPUT_DIR.relative_to(REPO_ROOT)})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for output CSVs (default: {DEFAULT_OUTPUT_DIR.relative_to(REPO_ROOT)})",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Marker counts > threshold are 'expressing' (default: 0, i.e. binary).",
    )
    parser.add_argument(
        "--use-x",
        action="store_true",
        help="Use adata.X instead of adata.raw. Default uses raw for full-panel coverage.",
    )
    parser.add_argument(
        "--strategy",
        choices=[
            "stringent",
            "stringent_p4",
            "stringent_p5",
            "brain-ref",
            "strict",
            "legacy",
        ],
        default="stringent",
        help=(
            "Classification strategy. "
            "'stringent' (P1): the reference cell-typing criteria literally — "
            "Olig2+Sox10+ + AT LEAST ONE of the cell-type-specific markers. "
            "Resolution: Intermediate > OPC > Mature (Mbp-diffusion-tolerant precedence). "
            "'stringent_p4' (P4 — diffusion-corrected variant): same as stringent, but the "
            "Mature OR-set drops Mbp (Mature = Olig2+ Sox10+ AND (Mog+ OR Opalin+)); "
            "precedence Intermediate > OPC > Mature unchanged. "
            "'stringent_p5' (P5 — diffusion-corrected variant with maturation-priority precedence): "
            "same Mature OR-set as P4 (Mog+ OR Opalin+, no Mbp), but precedence is "
            "Intermediate > Mature > OPC: a cell with both OPC criteria and verifiable mature "
            "markers (Mog or Opalin, not Mbp) is called Mature, on the biology that genuine "
            "myelin gene expression marks commitment to the mature programme. "
            "Targets the diffusion problem at its source (Mbp leaks ~5x more than Mog and "
            "~13x more than Opalin in our data — see methodology audit in lesion_analysis_plan.md). "
            "'brain-ref' (recall complement): trust brain-reference celltype_l1 for OPC vs OLG, "
            "apply Enpp6 / Mbp / Opalin / Mog markers within OLG. "
            "'strict' (sensitivity): same docx 'at least one' criteria but with mutual exclusion "
            "(cells with both OPC and mature markers go to 'other'). "
            "'legacy' (sensitivity): docx criteria with Mature > Intermediate > OPC precedence."
        ),
    )
    return parser.parse_args()


def marker_matrix(
    adata: ad.AnnData, use_raw: bool
) -> tuple[np.ndarray, dict[str, int]]:
    """Return (n_cells x len(ALL_MARKERS)) dense float matrix and gene->col index."""
    src = adata.raw.to_adata() if (use_raw and adata.raw is not None) else adata
    missing = [g for g in ALL_MARKERS if g not in src.var_names]
    if missing:
        raise SystemExit(f"Markers missing from matrix: {missing}")
    sub = src[:, list(ALL_MARKERS)]
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    return np.asarray(X), {g: i for i, g in enumerate(ALL_MARKERS)}


def classify(
    counts: np.ndarray,
    idx: dict[str, int],
    threshold: float,
    strategy: str,
    brain_ref_labels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    expr = counts > threshold
    core = expr[:, idx["Olig2"]] & expr[:, idx["Sox10"]]
    opc_extras = (
        expr[:, idx["Pdgfra"]] | expr[:, idx["Ptprz1"]] | expr[:, idx["Pcdh15"]]
    )
    inter_extras = expr[:, idx["Enpp6"]]
    mature_extras = expr[:, idx["Mbp"]] | expr[:, idx["Opalin"]] | expr[:, idx["Mog"]]

    # Per-rule boolean flags (kept in output for downstream re-classification).
    is_opc = core & opc_extras
    is_inter = core & inter_extras
    is_mature = core & mature_extras

    label = np.full(counts.shape[0], "other", dtype=object)

    if strategy == "brain-ref":
        # Trust brain-ref celltype_l1 for OPC and OLG, then use Enpp6 / mature markers
        # to split OLG into intermediate vs mature. Avoids the Mbp diffusion problem
        # where mature-marker contamination in real OPCs would otherwise mis-classify them.
        is_opc_ref = brain_ref_labels == "OPC"
        is_olg_ref = brain_ref_labels == "OLG"
        label[is_opc_ref] = "OPC"
        # Mature within OLG: any mature marker AND no Enpp6.
        label[is_olg_ref & mature_extras & ~inter_extras] = "Mature_Oligo"
        # Intermediate within OLG: Enpp6+ (overrides Mature if both — Enpp6 is the
        # transitional marker so its presence defines intermediate stage).
        label[is_olg_ref & inter_extras] = "Intermediate_Oligo"
    elif strategy == "stringent":
        # the reference cell-typing criteria, taken literally:
        #   OPC                = Olig2+ Sox10+ AND (Pdgfra+ OR Ptprz1+ OR Pcdh15+)
        #   Intermediate_Oligo = Olig2+ Sox10+ AND Enpp6+
        #   Mature_Oligo       = Olig2+ Sox10+ AND (Mbp+ OR Opalin+ OR Mog+)
        #
        # i.e. is_opc / is_inter / is_mature as already computed above match
        # the docx rules exactly (the marker UNION, not intersection).
        #
        # Resolution for cells satisfying multiple rules: Intermediate > OPC > Mature.
        # Rationale:
        #   - Enpp6 is the gating transitional marker, so its presence defines
        #     the intermediate state (Intermediate beats both OPC and Mature).
        #   - Among Enpp6-negative cells, positive OPC-marker evidence beats
        #     positive mature-marker evidence — Mbp probe transcripts diffuse
        #     beyond their cell of origin in Xenium (~89 % of non-oligo cells
        #     in our data have Mbp+ vs ~18 % Mog+ and ~7 % Opalin+ — see the
        #     methodology audit in docs/analysis/lesion_analysis_plan.md), so
        #     a Mature-precedence rule would mis-classify true OPCs that picked
        #     up neighbouring Mbp signal. Pdgfra / Ptprz1 / Pcdh15 don't have a
        #     known diffusion artefact, so OPC-marker positivity is the more
        #     reliable evidence when both marker sets are present.
        #
        # Implemented as reverse-precedence assignment (last write wins):
        label[is_mature] = "Mature_Oligo"
        label[is_opc] = "OPC"  # OPC overrides Mature (Mbp-diffusion-tolerant)
        label[is_inter] = "Intermediate_Oligo"  # Enpp6+ overrides everything
    elif strategy == "stringent_p4":
        # Diffusion-corrected variant (P4) — the methodology audit found that
        # ~89 % of cells outside the oligo lineage in spinal cord have Mbp+,
        # vs ~18 % Mog+ and ~7 % Opalin+. Mbp leaks ~5 x more than Mog and
        # ~13 x more than Opalin (Ainger et al. 1997; Müller et al. 2013;
        # Yang et al. 2025 MisTIC — citations in lesion_analysis_plan.md).
        # P4 targets the diffusion problem at its source: drop Mbp from the
        # Mature OR-set so the Mature criterion is:
        #   Mature_Oligo = Olig2+ Sox10+ AND (Mog+ OR Opalin+)
        # OPC and Intermediate criteria are unchanged. Precedence is
        # Intermediate > OPC > Mature (same as stringent). Cells that under
        # stringent were Mature only because of Mbp+ alone (no Mog+/Opalin+)
        # become "other" under P4 — exactly the population the diagnostics
        # identified as most likely Mbp-diffusion artefacts.
        is_mature_p4 = core & (expr[:, idx["Mog"]] | expr[:, idx["Opalin"]])
        label[is_mature_p4] = "Mature_Oligo"
        label[is_opc] = "OPC"
        label[is_inter] = "Intermediate_Oligo"
        # Update the per-rule flag so downstream consumers can re-derive:
        is_mature = is_mature_p4
    elif strategy == "stringent_p5":
        # P5 — same Mature OR-set as P4 (Mog+ OR Opalin+, Mbp dropped) but
        # precedence Intermediate > Mature > OPC: a cell satisfying both OPC
        # criteria and a verifiable mature marker (Mog or Opalin) is called
        # Mature, on the biology that genuine myelin gene expression marks
        # commitment to the mature programme. Differs from P4 only in
        # precedence: cells with both OPC and Mature criteria go to Mature
        # under P5 (vs OPC under P4). Intermediate (Enpp6+) still wins
        # highest priority.
        is_mature_p5 = core & (expr[:, idx["Mog"]] | expr[:, idx["Opalin"]])
        label[is_opc] = "OPC"
        label[is_mature_p5] = "Mature_Oligo"  # Mature beats OPC under P5
        label[is_inter] = "Intermediate_Oligo"  # Enpp6+ beats everything
        is_mature = is_mature_p5
    elif strategy == "strict":
        # Strict mutual exclusion: cells with both OPC and mature markers (no Enpp6)
        # go to 'other'; Enpp6 wins over the strict definitions if there's overlap.
        strict_opc = is_opc & ~mature_extras
        strict_mature = is_mature & ~opc_extras
        label[strict_opc] = "OPC"
        label[strict_mature] = "Mature_Oligo"
        label[is_inter] = "Intermediate_Oligo"
    elif strategy == "legacy":
        # Mature > Intermediate > OPC precedence via assignment order.
        label[is_opc] = "OPC"
        label[is_inter] = "Intermediate_Oligo"
        label[is_mature] = "Mature_Oligo"
    else:
        raise SystemExit(f"Unknown strategy: {strategy}")

    return is_opc, is_inter, is_mature, label


def process(
    path: Path, output_dir: Path, threshold: float, use_raw: bool, strategy: str
) -> dict:
    sample_id = path.name.replace("_spatial_with_annotations.h5ad", "")
    adata = ad.read_h5ad(path)
    counts, idx = marker_matrix(adata, use_raw)
    brain_ref_labels = adata.obs["celltype_l1"].astype(str).to_numpy()
    is_opc, is_inter, is_mature, label = classify(
        counts, idx, threshold, strategy, brain_ref_labels
    )

    df = pd.DataFrame(
        {
            "cell_id": np.asarray(adata.obs_names),
            "celltype_l1_brain_ref": adata.obs["celltype_l1"].astype(str).to_numpy(),
            "celltype_l1_confidence": adata.obs["celltype_l1_confidence"].to_numpy(),
            "celltype_refined": label,
            "is_opc_by_markers": is_opc,
            "is_intermediate_by_markers": is_inter,
            "is_mature_by_markers": is_mature,
            "x_centroid": adata.obs["x_centroid"].to_numpy(),
            "y_centroid": adata.obs["y_centroid"].to_numpy(),
        }
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{sample_id}_celltype_refined_{strategy}.csv"
    df.to_csv(out_path, index=False)

    counts_summary = Counter(label)
    return {
        "sample_id": sample_id,
        "n_cells": df.shape[0],
        "n_OPC": counts_summary.get("OPC", 0),
        "n_Intermediate_Oligo": counts_summary.get("Intermediate_Oligo", 0),
        "n_Mature_Oligo": counts_summary.get("Mature_Oligo", 0),
        "n_other": counts_summary.get("other", 0),
        "_df": df,
        "_path": out_path,
    }


def main() -> int:
    args = parse_args()
    use_raw = not args.use_x
    if not args.input_dir.is_dir():
        raise SystemExit(f"Input dir not found: {args.input_dir}")
    paths = sorted(args.input_dir.glob("*_spatial_with_annotations.h5ad"))
    if not paths:
        raise SystemExit(f"No .h5ad files in {args.input_dir}")

    print(f"Input dir  : {args.input_dir}")
    print(f"Output dir : {args.output_dir}")
    print(f"Threshold  : counts > {args.threshold}")
    print(f"Source     : adata.{'raw' if use_raw else 'X'}")
    print(f"Strategy   : {args.strategy}")
    print(f"Samples    : {len(paths)}")
    print()

    results = []
    for path in paths:
        r = process(path, args.output_dir, args.threshold, use_raw, args.strategy)
        n = r["n_cells"]
        pct_opc = r["n_OPC"] / n * 100 if n else 0.0
        pct_int = r["n_Intermediate_Oligo"] / n * 100 if n else 0.0
        pct_mat = r["n_Mature_Oligo"] / n * 100 if n else 0.0
        pct_oth = r["n_other"] / n * 100 if n else 0.0
        print(
            f"  {r['sample_id']}: n={n:>6}  "
            f"OPC={r['n_OPC']:>5} ({pct_opc:4.1f}%)  "
            f"Inter={r['n_Intermediate_Oligo']:>5} ({pct_int:4.1f}%)  "
            f"Mature={r['n_Mature_Oligo']:>5} ({pct_mat:4.1f}%)  "
            f"other={r['n_other']:>5} ({pct_oth:4.1f}%)"
        )
        results.append(r)

    summary = pd.DataFrame(
        [
            {
                k: r[k]
                for k in (
                    "sample_id",
                    "n_cells",
                    "n_OPC",
                    "n_Intermediate_Oligo",
                    "n_Mature_Oligo",
                    "n_other",
                )
            }
            for r in results
        ]
    )
    for col in ("OPC", "Intermediate_Oligo", "Mature_Oligo", "other"):
        summary[f"pct_{col}"] = (summary[f"n_{col}"] / summary["n_cells"] * 100).round(
            2
        )
    summary_path = args.output_dir / f"_summary_{args.strategy}.csv"
    summary.to_csv(summary_path, index=False)

    combined = pd.concat([r["_df"] for r in results], ignore_index=True)
    overlap = pd.crosstab(
        combined["celltype_refined"], combined["celltype_l1_brain_ref"]
    )
    overlap_path = (
        args.output_dir / f"_overlap_with_brain_reference_{args.strategy}.csv"
    )
    overlap.to_csv(overlap_path)

    print()
    print("=== Cohort totals ===")
    totals = summary[
        ["n_cells", "n_OPC", "n_Intermediate_Oligo", "n_Mature_Oligo", "n_other"]
    ].sum()
    print(totals.to_string())
    print()
    print(
        f"=== Refined label x brain-reference celltype_l1 (cohort, {combined.shape[0]} cells) ==="
    )
    print(overlap.to_string())
    print()
    print(
        f"Wrote per-sample CSVs and 2 summary tables to {args.output_dir.relative_to(REPO_ROOT)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
