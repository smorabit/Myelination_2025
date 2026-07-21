#!/usr/bin/env python3
"""Compare cell-type proportions inside the lesion vs the whole section.

Joins per-sample Xenium Explorer lesion exports
(data/cells_in_lesions/output-XETG<serial>__*_lesion[_check].csv)
with the celltype-refined labels for the chosen --strategy
(data/celltype_refined/{sample_id}_celltype_refined_{strategy}.csv)
to produce per-sample and cohort-level comparisons of OPC /
Intermediate_Oligo / Mature_Oligo / other proportions inside vs outside
the lesion mask.

Outputs (under data/lesion_celltype_proportions/):
  _per_sample.csv          long form: one row per (sample, celltype)
  _pct_pivot.csv           wide form: pct_in_lesion vs pct_overall side by side
  _enrichment_pivot.csv    wide form: enrichment ratio per (sample, celltype)
  _cohort_summary.csv      cohort-level totals + enrichment
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
import cohort
from spatial_paths import repo_root


import argparse
import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = repo_root()
LESION_DIR = REPO_ROOT / "data" / "cells_in_lesions"
REFINED_DIR = REPO_ROOT / "data" / "celltype_refined"
OUTPUT_DIR = REPO_ROOT / "data" / "lesion_celltype_proportions"
STRATEGY_CHOICES = [
    "stringent",
    "stringent_p4",
    "stringent_p5",
    "brain-ref",
    "strict",
    "legacy",
]
DEFAULT_STRATEGY = "stringent"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--strategy",
        choices=STRATEGY_CHOICES,
        default=DEFAULT_STRATEGY,
        help=(
            f"Cell-type annotation strategy to read from. Default: {DEFAULT_STRATEGY}. "
            "Joins against `data/celltype_refined/{sample_id}_celltype_refined_{strategy}.csv`. "
            "Outputs are also suffixed with the strategy so multiple runs coexist."
        ),
    )
    return parser.parse_args()


# Filename pattern: output-XETG<serial>__{sample_id}__{date}__{time}_lesion[_check].csv
SAMPLE_RE = re.compile(
    r"output-XETG\w+__(\d+__R\d+_\d+)__\d+__\d+_lesion(?:_check)?\.csv$"
)

CELL_TYPES = ["OPC", "Intermediate_Oligo", "Mature_Oligo", "other"]

# Young vs old assignment from the project description.
# Cohort from the samplesheet (shared/cohort.py); no real sample IDs stored here.
YOUNG_SAMPLES = cohort.young_samples()
OLD_SAMPLES = cohort.old_samples()


def parse_lesion_csv(path: Path) -> tuple[str | None, float | None, list[str], bool]:
    """Return (sample_id, lesion_area_um2, cell_ids, has_check_suffix)."""
    m = SAMPLE_RE.search(path.name)
    if not m:
        return None, None, [], False
    sample_id = m.group(1)
    has_check = "_lesion_check.csv" in path.name

    area = None
    with path.open() as fh:
        for line in fh:
            if line.startswith("#Area"):
                area = float(line.split(":", 1)[1].strip())
            elif not line.startswith("#"):
                break

    df = pd.read_csv(path, comment="#")
    if "Cell ID" not in df.columns:
        raise SystemExit(f"`Cell ID` column not found in {path.name}")
    cell_ids = df["Cell ID"].astype(str).tolist()
    return sample_id, area, cell_ids, has_check


def group_for(sample_id: str) -> str:
    return cohort.age_for(sample_id)


def main() -> int:
    args = parse_args()
    strategy = args.strategy
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    lesion_files = sorted(LESION_DIR.glob("output-XETG*__*_lesion*.csv"))
    if not lesion_files:
        raise SystemExit(f"No lesion CSVs in {LESION_DIR}")

    print(f"Lesion dir   : {LESION_DIR}")
    print(f"Refined dir  : {REFINED_DIR} (strategy: {strategy})")
    print(f"Output dir   : {OUTPUT_DIR}")
    print(f"Lesion files : {len(lesion_files)}")
    print()

    rows = []
    for path in lesion_files:
        sample_id, area, lesion_cell_ids, has_check = parse_lesion_csv(path)
        if sample_id is None:
            print(f"  [skip] cannot parse sample from {path.name}", file=sys.stderr)
            continue

        n_lesion = len(lesion_cell_ids)
        if n_lesion == 0:
            print(
                f"  [warn] {sample_id} has 0 lesion cells (empty export from {path.name})"
            )
            continue

        refined_path = REFINED_DIR / f"{sample_id}_celltype_refined_{strategy}.csv"
        if not refined_path.exists():
            print(
                f"  [skip] no refined CSV at {refined_path.relative_to(REPO_ROOT)}",
                file=sys.stderr,
            )
            continue

        refined = pd.read_csv(refined_path)
        n_total = len(refined)
        in_lesion = refined["cell_id"].astype(str).isin(set(lesion_cell_ids))
        n_match = int(in_lesion.sum())
        n_unmatched = n_lesion - n_match

        if n_unmatched > 0:
            print(
                f"  [info] {sample_id}: {n_unmatched} lesion cell IDs not in refined CSV "
                f"(probably segmentation diff between Xenium Explorer and pipeline)"
            )

        for ct in CELL_TYPES:
            mask_ct = refined["celltype_refined"] == ct
            n_in_lesion = int((mask_ct & in_lesion).sum())
            n_overall = int(mask_ct.sum())
            pct_in_lesion = n_in_lesion / n_match * 100 if n_match else 0
            pct_overall = n_overall / n_total * 100 if n_total else 0
            enrichment = pct_in_lesion / pct_overall if pct_overall else float("nan")
            rows.append(
                {
                    "sample_id": sample_id,
                    "group": group_for(sample_id),
                    "celltype": ct,
                    "n_in_lesion": n_in_lesion,
                    "n_overall": n_overall,
                    "n_lesion_matched": n_match,
                    "n_total": n_total,
                    "pct_in_lesion": round(pct_in_lesion, 2),
                    "pct_overall": round(pct_overall, 2),
                    "enrichment_ratio": round(enrichment, 2)
                    if pct_overall
                    else float("nan"),
                    "lesion_area_um2": area,
                    "n_lesion_unmatched": n_unmatched,
                    "lesion_export_flagged": has_check,
                }
            )

    if not rows:
        raise SystemExit(
            "No usable lesion / refined-celltype joins; nothing to report."
        )

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_DIR / f"_per_sample_{strategy}.csv", index=False)

    # Wide pivots for human inspection
    pct_pivot = df.pivot(
        index="sample_id", columns="celltype", values=["pct_in_lesion", "pct_overall"]
    )
    pct_pivot.to_csv(OUTPUT_DIR / f"_pct_pivot_{strategy}.csv")

    enrichment_pivot = df.pivot(
        index="sample_id", columns="celltype", values="enrichment_ratio"
    )
    enrichment_pivot.to_csv(OUTPUT_DIR / f"_enrichment_pivot_{strategy}.csv")

    # Cohort summary (sums weighted by per-sample matched cells)
    cohort_rows = []
    for ct in CELL_TYPES:
        sub = df[df["celltype"] == ct]
        in_lesion = sub["n_in_lesion"].sum()
        overall = sub["n_overall"].sum()
        lesion_total = sub["n_lesion_matched"].sum()
        cohort_total = sub["n_total"].sum()
        pct_l = in_lesion / lesion_total * 100 if lesion_total else 0
        pct_o = overall / cohort_total * 100 if cohort_total else 0
        cohort_rows.append(
            {
                "celltype": ct,
                "n_in_lesion": int(in_lesion),
                "n_overall": int(overall),
                "n_lesion_matched_cohort": int(lesion_total),
                "n_total_cohort": int(cohort_total),
                "pct_in_lesion": round(pct_l, 2),
                "pct_overall": round(pct_o, 2),
                "enrichment_ratio": round(pct_l / pct_o, 2) if pct_o else float("nan"),
            }
        )
    cohort = pd.DataFrame(cohort_rows)
    cohort.to_csv(OUTPUT_DIR / f"_cohort_summary_{strategy}.csv", index=False)

    # Console report
    print()
    print("=== Per-sample percentages ===")
    print()
    for sid in sorted(df["sample_id"].unique()):
        sub = df[df["sample_id"] == sid].set_index("celltype")
        grp = sub["group"].iloc[0]
        n_lm = sub["n_lesion_matched"].iloc[0]
        n_t = sub["n_total"].iloc[0]
        print(f"  {sid} ({grp}): lesion {n_lm}/{n_t} cells")
        for ct in CELL_TYPES:
            r = sub.loc[ct]
            arrow = (
                "↑"
                if r["enrichment_ratio"] > 1.5
                else ("↓" if r["enrichment_ratio"] < 0.66 else "·")
            )
            print(
                f"      {ct:<20s}  in_lesion={r['pct_in_lesion']:5.2f}%  "
                f"overall={r['pct_overall']:5.2f}%  enrichment={r['enrichment_ratio']:.2f}× {arrow}"
            )

    print()
    print("=== Cohort summary ===")
    print(cohort.to_string(index=False))
    print()
    print(
        f"Wrote per-sample, pivot, and cohort tables to {OUTPUT_DIR.relative_to(REPO_ROOT)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
