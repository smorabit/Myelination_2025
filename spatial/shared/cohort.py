"""Cohort / samplesheet loader -- the single source of sample identifiers.

Real sample identifiers are NOT stored in this repository. They live in a samplesheet
CSV that you provide and keep out of git. Copy ``samplesheet.template.csv`` to
``samplesheet.csv`` at your data root (or point ``$SPATIAL_SAMPLESHEET`` at it) and fill
in the real IDs. If no real sheet is found, the bundled template (placeholder IDs) is
used with a warning, so ``--help`` and import checks still work but real runs will not.

Columns:
    species    'mouse' | 'human'
    sample_id  mouse: the per-sample ID embedded in filenames; human: the Xenium output dir name
    label      short human-readable label (e.g. 'Young 1', 'region-1')
    age_group  mouse only: 'Young' | 'Old'
    status     mouse: 'cohort9' | 'outlier' | 'excluded';  human: 'used'
    slide      human only
    flagged    human only: 'true' | 'false' (dim-capture flag)
"""

from __future__ import annotations

import csv
import os
import sys
from functools import lru_cache
from pathlib import Path

from spatial_paths import repo_root

_TEMPLATE = Path(__file__).resolve().parent / "samplesheet.template.csv"


def _samplesheet_path() -> Path:
    env = os.environ.get("SPATIAL_SAMPLESHEET")
    if env:
        return Path(env).expanduser().resolve()
    candidate = repo_root() / "samplesheet.csv"
    if candidate.is_file():
        return candidate
    print(
        f"[cohort] WARNING: no samplesheet.csv found (looked for $SPATIAL_SAMPLESHEET and "
        f"{candidate}); falling back to the placeholder template {_TEMPLATE.name}. "
        f"Copy the template to samplesheet.csv and fill in real IDs for a real run.",
        file=sys.stderr,
    )
    return _TEMPLATE


@lru_cache(maxsize=1)
def _rows() -> tuple[dict, ...]:
    path = _samplesheet_path()
    with open(path, newline="") as fh:
        rows = [
            {k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()}
            for r in csv.DictReader(fh)
        ]
    if not rows:
        raise SystemExit(f"[cohort] samplesheet {path} has no rows")
    return tuple(rows)


def _mouse() -> list[dict]:
    return [r for r in _rows() if r.get("species") == "mouse"]


def young_samples() -> set[str]:
    return {
        r["sample_id"]
        for r in _mouse()
        if r.get("age_group") == "Young" and r.get("status") in ("cohort9", "outlier")
    }


def old_samples() -> set[str]:
    return {
        r["sample_id"]
        for r in _mouse()
        if r.get("age_group") == "Old" and r.get("status") in ("cohort9", "outlier")
    }


def outlier_sample() -> str:
    hits = [r["sample_id"] for r in _mouse() if r.get("status") == "outlier"]
    return hits[0] if hits else ""


def cohort9() -> list[str]:
    return sorted(r["sample_id"] for r in _mouse() if r.get("status") == "cohort9")


def cohort10() -> list[str]:
    return sorted(
        r["sample_id"] for r in _mouse() if r.get("status") in ("cohort9", "outlier")
    )


def default_exclude() -> set[str]:
    return {r["sample_id"] for r in _mouse() if r.get("status") == "excluded"}


def age_for(sample_id: str) -> str:
    if sample_id in young_samples():
        return "Young"
    if sample_id in old_samples():
        return "Old"
    return "Unknown"


def sample_labels() -> dict[str, str]:
    """{sample_id: friendly label} for every mouse cohort9 + outlier sample."""
    return {
        r["sample_id"]: (r.get("label") or r["sample_id"])
        for r in _mouse()
        if r.get("status") in ("cohort9", "outlier")
    }


def human_samples() -> list[dict]:
    """List of {dir, label, slide, flagged} for the human sections."""
    return [
        {
            "dir": r["sample_id"],
            "label": r.get("label", ""),
            "slide": r.get("slide", ""),
            "flagged": str(r.get("flagged", "")).lower() == "true",
        }
        for r in _rows()
        if r.get("species") == "human"
    ]
