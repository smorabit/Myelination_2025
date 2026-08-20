#!/usr/bin/env python3
"""Extract the custom Xenium gene-panel composition from a bundle `gene_panel.json`.

Writes, per panel, a machine-readable JSON summary and a flat TSV of every target
(gene symbol + control type), so the exact panel used for each species is a committed,
citable artifact for the manuscript (Methods / supplementary).

Usage:
    python scripts/extract_panel_composition.py <gene_panel.json> <out_prefix>
e.g.
    python scripts/extract_panel_composition.py data/gene_panel.json \
        Manuscript/xenium_panels/mouse_mBrain_50g
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    src = Path(sys.argv[1])
    out_prefix = Path(sys.argv[2])
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    doc = json.loads(src.read_text())
    payload = doc.get("payload", doc)
    panel = payload.get("panel", {})
    identity = panel.get("identity", {})
    targets = payload.get("targets", [])

    rows = []
    for t in targets:
        ttype = t.get("type", {})
        descriptor = ttype.get("descriptor") if isinstance(ttype, dict) else ttype
        # gene symbol lives under type.data.name (controls carry an id in type.data.id)
        data = ttype.get("data", {}) if isinstance(ttype, dict) else {}
        name = data.get("name") or data.get("id") or ""
        rows.append({"descriptor": descriptor, "name": name})

    by_type = Counter(r["descriptor"] for r in rows)
    genes = sorted(r["name"] for r in rows if r["descriptor"] == "gene")

    summary = {
        "source": str(src),
        "panel_name": identity.get("name"),
        "design_id": identity.get("design_id"),
        "version": identity.get("version"),
        "species": panel.get("species"),
        "tissue": panel.get("tissue"),
        "description": panel.get("description"),
        "num_gene_targets_declared": panel.get("num_gene_targets"),
        "num_targets_total": len(rows),
        "counts_by_type": dict(by_type),
        "num_genes": len(genes),
        "genes": genes,
    }
    (out_prefix.with_suffix(".json")).write_text(json.dumps(summary, indent=2))
    with (out_prefix.with_name(out_prefix.name + "_targets.tsv")).open(
        "w", newline=""
    ) as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["descriptor", "name"])
        for r in sorted(rows, key=lambda x: (x["descriptor"], x["name"])):
            w.writerow([r["descriptor"], r["name"]])

    print(
        f"{summary['panel_name']} ({summary['design_id']}, {summary['species']}): "
        f"{summary['num_genes']} genes, types={summary['counts_by_type']}"
    )
    print(f"  -> {out_prefix.with_suffix('.json')}")
    print(f"  -> {out_prefix.with_name(out_prefix.name + '_targets.tsv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
