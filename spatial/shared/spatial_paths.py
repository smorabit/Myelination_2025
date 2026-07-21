"""Repo-root / data-root resolution for the deposited spatial-Xenium scripts.

The analysis scripts read staged inputs from ``<root>/data/...`` and write
figures to ``<root>/docs/images`` (or the manuscript figure dir). In the original
analysis repo ``<root>`` was the repo root, derived from ``__file__``. In this
deposited layout the scripts live under ``spatial/{mouse,human,...}/`` at varying
depths, so ``__file__``-relative derivation is no longer reliable.

Resolution order:

1. ``$SPATIAL_REPO_ROOT`` if set. **Set this to your data root when reproducing**
   (the directory that contains ``data/`` with the staged GEO objects). This is
   the intended mechanism, matching the repo convention that paths are updated
   per machine.
2. Otherwise, the first ancestor of this file containing a ``.git`` directory
   (i.e. the checked-out repo root). Useful only if you stage ``data/`` inside
   the repo; otherwise data reads will fail loudly with a clear FileNotFound,
   which is the intended behaviour rather than silently wrong paths.
3. Otherwise, the ``spatial/`` directory (last-resort, import-time safe default).
"""

from __future__ import annotations

import os
from pathlib import Path


def repo_root() -> Path:
    """Return the data/figure root for the spatial scripts (see module docstring)."""
    env = os.environ.get("SPATIAL_REPO_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / ".git").exists():
            return parent
    # here == <root>/spatial/shared/spatial_paths.py -> parents[1] == <root>/spatial
    return here.parents[1]
