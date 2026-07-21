#!/usr/bin/env python3
"""Plot and test oligodendrocyte-lineage proportions in the lesion (Young vs Old).

For each strategy in --strategies (default: stringent + brain-ref):
  1. Loads per-sample lesion counts from data/lesion_celltype_proportions/_per_sample_{strategy}.csv
  2. Computes WITHIN-LINEAGE proportions: OPC / (OPC + Intermediate + Mature) etc.
     (so the three columns sum to 1 for each sample). This isolates the lineage
     composition question from the overall lesion-vs-section enrichment.
  3. Runs a two-sided Mann-Whitney U (Wilcoxon rank-sum) test for Young vs Old
     per cell type. Non-parametric is appropriate for n = 5 vs 5.
  4. Produces a strip + box plot (one panel per cell type).

Outputs (under data/lesion_celltype_proportions/):
  _lineage_proportions_{strategy}.csv   per-sample within-lineage proportions
  _young_vs_old_test_{strategy}.csv     Mann-Whitney U + medians per cell type
  _composition_permanova_{strategy}.csv global CLR + PERMANOVA test (one p, closure-aware)
  figures/lineage_young_vs_old_{strategy}.png

The per-cell-type Mann-Whitney tests treat the three within-lineage proportions as
separate outcomes, but they sum to 1 (compositional closure), so they are not
independent. The PERMANOVA is the closure-aware global companion: it centre-log-ratio
(CLR) transforms each sample's OPC:Intermediate:Mature composition, then tests whether
the Young and Old groups differ with a distance-based pseudo-F (Anderson 2001) on the
Aitchison distances, significance from an exact label permutation. It answers "does the
lineage composition differ with age?" with a single p that respects the closure.
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
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import plot_style

REPO_ROOT = repo_root()
PROP_DIR = REPO_ROOT / "data" / "lesion_celltype_proportions"
FIG_DIR = PROP_DIR / "figures"

OLIGO_TYPES = ["OPC", "Intermediate_Oligo", "Mature_Oligo"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--strategies",
        nargs="+",
        default=["stringent", "brain-ref"],
        help="Strategies to process (default: stringent brain-ref).",
    )
    return parser.parse_args()


def compute_lineage_proportions(strategy: str) -> pd.DataFrame:
    df = pd.read_csv(PROP_DIR / f"_per_sample_{strategy}.csv")
    sub = df[df["celltype"].isin(OLIGO_TYPES)].copy()
    pivot = sub.pivot(index="sample_id", columns="celltype", values="n_in_lesion")
    pivot = pivot[OLIGO_TYPES]  # ensure column order
    lineage_total = pivot.sum(axis=1)
    proportions = pivot.div(lineage_total.replace(0, np.nan), axis=0)
    group_lookup = sub.drop_duplicates("sample_id").set_index("sample_id")["group"]
    proportions["group"] = group_lookup
    proportions["lineage_total_in_lesion"] = lineage_total
    return proportions.reset_index()


def test_young_vs_old(prop_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ct in OLIGO_TYPES:
        young = prop_df.loc[prop_df["group"] == "Young", ct].dropna().values
        old = prop_df.loc[prop_df["group"] == "Old", ct].dropna().values
        u, p = mannwhitneyu(young, old, alternative="two-sided")
        rows.append(
            {
                "celltype": ct,
                "n_young": len(young),
                "n_old": len(old),
                "median_young": float(np.median(young)) if len(young) else float("nan"),
                "median_old": float(np.median(old)) if len(old) else float("nan"),
                "median_diff_old_minus_young": (
                    float(np.median(old) - np.median(young))
                    if len(young) and len(old)
                    else float("nan")
                ),
                "U_statistic": float(u),
                "p_value_two_sided": float(p),
            }
        )
    return pd.DataFrame(rows)


def _clr(counts: np.ndarray, pseudo: float = 0.5) -> np.ndarray:
    """Centre-log-ratio transform of each row (composition). Scale-invariant, so
    counts and proportions give identical coordinates; a small pseudocount keeps
    log() finite if any part is zero."""
    logx = np.log(counts + pseudo)
    return logx - logx.mean(axis=1, keepdims=True)


def permanova_young_vs_old(strategy: str) -> pd.DataFrame:
    """Closure-aware global test of the OPC:Intermediate:Mature composition, Young
    vs Old. Distance-based PERMANOVA (Anderson 2001): CLR-transform each sample's
    three-part composition, take Euclidean (Aitchison) distances, and compute the
    pseudo-F ratio of between- to within-group dispersion. Significance is an exact
    label permutation over all C(N, n_old) relabellings (252 at 5 vs 5), one-sided
    (larger F = more separation). Returns a one-row summary."""
    df = pd.read_csv(PROP_DIR / f"_per_sample_{strategy}.csv")
    sub = df[df["celltype"].isin(OLIGO_TYPES)].copy()
    pivot = sub.pivot(index="sample_id", columns="celltype", values="n_in_lesion")[
        OLIGO_TYPES
    ]
    groups = (
        sub.drop_duplicates("sample_id")
        .set_index("sample_id")["group"]
        .reindex(pivot.index)
    )
    keep = groups.isin(["Young", "Old"])
    pivot, groups = pivot[keep], groups[keep]
    counts = pivot.to_numpy(float)
    labels = (groups.to_numpy() == "Old").astype(int)
    n = len(labels)
    n_old = int(labels.sum())
    n_young = n - n_old

    clr = _clr(counts)
    diff = clr[:, None, :] - clr[None, :, :]
    d2 = (diff**2).sum(-1)  # NxN squared Aitchison distances
    triu = np.triu_indices(n, 1)
    ss_total = d2[triu].sum() / n

    def pseudo_f(lab: np.ndarray) -> float:
        ss_within = 0.0
        for g in (0, 1):
            idx = np.where(lab == g)[0]
            if len(idx) > 1:
                iu = np.triu_indices(len(idx), 1)
                ss_within += d2[np.ix_(idx, idx)][iu].sum() / len(idx)
        ss_between = ss_total - ss_within
        return (ss_between / 1.0) / (ss_within / (n - 2))  # a = 2 groups

    f_obs = pseudo_f(labels)
    ge = tot = 0
    for old_idx in combinations(range(n), n_old):
        lab = np.zeros(n, int)
        lab[list(old_idx)] = 1
        tot += 1
        ge += pseudo_f(lab) >= f_obs - 1e-12
    return pd.DataFrame(
        [
            {
                "test": "PERMANOVA (CLR / Aitchison, exact label permutation)",
                "n_young": n_young,
                "n_old": n_old,
                "pseudo_F": round(f_obs, 4),
                "p_value": round(ge / tot, 4),
                "n_permutations": tot,
            }
        ]
    )


def plot_proportions(
    prop_df: pd.DataFrame, test_df: pd.DataFrame, strategy: str, out_path: Path
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
    rng = np.random.default_rng(42)

    for ax, ct in zip(axes, OLIGO_TYPES):
        young = prop_df.loc[prop_df["group"] == "Young", ct].dropna().values
        old = prop_df.loc[prop_df["group"] == "Old", ct].dropna().values

        ax.boxplot(
            [young, old],
            positions=[0, 1],
            widths=0.55,
            showfliers=False,
            patch_artist=True,
            boxprops=dict(facecolor="#e8e8e8", edgecolor="black", alpha=0.6),
            medianprops=dict(color="black", linewidth=2),
            whiskerprops=dict(color="black"),
            capprops=dict(color="black"),
        )

        for vals, x, color in [
            (young, 0, plot_style.GROUP_COLORS["Young"]),
            (old, 1, plot_style.GROUP_COLORS["Old"]),
        ]:
            jitter = rng.uniform(-0.12, 0.12, len(vals))
            ax.scatter(
                x + jitter,
                vals,
                color=color,
                s=70,
                alpha=0.85,
                edgecolor="black",
                linewidth=0.6,
                zorder=3,
            )

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Young (n=5)", "Old (n=5)"])

        row = test_df.loc[test_df["celltype"] == ct].iloc[0]
        p = row["p_value_two_sided"]
        med_y = row["median_young"]
        med_o = row["median_old"]
        ax.set_title(
            f"{ct}\nmedians: Y={med_y:.3f}, O={med_o:.3f}\nMann-Whitney p = {p:.3g}",
            fontsize=10,
        )
        ax.set_ylabel("proportion within oligo lineage in lesion")
        ax.grid(True, axis="y", linestyle="--", alpha=0.3)
        ax.set_ylim(0, max(np.concatenate([young, old, [0.01]])) * 1.15)

    fig.suptitle(f"Lesion oligo-lineage composition: Young vs Old ({strategy})")
    plt.tight_layout()
    plot_style.save(fig, out_path.stem, fig_dir=out_path.parent)


def main() -> int:
    args = parse_args()
    plot_style.apply()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    for strategy in args.strategies:
        print(f"\n=== {strategy} ===")
        prop = compute_lineage_proportions(strategy)
        prop_path = PROP_DIR / f"_lineage_proportions_{strategy}.csv"
        prop.to_csv(prop_path, index=False)

        test = test_young_vs_old(prop)
        test_path = PROP_DIR / f"_young_vs_old_test_{strategy}.csv"
        test.to_csv(test_path, index=False)

        print("Per-sample within-lineage proportions:")
        print(prop.to_string(index=False))
        print()
        print("Mann-Whitney U test (Young vs Old, two-sided):")
        print(test.to_string(index=False))

        permanova = permanova_young_vs_old(strategy)
        permanova_path = PROP_DIR / f"_composition_permanova_{strategy}.csv"
        permanova.to_csv(permanova_path, index=False)
        print("\nGlobal composition test (CLR + PERMANOVA, closure-aware):")
        print(permanova.to_string(index=False))

        out_path = FIG_DIR / f"lineage_young_vs_old_{strategy}.png"
        plot_proportions(prop, test, strategy, out_path)
        print(f"  plot -> {out_path.relative_to(REPO_ROOT)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
