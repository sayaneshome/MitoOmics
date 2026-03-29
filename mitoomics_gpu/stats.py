"""
stats.py
Differential MHI analysis between groups (e.g. tumor vs normal).

Provides:
  - mhi_group_test()   : Mann-Whitney U + effect size (rank-biserial r) per group pair
  - mhi_correlation()  : Spearman correlation of MHI with a continuous variable
  - append_stats_to_report() : writes a stats section into the markdown report
"""
from __future__ import annotations
import itertools
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats


# ---------------------------------------------------------------------------
# Core tests
# ---------------------------------------------------------------------------

def mhi_group_test(
    results_df: pd.DataFrame,
    group_col: str,
    mhi_col: str = "MHI",
    fdr_method: str = "bh",
) -> pd.DataFrame:
    """
    Pairwise Mann-Whitney U tests for MHI between all groups in `group_col`.

    Returns a DataFrame with columns:
        group_a, group_b, n_a, n_b, median_a, median_b,
        U_stat, p_value, p_adj (BH-corrected), effect_r
    where effect_r is the rank-biserial correlation (∈ [-1, 1]).
    """
    if group_col not in results_df.columns:
        raise ValueError(f"Column '{group_col}' not in results DataFrame.")

    groups = results_df[group_col].dropna().unique()
    if len(groups) < 2:
        return pd.DataFrame()

    rows = []
    for a, b in itertools.combinations(groups, 2):
        xa = results_df.loc[results_df[group_col] == a, mhi_col].dropna().values
        xb = results_df.loc[results_df[group_col] == b, mhi_col].dropna().values
        if len(xa) < 2 or len(xb) < 2:
            continue
        U, p = scipy_stats.mannwhitneyu(xa, xb, alternative="two-sided")
        # Rank-biserial correlation as effect size
        r = 1 - (2 * U) / (len(xa) * len(xb))
        rows.append({
            "group_a": a, "group_b": b,
            "n_a": len(xa), "n_b": len(xb),
            "median_a": float(np.median(xa)),
            "median_b": float(np.median(xb)),
            "U_stat": float(U),
            "p_value": float(p),
            "effect_r": float(r),
        })

    if not rows:
        return pd.DataFrame()

    out = pd.DataFrame(rows)
    out["p_adj"] = _bh_correct(out["p_value"].values)
    out["significant"] = out["p_adj"] < 0.05
    return out.sort_values("p_adj").reset_index(drop=True)


def mhi_correlation(
    results_df: pd.DataFrame,
    covariate_col: str,
    mhi_col: str = "MHI",
) -> dict:
    """
    Spearman correlation between MHI and a continuous covariate.

    Returns dict with keys: rho, p_value, n.
    """
    sub = results_df[[mhi_col, covariate_col]].dropna()
    if len(sub) < 3:
        return {"rho": None, "p_value": None, "n": len(sub)}
    rho, p = scipy_stats.spearmanr(sub[mhi_col], sub[covariate_col])
    return {"rho": float(rho), "p_value": float(p), "n": len(sub)}


# ---------------------------------------------------------------------------
# Report integration
# ---------------------------------------------------------------------------

def append_stats_to_report(
    stats_df: pd.DataFrame,
    report_path: str,
    group_col: str,
) -> None:
    """Append a statistics section to an existing markdown report."""
    section = f"\n\n## Differential MHI Analysis (by {group_col})\n\n"
    if stats_df.empty:
        section += "_Insufficient group sizes for testing._\n"
    else:
        sig = stats_df[stats_df["significant"]]
        section += (
            f"Pairwise Mann-Whitney U tests across "
            f"{stats_df[['group_a','group_b']].apply(tuple,axis=1).nunique()} group pairs. "
            f"**{len(sig)} significant** after BH correction (α = 0.05).\n\n"
        )
        section += stats_df.to_markdown(index=False)
        section += "\n\n_Effect size r: rank-biserial correlation (|r| > 0.3 = medium, > 0.5 = large)._\n"

    with open(report_path, "a", encoding="utf-8") as f:
        f.write(section)


# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

def _bh_correct(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return p_values
    order = np.argsort(p_values)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    adj = p_values * n / ranked
    # Enforce monotonicity (cumulative minimum from the right)
    adj_sorted = adj[order]
    for i in range(n - 2, -1, -1):
        adj_sorted[i] = min(adj_sorted[i], adj_sorted[i + 1])
    adj[order] = adj_sorted
    return np.clip(adj, 0, 1)
