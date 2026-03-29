"""Tests for stats.py — differential MHI testing."""
import numpy as np
import pandas as pd
import pytest

from mitoomics_gpu.stats import mhi_group_test, mhi_correlation, _bh_correct


def _results_df(seed=0):
    rng = np.random.default_rng(seed)
    n = 30
    groups = np.repeat(["tumor", "normal", "benign"], 10)
    # tumor has noticeably higher MHI
    mhi = np.concatenate([
        rng.uniform(0.6, 1.0, 10),
        rng.uniform(0.2, 0.6, 10),
        rng.uniform(0.3, 0.7, 10),
    ])
    return pd.DataFrame({"subject_id": [f"S{i}" for i in range(n)],
                         "MHI": mhi, "disease": groups})


def test_mhi_group_test_returns_df():
    df = _results_df()
    out = mhi_group_test(df, group_col="disease")
    assert isinstance(out, pd.DataFrame)
    assert "p_value" in out.columns
    assert "p_adj" in out.columns


def test_mhi_group_test_n_pairs():
    df = _results_df()
    out = mhi_group_test(df, group_col="disease")
    # C(3,2) = 3 pairs
    assert len(out) == 3


def test_mhi_group_test_missing_col():
    df = _results_df()
    with pytest.raises(ValueError):
        mhi_group_test(df, group_col="nonexistent")


def test_mhi_group_test_single_group_empty():
    df = pd.DataFrame({"subject_id": ["S1","S2"], "MHI": [0.5,0.6], "disease": ["tumor","tumor"]})
    out = mhi_group_test(df, group_col="disease")
    assert out.empty


def test_mhi_correlation():
    df = _results_df()
    df["age"] = np.random.default_rng(1).uniform(30, 80, len(df))
    res = mhi_correlation(df, covariate_col="age")
    assert "rho" in res
    assert "p_value" in res


def test_bh_correct_monotone():
    p = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
    adj = _bh_correct(p)
    assert np.all(adj[:-1] <= adj[1:] + 1e-12)
    assert np.all(adj <= 1.0)
