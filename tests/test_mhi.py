"""Tests for mhi.py — MHI computation."""
import pandas as pd
import numpy as np
import pytest

from mitoomics_gpu.mhi import combine_components


def _subject_df(n=10, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "subject_id": [f"S{i}" for i in range(n)],
        "copy_number":   rng.uniform(0, 1, n),
        "fusion_fission":rng.uniform(0, 1, n),
        "mitophagy":     rng.uniform(0, 1, n),
        "heterogeneity": rng.uniform(0, 1, n),
    })


def test_mhi_in_range():
    df = _subject_df()
    out = combine_components(df)
    assert "MHI" in out.columns
    assert out["MHI"].between(0, 1).all()


def test_mhi_subject_count():
    df = _subject_df(n=15)
    out = combine_components(df)
    assert len(out) == 15


def test_missing_component_still_runs():
    df = _subject_df()[["subject_id", "copy_number", "mitophagy"]]
    out = combine_components(df)
    assert "MHI" in out.columns


def test_no_components_raises():
    df = pd.DataFrame({"subject_id": ["S1", "S2"]})
    with pytest.raises(ValueError, match="No components"):
        combine_components(df)


def test_ridge_with_targets():
    df = _subject_df(n=20)
    targets = pd.Series(
        np.random.default_rng(42).uniform(0, 1, 20),
        index=df["subject_id"],
        name="target",
    )
    out = combine_components(df, fit_ridge=True, targets=targets)
    assert "MHI" in out.columns
    # Weight columns should appear
    assert any(c.startswith("weight_") for c in out.columns)
