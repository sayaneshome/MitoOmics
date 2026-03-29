"""Tests for scoring.py — copy number proxy and program scores."""
import numpy as np
import pandas as pd
import pytest
import anndata as ad

from mitoomics_gpu.scoring import copy_number_proxy, program_scores, heterogeneity_scores


def _make_adata(n_cells=30, mt_frac=0.1, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(80)] + [f"MT-CO{i}" for i in range(10)]
    X = rng.poisson(1.0, size=(n_cells, len(genes))).astype(float)
    adata = ad.AnnData(X=X)
    adata.var.index = genes
    adata.obs["subject_id"] = np.repeat(["S1", "S2", "S3"], n_cells // 3)
    adata.obs["cell_type"]  = rng.choice(["T", "B", "NK"], size=n_cells)
    adata.obs["batch"]      = rng.choice(["b1", "b2"], size=n_cells)
    return adata


def test_copy_number_proxy_shape():
    adata = _make_adata()
    result = copy_number_proxy(adata)
    assert len(result) == adata.n_obs
    assert result.name == "copy_number_proxy"


def test_copy_number_proxy_detects_mt_genes():
    """Cells with higher MT expression should have higher (positive) proxy values."""
    adata = _make_adata()
    mt_mask = adata.var.index.str.startswith("MT-")
    # Boost MT expression for first 10 cells
    adata.X[:10, mt_mask] += 50
    result = copy_number_proxy(adata)
    # boosted cells should tend to be above zero
    assert result[:10].mean() > result[10:].mean()


def test_program_scores_returns_correct_shape():
    adata = _make_adata()
    pathways = {"fusion": ["G0", "G1", "G2"], "mitophagy": ["G5", "G6"]}
    df = program_scores(adata, pathways=pathways)
    assert df.shape == (adata.n_obs, 2)
    assert set(df.columns) == {"fusion", "mitophagy"}


def test_program_scores_missing_genes_returns_zeros():
    adata = _make_adata()
    pathways = {"phantom": ["NOTEXIST1", "NOTEXIST2"]}
    df = program_scores(adata, pathways=pathways)
    assert (df["phantom"] == 0).all()


def test_program_scores_empty_pathways():
    adata = _make_adata()
    df = program_scores(adata, pathways={})
    assert df.empty


def test_heterogeneity_scores_range():
    adata = _make_adata(n_cells=30)
    s = heterogeneity_scores(adata)
    assert s.min() >= 0.0
    assert s.max() <= 1.0 + 1e-6
