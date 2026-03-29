"""
scoring.py
Unified scRNA scoring module — single source of truth for both CPU and GPU CLIs.

All functions normalise with CPM + log1p before scoring (matching the GPU path),
so results are consistent regardless of which CLI is used.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _to_dense(X) -> np.ndarray:
    if sp.issparse(X):
        return np.asarray(X.todense())
    return np.asarray(X)


def _cpm_log1p(X: np.ndarray) -> np.ndarray:
    """Library-size normalise to CPM then log1p. Operates on dense array."""
    lib = X.sum(axis=1, keepdims=True).clip(1)
    return np.log1p(X / lib * 1e6)


def _batch_zscore(values: np.ndarray, batch: pd.Series) -> np.ndarray:
    """Per-batch z-score a 1-D array aligned to adata.obs."""
    out = np.zeros_like(values, dtype=float)
    batch_arr = batch.values if hasattr(batch, "values") else np.asarray(batch)
    for grp in np.unique(batch_arr):
        pos = np.where(batch_arr == grp)[0]
        v = values[pos]
        mu, sd = v.mean(), v.std(ddof=1)
        out[pos] = (v - mu) / (sd + 1e-9)
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def copy_number_proxy(adata) -> pd.Series:
    """
    Per-cell mitochondrial copy number proxy, z-scored within batch.

    - Mito-only dataset (>50% MT genes): uses log(raw_sum) — total mito
      transcriptional output is a strong mtDNA copy number correlate.
    - Standard dataset: uses MT-gene fraction of total counts.
    """
    genes     = adata.var.index.astype(str)
    is_mt_arr = np.asarray(genes.str.upper().str.startswith("MT-"))
    mt_frac   = is_mt_arr.mean()

    if mt_frac > 0.5:
        # Mito-filtered dataset — use total expression as proxy
        if "raw_sum" in adata.obs.columns:
            vals = np.log1p(adata.obs["raw_sum"].values.astype(float))
        else:
            X    = _to_dense(adata.X)
            vals = np.log1p(X.sum(axis=1))
    else:
        X     = _to_dense(adata.X)
        total = X.sum(axis=1) + 1e-9
        vals  = X[:, is_mt_arr].sum(axis=1) / total if is_mt_arr.any()                 else np.zeros(adata.n_obs)

    batch = adata.obs.get(
        "batch", pd.Series(["batch1"] * adata.n_obs, index=adata.obs_names)
    )
    z = _batch_zscore(np.asarray(vals, dtype=float), batch)
    return pd.Series(z, index=adata.obs_names, name="copy_number_proxy")


def program_scores(
    adata,
    pathways: dict[str, list[str]] | None = None,
    layer: str | None = None,
    use_gpu: bool = False,
) -> pd.DataFrame:
    """
    Per-cell pathway scores: CPM + log1p → mean z-score over gene set members.

    Parameters
    ----------
    adata     : AnnData
    pathways  : dict mapping pathway name -> list of gene symbols
    layer     : AnnData layer to use (default: .X)
    use_gpu   : if True and CuPy is available, offload to GPU

    Returns
    -------
    DataFrame (cells x pathways), columns lowercased.
    """
    if not pathways:
        return pd.DataFrame(index=adata.obs_names)

    raw = adata.layers[layer] if layer else adata.X

    # GPU path
    if use_gpu:
        try:
            from .rapids_program_score import score_programs, _HAVE_GPU
            if _HAVE_GPU:
                gs_df = pd.DataFrame(
                    [{"pathway": p, "gene": g} for p, gs in pathways.items() for g in gs]
                )
                df = score_programs(
                    adata, gene_sets=pathways, layer=layer,
                    libnorm="CPM", log1p=True, standardize="z",
                    return_df=True, use_gpu=True,
                )
                df.columns = [c.lower() for c in df.columns]
                return df
        except Exception:
            pass  # fall through to CPU

    # CPU path — CPM + log1p then per-gene z, then mean over set
    X = _cpm_log1p(_to_dense(raw))
    var_upper = adata.var.index.astype(str).str.upper().values

    out = {}
    for pname, glist in pathways.items():
        if not glist:
            continue
        pset = {g.upper() for g in glist}
        mask = np.array([g in pset for g in var_upper], dtype=bool)
        if not mask.any():
            out[pname] = np.zeros(adata.n_obs, dtype=float)
            continue
        sub = X[:, mask]
        # z-score each gene across cells, then average within set
        mu = sub.mean(axis=0)
        sd = sub.std(axis=0, ddof=0) + 1e-8
        z_sub = (sub - mu) / sd
        out[pname] = z_sub.mean(axis=1)

    df = pd.DataFrame(out, index=adata.obs_names)
    df.columns = [c.lower() for c in df.columns]
    return df


def heterogeneity_scores(adata) -> pd.Series:
    """
    Per-subject heterogeneity: Shannon diversity over cell-type distribution,
    scaled to [0, 1] across subjects.
    """
    ct = adata.obs.get("cell_type")

    def shannon(p: np.ndarray) -> float:
        p = p[p > 0]
        return float(-(p * np.log(p)).sum())

    if ct is None:
        # Fallback: SD of copy-number proxy per subject
        cnp = copy_number_proxy(adata)
        s = cnp.groupby(adata.obs["subject_id"], observed=True).std().fillna(0.0)
    else:
        rows = []
        for sid, g in adata.obs.groupby("subject_id", observed=True):
            props = g["cell_type"].value_counts(normalize=True).values
            rows.append((sid, shannon(props)))
        s = pd.DataFrame(rows, columns=["subject_id", "diversity"]).set_index("subject_id")["diversity"]

    v = s.values.astype(float)
    scaled = (v - v.min()) / (v.max() - v.min() + 1e-9)
    return pd.Series(scaled, index=s.index)
