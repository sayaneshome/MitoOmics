from __future__ import annotations
import pandas as pd
import numpy as np
from sklearn.linear_model import RidgeCV
from .config import DEFAULT_WEIGHTS


def _scale01(s: pd.Series) -> pd.Series:
    v = s.values.astype(float)
    lo, hi = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
    scaled = (np.clip(v, lo, hi) - lo) / (hi - lo + 1e-9)
    return pd.Series(scaled, index=s.index)


def combine_components(
    subject_df: pd.DataFrame,
    weights=DEFAULT_WEIGHTS,
    fit_ridge: bool = False,
    targets: pd.Series | None = None,
) -> pd.DataFrame:
    df = subject_df.copy().set_index("subject_id")

    cols = [
        c for c in ["copy_number", "fusion_fission", "mitophagy", "heterogeneity"]
        if c in df.columns
    ]
    if not cols:
        raise ValueError(
            "No components to combine. Need at least one of: "
            "copy_number, fusion_fission, mitophagy, heterogeneity."
        )

    X = df[cols].apply(_scale01)

    if fit_ridge and targets is not None:
        y = targets.reindex(X.index)
        # pandas 2.x: use .ffill()/.bfill() not fillna(method=...)
        y = y.ffill().bfill().fillna(y.mean())
        model = RidgeCV(alphas=[0.1, 1.0, 10.0], cv=min(5, len(X)))
        model.fit(X.values, y.values)
        w = pd.Series(model.coef_, index=cols)
    else:
        cfg = {
            "copy_number": weights.copy_number,
            "fusion_fission": weights.fusion_fission,
            "mitophagy": weights.mitophagy,
            "heterogeneity": weights.heterogeneity,
        }
        w = pd.Series({c: cfg[c] for c in cols})
        w = w / (w.sum() + 1e-9)

    mhi = (X * w).sum(axis=1).rename("MHI")
    out = df.copy()
    out["MHI"] = mhi
    out = out.reset_index()

    if fit_ridge and targets is not None:
        for c in w.index:
            out[f"weight_{c}"] = w.loc[c]

    return out
