"""
celltype_mhi.py
Compute MHI per cell type per subject and produce summary tables + figures.

Usage:
  python -m mitoomics_gpu.tools.celltype_mhi \
    --scrna           mitoomics_gpu/data/scrna.mito.h5ad \
    --mitocarta-table mitoomics_gpu/data/mitocarta3_table.csv \
    --outdir          results/
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import anndata
from ..io import load_mitocarta_pathways
from ..scoring import copy_number_proxy, program_scores, heterogeneity_scores
from ..mhi import combine_components


PALETTE = [
    "#1D9E75","#3B8BD4","#EF9F27","#D85A30","#7F77DD",
    "#D4537E","#0F6E56","#854F0B","#185FA5","#534AB7",
]


def compute_celltype_mhi(adata, pathways: dict) -> pd.DataFrame:
    """
    For each (subject_id, cell_type) group compute component scores and MHI.
    Returns a long DataFrame: subject_id, cell_type, copy_number,
    fusion_fission, mitophagy, heterogeneity, MHI.
    """
    cnp = copy_number_proxy(adata)
    adata.obs["copy_number_proxy"] = cnp.values
    prog = program_scores(adata, pathways=pathways)
    adata.obs = adata.obs.join(prog, rsuffix="_prog")

    def _pick(df, kws):
        return [c for c in df.columns if any(k in c for k in kws)]

    ff_cols   = _pick(prog, ["fusion", "fission"])
    mito_cols = _pick(prog, ["mitophagy"])

    rows = []
    for (sid, ct), grp in adata.obs.groupby(
        ["subject_id", "cell_type"], observed=True
    ):
        if len(grp) < 3:   # skip tiny groups
            continue
        cnp_val  = grp["copy_number_proxy"].mean()
        ff_val   = grp[ff_cols].mean().mean() if ff_cols else np.nan
        mito_val = grp[mito_cols].mean().mean() if mito_cols else np.nan
        rows.append({
            "subject_id":   sid,
            "cell_type":    ct,
            "disease":      grp["disease"].iloc[0] if "disease" in grp.columns else "unknown",
            "n_cells":      len(grp),
            "copy_number":  cnp_val,
            "fusion_fission": ff_val,
            "mitophagy":    mito_val,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Compute MHI per row
    mhi_rows = []
    for _, row in df.iterrows():
        sub = pd.DataFrame([{
            "subject_id":    f"{row['subject_id']}|{row['cell_type']}",
            "copy_number":   row["copy_number"],
            "fusion_fission":row["fusion_fission"],
            "mitophagy":     row["mitophagy"],
        }])
        try:
            out = combine_components(sub)
            mhi_rows.append(out["MHI"].iloc[0])
        except Exception:
            mhi_rows.append(np.nan)
    df["MHI"] = mhi_rows
    return df


def plot_celltype_heatmap(ct_df: pd.DataFrame, outdir: Path) -> str:
    """Heatmap: cell types (rows) × subjects (cols), color = MHI."""
    pivot = ct_df.pivot_table(
        index="cell_type", columns="subject_id", values="MHI"
    )
    # Keep top 15 most represented cell types
    top_ct = ct_df.groupby("cell_type", observed=True)["n_cells"].sum().nlargest(15).index
    pivot = pivot.loc[pivot.index.intersection(top_ct)]

    fig, ax = plt.subplots(figsize=(min(20, pivot.shape[1] * 0.35 + 2), 6))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=90, fontsize=6)
    ax.set_xlabel("Subject")
    ax.set_ylabel("Cell type")
    ax.set_title("MHI by cell type × subject")
    plt.colorbar(im, ax=ax, label="MHI", shrink=0.6)
    fig.tight_layout()
    path = outdir / "celltype_mhi_heatmap.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    return str(path)


def plot_celltype_boxplot(ct_df: pd.DataFrame, outdir: Path) -> str:
    """Boxplot: MHI distribution per cell type, colored by disease."""
    top_ct = (
        ct_df.groupby("cell_type", observed=True)["n_cells"]
        .sum().nlargest(12).index.tolist()
    )
    sub = ct_df[ct_df["cell_type"].isin(top_ct)]
    order = (
        sub.groupby("cell_type", observed=True)["MHI"]
        .median().sort_values(ascending=False).index.tolist()
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    data = [sub[sub["cell_type"] == ct]["MHI"].dropna().values for ct in order]
    bp = ax.boxplot(data, patch_artist=True, widths=0.5)
    for patch, color in zip(bp["boxes"], PALETTE):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("MHI")
    ax.set_title("MHI distribution by cell type")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = outdir / "celltype_mhi_boxplot.png"
    fig.savefig(path, dpi=180)
    plt.close()
    return str(path)


def plot_disease_celltype(ct_df: pd.DataFrame, outdir: Path) -> str:
    """Grouped bar: mean MHI per (disease, cell_type) for top cell types."""
    top_ct = (
        ct_df.groupby("cell_type", observed=True)["n_cells"]
        .sum().nlargest(6).index.tolist()
    )
    sub = ct_df[ct_df["cell_type"].isin(top_ct)]
    pivot = sub.pivot_table(
        index="disease", columns="cell_type", values="MHI", aggfunc="mean"
    )[top_ct]

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(pivot.index))
    w = 0.12
    for i, ct in enumerate(top_ct):
        ax.bar(x + i * w, pivot[ct].values, width=w,
               label=ct, color=PALETTE[i % len(PALETTE)], alpha=0.85)
    ax.set_xticks(x + w * (len(top_ct) - 1) / 2)
    ax.set_xticklabels(pivot.index, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Mean MHI")
    ax.set_title("Mean MHI by disease subtype and cell type")
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = outdir / "disease_celltype_mhi.png"
    fig.savefig(path, dpi=180)
    plt.close()
    return str(path)


def main():
    ap = argparse.ArgumentParser(description="Cell-type level MHI analysis")
    ap.add_argument("--scrna",           required=True)
    ap.add_argument("--mitocarta-table", required=True)
    ap.add_argument("--outdir",          default="results/")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] Loading data...")
    adata    = anndata.read_h5ad(args.scrna)
    pathways = load_mitocarta_pathways(args.mitocarta_table)

    print("[INFO] Computing cell-type MHI...")
    ct_df = compute_celltype_mhi(adata, pathways)

    out_csv = outdir / "celltype_mhi.csv"
    ct_df.to_csv(out_csv, index=False)
    print(f"[OK] Cell-type MHI → {out_csv}  ({len(ct_df)} rows)")

    print("[INFO] Generating figures...")
    plot_celltype_heatmap(ct_df, outdir)
    print(f"[OK] → {outdir}/celltype_mhi_heatmap.png")
    plot_celltype_boxplot(ct_df, outdir)
    print(f"[OK] → {outdir}/celltype_mhi_boxplot.png")
    plot_disease_celltype(ct_df, outdir)
    print(f"[OK] → {outdir}/disease_celltype_mhi.png")

    # Summary table
    summary = ct_df.groupby(["cell_type", "disease"], observed=True)["MHI"].agg(
        ["mean","median","std","count"]
    ).round(4).reset_index()
    summary.to_csv(outdir / "celltype_mhi_summary.csv", index=False)
    print(f"[OK] Summary table  → {outdir}/celltype_mhi_summary.csv")

    print("\n── Top 10 cell types by mean MHI ──")
    top = ct_df.groupby("cell_type", observed=True)["MHI"].mean().sort_values(ascending=False).head(10)
    print(top.round(4).to_string())

    print("\n── Mean MHI by disease subtype ──")
    dis = ct_df.groupby("disease", observed=True)["MHI"].mean().sort_values(ascending=False)
    print(dis.round(4).to_string())


if __name__ == "__main__":
    main()
