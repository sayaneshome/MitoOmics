"""
Case Study 2: Full MitoCarta pathway profiling (141 pathways).
Computes per-subject scores for all 141 MitoCarta pathways,
clusters subjects, and identifies pathways discriminating disease subtypes.

Usage:
  python mitoomics_gpu/tools/casestudy_pathway_profiling.py \
    --scrna           mitoomics_gpu/data/scrna.mito.h5ad \
    --mitocarta-table mitoomics_gpu/data/mitocarta3_table.csv \
    --outdir          results/casestudy_pathways/
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import anndata
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from mitoomics_gpu.io import load_mitocarta_pathways
from mitoomics_gpu.scoring import program_scores


DISEASE_COLORS = {
    "triple-negative breast carcinoma":          "#E24B4A",
    "estrogen-receptor positive breast cancer":  "#EF9F27",
    "invasive ductal breast carcinoma":          "#3B8BD4",
    "HER2 positive breast carcinoma":            "#7F77DD",
    "breast cancer":                             "#1D9E75",
    "breast carcinoma":                          "#0F6E56",
    "invasive lobular breast carcinoma":         "#D4537E",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scrna",           required=True)
    ap.add_argument("--mitocarta-table", required=True)
    ap.add_argument("--outdir",          default="results/casestudy_pathways/")
    ap.add_argument("--n-top-pathways",  type=int, default=30)
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] Loading data...")
    adata    = anndata.read_h5ad(args.scrna)
    pathways = load_mitocarta_pathways(args.mitocarta_table)
    print(f"[INFO] {len(pathways)} pathways loaded, {adata.n_obs} cells, {adata.n_vars} genes")

    print("[INFO] Scoring all 141 MitoCarta pathways per cell...")
    prog = program_scores(adata, pathways=pathways)
    adata.obs = adata.obs.join(prog, rsuffix="_path")
    path_cols = prog.columns.tolist()

    print("[INFO] Aggregating to subject level...")
    subj = adata.obs.groupby("subject_id", observed=True)[path_cols].mean()
    # Attach disease label
    disease = adata.obs.groupby("subject_id", observed=True)["disease"].agg(
        lambda x: x.value_counts().index[0]
    )
    subj["disease"] = disease
    subj.to_csv(outdir / "pathway_profiles.csv")
    print(f"[OK] Pathway profiles → {outdir}/pathway_profiles.csv  ({len(subj)} subjects × {len(path_cols)} pathways)")

    # ── Top discriminating pathways (ANOVA-like variance) ──
    X = subj[path_cols].values.astype(float)
    pathway_variance = np.nanvar(X, axis=0)
    top_idx = np.argsort(pathway_variance)[::-1][:args.n_top_pathways]
    top_paths = [path_cols[i] for i in top_idx]

    top_df = subj[top_paths + ["disease"]]
    top_df.to_csv(outdir / "top_pathways.csv")

    # ── Heatmap of top pathways ──
    disease_order = [
        "triple-negative breast carcinoma",
        "estrogen-receptor positive breast cancer",
        "invasive ductal breast carcinoma",
        "HER2 positive breast carcinoma",
        "breast cancer",
        "breast carcinoma",
        "invasive lobular breast carcinoma",
    ]
    subj_sorted = pd.concat([
        subj[subj["disease"] == d] for d in disease_order if d in subj["disease"].values
    ])
    mat = subj_sorted[top_paths].values
    mat_scaled = StandardScaler().fit_transform(mat)

    row_colors = [DISEASE_COLORS.get(d, "#888780") for d in subj_sorted["disease"]]
    short_paths = [p.split(">")[-1].strip()[:35] for p in top_paths]

    fig, (ax_cb, ax_hm) = plt.subplots(
        1, 2, figsize=(16, 8),
        gridspec_kw={"width_ratios": [0.025, 1]}
    )
    im = ax_hm.imshow(mat_scaled.T, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax_hm.set_yticks(range(len(short_paths))); ax_hm.set_yticklabels(short_paths, fontsize=7)
    ax_hm.set_xticks([]); ax_hm.set_xlabel(f"Subjects (n={len(subj_sorted)}, sorted by disease)")
    ax_hm.set_title(f"Top {args.n_top_pathways} MitoCarta pathways by variance across subjects")
    # Disease color bar on left
    for i, c in enumerate(row_colors):
        ax_cb.barh(i, 1, color=c, edgecolor="none")
    ax_cb.set_xlim(0,1); ax_cb.set_ylim(-0.5, len(row_colors)-0.5)
    ax_cb.axis("off")
    plt.colorbar(im, ax=ax_hm, label="z-score", shrink=0.4, pad=0.01)
    # Legend
    handles = [mpatches.Patch(color=DISEASE_COLORS.get(d,"#888780"), label=d[:35])
               for d in disease_order if d in subj["disease"].values]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=7,
               bbox_to_anchor=(0.55, -0.05))
    fig.tight_layout()
    fig.savefig(outdir / "pathway_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[OK] → {outdir}/pathway_heatmap.png")

    # ── PCA of pathway profiles ──
    X_top = subj[top_paths].values.astype(float)
    X_top = StandardScaler().fit_transform(X_top)
    pca   = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_top)

    fig, ax = plt.subplots(figsize=(8, 6))
    for d in disease_order:
        mask = subj["disease"].values == d
        if mask.sum() == 0: continue
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=DISEASE_COLORS.get(d,"#888780"), label=d[:35],
                   s=60, alpha=0.85, edgecolors="white", linewidths=0.5)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("PCA of MitoCarta pathway profiles")
    ax.legend(fontsize=7, loc="best")
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "pathway_pca.png", dpi=180)
    plt.close()
    print(f"[OK] → {outdir}/pathway_pca.png")

    # ── Top 10 most variable pathways summary ──
    print(f"\n── Top 10 most variable pathways across subjects ──")
    for i, col in enumerate(top_paths[:10]):
        pretty = col.replace("_"," ").replace(">"," > ")
        print(f"  {i+1:2d}. {pretty[:60]}  (variance={pathway_variance[top_idx[i]]:.4f})")

    # ── Per-disease pathway means ──
    disease_path = subj.groupby("disease", observed=True)[top_paths[:15]].mean().T
    disease_path.index = [p.split(">")[-1].strip()[:30] for p in disease_path.index]
    disease_path.to_csv(outdir / "disease_pathway_means.csv")
    print(f"[OK] Disease pathway means → {outdir}/disease_pathway_means.csv")


if __name__ == "__main__":
    main()
