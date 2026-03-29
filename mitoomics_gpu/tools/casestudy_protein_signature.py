"""
Case Study 3: EV Mitochondrial Protein Signature Analysis.
Identifies which mito proteins are most differentially abundant
across cancer types and builds a minimal diagnostic signature.

Usage:
  python mitoomics_gpu/tools/casestudy_protein_signature.py \
    --proteomics mitoomics_gpu/data/ev_human.csv \
    --genesets   mitoomics_gpu/data/genesets_curated.csv \
    --outdir     results/casestudy_proteins/
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats
from itertools import combinations

PALETTE = ["#1D9E75","#3B8BD4","#EF9F27","#D85A30","#7F77DD",
           "#D4537E","#0F6E56","#854F0B","#185FA5","#534AB7"]


def parse_cancer_type(s: str) -> str:
    s = s.lower()
    if "breast" in s or "mda" in s or "mcf" in s:                    return "Breast"
    if "pancreatic" in s or "bxpc" in s or "hpaf" in s or "panc" in s: return "Pancreatic"
    if "lung" in s or "h358" in s or "h1299" in s:                   return "Lung"
    if "ovarian" in s or "skov" in s:                                 return "Ovarian"
    if "prostate" in s:                                               return "Prostate"
    if "plasma" in s:                                                  return "Patient plasma"
    if "colon" in s or "hct" in s:                                    return "Colorectal"
    if "glioma" in s or "gbm" in s or "u87" in s:                    return "Glioma"
    if "normal" in s or "healthy" in s:                               return "Normal/Control"
    return "Other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteomics", required=True)
    ap.add_argument("--genesets",   required=True)
    ap.add_argument("--outdir",     default="results/casestudy_proteins/")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] Loading data...")
    ev = pd.read_csv(args.proteomics)
    gs = pd.read_csv(args.genesets)
    gs_prots = set(gs["gene"].str.upper())

    ev["cancer_type"]    = ev["subject_id"].apply(parse_cancer_type)
    ev["protein_upper"]  = ev["protein"].str.upper()

    # Filter to mito pathway proteins only
    ev_mito = ev[ev["protein_upper"].isin(gs_prots)].copy()
    ev_mito["log_abundance"] = np.log1p(ev_mito["abundance"])
    print(f"[INFO] {ev_mito['protein_upper'].nunique()} mito proteins found across {ev_mito['subject_id'].nunique()} subjects")

    # ── Subject-level protein matrix ──
    prot_matrix = ev_mito.pivot_table(
        index="subject_id", columns="protein_upper",
        values="log_abundance", aggfunc="mean"
    ).fillna(0)
    prot_matrix["cancer_type"] = prot_matrix.index.map(
        ev[["subject_id","cancer_type"]].drop_duplicates().set_index("subject_id")["cancer_type"]
    )
    prot_matrix.to_csv(outdir / "protein_matrix.csv")
    print(f"[OK] Protein matrix → {outdir}/protein_matrix.csv  ({prot_matrix.shape})")

    # ── Mean abundance per protein per cancer type ──
    prot_cols = [c for c in prot_matrix.columns if c != "cancer_type"]
    means = prot_matrix.groupby("cancer_type")[prot_cols].mean()
    means.to_csv(outdir / "protein_means_by_cancer.csv")

    # ── Heatmap ──
    means_norm = (means - means.mean()) / (means.std() + 1e-8)
    fig, ax = plt.subplots(figsize=(max(10, len(prot_cols)*0.45), 5))
    im = ax.imshow(means_norm.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(means_norm.index))); ax.set_yticklabels(means_norm.index, fontsize=9)
    ax.set_xticks(range(len(prot_cols))); ax.set_xticklabels(prot_cols, rotation=45, ha="right", fontsize=8)
    ax.set_title("Mito protein abundance across cancer types (z-scored)")
    plt.colorbar(im, ax=ax, label="z-score", shrink=0.6)
    fig.tight_layout()
    fig.savefig(outdir / "protein_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[OK] → {outdir}/protein_heatmap.png")

    # ── Differential proteins: top cancer type vs others ──
    cancer_types = [ct for ct in prot_matrix["cancer_type"].unique()
                    if (prot_matrix["cancer_type"]==ct).sum() >= 3]
    diff_rows = []
    for ct in cancer_types:
        mask = prot_matrix["cancer_type"] == ct
        X_in  = prot_matrix.loc[mask,  prot_cols].values
        X_out = prot_matrix.loc[~mask, prot_cols].values
        for j, prot in enumerate(prot_cols):
            a, b = X_in[:, j], X_out[:, j]
            if len(a) < 2 or len(b) < 2: continue
            U, p = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
            fc = np.mean(a) - np.mean(b)
            diff_rows.append({"cancer_type":ct,"protein":prot,
                               "mean_in":round(float(np.mean(a)),3),
                               "mean_out":round(float(np.mean(b)),3),
                               "fold_change":round(float(fc),3),
                               "p_value":float(p)})
    diff_df = pd.DataFrame(diff_rows).sort_values("p_value")
    diff_df.to_csv(outdir / "differential_proteins.csv", index=False)
    print(f"[OK] Differential proteins → {outdir}/differential_proteins.csv")

    # ── Top 5 differential proteins per cancer type ──
    print("\n── Top differential mito proteins per cancer type ──")
    for ct in cancer_types:
        top5 = diff_df[diff_df["cancer_type"]==ct].nsmallest(5,"p_value")
        print(f"\n  {ct}:")
        for _, row in top5.iterrows():
            direction = "↑" if row["fold_change"] > 0 else "↓"
            print(f"    {direction} {row['protein']:12s}  fc={row['fold_change']:+.3f}  p={row['p_value']:.3f}")

    # ── Bar plot: key mito proteins across cancer types ──
    key_prots = diff_df.nsmallest(12,"p_value")["protein"].unique()[:8]
    sub = prot_matrix[list(key_prots) + ["cancer_type"]]
    ct_order = means[list(key_prots)].mean(axis=1).sort_values(ascending=False).index.tolist()

    fig, axes = plt.subplots(2, 4, figsize=(14, 6), sharey=False)
    axes = axes.flatten()
    for i, prot in enumerate(key_prots):
        ax = axes[i]
        vals  = [sub[sub["cancer_type"]==ct][prot].dropna().values for ct in ct_order]
        meds  = [np.median(v) if len(v)>0 else 0 for v in vals]
        colors= [PALETTE[j%len(PALETTE)] for j in range(len(ct_order))]
        bars  = ax.bar(range(len(ct_order)), meds, color=colors, alpha=0.85)
        ax.set_xticks(range(len(ct_order)))
        ax.set_xticklabels(ct_order, rotation=40, ha="right", fontsize=6)
        ax.set_title(prot, fontsize=9, fontweight="bold")
        ax.set_ylabel("log abundance" if i%4==0 else "")
        ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("Key mito EV proteins across cancer types", fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / "key_proteins_by_cancer.png", dpi=180)
    plt.close()
    print(f"[OK] → {outdir}/key_proteins_by_cancer.png")


if __name__ == "__main__":
    main()
