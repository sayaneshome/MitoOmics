"""
Case Study 1: Multi-cancer EV proteomics MHI comparison.
Uses ev_human.csv to compute MHI across 9 cancer types
from EV/MDV proteomic profiles.

Usage:
  python mitoomics_gpu/tools/casestudy_ev_multicancer.py \
    --proteomics mitoomics_gpu/data/ev_human.csv \
    --genesets   mitoomics_gpu/data/genesets_curated.csv \
    --whitelist  mitoomics_gpu/data/ev_whitelist.csv \
    --outdir     results/casestudy_ev/
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats as scipy_stats

PALETTE = ["#1D9E75","#3B8BD4","#EF9F27","#D85A30","#7F77DD",
           "#D4537E","#0F6E56","#854F0B","#185FA5","#534AB7","#888780"]

def parse_cancer_type(subject_id: str) -> str:
    s = subject_id.lower()
    if "breast" in s or "mda" in s or "mcf" in s or "t47d" in s:    return "Breast"
    if "pancreatic" in s or "bxpc" in s or "hpaf" in s or "panc" in s: return "Pancreatic"
    if "lung" in s or "h358" in s or "h1299" in s or "a549" in s:    return "Lung"
    if "ovarian" in s or "skov" in s or "ovary" in s:                return "Ovarian"
    if "prostate" in s or "lncap" in s or "pc3" in s:                return "Prostate"
    if "plasma" in s:                                                  return "Patient plasma"
    if "colon" in s or "hct" in s or "sw480" in s:                   return "Colorectal"
    if "glioma" in s or "gbm" in s or "u87" in s or "glio" in s:    return "Glioma"
    if "normal" in s or "healthy" in s:                               return "Normal/Control"
    if "exo" in s:                                                     return "Exosome control"
    return "Other"

def compute_ev_mhi(ev_df: pd.DataFrame, pathways: dict,
                   whitelist: set | None = None) -> pd.DataFrame:
    """Compute subject-level MHI from EV proteomics."""
    df = ev_df.copy()
    df["protein_upper"] = df["protein"].str.upper()
    if whitelist:
        df = df[df["protein_upper"].isin({p.upper() for p in whitelist})]

    # Normalize abundance per subject
    df["abundance_norm"] = df.groupby("subject_id")["abundance"].transform(
        lambda x: x / (x.sum() + 1e-9)
    )

    # Build protein -> pathway map
    prot_path = pd.DataFrame([
        {"protein_upper": g.upper(), "pathway": p}
        for p, genes in pathways.items() for g in genes
    ]).drop_duplicates()

    merged = df.merge(prot_path, on="protein_upper", how="inner")
    if merged.empty:
        return pd.DataFrame()

    piv = merged.pivot_table(
        index="subject_id", columns="pathway",
        values="abundance_norm", aggfunc="sum", fill_value=0.0
    ).reset_index()
    piv.columns.name = None
    piv.rename(columns={c: c.lower() for c in piv.columns if c != "subject_id"}, inplace=True)

    # Rename to standard MHI component names
    if "fusion" in piv.columns and "fission" in piv.columns:
        piv["fusion_fission"] = piv[["fusion","fission"]].mean(axis=1)
        piv = piv.drop(columns=["fusion","fission"], errors="ignore")
    if "biogenesis" in piv.columns:
        piv = piv.rename(columns={"biogenesis": "copy_number"})

    # Scale to [0,1] per component
    comp_cols = [c for c in ["copy_number","fusion_fission","mitophagy"]
                 if c in piv.columns]
    if not comp_cols:
        return pd.DataFrame()

    for c in comp_cols:
        v = piv[c].values.astype(float)
        lo, hi = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
        piv[c] = (np.clip(v, lo, hi) - lo) / (hi - lo + 1e-9)

    weights = {"copy_number": 0.4, "fusion_fission": 0.35, "mitophagy": 0.25}
    w = pd.Series({c: weights.get(c, 0.33) for c in comp_cols})
    w /= w.sum()
    piv["MHI"] = (piv[comp_cols] * w).sum(axis=1)
    return piv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteomics", required=True)
    ap.add_argument("--genesets",   required=True)
    ap.add_argument("--whitelist",  default=None)
    ap.add_argument("--outdir",     default="results/casestudy_ev/")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] Loading data...")
    ev  = pd.read_csv(args.proteomics)
    gs  = pd.read_csv(args.genesets)
    pathways  = {p: gs[gs["pathway"]==p]["gene"].tolist() for p in gs["pathway"].unique()}
    whitelist = set(pd.read_csv(args.whitelist).iloc[:,0].str.upper()) if args.whitelist else None

    ev["cancer_type"] = ev["subject_id"].apply(parse_cancer_type)
    print(f"[INFO] Cancer types:\n{ev.groupby('cancer_type')['subject_id'].nunique().to_string()}")

    print("[INFO] Computing EV MHI per subject...")
    mhi_df = compute_ev_mhi(ev, pathways, whitelist)
    if mhi_df.empty:
        print("[ERROR] No MHI computed — check pathway gene overlap")
        return

    mhi_df["cancer_type"] = mhi_df["subject_id"].apply(parse_cancer_type)
    mhi_df.to_csv(outdir / "ev_mhi_multicancer.csv", index=False)
    print(f"[OK] MHI computed for {len(mhi_df)} subjects → {outdir}/ev_mhi_multicancer.csv")

    # ── Summary table ──
    summary = mhi_df.groupby("cancer_type")["MHI"].agg(
        n="count", median="median", mean="mean", std="std"
    ).round(3).sort_values("median", ascending=False)
    summary.to_csv(outdir / "ev_mhi_summary.csv")
    print(f"\n── MHI by cancer type ──\n{summary.to_string()}")

    # ── Boxplot ──
    order = summary.index.tolist()
    color_map = {ct: PALETTE[i % len(PALETTE)] for i, ct in enumerate(order)}
    data  = [mhi_df[mhi_df["cancer_type"]==ct]["MHI"].dropna().values for ct in order]
    counts= [len(d) for d in data]

    fig, ax = plt.subplots(figsize=(12, 5))
    bp = ax.boxplot(data, patch_artist=True, widths=0.55, notch=False)
    for patch, ct in zip(bp["boxes"], order):
        patch.set_facecolor(color_map[ct]); patch.set_alpha(0.8)
    ax.set_xticks(range(1, len(order)+1))
    ax.set_xticklabels([f"{ct}\n(n={n})" for ct, n in zip(order, counts)],
                       rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("MHI"); ax.set_title("EV proteomics MHI by cancer type")
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "ev_mhi_boxplot.png", dpi=180)
    plt.close()
    print(f"[OK] → {outdir}/ev_mhi_boxplot.png")

    # ── Component radar / bar ──
    comp_cols = [c for c in ["copy_number","fusion_fission","mitophagy"] if c in mhi_df.columns]
    comp_means = mhi_df.groupby("cancer_type")[comp_cols].mean().loc[order]
    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(order)); w = 0.22
    colors_comp = ["#1D9E75","#3B8BD4","#EF9F27"]
    for i, c in enumerate(comp_cols):
        ax.bar(x + i*w, comp_means[c].values, width=w,
               label=c.replace("_"," ").title(), color=colors_comp[i], alpha=0.85)
    ax.set_xticks(x + w*(len(comp_cols)-1)/2)
    ax.set_xticklabels(order, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Mean score (0–1)")
    ax.set_title("EV MHI components by cancer type")
    ax.legend(fontsize=8); ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "ev_components_by_cancer.png", dpi=180)
    plt.close()
    print(f"[OK] → {outdir}/ev_components_by_cancer.png")

    # ── Top mito proteins per cancer type ──
    ev_mito = ev.copy()
    ev_mito["protein_upper"] = ev_mito["protein"].str.upper()
    gs_prots = set(gs["gene"].str.upper())
    ev_mito = ev_mito[ev_mito["protein_upper"].isin(gs_prots)]
    if not ev_mito.empty:
        prot_pivot = ev_mito.groupby(["cancer_type","protein_upper"])["abundance"].mean().reset_index()
        prot_heat  = prot_pivot.pivot_table(
            index="protein_upper", columns="cancer_type", values="abundance", fill_value=0
        )
        # log scale, keep top 30 proteins by variance
        prot_log  = np.log1p(prot_heat)
        top30     = prot_log.var(axis=1).nlargest(30).index
        prot_plot = prot_log.loc[top30]
        fig, ax = plt.subplots(figsize=(12, 8))
        im = ax.imshow(prot_plot.values, aspect="auto", cmap="YlOrRd")
        ax.set_yticks(range(len(prot_plot.index)));  ax.set_yticklabels(prot_plot.index, fontsize=8)
        ax.set_xticks(range(len(prot_plot.columns))); ax.set_xticklabels(prot_plot.columns, rotation=30, ha="right", fontsize=8)
        ax.set_title("Top 30 mito proteins by variance across cancer types\n(log mean abundance)")
        plt.colorbar(im, ax=ax, label="log abundance", shrink=0.6)
        fig.tight_layout()
        fig.savefig(outdir / "ev_protein_heatmap.png", dpi=180, bbox_inches="tight")
        plt.close()
        print(f"[OK] → {outdir}/ev_protein_heatmap.png")

    # ── Differential stats ──
    types = [ct for ct in order if len(mhi_df[mhi_df["cancer_type"]==ct]) >= 3]
    rows  = []
    from itertools import combinations
    for a, b in combinations(types, 2):
        xa = mhi_df[mhi_df["cancer_type"]==a]["MHI"].dropna().values
        xb = mhi_df[mhi_df["cancer_type"]==b]["MHI"].dropna().values
        if len(xa) < 2 or len(xb) < 2: continue
        U, p = scipy_stats.mannwhitneyu(xa, xb, alternative="two-sided")
        r = 1 - 2*U/(len(xa)*len(xb))
        rows.append({"group_a":a,"group_b":b,"n_a":len(xa),"n_b":len(xb),
                     "median_a":round(float(np.median(xa)),3),
                     "median_b":round(float(np.median(xb)),3),
                     "p_value":float(p),"effect_r":round(float(r),3)})
    if rows:
        st_df = pd.DataFrame(rows).sort_values("p_value")
        st_df.to_csv(outdir / "ev_differential_stats.csv", index=False)
        print(f"\n[OK] Stats → {outdir}/ev_differential_stats.csv")
        print(st_df[["group_a","group_b","n_a","n_b","median_a","median_b","p_value","effect_r"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
