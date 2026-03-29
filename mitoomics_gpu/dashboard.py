"""
dashboard.py
Run with:  streamlit run -m mitoomics_gpu.dashboard -- --results results/results_summary.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ── page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="MitoOmics-GPU · MHI Dashboard",
    page_icon="🧬",
    layout="wide",
)


# ── helpers ──────────────────────────────────────────────────────────────────
@st.cache_data
def load(path: str) -> pd.DataFrame:
    return pd.read_csv(path)

COMPONENTS = ["copy_number", "fusion_fission", "mitophagy", "heterogeneity"]
COMP_LABELS = {
    "copy_number":   "mtDNA copy number",
    "fusion_fission":"Fusion / Fission",
    "mitophagy":     "Mitophagy",
    "heterogeneity": "Cell heterogeneity",
}
PALETTE = ["#1D9E75", "#3B8BD4", "#EF9F27", "#D85A30", "#7F77DD", "#D4537E"]


def radar_chart(values: dict[str, float], subject: str) -> plt.Figure:
    cats = list(values.keys())
    N = len(cats)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]
    vals = list(values.values()) + [list(values.values())[0]]

    fig, ax = plt.subplots(figsize=(3, 3), subplot_kw=dict(polar=True))
    ax.plot(angles, vals, "o-", linewidth=1.5, color="#1D9E75")
    ax.fill(angles, vals, alpha=0.2, color="#1D9E75")
    ax.set_thetagrids(np.degrees(angles[:-1]), cats, fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_title(subject, size=9, pad=8)
    ax.set_yticklabels([])
    ax.spines["polar"].set_visible(False)
    fig.tight_layout()
    return fig


# ── sidebar ───────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--results", default="results/results_summary.csv")
parser.add_argument("--stats",   default=None, help="mhi_differential_stats.csv (optional)")
args, _ = parser.parse_known_args()

results_path = st.sidebar.text_input("Results CSV", value=args.results)
stats_path   = st.sidebar.text_input("Stats CSV (optional)", value=args.stats or "")

if not Path(results_path).exists():
    st.error(f"File not found: {results_path}")
    st.stop()

df = load(results_path)
avail_comps = [c for c in COMPONENTS if c in df.columns]

st.sidebar.markdown("---")
st.sidebar.markdown("**Filter subjects**")
top_n = st.sidebar.slider("Show top N by MHI", 5, min(100, len(df)), min(30, len(df)))
group_col = st.sidebar.selectbox(
    "Colour by group",
    ["(none)"] + [c for c in df.columns if c not in ["subject_id", "MHI"] + COMPONENTS],
)

# ── header ────────────────────────────────────────────────────────────────────
st.title("🧬 MitoOmics-GPU — Mitochondrial Health Index")
st.caption("Research use only · Not for clinical diagnosis")

# ── KPI row ───────────────────────────────────────────────────────────────────
c1, c2, c3, c4 = st.columns(4)
c1.metric("Subjects", len(df))
c2.metric("Mean MHI", f"{df['MHI'].mean():.3f}")
c3.metric("Max MHI",  f"{df['MHI'].max():.3f}")
c4.metric("Components", len(avail_comps))

st.divider()

# ── MHI bar chart ─────────────────────────────────────────────────────────────
col_a, col_b = st.columns([2, 1])

with col_a:
    st.subheader(f"Top {top_n} subjects by MHI")
    plot_df = df.sort_values("MHI", ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=(10, 3.5))
    colors = PALETTE[0]
    if group_col != "(none)" and group_col in plot_df.columns:
        unique_groups = plot_df[group_col].unique()
        color_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(unique_groups)}
        colors = [color_map[g] for g in plot_df[group_col]]
        handles = [plt.Rectangle((0,0),1,1, color=color_map[g]) for g in unique_groups]
        ax.legend(handles, unique_groups, fontsize=8, loc="upper right")

    ax.bar(range(len(plot_df)), plot_df["MHI"].values, color=colors, width=0.7)
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(plot_df["subject_id"].astype(str), rotation=90, fontsize=7)
    ax.set_ylabel("MHI")
    ax.set_ylim(0, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    st.pyplot(fig)
    plt.close()

with col_b:
    st.subheader("MHI distribution")
    fig2, ax2 = plt.subplots(figsize=(4, 3.5))
    ax2.hist(df["MHI"].dropna(), bins=20, color="#1D9E75", edgecolor="white")
    ax2.set_xlabel("MHI")
    ax2.set_ylabel("Subjects")
    ax2.spines[["top", "right"]].set_visible(False)
    fig2.tight_layout()
    st.pyplot(fig2)
    plt.close()

st.divider()

# ── Component scatter grid ────────────────────────────────────────────────────
if avail_comps:
    st.subheader("Component vs MHI")
    cols = st.columns(len(avail_comps))
    for i, comp in enumerate(avail_comps):
        with cols[i]:
            fig3, ax3 = plt.subplots(figsize=(3, 3))
            ax3.scatter(df[comp], df["MHI"], s=12, alpha=0.7, color=PALETTE[i % len(PALETTE)])
            ax3.set_xlabel(COMP_LABELS.get(comp, comp), fontsize=8)
            ax3.set_ylabel("MHI" if i == 0 else "", fontsize=8)
            ax3.spines[["top", "right"]].set_visible(False)
            corr = df[[comp, "MHI"]].dropna().corr().iloc[0, 1]
            ax3.set_title(f"r = {corr:.2f}", fontsize=8)
            fig3.tight_layout()
            st.pyplot(fig3)
            plt.close()

st.divider()

# ── Subject deep-dive ─────────────────────────────────────────────────────────
st.subheader("Subject deep-dive")
selected = st.selectbox("Select subject", df["subject_id"].astype(str).tolist())
row = df[df["subject_id"].astype(str) == selected].iloc[0]

d1, d2 = st.columns([1, 2])
with d1:
    comp_vals = {COMP_LABELS.get(c, c): float(row[c]) for c in avail_comps if c in row and pd.notna(row[c])}
    if comp_vals:
        fig4 = radar_chart(comp_vals, selected)
        st.pyplot(fig4)
        plt.close()

with d2:
    st.markdown(f"**MHI: {row['MHI']:.4f}**")
    detail = {COMP_LABELS.get(c, c): f"{row[c]:.4f}" for c in avail_comps if c in row}
    st.table(pd.DataFrame(detail, index=["score"]).T.rename(columns={"score": "Value"}))

st.divider()

# ── Differential stats (optional) ─────────────────────────────────────────────
if stats_path and Path(stats_path).exists():
    st.subheader("Differential MHI — group comparison")
    sdf = load(stats_path)
    sig = sdf[sdf["p_adj"] < 0.05] if "p_adj" in sdf.columns else sdf
    st.markdown(f"**{len(sig)}** significant pairwise comparisons (BH-adjusted p < 0.05)")

    def highlight_sig(row):
        color = "background-color: #E1F5EE" if row.get("significant", False) else ""
        return [color] * len(row)

    st.dataframe(
        sdf.style.apply(highlight_sig, axis=1).format(
            {"p_value": "{:.4f}", "p_adj": "{:.4f}", "effect_r": "{:.3f}",
             "median_a": "{:.3f}", "median_b": "{:.3f}"}
        ),
        use_container_width=True,
    )

st.divider()

# ── Full results table ────────────────────────────────────────────────────────
with st.expander("Full results table"):
    st.dataframe(df.sort_values("MHI", ascending=False), use_container_width=True)
