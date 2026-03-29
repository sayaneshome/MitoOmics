"""
cli.py — MitoOmics-GPU end-to-end pipeline.

Handles two analysis modes automatically:
  - scRNA-only subjects  : scored from copy_number + pathway programs + heterogeneity
  - EV-only subjects     : scored from EV proteomics pathway scores
  - Overlapping subjects : full integration (if subject IDs match)
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

from .io import (
    load_scrna, load_proteomics, load_imaging,
    load_mitocarta_pathways, load_ev_whitelist,
)
from .scoring import copy_number_proxy, program_scores, heterogeneity_scores
from .ev_integration import ev_pathway_scores
from .mhi import combine_components
from .report import save_figures, write_report
from .stats import mhi_group_test, append_stats_to_report


def main():
    ap = argparse.ArgumentParser(description="MitoOmics-GPU pipeline")
    ap.add_argument("--scrna",           required=True)
    ap.add_argument("--proteomics",      required=True)
    ap.add_argument("--imaging",         default=None)
    ap.add_argument("--outdir",          required=True)
    ap.add_argument("--mitocarta-table", required=True)
    ap.add_argument("--ev-whitelist",    default=None)
    ap.add_argument("--fit-ridge",       action="store_true")
    ap.add_argument("--targets",         default=None)
    ap.add_argument("--genesets",        default=None,
                    help="genesets_curated.csv (pathway,gene) — used as primary scorer. Falls back to MitoCarta if not provided.")
    ap.add_argument("--group-col",       default=None)
    ap.add_argument("--use-gpu",         action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    adata    = load_scrna(args.scrna)
    ev_df    = load_proteomics(args.proteomics)
    # Use curated gene sets as primary (fusion/fission/mitophagy/biogenesis)
    # Fall back to full MitoCarta if not provided
    if args.genesets and __import__('pathlib').Path(args.genesets).exists():
        import pandas as _pd
        gs = _pd.read_csv(args.genesets)
        pathways = {p: gs[gs["pathway"]==p]["gene"].tolist()
                    for p in gs["pathway"].unique()}
        print(f"[INFO] Using curated gene sets: {list(pathways.keys())}")
    else:
        pathways = load_mitocarta_pathways(args.mitocarta_table)
    if not pathways:
        raise ValueError("No pathways found.")
    whitelist = load_ev_whitelist(args.ev_whitelist) if args.ev_whitelist else None

    scrna_subjects = set(adata.obs["subject_id"].unique())
    ev_subjects    = set(ev_df["subject_id"].unique())
    overlap        = scrna_subjects & ev_subjects

    print(f"[INFO] scRNA subjects : {len(scrna_subjects)}")
    print(f"[INFO] EV subjects    : {len(ev_subjects)}")
    print(f"[INFO] Overlap        : {len(overlap)}")

    # ── scRNA components (for scRNA subjects) ────────────────────────────────
    cnp_cell = copy_number_proxy(adata)
    adata.obs["copy_number_proxy"] = cnp_cell.values
    cnp_subject = (
        adata.obs.groupby("subject_id", observed=True)["copy_number_proxy"]
        .mean().rename("copy_number")
    )

    prog_cell = program_scores(adata, pathways=pathways, use_gpu=args.use_gpu)
    adata.obs = adata.obs.join(prog_cell)

    def _pick(df, kws):
        return [c for c in df.columns if any(k in c for k in kws)]

    ff_cols    = _pick(prog_cell, ["fusion", "fission"])
    mito_cols  = _pick(prog_cell, ["mitophagy"])
    bio_cols   = _pick(prog_cell, ["biogenesis"])

    ff_subject = (
        adata.obs.groupby("subject_id", observed=True)[ff_cols]
        .mean().mean(axis=1).rename("fusion_fission")
    ) if ff_cols else pd.Series(dtype=float, name="fusion_fission")

    mitophagy_subject = (
        adata.obs.groupby("subject_id", observed=True)[mito_cols]
        .mean().mean(axis=1).rename("mitophagy")
    ) if mito_cols else pd.Series(dtype=float, name="mitophagy")

    biogenesis_subject = (
        adata.obs.groupby("subject_id", observed=True)[bio_cols]
        .mean().mean(axis=1).rename("biogenesis")
    ) if bio_cols else pd.Series(dtype=float, name="biogenesis")

    het_subject = heterogeneity_scores(adata).rename("heterogeneity")

    # ── EV proteomics components ──────────────────────────────────────────────
    ev_path = ev_pathway_scores(ev_df, pathways=pathways, whitelist=whitelist)

    # ── Attach subject-level metadata from adata.obs ─────────────────────────
    meta_cols = ["disease", "sex", "tissue_general", "cell_type"]
    meta_cols = [c for c in meta_cols if c in adata.obs.columns]
    subject_meta = pd.DataFrame({"subject_id": list(scrna_subjects)})
    for col in meta_cols:
        agg = (
            adata.obs.groupby("subject_id", observed=True)[col]
            .agg(lambda x: x.value_counts().index[0])
            .reset_index()
            .rename(columns={col: col})
        )
        subject_meta = subject_meta.merge(agg, on="subject_id", how="left")

    # ── Build scRNA subject DataFrame ─────────────────────────────────────────
    scrna_df = pd.DataFrame({"subject_id": sorted(scrna_subjects)})
    for s in [cnp_subject, ff_subject, mitophagy_subject, biogenesis_subject, het_subject]:
        if s is not None and len(s) > 0:
            scrna_df = scrna_df.merge(s.to_frame(), left_on="subject_id",
                                      right_index=True, how="left")
    scrna_df = scrna_df.merge(subject_meta, on="subject_id", how="left")

    # ── Build EV subject DataFrame ────────────────────────────────────────────
    if not ev_path.empty:
        ev_df2 = ev_path.copy()
    else:
        ev_df2 = pd.DataFrame({"subject_id": sorted(ev_subjects)})

    # ── Compute MHI separately for each cohort ────────────────────────────────
    ridge_targets = None
    if args.targets:
        tgt = pd.read_csv(args.targets)
        ridge_targets = tgt.set_index("subject_id")["target"]

    # scRNA MHI
    scrna_mhi = combine_components(
        scrna_df,
        fit_ridge=args.fit_ridge and ridge_targets is not None,
        targets=ridge_targets,
    )
    scrna_mhi["cohort"] = "scRNA"

    # EV MHI — rename pathway cols to standard component names where possible
    ev_rename = {"fusion": "fusion_fission", "fission": "fusion_fission",
                 "biogenesis": "copy_number"}
    if not ev_df2.empty and "subject_id" in ev_df2.columns:
        ev_agg = ev_df2.copy()
        # Merge fusion + fission into fusion_fission
        if "fusion" in ev_agg.columns and "fission" in ev_agg.columns:
            ev_agg["fusion_fission"] = ev_agg[["fusion","fission"]].mean(axis=1)
            ev_agg = ev_agg.drop(columns=["fusion","fission"], errors="ignore")
        elif "fusion" in ev_agg.columns:
            ev_agg = ev_agg.rename(columns={"fusion": "fusion_fission"})
        elif "fission" in ev_agg.columns:
            ev_agg = ev_agg.rename(columns={"fission": "fusion_fission"})
        if "biogenesis" in ev_agg.columns:
            ev_agg = ev_agg.rename(columns={"biogenesis": "copy_number"})
    else:
        ev_agg = ev_df2

    ev_comp_cols = [c for c in ["copy_number","fusion_fission","mitophagy","heterogeneity"]
                    if c in ev_agg.columns and ev_agg[c].notna().any()]
    if ev_comp_cols:
        try:
            ev_mhi = combine_components(ev_agg)
            ev_mhi["cohort"] = "EV_proteomics"
        except Exception as e:
            print(f"[WARN] EV MHI skipped: {e}")
            ev_mhi = pd.DataFrame()
    else:
        ev_mhi = pd.DataFrame()

    # ── Save results ──────────────────────────────────────────────────────────
    scrna_mhi.to_csv(outdir / "results_scrna.csv", index=False)
    print(f"[OK] scRNA results  → {outdir}/results_scrna.csv  ({len(scrna_mhi)} subjects)")

    if not ev_mhi.empty:
        ev_mhi.to_csv(outdir / "results_ev.csv", index=False)
        print(f"[OK] EV results     → {outdir}/results_ev.csv  ({len(ev_mhi)} subjects)")

    # Combined for downstream tools
    all_results = pd.concat([scrna_mhi, ev_mhi], ignore_index=True) if not ev_mhi.empty else scrna_mhi
    all_results.to_csv(outdir / "results_summary.csv", index=False)
    print(f"[OK] Combined       → {outdir}/results_summary.csv  ({len(all_results)} subjects)")

    # ── Figures + report ──────────────────────────────────────────────────────
    figs = save_figures(scrna_mhi, outdir)
    write_report(scrna_mhi, figs, outdir)

    # ── Differential stats ────────────────────────────────────────────────────
    group_col = args.group_col
    if group_col:
        # Try scRNA results first (has disease labels)
        target_df = scrna_mhi if group_col in scrna_mhi.columns else all_results
        if group_col in target_df.columns:
            stats_df = mhi_group_test(target_df, group_col=group_col)
            if not stats_df.empty:
                stats_df.to_csv(outdir / "mhi_differential_stats.csv", index=False)
                append_stats_to_report(stats_df, str(outdir / "report.md"), group_col)
                print(f"[OK] Stats          → {outdir}/mhi_differential_stats.csv")
                print(f"\n{'─'*80}")
                print(f"Differential MHI by {group_col}:")
                print(f"{'─'*80}")
                print(stats_df[["group_a","group_b","n_a","n_b",
                                 "median_a","median_b","p_value","p_adj","effect_r"]]
                      .to_string(index=False))
            else:
                print(f"[WARN] Not enough subjects per group for stats.")
        else:
            print(f"[WARN] '{group_col}' not found. Available: {list(scrna_mhi.columns)}")

if __name__ == "__main__":
    main()
