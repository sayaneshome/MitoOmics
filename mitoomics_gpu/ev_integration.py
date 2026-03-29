from __future__ import annotations
import pandas as pd


def ev_pathway_scores(
    ev_df: pd.DataFrame,
    pathways: dict[str, list[str]],
    whitelist: set[str] | None = None,
) -> pd.DataFrame:
    """
    Aggregate EV/MDV proteomics to subject-level pathway scores.

    - Optionally filter to whitelist proteins (uppercased symbols)
    - Map proteins -> pathways via MitoCarta-derived pathway dict
    - Sum normalized abundance per (subject, pathway)

    Fully vectorized — no iterrows over the full dataset.
    """
    if ev_df is None or ev_df.empty:
        return pd.DataFrame(columns=["subject_id"])

    df = ev_df.copy()
    df["protein_upper"] = df["protein"].astype(str).str.upper()

    if whitelist:
        wl_upper = {p.upper() for p in whitelist}
        df = df[df["protein_upper"].isin(wl_upper)]
        if df.empty:
            return pd.DataFrame({"subject_id": ev_df["subject_id"].unique()})

    # Build protein -> pathway mapping as a DataFrame for a vectorized merge
    prot_path_rows = [
        {"protein_upper": str(sym).upper(), "pathway": pth}
        for pth, syms in pathways.items()
        for sym in syms
    ]
    if not prot_path_rows:
        return pd.DataFrame({"subject_id": df["subject_id"].unique()})

    prot_path = pd.DataFrame(prot_path_rows).drop_duplicates()

    # Vectorized: merge abundance rows with pathway membership, then pivot
    merged = df.merge(prot_path, on="protein_upper", how="inner")

    if merged.empty:
        return pd.DataFrame({"subject_id": df["subject_id"].unique()})

    piv = merged.pivot_table(
        index="subject_id",
        columns="pathway",
        values="abundance_norm",
        aggfunc="sum",
        fill_value=0.0,
    ).reset_index()

    piv.columns.name = None
    piv.rename(
        columns={c: c.lower() for c in piv.columns if c != "subject_id"},
        inplace=True,
    )
    return piv
