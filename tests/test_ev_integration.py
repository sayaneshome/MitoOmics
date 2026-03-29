"""Tests for ev_integration.py — vectorized pathway aggregation."""
import pandas as pd
import pytest

from mitoomics_gpu.ev_integration import ev_pathway_scores


def _make_ev(subjects=("S1","S2","S3"), proteins=("MFN1","DNM1L","PINK1","OTHER")):
    rows = []
    for s in subjects:
        for p in proteins:
            rows.append({"subject_id": s, "protein": p, "abundance": 1.0, "abundance_norm": 0.25})
    return pd.DataFrame(rows)


PATHWAYS = {
    "fusion":    ["MFN1", "MFN2", "OPA1"],
    "fission":   ["DNM1L", "FIS1"],
    "mitophagy": ["PINK1", "PRKN"],
}


def test_basic_pivot_shape():
    ev = _make_ev()
    out = ev_pathway_scores(ev, PATHWAYS)
    assert "subject_id" in out.columns
    assert len(out) == 3  # 3 subjects


def test_whitelist_filters_proteins():
    ev = _make_ev()
    # Only keep mitophagy proteins in whitelist
    out = ev_pathway_scores(ev, PATHWAYS, whitelist={"PINK1", "PRKN"})
    assert "mitophagy" in out.columns
    # fusion and fission proteins not in whitelist → columns absent or zero
    if "fusion" in out.columns:
        assert (out["fusion"] == 0).all()


def test_empty_ev_returns_empty():
    ev = pd.DataFrame(columns=["subject_id","protein","abundance","abundance_norm"])
    out = ev_pathway_scores(ev, PATHWAYS)
    assert "subject_id" in out.columns
    assert len(out) == 0


def test_no_pathway_overlap_returns_subjects():
    ev = _make_ev(proteins=("NOCALL",))
    out = ev_pathway_scores(ev, PATHWAYS)
    # Should still have subject_id column, no pathway columns
    assert "subject_id" in out.columns
    assert set(out.columns) == {"subject_id"}


def test_uppercase_insensitive():
    ev = _make_ev(proteins=("mfn1", "dnm1l"))
    out = ev_pathway_scores(ev, PATHWAYS)
    assert "fusion" in out.columns or "fission" in out.columns
