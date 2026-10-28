"""Regression tests: the framework reproduces the published results for the two example instruments.

The reference tables (tests/reference/paper_*.csv) hold aggregate values of the original study (N = 484).
Python-only tests run everywhere; tests that fit lavaan models need Rscript and are marked ``r``.
"""
import numpy as np
import pandas as pd
import pytest

from twinval import item, prep, scale, structural as S
from conftest import REF


def cmp(a, b, tol, what):
    m = a.merge(b, on="scale_id" if "scale_id" in a else "item_id", suffixes=("", "_ref"))
    assert len(m) == len(b), f"{what}: rows missing"
    return m


def test_scale_level_reproduces_paper(study, cfg):
    ref = pd.read_csv(REF / "paper_scale_level.csv")
    core = scale.core_table(study, cfg, np.random.default_rng(0))
    m = core.merge(ref, on="scale_id", suffixes=("", "_ref"))
    assert len(m) == 24
    assert (m["pearson_r"] - m["pearson_r_ref"]).abs().max() < 1e-3
    assert (m["ccc"] - m["ccc_ref"]).abs().max() < 1e-3
    assert (m["cb"] - m["cb_ref"]).abs().max() < 1e-3
    assert (m["auc"] - m["auc_ref"]).abs().max() < 0.01          # published scores were rounded to 3 decimals (ties in AUC)
    assert (m["typology_r_cb"] == m["typology_r_cb_ref"]).all()


def test_recalibration_and_disattenuation_reproduce_paper(study, cfg):
    ref = pd.read_csv(REF / "paper_scale_level.csv")
    rng = np.random.default_rng(0)
    core = scale.core_table(study, cfg, rng)
    _, per_scale = scale.recalibration(study, cfg, rng)
    m = per_scale.merge(ref, on="scale_id", suffixes=("", "_ref"))
    assert (m["delta_ccc_ls"] - m["delta_ccc_ls_ref"]).abs().max() < 5e-3
    assert (m["delta_ccc_ols"] - m["delta_ccc_ols_ref"]).abs().max() < 5e-3
    dis = scale.disattenuation(study, cfg, rng, core).merge(ref, on="scale_id", suffixes=("", "_ref"))
    ok = dis.dropna(subset=["r_dis_alpha"])
    assert (ok["r_dis"] - ok["r_dis_alpha"]).abs().max() < 1e-3


def test_item_level_reproduces_paper_with_original_categories(study, cfg):
    ref = pd.read_csv(REF / "paper_item_level.csv")
    ag, bias, _ = item.agreement_and_bias(item.build_panel(study), cfg, np.random.default_rng(0))
    m = ag.merge(bias[["item_id", "sd_ratio", "bias_range", "entropy_difference"]], on="item_id").merge(ref, on="item_id", suffixes=("", "_ref"))
    assert len(m) == 75
    for col in ("exact_agreement", "mae_range", "spearman_rho", "ccc", "sd_ratio", "bias_range", "entropy_difference"):
        assert (m[col] - m[col + "_ref"]).abs().max() < 1e-6, col


def test_variance_decomposition_reproduces_paper(study, cfg):
    ref = pd.read_csv(REF / "paper_variance_decomposition.csv")
    dec, _ = item.variance_decomposition(item.build_panel(study), cfg, np.random.default_rng(0))
    m = dec.merge(ref, on=["instrument_id", "side"], suffixes=("", "_ref"))
    assert len(m) == 6
    for col in ("pct_person", "pct_item", "pct_person_item_residual"):
        assert (m[col] - m[col + "_ref"]).abs().max() < 1e-6


def test_harmonisation_status_matches_paper(study):
    hh, ht, log, bad = prep.harmonize(study)
    st = prep.harmonization_status(study, log, bad).set_index("scale_id")["category_status"]
    assert st["148_Extraversion"] == "harmonized" and st["148_Depression"] == "clean" and len(bad) == 0


@pytest.mark.r
def test_structural_chain_reproduces_paper(study, cfg, tmp_path):
    from twinval.structural import r_available
    ok, msg = r_available()
    if not ok:
        pytest.skip(msg)
    wd = tmp_path / "w"
    inp = S.prepare_inputs(study, cfg, wd)
    S.run_r("esem_fit.R", [wd], wd / "esem.log")
    diag = pd.read_csv(wd / "esem_model_diagnostics.csv")
    pts = S.tucker_point_estimates(wd, diag, inp["esem_instruments"])["tucker_by_factor"]
    ref = pd.read_csv(REF / "paper_tucker_phi.csv")
    m = pts.rename(columns={"factor": "scale_id"}).merge(ref, on="scale_id", suffixes=("", "_ref"))
    assert len(m) == 8 and (m["phi_full_column"] - m["phi_full_column_ref"]).abs().max() < 1e-5
    assert S.admissibility(diag, [148, 161])["paired_esem_admissible"].all()


@pytest.mark.r
def test_full_pipeline_types_and_cfa_reproduce_paper(tmp_path):
    """CFA verdicts and the final type of every example scale agree with the article (needs R; ~2 min)."""
    from twinval.config import load_config
    from twinval.pipeline import Pipeline
    from twinval.structural import r_available
    from conftest import EXAMPLE, FAST
    ok, msg = r_available()
    if not ok:
        pytest.skip(msg)
    cfg = load_config(overrides={**FAST, "data": {"data_dir": str(EXAMPLE)}, "structural": {"tucker_bootstrap": False}})
    Pipeline(cfg, tmp_path / "o").run(["structural", "scale", "item", "synthesis"])
    out = tmp_path / "o"
    cfa = pd.read_csv(out / "03_structural" / "invariance_decisions.csv").merge(
        pd.read_csv(REF / "paper_cfa_status.csv"), on="scale_id", suffixes=("", "_ref"))
    assert len(cfa) == 24
    for col in ("wu_threshold_status", "ms_loading_status", "wu_conditional_scalar_status"):
        assert (cfa[col].astype(str) == cfa[col + "_ref"].astype(str)).all(), col
    typ = pd.read_csv(out / "05_synthesis" / "scale_typology.csv").merge(pd.read_csv(REF / "paper_typology.csv"), on="scale_id", suffixes=("", "_ref"))
    assert (typ["final_type"] == typ["final_type_ref"]).all()
    assert (typ["structural_evidence_level"] == typ["structural_evidence_level_ref"]).all()
    assert set(typ["final_type"]) >= {"A", "B", "C", "D", "E"}          # the example shows all five substantive types
