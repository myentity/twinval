"""Rule-based typology and the structural admissibility gate (pure logic, no data files)."""
import numpy as np
import pandas as pd
import pytest

from twinval import structural as S, synthesis as synth
from twinval.config import load_config


def row(r, auc_lo, cb, ccc, level, **kw):
    d = {"scale_id": "x", "scale_name": "x", "instrument_name": "i", "r": r, "AUC": 0.6, "AUC_ci_low": auc_lo, "CCC": ccc, "Cb": cb,
         "pearson_q_bh_sig": True, "structural_evidence_level": level, "paired_esem_admissible": False,
         "tucker_phi_primary": np.nan, "tucker_phi_exploratory": np.nan, "metric_invariance": "no", "scalar_invariance": "no"}
    d.update(kw)
    return d


@pytest.mark.parametrize("args,expected", [
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="metric_invariance"), "A"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="scalar_invariance"), "A"),
    (dict(r=.6, auc_lo=.55, cb=.7, ccc=.4, level="metric_invariance"), "B"),
    (dict(r=.6, auc_lo=.55, cb=.7, ccc=.4, level="structurally_discrepant"), "B"),           # B is evaluated before D
    (dict(r=.3, auc_lo=.55, cb=.9, ccc=.3, level="configural_only"), "C"),
    (dict(r=.3, auc_lo=.55, cb=.9, ccc=.3, level="congruence_supported"), "C"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="structurally_discrepant"), "D"),
    (dict(r=.3, auc_lo=.55, cb=.7, ccc=.45, level="structurally_discrepant"), "D"),          # branch (b): CCC >= .40
    (dict(r=.3, auc_lo=.55, cb=.7, ccc=.30, level="structurally_discrepant"), "F"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="configural_only"), "E"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="model_inadmissible"), "model_inadmissible"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="not_tested"), "not_tested"),
    (dict(r=.5, auc_lo=.501, cb=.8, ccc=.5, level="metric_invariance"), "A"),                # thresholds are inclusive (>=)
    (dict(r=.499, auc_lo=.55, cb=.9, ccc=.5, level="metric_invariance"), "C"),
    (dict(r=.6, auc_lo=.50, cb=.9, ccc=.5, level="metric_invariance"), "C"),                  # AUC CI must be strictly above .50
])
def test_typology_rules(args, expected):
    typ, summ, rules = synth.typology(pd.DataFrame([row(**args)]), load_config())
    assert typ.loc[0, "final_type"] == expected
    assert summ["n"].sum() == 1 and "decision_order" in rules


def test_typology_labels_come_from_config():
    typ, *_ = synth.typology(pd.DataFrame([row(.6, .55, .9, .5, "metric_invariance")]),
                             load_config(overrides={"synthesis": {"type_labels": {"A": "My label"}}}))
    assert typ.loc[0, "final_type_label"] == "My label"


@pytest.mark.parametrize("conf,metric,scalar,applic,paired,phi,expected", [
    (True, "yes", "yes", True, True, .9, "scalar_invariance"),
    (True, "yes", "no", True, True, .9, "metric_invariance"),
    (True, "no", "not_tested", True, True, .3, "configural_only"),               # CFA outranks ESEM/Tucker
    (False, "not_tested", "not_tested", True, True, .90, "congruence_supported"),
    (False, "not_tested", "not_tested", True, True, .80, "structurally_discrepant"),
    (False, "not_tested", "not_tested", True, False, .20, "model_inadmissible"),    # low phi from an inadmissible ESEM is NOT a discrepancy
    (False, "not_tested", "not_tested", False, False, np.nan, "not_tested"),
])
def test_structural_evidence_level(conf, metric, scalar, applic, paired, phi, expected):
    assert S.structural_evidence_level(conf, metric, scalar, applic, paired, phi) == expected


def test_admissibility_requires_both_sides():
    diag = pd.DataFrame({"instrument_id": [1, 1, 2, 2, 3, 3], "side": ["human", "twin"] * 3,
                         "status": ["ok", "ok", "ok", "inadmissible_heywood", "nonconverged", "ok"]})
    g = S.admissibility(diag, [1, 2, 3, 4]).set_index("instrument_id")
    assert list(g["paired_esem_admissible"]) == [True, False, False, False]


def test_cfa_scale_status_conditional_scalar():
    res = pd.DataFrame([
        {"scale_id": "a", "wu_configural_converged": True, "wu_configural_heywood": False, "ms_loading_status": "supported",
         "wu_threshold_status": "supported", "wu_conditional_scalar_status": "supported"},
        {"scale_id": "b", "wu_configural_converged": True, "wu_configural_heywood": False, "ms_loading_status": "rejected",
         "wu_threshold_status": "rejected", "wu_conditional_scalar_status": "supported"},   # thresholds not supported -> scalar not interpretable
        {"scale_id": "c", "wu_configural_converged": True, "wu_configural_heywood": True, "ms_loading_status": "supported",
         "wu_threshold_status": "supported", "wu_conditional_scalar_status": "supported"}])
    s = S.cfa_scale_status(res).set_index("scale_id")
    assert s.loc["a", "metric_invariance"] == "yes" and s.loc["a", "scalar_invariance"] == "yes"
    assert s.loc["b", "metric_invariance"] == "no" and s.loc["b", "scalar_invariance"] == "not_tested"
    assert s.loc["c", "metric_invariance"] == "not_tested" and not s.loc["c", "configural_cfa_admissible"]


def test_tucker_phi_and_bootstrap_merge():
    assert np.isclose(S.tucker_phi(np.array([1., 2, 3]), np.array([2., 4, 6])), 1.0)
    boot = pd.DataFrame({"instrument_id": 1, "factor": ["f"] * 20, "rep": range(20),
                         "phi": [0.9] * 15 + [np.nan] * 5, "valid": [True] * 15 + [False] * 5})
    by_factor, mean_ci = S.merge_tucker_bootstrap(pd.DataFrame(), boot, 0.80)
    assert by_factor.loc[0, "frac_valid_bootstrap"] == 0.75 and not by_factor.loc[0, "ci_stable"]


def test_order_sensitivity_returns_empty_when_all_models_are_clean(tmp_path):
    res = pd.DataFrame({"scale_id": ["a"], "wu_configural_converged": [True], "wu_configural_heywood": [False],
                        "wu_threshold_status": ["supported"]})
    out = S.cfa_order_sensitivity(res, pd.DataFrame(), tmp_path, 3, 1, load_config(), tmp_path)
    assert out.empty
