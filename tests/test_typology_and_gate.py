"""Rule-based typology and the structural admissibility gate (pure logic, no data files)."""
import numpy as np
import pandas as pd
import pytest

from twinval import structural as S, synthesis as synth
from twinval.config import load_config

V1 = load_config(overrides={"structural": {"rules_version": 1}})


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
def test_typology_rules_v1(args, expected):
    typ, summ, rules = synth.typology(pd.DataFrame([row(**args)]), V1)
    assert typ.loc[0, "final_type"] == expected
    assert summ["n"].sum() == 1 and "decision_order" in rules


def test_typology_labels_come_from_config():
    typ, *_ = synth.typology(pd.DataFrame([row(.6, .55, .9, .5, "metric_invariance")]),
                             load_config(overrides={"synthesis": {"type_labels": {"A": "My label"}}}))
    assert typ.loc[0, "final_type_label"] == "My label"
    typ, *_ = synth.typology(pd.DataFrame([row(.6, .55, .9, .5, "metric_invariance")]),
                             load_config(overrides={"structural": {"rules_version": 1},
                                                    "synthesis": {"type_labels_v1": {"A": "My label"}}}))
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


def test_v1_order_sensitivity_returns_empty_when_all_models_are_clean(tmp_path):
    res = pd.DataFrame({"scale_id": ["a"], "wu_configural_converged": [True], "wu_configural_heywood": [False],
                        "wu_threshold_status": ["supported"]})
    out = S.cfa_order_sensitivity(res, pd.DataFrame(), tmp_path, 3, 1, load_config(overrides={"structural": {"rules_version": 1}}), tmp_path)
    assert out.empty


# =============================================================================== rules_version 2
V2 = load_config()


@pytest.mark.parametrize("args,expected", [
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="metric_invariance"), "A"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="scalar_invariance"), "A"),
    (dict(r=.6, auc_lo=.55, cb=.7, ccc=.4, level="metric_invariance"), "B"),
    (dict(r=.6, auc_lo=.55, cb=.7, ccc=.4, level="structurally_discrepant"), "B"),           # B is evaluated before D
    (dict(r=.3, auc_lo=.55, cb=.9, ccc=.3, level="configural_fit"), "C"),
    (dict(r=.3, auc_lo=.55, cb=.9, ccc=.3, level="congruence_supported"), "C"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="structurally_discrepant"), "D"),
    (dict(r=.3, auc_lo=.55, cb=.7, ccc=.45, level="structurally_discrepant"), "F"),          # no CCC >= .40 branch any more
    (dict(r=.3, auc_lo=.55, cb=.9, ccc=.45, level="structurally_discrepant"), "F"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="configural_fit"), "E"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="cfa_misfit"), "cfa_misfit"),               # analysis statuses are not types
    (dict(r=.3, auc_lo=.55, cb=.9, ccc=.3, level="cfa_fit_not_testable"), "cfa_fit_not_testable"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="model_not_estimable"), "model_not_estimable"),
    (dict(r=.6, auc_lo=.55, cb=.9, ccc=.5, level="not_tested"), "not_tested"),
])
def test_typology_rules_v2(args, expected):
    typ, summ, rules = synth.typology(pd.DataFrame([row(**args)]), V2)
    assert typ.loc[0, "final_type"] == expected
    assert rules["rules_version"] == 2 and summ["n"].sum() == 1


def test_evidence_quality_reports_flags_without_changing_the_type():
    r = row(.6, .55, .9, .5, "metric_invariance",
            structural_flags="near_boundary_loading;identification_scheme_sensitive;baseline_fit_not_testable")
    typ, *_ = synth.typology(pd.DataFrame([r]), V2)
    assert typ.loc[0, "final_type"] == "A"
    q = typ.loc[0, "evidence_quality"]
    assert "baseline fit not testable" in q and "identification-sensitive" in q and "near-boundary" in q


def model(prefix, **kw):
    """Diagnostics of one clean admissible model in the column layout of paired_cfa_invariance.R."""
    d = dict(converged=1, df=15, heywood=0, vcov_min_eigen=1e-3, theta_min_eigen=0.2, psi_min_eigen=0.3,
             max_std_loading=0.8, max_abs_resid_cor=0.3, abs_factor_cor=0.6, srmr=0.04, cfi=0.99, rmsea_ci_lo=0.02)
    d.update(kw)
    return {f"{prefix}_{k}": v for k, v in d.items()}


@pytest.mark.parametrize("kw,status,flags", [
    (dict(), "admissible", []),
    (dict(converged=0), "nonconverged", []),
    (dict(df=-1), "identification_not_established", []),
    (dict(vcov_min_eigen=np.nan), "identification_not_established", []),        # vcov not computable
    (dict(heywood=1), "inadmissible", []),
    (dict(theta_min_eigen=-0.5), "inadmissible", []),                           # residual matrix not PD
    (dict(theta_min_eigen=-1e-9), "admissible", []),                            # numerical deviation within tolerance
    (dict(max_abs_resid_cor=1.2), "inadmissible", []),
    (dict(abs_factor_cor=1.01), "inadmissible", []),
    (dict(max_std_loading=0.97), "admissible", ["near_boundary_loading"]),      # warning, not exclusion
    (dict(max_abs_resid_cor=0.96), "admissible", ["near_boundary_residual_correlation"]),
])
def test_model_status(kw, status, flags):
    st, fl = S.model_status(model("m", **kw), "m", V2.thresholds)
    assert st == status and fl == flags


def test_model_status_not_fitted():
    assert S.model_status({}, "m", V2.thresholds) == ("not_fitted", [])


@pytest.mark.parametrize("kw,fit,flags", [
    (dict(), "supported", []),
    (dict(srmr=0.09), "not_supported", []),
    (dict(df=0), "not_testable", ["baseline_fit_not_testable"]),
    (dict(cfi=0.85), "supported", ["conflicting_fit_diagnostics"]),
    (dict(rmsea_ci_lo=0.12), "supported", ["conflicting_fit_diagnostics"]),
])
def test_baseline_fit(kw, fit, flags):
    assert S.baseline_fit(model("m", **kw), "m", V2.thresholds) == (fit, flags)


def wu_row(w0=None, w1=None, w2=None, w3=None, d1=(0, 0.01, 0), d2=(0, 0.01, 0.01), d3=(0, 0.01, 0.005), ms=None, dm=(0, 0.01, 0.01)):
    """One scale of invariance_results.csv: Wu W0-W3 and Millsap M0-M1 with per-step deltas (CFI, RMSEA, SRMR)."""
    r = {"scale_id": "s"}
    for pref, kw, df in (("wu_configural", w0, 15), ("wu_thresholds", w1, 20), ("wu_conditional_loading", w2, 23),
                         ("wu_conditional_scalar", w3, 26), ("ms_baseline", w0, 15), ("ms_loading", ms, 18)):
        if kw is not False:
            r.update(model(pref, **{"df": df, **(kw or {})}))
    for dec, d in (("wu_threshold", d1), ("wu_conditional_loading", d2), ("wu_conditional_scalar", d3), ("ms_loading", dm)):
        r.update({f"{dec}_delta_cfi": d[0], f"{dec}_delta_rmsea": d[1], f"{dec}_delta_srmr": d[2]})
    return r


def status(r, **struct):
    cfg = load_config(overrides={"structural": struct}) if struct else V2
    return S.cfa_scale_status_v2(pd.DataFrame([r]), cfg).iloc[0]


def test_wu_sequence_scalar():
    s = status(wu_row())
    assert (s.threshold_decision, s.metric_decision, s.scalar_decision) == ("supported", "supported", "supported")
    assert s.cfa_conclusion == "scalar" and "threshold_supported" in s.cfa_flags


def test_wu_threshold_rejected_blocks_later_steps_and_millsap_does_not_replace_it():
    s = status(wu_row(d1=(-0.05, 0.01, 0), w2=False, w3=False))
    assert s.threshold_decision == "rejected" and s.metric_decision == "not_interpreted_prior_not_supported"
    assert s.cfa_conclusion == "configural_fit"                                    # Millsap M1 supported, but not used
    # the Wu loading step is not interpreted, so the tests do not contradict each other: route dependence, not a conflict
    assert "identification_route_dependent" in s.cfa_flags and "identification_scheme_sensitive" not in s.cfa_flags


def test_identification_flags_need_comparable_results():
    s = status(wu_row(dm=(-0.05, 0.03, 0.05)))                                     # Wu W2 supported, Millsap M1 rejected
    assert s.metric_decision == "supported" and "identification_scheme_sensitive" in s.cfa_flags
    s = status(wu_row(w0=dict(srmr=0.12)))                                         # baseline misfit, both loading tests supported
    assert s.cfa_conclusion == "none" and s.sensitivity_millsap_metric_decision == "supported"
    assert not any(f.startswith("identification_") for f in s.cfa_flags.split(";"))
    s = status(wu_row(w0=dict(srmr=0.12), d1=(-0.05, 0.01, 0), w2=False, w3=False))  # misfit + rejected thresholds
    assert not any(f.startswith("identification_") for f in s.cfa_flags.split(";"))


def test_delta_srmr_is_part_of_the_decision():
    r = wu_row(d2=(0, 0.01, 0.05))                                                 # loadings: dSRMR .05 > .030
    assert status(r).metric_decision == "rejected"
    assert status(r, use_delta_srmr=False).metric_decision == "supported"         # sensitivity without dSRMR
    assert status(wu_row(d2=(0, 0.01, np.nan))).metric_decision == "decision_not_evaluable"


def test_unidentified_baseline_gives_no_cfa_evidence():
    s = status(wu_row(w0=dict(df=-1, vcov_min_eigen=np.nan), w1=dict(df=5)))
    assert s.baseline_status == "identification_not_established"
    assert s.threshold_decision.startswith("not_evaluable") and s.cfa_conclusion == "none"
    assert "baseline_not_identified" in s.cfa_flags


def test_inadmissible_constrained_model_makes_the_step_not_evaluable_not_rejected():
    s = status(wu_row(w2=dict(theta_min_eigen=-0.3), w3=False))
    assert s.metric_decision == "not_evaluable_inadmissible" and s.cfa_conclusion == "configural_fit"


def test_order_sensitivity_candidates_follow_v2_statuses(tmp_path):
    clean = pd.DataFrame([wu_row()])
    assert S.cfa_order_sensitivity(clean, pd.DataFrame(), tmp_path, 3, 1, V2, tmp_path).empty
    # admissible-by-v1 (converged, no Heywood) but inadmissible under v2 (non-PD residual matrix): a candidate
    bad = pd.DataFrame([wu_row(w0=dict(theta_min_eigen=-0.2))])
    with pytest.raises(Exception):                                                            # the refit needs real data
        S.cfa_order_sensitivity(bad, pd.DataFrame(), tmp_path, 1, 1, V2, tmp_path)


def test_baseline_misfit_and_untestable_fit():
    s = status(wu_row(w0=dict(srmr=0.12)))
    assert s.cfa_conclusion == "none" and s.cfa_no_conclusion_reason == "baseline_misfit"
    assert "invariance_supported_under_baseline_misfit" in s.cfa_flags
    assert s.metric_decision == "supported" and s.metric_invariance == "not_tested" and not s.metric_test_available
    s = status(wu_row(w0=dict(df=0), w1=dict(df=5), w2=dict(df=8), w3=dict(df=11)))
    assert s.cfa_conclusion == "scalar" and "baseline_fit_not_testable" in s.cfa_flags     # df = 0: later steps testable
    s = status(wu_row(w0=dict(df=0), w1=dict(df=5), d1=(-0.05, 0.01, 0), w2=False, w3=False))
    assert s.cfa_conclusion == "none" and s.cfa_no_conclusion_reason == "baseline_fit_not_testable"


def test_near_boundary_models_are_kept_but_can_be_excluded_in_sensitivity():
    r = wu_row(w0=dict(max_std_loading=0.99))
    assert status(r).cfa_conclusion == "scalar" and "near_boundary_loading" in status(r).cfa_flags
    assert status(r, exclude_near_boundary=True).cfa_conclusion == "none"


def test_millsap_as_primary_branch():
    s = status(wu_row(d1=(-0.05, 0.01, 0), w2=False, w3=False), primary_branch="millsap")
    assert s.metric_decision == "supported" and s.cfa_conclusion == "metric"


@pytest.mark.parametrize("conclusion,reason,cfa_run,applic,paired,phi,priority,expected", [
    ("scalar", "", True, True, True, .9, "cfa", "scalar_invariance"),
    ("configural_fit", "", True, True, True, .80, "cfa", "configural_fit"),                 # CFA first; conflict flagged elsewhere
    ("configural_fit", "", True, True, True, .80, "esem", "structurally_discrepant"),       # sensitivity: ESEM first
    ("none", "baseline_misfit", True, True, True, .90, "cfa", "congruence_supported"),      # ESEM used when CFA gives no conclusion
    ("none", "baseline_misfit", True, True, False, .90, "cfa", "cfa_misfit"),
    ("none", "baseline_fit_not_testable", True, False, False, np.nan, "cfa", "cfa_fit_not_testable"),
    ("none", "inadmissible", True, True, False, .20, "cfa", "model_not_estimable"),        # inadmissible ESEM is never evidence
    ("none", "nonconverged", True, False, False, np.nan, "cfa", "model_not_estimable"),
    ("none", "not_fitted", False, False, False, np.nan, "cfa", "not_tested"),
])
def test_structural_evidence_level_v2(conclusion, reason, cfa_run, applic, paired, phi, priority, expected):
    assert S.structural_evidence_level_v2(conclusion, reason, cfa_run, applic, paired, phi, 0.85, priority) == expected


def test_instrument_summary_follows_cfa_conclusion():
    from types import SimpleNamespace as NS
    rows = {"good": wu_row(), "misfit": wu_row(w0=dict(srmr=0.12)), "noconv": wu_row(w0=dict(converged=False))}
    res = pd.DataFrame([{**r, "scale_id": k} for k, r in rows.items()])
    study = NS(registry=NS(scales={k: NS(instrument_id=1) for k in rows}, instruments={1: NS(name_en="I1"), 2: NS(name_en="I2")}))
    st = S.cfa_scale_status_v2(res, V2)
    assert st.set_index("scale_id").loc["misfit", "metric_decision"] == "supported"        # supported step, misfitting baseline
    inst = S.cfa_by_instrument(res, study, st).iloc[0]
    assert inst.configural_status == "mixed" and inst.n_scales_configural_ok == 1
    assert (inst.n_scales_cfa_misfit, inst.n_scales_not_estimable) == (1, 1)
    # only the scale with a CFA conclusion counts; the instrument level is not established
    assert inst.n_scales_metric_supported == 1 and inst.n_scales_scalar_supported == 1
    assert inst.metric_status == "not_established" and inst.scalar_status == "not_established"
    # every scale concluded -> supported; a rejected step of a concluded scale -> not_supported
    two = pd.DataFrame([{**wu_row(), "scale_id": "good"}, {**wu_row(d2=(-0.05, 0.01, 0.01), w3=False), "scale_id": "noconv"}])
    inst = S.cfa_by_instrument(two, study, S.cfa_scale_status_v2(two, V2)).iloc[0]
    assert inst.configural_status == "supported" and inst.metric_status == "not_supported" and inst.n_scales_metric_rejected == 1
    one = pd.DataFrame([{**wu_row(w0=dict(srmr=0.12)), "scale_id": "misfit"}])
    assert S.cfa_by_instrument(one, study, S.cfa_scale_status_v2(one, V2)).iloc[0].configural_status == "misfit"


def test_step_diagnostics_long_table():
    d = S.cfa_step_diagnostics(pd.DataFrame([wu_row(d1=(-0.05, 0.01, 0), w2=False, w3=False)]), V2)
    wu = d[d.sequence == "wu"].set_index("step")
    assert list(wu.index) == ["baseline", "threshold", "metric", "scalar"]
    assert wu.loc["baseline", "decision"] == "supported" and wu.loc["threshold", "decision"] == "rejected"
    assert wu.loc["threshold", "delta_df"] == 5 and wu.loc["threshold", "delta_cfi"] == -0.05
    assert wu.loc["metric", "status"] == "not_fitted_prior_step_not_supported"
