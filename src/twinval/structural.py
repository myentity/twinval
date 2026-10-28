"""Level 3 - structural (psychometric) equivalence of Human and Twin responses.

Three complementary procedures (the model fitting itself is done in R with lavaan, see ``twinval/r``):

1. **Paired ordinal CFA / measurement invariance** (per scale). Human_i and Twin_i are two measurements of the same
   person, so a *dependent-data* model is fitted. Two non-nested schemes are kept apart:
   Wu & Estabrook (2016) - thresholds tested first, loadings/intercepts conditional on them - and Millsap & Tein (2004) /
   Liu et al. (2017) - loadings tested first (unconditional loading invariance).
2. **ESEM with target rotation** (per multidimensional instrument), separately for Human and Twin, followed by an
   orthogonal Procrustes rotation of the Twin loadings onto the Human loadings and Tucker's congruence
   coefficient phi per factor (with a person-level bootstrap).
3. **Admissibility gate.** Tucker's phi is interpreted only when BOTH ESEM solutions of an instrument are admissible
   (converged, no Heywood case, no factor correlation > 1). Estimates from inadmissible solutions are kept as
   exploratory point estimates and never used to declare congruence or discrepancy.
"""
from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import orthogonal_procrustes
from scipy.optimize import linear_sum_assignment

from . import prep
from .io import Study

log = logging.getLogger("twinval")
R_DIR = Path(__file__).parent / "r"


class RError(RuntimeError):
    pass


# =============================================================================== R plumbing
def r_available(rscript: str = "Rscript") -> tuple[bool, str]:
    """Check that Rscript exists and lavaan + semTools can be loaded."""
    if shutil.which(rscript) is None:
        return False, f"'{rscript}' not found on PATH"
    p = subprocess.run([rscript, "-e", "suppressMessages({library(lavaan); library(semTools)}); cat('ok')"],
                       capture_output=True, text=True)
    if p.returncode != 0 or "ok" not in p.stdout:
        return False, "R packages lavaan and semTools are required (install.packages(c('lavaan','semTools')))"
    return True, ""


def run_r(script: str, args: list, log_path: Path, rscript: str = "Rscript") -> None:
    """Run one of the packaged R scripts, streaming stdout/stderr to ``log_path``."""
    cmd = [rscript, str(R_DIR / script)] + [str(a) for a in args]
    log.info("running %s", " ".join(cmd))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        raise RError(f"{script} failed (exit {p.returncode}); see {log_path}")


# =============================================================================== input preparation
def prepare_inputs(study: Study, cfg, workdir: Path) -> dict:
    """Harmonise categories and write everything the R scripts need under ``workdir``.

    Returns a dict with the harmonised item arrays, the harmonisation logs, the CFA manifest path and the list of
    ESEM-applicable instruments.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    hh, ht, hlog, hbad = prep.harmonize(study)
    status = prep.harmonization_status(study, hlog, hbad)
    unusable = set(hbad["scale_id"]) if len(hbad) else set()

    # ---- CFA input: one harmonised paired_wide file per usable scale
    data_dir = workdir / "cfa_input"
    data_dir.mkdir(exist_ok=True)
    rows = []
    for sid in study.scale_ids:
        sc = study.registry.scales[sid]
        if sid in unusable or (sc.redundant and cfg.structural.cfa_scales == "nonredundant"):
            continue
        k = sc.n_items
        df = pd.DataFrame({"ID": study.ids})
        for i in range(k):
            df[f"it{i + 1}_h"] = hh[sid][:, i].astype(int)
        for i in range(k):
            df[f"it{i + 1}_t"] = ht[sid][:, i].astype(int)
        path = data_dir / f"{sid}_paired_wide.csv"
        df.to_csv(path, index=False)
        rows.append({"scale_id": sid, "file": str(path.resolve()), "n_items": k, "n_pairs": study.n,
                     "instrument": study.registry.instruments[sc.instrument_id].name_en, "scale_name": sc.name_en,
                     "method_count": sc.scoring,
                     "category_status": status.set_index("scale_id").loc[sid, "category_status"]})
    manifest = pd.DataFrame(rows)
    manifest.to_csv(workdir / "cfa_manifest.csv", index=False)

    # ---- ESEM input: instruments with > 1 non-redundant scale, none of them unusable
    esem_instr, qc_frames = [], []
    (workdir / "item_battery").mkdir(exist_ok=True)
    (workdir / "target_matrices").mkdir(exist_ok=True)
    for iid in study.instrument_ids:
        ins = study.registry.instruments[iid]
        nr = [s for s in ins.scales if not s.redundant]
        if len(nr) < 2 or any(s.id in unusable for s in nr):
            continue
        battery, target = prep.item_battery(study, iid, hh, ht)
        battery.to_csv(workdir / "item_battery" / f"{iid}_item_battery.csv", index=False)
        target.to_csv(workdir / "target_matrices" / f"{iid}_target_matrix.csv")
        qc_frames.append(prep.item_ordinal_qc(battery, iid))
        esem_instr.append({"instrument_id": iid, "instrument_name": ins.name_en})
    pd.DataFrame(esem_instr, columns=["instrument_id", "instrument_name"]).to_csv(workdir / "esem_registry.csv", index=False)
    return {"items_h": hh, "items_t": ht, "harmonization_log": hlog, "harmonization_unusable": hbad,
            "harmonization_status": status, "manifest": manifest, "esem_instruments": esem_instr,
            "item_qc": pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame()}


# =============================================================================== CFA / invariance results
_STEPS = {
    "wu_configural": ("wu", "configural"), "wu_thresholds": ("wu", "thresholds"),
    "wu_conditional_loading": ("wu", "conditional_loading (thresholds+loadings)"),
    "wu_conditional_scalar": ("wu", "conditional_scalar (thresholds+loadings+intercepts)"),
    "wu_conditional_strict": ("wu", "conditional_strict (+residuals)"),
    "ms_baseline": ("millsap", "baseline"), "ms_loading": ("millsap", "loading"),
    "ms_threshold_after_loading": ("millsap", "threshold_after_loading"),
    "ms_unique_after_threshold": ("millsap", "unique_after_threshold"),
}
_DECISION_COLS = [
    "wu_threshold_status", "wu_threshold_delta_cfi", "wu_threshold_delta_rmsea", "wu_threshold_delta_srmr", "wu_threshold_lrt_p",
    "wu_conditional_loading_status", "wu_conditional_loading_delta_cfi", "wu_conditional_loading_delta_rmsea",
    "wu_conditional_scalar_status", "wu_conditional_scalar_delta_cfi", "wu_conditional_scalar_delta_rmsea",
    "wu_conditional_strict_status", "ms_loading_status", "ms_loading_delta_cfi", "ms_loading_delta_rmsea",
    "ms_loading_delta_srmr", "ms_loading_lrt_p", "ms_threshold_after_loading_status",
    "ms_threshold_after_loading_delta_cfi", "ms_unique_after_threshold_status"]


def parse_cfa(results: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Turn ``invariance_results.csv`` into long fit table, decision table and reliability table."""
    rows = []
    for r in results.to_dict("records"):
        for prefix, (scheme, label) in _STEPS.items():
            conv, cfi = r.get(f"{prefix}_converged"), r.get(f"{prefix}_cfi")
            if pd.isna(conv) and pd.isna(cfi):
                continue
            rows.append({"scale_id": r["scale_id"], "scheme": scheme, "step": label, "converged": conv, "cfi": cfi,
                         "rmsea": r.get(f"{prefix}_rmsea"), "srmr": r.get(f"{prefix}_srmr"), "df": r.get(f"{prefix}_df"),
                         "heywood": r.get(f"{prefix}_heywood")})
    fit = pd.DataFrame(rows)
    dec_cols = ["scale_id"] + [c for c in _DECISION_COLS if c in results.columns]
    decisions = results[dec_cols].copy()
    rel_cols = [c for c in ("wu_cr_h", "wu_cr_t", "wu_ave_h", "wu_ave_t", "wu_tucker_phi", "wu_latent_cor_oriented",
                            "wu_negative_loadings_h", "wu_negative_loadings_t") if c in results.columns]
    rel = results[["scale_id"] + rel_cols].rename(columns={"wu_cr_h": "cr_h", "wu_cr_t": "cr_t", "wu_ave_h": "ave_h",
                                                            "wu_ave_t": "ave_t"})
    return {"paired_cfa_invariance_fit": fit, "invariance_decisions": decisions, "cfa_reliability": rel}


# =============================================================================== CFA status
# Three properties are kept apart for EVERY model that a conclusion uses (DECISIONS.md):
#   estimability  - nonconverged / identification_not_established / inadmissible / admissible (+ warning flags)
#   absolute fit  - only for the baseline model with df > 0 (SRMR primary, CFI/TLI/RMSEA CI as diagnostics)
#   invariance    - support of the added equality constraints in one pre-specified sequence (Wu primary)
#
# branch -> baseline model, then (step, model, decision prefix in the R output, constraint kind for delta-SRMR)
_BRANCHES = {
    "wu": ("wu_configural", [("threshold", "wu_thresholds", "thresholds"),
                             ("metric", "wu_conditional_loading", "loadings"),
                             ("scalar", "wu_conditional_scalar", "thresholds")]),
    "millsap": ("ms_baseline", [("metric", "ms_loading", "loadings"),
                                ("scalar", "ms_threshold_after_loading", "thresholds")]),
}


def _num(r, key) -> float:
    v = r.get(key, np.nan)
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def model_status(r, prefix: str, t) -> tuple[str, list[str]]:
    """Estimability of one fitted lavaan model from the per-model diagnostics of ``paired_cfa_invariance.R``.

    Returns ``(status, flags)`` with status ``not_fitted`` | ``nonconverged`` | ``identification_not_established`` |
    ``inadmissible`` | ``admissible``. Near-boundary estimates are warning flags, never exclusions.
    """
    conv = r.get(f"{prefix}_converged", np.nan)
    if pd.isna(conv):
        return "not_fitted", []
    if not bool(conv):
        return "nonconverged", []
    df, vmin = _num(r, f"{prefix}_df"), _num(r, f"{prefix}_vcov_min_eigen")
    if np.isnan(df) or df < 0 or np.isnan(vmin) or vmin < -t.cfa_pd_tolerance:
        return "identification_not_established", []
    th, ps = _num(r, f"{prefix}_theta_min_eigen"), _num(r, f"{prefix}_psi_min_eigen")
    rc, fc = _num(r, f"{prefix}_max_abs_resid_cor"), _num(r, f"{prefix}_abs_factor_cor")
    bad = bool(r.get(f"{prefix}_heywood") == True) or th < -t.cfa_pd_tolerance or ps < -t.cfa_pd_tolerance \
        or rc > 1 + t.cfa_bound_tolerance or fc > 1 + t.cfa_bound_tolerance
    if bad:
        return "inadmissible", []
    flags = []
    if _num(r, f"{prefix}_max_std_loading") >= t.near_boundary_loading:
        flags.append("near_boundary_loading")
    if rc >= t.near_boundary_correlation:
        flags.append("near_boundary_residual_correlation")
    if fc >= t.near_boundary_correlation:
        flags.append("near_boundary_latent_correlation")
    return "admissible", flags


def baseline_fit(r, prefix: str, t) -> tuple[str, list[str]]:
    """Absolute fit of an admissible baseline model: ``supported`` | ``not_supported`` | ``not_testable`` | ``index_missing``."""
    df, srmr = _num(r, f"{prefix}_df"), _num(r, f"{prefix}_srmr")
    if df == 0:
        return "not_testable", ["baseline_fit_not_testable"]
    if np.isnan(srmr):
        return "index_missing", []
    if srmr > t.cfa_srmr_fit_max:
        return "not_supported", []
    cfi, rlo = _num(r, f"{prefix}_cfi"), _num(r, f"{prefix}_rmsea_ci_lo")
    conflict = (not np.isnan(cfi) and cfi < t.cfa_fit_conflict_cfi_min) or (not np.isnan(rlo) and rlo > t.cfa_fit_conflict_rmsea_ci_lo_max)
    return "supported", (["conflicting_fit_diagnostics"] if conflict else [])


def step_decision(r, less: str, more: str, dec: str, kind: str, prior_supported: bool, usable: dict, cfg) -> str:
    """Invariance decision for the comparison ``less`` -> ``more`` (models given by their R-output prefixes).

    ``supported`` / ``rejected`` only when both models are usable, the step adds testable constraints (df_diff > 0)
    and every change-in-fit index is available; otherwise the reason is returned. A step whose preceding step is not
    supported is not interpreted.
    """
    t = cfg.thresholds
    if not prior_supported:
        return "not_interpreted_prior_not_supported"
    for p in (less, more):
        if usable[p] != "admissible":
            return f"not_evaluable_{usable[p]}"
    df_diff = _num(r, f"{more}_df") - _num(r, f"{less}_df")
    if np.isnan(df_diff) or df_diff <= 0:
        return "not_evaluable_no_added_constraints"
    d_cfi, d_rmsea, d_srmr = (_num(r, f"{dec}_delta_{k}") for k in ("cfi", "rmsea", "srmr"))
    use_srmr = cfg.structural.use_delta_srmr
    if np.isnan(d_cfi) or np.isnan(d_rmsea) or (use_srmr and np.isnan(d_srmr)):
        return "decision_not_evaluable"
    srmr_max = t.cfa_delta_srmr_loadings_max if kind == "loadings" else t.cfa_delta_srmr_thresholds_max
    ok = d_cfi >= t.cfa_delta_cfi_min and d_rmsea <= t.cfa_delta_rmsea_max and (not use_srmr or d_srmr <= srmr_max)
    return "supported" if ok else "rejected"


def _branch_status(r, branch: str, cfg) -> dict:
    """Baseline estimability and fit, step decisions and the CFA conclusion of one scale in one branch."""
    t, s = cfg.thresholds, cfg.structural
    base, steps = _BRANCHES[branch]
    models = [base] + [m for _, m, _ in steps]
    status, flags = {}, {}
    for m in models:
        status[m], flags[m] = model_status(r, m, t)
        if s.exclude_near_boundary and status[m] == "admissible" and flags[m]:
            status[m] = "excluded_near_boundary"
    out = {"baseline_status": status[base], "baseline_fit": "not_evaluated"}
    used_flags: list[str] = []
    if status[base] == "admissible":
        out["baseline_fit"], fit_flags = baseline_fit(r, base, t)
        used_flags += flags[base] + fit_flags
    elif status[base] == "identification_not_established" and _num(r, f"{base}_df") < 0:
        used_flags.append("baseline_not_identified")
    prior, less, level = True, base, "none"
    for name, more, kind in steps:
        dec = {"threshold": "wu_threshold", "metric": "wu_conditional_loading" if branch == "wu" else "ms_loading",
               "scalar": "wu_conditional_scalar" if branch == "wu" else "ms_threshold_after_loading"}[name]
        out[f"{name}_decision"] = step_decision(r, less, more, dec, kind, prior, status, cfg)
        prior = out[f"{name}_decision"] == "supported"
        if prior:
            used_flags += flags[more]
            if name in ("metric", "scalar"):
                level = name
        less = more
    if out.get("threshold_decision") == "supported":
        used_flags.append("threshold_supported")
    # the CFA gives a conclusion only from an admissible baseline whose fit is supported or (df = 0) not testable
    fit = out["baseline_fit"]
    if status[base] != "admissible":
        conclusion, reason = "none", status[base]
    elif fit == "supported" or (fit == "not_testable" and level != "none"):
        conclusion, reason = (level if level != "none" else "configural_fit"), ""
    elif fit == "not_testable":
        conclusion, reason = "none", "baseline_fit_not_testable"
    else:
        conclusion, reason = "none", "baseline_misfit"
        if level != "none":
            used_flags.append("invariance_supported_under_baseline_misfit")
    out.update(cfa_conclusion=conclusion, cfa_no_conclusion_reason=reason, flags=used_flags)
    return out


def cfa_scale_status(results: pd.DataFrame, cfg) -> pd.DataFrame:
    """Per scale: estimability and fit of the baseline, every step decision of the primary branch,
    the CFA conclusion, and the secondary branch as a sensitivity result (flags ``identification_scheme_sensitive`` and
    ``identification_route_dependent``)."""
    primary = cfg.structural.primary_branch
    secondary = "millsap" if primary == "wu" else "wu"
    rows = []
    for r in results.to_dict("records"):
        p, q = _branch_status(r, primary, cfg), _branch_status(r, secondary, cfg)
        flags = list(dict.fromkeys(p["flags"]))
        p_metric, q_metric = p.get("metric_decision"), q.get("metric_decision")
        # the secondary branch never replaces a primary step; it only adds one of two flags:
        #  identification_scheme_sensitive  both loading tests were evaluated and disagree (direct conflict)
        #  identification_route_dependent   the primary sequence gives a CFA conclusion but does not interpret the loading
        #                                   step (its threshold step was rejected), while the secondary loading step is
        #                                   supported - the available conclusion depends on the sequence, the tests do not
        #                                   contradict each other
        # a missing CFA conclusion because of baseline misfit or estimability never triggers either flag
        ev = ("supported", "rejected")
        if p_metric in ev and q_metric in ev and p_metric != q_metric:
            flags.append("identification_scheme_sensitive")
        elif p["cfa_conclusion"] != "none" and p_metric == "not_interpreted_prior_not_supported" and q_metric == "supported":
            flags.append("identification_route_dependent")
        # summary columns follow the CFA conclusion: a step supported under a misfitting baseline is only a flag
        concl = p["cfa_conclusion"]
        metric_yn = "yes" if concl in ("metric", "scalar") else "no" if concl != "none" and p_metric == "rejected" else "not_tested"
        scalar_yn = "yes" if concl == "scalar" else "no" if concl == "metric" and p.get("scalar_decision") == "rejected" else "not_tested"
        rows.append({
            "scale_id": r["scale_id"], "primary_branch": primary,
            "baseline_status": p["baseline_status"], "baseline_fit": p["baseline_fit"],
            "threshold_decision": p.get("threshold_decision", "not_applicable"),
            "metric_decision": p_metric, "scalar_decision": p.get("scalar_decision"),
            "cfa_conclusion": p["cfa_conclusion"], "cfa_no_conclusion_reason": p["cfa_no_conclusion_reason"],
            "configural_cfa_admissible": p["baseline_status"] == "admissible",
            "metric_invariance": metric_yn, "scalar_invariance": scalar_yn,
            "metric_test_available": metric_yn != "not_tested", "scalar_test_available": scalar_yn != "not_tested",
            f"sensitivity_{secondary}_baseline_status": q["baseline_status"],
            f"sensitivity_{secondary}_metric_decision": q_metric,
            f"sensitivity_{secondary}_scalar_decision": q.get("scalar_decision"),
            f"sensitivity_{secondary}_cfa_conclusion": q["cfa_conclusion"],
            "cfa_flags": ";".join(flags)})
    return pd.DataFrame(rows)


def cfa_step_diagnostics(results: pd.DataFrame, cfg) -> pd.DataFrame:
    """Machine-readable basis of every CFA decision: one row per scale and fitted model of both sequences.

    Columns: the model and the model it is compared with, its estimability status and near-boundary flags, the numerical
    diagnostics used for the status, absolute fit (df, scaled chi-square, CFI, TLI, RMSEA with 90% CI, SRMR), and for
    comparisons the change in df and in CFI / RMSEA / SRMR with the step decision of ``cfa_scale_status``.
    """
    t = cfg.thresholds
    diag_cols = ("converged", "heywood", "post_check", "vcov_min_eigen", "theta_min_eigen", "psi_min_eigen", "max_std_loading",
                 "min_std_resid_var", "max_abs_resid_cor", "abs_factor_cor")
    fit_cols = ("df", "chisq", "cfi", "tli", "rmsea", "rmsea_ci_lo", "rmsea_ci_hi", "srmr")
    dec_key = {"wu": {"threshold": "wu_threshold", "metric": "wu_conditional_loading", "scalar": "wu_conditional_scalar"},
               "millsap": {"metric": "ms_loading", "scalar": "ms_threshold_after_loading"}}
    rows = []
    for r in results.to_dict("records"):
        for branch, (base, steps) in _BRANCHES.items():
            dec = _branch_status(r, branch, cfg)
            seq = [("baseline", base, None)] + [(name, m, prev) for (name, m, _), prev in zip(steps, [base] + [m for _, m, _ in steps])]
            for name, m, prev in seq:
                st, fl = model_status(r, m, t)
                row = {"scale_id": r["scale_id"], "sequence": branch, "step": name, "model": m, "compared_with": prev or "",
                       "status": st, "near_boundary_flags": ";".join(fl)}
                row.update({c: r.get(f"{m}_{c}", np.nan) for c in fit_cols + diag_cols})
                if prev:
                    row["delta_df"] = _num(r, f"{m}_df") - _num(r, f"{prev}_df")
                    k = dec_key[branch][name]
                    row.update({f"delta_{x}": r.get(f"{k}_delta_{x}", np.nan) for x in ("cfi", "rmsea", "srmr")})
                    row["decision"] = dec.get(f"{name}_decision")
                    if row["decision"] == "not_interpreted_prior_not_supported" and np.isnan(_num(r, f"{m}_df")):
                        row["status"] = "not_fitted_prior_step_not_supported"    # the R script skips the model
                else:
                    row["decision"] = dec.get("baseline_fit")
                rows.append(row)
    return pd.DataFrame(rows)


def structural_evidence_level(conclusion: str, reason: str, cfa_run: bool, esem_applicable: bool,
                              paired_admissible: bool, phi: float, phi_congruence: float = 0.85,
                              priority: str = "cfa") -> str:
    """Structural status S.

    1-3 ``scalar_invariance`` / ``metric_invariance`` / ``configural_fit`` - CFA conclusion (admissible baseline, fit
    supported or not testable with supported invariance, primary sequence);
    4-5 ``congruence_supported`` / ``structurally_discrepant`` - admissible paired ESEM when the CFA gives no conclusion;
    analysis statuses (no structural evidence): ``cfa_misfit``, ``cfa_fit_not_testable``, ``model_not_estimable``,
    ``not_tested``. With ``priority='esem'`` (sensitivity) an admissible paired ESEM is used first.
    """
    esem = None
    if esem_applicable and paired_admissible and pd.notna(phi):
        esem = "congruence_supported" if phi >= phi_congruence else "structurally_discrepant"
    cfa = {"scalar": "scalar_invariance", "metric": "metric_invariance", "configural_fit": "configural_fit"}.get(conclusion)
    if priority == "esem" and esem:
        return esem
    if cfa:
        return cfa
    if esem:
        return esem
    if reason == "baseline_misfit":
        return "cfa_misfit"
    if reason == "baseline_fit_not_testable":
        return "cfa_fit_not_testable"
    if cfa_run or esem_applicable:
        return "model_not_estimable"
    return "not_tested"


def cfa_by_instrument(study: Study, status: pd.DataFrame) -> pd.DataFrame:
    """Instrument-level synthesis. Every verdict follows the CFA conclusion of the scales (primary branch); step
    decisions without a conclusion (e.g. supported under a misfitting baseline) are not counted:
      configural_status  ``supported`` (every scale has a CFA conclusion) | ``misfit`` / ``not_estimable`` /
                         ``fit_not_testable`` (every scale without a conclusion has this reason) | ``mixed``
      metric_status      ``supported`` (every scale metric or scalar) | ``not_supported`` (a scale with a conclusion
                         rejects the loading step) | ``not_established`` (otherwise)
      scalar_status      the same for the scalar level."""
    reg = study.registry.scales
    s = status.copy()
    s["instrument_id"] = s["scale_id"].map(lambda x: reg[x].instrument_id)
    reason_name = {"baseline_misfit": "misfit", "baseline_fit_not_testable": "fit_not_testable"}

    def level(yn: pd.Series) -> str:
        if (yn == "yes").all():
            return "supported"
        return "not_supported" if (yn == "no").any() else "not_established"

    rows = []
    for iid, sub in s.groupby("instrument_id"):
        ok = sub["cfa_conclusion"] != "none"
        reasons = {reason_name.get(r, "not_estimable") for r in sub.loc[~ok, "cfa_no_conclusion_reason"]}
        conf = "supported" if ok.all() else reasons.pop() if len(reasons) == 1 else "mixed"
        rows.append({
            "instrument_id": iid, "instrument_name": study.registry.instruments[iid].name_en, "n_scales": len(sub),
            "primary_branch": sub["primary_branch"].iloc[0], "configural_status": conf, "n_scales_configural_ok": int(ok.sum()),
            "n_scales_cfa_misfit": int((sub["cfa_no_conclusion_reason"] == "baseline_misfit").sum()),
            "n_scales_fit_not_testable": int((sub["cfa_no_conclusion_reason"] == "baseline_fit_not_testable").sum()),
            "n_scales_not_estimable": int((~ok & ~sub["cfa_no_conclusion_reason"].isin(list(reason_name))).sum()),
            "metric_status": level(sub["metric_invariance"]), "scalar_status": level(sub["scalar_invariance"]),
            "n_scales_metric_supported": int((sub["metric_invariance"] == "yes").sum()),
            "n_scales_metric_rejected": int((sub["metric_invariance"] == "no").sum()),
            "n_scales_scalar_supported": int((sub["scalar_invariance"] == "yes").sum()),
            "n_scales_scalar_rejected": int((sub["scalar_invariance"] == "no").sum())})
    return pd.DataFrame(rows)


def cfa_order_sensitivity(results: pd.DataFrame, manifest: pd.DataFrame, workdir: Path, n_perm: int, seed: int,
                          cfg, logs: Path) -> pd.DataFrame:
    """Numerical-stability check of the CFA verdicts.

    Weighted-least-squares fits of small, nearly saturated ordinal models can converge or not depending only on the
    ORDER of the rows (floating point differences). For every *candidate* scale - one whose configural model did not
    converge / is inadmissible, or whose threshold step is not testable - the CFA is refitted on ``n_perm`` random
    row orders and the scale is flagged ``order_sensitive`` if the admissibility of the baseline or the CFA conclusion
    is not the same in all fits (candidates are scales whose baseline is not admissible or whose first step is not
    evaluable).
    """
    clean_steps = {"supported", "rejected"}

    def verdicts(df):
        """Per scale: (baseline admissible, candidate for the check, verdict compared across row orders)."""
        st = cfa_scale_status(df, cfg)
        ok = st["baseline_status"] == "admissible"
        first = st["threshold_decision"] if cfg.structural.primary_branch == "wu" else st["metric_decision"]
        # the verdict is admissibility + CFA conclusion; the reason a baseline fails (nonconverged / not identified /
        # inadmissible) does not change the structural status and is only reported (column baseline_statuses)
        return ok.tolist(), (~(ok & first.isin(clean_steps))).tolist(), \
            list(zip(ok.astype(str) + "|" + st["cfa_conclusion"], st["baseline_status"]))

    ok0, cand0, v0 = verdicts(results)
    cand = [sid for sid, c in zip(results["scale_id"], cand0) if c]
    if not cand or n_perm < 1:
        return pd.DataFrame(columns=["scale_id", "n_fits", "share_configural_admissible", "n_distinct_verdicts", "baseline_statuses", "order_sensitive"])
    rng = np.random.default_rng(seed)
    ok_all = {sid: [o] for sid, o in zip(results["scale_id"], ok0) if sid in cand}
    vd_all = {sid: [v] for sid, v in zip(results["scale_id"], v0) if sid in cand}
    for j in range(n_perm):
        pdir = workdir / "order_check" / f"perm{j + 1}"
        pdir.mkdir(parents=True, exist_ok=True)
        rows = []
        for sid in cand:
            src = manifest[manifest.scale_id == sid].iloc[0].to_dict()
            d = pd.read_csv(src["file"]).iloc[rng.permutation(len(pd.read_csv(src["file"])))]
            path = pdir / f"{sid}_paired_wide.csv"
            d.to_csv(path, index=False)
            src["file"] = str(path.resolve())
            rows.append(src)
        mpath = pdir / "manifest.csv"
        res_path = pdir / "results" / "invariance_results.csv"
        prev = pd.read_csv(mpath) if mpath.exists() else None
        # resume: the permutation is drawn anyway (same random stream); a refit of the same candidates is reused
        reuse = bool(cfg.structural.get("resume", False)) and res_path.exists() and prev is not None \
            and sorted(prev["scale_id"]) == sorted(cand) and "wu_configural_vcov_min_eigen" in pd.read_csv(res_path, nrows=0)
        if not reuse:
            pd.DataFrame(rows).to_csv(mpath, index=False)
            run_r("paired_cfa_invariance.R", [mpath, pdir / "results", cfg.thresholds.cfa_delta_cfi_min,
                                              cfg.thresholds.cfa_delta_rmsea_max], logs / "cfa_order_check.log", cfg.structural.rscript)
        r = pd.read_csv(res_path)
        ok, _, vd = verdicts(r)
        for sid, o, v in zip(r["scale_id"], ok, vd):
            ok_all[sid].append(bool(o))
            vd_all[sid].append(v)
    # sensitive when the baseline status or the CFA conclusion differs between row orders
    out = [{"scale_id": sid, "n_fits": len(v), "share_configural_admissible": float(np.mean(v)),
            "n_distinct_verdicts": len({x[0] for x in vd_all[sid]}),
            "baseline_statuses": ";".join(sorted({x[1] for x in vd_all[sid]})),
            "order_sensitive": bool(len({x[0] for x in vd_all[sid]}) > 1)}
           for sid, v in ok_all.items()]
    return pd.DataFrame(out)


# =============================================================================== ESEM / Tucker
def tucker_phi(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.sqrt((a ** 2).sum()) * np.sqrt((b ** 2).sum()))
    return float((a * b).sum() / den) if den > 0 else float("nan")


def phi_category(phi: float) -> str:
    return "very_high" if phi >= 0.95 else "high" if phi >= 0.90 else "borderline" if phi >= 0.85 else "insufficient"


def tucker_point_estimates(workdir: Path, esem_diag: pd.DataFrame, instruments: list[dict]) -> dict[str, pd.DataFrame]:
    """Orthogonal Procrustes (Twin -> Human), Hungarian factor matching check and Tucker's phi per factor.

    The Human solution is the reference. The factor matching is *verified* against the theoretical order (both
    sides use the same target) and not assumed.
    """
    phi_rows, sum_rows, aligned_rows = [], [], []
    status = esem_diag.pivot(index="instrument_id", columns="side", values="status") if len(esem_diag) else pd.DataFrame()
    for ins in instruments:
        iid, name = ins["instrument_id"], ins["instrument_name"]
        lh, lt = workdir / "esem_loadings_matrices" / f"{iid}_human_loadings.csv", workdir / "esem_loadings_matrices" / f"{iid}_twin_loadings.csv"
        sh = status.loc[iid, "human"] if iid in status.index else "not_estimated"
        st = status.loc[iid, "twin"] if iid in status.index else "not_estimated"
        if not (lh.exists() and lt.exists()):
            sum_rows.append({"instrument_id": iid, "instrument_name": name, "status_human": sh, "status_twin": st,
                             "phi_mean": np.nan, "phi_min": np.nan, "n_factors": np.nan})
            continue
        Lh, Lt = pd.read_csv(lh, index_col=0), pd.read_csv(lt, index_col=0)
        if Lh.isna().any().any() or Lt.isna().any().any():
            sum_rows.append({"instrument_id": iid, "instrument_name": name, "status_human": sh, "status_twin": st,
                             "phi_mean": np.nan, "phi_min": np.nan, "n_factors": Lh.shape[1]})
            continue
        assert list(Lh.index) == list(Lt.index) and list(Lh.columns) == list(Lt.columns)
        H, T = Lh.values, Lt.values
        R, _ = orthogonal_procrustes(T, H)
        Ta = T @ R
        nf = H.shape[1]
        cong = np.array([[tucker_phi(H[:, i], Ta[:, j]) for j in range(nf)] for i in range(nf)])
        row_ind, col_ind = linear_sum_assignment(-cong)
        matched = bool(np.all(row_ind == col_ind))
        target = pd.read_csv(workdir / "target_matrices" / f"{iid}_target_matrix.csv", index_col=0).loc[Lh.index]
        phis = []
        for i, fac in enumerate(Lh.columns):
            j = int(col_ind[list(row_ind).index(i)])
            phi = tucker_phi(H[:, i], Ta[:, j])
            tgt = (target[fac] == 1).values
            phis.append(phi)
            phi_rows.append({"instrument_id": iid, "instrument_name": name, "factor": fac, "matched_twin_factor": Lh.columns[j],
                             "phi_full_column": phi, "phi_target_items_only": tucker_phi(H[tgt, i], Ta[tgt, j]),
                             "status_human": sh, "status_twin": st, "category": phi_category(phi)})
            for it in Lh.index:
                aligned_rows.append({"instrument_id": iid, "instrument_name": name, "factor": fac, "item": it,
                                     "loading_human": Lh.loc[it, fac], "loading_twin_aligned": float(Ta[list(Lh.index).index(it), j]),
                                     "is_target": int(target.loc[it, fac])})
        phis = np.array(phis)
        sum_rows.append({"instrument_id": iid, "instrument_name": name, "n_factors": nf, "status_human": sh, "status_twin": st,
                         "matched_theoretical_order": matched, "procrustes_is_reflection": bool(np.linalg.det(R) < 0),
                         "phi_mean": float(phis.mean()), "phi_median": float(np.median(phis)), "phi_min": float(phis.min()),
                         "n_phi_ge_095": int((phis >= 0.95).sum()), "n_phi_090_095": int(((phis >= 0.90) & (phis < 0.95)).sum()),
                         "n_phi_085_090": int(((phis >= 0.85) & (phis < 0.90)).sum()), "n_phi_lt_085": int((phis < 0.85).sum())})
    return {"tucker_by_factor": pd.DataFrame(phi_rows), "tucker_summary": pd.DataFrame(sum_rows),
            "esem_aligned_loadings": pd.DataFrame(aligned_rows)}


def merge_tucker_bootstrap(tucker_by_factor: pd.DataFrame, boot: pd.DataFrame, valid_min: float = 0.80):
    """Percentile CIs of phi per factor, share of valid replicates and instrument-level CI of the *mean* phi.

    The instrument-level CI is built from the per-replicate mean of phi over the factors (not from the mean of the
    factor-level CI limits, which is not a confidence interval of the mean).
    """
    rows = []
    for (iid, fac), sub in boot.groupby(["instrument_id", "factor"]):
        n_valid = int(sub["valid"].sum())
        frac = n_valid / len(sub)
        v = sub.loc[sub["valid"], "phi"].dropna()
        lo, hi = (np.percentile(v, [2.5, 97.5]) if len(v) >= 10 else (np.nan, np.nan))
        rows.append({"instrument_id": iid, "factor": fac, "n_valid_bootstrap": n_valid, "n_total_bootstrap": len(sub),
                     "frac_valid_bootstrap": frac, "phi_boot_mean": float(v.mean()) if len(v) >= 10 else np.nan,
                     "phi_ci_lo": lo, "phi_ci_hi": hi, "ci_stable": bool(frac >= valid_min)})
    ci = pd.DataFrame(rows)
    by_factor = tucker_by_factor.merge(ci, on=["instrument_id", "factor"], how="left") if len(tucker_by_factor) else ci
    mean_rows = []
    for iid, sub in boot.groupby("instrument_id"):
        per_rep = sub[sub["valid"]].groupby("rep")["phi"].mean()
        lo, hi = (np.percentile(per_rep, [2.5, 97.5]) if len(per_rep) >= 10 else (np.nan, np.nan))
        mean_rows.append({"instrument_id": iid, "n_valid_reps_for_mean_phi": len(per_rep), "phi_mean_ci_lo": lo, "phi_mean_ci_hi": hi})
    return by_factor, pd.DataFrame(mean_rows)


# =============================================================================== admissibility gate & evidence level
def admissibility(esem_diag: pd.DataFrame, instrument_ids: list[int]) -> pd.DataFrame:
    """Per instrument: is the ESEM of each side admissible, and are BOTH admissible (paired)?"""
    rows = []
    piv = esem_diag.pivot(index="instrument_id", columns="side", values="status") if len(esem_diag) else pd.DataFrame()
    for iid in instrument_ids:
        h = bool(iid in piv.index and piv.loc[iid].get("human") == "ok")
        t = bool(iid in piv.index and piv.loc[iid].get("twin") == "ok")
        rows.append({"instrument_id": iid, "human_esem_admissible": h, "twin_esem_admissible": t,
                     "paired_esem_admissible": h and t})
    return pd.DataFrame(rows)


def scale_structural_table(study: Study, cfg, cfa_status: pd.DataFrame | None, tucker_by_factor: pd.DataFrame | None,
                           gate: pd.DataFrame | None, esem_instruments: list[int],
                           order_check: pd.DataFrame | None = None) -> pd.DataFrame:
    """One row per scale with every structural flag used by the cross-level synthesis."""
    reg = study.registry.scales
    cfa_map = cfa_status.set_index("scale_id") if cfa_status is not None and len(cfa_status) else pd.DataFrame()
    gate_map = gate.set_index("instrument_id") if gate is not None and len(gate) else pd.DataFrame()
    tk = tucker_by_factor
    rows = []
    for sid in study.scale_ids:
        sc = reg[sid]
        c = cfa_map.loc[sid] if sid in cfa_map.index else None
        conf = bool(c["configural_cfa_admissible"]) if c is not None else False
        metric = c["metric_invariance"] if c is not None else "not_tested"
        scalar = c["scalar_invariance"] if c is not None else "not_tested"
        applicable = (not sc.redundant) and sc.instrument_id in esem_instruments
        phi = frac = np.nan
        ci_lo = ci_hi = np.nan
        stable = False
        if applicable and tk is not None and len(tk):
            m = tk[(tk.instrument_id == sc.instrument_id) & (tk.factor == sid)]
            if len(m):
                phi = float(m["phi_full_column"].iloc[0])
                frac = float(m["frac_valid_bootstrap"].iloc[0]) if "frac_valid_bootstrap" in m else np.nan
                ci_lo = float(m["phi_ci_lo"].iloc[0]) if "phi_ci_lo" in m else np.nan
                ci_hi = float(m["phi_ci_hi"].iloc[0]) if "phi_ci_hi" in m else np.nan
                stable = bool(pd.notna(frac) and frac >= cfg.thresholds.bootstrap_valid_min)
        g = gate_map.loc[sc.instrument_id] if sc.instrument_id in gate_map.index else None
        h_ok = bool(g["human_esem_admissible"]) if g is not None else False
        t_ok = bool(g["twin_esem_admissible"]) if g is not None else False
        paired = h_ok and t_ok
        interpretable = bool(applicable and paired and pd.notna(phi))
        conclusion = c["cfa_conclusion"] if c is not None else "none"
        reason = c["cfa_no_conclusion_reason"] if c is not None else "not_fitted"
        level = structural_evidence_level(conclusion, reason, c is not None, applicable, paired, phi,
                                             cfg.thresholds.phi_congruence, cfg.structural.structural_priority)
        flags = [f for f in (str(c["cfa_flags"]).split(";") if c is not None and pd.notna(c["cfa_flags"]) else []) if f]
        if interpretable and level in ("scalar_invariance", "metric_invariance", "configural_fit") \
                and phi < cfg.thresholds.phi_congruence:
            flags.append("cfa_esem_conflict")
        if level in ("congruence_supported", "structurally_discrepant") and not stable:
            flags.append("bootstrap_unstable")
        rows.append({
            "scale_id": sid, "instrument_id": sc.instrument_id, "configural_cfa_admissible": conf,
            "metric_invariance": metric, "scalar_invariance": scalar,
            "metric_test_available": bool(c["metric_test_available"]) if c is not None else False,
            "scalar_test_available": bool(c["scalar_test_available"]) if c is not None else False,
            "esem_applicable": applicable, "human_esem_admissible": h_ok, "twin_esem_admissible": t_ok,
            "paired_esem_admissible": paired, "tucker_point_available": bool(pd.notna(phi)),
            "tucker_interpretable": interpretable, "tucker_bootstrap_stable": stable,
            "tucker_phi_exploratory": phi, "tucker_phi_primary": phi if interpretable else np.nan,
            "tucker_phi_ci_lo": ci_lo, "tucker_phi_ci_hi": ci_hi, "tucker_frac_valid_bootstrap": frac,
            "structural_evidence_level": level, "structural_flags": ";".join(dict.fromkeys(flags))})
        if c is not None:
            rows[-1].update({k: c[k] for k in ("baseline_status", "baseline_fit", "threshold_decision", "metric_decision",
                                              "scalar_decision", "cfa_conclusion", "cfa_no_conclusion_reason")})
    out = pd.DataFrame(rows)
    out["cfa_order_sensitive"] = False
    if order_check is not None and len(order_check):
        out["cfa_order_sensitive"] = out["scale_id"].isin(order_check.loc[order_check["order_sensitive"], "scale_id"])
        sens = out["cfa_order_sensitive"]
        out.loc[sens, "structural_flags"] = out.loc[sens, "structural_flags"].map(lambda s: ";".join(filter(None, [s, "order_sensitive"])))
    return out


def instrument_master(study: Study, cfa_instr: pd.DataFrame | None, esem_fit: pd.DataFrame | None,
                      tucker_summary: pd.DataFrame | None, mean_ci: pd.DataFrame | None, gate: pd.DataFrame | None,
                      by_factor: pd.DataFrame | None) -> pd.DataFrame:
    """Instrument-level structural synthesis including a conservative 'validation trap' flag.

    *Validation trap* = the structure looks similar (high phi from admissible, bootstrap-stable ESEM) while strict
    measurement equivalence (scalar) is not established. It is assigned only when both ESEM solutions are
    admissible, the configural CFA is supported and the bootstrap CI of phi is stable.
    """
    rows = []
    for iid in study.instrument_ids:
        ins = study.registry.instruments[iid]
        applicable = study.registry.esem_applicable(iid)
        c = cfa_instr[cfa_instr.instrument_id == iid].iloc[0] if cfa_instr is not None and (cfa_instr.instrument_id == iid).any() else None
        t = tucker_summary[tucker_summary.instrument_id == iid].iloc[0] if tucker_summary is not None and len(tucker_summary) and (tucker_summary.instrument_id == iid).any() else None
        g = gate[gate.instrument_id == iid].iloc[0] if gate is not None and (gate.instrument_id == iid).any() else None
        bf = by_factor[by_factor.instrument_id == iid] if by_factor is not None and len(by_factor) else pd.DataFrame()
        phi_mean = float(t["phi_mean"]) if t is not None and pd.notna(t.get("phi_mean")) else np.nan
        phi_min = float(t["phi_min"]) if t is not None and pd.notna(t.get("phi_min")) else np.nan
        paired = bool(g["paired_esem_admissible"]) if g is not None else False
        stable = bool(len(bf) and "ci_stable" in bf and bf["ci_stable"].all())
        conf = c["configural_status"] if c is not None else "not_estimable"
        scalar = c["scalar_status"] if c is not None else "not_estimable"
        trap = "not_applicable"
        if applicable:
            if conf != "supported":
                trap = "not_applicable_configural_not_supported"
            elif not paired:
                trap = "not_assessable_esem_inadmissible"
            elif pd.isna(phi_mean):
                trap = "not_assessable"
            else:
                strict = phi_mean >= 0.95 and phi_min >= 0.90
                trap = bool(strict and scalar != "supported" and stable) if strict else False
        rows.append({"instrument_id": iid, "instrument_name": ins.name_en, "esem_applicable": applicable,
                     "paired_esem_admissible": paired, "phi_mean": phi_mean, "phi_min": phi_min,
                     "phi_ci_stable": stable, "configural_status": conf,
                     "metric_status": c["metric_status"] if c is not None else None, "scalar_status": scalar,
                     "validation_trap_strict": trap,
                     "n_scales_metric_supported": c["n_scales_metric_supported"] if c is not None else np.nan,
                     "n_scales_metric_rejected": c["n_scales_metric_rejected"] if c is not None else np.nan})
    out = pd.DataFrame(rows)
    if mean_ci is not None and len(mean_ci):
        out = out.merge(mean_ci, on="instrument_id", how="left")
    return out
