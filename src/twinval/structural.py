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


def cfa_scale_status(results: pd.DataFrame) -> pd.DataFrame:
    """Per scale: admissibility of the configural model and the metric / scalar invariance verdicts.

    * ``configural_cfa_admissible`` - the configural Wu model converged without a Heywood case;
    * ``metric``  - loading invariance from the Millsap/Liu scheme (unconditional);
    * ``scalar``  - Wu chain thresholds + loadings + intercepts, interpretable only if thresholds are supported.
    Scales whose configural model is inadmissible get ``not_tested`` for both.
    """
    out = []
    for r in results.to_dict("records"):
        ok = bool(r.get("wu_configural_converged") is True or r.get("wu_configural_converged") == True) and \
            not (r.get("wu_configural_heywood") is True or r.get("wu_configural_heywood") == True)
        metric, scalar = "not_tested", "not_tested"
        if ok:
            s = r.get("ms_loading_status")
            metric = "yes" if s == "supported" else "no" if s == "rejected" else "not_tested"
            if r.get("wu_threshold_status") == "supported":
                s = r.get("wu_conditional_scalar_status")
                scalar = "yes" if s == "supported" else "no" if s == "rejected" else "not_tested"
        out.append({"scale_id": r["scale_id"], "configural_cfa_admissible": ok, "metric_invariance": metric,
                    "scalar_invariance": scalar, "wu_threshold_status": r.get("wu_threshold_status"),
                    "metric_test_available": ok and r.get("ms_loading_status") in ("supported", "rejected"),
                    "scalar_test_available": ok and r.get("wu_threshold_status") == "supported"})
    return pd.DataFrame(out)


def cfa_by_instrument(results: pd.DataFrame, study: Study) -> pd.DataFrame:
    """Instrument-level synthesis. Metric (Millsap/Liu) and scalar (Wu chain) verdicts come from DIFFERENT
    identification schemes and are reported in separate columns - they do not form one nested ladder."""
    def collapse(vals):
        vals = vals.dropna()
        if len(vals) == 0:
            return "not_estimable"
        if (vals == "rejected").any():
            return "not_supported"
        return "supported" if (vals == "supported").all() else "not_estimable"

    reg = study.registry.scales
    res = results.copy()
    res["instrument_id"] = res["scale_id"].map(lambda s: reg[s].instrument_id)
    rows = []
    for iid, sub in res.groupby("instrument_id"):
        conf_ok = (sub["wu_configural_converged"] == True) & (sub["wu_configural_heywood"] != True)
        if conf_ok.all():
            conf = "supported"
        elif (sub["wu_configural_converged"] == False).any():
            conf = "nonconverged"
        elif (sub["wu_configural_heywood"] == True).any():
            conf = "inadmissible_solution"
        else:
            conf = "not_estimable"
        interp = sub[sub["wu_threshold_status"] == "supported"]
        rows.append({
            "instrument_id": iid, "instrument_name": study.registry.instruments[iid].name_en, "n_scales": len(sub),
            "configural_status": conf, "n_scales_configural_ok": int(conf_ok.sum()),
            "metric_status_millsap_primary": collapse(sub["ms_loading_status"]),
            "metric_status_wu_conditional": collapse(interp["wu_conditional_loading_status"]) if len(interp) else "not_interpreted_threshold_not_supported",
            "scalar_status_wu_chain": collapse(interp["wu_conditional_scalar_status"]) if len(interp) else "not_interpreted_threshold_not_supported",
            "n_scales_metric_supported": int((sub["ms_loading_status"] == "supported").sum()),
            "n_scales_metric_rejected": int((sub["ms_loading_status"] == "rejected").sum()),
            "n_scales_scalar_supported_of_interpretable": int((interp["wu_conditional_scalar_status"] == "supported").sum()),
            "n_scales_scalar_rejected_of_interpretable": int((interp["wu_conditional_scalar_status"] == "rejected").sum()),
            "n_scales_scalar_not_interpretable": int(len(sub) - len(interp))})
    return pd.DataFrame(rows)


def cfa_order_sensitivity(results: pd.DataFrame, manifest: pd.DataFrame, workdir: Path, n_perm: int, seed: int,
                          cfg, logs: Path) -> pd.DataFrame:
    """Numerical-stability check of the CFA verdicts.

    Weighted-least-squares fits of small, nearly saturated ordinal models can converge or not depending only on the
    ORDER of the rows (floating point differences). For every *candidate* scale - one whose configural model did not
    converge / is inadmissible, or whose threshold step is not testable - the CFA is refitted on ``n_perm`` random
    row orders and the scale is flagged ``order_sensitive`` if the admissibility of the configural model is not the
    same in all fits.
    """
    def configural_ok(df):
        return (df["wu_configural_converged"] == True) & (df["wu_configural_heywood"] != True)

    clean_steps = {"supported", "rejected"}
    cand = results[~(configural_ok(results) & results["wu_threshold_status"].isin(clean_steps))]["scale_id"].tolist()
    if not cand or n_perm < 1:
        return pd.DataFrame(columns=["scale_id", "n_fits", "share_configural_admissible", "order_sensitive"])
    rng = np.random.default_rng(seed)
    ok_all = {sid: [bool(configural_ok(results.set_index("scale_id").loc[[sid]]).iloc[0])] for sid in cand}
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
        pd.DataFrame(rows).to_csv(mpath, index=False)
        run_r("paired_cfa_invariance.R", [mpath, pdir / "results", cfg.thresholds.cfa_delta_cfi_min,
                                          cfg.thresholds.cfa_delta_rmsea_max], logs / "cfa_order_check.log", cfg.structural.rscript)
        r = pd.read_csv(pdir / "results" / "invariance_results.csv")
        ok = configural_ok(r).tolist()
        for sid, o in zip(r["scale_id"], ok):
            ok_all[sid].append(bool(o))
    out = [{"scale_id": sid, "n_fits": len(v), "share_configural_admissible": float(np.mean(v)),
            "order_sensitive": bool(0 < np.mean(v) < 1)} for sid, v in ok_all.items()]
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


def structural_evidence_level(configural_ok: bool, metric: str, scalar: str, esem_applicable: bool,
                              paired_admissible: bool, phi: float, phi_congruence: float = 0.85) -> str:
    """Seven-step, priority-ordered structural status of a scale.

    Formal CFA evidence (a test with a decision rule) takes precedence over ESEM/Tucker evidence (descriptive):

    1 ``scalar_invariance`` 2 ``metric_invariance`` 3 ``configural_only`` (admissible configural CFA, metric not
    confirmed) 4 ``congruence_supported`` / 5 ``structurally_discrepant`` (only if the CFA configural model is
    inadmissible but BOTH ESEM solutions are admissible; phi >= / < threshold) 6 ``model_inadmissible`` (ESEM applies but
    at least one solution is inadmissible and there is no stronger CFA evidence) 7 ``not_tested``.
    """
    if scalar == "yes":
        return "scalar_invariance"
    if metric == "yes":
        return "metric_invariance"
    if configural_ok:
        return "configural_only"
    if esem_applicable and paired_admissible:
        return "congruence_supported" if (pd.notna(phi) and phi >= phi_congruence) else "structurally_discrepant"
    if esem_applicable:
        return "model_inadmissible"
    return "not_tested"


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
            "structural_evidence_level": structural_evidence_level(conf, metric, scalar, applicable, paired, phi,
                                                                   cfg.thresholds.phi_congruence)})
    out = pd.DataFrame(rows)
    out["cfa_order_sensitive"] = False
    if order_check is not None and len(order_check):
        out["cfa_order_sensitive"] = out["scale_id"].isin(order_check.loc[order_check["order_sensitive"], "scale_id"])
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
        scalar = c["scalar_status_wu_chain"] if c is not None else "not_estimable"
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
                     "metric_status_millsap_primary": c["metric_status_millsap_primary"] if c is not None else None,
                     "scalar_status_wu_chain": scalar, "validation_trap_strict": trap,
                     "n_scales_metric_supported": c["n_scales_metric_supported"] if c is not None else np.nan,
                     "n_scales_metric_rejected": c["n_scales_metric_rejected"] if c is not None else np.nan})
    out = pd.DataFrame(rows)
    if mean_ci is not None and len(mean_ci):
        out = out.merge(mean_ci, on="instrument_id", how="left")
    return out
