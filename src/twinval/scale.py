"""Level 2 - scale-level fidelity: individuation, calibration, recalibration, reliability correction.

Per scale ``s`` (score = mean or sum of items, computed from the paired item data):

* **Individuation** - Pearson r, Spearman rho and a scale-level identification AUC
  (does the twin's score of person i lie closer to person i's human score than the twins of other people?).
* **Calibration / agreement** - Lin's CCC decomposed as ``CCC = r * C_b`` (C_b = accuracy component).
* **Cross-validated recalibration** - location-scale (LS) and OLS recalibration of twin scores, evaluated
  out-of-fold, with a nested person-level bootstrap whose folds keep duplicated respondents together.
* **Disattenuation** - r corrected for the unreliability of both sides (Cronbach's alpha for k >= 3,
  Spearman-Brown for k = 2 items, composite reliability from the CFA when available).
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.model_selection import KFold

from .io import Study
from .stats import (bh_adjust, boot_indices, ccc_components, cohens_dz, grouped_fold_ids, human_z,
                    own_cross_stats, pct_ci, permutation_auc_p)

log = logging.getLogger("twinval")


# =============================================================================== core table
def _boot_r_ccc_cb(h: np.ndarray, t: np.ndarray, idx: np.ndarray):
    """Vectorised bootstrap of r, CCC and C_b for one scale; ``idx`` has shape (B, N)."""
    hb, tb = h[idx], t[idx]
    mh, mt = hb.mean(1), tb.mean(1)
    sh, st = hb.std(1, ddof=1), tb.std(1, ddof=1)
    cov = ((hb - mh[:, None]) * (tb - mt[:, None])).sum(1) / (h.size - 1)
    valid = (sh > 0) & (st > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = cov / (sh * st)
        v = st / sh
        u = (mt - mh) / np.sqrt(sh * st)
        cb = 2 / (v + 1 / v + u ** 2)
    return r[valid], (r * cb)[valid], cb[valid], int(valid.sum())


def typology_r_cb(r: float, cb: float, r_thr: float = 0.50, cb_thr: float = 0.80) -> str:
    """Descriptive two-component classification of a scale by r and C_b (not the final fidelity typology)."""
    if r >= r_thr and cb >= cb_thr:
        return "Good reconstruction"
    if r >= r_thr:
        return "Ordering preserved / inadequate calibration"
    if cb >= cb_thr:
        return "Limited individuation"
    return "Combined limitations"


def core_table(study: Study, cfg, rng: np.random.Generator) -> pd.DataFrame:
    """Pearson r, scale-level AUC, CCC/C_b/u/v with person-level bootstrap CIs (shared resamples)."""
    n_boot, n_perm = cfg.scale.n_boot, cfg.scale.n_perm
    H, T = study.matrices()
    zH, zT = human_z(H, T)
    N = study.n
    idx = boot_indices(rng, n_boot, N)
    rows = []
    for k, sid in enumerate(study.scale_ids):
        h, t, zh, zt = H[:, k], T[:, k], zH[:, k], zT[:, k]
        r_pt, p_pt = sps.pearsonr(h, t)
        rho, p_rho = sps.spearmanr(h, t)
        rb, cb_, cbb, n_valid = _boot_r_ccc_cb(h, t, idx)
        comp = ccc_components(h, t)
        D = np.abs(zh[:, None] - zt[None, :])
        st = own_cross_stats(D, lower_is_better=True)
        auc = float(st["auc"].mean())
        auc_boot = st["auc"][idx].mean(1)
        diff = st["cross_mean"] - st["own"]
        p_perm, null_m, null_sd = permutation_auc_p(D, auc, rng, n_perm)
        r_lo, r_hi = np.percentile(rb, [2.5, 97.5])
        c_lo, c_hi = np.percentile(cb_, [2.5, 97.5])
        b_lo, b_hi = np.percentile(cbb, [2.5, 97.5])
        a_lo, a_hi = np.percentile(auc_boot, [2.5, 97.5])
        sc = study.registry.scales[sid]
        rows.append({
            "scale_id": sid, "instrument_id": sc.instrument_id, "n": N, "n_items": sc.n_items,
            "is_nonredundant": not sc.redundant,
            "pearson_r": r_pt, "pearson_r_ci_lo": r_lo, "pearson_r_ci_hi": r_hi, "pearson_p": p_pt,
            "spearman_rho": rho, "spearman_p": p_rho,
            "auc": auc, "auc_ci_lo": a_lo, "auc_ci_hi": a_hi, "auc_perm_p": p_perm,
            "auc_perm_null_mean": null_m, "auc_perm_null_sd": null_sd,
            "d_own_mean": float(st["own"].mean()), "d_cross_mean": float(st["cross_mean"].mean()),
            "auc_cohens_dz": cohens_dz(diff), "tie_share": st["tie_share"],
            "ccc": float(comp["ccc"][0]), "ccc_ci_lo": c_lo, "ccc_ci_hi": c_hi,
            "cb": float(comp["cb"][0]), "cb_ci_lo": b_lo, "cb_ci_hi": b_hi,
            "u": float(comp["u"][0]), "v": float(comp["v"][0]),
            "n_valid_bootstrap": n_valid,
        })
    res = pd.DataFrame(rows)
    res["pearson_q_bh"] = bh_adjust(res["pearson_p"])
    res["pearson_q_bh_sig"] = res["pearson_q_bh"] < cfg.thresholds.alpha_sig
    res["auc_perm_q_bh"] = bh_adjust(res["auc_perm_p"])
    res["auc_perm_q_bh_sig"] = res["auc_perm_q_bh"] < cfg.thresholds.alpha_sig
    res["typology_r_cb"] = [typology_r_cb(r, c, cfg.thresholds.r_practical, cfg.thresholds.cb_adequate)
                            for r, c in zip(res["pearson_r"], res["cb"])]
    return res


# =============================================================================== recalibration
def _fold_calibrate(Htr, Ttr, Tte):
    """LS and OLS recalibration fitted on the training fold, applied to the test fold (vectorised over scales)."""
    mh, mt = Htr.mean(0), Ttr.mean(0)
    sh, st = Htr.std(0, ddof=1), Ttr.std(0, ddof=1)
    b_ls = sh / st
    cal_ls = (mh - b_ls * mt)[None, :] + b_ls[None, :] * Tte
    cov = ((Htr - mh) * (Ttr - mt)).sum(0) / (Htr.shape[0] - 1)
    b_ols = cov / st ** 2
    cal_ols = (mh - b_ols * mt)[None, :] + b_ols[None, :] * Tte
    return cal_ls, cal_ols


def _agreement(H, X):
    c = ccc_components(H, X)
    c["mae"] = np.abs(H - X).mean(0)
    c["rmse"] = np.sqrt(((H - X) ** 2).mean(0))
    return c


def recalibration(study: Study, cfg, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Repeated K-fold CV recalibration plus a nested, leakage-free person-level bootstrap.

    * ``LS``  b = SD_H,train / SD_T,train  (corrects location and scale only; keeps the ranking of twin scores)
    * ``OLS`` b = cov / var(T_train)       (regression to the mean: also shrinks the twin scores by r)

    Returns ``(per_repeat, per_scale)``. MAE/RMSE are additionally expressed in human SD units so that
    scales scored as sums and as means can be aggregated.
    """
    H, T = study.matrices()
    N, S = H.shape
    folds, reps, B = cfg.scale.cv_folds, cfg.scale.cv_repeats, cfg.scale.n_nested_boot
    raw = _agreement(H, T)
    sd_h = H.std(0, ddof=1)

    per_rep = []
    for rep in range(reps):
        kf = KFold(n_splits=folds, shuffle=True, random_state=1000 + rep)
        oof_ls, oof_ols = np.empty((N, S)), np.empty((N, S))
        for tr, te in kf.split(H):
            oof_ls[te], oof_ols[te] = _fold_calibrate(H[tr], T[tr], T[te])
        m_ls, m_ols = _agreement(H, oof_ls), _agreement(H, oof_ols)
        for k, sid in enumerate(study.scale_ids):
            rec = {"scale_id": sid, "repeat": rep}
            for name, m in (("raw", raw), ("ls", m_ls), ("ols", m_ols)):
                for met in ("ccc", "cb", "r", "mae", "rmse"):
                    rec[f"{met}_{name}"] = float(m[met][k])
                rec[f"mae_{name}_std"] = float(m["mae"][k] / sd_h[k])
                rec[f"rmse_{name}_std"] = float(m["rmse"][k] / sd_h[k])
            rec["delta_ccc_ls"] = rec["ccc_ls"] - rec["ccc_raw"]
            rec["delta_ccc_ols"] = rec["ccc_ols"] - rec["ccc_raw"]
            per_rep.append(rec)
    per_repeat = pd.DataFrame(per_rep)

    # nested person-level bootstrap; all copies of a respondent stay in one fold
    d_ls, d_ols = np.empty((B, S)), np.empty((B, S))
    for b in range(B):
        idx = rng.integers(0, N, size=N)
        Hb, Tb = H[idx], T[idx]
        fid = grouped_fold_ids(idx, folds, rng)
        o_ls, o_ols = np.empty((N, S)), np.empty((N, S))
        for f in range(folds):
            te = fid == f
            tr = ~te
            if tr.sum() < 10 or te.sum() == 0:
                o_ls[te], o_ols[te] = Tb[te], Tb[te]
                continue
            o_ls[te], o_ols[te] = _fold_calibrate(Hb[tr], Tb[tr], Tb[te])
        c_raw = ccc_components(Hb, Tb)["ccc"]
        d_ls[b] = ccc_components(Hb, o_ls)["ccc"] - c_raw
        d_ols[b] = ccc_components(Hb, o_ols)["ccc"] - c_raw

    rows = []
    g = per_repeat.groupby("scale_id", sort=False)
    means = g.mean(numeric_only=True)
    lo_ls, hi_ls = pct_ci(d_ls)
    lo_ols, hi_ols = pct_ci(d_ols)
    for k, sid in enumerate(study.scale_ids):
        sub = per_repeat[per_repeat.scale_id == sid]
        rec = {"scale_id": sid}
        for proc, lo, hi, boot in (("ls", lo_ls, hi_ls, d_ls), ("ols", lo_ols, hi_ols, d_ols)):
            d = sub[f"delta_ccc_{proc}"].values
            rec[f"delta_ccc_{proc}"] = float(d.mean())
            rec[f"delta_ccc_{proc}_sd_repeats"] = float(d.std(ddof=1)) if len(d) > 1 else np.nan
            rec[f"delta_ccc_{proc}_n_positive_repeats"] = int((d > 0).sum())
            rec[f"delta_ccc_{proc}_boot_mean"] = float(boot[:, k].mean())
            rec[f"delta_ccc_{proc}_boot_ci_lo"] = float(lo[k])
            rec[f"delta_ccc_{proc}_boot_ci_hi"] = float(hi[k])
            rec[f"delta_ccc_{proc}_ci_excludes_zero"] = bool(lo[k] > 0 or hi[k] < 0)
            rec[f"ccc_{proc}"] = float(means.loc[sid, f"ccc_{proc}"])
            rec[f"cb_{proc}"] = float(means.loc[sid, f"cb_{proc}"])
            for met in ("mae", "rmse"):
                rec[f"{met}_{proc}_std"] = float(means.loc[sid, f"{met}_{proc}_std"])
        rec["cb_raw"] = float(means.loc[sid, "cb_raw"])
        rec["mae_raw_std"], rec["rmse_raw_std"] = float(means.loc[sid, "mae_raw_std"]), float(means.loc[sid, "rmse_raw_std"])
        rows.append(rec)
    per_scale = pd.DataFrame(rows)
    return per_repeat, per_scale


# =============================================================================== reliability / disattenuation
def cronbach_alpha_boot(items: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Cronbach's alpha of ``items[idx]`` for every bootstrap row (vectorised over the B axis)."""
    k = items.shape[1]
    ib = items[idx]                                   # (B, N, K)
    item_var = ib.var(axis=1, ddof=1).sum(1)
    tot_var = ib.sum(2).var(axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (k / (k - 1)) * (1 - item_var / tot_var)


def cronbach_alpha(items: np.ndarray) -> float:
    k = items.shape[1]
    tv = items.sum(1).var(ddof=1)
    return float((k / (k - 1)) * (1 - items.var(0, ddof=1).sum() / tv)) if tv > 0 else float("nan")


def disattenuation(study: Study, cfg, rng: np.random.Generator, core: pd.DataFrame,
                   cfa: pd.DataFrame | None = None) -> pd.DataFrame:
    """Correlation corrected for attenuation ``r_dis = r / sqrt(rel_H * rel_T)``.

    Primary estimate uses Cronbach's alpha (k >= 3 items). For two-item scales the Spearman-Brown
    coefficient of the inter-item correlation is used. If CFA results are available (``cfa`` with columns
    ``scale_id, cr_h, cr_t``), a composite-reliability variant is added as a sensitivity analysis.
    ``r_dis`` values > 1 are reported as they are (never clipped).
    """
    B = cfg.scale.n_boot_disattenuation
    r_map = core.set_index("scale_id")["pearson_r"]
    cr_map = cfa.set_index("scale_id") if cfa is not None and len(cfa) else None
    rows = []
    N = study.n
    for sid in study.scale_ids:
        h, t = study.items_h[sid], study.items_t[sid]
        k = h.shape[1]
        r = float(r_map[sid])
        rec = {"scale_id": sid, "n_items": k, "pearson_r": r, "reliability_type": None,
               "rel_human": np.nan, "rel_twin": np.nan, "reliability_ceiling": np.nan,
               "r_dis": np.nan, "r_dis_ci_lo": np.nan, "r_dis_ci_hi": np.nan, "undefined_reason": None}
        idx = boot_indices(rng, B, N)
        hs, ts = h.mean(1), t.mean(1)
        if k >= 3:
            rec["reliability_type"] = "alpha"
            ah, at = cronbach_alpha(h), cronbach_alpha(t)
            rec["rel_human"], rec["rel_twin"] = ah, at
            if ah > 0 and at > 0:
                rec["reliability_ceiling"] = float(np.sqrt(ah * at))
                rec["r_dis"] = r / rec["reliability_ceiling"]
                a_h, a_t = cronbach_alpha_boot(h, idx), cronbach_alpha_boot(t, idx)
                rb = _boot_corr(hs, ts, idx)
                with np.errstate(invalid="ignore", divide="ignore"):
                    vals = rb / np.sqrt(a_h * a_t)
                vals = vals[np.isfinite(vals) & (a_h > 0) & (a_t > 0)]
                if len(vals) > 100:
                    rec["r_dis_ci_lo"], rec["r_dis_ci_hi"] = np.percentile(vals, [2.5, 97.5])
            else:
                rec["undefined_reason"] = "alpha <= 0 on at least one side"
        elif k == 2:
            rec["reliability_type"] = "spearman_brown"
            rh = float(np.corrcoef(h[:, 0], h[:, 1])[0, 1])
            rt = float(np.corrcoef(t[:, 0], t[:, 1])[0, 1])
            sbh, sbt = 2 * rh / (1 + rh), 2 * rt / (1 + rt)
            rec["rel_human"], rec["rel_twin"] = sbh, sbt
            if sbh > 0 and sbt > 0:
                rec["reliability_ceiling"] = float(np.sqrt(sbh * sbt))
                rec["r_dis"] = r / rec["reliability_ceiling"]
                rb = _boot_corr(hs, ts, idx)
                h0, h1, t0, t1 = h[:, 0], h[:, 1], t[:, 0], t[:, 1]
                rh_b, rt_b = _boot_corr(h0, h1, idx), _boot_corr(t0, t1, idx)
                with np.errstate(invalid="ignore", divide="ignore"):
                    vals = rb / np.sqrt((2 * rh_b / (1 + rh_b)) * (2 * rt_b / (1 + rt_b)))
                vals = vals[np.isfinite(vals)]
                if len(vals) > 100:
                    rec["r_dis_ci_lo"], rec["r_dis_ci_hi"] = np.percentile(vals, [2.5, 97.5])
            else:
                rec["undefined_reason"] = "Spearman-Brown <= 0 on at least one side"
        else:
            rec["undefined_reason"] = "single-item scale"
        rec["r_dis_alpha"] = rec["r_dis"] if rec["reliability_type"] == "alpha" else np.nan      # primary (k >= 3)
        rec["r_dis_sb"] = rec["r_dis"] if rec["reliability_type"] == "spearman_brown" else np.nan   # two-item scales
        rec["flag_abs_r_dis_gt_1"] = bool(abs(rec["r_dis"]) > 1) if pd.notna(rec["r_dis"]) else False
        if cr_map is not None and sid in cr_map.index:
            crh, crt = cr_map.loc[sid, "cr_h"], cr_map.loc[sid, "cr_t"]
            if pd.notna(crh) and pd.notna(crt) and crh > 0 and crt > 0:
                rec["cr_human"], rec["cr_twin"] = float(crh), float(crt)
                rec["r_dis_cr"] = r / float(np.sqrt(crh * crt))
        rows.append(rec)
    return pd.DataFrame(rows)


def _boot_corr(x: np.ndarray, y: np.ndarray, idx: np.ndarray) -> np.ndarray:
    xb, yb = x[idx], y[idx]
    xm, ym = xb.mean(1, keepdims=True), yb.mean(1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return ((xb - xm) * (yb - ym)).sum(1) / np.sqrt(((xb - xm) ** 2).sum(1) * ((yb - ym) ** 2).sum(1))


# =============================================================================== assembly
def build_master(core: pd.DataFrame, recal: pd.DataFrame, disatt: pd.DataFrame, study: Study) -> pd.DataFrame:
    """Scale-level master table: one row per scale with every scale-level result."""
    keep_dis = ["scale_id", "reliability_type", "rel_human", "rel_twin", "reliability_ceiling", "r_dis", "r_dis_alpha", "r_dis_sb",
                "r_dis_ci_lo", "r_dis_ci_hi", "flag_abs_r_dis_gt_1"] + \
               [c for c in ("cr_human", "cr_twin", "r_dis_cr") if c in disatt.columns]
    m = core.merge(recal, on="scale_id", how="left").merge(disatt[keep_dis], on="scale_id", how="left")
    reg = study.registry.scales
    m.insert(1, "scale_name", [reg[s].name_en for s in m["scale_id"]])
    m.insert(2, "instrument_name", [study.registry.instruments[reg[s].instrument_id].name_en for s in m["scale_id"]])
    return m


def summary_tables(master: pd.DataFrame) -> pd.DataFrame:
    """Descriptive summary (mean, median, IQR, range) of the main scale-level metrics, three weighting schemes."""
    from .stats import fisher_z_mean

    rows = []
    tasks = sorted(master["instrument_id"].unique())

    def ew(col):
        return float(np.mean([master.loc[master.instrument_id == t, col].mean() for t in tasks]))

    for col, label in (("pearson_r", "Pearson r"), ("auc", "Scale-level AUC"), ("ccc", "CCC"), ("cb", "C_b"),
                       ("delta_ccc_ls", "Delta CCC (LS recalibration)"), ("delta_ccc_ols", "Delta CCC (OLS recalibration)")):
        if col not in master:
            continue
        v = master[col].dropna()
        rows.append({"metric": label, "n_scales": len(v), "mean": v.mean(), "median": v.median(),
                     "iqr_lo": v.quantile(.25), "iqr_hi": v.quantile(.75), "min": v.min(), "max": v.max(),
                     "fisher_z_mean": fisher_z_mean(v.values) if col == "pearson_r" else np.nan,
                     "instrument_equal_weight_mean": ew(col),
                     "nonredundant_mean": master.loc[master.is_nonredundant, col].mean()})
    return pd.DataFrame(rows)
