"""Level 1 - profile-level fidelity: can the twin of person i be told apart from the twins of other people,
and how close is a twin profile to the profile of its own human?

Profile = vector of human-standardised scores over the 60 non-redundant scales (scale ``s`` of instrument
``m``: ``z = (X - mean_H) / SD_H``). Every instrument gets the same vote regardless of how many scales it has:

* **Equal-weight RMS distance** (final identification model)::

      D_ij = (1 / M) * sum_m sqrt( (1 / k_m) * sum_{s in m} (zH_is - zT_js)^2 )

* **Full** profiles use the z-scores as they are; **distinctive** profiles subtract, for every scale, the sample
  mean of humans and of twins separately, so that a systematic twin-wide shift cannot help or hurt identification.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps
from scipy.spatial import procrustes
from sklearn.covariance import LedoitWolf

from .io import Study
from .stats import (boot_indices, ccc_1d, cohens_dz, human_z, own_cross_stats, pct_ci)


# =============================================================================== distance / similarity matrices
def rms_distance_matrix(zH: np.ndarray, zT: np.ndarray, ids: list[str], groups: dict[int, list[str]]) -> np.ndarray:
    """Equal-instrument-weight RMS distance matrix (rows: human i, columns: twin j)."""
    col = {c: i for i, c in enumerate(ids)}
    n = zH.shape[0]
    acc = np.zeros((n, n))
    for cols in groups.values():
        ix = [col[c] for c in cols]
        Hm, Tm = zH[:, ix], zT[:, ix]
        sq = (Hm ** 2).sum(1, keepdims=True) + (Tm ** 2).sum(1, keepdims=True).T - 2 * Hm @ Tm.T
        acc += np.sqrt(np.clip(sq, 0, None) / Hm.shape[1])
    return acc / len(groups)


def _row_normalize(A: np.ndarray) -> np.ndarray:
    Ac = A - A.mean(1, keepdims=True)
    nrm = np.linalg.norm(Ac, axis=1, keepdims=True)
    return Ac / np.where(nrm == 0, 1.0, nrm)


def pearson_similarity(zH: np.ndarray, zT: np.ndarray) -> np.ndarray:
    return _row_normalize(zH) @ _row_normalize(zT).T


def eqw_fisher_similarity(zH, zT, ids, groups) -> np.ndarray:
    """Mean over instruments (with >= 2 scales) of Fisher-z transformed profile correlations."""
    col = {c: i for i, c in enumerate(ids)}
    acc, m = np.zeros((zH.shape[0],) * 2), 0
    for cols in groups.values():
        if len(cols) < 2:
            continue
        ix = [col[c] for c in cols]
        r = np.clip(_row_normalize(zH[:, ix]) @ _row_normalize(zT[:, ix]).T, -0.999999, 0.999999)
        acc += np.arctanh(r)
        m += 1
    return np.tanh(acc / max(m, 1))


def _distinctive(zH, zT):
    return zH - zH.mean(0, keepdims=True), zT - zT.mean(0, keepdims=True)


# =============================================================================== identification
def _summarize(M: np.ndarray, lower_is_better: bool, boot_idx: np.ndarray, label: str) -> tuple[dict, dict]:
    s = own_cross_stats(M, lower_is_better)
    n = M.shape[0]
    B = boot_idx.shape[0]
    top = {"top1": (s["rank"] == 1), "top5": (s["rank"] <= 5), "top10": (s["rank"] <= 10)}
    out = {"spec": label, "n": n}
    for name, flag in top.items():
        b = flag.astype(float)[boot_idx].mean(1)
        lo, hi = np.percentile(b, [2.5, 97.5])
        out[f"{name}_pct"], out[f"{name}_ci_lo"], out[f"{name}_ci_hi"] = flag.mean() * 100, lo * 100, hi * 100
    for name, arr in (("auc", s["auc"]), ("mean_rank", s["rank"].astype(float)), ("own_mean", s["own"]),
                      ("cross_mean", s["cross_mean"])):
        b = arr[boot_idx].mean(1)
        lo, hi = np.percentile(b, [2.5, 97.5])
        out[name], out[f"{name}_ci_lo"], out[f"{name}_ci_hi"] = float(arr.mean()), lo, hi
    out["median_rank"] = float(np.median(s["rank"]))
    out["rank_iqr_lo"], out["rank_iqr_hi"] = np.percentile(s["rank"], [25, 75])
    diff = (s["cross_mean"] - s["own"]) if lower_is_better else (s["own"] - s["cross_mean"])
    db = diff[boot_idx].mean(1)
    dzb = diff[boot_idx].mean(1) / diff[boot_idx].std(1, ddof=1)
    out["diff_mean"], out["diff_ci_lo"], out["diff_ci_hi"] = float(diff.mean()), *np.percentile(db, [2.5, 97.5])
    out["dz"], out["dz_ci_lo"], out["dz_ci_hi"] = cohens_dz(diff), *np.percentile(dzb, [2.5, 97.5])
    t, p = (sps.ttest_rel(s["cross_mean"], s["own"]) if lower_is_better else sps.ttest_rel(s["own"], s["cross_mean"]))
    out["paired_t"], out["paired_t_p"] = float(t), float(p)
    return out, s


def _permutation(M: np.ndarray, lower_is_better: bool, rng: np.random.Generator, n_perm: int) -> dict:
    n = M.shape[0]
    best = np.argmin(M, 1) if lower_is_better else np.argmax(M, 1)
    order = np.argsort(M if lower_is_better else -M, 1)
    top5 = order[:, :5]
    obs_top1 = float(np.mean(best == np.arange(n)))
    obs_top5 = float(np.mean([i in top5[i] for i in range(n)]))
    obs_auc = float(own_cross_stats(M, lower_is_better)["auc"].mean())
    idx = np.arange(n)
    off = ~np.eye(n, dtype=bool)
    t1, t5, ta = np.empty(n_perm), np.empty(n_perm), np.empty(n_perm)
    for p in range(n_perm):
        perm = rng.permutation(n)
        t1[p] = np.mean(best == perm)
        t5[p] = np.mean([perm[i] in top5[i] for i in range(n)])
        own_p = M[idx, perm]
        better = (own_p[:, None] < M) if lower_is_better else (own_p[:, None] > M)
        eq = own_p[:, None] == M
        mask = off.copy()
        mask[idx, perm] = False
        ta[p] = (((better & mask).sum(1) + 0.5 * (eq & mask).sum(1)) / mask.sum(1)).mean()
    pv = lambda obs, null: float((1 + np.sum(null >= obs)) / (n_perm + 1))
    return {"n_permutations": n_perm, "obs_top1": obs_top1, "obs_top5": obs_top5, "obs_auc": obs_auc,
            "null_top1_mean": float(t1.mean()), "null_auc_mean": float(ta.mean()),
            "p_top1": pv(obs_top1, t1), "p_top5": pv(obs_top5, t5), "p_auc": pv(obs_auc, ta)}


def identification(study: Study, cfg, rng: np.random.Generator) -> dict[str, pd.DataFrame]:
    """Own-twin identification with the final equal-weight RMS-distance model (full and distinctive profiles),
    permutation tests, a robustness table over alternative specifications and leave-one-instrument-out."""
    ids = study.nonredundant_ids
    groups = study.groups(ids)
    H, T = study.matrices(ids)
    zH, zT = human_z(H, T)
    dH, dT = _distinctive(zH, zT)
    n = study.n
    bidx = boot_indices(rng, cfg.profile.n_boot, n)

    D_full = rms_distance_matrix(zH, zT, ids, groups)
    D_dist = rms_distance_matrix(dH, dT, ids, groups)
    sum_full, s_full = _summarize(D_full, True, bidx, "full")
    sum_dist, s_dist = _summarize(D_dist, True, bidx, "distinctive")
    for d, lab in ((sum_full, "full"), (sum_dist, "distinctive")):
        d["n_instruments"], d["n_scales"] = len(groups), len(ids)
    summary = pd.DataFrame([sum_full, sum_dist])

    perm_rows = []
    for lab, M in (("full", D_full), ("distinctive", D_dist)):
        perm_rows.append({"spec": lab, **_permutation(M, True, rng, cfg.profile.n_perm)})
    permutation = pd.DataFrame(perm_rows)

    per_person = pd.DataFrame({
        "id": study.ids,
        "own_distance_full": s_full["own"], "cross_distance_full": s_full["cross_mean"], "rank_full": s_full["rank"],
        "auc_full": s_full["auc"], "percentile_full": s_full["percentile"],
        "own_distance_distinctive": s_dist["own"], "cross_distance_distinctive": s_dist["cross_mean"],
        "rank_distinctive": s_dist["rank"], "auc_distinctive": s_dist["auc"], "percentile_distinctive": s_dist["percentile"]})

    out = {"identification_summary": summary, "identification_permutation": permutation,
           "identification_per_person": per_person}

    if cfg.profile.robustness_specs:
        rows = []
        all_ids = study.scale_ids
        Ha, Ta = study.matrices(all_ids)
        zHa, zTa = human_z(Ha, Ta)
        dHa, dTa = _distinctive(zHa, zTa)
        g_all, g_nr = study.groups(all_ids), groups

        def add(name, n_coord, n_instr, metric, centred, M, lower):
            s = own_cross_stats(M, lower)
            diff = (s["cross_mean"] - s["own"]) if lower else (s["own"] - s["cross_mean"])
            rows.append({"spec": name, "n_scales": n_coord, "n_instruments": n_instr, "metric": metric,
                         "centred": centred, "top1_pct": float(np.mean(s["rank"] == 1) * 100),
                         "top5_pct": float(np.mean(s["rank"] <= 5) * 100), "auc": float(s["auc"].mean()),
                         "dz": cohens_dz(diff)})
        add("Pearson, all scales", len(all_ids), len(g_all), "Pearson", "no", pearson_similarity(zHa, zTa), False)
        add("Pearson, all scales", len(all_ids), len(g_all), "Pearson", "yes", pearson_similarity(dHa, dTa), False)
        add("Pearson, non-redundant scales", len(ids), len(g_nr), "Pearson", "no", pearson_similarity(zH, zT), False)
        add("Pearson, non-redundant scales", len(ids), len(g_nr), "Pearson", "yes", pearson_similarity(dH, dT), False)
        n_multi = sum(1 for c in g_nr.values() if len(c) >= 2)
        if n_multi:
            add("Equal-weight Fisher-z Pearson, non-redundant", len(ids), n_multi, "Fisher-z mean", "no",
                eqw_fisher_similarity(zH, zT, ids, g_nr), False)
            add("Equal-weight Fisher-z Pearson, non-redundant", len(ids), n_multi, "Fisher-z mean", "yes",
                eqw_fisher_similarity(dH, dT, ids, g_nr), False)
        add("Equal-weight RMS distance (final)", len(ids), len(g_nr), "RMS distance", "no", D_full, True)
        add("Equal-weight RMS distance (final)", len(ids), len(g_nr), "RMS distance", "yes", D_dist, True)
        out["identification_robustness"] = pd.DataFrame(rows)

    if cfg.profile.leave_one_instrument_out and len(groups) >= 3:
        rows = []
        base = {"full": summary.loc[0], "distinctive": summary.loc[1]}
        for tid in groups:
            g = {t: c for t, c in groups.items() if t != tid}
            rec = {"excluded_instrument_id": tid, "excluded_instrument": study.registry.instruments[tid].name_en}
            for lab, (a, b) in (("full", (zH, zT)), ("distinctive", (dH, dT))):
                s = own_cross_stats(rms_distance_matrix(a, b, ids, g), True)
                rec[f"top1_pct_{lab}"] = float(np.mean(s["rank"] == 1) * 100)
                rec[f"auc_{lab}"] = float(s["auc"].mean())
                rec[f"delta_top1_pp_{lab}"] = rec[f"top1_pct_{lab}"] - float(base[lab]["top1_pct"])
                rec[f"delta_auc_{lab}"] = rec[f"auc_{lab}"] - float(base[lab]["auc"])
            rows.append(rec)
        out["identification_leave_one_instrument_out"] = pd.DataFrame(rows)
    return out


# =============================================================================== calibration / agreement / error
def _weights(ids, groups) -> np.ndarray:
    """Per-scale weight ``1 / (k_m * M)``: every instrument contributes 1/M in total."""
    col = {c: i for i, c in enumerate(ids)}
    w = np.zeros(len(ids))
    for cols in groups.values():
        for c in cols:
            w[col[c]] = 1.0 / (len(cols) * len(groups))
    return w


def _weighted_ccc_matrix(zH, zT, w):
    muH, muT = zH @ w, zT @ w
    Hc, Tc = zH - muH[:, None], zT - muT[:, None]
    varH, varT = (Hc ** 2) @ w, (Tc ** 2) @ w
    cov = (Hc * w[None, :]) @ Tc.T
    return 2 * cov / (varH[:, None] + varT[None, :] + (muH[:, None] - muT[None, :]) ** 2)


def calibration(study: Study, cfg, rng: np.random.Generator) -> dict[str, pd.DataFrame]:
    """Profile calibration in the SAME specification as the identification model: non-redundant scales,
    equal instrument weights, human-based z. Reported for full and distinctive (= bias-corrected) profiles.

    * profile CCC of a person: Lin's CCC between the person's human and twin profiles (own) and, as a baseline,
      averaged over all other twins (cross);
    * MAE / RMSE of the profile, full and after removing the mean twin bias of each scale;
    * mean absolute bias and RMS bias across scales (equal instrument weights).
    """
    ids = study.nonredundant_ids
    groups = study.groups(ids)
    H, T = study.matrices(ids)
    zH, zT = human_z(H, T)
    dH, dT = _distinctive(zH, zT)
    w = _weights(ids, groups)
    n = study.n

    def own_cross(M):
        own = np.diag(M).copy()
        off = ~np.eye(n, dtype=bool)
        return own, np.array([M[i, off[i]].mean() for i in range(n)])

    own_f, cross_f = own_cross(_weighted_ccc_matrix(zH, zT, w))
    own_d, cross_d = own_cross(_weighted_ccc_matrix(dH, dT, w))
    mae_f, rmse_f = np.abs(zH - zT) @ w, np.sqrt(((zH - zT) ** 2) @ w)
    mae_b, rmse_b = np.abs(dH - dT) @ w, np.sqrt(((dH - dT) ** 2) @ w)
    bias = (zT - zH).mean(0)

    bias_df = pd.DataFrame({
        "scale_id": ids, "instrument_id": [study.registry.scales[s].instrument_id for s in ids], "weight": w,
        "bias_twin_minus_human_sd": bias,
        "bias_ci_lo": np.nan, "bias_ci_hi": np.nan})
    bidx = boot_indices(rng, cfg.profile.n_boot, n)
    bb = (zT - zH)[bidx].mean(1)                       # (B, S)
    bias_df["bias_ci_lo"], bias_df["bias_ci_hi"] = pct_ci(bb)
    bias_df["scale_name"] = [study.registry.scales[s].name_en for s in ids]
    bias_df = bias_df.sort_values("bias_twin_minus_human_sd", key=np.abs, ascending=False)

    per_person = pd.DataFrame({"id": study.ids, "own_ccc_full": own_f, "cross_ccc_full_mean": cross_f,
                               "own_ccc_distinctive": own_d, "cross_ccc_distinctive_mean": cross_d,
                               "mae_full": mae_f, "rmse_full": rmse_f, "mae_bias_corrected": mae_b,
                               "rmse_bias_corrected": rmse_b})

    def add(rows, name, arr):
        b = arr[bidx].mean(1)
        lo, hi = np.percentile(b, [2.5, 97.5])
        rows.append({"metric": name, "statistic": "mean", "value": float(arr.mean()), "ci_lo": lo, "ci_hi": hi,
                     "median": float(np.median(arr)), "iqr_lo": np.percentile(arr, 25), "iqr_hi": np.percentile(arr, 75)})

    rows = []
    for nm, a in (("own_ccc_full", own_f), ("cross_ccc_full_mean", cross_f), ("own_ccc_distinctive", own_d),
                  ("cross_ccc_distinctive_mean", cross_d), ("mae_full", mae_f), ("rmse_full", rmse_f),
                  ("mae_bias_corrected", mae_b), ("rmse_bias_corrected", rmse_b)):
        add(rows, nm, a)
    for lab, o, c in (("full", own_f, cross_f), ("distinctive", own_d, cross_d)):
        diff = o - c
        b = diff[bidx].mean(1)
        dzb = b / diff[bidx].std(1, ddof=1)
        rows.append({"metric": f"own_minus_cross_ccc_{lab}", "statistic": "mean_diff", "value": float(diff.mean()),
                     "ci_lo": np.percentile(b, 2.5), "ci_hi": np.percentile(b, 97.5)})
        rows.append({"metric": f"own_minus_cross_ccc_{lab}", "statistic": "cohens_dz", "value": cohens_dz(diff),
                     "ci_lo": np.percentile(dzb, 2.5), "ci_hi": np.percentile(dzb, 97.5)})
    rows.append({"metric": "mean_absolute_bias_by_scale", "statistic": "point", "value": float(np.abs(bias) @ w)})
    rows.append({"metric": "rms_bias_by_scale", "statistic": "point", "value": float(np.sqrt((bias ** 2) @ w))})
    rows.append({"metric": "share_scales_abs_bias_gt_0.5", "statistic": "point", "value": float((np.abs(bias) > 0.5).mean())})
    summary = pd.DataFrame(rows)
    summary["n_scales"], summary["n_instruments"] = len(ids), len(groups)
    summary["weighting"] = "equal instrument weight, per-scale weight 1/(k_m*M)"
    summary["standardization"] = "human-based z: (X - mean_H) / SD_H"
    return {"profile_calibration_summary": summary, "profile_calibration_per_person": per_person,
            "profile_bias_by_scale": bias_df}


# =============================================================================== Mahalanobis / mean profiles
def mahalanobis_and_mean_profile(study: Study, cfg, rng: np.random.Generator) -> dict[str, pd.DataFrame]:
    """Covariance-aware own-vs-other distance (Ledoit-Wolf shrunk covariance of human z-profiles) and a
    comparison of the *mean* human and twin profiles (Procrustes disparity, Euclidean distance, pattern r)."""
    ids = study.nonredundant_ids if cfg.profile.mahalanobis_scales == "nonredundant" else study.scale_ids
    H, T = study.matrices(ids)
    zH, zT = human_z(H, T)
    lw = LedoitWolf().fit(zH)
    L = np.linalg.cholesky(np.linalg.inv(lw.covariance_))
    A, B = zH @ L, zT @ L
    sq = (A ** 2).sum(1, keepdims=True) + (B ** 2).sum(1, keepdims=True).T - 2 * A @ B.T
    M = np.sqrt(np.clip(sq, 0, None))
    s = own_cross_stats(M, True)
    bidx = boot_indices(rng, cfg.profile.n_boot, study.n)
    b = s["own"][bidx].mean(1)
    diff = s["cross_mean"] - s["own"]
    mh, mt = H.mean(0), T.mean(0)
    _, _, disparity = procrustes(mh.reshape(-1, 1), mt.reshape(-1, 1))
    tab = pd.DataFrame([{
        "n_scales": len(ids), "ledoit_wolf_shrinkage": float(lw.shrinkage_),
        "mahalanobis_own_mean": float(s["own"].mean()), "mahalanobis_own_ci_lo": np.percentile(b, 2.5),
        "mahalanobis_own_ci_hi": np.percentile(b, 97.5), "mahalanobis_cross_mean": float(s["cross_mean"].mean()),
        "mahalanobis_dz": cohens_dz(diff), "mahalanobis_auc": float(s["auc"].mean()),
        "mean_profile_procrustes_disparity": float(disparity),
        "mean_profile_euclidean": float(np.linalg.norm(mh - mt)),
        "mean_profile_pearson_r": float(np.corrcoef(mh, mt)[0, 1]) if len(ids) > 2 else np.nan,
        "mean_profile_cosine": float(mh @ mt / (np.linalg.norm(mh) * np.linalg.norm(mt))),
        "note": "mean-profile statistics use raw score units (sum- and mean-scored scales differ in unit)"}])
    return {"profile_mahalanobis_mean_profile": tab}
