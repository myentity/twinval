"""Level 4 - item-level fidelity and mechanism.

Unit of analysis: one item of a *non-redundant* scale (nested facets and total scores re-use items and would
double-count them). For each item the twin's answer is compared with the answer of the same human:

* agreement and error (exact / adjacent agreement, MAE, RMSE in units of the response range, Spearman rho, CCC),
* bias and variance compression (SD ratio, normalised entropy, floor/ceiling use, calibration slope),
* an ANOVA-style crossed decomposition of variance (person / item / person x item + residual) for the human
  data, the twin data and their difference, with paired bootstrap contrasts,
* the link between item fidelity and the discrepancy of ESEM loadings (only for instruments whose Human and Twin
  ESEM solutions are both admissible),
* reverse-keyed versus direct items.

All range-normalised quantities use the *theoretical* response range of the instrument (registry), never the
range that happened to be observed in the sample.

DATA VERSION. ``build_panel`` accepts the raw item arrays (``study.items_h/items_t``) or the harmonised ones from
``prep.harmonize`` (rare response categories merged identically for Human and Twin when their observed category sets
differ). The pipeline uses the raw arrays by default (``item.data_version: raw``, the study's primary analysis);
``harmonized`` is a sensitivity analysis. Harmonisation re-numbers merged categories 1..k' while the range-normalised and
category-based quantities keep the instrument's response format, so for merged items MAE, bias, SD ratio, entropy and
extreme-response shares refer to a coarser coding (the highest category becomes unreachable).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.metrics import cohen_kappa_score

from .io import Study
from .stats import bh_adjust, boot_indices, cluster_bootstrap_spearman, pct_ci


# =============================================================================== panel
@dataclass
class ItemPanel:
    items: list[str]
    H: np.ndarray            # (N, K)
    T: np.ndarray
    meta: pd.DataFrame       # item_id, instrument_id, instrument_name, scale_id, scale_name, reverse_keyed, item_min, item_max, n_categories
    ids: np.ndarray

    @property
    def K(self):
        return len(self.items)


def build_panel(study: Study, items_h: dict[str, np.ndarray] | None = None,
                items_t: dict[str, np.ndarray] | None = None) -> ItemPanel:
    """Assemble the (N x K) item matrices of the non-redundant scales (harmonised arrays may be supplied)."""
    items_h = items_h or study.items_h
    items_t = items_t or study.items_t
    cols_h, cols_t, meta = [], [], []
    for sid in study.nonredundant_ids:
        sc = study.registry.scales[sid]
        ins = study.registry.instruments[sc.instrument_id]
        for k, item in enumerate(sc.items):
            cols_h.append(items_h[sid][:, k])
            cols_t.append(items_t[sid][:, k])
            meta.append({"item_id": item, "instrument_id": ins.id, "instrument_name": ins.name_en, "scale_id": sid,
                         "scale_name": sc.name_en, "reverse_keyed": item in sc.reverse_items,
                         "item_min": float(ins.response_min), "item_max": float(ins.response_max),
                         "n_categories": ins.n_categories})
    meta = pd.DataFrame(meta)
    if meta["item_id"].duplicated().any():
        raise ValueError("an item belongs to two non-redundant scales")
    return ItemPanel(meta["item_id"].tolist(), np.column_stack(cols_h), np.column_stack(cols_t), meta, study.ids)


# =============================================================================== agreement / bias / compression
def _corr_cols(A, B):
    Am, Bm = A - A.mean(0), B - B.mean(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (Am * Bm).sum(0) / np.sqrt((Am ** 2).sum(0) * (Bm ** 2).sum(0))


def _ccc_cols(H, T):
    mh, mt = H.mean(0), T.mean(0)
    sh, st = H.std(0, ddof=1), T.std(0, ddof=1)
    r = _corr_cols(H, T)
    with np.errstate(invalid="ignore", divide="ignore"):
        v, u = st / sh, (mt - mh) / np.sqrt(sh * st)
        c = r * 2 / (v + 1 / v + u ** 2)
    c[(sh == 0) | (st == 0)] = np.nan
    return c


def _rank_cols(A):
    return sps.rankdata(A, axis=0)


def _entropy_cols(A, n_cat, lo):
    """Shannon entropy of the response distribution normalised by log(number of theoretical categories)."""
    N, K = A.shape
    out = np.zeros(K)
    for c in range(int(n_cat.max())):
        p = (A == (lo + c)).mean(0)
        with np.errstate(divide="ignore", invalid="ignore"):
            out -= np.where(p > 0, p * np.log(p), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(n_cat > 1, out / np.log(n_cat), np.nan)


def agreement_and_bias(panel: ItemPanel, cfg, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Item-level agreement/error table, bias/compression table and a data audit table."""
    H, T, meta = panel.H, panel.T, panel.meta
    N, K = H.shape
    lo, hi = meta["item_min"].to_numpy(float), meta["item_max"].to_numpy(float)
    rng_ = hi - lo
    ncat = meta["n_categories"].to_numpy(float)
    err = T - H
    mh, mt, sh, st = H.mean(0), T.mean(0), H.std(0, ddof=1), T.std(0, ddof=1)
    r = _corr_cols(H, T)
    rho = _corr_cols(_rank_cols(H), _rank_cols(T))
    ccc = _ccc_cols(H, T)

    kappa, cal_slope, cal_int = np.full(K, np.nan), np.full(K, np.nan), np.full(K, np.nan)
    for j in range(K):
        if sh[j] > 0 and st[j] > 0:
            try:
                kappa[j] = cohen_kappa_score(H[:, j].astype(int), T[:, j].astype(int), weights="quadratic")
            except ValueError:
                pass
            cal_slope[j], cal_int[j] = np.polyfit(H[:, j], T[:, j], 1)
    ent_h, ent_t = _entropy_cols(H, ncat, lo), _entropy_cols(T, ncat, lo)

    # ------------------------------------------------------------ person-level bootstrap
    B = cfg.item.n_boot
    bidx = boot_indices(rng, B, N)
    b_exact, b_mae, b_bias, b_sdr = (np.empty((B, K)) for _ in range(4))
    b_rho, b_ccc, b_ent = (np.empty((B, K)) for _ in range(3))
    for b in range(B):
        ix = bidx[b]
        Hb, Tb = H[ix], T[ix]
        e = Tb - Hb
        b_exact[b] = (Hb == Tb).mean(0)
        b_mae[b] = np.abs(e).mean(0) / rng_
        b_bias[b] = e.mean(0) / rng_
        shb, stb = Hb.std(0, ddof=1), Tb.std(0, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            b_sdr[b] = np.where(shb > 0, stb / shb, np.nan)
        b_rho[b] = _corr_cols(_rank_cols(Hb), _rank_cols(Tb))
        b_ccc[b] = _ccc_cols(Hb, Tb)
        b_ent[b] = _entropy_cols(Tb, ncat, lo) - _entropy_cols(Hb, ncat, lo)
    ci = lambda a: pct_ci(a)

    ag = pd.DataFrame({"item_id": panel.items, "n_pairs": N,
                       "exact_agreement": (H == T).mean(0), "adjacent_agreement": (np.abs(err) <= 1).mean(0),
                       "mae_raw": np.abs(err).mean(0), "mae_range": np.abs(err).mean(0) / rng_,
                       "rmse_raw": np.sqrt((err ** 2).mean(0)), "rmse_range": np.sqrt((err ** 2).mean(0)) / rng_,
                       "spearman_rho": rho, "pearson_r": r, "ccc": ccc, "weighted_kappa_quadratic": kappa})
    for name, arr in (("exact_agreement", b_exact), ("mae_range", b_mae), ("spearman_rho", b_rho), ("ccc", b_ccc)):
        ag[f"{name}_ci_lo"], ag[f"{name}_ci_hi"] = ci(arr)
    ag = ag.merge(meta, on="item_id")

    # Wilcoxon signed-rank test of bias (H0: no systematic difference), BH-adjusted
    p_bias = np.full(K, np.nan)
    for j in range(K):
        if np.any(err[:, j] != 0):
            p_bias[j] = sps.wilcoxon(err[:, j]).pvalue
        else:
            p_bias[j] = 1.0
    modal_t = np.array([np.bincount((T[:, j] - lo[j]).astype(int)).max() / N for j in range(K)])
    bias = pd.DataFrame({
        "item_id": panel.items, "mean_human": mh, "mean_twin": mt, "bias_raw": mt - mh, "bias_range": (mt - mh) / rng_,
        "bias_p_wilcoxon": p_bias, "bias_q_bh": bh_adjust(p_bias),
        "sd_human": sh, "sd_twin": st, "sd_ratio": st / sh, "log_sd_ratio": np.where(st > 0, np.log(st / sh), -np.inf),
        "var_ratio": (st / sh) ** 2, "entropy_human_norm": ent_h, "entropy_twin_norm": ent_t,
        "entropy_difference": ent_t - ent_h,
        "n_categories_used_human": [len(np.unique(H[:, j])) for j in range(K)],
        "n_categories_used_twin": [len(np.unique(T[:, j])) for j in range(K)],
        "modal_category_share_twin": modal_t,
        "floor_share_twin": [(T[:, j] == lo[j]).mean() for j in range(K)],
        "ceiling_share_twin": [(T[:, j] == hi[j]).mean() for j in range(K)],
        "calibration_slope": cal_slope, "calibration_intercept": cal_int})
    # extreme categories of the allowable range on both sides: SD ratio and entropy describe compression, not its form
    fh, ch = (H == lo[None, :]).mean(0), (H == hi[None, :]).mean(0)
    ft, ct = bias["floor_share_twin"].to_numpy(), bias["ceiling_share_twin"].to_numpy()
    bias["floor_share_human"], bias["ceiling_share_human"] = fh, ch
    bias["extreme_share_human"], bias["extreme_share_twin"] = fh + ch, ft + ct
    bias["extreme_share_difference"] = (ft + ct) - (fh + ch)
    bias["category_coverage_ratio"] = bias["n_categories_used_twin"] / bias["n_categories_used_human"]
    bias["bias_range_ci_lo"], bias["bias_range_ci_hi"] = ci(b_bias)
    bias["sd_ratio_ci_lo"], bias["sd_ratio_ci_hi"] = ci(b_sdr)
    bias["entropy_difference_ci_lo"], bias["entropy_difference_ci_hi"] = ci(b_ent)
    bias = bias.merge(meta, on="item_id")

    top2 = [pd.Series(T[:, j]).value_counts(normalize=True).nlargest(2).sum() for j in range(K)]
    audit = pd.DataFrame({
        "item_id": panel.items,
        "n_categories_human": bias["n_categories_used_human"], "n_categories_twin": bias["n_categories_used_twin"],
        "zero_variance_human": sh == 0, "zero_variance_twin": st == 0,
        "twin_top2_categories_share": top2, "twin_near_degenerate_flag": np.array(top2) > 0.95}).merge(meta, on="item_id")
    return ag, bias, audit


def range_codebook(panel: ItemPanel) -> pd.DataFrame:
    """Theoretical response range per item (from the registry) next to the range observed in the sample."""
    obs_min = np.minimum(panel.H.min(0), panel.T.min(0))
    obs_max = np.maximum(panel.H.max(0), panel.T.max(0))
    cb = panel.meta[["item_id", "instrument_id", "instrument_name", "scale_id", "reverse_keyed"]].copy()
    cb["theoretical_min"] = panel.meta["item_min"].values
    cb["theoretical_max"] = panel.meta["item_max"].values
    cb["theoretical_n_categories"] = panel.meta["n_categories"].values
    cb["observed_min"], cb["observed_max"] = obs_min, obs_max
    cb["observed_range_narrower_than_theoretical"] = (obs_min > cb["theoretical_min"]) | (obs_max < cb["theoretical_max"])
    return cb


# =============================================================================== variance decomposition
def _decompose(Y: np.ndarray):
    mu = Y.mean()
    pm, im = Y.mean(1), Y.mean(0)
    resid = Y - pm[:, None] - im[None, :] + mu
    vp, vi, vr = np.var(pm - mu, ddof=1), np.var(im - mu, ddof=1), np.var(resid, ddof=1)
    return vp, vi, vr, vp + vi + vr


def variance_decomposition(panel: ItemPanel, cfg, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ANOVA-style crossed decomposition Y_ij = mu + P_i + I_j + (PI)_ij for each instrument.

    Items are rescaled to [0, 1] with the theoretical range. Absolute variances (``var_*``) are reported next to the
    percentages, and ``person_effect_r`` (contrasts table) is the correlation of the human and twin person effects.
    With one observation per person x item cell the
    interaction cannot be separated from measurement error, so it is reported as the *person x item / residual*
    component. Decomposed are the human data, the twin data and the error (twin - human). Bootstrap contrasts
    (twin - human, in percentage points of variance) are paired within every replicate.
    """
    B = cfg.item.n_boot_variance
    meta = panel.meta
    rows, contrasts, data = [], [], {}
    for iid, sub in meta.groupby("instrument_id", sort=False):
        ix = sub.index.to_numpy()
        if len(ix) < 2:
            continue
        lo, hi = sub["item_min"].to_numpy(), sub["item_max"].to_numpy()
        Yh, Yt = (panel.H[:, ix] - lo) / (hi - lo), (panel.T[:, ix] - lo) / (hi - lo)
        data[iid] = (Yh, Yt)
        name = sub["instrument_name"].iloc[0]
        pct = {}
        for side, Y in (("human", Yh), ("twin", Yt), ("error", Yt - Yh)):
            vp, vi, vr, vt = _decompose(Y)
            pct[side] = {"pct_person": vp / vt * 100, "pct_item": vi / vt * 100, "pct_person_item_residual": vr / vt * 100}
            rows.append({"instrument_id": iid, "instrument_name": name, "side": side, "n_persons": Y.shape[0],
                         "n_items": Y.shape[1], "var_person": vp, "var_item": vi, "var_person_item_residual": vr,
                         "var_total": vt, **pct[side]})
        bidx = boot_indices(rng, B, Yh.shape[0])
        d_item, d_res = np.empty(B), np.empty(B)
        for b in range(B):
            _, vih, vrh, vth = _decompose(Yh[bidx[b]])
            _, vit, vrt, vtt = _decompose(Yt[bidx[b]])
            d_item[b], d_res[b] = vit / vtt * 100 - vih / vth * 100, vrt / vtt * 100 - vrh / vth * 100
        # correspondence of the person effects (person means over the instrument's items): similar person SHARES of the
        # variance do not imply that the twin keeps the same persons in the same positions
        person_r = float(np.corrcoef(Yh.mean(1), Yt.mean(1))[0, 1])
        contrasts.append({"instrument_id": iid, "instrument_name": name, "n_valid_reps": int(np.isfinite(d_item).sum()),
                          "person_effect_r": person_r,
                          "delta_item_point": pct["twin"]["pct_item"] - pct["human"]["pct_item"],
                          "delta_item_ci_lo": np.nanpercentile(d_item, 2.5), "delta_item_ci_hi": np.nanpercentile(d_item, 97.5),
                          "delta_resid_point": pct["twin"]["pct_person_item_residual"] - pct["human"]["pct_person_item_residual"],
                          "delta_resid_ci_lo": np.nanpercentile(d_res, 2.5), "delta_resid_ci_hi": np.nanpercentile(d_res, 97.5)})
    decomp, contr = pd.DataFrame(rows), pd.DataFrame(contrasts)

    # instrument-weighted contrast: every instrument is resampled separately, then averaged with equal weight
    if len(data) >= 1:
        oi, orr = np.empty(B), np.empty(B)
        for b in range(B):
            di, dr = [], []
            for Yh, Yt in data.values():
                ix = rng.integers(0, Yh.shape[0], size=Yh.shape[0])
                _, vih, vrh, vth = _decompose(Yh[ix])
                _, vit, vrt, vtt = _decompose(Yt[ix])
                di.append(vit / vtt * 100 - vih / vth * 100)
                dr.append(vrt / vtt * 100 - vrh / vth * 100)
            oi[b], orr[b] = np.mean(di), np.mean(dr)
        tw, hu = decomp[decomp.side == "twin"], decomp[decomp.side == "human"]
        contr = pd.concat([contr, pd.DataFrame([{
            "instrument_id": "ALL_INSTRUMENT_WEIGHTED", "instrument_name": f"Instrument-weighted mean ({len(data)} instruments)",
            "n_valid_reps": B, "person_effect_r": float(np.mean([c["person_effect_r"] for c in contrasts])),
            "delta_item_point": float((tw.pct_item.values - hu.pct_item.values).mean()),
            "delta_item_ci_lo": np.percentile(oi, 2.5), "delta_item_ci_hi": np.percentile(oi, 97.5),
            "delta_resid_point": float((tw.pct_person_item_residual.values - hu.pct_person_item_residual.values).mean()),
            "delta_resid_ci_lo": np.percentile(orr, 2.5), "delta_resid_ci_hi": np.percentile(orr, 97.5)}])], ignore_index=True)
    return decomp, contr


# =============================================================================== loadings, reverse keying, summaries
def loading_discrepancy(agreement: pd.DataFrame, bias: pd.DataFrame, loadings: pd.DataFrame,
                        admissible_instruments: list[int], cfg, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Link item fidelity to the discrepancy between human and (Procrustes-aligned) twin ESEM loadings.

    ``loadings`` has columns instrument_id, factor, item, loading_human, loading_twin_aligned, is_target.
    Only instruments whose Human AND Twin ESEM are admissible enter the association table; other instruments are
    kept in the discrepancy table but flagged ``esem_admissible_both_sides = False``.
    """
    rows = []
    for (iid, item), g in loadings.groupby(["instrument_id", "item"]):
        tgt = g[g["is_target"] == 1]
        rows.append({"instrument_id": iid, "item_id": item,
                     "abs_delta_loading_target": float((tgt.loading_human - tgt.loading_twin_aligned).abs().iloc[0]) if len(tgt) else np.nan,
                     "loading_vector_distance": float(np.sqrt(((g.loading_human - g.loading_twin_aligned) ** 2).sum())),
                     "n_factors": len(g), "esem_admissible_both_sides": iid in admissible_instruments})
    disc = pd.DataFrame(rows, columns=["instrument_id", "item_id", "abs_delta_loading_target", "loading_vector_distance",
                                       "n_factors", "esem_admissible_both_sides"])
    if disc.empty:
        return disc, pd.DataFrame()
    m = disc.merge(agreement[["item_id", "spearman_rho", "ccc", "mae_range", "exact_agreement"]], on="item_id", how="left") \
            .merge(bias[["item_id", "bias_range", "log_sd_ratio", "entropy_difference"]], on="item_id", how="left")
    m["abs_bias_range"] = m["bias_range"].abs()
    m["abs_log_sd_ratio"] = m["log_sd_ratio"].abs().replace(np.inf, np.nan)
    adm = m[m["esem_admissible_both_sides"]].dropna(subset=["loading_vector_distance"])
    rows = []
    for metric in ("spearman_rho", "ccc", "mae_range", "abs_bias_range", "abs_log_sd_ratio", "entropy_difference", "exact_agreement"):
        sub = adm.dropna(subset=[metric])
        if len(sub) < 5:
            continue
        res = cluster_bootstrap_spearman(sub, "loading_vector_distance", metric, "instrument_id", cfg.synthesis.n_boot_cluster,
                                         rng, cfg.synthesis.min_clusters_for_cluster_bootstrap)
        per = {int(t): float(sps.spearmanr(g["loading_vector_distance"], g[metric])[0])
               for t, g in sub.groupby("instrument_id") if len(g) >= 5 and g[metric].nunique() > 1}
        rows.append({"metric": metric, **res, "per_instrument_rho": per})
    assoc = pd.DataFrame(rows)
    if len(assoc):
        assoc["q_bh"] = bh_adjust(assoc["p_naive"])
    return disc, assoc


def reverse_keyed_analysis(agreement: pd.DataFrame, bias: pd.DataFrame, cfg, rng: np.random.Generator) -> pd.DataFrame:
    """Reverse-keyed versus direct items on the core item metrics (cluster bootstrap over instruments)."""
    df = agreement.merge(bias[["item_id", "bias_range", "sd_ratio", "log_sd_ratio", "entropy_difference"]], on="item_id")
    df["abs_bias_range"] = df["bias_range"].abs()
    df["abs_log_sd_ratio"] = df["log_sd_ratio"].abs().replace(np.inf, np.nan)
    rows = []
    for metric in ("exact_agreement", "mae_range", "spearman_rho", "ccc", "bias_range", "abs_bias_range", "sd_ratio",
                   "abs_log_sd_ratio", "entropy_difference"):
        sub = df.dropna(subset=[metric])
        rev, dr = sub[sub.reverse_keyed][metric].values, sub[~sub.reverse_keyed][metric].values
        if len(rev) < 3 or len(dr) < 3:
            continue
        p_mw = sps.mannwhitneyu(rev, dr, alternative="two-sided").pvalue
        clusters = sub["instrument_id"].unique()
        lo = hi = np.nan
        if len(clusters) >= cfg.synthesis.min_clusters_for_cluster_bootstrap:
            groups = {c: sub[sub.instrument_id == c] for c in clusters}
            boots = []
            for _ in range(cfg.synthesis.n_boot_cluster):
                cat = pd.concat([groups[c] for c in rng.choice(clusters, len(clusters))], ignore_index=True)
                a, b = cat[cat.reverse_keyed][metric], cat[~cat.reverse_keyed][metric]
                if len(a) >= 2 and len(b) >= 2:
                    boots.append(a.mean() - b.mean())
            if len(boots) > 100:
                lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append({"metric": metric, "n_reverse": len(rev), "n_direct": len(dr), "mean_reverse": rev.mean(),
                     "mean_direct": dr.mean(), "diff_reverse_minus_direct": rev.mean() - dr.mean(),
                     "cluster_bootstrap_ci_lo": lo, "cluster_bootstrap_ci_hi": hi, "mannwhitney_p": p_mw})
    out = pd.DataFrame(rows)
    if len(out):
        out["q_bh"] = bh_adjust(out["mannwhitney_p"])
    return out


def summarize(master: pd.DataFrame) -> pd.DataFrame:
    """Item-weighted and instrument-weighted summaries of the main item metrics."""
    rows = []
    for col in ("exact_agreement", "mae_range", "spearman_rho", "ccc", "bias_range", "sd_ratio", "entropy_difference"):
        v = master[col].replace([np.inf, -np.inf], np.nan).dropna()
        per_ins = master.groupby("instrument_id")[col].mean()
        rows.append({"metric": col, "n_items": len(v), "item_weighted_mean": v.mean(), "item_weighted_median": v.median(),
                     "item_weighted_iqr_lo": v.quantile(.25), "item_weighted_iqr_hi": v.quantile(.75),
                     "item_weighted_min": v.min(), "item_weighted_max": v.max(), "n_instruments": len(per_ins),
                     "instrument_weighted_mean": per_ins.mean(), "instrument_weighted_median": per_ins.median()})
    return pd.DataFrame(rows)
