"""Cross-level synthesis: one row per scale, links between levels, rule-based typology, results matrix.

The synthesis does **not** average different units of analysis (respondent, scale, factor, item) into a single
fidelity index. Instead it

1. joins the results of all levels in a scale-level master table (:func:`build_master`),
2. tests pre-defined families of associations between levels with a cluster (per-instrument) bootstrap
   (:func:`associations`), and
3. assigns every scale to a type with explicit, pre-specified rules on three axes - individuation (I),
   calibration (C) and structural equivalence (S) - without clustering (:func:`typology`).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from .io import Study
from .stats import bh_adjust, cluster_bootstrap_spearman


# =============================================================================== master table
def _item_aggregates(item_master: pd.DataFrame) -> pd.DataFrame:
    g = item_master.groupby("scale_id")
    agg = pd.DataFrame({
        "n_items_in_item_level": g.size(),
        "item_median_spearman_rho": g["spearman_rho"].median(),
        "item_median_ccc": g["ccc"].median(),
        "item_median_mae_range": g["mae_range"].median(),
        "item_median_abs_bias_range": g["bias_range"].apply(lambda x: x.abs().median()),
        "item_median_sd_ratio": g["sd_ratio"].median(),
        "item_pct_sd_ratio_lt1": g["sd_ratio"].apply(lambda x: (x < 1).mean()),
        "item_median_entropy_difference": g["entropy_difference"].median(),
        "item_pct_near_degenerate_twin": g["twin_near_degenerate_flag"].mean() if "twin_near_degenerate_flag" in item_master else np.nan,
    })
    if "loading_vector_distance" in item_master:
        agg["item_median_loading_vector_distance"] = g["loading_vector_distance"].median()
        agg["item_median_abs_delta_loading"] = g["abs_delta_loading_target"].median()
    return agg.reset_index()


def build_master(study: Study, scale_master: pd.DataFrame | None, structural: pd.DataFrame | None,
                 item_master: pd.DataFrame | None) -> pd.DataFrame:
    """Join scale-level, structural and item-level results (whatever is available) on ``scale_id``."""
    reg = study.registry
    base = pd.DataFrame({"scale_id": study.scale_ids})
    base["scale_name"] = [reg.scales[s].name_en for s in base.scale_id]
    base["instrument_id"] = [reg.scales[s].instrument_id for s in base.scale_id]
    base["instrument_name"] = [reg.instruments[i].name_en for i in base.instrument_id]
    base["n_items_scale"] = [reg.scales[s].n_items for s in base.scale_id]
    base["is_nonredundant"] = [not reg.scales[s].redundant for s in base.scale_id]
    m = base
    if scale_master is not None:
        keep = ["scale_id", "pearson_r", "pearson_r_ci_lo", "pearson_r_ci_hi", "pearson_q_bh_sig", "auc", "auc_ci_lo",
                "auc_ci_hi", "ccc", "cb", "u", "v", "r_dis_alpha", "delta_ccc_ls", "delta_ccc_ols", "typology_r_cb"]
        m = m.merge(scale_master[[c for c in keep if c in scale_master.columns]], on="scale_id", how="left")
        m = m.rename(columns={"pearson_r": "r", "auc": "AUC", "ccc": "CCC", "cb": "Cb", "auc_ci_lo": "AUC_ci_low",
                              "auc_ci_hi": "AUC_ci_high"})
    if structural is not None and len(structural):
        m = m.merge(structural.drop(columns=["instrument_id"]), on="scale_id", how="left")
    if item_master is not None and len(item_master):
        m = m.merge(_item_aggregates(item_master), on="scale_id", how="left")
    return m


# =============================================================================== associations
BLOCK_A = [("AUC", "r", "AUC vs Pearson r"), ("AUC", "CCC", "AUC vs CCC"), ("r", "CCC", "Pearson r vs CCC"),
           ("r_dis_alpha", "r", "Disattenuated r vs raw r"),
           ("Cb", "delta_ccc_ols", "C_b vs change of CCC after OLS recalibration")]
BLOCK_B = [("item_median_ccc", "CCC", "Median item CCC vs scale CCC"),
           ("item_median_spearman_rho", "r", "Median item rho vs scale r"),
           ("item_median_mae_range", "CCC", "Median normalised item MAE vs scale CCC"),
           ("item_median_abs_bias_range", "Cb", "Median |item bias| vs C_b"),
           ("item_median_sd_ratio", "Cb", "Median item SD ratio vs C_b"),
           ("item_pct_sd_ratio_lt1", "AUC", "Share of compressed items vs scale AUC"),
           ("item_pct_sd_ratio_lt1", "r", "Share of compressed items vs scale r"),
           ("item_pct_sd_ratio_lt1", "CCC", "Share of compressed items vs scale CCC"),
           ("item_median_entropy_difference", "AUC", "Median entropy difference vs scale AUC"),
           ("item_median_entropy_difference", "r", "Median entropy difference vs scale r"),
           ("item_median_entropy_difference", "CCC", "Median entropy difference vs scale CCC")]
BLOCK_C = [("r", "Tucker phi vs scale r"), ("AUC", "Tucker phi vs scale AUC"), ("CCC", "Tucker phi vs scale CCC"),
           ("item_median_loading_vector_distance", "Tucker phi vs item loading discrepancy"),
           ("item_median_ccc", "Tucker phi vs item CCC")]


def associations(master: pd.DataFrame, cfg, rng: np.random.Generator) -> dict[str, pd.DataFrame]:
    """Spearman associations between levels in three pre-defined families.

    * **A** within the scale level, **B** scale x item level (both with a cluster bootstrap over instruments and
      bootstrap p-values, BH-adjusted within the family);
    * **C1** (gated, primary) Tucker phi vs accuracy, restricted to factors whose Human AND Twin ESEM are admissible.
      Descriptive only when there are few instruments (no CI / BH);
    * **C2** (exploratory) the same with ALL point estimates of phi, including those from inadmissible solutions.
    """
    B = cfg.synthesis.n_boot_cluster
    mc = cfg.synthesis.min_clusters_for_cluster_bootstrap
    rows = []
    fams = [("A", [(x, y, l) for x, y, l in BLOCK_A]), ("B", [(x, y, l) for x, y, l in BLOCK_B])]
    if "tucker_phi_exploratory" in master and master["tucker_phi_exploratory"].notna().any():
        fams.append(("C2_exploratory", [("tucker_phi_exploratory", y, f"[EXPLORATORY] {l}") for y, l in BLOCK_C]))
    for name, pairs in fams:
        fam_rows = []
        for x, y, label in pairs:
            if x not in master or y not in master:
                continue
            r = cluster_bootstrap_spearman(master, x, y, "instrument_id", B, rng, mc)
            if r is not None:
                fam_rows.append({"family": name, "pair": label, "x": x, "y": y, **r})
        if fam_rows:
            df = pd.DataFrame(fam_rows)
            df["p_BH"] = bh_adjust(df["p_boot"]) if df["p_boot"].notna().any() else np.nan
            rows.append(df)
    assoc = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    out = {"associations": assoc}
    if "tucker_phi_primary" in master and master["tucker_phi_primary"].notna().sum() >= 3:
        gated = master[master["tucker_interpretable"] == True]
        c1 = []
        for y, label in BLOCK_C:
            if y not in gated:
                continue
            sub = gated.dropna(subset=["tucker_phi_primary", y])
            if len(sub) >= 3:
                rho, p = sps.spearmanr(sub["tucker_phi_primary"], sub[y])
                c1.append({"family": "C1_gated_descriptive", "pair": label.replace("Tucker phi", "Tucker phi (gated)"),
                           "x": "tucker_phi_primary", "y": y, "rho": rho, "p_naive": p, "n": len(sub),
                           "n_clusters": sub["instrument_id"].nunique()})
        out["associations_gated_descriptive"] = pd.DataFrame(c1)
        cols = ["scale_id", "scale_name", "instrument_name", "tucker_phi_primary", "tucker_phi_ci_lo", "tucker_phi_ci_hi",
                "tucker_bootstrap_stable", "r", "AUC", "CCC"]
        out["gated_scale_detail"] = gated[[c for c in cols if c in gated.columns]]

    if "metric_invariance" in master:
        rows = []
        for col, label in (("metric_invariance", "Metric invariance"), ("scalar_invariance", "Scalar invariance")):
            for met in ("r", "AUC", "CCC"):
                yes = master[master[col] == "yes"][met].dropna()
                no = master[master[col] == "no"][met].dropna()
                rec = {"status": label, "metric": met, "n_yes": len(yes), "n_no": len(no),
                       "median_yes": yes.median() if len(yes) else np.nan, "median_no": no.median() if len(no) else np.nan}
                if len(yes) >= 3 and len(no) >= 3:
                    rec["mannwhitney_p"] = float(sps.mannwhitneyu(yes, no, alternative="two-sided").pvalue)
                rows.append(rec)
        out["invariance_group_comparison"] = pd.DataFrame(rows)

    num = [c for c in ("r", "AUC", "CCC", "Cb", "r_dis_alpha", "delta_ccc_ls", "delta_ccc_ols", "tucker_phi_exploratory",
                       "tucker_phi_primary", "item_median_spearman_rho", "item_median_ccc", "item_median_mae_range",
                       "item_median_abs_bias_range", "item_median_sd_ratio", "item_median_entropy_difference",
                       "item_median_loading_vector_distance") if c in master]
    out["correlation_matrix_exploratory"] = master[num].corr(method="spearman")
    return out


# =============================================================================== typology
TYPE_ORDER = ["A", "B", "C", "D", "E", "F", "model_inadmissible", "not_tested"]


def _axis_i_detectable(row) -> str:
    if pd.isna(row["AUC"]) or pd.isna(row.get("AUC_ci_low")) or pd.isna(row["r"]) or pd.isna(row.get("pearson_q_bh_sig")):
        return "unavailable"
    ok = row["AUC"] > 0.5 and row["AUC_ci_low"] > 0.5 and row["r"] > 0 and bool(row["pearson_q_bh_sig"])
    return "positive" if ok else "weak_or_not_confirmed"


def typology(master: pd.DataFrame, cfg) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Rule-based scale types (no clustering, thresholds fixed in the configuration).

    Axes: **I** practically adequate individuation (``r >= r_practical`` and lower AUC CI > ``auc_ci_low``; the
    statistical detectability criterion is kept in a separate column), **C** calibration (``C_b >= cb_adequate``),
    **S** structural status (``structural_evidence_level``). Rules are applied in this order (first match wins)::

        model_inadmissible / not_tested                                     S is unknown - kept as separate categories
        A  I adequate  & C adequate    & S in {scalar, metric}
        B  I adequate  & C distorted                                         (evaluated before D)
        C  I inadequate               & S in {scalar, metric, configural, congruence}
        D  (I adequate or CCC >= ccc_acceptable) & S == structurally_discrepant
        E  I adequate  & C adequate    & S in {configural, congruence}
        F  everything else (defensive residual)
    """
    t = cfg.thresholds
    df = master.copy()
    need = ["r", "AUC", "AUC_ci_low", "Cb", "CCC", "structural_evidence_level"]
    missing = [c for c in need if c not in df]
    if missing:
        raise ValueError(f"typology needs scale-level and structural results (missing: {missing})")
    df["axis_I_detectable"] = df.apply(_axis_i_detectable, axis=1)
    df["axis_I_practical"] = np.where(df["r"].isna() | df["AUC_ci_low"].isna(), "unavailable",
                                      np.where((df["r"] >= t.r_practical) & (df["AUC_ci_low"] > t.auc_ci_low), "adequate", "inadequate"))
    df["axis_C"] = np.where(df["Cb"].isna(), "uncertain", np.where(df["Cb"] >= t.cb_adequate, "adequate", "distorted"))
    df["axis_S"] = df["structural_evidence_level"]

    def final(row):
        i, c, s = row["axis_I_practical"], row["axis_C"], row["axis_S"]
        if s == "model_inadmissible":
            return "model_inadmissible"
        if s == "not_tested":
            return "not_tested"
        if i == "adequate" and c == "adequate" and s in ("scalar_invariance", "metric_invariance"):
            return "A"
        if i == "adequate" and c == "distorted":
            return "B"
        if i != "adequate" and s in ("scalar_invariance", "metric_invariance", "configural_only", "congruence_supported"):
            return "C"
        if (i == "adequate" or (pd.notna(row["CCC"]) and row["CCC"] >= t.ccc_acceptable)) and s == "structurally_discrepant":
            return "D"
        if i == "adequate" and c == "adequate" and s in ("configural_only", "congruence_supported"):
            return "E"
        return "F"

    df["final_type"] = df.apply(final, axis=1)
    labels = cfg.synthesis.type_labels
    df["final_type_label"] = df["final_type"].map(lambda k: labels.get(k, k))
    df["classification_reason"] = ("I_practical=" + df["axis_I_practical"] + " (I_detectable=" + df["axis_I_detectable"] +
                                   "), C=" + df["axis_C"] + ", S=" + df["axis_S"] + " -> " + df["final_type"])
    cols = ["scale_id", "scale_name", "instrument_name", "r", "AUC", "AUC_ci_low", "CCC", "Cb", "structural_evidence_level",
            "paired_esem_admissible", "tucker_phi_primary", "tucker_phi_exploratory", "metric_invariance", "scalar_invariance",
            "axis_I_detectable", "axis_I_practical", "axis_C", "axis_S", "final_type", "final_type_label", "classification_reason"]
    df = df[[c for c in cols if c in df.columns]]
    summ = df["final_type"].value_counts().reindex(TYPE_ORDER).dropna().astype(int).rename("n").reset_index()
    summ.columns = ["final_type", "n"]
    summ["pct"] = (summ["n"] / len(df) * 100).round(1)
    summ["label"] = summ["final_type"].map(lambda k: labels.get(k, k))
    summ["example_scales"] = summ["final_type"].map(lambda k: "; ".join(df.loc[df.final_type == k, "scale_name"].head(3)))
    rules = {"version": "2.0", "thresholds": {"r_practical": t.r_practical, "auc_ci_low": t.auc_ci_low,
                                              "cb_adequate": t.cb_adequate, "ccc_acceptable": t.ccc_acceptable,
                                              "phi_congruence": t.phi_congruence},
             "decision_order": ["model_inadmissible / not_tested (S unknown; separate categories)",
                                "A: I adequate & C adequate & S in {scalar, metric}",
                                "B: I adequate & C distorted (evaluated before D)",
                                "C: I inadequate & S in {scalar, metric, configural, congruence}",
                                "D: (I adequate or CCC >= ccc_acceptable) & S == structurally_discrepant",
                                "E: I adequate & C adequate & S in {configural, congruence}",
                                "F: residual"],
             "labels": dict(labels)}
    return df, summ, rules


# =============================================================================== results matrix
def results_matrix(study: Study, profile_res: dict | None, scale_master: pd.DataFrame | None,
                   structural: pd.DataFrame | None, item_summary: pd.DataFrame | None,
                   contrasts: pd.DataFrame | None) -> pd.DataFrame:
    """Navigation table: the main result of every level with its unit of analysis (values are NOT combined)."""
    rows = []

    def add(level, unit, n, comp, metric, value, lo=np.nan, hi=np.nan, detail="", status="main"):
        rows.append({"level": level, "unit": unit, "n": n, "component": comp, "metric": metric, "value": value,
                     "ci_lo": lo, "ci_hi": hi, "detail": detail, "status": status})

    if profile_res:
        s = profile_res["identification_summary"].set_index("spec")
        if "distinctive" in s.index:
            d = s.loc["distinctive"]
            add("profile", "respondent", f"{int(d['n'])} (60 non-redundant scales)", "individuation",
                "Own-twin identification AUC (distinctive)", d["auc"], d["auc_ci_lo"], d["auc_ci_hi"])
            add("profile", "respondent", int(d["n"]), "individuation", "Top-1 identification (%)", d["top1_pct"],
                d["top1_ci_lo"], d["top1_ci_hi"], status="supplementary")
        if "profile_calibration_summary" in profile_res:
            c = profile_res["profile_calibration_summary"].set_index("metric")
            for name, lab in (("own_ccc_full", "Profile CCC, own (full)"), ("own_ccc_distinctive", "Profile CCC, own (distinctive / bias-corrected)")):
                if name in c.index:
                    add("profile", "respondent", int(study.n), "calibration/agreement", lab, c.loc[name, "value"],
                        c.loc[name, "ci_lo"], c.loc[name, "ci_hi"],
                        "same specification as the identification model (non-redundant scales, equal instrument weights)")
    if scale_master is not None:
        from .stats import fisher_z_mean
        m = scale_master
        add("scale", "scale", f"{len(m)} scales", "individuation", "Pearson r, Fisher-z mean", fisher_z_mean(m["pearson_r"].values),
            detail=f"raw arithmetic mean = {m['pearson_r'].mean():.3f}; median = {m['pearson_r'].median():.3f} "
                   f"[IQR {m['pearson_r'].quantile(.25):.3f}, {m['pearson_r'].quantile(.75):.3f}]")
        add("scale", "scale", len(m), "individuation", "Scale-level AUC, mean", m["auc"].mean(),
            detail=f"median = {m['auc'].median():.3f}")
        if "r_dis_alpha" in m:
            add("scale", "scale", int(m["r_dis_alpha"].notna().sum()), "individuation", "Disattenuated r (alpha-based, k >= 3), median",
                m["r_dis_alpha"].median(), detail=f"IQR [{m['r_dis_alpha'].quantile(.25):.3f}, {m['r_dis_alpha'].quantile(.75):.3f}]", status="supplementary")
        add("scale", "scale", len(m), "calibration/agreement", "CCC, mean", m["ccc"].mean(), detail=f"median = {m['ccc'].median():.3f}")
        add("scale", "scale", len(m), "calibration/agreement", "C_b, mean", m["cb"].mean(), detail=f"median = {m['cb'].median():.3f}")
        if "delta_ccc_ls" in m:
            add("scale", "scale", len(m), "calibration/agreement", "Delta CCC after LS recalibration, mean", m["delta_ccc_ls"].mean(),
                detail=f"median = {m['delta_ccc_ls'].median():+.3f}", status="supplementary")
    if structural is not None and len(structural):
        s = structural
        g = s["tucker_phi_primary"].dropna()
        add("structure", "factor", f"{len(g)} factors, {s.loc[s.tucker_interpretable == True, 'instrument_id'].nunique()} instruments (gated)",
            "structural fidelity", "Tucker phi, GATED (both ESEM admissible), median", g.median() if len(g) else np.nan,
            detail="descriptive; bootstrap-stable only where >= 80% replicates are admissible")
        e = s["tucker_phi_exploratory"].dropna()
        add("structure", "factor", f"{len(e)} factors (all point estimates)", "structural fidelity",
            "Tucker phi, EXPLORATORY (incl. inadmissible ESEM), median", e.median() if len(e) else np.nan,
            detail="NOT used for typology or main conclusions", status="supplementary")
        for col, lab in (("metric_invariance", "Metric invariance (Millsap/Liu)"), ("scalar_invariance", "Scalar invariance (Wu chain)")):
            vc = s[col].value_counts()
            add("structure", "scale", len(s), "structural fidelity", lab, float(vc.get("yes", 0)),
                detail=f"yes: {vc.get('yes', 0)}, no: {vc.get('no', 0)}, not_tested: {vc.get('not_tested', 0)}")
    if item_summary is not None and len(item_summary):
        i = item_summary.set_index("metric")
        for met, lab, comp in (("spearman_rho", "Item-level Spearman rho, median", "individuation"),
                               ("ccc", "Item-level CCC, median", "calibration/agreement"),
                               ("sd_ratio", "Item-level SD ratio, median", "calibration/agreement")):
            if met in i.index:
                r = i.loc[met]
                add("item", "item", int(r["n_items"]), comp, lab, r["item_weighted_median"],
                    detail=f"mean = {r['item_weighted_mean']:.3f}; IQR [{r['item_weighted_iqr_lo']:.3f}, {r['item_weighted_iqr_hi']:.3f}]")
    if contrasts is not None and len(contrasts):
        a = contrasts[contrasts["instrument_id"].astype(str) == "ALL_INSTRUMENT_WEIGHTED"]
        if len(a):
            r = a.iloc[0]
            add("item", "instrument", r["instrument_name"], "mechanism", "Delta %var item (twin - human), pp",
                r["delta_item_point"], r["delta_item_ci_lo"], r["delta_item_ci_hi"], "ANOVA-style decomposition, paired bootstrap contrast")
    return pd.DataFrame(rows)
