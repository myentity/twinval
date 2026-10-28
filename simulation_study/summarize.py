"""Step 3: performance measures of the Monte Carlo study (DESIGN.md, section P).

Reads replications/*.csv and population/population_truth.csv; writes summary/*.csv and summary/summary.md.

Conventions
* A replication counts if it finished; failed replications (error rows without a later success) are counted and reported.
* Non-estimable models stay in every main denominator; accuracy among estimable scales is a separate, labelled measure.
* Monte Carlo standard error (MCSE): the three scales of a replication share one data set, so a proportion is first
  averaged within the replication and MCSE = SD(replication means) / sqrt(R). Pooled proportions over several conditions
  are ratio estimates with replications as clusters.
* Conditions at the I cut-off (population r about .50) have no defined population I axis; they enter no sensitivity /
  specificity and are reported as P(I adequate).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from conditions import conditions  # noqa: E402

REP = HERE / "replications"
OUT = HERE / "summary"
SUBSTANTIVE = {"scalar_invariance", "metric_invariance", "configural_fit", "congruence_supported", "structurally_discrepant"}
BOUNDARY_R = 0.005                                      # |r_pop - .50| below this: population I axis undefined


def cb_is_one(c: dict) -> bool:
    """Generator gives identical Human and Twin score distributions, so the true C_b is exactly 1 (the boundary)."""
    return c["mu"] == 0 and c["sigma"] == 1 and not (c["dlambda"] or c["shift"] or c["rare"] or c["dup"])


CB_ONE = {c["cond_id"] for c in conditions() if cb_is_one(c)}


def step_class(x: str) -> str:
    if x in ("supported", "rejected"):
        return x
    if isinstance(x, str) and x.startswith("not_evaluable"):
        return "not_evaluable"
    if x == "decision_not_evaluable":
        return "not_evaluable"
    return "not_interpreted"


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = [pd.read_csv(f) for f in sorted(REP.glob("*.csv"))]
    d = pd.concat(frames, ignore_index=True)
    ok = d[d.error.isna() | (d.error == "")].copy()
    ok = ok.drop_duplicates(["cond_id", "rep", "scale"], keep="last")
    failed = d[~(d.error.isna() | (d.error == ""))][["cond_id", "rep", "error"]]
    failed = failed[~failed.set_index(["cond_id", "rep"]).index.isin(ok.set_index(["cond_id", "rep"]).index)]
    cond = pd.DataFrame(conditions())
    truth = pd.read_csv(HERE / "population" / "population_truth.csv")
    ok = ok.merge(cond[["cond_id", "param_id", "block"]], on="cond_id").merge(truth, on=["param_id", "scale"], how="left",
                                                                             suffixes=("", "_pop"))
    ok["boundary"] = (ok.r_pop - 0.50).abs() < BOUNDARY_R
    return ok, failed


def mc_prop(df: pd.DataFrame, ind: pd.Series) -> tuple[float, float, int]:
    """Proportion with MCSE from replication means (replications as clusters, scales averaged within)."""
    g = ind.groupby([df.cond_id, df.rep]).mean()
    R = len(g)
    return float(g.mean()), float(g.std(ddof=1) / np.sqrt(R)) if R > 1 else np.nan, R


def ratio(df: pd.DataFrame, num: pd.Series, den: pd.Series) -> tuple[float, float, int]:
    """Pooled ratio sum(num)/sum(den) with a cluster (replication) standard error."""
    keys = [df.cond_id, df.rep]
    y, n = num.astype(float).groupby(keys).sum(), den.astype(float).groupby(keys).sum()
    y, n = y[n > 0], n[n > 0]
    if n.sum() == 0:
        return np.nan, np.nan, 0
    p = y.sum() / n.sum()
    k = len(n)
    se = np.sqrt(((y - p * n) ** 2).sum() * k / (k - 1)) / n.sum() if k > 1 else np.nan
    return float(p), float(se), int(n.sum())


def per_condition(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cid, g in d.groupby("cond_id", sort=False):
        r = {"cond_id": cid, "block": g.block.iloc[0], "N": g.N.iloc[0], "R": g.rep.nunique(),
             "r_pop": g.r_pop.mean(), "cb_pop": g.cb_pop.mean(), "S_design": g.S_design.iloc[0],
             "S_ref": "/".join(sorted(set(g.S_ref))), "type_ref": "/".join(sorted(set(g.type_ref)))}
        for v, pop in (("r", "r_pop"), ("ccc", "ccc_pop"), ("cb", "cb_pop"), ("auc", "auc_pop")):
            e = g[v] - g[pop]
            r[f"{v}_bias"], r[f"{v}_rmse"] = e.mean(), np.sqrt((e ** 2).mean())
            if v == "cb" and cid in CB_ONE:
                # true C_b = 1 exactly: the estimate (and every bootstrap value) is < 1 unless means and SDs coincide, so
                # a percentile interval cannot contain the true value; coverage is not reported (NaN). Conditions with a
                # true C_b strictly below 1, however close, keep their coverage.
                r[f"{v}_coverage"] = r[f"{v}_coverage_mcse"] = np.nan
                continue
            r[f"{v}_coverage"], r[f"{v}_coverage_mcse"], _ = mc_prop(g, (g[f"{v}_ci_lo"] <= g[pop]) & (g[pop] <= g[f"{v}_ci_hi"]))
        if g.boundary.all():
            r["P_I_adequate"], r["P_I_adequate_mcse"], _ = mc_prop(g, g.I_rec == "adequate")
        else:
            r["I_correct"], r["I_correct_mcse"], _ = mc_prop(g, g.I_rec == g.I_pop)
        r["C_correct"], r["C_correct_mcse"], _ = mc_prop(g, g.C_rec == g.C_pop)
        for st in ("admissible", "nonconverged", "identification_not_established", "inadmissible"):
            r[f"W0_{st}"] = (g.baseline_status == st).mean()
        r["W0_misfit"] = (g.baseline_fit == "not_supported").mean()
        r["near_boundary_flag"] = g.cfa_flags.fillna("").str.contains("near_boundary").mean()
        r["paired_esem_admissible"] = g.paired_esem_admissible.mean()
        for stp in ("W1", "W2", "W3", "M1"):
            cls = g[stp].map(step_class)
            for c in ("supported", "rejected", "not_evaluable", "not_interpreted"):
                r[f"{stp}_{c}"] = (cls == c).mean()
        r["S_available"], r["S_available_mcse"], _ = mc_prop(g, g.S_rec.isin(SUBSTANTIVE))
        if g.S_design.notna().all():
            r["S_eq_design"], r["S_eq_design_mcse"], _ = mc_prop(g, g.S_rec == g.S_design)
        r["S_eq_ref"], r["S_eq_ref_mcse"], _ = mc_prop(g, g.S_rec == g.S_ref)
        if not g.boundary.all():
            if g.type_design.notna().all():
                r["type_eq_design"], r["type_eq_design_mcse"], _ = mc_prop(g, g.type_rec == g.type_design)
            r["type_eq_ref"], r["type_eq_ref_mcse"], _ = mc_prop(g, g.type_rec == g.type_ref)
            est = g.baseline_status == "admissible"
            if est.any():
                r["type_eq_ref_given_estimable"] = (g.type_rec[est] == g.type_ref[est]).mean()
        r["D_rate"], r["F_rate"] = (g.type_rec == "D").mean(), (g.type_rec == "F").mean()
        r["harmonised"] = (g.n_items_harmonised > 0).mean()
        r["min_categories_after"] = g.min_categories_after.min()
        r["seconds_per_rep"] = g.groupby("rep").seconds.first().mean()
        rows.append(r)
    return pd.DataFrame(rows)


def d_f_causes(d: pd.DataFrame) -> pd.DataFrame:
    x = d[d.type_rec.isin(["D", "F"])].copy()
    if not len(x):
        return pd.DataFrame(columns=["cond_id", "type_rec", "cause", "n"])

    def cause(r):
        if r.S_rec == "structurally_discrepant":
            return "ESEM discrepancy after loss of CFA conclusion (" + str(r.cfa_no_conclusion_reason) + ")"
        if r.type_rec == "F":
            return "combination outside the rules (I " + r.I_rec + ", C " + r.C_rec + ", S " + r.S_rec + ")"
        return "other"
    x["cause"] = x.apply(cause, axis=1)
    return x.groupby(["cond_id", "type_rec", "cause"]).size().rename("n").reset_index()


def sens_spec(d: pd.DataFrame, ref_col: str) -> pd.DataFrame:
    x = d[~d.boundary & d[ref_col].notna()]
    rows = []
    for k in sorted(set(x[ref_col]) | set(x.type_rec)):
        for n in sorted(x.N.unique()) + ["all"]:
            g = x if n == "all" else x[x.N == n]
            pos, neg = g[ref_col] == k, g[ref_col] != k
            se, se_mc, n_pos = ratio(g, pos & (g.type_rec == k), pos)
            sp, sp_mc, n_neg = ratio(g, neg & (g.type_rec != k), neg)
            rows.append({"type": k, "N": n, "sensitivity": se, "sensitivity_mcse": se_mc, "n_scales_pos": n_pos,
                         "specificity": sp, "specificity_mcse": sp_mc, "n_scales_neg": n_neg})
    return pd.DataFrame(rows)


STEP_GROUPS = {                                         # generator -> conditions pooled with equal weight, per N
    "no violation, k = 5": ("I_tracking", "I_threshold", "C_calibration"),
    "intercept shift": ("S_thresholds",),
    "loading difference": ("S_loadings",),
    "no violation, k = 3": ("K_short_k3",),
    "no violation, k = 2": ("K_short_k2",),
}


def step_rates(d: pd.DataFrame) -> pd.DataFrame:
    """Wu step decisions by generator and N: false rejections of true equalities and detection of the built-in ones.

    Conditions are pooled with equal weight; MCSE from replication means (replications as clusters)."""
    blk = d.block.where(d.block != "K_short", "K_short_k" + d.k.astype(str)) if "k" in d else d.block
    rows = []
    for name, blocks in STEP_GROUPS.items():
        for n in sorted(d.N.unique()):
            g = d[blk.isin(blocks) & (d.N == n)]
            if not len(g):
                continue
            r = {"generator": name, "N": n, "n_conditions": g.cond_id.nunique()}
            for lab, ind in (("W0_admissible", g.baseline_status == "admissible"),
                             ("W1_rejected", g.W1.map(step_class) == "rejected"),
                             ("W2_rejected", g.W2.map(step_class) == "rejected"),
                             ("W3_rejected", g.W3.map(step_class) == "rejected"),
                             ("S_eq_design", g.S_rec == g.S_design)):
                r[lab], r[f"{lab}_mcse"], _ = mc_prop(g, ind)
            rows.append(r)
    return pd.DataFrame(rows)


def order_check() -> pd.DataFrame:
    """Row-order subset: refit candidates and order-sensitive flags (false alarms), per scale and per replication.

    The three scales of a replication share one data set, so the upper bound is computed at the replication level: event
    "at least one scale of the replication flagged", among replications with at least one refit candidate (exact one-sided
    95% Clopper-Pearson bound for zero events, 1 - .05^(1/n)).
    """
    files = sorted((HERE / "replications_order_check").glob("*.csv"))
    if not files:
        return pd.DataFrame()
    o = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    o = o[o.error.isna() | (o.error == "")]
    rep = o.groupby(["cond_id", "rep"]).agg(c=("order_candidate", "sum"), s=("order_sensitive", "sum"))
    g, gr = o.groupby("cond_id"), rep.groupby("cond_id")
    out = pd.DataFrame({"R": g.rep.nunique(), "scales": g.size(), "candidates": g.order_candidate.sum(),
                        "order_sensitive": g.order_sensitive.sum(), "reps_with_candidates": gr.c.apply(lambda x: (x > 0).sum()),
                        "reps_flagged": gr.s.apply(lambda x: (x > 0).sum())}).reset_index()
    tot = out.drop(columns="cond_id").sum()
    out = pd.concat([out, pd.DataFrame([{"cond_id": "ALL", **tot.to_dict()}])], ignore_index=True)
    n, x = out.reps_with_candidates, out.reps_flagged
    out["rep_rate_upper95"] = np.where((x == 0) & (n > 0), 1 - 0.05 ** (1 / n.clip(lower=1)), np.nan)
    return out


def main():
    OUT.mkdir(exist_ok=True)
    d, failed = load()
    if "k" not in d:
        d = d.merge(pd.DataFrame(conditions())[["cond_id", "k"]], on="cond_id")
    sr = step_rates(d)
    sr.to_csv(OUT / "step_rates_by_generator.csv", index=False)
    oc = order_check()
    oc.to_csv(OUT / "order_check.csv", index=False)
    pc = per_condition(d)
    pc.to_csv(OUT / "condition_measures.csv", index=False)
    tr_ref = pd.crosstab([d.block, d.type_ref], d.type_rec)
    tr_ref.to_csv(OUT / "type_transitions_reference.csv")
    tr_des = pd.crosstab([d.block, d.type_design.fillna("no expectation")], d.type_rec)
    tr_des.to_csv(OUT / "type_transitions_design.csv")
    s_ref = pd.crosstab([d.block, d.S_ref], d.S_rec)
    s_ref.to_csv(OUT / "S_transitions_reference.csv")
    ss_ref, ss_des = sens_spec(d, "type_ref"), sens_spec(d, "type_design")
    ss_ref.to_csv(OUT / "sens_spec_reference.csv", index=False)
    ss_des.to_csv(OUT / "sens_spec_design.csv", index=False)
    # the N = 500 grid has extra conditions (rare categories, combinations, duplicates): compare N on the common set too
    n_by_param = d.groupby("param_id").N.nunique()
    common = d.param_id.isin(n_by_param[n_by_param == d.N.nunique()].index)
    sens_spec(d[common], "type_design").to_csv(OUT / "sens_spec_design_common.csv", index=False)
    te = pc[pc.type_eq_design.notna()].assign(common=lambda x: x.cond_id.isin(set(d.cond_id[common])))
    pd.concat([te.groupby("N").type_eq_design.agg(["mean", "count"]).assign(set="all"),
               te[te.common].groupby("N").type_eq_design.agg(["mean", "count"]).assign(set="common")]).reset_index() \
        .rename(columns={"mean": "type_eq_design_mean", "count": "n_conditions"}).to_csv(OUT / "type_agreement_by_N.csv", index=False)
    df_c = d_f_causes(d)
    df_c.to_csv(OUT / "D_F_causes.csv", index=False)
    failed.to_csv(OUT / "failed_replications.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["cond_id", "R", "r_pop", "cb_pop", "I_correct", "P_I_adequate", "C_correct", "W0_admissible", "S_available",
            "S_eq_design", "S_eq_ref", "type_eq_design", "type_eq_ref", "D_rate", "F_rate", "seconds_per_rep"]
    lines = ["# Monte Carlo summary", "",
             f"Replications: {d.groupby('cond_id').rep.nunique().sum()} finished in {d.cond_id.nunique()} conditions; "
             f"failed (no later success): {failed[['cond_id', 'rep']].drop_duplicates().shape[0]}", "",
             "## Per condition (proportions over scales x replications)", "", pc[[c for c in cols if c in pc]].round(3).to_markdown(index=False), "",
             "## Final type: reference -> recovered", "", tr_ref.to_markdown(), "",
             "## Sensitivity / specificity vs reference type (boundary conditions excluded)", "",
             ss_ref.round(3).to_markdown(index=False), "",
             "## Wu step decisions by generator (conditions pooled with equal weight)", "", sr.round(3).to_markdown(index=False), "",
             "## Row-order check (false alarms)", "", oc.to_markdown(index=False) if len(oc) else "not run", "",
             "## D / F with cause", "", df_c.to_markdown(index=False) if len(df_c) else "none", ""]
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[:8]))


if __name__ == "__main__":
    main()
