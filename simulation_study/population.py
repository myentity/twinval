"""Step 1: the three bases of comparison for every parameter set of the Monte Carlo grid.

1. Generator properties (known by construction): equality of loadings, threshold spacing and item intercepts between human
   and twin, local dependence (duplicated answers), number of items. The expected structural status of the PRIMARY
   Wu/Estabrook sequence is derived from them by an explicit mapping (``design_expectation``).
2. Population score values: r, C_b, CCC, AUC of the scale scores on a very large independent sample.
3. Large-sample reference classification: twinval's rules applied to one sample of N = 20,000. This is the
   operational verdict of the procedure itself with nearly unlimited data, NOT an independent truth.

Scale level (I and C axes): Pearson r, C_b, CCC and the scale-level AUC of the scale scores from one very large sample
(N = 200,000; the AUC from 2 million random human-twin pairs). With r > 0 the population AUC exceeds .5, so the population
I axis is r >= .50 and the C axis is C_b >= .80.

Structural level (S axis): twinval's structural decisions (paired CFA, ESEM, gate) on one large sample
(N = 20,000). With this N the change-in-fit indices and the absolute fit are essentially their population values, so the
resulting S status is the verdict the rules give with unlimited data ("population verdict"). The design intent (invariant,
loadings differ, thresholds shifted) is kept next to it.

Writes population/population_truth.csv (one row per parameter set and scale) and population/population_by_condition.csv.
"""
from __future__ import annotations

import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from conditions import conditions  # noqa: E402
import mc_core as M  # noqa: E402

OUT = HERE / "population"
WORK = HERE / "work_population"
N_SCORES, N_STRUCT, N_PAIRS = 200_000, 20_000, 2_000_000
SEED_SCORES, SEED_STRUCT = 7_000_001, 7_000_002


def design_expectation(c) -> dict:
    """Generator properties and the Wu/Estabrook status they imply.

    Wu/Estabrook steps as fitted (theta parameterisation, ID.cat = Wu.Estabrook.2016):
      W0 configural  - thresholds free on both sides; item intercepts fixed at 0 and residual variances at 1 on both sides
      W1 thresholds  - thresholds equal; twin intercepts and residual variances freed, so W1 tests only the SPACING of the
                       thresholds (with 5 categories: 4 equalities - 2 freed parameters = 2 constraints per item)
      W2 loadings    - W1 + equal loadings
      W3 scalar      - W2 + equal intercepts
    Generator -> expected decisions:
      no violation                          W1, W2, W3 supported            -> scalar_invariance
      item_shift (+/-.5 on items 1 and 5)   location of x* moves; absorbed by the freed twin intercept at W1, so W1 and W2
                                            hold and W3 (equal intercepts) fails -> metric_invariance
      twin loadings differ                  x* keeps unit variance, so threshold spacing is unchanged (W1 holds) and W2
                                            fails -> configural_fit
      2-item scales                         paired W0 has negative df -> identification not established -> model_not_estimable
      duplicated twin answers               local dependence: the one-factor twin model is misspecified; no Wu status is
                                            implied (estimability depends on the data) -> no expectation
    Whether a built-in difference also exceeds the change-in-fit cut-offs is an empirical question; the reference
    classification (N = 20,000) answers it for the procedure as implemented.
    """
    props = {"loadings_equal": not c["dlambda"], "intercepts_equal": not c["shift"],
             "threshold_spacing_equal": True, "local_dependence": bool(c["dup"]), "n_items": int(c["k"]),
             "twin_never_highest": bool(c["rare"])}
    if c["dup"]:
        exp = None
    elif c["k"] <= 2:
        exp = "model_not_estimable"
    elif c["dlambda"]:
        exp = "configural_fit"
    elif c["shift"]:
        exp = "metric_invariance"
    else:
        exp = "scalar_invariance"
    return {**props, "S_design": exp}


def score_truth(c) -> pd.DataFrame:
    wd = WORK / f"{c['param_id']}_scores"
    st = M.generate(c, SEED_SCORES, wd, n=N_SCORES)
    H, T = st.matrices()
    rng = np.random.default_rng(1)
    rows = []
    for k, sid in enumerate(st.scale_ids):
        h, t = H[:, k], T[:, k]
        r = float(np.corrcoef(h, t)[0, 1])
        v = t.std(ddof=1) / h.std(ddof=1)
        u = (t.mean() - h.mean()) / np.sqrt(h.std(ddof=1) * t.std(ddof=1))
        cb = 2 / (v + 1 / v + u ** 2)
        zh, zt = (h - h.mean()) / h.std(ddof=1), (t - h.mean()) / h.std(ddof=1)
        i = rng.integers(0, len(h), N_PAIRS)
        j = (i + rng.integers(1, len(h), N_PAIRS)) % len(h)
        own, other = np.abs(zh[i] - zt[i]), np.abs(zh[i] - zt[j])
        auc = float(np.mean(own < other) + 0.5 * np.mean(own == other))
        rows.append({"param_id": c["param_id"], "scale": sid.split("_")[1], "r_pop": r, "cb_pop": cb, "ccc_pop": r * cb,
                     "auc_pop": auc, "sd_ratio_pop": v, "u_pop": u})
    M.cleanup(wd)
    return pd.DataFrame(rows)


def structural_truth(c) -> pd.DataFrame:
    wd = WORK / f"{c['param_id']}_struct"
    M.cleanup(wd)
    t0 = time.time()
    st = M.generate(c, SEED_STRUCT, wd / "sim", n=N_STRUCT)
    res = M.structural(st, wd / "r")
    tab = res["table"].set_index("scale_id")
    s = res["status"].set_index("scale_id")
    rows = []
    for sid in st.scale_ids:
        rows.append({"param_id": c["param_id"], "scale": sid.split("_")[1], "S_ref": tab.loc[sid, "structural_evidence_level"],
                     "cfa_conclusion_ref": s.loc[sid, "cfa_conclusion"], "baseline_status_ref": s.loc[sid, "baseline_status"],
                     "baseline_fit_ref": s.loc[sid, "baseline_fit"], "threshold_ref": s.loc[sid, "threshold_decision"],
                     "metric_ref": s.loc[sid, "metric_decision"], "scalar_ref": s.loc[sid, "scalar_decision"],
                     "millsap_metric_ref": s.loc[sid, "sensitivity_millsap_metric_decision"],
                     "phi_ref": tab.loc[sid, "tucker_phi_exploratory"], "esem_admissible_ref": tab.loc[sid, "paired_esem_admissible"],
                     "flags_ref": s.loc[sid, "cfa_flags"], "seconds": round(time.time() - t0, 1)})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)
    conds = pd.DataFrame(conditions())
    psets = conds.drop_duplicates("param_id").to_dict("records")
    with ProcessPoolExecutor(max_workers=7) as ex:
        sc = pd.concat(ex.map(score_truth, psets))
        stt = pd.concat(ex.map(structural_truth, psets))
    truth = sc.merge(stt, on=["param_id", "scale"])
    pmeta = conds.drop_duplicates("param_id").set_index("param_id")[["rho", "mu", "sigma", "dlambda", "shift", "rare", "dup", "k"]]
    truth = truth.merge(pmeta, left_on="param_id", right_index=True)
    design = pd.DataFrame([{"param_id": p, **design_expectation(pmeta.loc[p])} for p in pmeta.index])
    truth = truth.merge(design, on="param_id")
    # population score axes (basis 2)
    truth["I_pop"] = np.where(truth.r_pop >= M.CFG.thresholds.r_practical, "adequate", "inadequate")
    truth["C_pop"] = truth.cb_pop.map(M.axis_C)
    # expected type: population I and C with the generator-implied S (basis 1 + 2); reference type: with S_ref (basis 2 + 3)
    truth["type_design"] = [M.final_type(i, c, s) if isinstance(s, str) else None
                            for i, c, s in zip(truth.I_pop, truth.C_pop, truth.S_design)]
    truth["type_ref"] = [M.final_type(i, c, s) for i, c, s in zip(truth.I_pop, truth.C_pop, truth.S_ref)]
    truth["S_design_matches_ref"] = truth.S_design == truth.S_ref
    truth.to_csv(OUT / "population_truth.csv", index=False)
    agg = truth.groupby("param_id").agg(
        r_pop=("r_pop", "mean"), cb_pop=("cb_pop", "mean"), auc_pop=("auc_pop", "mean"), ccc_pop=("ccc_pop", "mean"),
        S_design=("S_design", "first"), S_ref=("S_ref", lambda x: "/".join(sorted(set(x)))),
        type_design=("type_design", "first"), type_ref=("type_ref", lambda x: "/".join(sorted(set(x)))))
    by = conds.merge(agg, on="param_id")
    by.to_csv(OUT / "population_by_condition.csv", index=False)
    pd.set_option("display.width", 250)
    print(agg.round(3).to_string())


if __name__ == "__main__":
    main()
