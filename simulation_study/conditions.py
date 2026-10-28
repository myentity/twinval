"""Monte Carlo design (ADEMP, Morris, White & Crowther 2019): the grid of data-generating conditions.

Every condition is ONE instrument with three 5-item scales (5 categories, factor correlation .30). The three scales share the
twin parameters of the condition, so the scale is the unit of classification and the replication is the unit of resampling.
Item parameters (human loadings, thresholds) are FIXED across conditions and replications (``item_seed``); only the twin
process, the sample size and the random seed of the replication vary.

Factors (fractional design: base condition plus one- and two-factor variations)
  N          250 / 500 / 1000
  rho        twin-human latent correlation (tracking): .30 / .50 / .70 / .80 / .90
  mu         twin mean shift (latent SD units, bias): 0 / .5 / .6 / 1
  sigma      twin factor scale (variance compression): 1 / .6
  dlambda    twin loadings differ from the human ones (metric non-invariance): no / yes
  shift      two twin items shifted by +/-.5 latent SD (scalar non-invariance): no / yes
  rare       twin never uses the highest category (rare categories, harmonisation): no / yes
  dup        twin answers items 2 and 3 almost exactly like item 1 (r = .98 before categorisation): no / yes
  k          items per scale: 5 (default) / 3 / 2 (short scales; 2-item paired models are not identified)
"""
from __future__ import annotations

ITEM_SEED = 20261027
K_ITEMS = 5
HUMAN_LOADINGS = [0.78, 0.70, 0.62, 0.74, 0.66]
TWIN_LOADINGS_DIFF = [0.45, 0.90, 0.40, 0.90, 0.45]        # |difference| about .25-.33 per item
ITEM_SHIFT = [0.5, 0.0, 0.0, 0.0, -0.5]
SCALES = ("S1", "S2", "S3")

BASE = dict(rho=0.9, mu=0.0, sigma=1.0, dlambda=False, shift=False, rare=False, dup=False, k=K_ITEMS)


def _cond(block, n, **kw):
    p = {**BASE, **kw}
    pid = "rho{rho:.3f}_mu{mu:.1f}_sg{sigma:.1f}".format(**p).replace("0_mu", "_mu") + ("_dl" if p["dlambda"] else "") + \
          ("_sh" if p["shift"] else "") + ("_rare" if p["rare"] else "") + ("_dup" if p["dup"] else "") + \
          (f"_k{p['k']}" if p["k"] != K_ITEMS else "")
    return {"block": block, "N": n, "param_id": pid, **p, "cond_id": f"{pid}_N{n}"}


def conditions() -> list[dict]:
    out = []
    for rho in (0.3, 0.5, 0.7, 0.8, 0.9):                                      # individuation (incl. base rho .9)
        for n in (250, 500, 1000):
            out.append(_cond("I_tracking", n, rho=rho))
    for rho in (0.60, 0.625, 0.65):                         # observed population r about .48 / .50 / .52 (the I cut-off)
        for n in (250, 500, 1000):
            out.append(_cond("I_threshold", n, rho=rho))
    # calibration: population C_b about .91 (bias only), .82 and .77 (just above / below .80) and .71 (strong distortion)
    for mu, sg in ((0.5, 1.0), (0.5, 0.6), (0.6, 0.6), (1.0, 1.0)):
        for n in (250, 500, 1000):
            out.append(_cond("C_calibration", n, mu=mu, sigma=sg))
    for n in (250, 500, 1000):                                                  # structure
        out.append(_cond("S_loadings", n, dlambda=True))
        out.append(_cond("S_thresholds", n, shift=True))
    out.append(_cond("R_rare", 500, rare=True))                                 # rare categories
    out.append(_cond("R_rare", 500, rare=True, mu=0.5))
    out.append(_cond("X_combined", 500, rho=0.5, dlambda=True))                 # combinations
    out.append(_cond("X_combined", 500, sigma=0.6, dlambda=True))
    out.append(_cond("P_duplicates", 500, dup=True))                           # pathological twin answers
    for k in (3, 2):                                                            # short scales
        for n in (250, 500, 1000):
            out.append(_cond("K_short", n, k=k))
    return out


def spec(cond: dict, seed: int, n: int | None = None) -> dict:
    """Scenario dictionary for ``twinval.synth.simulate`` (one instrument, three scales)."""
    k = int(cond.get("k", K_ITEMS))
    tw = {"rho": cond["rho"], "mu": cond["mu"], "sigma": cond["sigma"]}
    if cond["dlambda"]:
        tw["loading"] = TWIN_LOADINGS_DIFF
    if cond["shift"]:
        tw["item_shift"] = ITEM_SHIFT
    if cond["rare"]:
        tw["never_highest"] = True
    if cond["dup"]:
        tw["duplicate"] = [[1, 2], [1, 3]]
        tw["duplicate_corr"] = 0.98
    return {"seed": int(seed), "n_respondents": int(n or cond["N"]), "factor_correlation": 0.3, "id_prefix": "MC",
            "instruments": [{"id": 950, "key": "MC", "name": "Monte Carlo instrument", "response": [1, 5],
                             "item_seed": ITEM_SEED,
                             "scales": [{"name": s, "n_items": k, "human": {"loading": HUMAN_LOADINGS[:k]}, "twin": dict(tw)}
                                        for s in SCALES]}]}


if __name__ == "__main__":
    import pandas as pd
    c = pd.DataFrame(conditions())
    print(c.to_string())
    print(len(c), "conditions,", c.param_id.nunique(), "distinct parameter sets")
