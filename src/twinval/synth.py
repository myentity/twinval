"""Synthetic paired data with a KNOWN twin process ("twinval simulate").

Purpose: show every kind of result the framework can produce (all types of the typology, model-inadmissible and
not-tested scales, degenerate/robust cluster bootstrap, leave-one-instrument-out, harmonisation, ...) on data that carry no
personal information, and let users check that the pipeline recovers what was built into the data.

Generative model (all latent variables standardised)

* **Human.** An instrument has one factor per (non-redundant) scale; factors correlate with ``factor_correlation``.
  Item response  x*_ij = lambda_j * eta_i,f + sqrt(1 - lambda_j^2) * e_ij,  cut into ``K`` ordered categories by fixed
  item thresholds (the same thresholds are used for the twin, so a shift of the twin is a real shift).
* **Twin.** The twin factor is  eta_T = mu + sigma * ( sqrt(1-h^2) * (rho * eta_H + sqrt(1-rho^2) * zeta) + h * G )  where

  ========== ==========================================================================================
  ``rho``    latent correlation with the human factor  -> *individuation* (Pearson r of the scale scores)
  ``mu``     mean shift (in SD units)                   -> *calibration* (C_b), bias
  ``sigma``  scale of the twin factor (<1: compression) -> *calibration* (C_b), SD ratio
  ``h``      share of an instrument-wide "halo" factor G shared by all scales of the instrument
  ========== ==========================================================================================

  Twin items use their own loadings ``loading`` (default: the human ones) - unequal loadings break *metric invariance*
  (structural evidence "configural only"), loadings close to 1 make all items of the scale nearly identical (twin
  "answers the same to every item"), ``duplicate`` makes the twin answer some items almost identically (provokes
  non-convergence / Heywood cases; three or more exactly identical twin items make the scale-level CFA and the twin ESEM inadmissible), ``resid_corr`` lets the twin reproduce the
  human's item-specific residuals, ``item_shift`` moves single twin items up/down (breaks *scalar* but not metric invariance), ``copy_items`` lets the twin reproduce the human's answer verbatim on some items (leakage), ``near_constant`` makes the twin give the same answer to (almost) everybody on an item. ``reverse_loading_scale`` shrinks the
  loadings of reverse-keyed items, ``never_lowest`` / ``never_highest`` remove an extreme category from the twin's answers
  (triggers harmonisation).
* Data are stored in the framework's format (reverse-keyed items already reverse-scored).

The intended type of every scale is written to ``scenario_truth.csv``; ``compare_with_truth`` compares it with the typology
recovered by the pipeline. Intended types are *designs*, not guarantees: a scale built as type "E" can end up "F" if
sampling noise pushes a statistic over a threshold - that is exactly the kind of variation the example shows.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import norm


# =============================================================================== helpers
def _thresholds(K: int, n_items: int, rng, centre=0.6, spread=0.25, jitter=0.25) -> np.ndarray:
    """(n_items, K-1) thresholds of a unimodal category distribution centred at ``centre`` of the range, jittered per item."""
    ks = np.arange(1, K + 1)
    w = np.exp(-0.5 * ((ks - (1 + centre * (K - 1))) / (spread * K)) ** 2)
    base = norm.ppf(np.cumsum(w / w.sum())[:-1])
    return base[None, :] + rng.normal(0, jitter, (n_items, 1))


def _categorise(x: np.ndarray, thr: np.ndarray, lo: int) -> np.ndarray:
    out = np.zeros(x.shape, dtype=int)
    for k in range(thr.shape[1]):
        out += (x > thr[:, k][None, :])
    return out + lo


def _vec(v, n, default):
    if v is None:
        return np.array(default, float)
    a = np.asarray(v, float)
    return np.full(n, float(a)) if a.ndim == 0 else a


# =============================================================================== generator
def simulate(spec: dict, out_dir: str | Path) -> dict:
    """Generate registry, paired data and scenario truth from a scenario dictionary; returns paths and the truth table."""
    out = Path(out_dir)
    data_dir = out / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    seed0 = int(spec.get("seed", 1))
    N = int(spec.get("n_respondents", 500))
    prefix = spec.get("id_prefix", "SY")
    ids = [f"{prefix}{i + 1:04d}" for i in range(N)]
    fcorr = float(spec.get("factor_correlation", 0.3))

    reg = {"schema_version": 1, "instruments": {}}
    truth = []
    for ins in spec["instruments"]:
        iid = int(ins["id"])
        # every instrument has its own random stream (base seed, instrument id, optional instrument "seed"), so changing one
        # instrument never changes the data of another
        rng = np.random.default_rng([seed0, iid, int(ins.get("seed", 0))])
        lo, hi = ins.get("response", [1, 5])
        K = hi - lo + 1
        scales = ins["scales"]
        F = len(scales)
        phi = (1 - fcorr) * np.eye(F) + fcorr * np.ones((F, F))
        eta_h = rng.multivariate_normal(np.zeros(F), phi, size=N)
        G = rng.normal(size=N)                                          # instrument-wide halo factor of the twin
        halo0 = float(ins.get("twin", {}).get("halo", 0.0))
        n_items_total = sum(int(s["n_items"]) for s in scales)
        thr_all = _thresholds(K, n_items_total, rng)
        node = {"key": ins.get("key", f"SYN{iid}"), "name_en": ins["name"], "name_ru": ins["name"],
                "item_prefix": ins.get("key", f"SYN{iid}"), "reference": "synthetic (twinval simulate)",
                "response_min": int(lo), "response_max": int(hi), "scales": []}
        item_no, sids, sums_h, sums_t = 0, [], [], []
        for f, s in enumerate(scales):
            k = int(s["n_items"])
            sid = f"{iid}_{s['name'].replace(' ', '')}"
            items = [f"{node['item_prefix']}_{item_no + j + 1}" for j in range(k)]
            rev_pos = [int(p) for p in s.get("reverse", [])]                 # 1-based positions inside the scale
            thr = thr_all[item_no:item_no + k]
            hu, tw = s.get("human", {}), s.get("twin", {})
            lam_h = np.clip(_vec(hu.get("loading"), k, rng.uniform(0.62, 0.78, k)), 0.05, 0.99)
            lam_t = _vec(tw.get("loading"), k, lam_h)
            for p in rev_pos:
                lam_t[p - 1] *= float(tw.get("reverse_loading_scale", 1.0))
            lam_t = np.clip(lam_t, 0.02, 0.985)
            rho, mu, sg = float(tw.get("rho", 0.85)), float(tw.get("mu", 0.0)), float(tw.get("sigma", 1.0))
            h = float(tw.get("halo", halo0))
            zeta = rng.normal(size=N)
            eta_t = mu + sg * (np.sqrt(1 - h ** 2) * (rho * eta_h[:, f] + np.sqrt(1 - rho ** 2) * zeta) + h * G)
            e_h = rng.normal(size=(N, k))
            re = float(tw.get("resid_corr", 0.0))                           # twin reproduces the human item-specific residual
            e_t = re * e_h + np.sqrt(1 - re ** 2) * rng.normal(size=(N, k))
            xh = eta_h[:, [f]] * lam_h + e_h * np.sqrt(1 - lam_h ** 2)
            xt = eta_t[:, None] * lam_t + e_t * np.sqrt(1 - lam_t ** 2) + _vec(tw.get("item_shift"), k, np.zeros(k))
            for i_, j_ in tw.get("duplicate", []):                          # twin answers item j almost like item i
                a = float(tw.get("duplicate_corr", 0.95))
                xt[:, j_ - 1] = a * xt[:, i_ - 1] + np.sqrt(1 - a ** 2) * rng.normal(size=N)
            H, T = _categorise(xh, thr, lo), _categorise(xt, thr, lo)
            for j_ in tw.get("copy_items", []):                              # twin reproduces the human answer exactly (leakage)
                T[:, j_ - 1] = H[:, j_ - 1]
            for j_ in tw.get("near_constant", []):                           # twin gives (almost) everybody the same answer
                keep = rng.random(N) < float(tw.get("near_constant_keep", 0.015))
                T[:, j_ - 1] = np.where(keep, T[:, j_ - 1], int(round(np.median(T[:, j_ - 1]))))
            if tw.get("never_lowest"):
                T[T == lo] = lo + 1
            if tw.get("never_highest"):
                T[T == hi] = hi - 1
            df = pd.DataFrame({"ID": ids})
            for j in range(k):
                df[f"it{j + 1}_h"] = H[:, j]
            for j in range(k):
                df[f"it{j + 1}_t"] = T[:, j]
            df.to_csv(data_dir / f"{sid}_paired_wide.csv", index=False)
            node["scales"].append({"id": sid, "name_en": s["name"], "name_ru": s["name"], "legacy_id": sid,
                                   "level": s.get("level", "scale"), "scoring": s.get("scoring", "mean"), "items": items,
                                   "reverse_items": [items[p - 1] for p in rev_pos], "parent": None, "redundant": False})
            truth.append({"scale_id": sid, "instrument_id": iid, "scale_name": s["name"], "intended_type": s.get("intended", ""),
                          "design": s.get("design", ""), "rho": rho, "mu": mu, "sigma": sg, "halo": h,
                          "twin_loading_mean": float(np.mean(lam_t)), "human_loading_mean": float(np.mean(lam_h)), "n_items": k})
            sids.append(sid)
            sums_h.append(H)
            sums_t.append(T)
            item_no += k
        if ins.get("total"):                                                # a total score re-using all items (redundant)
            t = ins["total"]
            H, T = np.column_stack(sums_h), np.column_stack(sums_t)
            sid = f"{iid}_{t['name'].replace(' ', '')}"
            df = pd.DataFrame({"ID": ids})
            for j in range(H.shape[1]):
                df[f"it{j + 1}_h"] = H[:, j]
            for j in range(T.shape[1]):
                df[f"it{j + 1}_t"] = T[:, j]
            df.to_csv(data_dir / f"{sid}_paired_wide.csv", index=False)
            allitems = [i for sc in node["scales"] for i in sc["items"]]
            allrev = [i for sc in node["scales"] for i in sc["reverse_items"]]
            node["scales"].append({"id": sid, "name_en": t["name"], "name_ru": t["name"], "legacy_id": sid, "level": "total",
                                   "scoring": t.get("scoring", "sum"), "items": allitems, "reverse_items": allrev,
                                   "parent": None, "redundant": True})
            truth.append({"scale_id": sid, "instrument_id": iid, "scale_name": t["name"], "intended_type": "redundant",
                          "design": "total score re-using all items (excluded from the analyses that need independence)",
                          "rho": np.nan, "mu": np.nan, "sigma": np.nan, "halo": np.nan, "twin_loading_mean": np.nan,
                          "human_loading_mean": np.nan, "n_items": len(allitems)})
        reg["instruments"][iid] = node

    with open(out / "registry.yaml", "w", encoding="utf-8") as fh:
        yaml.safe_dump(reg, fh, sort_keys=False, allow_unicode=True)
    truth = pd.DataFrame(truth)
    truth.to_csv(out / "scenario_truth.csv", index=False)
    return {"data_dir": data_dir, "registry": out / "registry.yaml", "truth": truth}


def compare_with_truth(truth: pd.DataFrame, typology: pd.DataFrame) -> pd.DataFrame:
    """Intended (designed) versus recovered type per scale."""
    t = truth[truth["intended_type"] != "redundant"][["scale_id", "scale_name", "intended_type", "design"]]
    r = typology[["scale_id", "final_type", "r", "Cb", "structural_evidence_level"]].rename(
        columns={"final_type": "recovered_type", "structural_evidence_level": "structural_level"})
    m = t.merge(r, on="scale_id", how="left")
    m["match"] = m["intended_type"].astype(str) == m["recovered_type"].astype(str)
    return m


def load_scenario(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)
