"""Step 2: independent Monte Carlo replications of the twinval classification.

One replication = one data set of a condition (three scales). Each finished replication appends one row per scale to
``replications/<cond_id>.csv``. Replications already present there are skipped, so an interrupted run resumes where it
stopped. Replication seed = 1,000,000 + 1000 * condition index + replication number.

Usage
  python run_replications.py --reps 300 --workers 10                  # main grid
  python run_replications.py --reps 100 --workers 10 --order-check    # row-order subset (7 conditions, N = 500)
  python run_replications.py --reps 2 --conditions rho0.90_mu0.0_sg1.0_N500 ...   # selected conditions
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")          # one thread per worker (numpy / R BLAS)
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from conditions import conditions  # noqa: E402
import mc_core as M  # noqa: E402
from twinval import scale as scale_mod  # noqa: E402
from twinval.config import load_config  # noqa: E402

OUT = HERE / "replications"
OUT_ORDER = HERE / "replications_order_check"
WORK = HERE / "work_replications"
ORDER_CONDITIONS = ["rho0.90_mu0.0_sg1.0_N500", "rho0.90_mu0.0_sg1.0_dl_N500", "rho0.90_mu0.0_sg1.0_sh_N500",
                    "rho0.90_mu0.0_sg1.0_rare_N500", "rho0.90_mu0.5_sg1.0_rare_N500", "rho0.90_mu0.0_sg1.0_dup_N500",
                    "rho0.90_mu0.0_sg1.0_k2_N500"]
N_BOOT = 2000
CFG_SCALE = load_config(overrides={"scale": {"n_boot": N_BOOT, "n_perm": 1}})   # the permutation p does not enter the type


COLS = ["cond_id", "rep", "seed", "scale", "N", "r", "r_ci_lo", "r_ci_hi", "auc", "auc_ci_lo", "auc_ci_hi", "ccc", "ccc_ci_lo",
        "ccc_ci_hi", "cb", "cb_ci_lo", "cb_ci_hi", "I_rec", "C_rec", "baseline_status", "baseline_fit", "W1", "W2", "W3", "M1",
        "cfa_conclusion", "cfa_no_conclusion_reason", "cfa_flags", "w0_df", "w0_srmr", "w0_max_std_loading", "esem_applicable",
        "human_esem_admissible", "twin_esem_admissible", "paired_esem_admissible", "phi", "S_rec", "structural_flags",
        "type_rec", "n_items_harmonised", "min_categories_after", "scale_unusable", "order_candidate", "order_sensitive",
        "order_baseline_statuses", "seconds", "error"]


def seed_of(ci: int, rep: int) -> int:
    return 1_000_000 + 1000 * ci + rep


def one_replication(ci: int, cond: dict, rep: int, order_perm: int) -> list[dict]:
    seed = seed_of(ci, rep)
    wd = WORK / f"{cond['cond_id']}_r{rep:04d}{'_oc' if order_perm else ''}"
    M.cleanup(wd)
    t0 = time.time()
    try:
        st = M.generate(cond, seed, wd / "sim")
        core = scale_mod.core_table(st, CFG_SCALE, np.random.default_rng(seed)).set_index("scale_id")
        res = M.structural(st, wd / "r", order_perm=order_perm, seed=seed)
        tab, s = res["table"].set_index("scale_id"), res["status"].set_index("scale_id")
        inv = res["invariance"].set_index("scale_id")
        hlog, hbad = res["harmonization_log"], res["harmonization_unusable"]
        order = res["order"].set_index("scale_id") if res["order"] is not None and len(res["order"]) else None
        rows = []
        for sid in st.scale_ids:
            c = core.loc[sid]
            i_ax, c_ax = M.axis_I(c.pearson_r, c.auc_ci_lo), M.axis_C(c.cb)
            S_rec = tab.loc[sid, "structural_evidence_level"]
            hl = hlog[hlog.scale_id == sid] if len(hlog) else hlog
            row = {"cond_id": cond["cond_id"], "rep": rep, "seed": seed, "scale": sid.split("_")[1], "N": cond["N"],
                   "r": c.pearson_r, "r_ci_lo": c.pearson_r_ci_lo, "r_ci_hi": c.pearson_r_ci_hi,
                   "auc": c.auc, "auc_ci_lo": c.auc_ci_lo, "auc_ci_hi": c.auc_ci_hi,
                   "ccc": c.ccc, "ccc_ci_lo": c.ccc_ci_lo, "ccc_ci_hi": c.ccc_ci_hi,
                   "cb": c.cb, "cb_ci_lo": c.cb_ci_lo, "cb_ci_hi": c.cb_ci_hi,
                   "I_rec": i_ax, "C_rec": c_ax,
                   "baseline_status": s.loc[sid, "baseline_status"], "baseline_fit": s.loc[sid, "baseline_fit"],
                   "W1": s.loc[sid, "threshold_decision"], "W2": s.loc[sid, "metric_decision"], "W3": s.loc[sid, "scalar_decision"],
                   "M1": s.loc[sid, "sensitivity_millsap_metric_decision"],
                   "cfa_conclusion": s.loc[sid, "cfa_conclusion"], "cfa_no_conclusion_reason": s.loc[sid, "cfa_no_conclusion_reason"],
                   "cfa_flags": s.loc[sid, "cfa_flags"], "w0_df": inv.loc[sid, "wu_configural_df"],
                   "w0_srmr": inv.loc[sid, "wu_configural_srmr"], "w0_max_std_loading": inv.loc[sid, "wu_configural_max_std_loading"],
                   "esem_applicable": bool(tab.loc[sid, "esem_applicable"]),
                   "human_esem_admissible": tab.loc[sid, "human_esem_admissible"],
                   "twin_esem_admissible": tab.loc[sid, "twin_esem_admissible"],
                   "paired_esem_admissible": bool(tab.loc[sid, "paired_esem_admissible"]),
                   "phi": tab.loc[sid, "tucker_phi_exploratory"], "S_rec": S_rec,
                   "structural_flags": tab.loc[sid, "structural_flags"],
                   "type_rec": M.final_type(i_ax, c_ax, S_rec),
                   "n_items_harmonised": int(len(hl)), "min_categories_after": int(hl.n_cat_after.min()) if len(hl) else 5,
                   "scale_unusable": bool(len(hbad) and sid in set(hbad.scale_id)),
                   "seconds": round(time.time() - t0, 1), "error": ""}
            if order_perm > 0:                               # no candidate at all -> nothing refitted, nothing flagged
                order = order if order is not None else pd.DataFrame(columns=["order_sensitive", "baseline_statuses"])
                row.update(order_candidate=sid in order.index,
                           order_sensitive=bool(order.loc[sid, "order_sensitive"]) if sid in order.index else False,
                           order_baseline_statuses=order.loc[sid, "baseline_statuses"] if sid in order.index else "")
            rows.append(row)
        return rows
    except Exception as e:                                  # a failed replication is recorded, never silently dropped
        return [{"cond_id": cond["cond_id"], "rep": rep, "seed": seed, "scale": "ALL", "N": cond["N"],
                 "error": f"{type(e).__name__}: {e} | {traceback.format_exc(limit=2)[-300:]}",
                 "seconds": round(time.time() - t0, 1)}]
    finally:
        M.cleanup(wd)


def done_reps(path: Path) -> set[int]:
    if not path.exists():
        return set()
    d = pd.read_csv(path, usecols=["rep", "error"])
    return set(d.loc[d.error.isna() | (d.error == ""), "rep"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--conditions", nargs="*", help="cond_id values (default: all, or the order-check subset)")
    ap.add_argument("--order-check", action="store_true", help="row-order check with 3 permutations")
    a = ap.parse_args()
    out = OUT_ORDER if a.order_check else OUT
    out.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)
    conds = conditions()
    wanted = set(a.conditions or (ORDER_CONDITIONS if a.order_check else [c["cond_id"] for c in conds]))
    tasks = []
    for ci, c in enumerate(conds):
        if c["cond_id"] in wanted:
            have = done_reps(out / f"{c['cond_id']}.csv")
            tasks += [(ci, c, r) for r in range(1, a.reps + 1) if r not in have]
    print(f"{len(tasks)} replications to run ({len(wanted)} conditions, {a.workers} workers)", flush=True)
    t0, n = time.time(), 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(one_replication, ci, c, r, 3 if a.order_check else 0): c["cond_id"] for ci, c, r in tasks}
        for f in as_completed(futs):
            rows = pd.DataFrame(f.result()).reindex(columns=COLS)        # one fixed column layout per file
            path = out / f"{futs[f]}.csv"
            rows.to_csv(path, mode="a", header=not path.exists(), index=False)
            n += 1
            if n % 50 == 0 or n == len(tasks):
                el = time.time() - t0
                print(f"{n}/{len(tasks)} done, {el / 60:.1f} min, ~{el / n * (len(tasks) - n) / 60:.0f} min left", flush=True)


if __name__ == "__main__":
    main()
