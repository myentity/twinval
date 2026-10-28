"""Shared pieces of the Monte Carlo study: data generation and the twinval decisions for one data set.

The structural part is the same code path as ``twinval run`` (stage ``structural``): category harmonisation, paired ordinal
CFA (r/paired_cfa_invariance.R), ESEM + Procrustes/Tucker point estimates, admissibility gate, structural status and
``scale_structural_table``. The Tucker bootstrap is not run: it does not enter the final type except through the
``bootstrap_unstable`` flag.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from twinval import structural as S
from twinval import synthesis as synth
from twinval.config import load_config
from twinval.io import load_study
from twinval.registry import load_registry
from twinval.synth import simulate

from conditions import spec

CFG = load_config(overrides={"structural": {"tucker_bootstrap": False, "order_check_permutations": 0, "resume": False}})


def generate(cond: dict, seed: int, out: Path, n: int | None = None):
    simulate(spec(cond, seed, n), out)
    return load_study(out / "data", load_registry(out / "registry.yaml"))


def structural(study, workdir: Path, cfg=CFG, order_perm: int = 0, seed: int = 1) -> dict:
    """Paired CFA + ESEM + gate -> structural table (one row per scale)."""
    workdir.mkdir(parents=True, exist_ok=True)
    logs = workdir / "logs"
    logs.mkdir(exist_ok=True)
    inp = S.prepare_inputs(study, cfg, workdir)
    S.run_r("paired_cfa_invariance.R", [workdir / "cfa_manifest.csv", workdir / "cfa_results", cfg.thresholds.cfa_delta_cfi_min,
                                        cfg.thresholds.cfa_delta_rmsea_max], logs / "cfa.log", cfg.structural.rscript)
    inv = pd.read_csv(workdir / "cfa_results" / "invariance_results.csv")
    status = S.cfa_scale_status(inv, cfg)
    order = None
    if order_perm > 0:
        order = S.cfa_order_sensitivity(inv, inp["manifest"], workdir, order_perm, seed, cfg, logs)
    esem_ids = [e["instrument_id"] for e in inp["esem_instruments"]]
    tucker = gate = None
    if esem_ids:
        S.run_r("esem_fit.R", [workdir], logs / "esem.log", cfg.structural.rscript)
        diag = pd.read_csv(workdir / "esem_model_diagnostics.csv")
        tucker = S.tucker_point_estimates(workdir, diag, inp["esem_instruments"])["tucker_by_factor"]
        gate = S.admissibility(diag, esem_ids)
    table = S.scale_structural_table(study, cfg, status, tucker, gate, esem_ids, order)
    return {"invariance": inv, "status": status, "table": table, "tucker": tucker, "gate": gate, "order": order,
            "harmonization_log": inp["harmonization_log"], "harmonization_unusable": inp["harmonization_unusable"]}


def axis_I(r: float, auc_ci_lo: float, t=CFG.thresholds) -> str:
    return "adequate" if (r >= t.r_practical and auc_ci_lo > t.auc_ci_low) else "inadequate"


def axis_C(cb: float, t=CFG.thresholds) -> str:
    return "adequate" if cb >= t.cb_adequate else "distorted"


def final_type(i: str, c: str, s: str) -> str:
    """The type decision of twinval (synthesis._final_type) for given axes."""
    return synth._final_type({"axis_I_practical": i, "axis_C": c, "axis_S": s})


def scale_scores(study) -> tuple[np.ndarray, np.ndarray]:
    H, T = study.matrices()
    return H, T


def cleanup(path: Path):
    shutil.rmtree(path, ignore_errors=True)
