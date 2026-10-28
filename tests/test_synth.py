"""Synthetic data generator (`twinval simulate`) and the branches that need many instruments."""
import numpy as np
import pandas as pd
import pytest

from twinval.config import load_config
from twinval.io import load_study
from twinval.pipeline import Pipeline
from twinval.registry import load_registry
from twinval.synth import compare_with_truth, load_scenario, simulate
from conftest import FAST, ROOT

SHIPPED = ROOT / "examples" / "synthetic"


def spec(**twin_a):
    return {"seed": 3, "n_respondents": 400, "instruments": [
        {"id": 1, "key": "S1", "name": "One", "response": [1, 5], "scales": [
            {"name": "A", "n_items": 5, "twin": twin_a}, {"name": "B", "n_items": 4, "twin": {"rho": 0.85}}]},
        {"id": 2, "key": "S2", "name": "Two", "response": [1, 7], "total": {"name": "Total"}, "scales": [
            {"name": "C", "n_items": 4, "reverse": [2], "twin": {"rho": 0.5}}, {"name": "D", "n_items": 3}]}]}


def load(out):
    reg = load_registry(out / "registry.yaml")
    return reg, load_study(out / "data", reg, None, "ID", True, "drop")


def score(out, sid):
    d = pd.read_csv(out / "data" / f"{sid}_paired_wide.csv")
    return d.filter(regex="_h$").mean(axis=1).to_numpy(), d.filter(regex="_t$").mean(axis=1).to_numpy(), d


def test_simulate_writes_valid_framework_input(tmp_path):
    res = simulate(spec(rho=0.85), tmp_path)
    reg, st = load(tmp_path)
    assert st.n == 400 and len(st.scale_ids) == 5 and st.instrument_ids == [1, 2]
    assert reg.scales["2_Total"].redundant and reg.scales["2_Total"].n_items == 7       # total score re-uses all items
    assert reg.scales["2_C"].reverse_items == ("S2_2",)
    assert len(res["truth"]) == 5
    for sid in st.scale_ids:
        d = pd.read_csv(tmp_path / "data" / f"{sid}_paired_wide.csv")
        lo, hi = (1, 5) if sid.startswith("1_") else (1, 7)
        assert d.drop(columns="ID").to_numpy().min() >= lo and d.drop(columns="ID").to_numpy().max() <= hi


def test_simulate_is_deterministic_and_instruments_are_independent(tmp_path):
    simulate(spec(rho=0.85), tmp_path / "a")
    simulate(spec(rho=0.85), tmp_path / "b")
    s2 = spec(rho=0.85)
    s2["instruments"][1]["scales"][0]["twin"] = {"rho": 0.1, "mu": 1.0}                # change ONLY instrument 2
    simulate(s2, tmp_path / "c")
    for f in (tmp_path / "a" / "data").glob("*.csv"):
        assert f.read_text() == (tmp_path / "b" / "data" / f.name).read_text()
    for f in (tmp_path / "a" / "data").glob("1_*.csv"):
        assert f.read_text() == (tmp_path / "c" / "data" / f.name).read_text()
    assert (tmp_path / "a" / "data" / "2_C_paired_wide.csv").read_text() != (tmp_path / "c" / "data" / "2_C_paired_wide.csv").read_text()


def test_twin_parameters_do_what_they_say(tmp_path):
    simulate(spec(rho=0.9), tmp_path / "hi")
    simulate(spec(rho=0.2), tmp_path / "lo")
    h, t, _ = score(tmp_path / "hi", "1_A")
    h2, t2, _ = score(tmp_path / "lo", "1_A")
    assert np.corrcoef(h, t)[0, 1] > 0.6 > 0.3 > np.corrcoef(h2, t2)[0, 1]              # rho -> individuation

    simulate(spec(rho=0.9, mu=0.8, sigma=0.5), tmp_path / "cal")
    h, t, _ = score(tmp_path / "cal", "1_A")
    assert t.mean() - h.mean() > 0.3 and t.std() < 0.8 * h.std()                          # mu, sigma -> bias, compression

    simulate(spec(rho=0.9, never_lowest=True, never_highest=True), tmp_path / "ext")
    _, _, d = score(tmp_path / "ext", "1_A")
    assert d.filter(regex="_t$").to_numpy().min() == 2 and d.filter(regex="_t$").to_numpy().max() == 4
    assert d.filter(regex="_h$").to_numpy().min() == 1

    simulate(spec(rho=0.9, duplicate=[[1, 2], [1, 3]], duplicate_corr=1.0), tmp_path / "dup")
    _, _, d = score(tmp_path / "dup", "1_A")
    assert d["it1_t"].corr(d["it2_t"], method="spearman") > 0.85 > d["it1_h"].corr(d["it2_h"], method="spearman")

    simulate(spec(rho=0.9, near_constant=[2], near_constant_keep=0.02), tmp_path / "nc")
    _, _, d = score(tmp_path / "nc", "1_A")
    assert d["it2_t"].value_counts(normalize=True).iloc[0] > 0.9

    simulate(spec(rho=0.9, copy_items=[1]), tmp_path / "cp")
    _, _, d = score(tmp_path / "cp", "1_A")
    assert (d["it1_t"] == d["it1_h"]).all()


def test_shipped_synthetic_example_is_consistent():
    truth = pd.read_csv(SHIPPED / "scenario_truth.csv")
    reg = load_registry(SHIPPED / "registry.yaml")
    files = {f.name.replace("_paired_wide.csv", "") for f in (SHIPPED / "data").glob("*_paired_wide.csv")}
    assert files == set(reg.scales) == set(truth["scale_id"])
    scen = load_scenario(SHIPPED / "scenario.yaml")
    n_scales = sum(len(i["scales"]) + bool(i.get("total")) for i in scen["instruments"])
    assert n_scales == len(truth) == 29
    st = load_study(SHIPPED / "data", reg, None, "ID", True, "drop")
    assert st.n == 500 and len(st.instrument_ids) == 10
    assert all(str(i).startswith("SY") for i in st.ids)                                   # synthetic ids, unrelated to the real example


def test_many_instrument_branches_without_r(tmp_path):
    """Leave-one-instrument-out and the cluster bootstrap need >= 3 / >= 4 instruments: exercised here on synthetic data."""
    sp = {"seed": 5, "n_respondents": 250, "instruments": [
        {"id": 10 + i, "key": f"M{i}", "name": f"Inst {i}", "response": [1, 5],
         "scales": [{"name": f"S{j}", "n_items": 4, "twin": {"rho": 0.4 + 0.1 * ((i + j) % 5)}} for j in range(2)]} for i in range(5)]}
    simulate(sp, tmp_path / "s")
    cfg = load_config(overrides={**FAST, "data": {"data_dir": str(tmp_path / "s" / "data")}, "structural": {"enabled": False}})
    p = Pipeline(cfg, tmp_path / "out", load_registry(tmp_path / "s" / "registry.yaml"))
    p.run(["audit", "profile", "scale", "item", "synthesis", "report"])
    loio = pd.read_csv(tmp_path / "out" / "01_profile" / "identification_leave_one_instrument_out.csv")
    assert len(loio) == 5 and loio["auc_full"].notna().all()
    assoc = pd.read_csv(tmp_path / "out" / "05_synthesis" / "cross_level_associations.csv")
    assert (assoc["n_clusters"] == 5).all() and assoc["ci_lo"].notna().any()             # cluster bootstrap is computed with 5 clusters
    man = (tmp_path / "out" / "report.md").read_text()
    assert "Only 2 instruments" not in man


@pytest.mark.r
def test_shipped_synthetic_example_shows_all_types(tmp_path):
    from twinval.structural import r_available
    ok, msg = r_available()
    if not ok:
        pytest.skip(msg)
    cfg = load_config(overrides={**FAST, "data": {"data_dir": str(SHIPPED / "data")}, "structural": {"n_boot_tucker": 3}})
    p = Pipeline(cfg, tmp_path / "out", load_registry(SHIPPED / "registry.yaml"))
    p.run()
    ty = pd.read_csv(tmp_path / "out" / "05_synthesis" / "scale_typology.csv")
    # rules_version 2: the redundant total score (906_Total) has a misfitting baseline CFA
    assert set(ty["final_type"]) == {"A", "B", "C", "D", "E", "F", "cfa_misfit", "model_not_estimable"}
    cmp = compare_with_truth(pd.read_csv(SHIPPED / "scenario_truth.csv"), ty, 2)
    assert cmp["match"].mean() >= 0.85
    assoc = pd.read_csv(tmp_path / "out" / "05_synthesis" / "cross_level_associations.csv")
    assert assoc["ci_lo"].notna().any()

    # figures: one figures/ subfolder per stage (including structural, which needs R), no flat top-level folder
    assert not (tmp_path / "out" / "figures").exists()
    for f in ["01_profile/figures/profile_ccc_own_vs_cross.png", "02_scale/figures/scale_r_vs_cb.png",
              "03_structural/figures/tucker_phi_admissibility.png", "04_item/figures/item_fidelity_map.png",
              "05_synthesis/figures/scale_typology_map.png", "05_synthesis/figures/cross_level_associations_heatmap.png"]:
        assert (tmp_path / "out" / f).exists(), f
