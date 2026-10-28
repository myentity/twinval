"""Level 1, 2 and 4 on synthetic data with a known twin process (no R needed)."""
import numpy as np
import pytest

from twinval import item, profile, scale
from twinval.io import load_study
from conftest import make_synthetic_dir


@pytest.fixture(scope="module")
def good_and_bad(tmp_path_factory, cfg):
    good = load_study(make_synthetic_dir(tmp_path_factory.mktemp("good"), n=200, seed=1, noise=0.3))
    bad = load_study(make_synthetic_dir(tmp_path_factory.mktemp("bad"), n=200, seed=1, noise=3.0))
    return good, bad


def test_scale_metrics_order_with_twin_quality(good_and_bad, cfg):
    good, bad = good_and_bad
    rg = scale.core_table(good, cfg, np.random.default_rng(0))
    rb = scale.core_table(bad, cfg, np.random.default_rng(0))
    assert rg["pearson_r"].mean() > rb["pearson_r"].mean() + 0.2
    assert rg["auc"].mean() > rb["auc"].mean()
    assert np.allclose(rg["ccc"], rg["pearson_r"] * rg["cb"])
    assert (rg["pearson_r_ci_lo"] < rg["pearson_r"]).all() and (rg["pearson_r"] < rg["pearson_r_ci_hi"]).all()


def test_recalibration_improves_biased_twin(good_and_bad, cfg):
    _, bad = good_and_bad
    # inject a constant bias into the twin scores: LS recalibration must remove it
    bad.T["148_Extraversion"] = bad.T["148_Extraversion"] + 1.5
    _, per_scale = scale.recalibration(bad, cfg, np.random.default_rng(0))
    row = per_scale.set_index("scale_id").loc["148_Extraversion"]
    assert row["delta_ccc_ls"] > 0.05 and row["cb_ls"] > row["cb_raw"]


def test_disattenuation_exceeds_raw_r(good_and_bad, cfg):
    good, _ = good_and_bad
    core = scale.core_table(good, cfg, np.random.default_rng(0))
    dis = scale.disattenuation(good, cfg, np.random.default_rng(0), core)
    ok = dis.dropna(subset=["r_dis"])
    assert (ok["r_dis"] >= ok["pearson_r"] - 1e-9).all()
    assert set(dis["reliability_type"].dropna()) <= {"alpha", "spearman_brown"}


def test_profile_identification_and_calibration(good_and_bad, cfg):
    good, bad = good_and_bad
    ig = profile.identification(good, cfg, np.random.default_rng(0))["identification_summary"]
    ib = profile.identification(bad, cfg, np.random.default_rng(0))["identification_summary"]
    assert ig["auc"].min() > ib["auc"].max() and ig["auc"].min() > 0.75
    cal = profile.calibration(good, cfg, np.random.default_rng(0))["profile_calibration_summary"].set_index("metric")
    assert cal.loc["own_ccc_full", "value"] > cal.loc["cross_ccc_full_mean", "value"] + 0.2
    assert abs(cal.loc["cross_ccc_full_mean", "value"]) < 0.1


def test_equal_weights_sum_to_one(good_and_bad):
    good, _ = good_and_bad
    ids = good.nonredundant_ids
    w = profile._weights(ids, good.groups(ids))
    assert np.isclose(w.sum(), 1.0)
    for cols in good.groups(ids).values():
        assert np.isclose(sum(w[ids.index(c)] for c in cols), 1 / len(good.groups(ids)))


def test_item_level_uses_theoretical_range(good_and_bad, cfg):
    good, _ = good_and_bad
    panel = item.build_panel(good)
    ag, bias, audit = item.agreement_and_bias(panel, cfg, np.random.default_rng(0))
    # normalised MAE is computed against the registry range, not the observed one
    j = 0
    rng_ = panel.meta.loc[j, "item_max"] - panel.meta.loc[j, "item_min"]
    assert np.isclose(ag.loc[j, "mae_range"], ag.loc[j, "mae_raw"] / rng_)
    assert ag["exact_agreement"].between(0, 1).all() and bias["sd_ratio"].notna().all()
    # extreme categories of the allowable range, both sides
    lo, hi = panel.meta.loc[j, "item_min"], panel.meta.loc[j, "item_max"]
    assert np.isclose(bias.loc[j, "extreme_share_human"], np.mean((panel.H[:, j] == lo) | (panel.H[:, j] == hi)))
    assert np.allclose(bias.extreme_share_difference, bias.extreme_share_twin - bias.extreme_share_human)


def test_variance_decomposition_sums_to_100(good_and_bad, cfg):
    good, _ = good_and_bad
    dec, con = item.variance_decomposition(item.build_panel(good), cfg, np.random.default_rng(0))
    tot = dec[["pct_person", "pct_item", "pct_person_item_residual"]].sum(axis=1)
    assert np.allclose(tot, 100)
    assert (con["instrument_id"].astype(str) == "ALL_INSTRUMENT_WEIGHTED").any()
    # absolute components add up and the person effects of a good twin correspond
    assert np.allclose(dec[["var_person", "var_item", "var_person_item_residual"]].sum(axis=1), dec["var_total"])
    per = con[con["instrument_id"].astype(str) != "ALL_INSTRUMENT_WEIGHTED"]["person_effect_r"]
    assert per.between(-1, 1).all() and per.mean() > 0.3
