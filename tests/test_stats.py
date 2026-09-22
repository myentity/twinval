import numpy as np
import pytest

from twinval import stats


def test_ccc_equals_r_times_cb():
    rng = np.random.default_rng(1)
    h = rng.normal(size=(300, 4))
    t = 0.5 * h + rng.normal(size=(300, 4)) + 0.7
    c = stats.ccc_components(h, t)
    assert np.allclose(c["ccc"], c["r"] * c["cb"])
    for j in range(4):
        assert np.isclose(c["ccc"][j], stats.ccc_1d(h[:, j], t[:, j]))


def test_cb_is_one_for_identical_distribution_shift_free():
    x = np.random.default_rng(2).normal(size=200)
    c = stats.ccc_components(x, x.copy())
    assert np.isclose(c["cb"][0], 1) and np.isclose(c["ccc"][0], 1)


def test_human_z_uses_human_reference_for_both_sides():
    rng = np.random.default_rng(3)
    H = rng.normal(5, 2, size=(100, 2))
    T = H + 3
    zH, zT = stats.human_z(H, T)
    assert np.allclose(zH.mean(0), 0) and np.allclose(zH.std(0, ddof=1), 1)
    assert np.allclose(zT.mean(0), 3 / H.std(0, ddof=1))     # twin bias stays visible


def test_bootstrap_p_floor_and_symmetry():
    p, floor = stats.bootstrap_p(np.linspace(0.2, 0.9, 5000))
    assert floor and np.isclose(p, 2 / 5001)                   # never exactly 0
    p2, floor2 = stats.bootstrap_p(np.random.default_rng(0).normal(0, 1, 5000))
    assert 0.5 < p2 <= 1 and not floor2
    assert np.isnan(stats.bootstrap_p(np.arange(10.0))[0])


def test_grouped_folds_keep_duplicates_together():
    rng = np.random.default_rng(4)
    idx = rng.integers(0, 60, size=60)
    fid = stats.grouped_fold_ids(idx, 5, rng)
    for orig in np.unique(idx):
        assert len(set(fid[idx == orig])) == 1
    assert set(fid) == set(range(5))


def test_own_cross_stats_perfect_identification():
    n = 40
    M = 1 + np.ones((n, n)) - np.eye(n)               # own distance 0, others 1
    s = stats.own_cross_stats(M, lower_is_better=True)
    assert np.all(s["rank"] == 1) and np.isclose(s["auc"].mean(), 1)


def test_bh_adjust_handles_nan():
    q = stats.bh_adjust(np.array([0.01, 0.04, np.nan, 0.03]))
    assert np.isnan(q[2]) and (q[[0, 1, 3]] >= [0.01, 0.04, 0.03]).all()


def test_cluster_bootstrap_needs_enough_clusters():
    import pandas as pd
    rng = np.random.default_rng(5)
    df = pd.DataFrame({"x": rng.normal(size=40), "cl": np.repeat([1, 2], 20)})
    df["y"] = df["x"] + rng.normal(size=40)
    r = stats.cluster_bootstrap_spearman(df, "x", "y", "cl", 300, rng, min_clusters=4)
    assert np.isnan(r["ci_lo"]) and r["n_clusters"] == 2 and not np.isnan(r["rho"])
