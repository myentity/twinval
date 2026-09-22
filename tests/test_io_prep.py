import numpy as np
import pandas as pd
import pytest

from twinval import prep
from twinval.io import DataError, load_study
from twinval.registry import load_registry


def test_registry_counts(registry):
    assert set(registry.instruments) == {148, 161} and len(registry.scales) == 24
    assert len(registry.nonredundant()) == 8                      # 5 BFI-2 domains + 3 SPH sub-scales
    assert registry.esem_applicable(148) and registry.esem_applicable(161)
    assert registry.scales["148_Trust"].parent == "148_Agreeableness" and registry.scales["148_Trust"].redundant


def test_example_loads(study):
    assert study.n == 484 and set(study.instrument_ids) == {148, 161}
    assert len(study.scale_ids) == 24 and len(study.nonredundant_ids) == 8
    # BFI-2 domains are means, SPH scores are sums
    assert study.H["148_Extraversion"].between(1, 5).all() and study.H["161_OverallScore"].max() > 12


def test_ids_aligned_across_scales(synthetic_dir):
    st = load_study(synthetic_dir)
    a = pd.read_csv(synthetic_dir / "148_Extraversion_paired_wide.csv", dtype={"ID": str})
    b = pd.read_csv(synthetic_dir / "161_EmotionalWellBeing_paired_wide.csv", dtype={"ID": str})
    assert list(a["ID"]) == list(b["ID"]) == list(st.ids)


def test_missing_column_and_range_errors(synthetic_dir, tmp_path):
    import shutil
    d = tmp_path / "bad"
    shutil.copytree(synthetic_dir, d)
    f = d / "161_EmotionalWellBeing_paired_wide.csv"
    df = pd.read_csv(f)
    df.drop(columns=["it3_t"]).to_csv(f, index=False)
    with pytest.raises(DataError, match="missing columns"):
        load_study(d)
    df.assign(it1_h=99).to_csv(f, index=False)
    with pytest.raises(DataError, match="outside the registry range"):
        load_study(d)


def test_incomplete_respondents_are_dropped(synthetic_dir, tmp_path):
    import shutil
    d = tmp_path / "gaps"
    shutil.copytree(synthetic_dir, d)
    f = d / "148_Anxiety_paired_wide.csv"
    df = pd.read_csv(f)
    df.loc[:4, "it1_h"] = np.nan
    df.to_csv(f, index=False)
    st = load_study(d)
    assert st.n == 145


def test_harmonize_merges_categories_symmetrically(synthetic_dir, tmp_path):
    import shutil
    d = tmp_path / "harm"
    shutil.copytree(synthetic_dir, d)
    f = d / "161_EmotionalWellBeing_paired_wide.csv"
    df = pd.read_csv(f)
    df.loc[df["it1_h"] == 1, "it1_h"] = 2            # human never uses category 1 -> categories 1|2 merged
    df.loc[df["it1_t"] == 6, "it1_t"] = 5            # twin never uses 6 -> 5|6 merged
    df.to_csv(f, index=False)
    st = load_study(d)
    hh, ht, log, bad = prep.harmonize(st)
    assert len(bad) == 0 and (log["scale_id"] == "161_EmotionalWellBeing").sum() == 1
    assert sorted(np.unique(hh["161_EmotionalWellBeing"][:, 0])) == sorted(np.unique(ht["161_EmotionalWellBeing"][:, 0]))
    # untouched items keep their values
    assert np.array_equal(hh["161_EmotionalWellBeing"][:, 1], st.items_h["161_EmotionalWellBeing"][:, 1])


def test_item_battery_excludes_redundant_scales(study):
    hh, ht, *_ = prep.harmonize(study)
    battery, target = prep.item_battery(study, 148, hh, ht)
    assert target.shape == (61, 5) and battery.shape[1] == 1 + 2 * 61      # 5 domains, 61 items
    assert (target.sum(1) == 1).all()                 # simple structure: each item loads on exactly one target factor
