import json
from pathlib import Path

import pandas as pd

from twinval.cli import main
from twinval.config import load_config
from twinval.pipeline import Pipeline
from conftest import EXAMPLE, FAST


def test_pipeline_without_r(tmp_path):
    cfg = load_config(overrides={**FAST, "data": {"data_dir": str(EXAMPLE)}, "structural": {"enabled": False}})
    p = Pipeline(cfg, tmp_path / "out")
    p.run(["audit", "profile", "scale", "item", "synthesis", "report"])
    out = tmp_path / "out"
    for f in ["00_audit/data_audit.csv", "01_profile/identification_summary.csv", "01_profile/profile_calibration_summary.csv",
              "02_scale/scale_level_master.csv", "04_item/item_level_master.csv", "05_synthesis/cross_level_master.csv",
              "05_synthesis/cross_level_results_matrix.csv", "report.md", "run_manifest.json", "config_used.yaml"]:
        assert (out / f).exists(), f
    assert not (out / "05_synthesis" / "scale_typology.csv").exists()          # typology needs the structural level
    man = json.loads((out / "run_manifest.json").read_text())
    assert man["n_respondents"] == 484 and any("instruments" in n for n in man["notes"])
    rm = pd.read_csv(out / "05_synthesis" / "cross_level_results_matrix.csv")
    assert {"profile", "scale", "item"} <= set(rm["level"])
    assert "typology skipped" in " ".join(man["notes"])

    # data versions are stated in the manifest, the report and every item table
    assert man["data_versions"] == {"profile_and_scale_level": "raw (input files)", "structural_models": "harmonized",
                                    "item_level": "raw"}
    assert "## Data versions" in (out / "report.md").read_text()
    assert set(pd.read_csv(out / "04_item" / "item_agreement_error.csv")["data_version"]) == {"raw"}

    # harmonised paired data are written in the input format; items with equal category sets are unchanged
    raw_dir, hdir = Path(EXAMPLE), out / "00_audit" / "harmonized_data"
    files = sorted(f.name for f in hdir.glob("*_paired_wide.csv"))
    assert files == sorted(f.name for f in raw_dir.glob("*_paired_wide.csv"))
    changed = 0
    for f in files:
        h = pd.read_csv(hdir / f)
        r = pd.read_csv(raw_dir / f).sort_values("ID").reset_index(drop=True)
        assert list(h.columns) == list(r.columns) and (h["ID"].values == r["ID"].values).all()
        changed += bool((h.drop(columns="ID").values != r.drop(columns="ID").values).any())
    log = pd.read_csv(out / "00_audit" / "harmonization_log.csv")
    assert changed == log["scale_id"].nunique() > 0

    # every figure lives in the figures/ subfolder of ITS OWN stage, never in one flat top-level folder
    assert not (out / "figures").exists()
    for f in ["01_profile/figures/profile_ccc_own_vs_cross.png", "01_profile/figures/profile_ccc_own_vs_cross.pdf",
              "02_scale/figures/scale_r_vs_cb.png", "04_item/figures/item_fidelity_map.png",
              "04_item/figures/item_variance_decomposition.png"]:
        assert (out / f).exists(), f
    assert not (out / "03_structural" / "figures").exists()                    # structural stage did not run


def test_item_level_harmonized_data_version(tmp_path):
    cfg = load_config(overrides={**FAST, "data": {"data_dir": str(EXAMPLE), "instruments": [161]}, "structural": {"enabled": False},
                                 "item": {**FAST["item"], "data_version": "harmonized"}})
    p = Pipeline(cfg, tmp_path / "out")
    p.run(["audit", "item"])
    assert set(pd.read_csv(tmp_path / "out" / "04_item" / "item_agreement_error.csv")["data_version"]) == {"harmonized"}


def test_cli_validate_and_scale_key(capsys):
    assert main(["validate", "--data", str(EXAMPLE)]) == 0
    assert "OK: 484 complete respondents" in capsys.readouterr().out
    assert main(["scale-key", "--format", "csv"]) == 0
    assert capsys.readouterr().out.count("\n") == 25
    assert main(["validate", "--data", "/nonexistent"]) == 2
