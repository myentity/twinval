from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from twinval.config import load_config
from twinval.io import load_study
from twinval.registry import load_registry

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "data"
REF = Path(__file__).parent / "reference"

FAST = {"profile": {"n_boot": 200, "n_perm": 30}, "scale": {"n_boot": 200, "n_perm": 30, "cv_repeats": 10, "n_nested_boot": 60,
                                                              "n_boot_disattenuation": 200},
        "item": {"n_boot": 100, "n_boot_variance": 100}, "synthesis": {"n_boot_cluster": 300}}


@pytest.fixture(scope="session")
def registry():
    return load_registry()


@pytest.fixture(scope="session")
def cfg():
    return load_config(overrides=FAST)


@pytest.fixture(scope="session")
def study():
    return load_study(EXAMPLE)


def make_synthetic_dir(path: Path, n=150, seed=0, instruments=(148, 161), noise=1.0):
    """Write paired_wide files for a few real registry scales with a known twin process (human + noise)."""
    reg = load_registry()
    rng = np.random.default_rng(seed)
    ids = [f"S{i:04d}" for i in range(n)]
    path.mkdir(parents=True, exist_ok=True)
    trait = {iid: rng.normal(size=n) for iid in instruments}
    for iid in instruments:
        ins = reg.instruments[iid]
        for sc in ins.scales:
            k = sc.n_items
            lat = trait[iid] + rng.normal(scale=0.3, size=n)
            def resp(l, sd):
                x = l[:, None] + rng.normal(scale=sd, size=(n, k))
                q = np.quantile(x, np.linspace(0, 1, ins.n_categories + 1)[1:-1])
                return np.digitize(x, q) + ins.response_min
            h, t = resp(lat, 1.0), resp(lat + rng.normal(scale=noise, size=n), 1.0)
            df = pd.DataFrame({"ID": ids})
            for i in range(k):
                df[f"it{i + 1}_h"] = h[:, i]
            for i in range(k):
                df[f"it{i + 1}_t"] = t[:, i]
            df.to_csv(path / f"{sc.id}_paired_wide.csv", index=False)
    return path


@pytest.fixture(scope="session")
def synthetic_dir(tmp_path_factory):
    return make_synthetic_dir(tmp_path_factory.mktemp("synthetic"))
