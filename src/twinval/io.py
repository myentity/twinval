"""Loading and validation of paired Human/Twin data.

Input convention (one file per scale or sub-scale)::

    <scale_id>_paired_wide.csv   columns: ID, it1_h ... itK_h, it1_t ... itK_t

* ``ID``    respondent identifier (the same across all files of a study)
* ``itk_h`` k-th item answered by the human (already reverse-scored; integer categories)
* ``itk_t`` the same item answered by the digital twin of that human
* the order ``it1..itK`` equals the order of ``items`` of the scale in the registry
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .registry import Registry, load_registry

log = logging.getLogger("twinval")


class DataError(ValueError):
    """Raised when the input files do not follow the documented format."""


@dataclass
class Study:
    """A complete-case paired data set: item responses and scale scores for N respondents."""

    registry: Registry
    ids: np.ndarray                        # (N,) respondent ids (sorted)
    items_h: dict[str, np.ndarray]         # scale_id -> (N, K) human item responses
    items_t: dict[str, np.ndarray]         # scale_id -> (N, K) twin item responses
    H: pd.DataFrame                        # (N, S) human scale scores
    T: pd.DataFrame                        # (N, S) twin scale scores
    audit: pd.DataFrame = field(default_factory=pd.DataFrame)

    @property
    def n(self) -> int:
        return len(self.ids)

    @property
    def scale_ids(self) -> list[str]:
        return list(self.H.columns)

    @property
    def nonredundant_ids(self) -> list[str]:
        return [s for s in self.scale_ids if not self.registry.scales[s].redundant]

    @property
    def instrument_ids(self) -> list[int]:
        return sorted({self.registry.scales[s].instrument_id for s in self.scale_ids})

    def groups(self, scale_ids=None) -> dict[int, list[str]]:
        return self.registry.instrument_groups(scale_ids if scale_ids is not None else self.scale_ids)

    def matrices(self, scale_ids=None) -> tuple[np.ndarray, np.ndarray]:
        cols = scale_ids if scale_ids is not None else self.scale_ids
        return self.H[cols].to_numpy(float), self.T[cols].to_numpy(float)


def _score(items: np.ndarray, scoring: str) -> np.ndarray:
    return items.mean(axis=1) if scoring == "mean" else items.sum(axis=1)


def _find_file(data_dir: Path, scale, legacy: bool) -> Path | None:
    cands = [data_dir / f"{scale.id}_paired_wide.csv"]
    if legacy:
        cands.append(data_dir / f"{scale.legacy_id}_paired_wide.csv")
    for c in cands:
        if c.exists():
            return c
    return None


def load_study(data_dir: str | Path, registry: Registry | None = None, instruments=None,
               id_column: str = "ID", legacy_file_names: bool = True, missing: str = "drop") -> Study:
    """Load all ``*_paired_wide.csv`` files found for the registry scales and build a :class:`Study`."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        raise DataError(f"data directory not found: {data_dir}")
    registry = registry or load_registry()
    reg = registry.subset(instrument_ids=set(instruments) if instruments else None)

    frames: dict[str, pd.DataFrame] = {}
    audit_rows = []
    for scale in reg.scales.values():
        path = _find_file(data_dir, scale, legacy_file_names)
        if path is None:
            audit_rows.append({"scale_id": scale.id, "status": "file_not_found"})
            continue
        df = pd.read_csv(path, dtype={id_column: str})
        need = [id_column] + [f"it{k}_h" for k in range(1, scale.n_items + 1)] + \
               [f"it{k}_t" for k in range(1, scale.n_items + 1)]
        absent = [c for c in need if c not in df.columns]
        if absent:
            raise DataError(f"{path.name}: missing columns {absent[:6]}{'...' if len(absent) > 6 else ''} "
                            f"(registry expects {scale.n_items} items)")
        extra = [c for c in df.columns if re.fullmatch(r"it\d+_[ht]", c) and c not in need]
        if extra:
            raise DataError(f"{path.name}: unexpected item columns {extra[:6]} (registry expects {scale.n_items} items)")
        if df[id_column].duplicated().any():
            raise DataError(f"{path.name}: duplicated respondent IDs")
        frames[scale.id] = df[need].set_index(id_column)
        audit_rows.append({"scale_id": scale.id, "status": "loaded", "file": path.name, "n_rows": len(df)})

    if not frames:
        raise DataError(f"no *_paired_wide.csv files matching the registry were found in {data_dir}")

    # ------------------------------------------------------------ alignment (complete case)
    ids_all = set.union(*[set(f.index) for f in frames.values()])
    ids_common = set.intersection(*[set(f.index) for f in frames.values()])
    complete = set(ids_common)
    for sid, f in frames.items():
        ok = f.notna().all(axis=1)
        complete &= set(f.index[ok])
    if missing == "error" and len(complete) != len(ids_all):
        raise DataError(f"{len(ids_all) - len(complete)} respondents are missing or incomplete in some scale file")
    ids = np.array(sorted(complete))
    if len(ids) < 30:
        raise DataError(f"only {len(ids)} complete respondents remain after alignment")
    if len(complete) < len(ids_all):
        log.warning("dropped %d of %d respondents (not present/complete in all scale files)",
                    len(ids_all) - len(complete), len(ids_all))

    # ------------------------------------------------------------ range check + scoring
    items_h, items_t, sc_h, sc_t = {}, {}, {}, {}
    for sid, f in frames.items():
        sc = reg.scales[sid]
        ins = reg.instruments[sc.instrument_id]
        f = f.loc[ids]
        k = sc.n_items
        h = f[[f"it{i}_h" for i in range(1, k + 1)]].to_numpy(float)
        t = f[[f"it{i}_t" for i in range(1, k + 1)]].to_numpy(float)
        if not (np.all(h == np.round(h)) and np.all(t == np.round(t))):
            raise DataError(f"{sid}: item responses must be integers")
        lo, hi = min(h.min(), t.min()), max(h.max(), t.max())
        if lo < ins.response_min or hi > ins.response_max:
            raise DataError(f"{sid}: responses in [{lo:g}, {hi:g}] fall outside the registry range "
                            f"[{ins.response_min}, {ins.response_max}] of {ins.name_en}")
        items_h[sid], items_t[sid] = h, t
        sc_h[sid], sc_t[sid] = _score(h, sc.scoring), _score(t, sc.scoring)

    loaded = [s for s in reg.scales if s in frames]
    H = pd.DataFrame(sc_h, index=ids)[loaded]
    T = pd.DataFrame(sc_t, index=ids)[loaded]
    reg = reg.subset(scale_ids=set(loaded))

    audit = pd.DataFrame(audit_rows)
    audit["n_complete_respondents"] = len(ids)
    n_missing_files = int((audit["status"] == "file_not_found").sum())
    if n_missing_files:
        log.warning("%d registry scales have no data file and are skipped", n_missing_files)
    return Study(reg, ids, items_h, items_t, H, T, audit)
