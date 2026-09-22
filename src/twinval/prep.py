"""Data preparation shared by the item-level and structural stages.

* :func:`harmonize` - symmetric merging of adjacent response categories when the observed category sets of
  Human and Twin differ for an item. Threshold-based models (ordinal CFA / ESEM) need identical category sets on
  both sides; one shared recoding is applied to Human and Twin so that ordinality and comparability are kept.
* :func:`item_battery` - the per-instrument item matrix (canonical item codes, ``_h`` / ``_t`` columns) and the
  target loading matrix used by ESEM.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .io import Study


def _levels(x: np.ndarray) -> list[int]:
    return sorted(int(v) for v in np.unique(x[~np.isnan(x)]))


def _merge_bins(lh: set, lt: set) -> list[list[int]]:
    """Greedy adjacent merge until every bin is observed on both sides."""
    bins = [[v] for v in sorted(lh | lt)]
    changed = True
    while changed:
        changed = False
        for j in range(len(bins)):
            ph, pt = any(v in lh for v in bins[j]), any(v in lt for v in bins[j])
            if ph != pt:
                if j + 1 < len(bins):
                    bins[j:j + 2] = [bins[j] + bins[j + 1]]
                else:
                    bins[j - 1:j + 1] = [bins[j - 1] + bins[j]]
                changed = True
                break
    return bins


def harmonize(study: Study) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], pd.DataFrame, pd.DataFrame]:
    """Return harmonised item arrays (same layout as ``study.items_h/items_t``) and two logs.

    Scales with an item that collapses to a single category cannot be harmonised; they are listed in the
    ``unusable`` log and returned unchanged (structural stages skip them).
    """
    out_h, out_t, log_rows, bad_rows = {}, {}, [], []
    for sid in study.scale_ids:
        h, t = study.items_h[sid].copy(), study.items_t[sid].copy()
        for i in range(h.shape[1]):
            lh, lt = _levels(h[:, i]), _levels(t[:, i])
            if lh == lt:
                continue
            bins = _merge_bins(set(lh), set(lt))
            if len(bins) < 2:
                bad_rows.append({"scale_id": sid, "item_index": i + 1, "human_levels": "|".join(map(str, lh)),
                                 "twin_levels": "|".join(map(str, lt)), "reason": "collapses to a single category"})
                continue
            m = {v: new for new, b in enumerate(bins, start=1) for v in b}
            h[:, i] = np.vectorize(m.get)(h[:, i])
            t[:, i] = np.vectorize(m.get)(t[:, i])
            log_rows.append({"scale_id": sid, "item_index": i + 1,
                             "human_levels_before": "|".join(map(str, lh)), "twin_levels_before": "|".join(map(str, lt)),
                             "n_cat_after": len(bins),
                             "merge_map": "; ".join(f"{k}={'-'.join(map(str, b))}" for k, b in enumerate(bins, 1))})
        out_h[sid], out_t[sid] = h, t
    log = pd.DataFrame(log_rows, columns=["scale_id", "item_index", "human_levels_before", "twin_levels_before",
                                          "n_cat_after", "merge_map"])
    bad = pd.DataFrame(bad_rows, columns=["scale_id", "item_index", "human_levels", "twin_levels", "reason"])
    return out_h, out_t, log, bad


def harmonization_status(study: Study, log: pd.DataFrame, bad: pd.DataFrame) -> pd.DataFrame:
    """Per-scale status: clean | harmonized | unusable."""
    merged = set(log["scale_id"]) if len(log) else set()
    unusable = set(bad["scale_id"]) if len(bad) else set()
    rows = [{"scale_id": s, "category_status": "unusable" if s in unusable else "harmonized" if s in merged else "clean"}
            for s in study.scale_ids]
    return pd.DataFrame(rows)


def item_battery(study: Study, instrument_id: int, items_h: dict, items_t: dict):
    """Instrument-level item matrix over the NON-REDUNDANT scales.

    Returns ``(battery, target)``: ``battery`` has column ``id`` and ``<item>_h`` / ``<item>_t`` for every item;
    ``target`` is the item x scale 0/1 matrix of the theoretical simple structure (target rotation for ESEM).
    Facets nested in a domain and total scores are excluded: they re-use the very same items and would make
    the joint model unidentified.
    """
    reg = study.registry
    scales = [s for s in reg.instruments[instrument_id].scales if not s.redundant and s.id in items_h]
    cols, item_scale = {}, {}
    for sc in scales:
        for k, item in enumerate(sc.items):
            if item in item_scale:
                raise ValueError(f"item {item} belongs to two non-redundant scales of instrument {instrument_id}")
            item_scale[item] = sc.id
            cols[f"{item}_h"] = items_h[sc.id][:, k]
            cols[f"{item}_t"] = items_t[sc.id][:, k]
    battery = pd.DataFrame(cols)
    battery.insert(0, "id", study.ids)
    target = pd.DataFrame(0, index=list(item_scale), columns=[s.id for s in scales])
    for item, sid in item_scale.items():
        target.loc[item, sid] = 1
    return battery, target


def item_ordinal_qc(battery: pd.DataFrame, instrument_id: int) -> pd.DataFrame:
    """Per-item, per-side descriptive QC for ordinal models: categories, rarest category, modal share."""
    rows = []
    for col in [c for c in battery.columns if c.endswith(("_h", "_t"))]:
        v = battery[col].dropna()
        counts = v.value_counts()
        full = set(range(int(v.min()), int(v.max()) + 1))
        rows.append({"instrument_id": instrument_id, "item": col[:-2], "side": "human" if col.endswith("_h") else "twin",
                     "n_categories_observed": int(counts.shape[0]), "min_category_freq": int(counts.min()),
                     "modal_category_share": float(counts.max() / len(v)),
                     "empty_categories_within_range": sorted(full - set(v.unique())),
                     "near_zero_variance": bool(v.std(ddof=1) < 0.05)})
    return pd.DataFrame(rows)
