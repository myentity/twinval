"""Instrument / scale registry: the key that links data files, items, scoring and redundancy rules."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

DEFAULT_REGISTRY_PATH = Path(__file__).parent / "registry" / "instruments.yaml"


@dataclass(frozen=True)
class Scale:
    id: str
    instrument_id: int
    name_en: str
    name_ru: str
    legacy_id: str
    level: str                 # scale | subscale | domain | facet | total
    scoring: str               # 'mean' | 'sum'
    items: tuple[str, ...]     # canonical item codes in the order of it1..itK columns
    reverse_items: tuple[str, ...]
    parent: str | None
    redundant: bool            # True for facets nested in a domain and for total scores

    @property
    def n_items(self) -> int:
        return len(self.items)


@dataclass
class Instrument:
    id: int
    key: str
    name_en: str
    name_ru: str
    item_prefix: str
    reference: str
    response_min: int
    response_max: int
    scales: list[Scale] = field(default_factory=list)

    @property
    def n_categories(self) -> int:
        return self.response_max - self.response_min + 1


class Registry:
    def __init__(self, instruments: dict[int, Instrument]):
        self.instruments = instruments
        self.scales: dict[str, Scale] = {s.id: s for i in instruments.values() for s in i.scales}

    # ---------------------------------------------------------------- lookups
    def instrument_of(self, scale_id: str) -> Instrument:
        return self.instruments[self.scales[scale_id].instrument_id]

    def subset(self, instrument_ids=None, scale_ids=None) -> "Registry":
        """Registry restricted to some instruments and/or scales (used for partial data sets)."""
        out = {}
        for iid, ins in self.instruments.items():
            if instrument_ids is not None and iid not in instrument_ids:
                continue
            keep = [s for s in ins.scales if scale_ids is None or s.id in scale_ids]
            if keep:
                out[iid] = Instrument(ins.id, ins.key, ins.name_en, ins.name_ru, ins.item_prefix,
                                      ins.reference, ins.response_min, ins.response_max, keep)
        return Registry(out)

    def nonredundant(self, scale_ids=None) -> list[str]:
        return [s.id for s in self.scales.values() if not s.redundant and (scale_ids is None or s.id in scale_ids)]

    def instrument_groups(self, scale_ids) -> dict[int, list[str]]:
        groups: dict[int, list[str]] = {}
        for sid in scale_ids:
            groups.setdefault(self.scales[sid].instrument_id, []).append(sid)
        return groups

    def esem_applicable(self, instrument_id: int) -> bool:
        """An instrument is multidimensional (ESEM applies) if it has > 1 non-redundant scale."""
        n = sum(1 for s in self.instruments[instrument_id].scales if not s.redundant)
        return n > 1


def load_registry(path: str | Path | None = None) -> Registry:
    with open(path or DEFAULT_REGISTRY_PATH, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    instruments: dict[int, Instrument] = {}
    for iid, node in raw["instruments"].items():
        iid = int(iid)
        ins = Instrument(iid, node["key"], node["name_en"], node["name_ru"], node["item_prefix"],
                         node.get("reference", ""), int(node["response_min"]), int(node["response_max"]))
        for s in node["scales"]:
            ins.scales.append(Scale(
                id=s["id"], instrument_id=iid, name_en=s["name_en"], name_ru=s["name_ru"],
                legacy_id=s["legacy_id"], level=s["level"], scoring=s["scoring"],
                items=tuple(s["items"]), reverse_items=tuple(s.get("reverse_items", [])),
                parent=s.get("parent"), redundant=bool(s["redundant"])))
        instruments[iid] = ins
    return Registry(instruments)
