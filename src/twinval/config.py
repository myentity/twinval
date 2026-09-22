"""Configuration handling: packaged defaults + user YAML overrides (deep merge)."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import yaml

DEFAULT_CONFIG_PATH = Path(__file__).parent / "configs" / "default.yaml"


class AttrDict(dict):
    """dict with attribute access; nested dicts are wrapped lazily."""

    def __getattr__(self, name):
        try:
            v = self[name]
        except KeyError as e:  # pragma: no cover
            raise AttributeError(name) from e
        return AttrDict(v) if isinstance(v, dict) and not isinstance(v, AttrDict) else v

    def __setattr__(self, name, value):
        self[name] = value


def _deep_update(base: dict, upd: dict) -> dict:
    for k, v in upd.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _deep_update(base[k], v)
        else:
            base[k] = copy.deepcopy(v)
    return base


def load_config(path: str | Path | None = None, overrides: dict | None = None) -> AttrDict:
    """Load the default configuration, then apply a user YAML file and explicit overrides."""
    with open(DEFAULT_CONFIG_PATH, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if path is not None:
        with open(path, encoding="utf-8") as fh:
            user = yaml.safe_load(fh) or {}
        _deep_update(cfg, user)
    if overrides:
        _deep_update(cfg, overrides)
    return _wrap(cfg)


def _wrap(d):
    if isinstance(d, dict):
        return AttrDict({k: _wrap(v) for k, v in d.items()})
    return d


def config_hash(cfg: dict) -> str:
    blob = json.dumps(cfg, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:12]
