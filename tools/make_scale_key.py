"""Regenerate docs/SCALE_KEY.md from the packaged registry."""
from pathlib import Path

from twinval.keydoc import scale_key_markdown
from twinval.registry import load_registry

root = Path(__file__).resolve().parents[1]
reg = load_registry()
text = ["# Scale and sub-scale key", "",
        "Generated from `src/twinval/registry/instruments.yaml` (`python tools/make_scale_key.py`).",
        "Column *items* lists the item codes **in the order of the `it1 ... itK` columns** of the paired data files.",
        "*Redundant* scales (facets nested in a domain, total scores) are analysed at the scale level but are excluded from the",
        "profile, ESEM and item-level analyses, where they would count the same items twice.", "",
        scale_key_markdown(reg)]
(root / "docs" / "SCALE_KEY.md").write_text("\n".join(text), encoding="utf-8")
print("written", len(text))
