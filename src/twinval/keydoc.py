"""Human-readable scale / sub-scale key generated from the registry (used by ``twinval scale-key`` and the docs)."""
from __future__ import annotations

from .registry import Registry


def scale_key_markdown(reg: Registry, instruments=None) -> str:
    out = []
    for iid, ins in reg.instruments.items():
        if instruments is not None and iid not in instruments:
            continue
        ref = f" - {ins.reference}" if ins.reference else ""
        out += [f"### {iid} - {ins.name_en}{ref}", "",
                f"Response scale {ins.response_min}-{ins.response_max}; item prefix `{ins.item_prefix}`.", "",
                "| scale_id | name | level | scoring | items (it1...itK) | reverse-keyed items | in profile / ESEM / item level |",
                "|---|---|---|---|---|---|---|"]
        for s in ins.scales:
            items = ", ".join(s.items)
            rev = ", ".join(s.reverse_items) or "-"
            note = "no (redundant" + (f": nested in `{s.parent}`)" if s.parent else ": total score)") if s.redundant else "yes"
            out.append(f"| `{s.id}` | {s.name_en} | {s.level} | {s.scoring} | {items} | {rev} | {note} |")
        out.append("")
    return "\n".join(out)
