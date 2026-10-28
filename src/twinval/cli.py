"""Command line interface: ``twinval <command>``."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import load_config
from .io import DataError, load_study
from .pipeline import STAGES, Pipeline
from .registry import load_registry


def _common(p):
    p.add_argument("--config", help="YAML file with settings that override the defaults")
    p.add_argument("--data", help="folder with <scale_id>_paired_wide.csv files (overrides data.data_dir)")
    p.add_argument("--registry", help="custom instruments.yaml (default: the packaged registry)")
    p.add_argument("-v", "--verbose", action="store_true")


def _cfg(args):
    over = {}
    if getattr(args, "data", None):
        over["data"] = {"data_dir": args.data}
    if getattr(args, "instruments", None):
        over.setdefault("data", {})["instruments"] = [int(x) for x in args.instruments.split(",")]
    if getattr(args, "seed", None) is not None:
        over["run"] = {"seed": args.seed}
    if getattr(args, "fast", False):
        over.update({"profile": {"n_boot": 300, "n_perm": 50}, "scale": {"n_boot": 300, "n_perm": 50, "cv_repeats": 3, "n_nested_boot": 100,
                                                                       "n_boot_disattenuation": 300},
                     "item": {"n_boot": 300, "n_boot_variance": 200}, "structural": {"n_boot_tucker": 10},
                     "synthesis": {"n_boot_cluster": 500}})
    return load_config(args.config, over)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="twinval", description="Multi-level validation of LLM-based digital twins")
    ap.add_argument("--version", action="version", version=f"twinval {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the pipeline (all stages or a subset)")
    _common(r)
    r.add_argument("--out", required=True, help="output folder")
    r.add_argument("--stages", default=",".join(STAGES), help=f"comma separated subset of: {', '.join(STAGES)}")
    r.add_argument("--instruments", help="comma separated instrument ids, e.g. 148,161")
    r.add_argument("--seed", type=int)
    r.add_argument("--fast", action="store_true", help="small resampling sizes for a quick smoke run (NOT for reporting)")

    v = sub.add_parser("validate", help="check the input files without running any analysis")
    _common(v)
    v.add_argument("--instruments")

    k = sub.add_parser("scale-key", help="print the scale / sub-scale key of the registry")
    k.add_argument("--registry")
    k.add_argument("--format", choices=["md", "csv"], default="md")

    sm = sub.add_parser("simulate", help="generate synthetic paired data (registry + data + scenario truth) from a scenario YAML")
    sm.add_argument("--scenario", required=True, help="scenario YAML, e.g. examples/synthetic/scenario.yaml")
    sm.add_argument("--out", required=True, help="output folder (gets registry.yaml, data/, scenario_truth.csv)")
    sm.add_argument("--seed", type=int, help="override the seed of the scenario")

    ct = sub.add_parser("compare-truth", help="compare the designed (intended) scale types of a synthetic scenario with a run's typology")
    ct.add_argument("--truth", required=True, help="scenario_truth.csv written by `simulate`")
    ct.add_argument("--results", required=True, help="results folder of `twinval run`")
    ct.add_argument("--out", help="write the comparison as CSV")
    ct.add_argument("--rules-version", type=int, help="typology rules of the run (default: read from its config_used.yaml)")

    sub.add_parser("check-r", help="check the optional R backend (Rscript + lavaan + semTools)")
    sub.add_parser("write-config", help="print the default configuration").add_argument("--out")

    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")

    if args.cmd == "simulate":
        from .synth import load_scenario, simulate
        spec = load_scenario(args.scenario)
        if args.seed is not None:
            spec["seed"] = args.seed
        res = simulate(spec, args.out)
        print(f"synthetic data written to {Path(args.out).resolve()}: {len(res['truth'])} scales; run with\n"
              f"  twinval run --data {Path(args.out) / 'data'} --registry {Path(args.out) / 'registry.yaml'} --out <results>")
        return 0
    if args.cmd == "compare-truth":
        import pandas as pd
        from .synth import compare_with_truth
        rv = args.rules_version
        if rv is None:
            import yaml
            used = Path(args.results) / "config_used.yaml"
            run_cfg = yaml.safe_load(used.read_text(encoding="utf-8")) if used.exists() else {}
            rv = int(((run_cfg or {}).get("structural") or {}).get("rules_version", 1))   # runs before rules v2 have no key
        cmp = compare_with_truth(pd.read_csv(args.truth), pd.read_csv(Path(args.results) / "05_synthesis" / "scale_typology.csv"), rv)
        if args.out:
            cmp.round(3).to_csv(args.out, index=False)
        print(cmp.drop(columns="design").round(2).to_string(index=False))
        print(f"\nintended type recovered for {int(cmp['match'].sum())} of {len(cmp)} scales")
        return 0
    if args.cmd == "check-r":
        from .structural import r_available
        ok, msg = r_available()
        print("R backend OK" if ok else f"R backend NOT available: {msg}")
        return 0 if ok else 1
    if args.cmd == "write-config":
        from .config import DEFAULT_CONFIG_PATH
        txt = Path(DEFAULT_CONFIG_PATH).read_text(encoding="utf-8")
        if args.out:
            Path(args.out).write_text(txt, encoding="utf-8")
        else:
            print(txt)
        return 0
    if args.cmd == "scale-key":
        reg = load_registry(args.registry)
        if args.format == "csv":
            print("instrument_id,instrument,scale_id,scale_name,level,scoring,n_items,redundant,parent,items")
            for i in reg.instruments.values():
                for s in i.scales:
                    print(f'{i.id},"{i.name_en}",{s.id},"{s.name_en}",{s.level},{s.scoring},{s.n_items},{s.redundant},{s.parent or ""},{" ".join(s.items)}')
        else:
            from .keydoc import scale_key_markdown
            print(scale_key_markdown(reg))
        return 0

    cfg = _cfg(args)
    reg = load_registry(getattr(args, "registry", None))
    if args.cmd == "validate":
        try:
            st = load_study(cfg.data.data_dir, reg, cfg.data.instruments, cfg.data.id_column, cfg.data.legacy_file_names, cfg.data.missing)
        except DataError as e:
            print(f"DATA ERROR: {e}")
            return 2
        print(f"OK: {st.n} complete respondents, {len(st.scale_ids)} scales, instruments {st.instrument_ids}")
        print(st.audit.to_string(index=False))
        return 0

    pipe = Pipeline(cfg, args.out, reg)
    try:
        pipe.run([s.strip() for s in args.stages.split(",") if s.strip()])
    except DataError as e:
        print(f"DATA ERROR: {e}", file=sys.stderr)
        return 2
    print(f"done - results in {Path(args.out).resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
