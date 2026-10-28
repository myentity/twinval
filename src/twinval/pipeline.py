"""Orchestration of the validation stages and of the output folder layout.

Stages (each one can be run on its own; results of earlier stages are re-read from disk when needed)::

    audit       data checks, harmonisation log
    profile     level 1 - identification, calibration, Mahalanobis / mean profile
    structural  level 3 - paired ordinal CFA, ESEM + Tucker phi, admissibility gate      (needs R + lavaan)
    scale       level 2 - r / AUC / CCC / C_b, CV recalibration, disattenuation
    item        level 4 - agreement, bias/compression, variance decomposition, loading discrepancy
    synthesis   level 5 - cross-level master, associations, typology, results matrix
    report      markdown summary

Every stage that produces a figure writes it to a `figures/` subfolder of ITS OWN result folder (e.g. `02_scale/figures/`),
not to one flat top-level folder - each level of the analysis is then self-contained.
"""
from __future__ import annotations

import json
import logging
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import __version__, figures, item as item_mod, prep, profile as profile_mod, scale as scale_mod
from . import structural as S, synthesis as synth
from .config import config_hash
from .io import Study, load_study
from .registry import Registry, load_registry

log = logging.getLogger("twinval")

STAGES = ["audit", "profile", "structural", "scale", "item", "synthesis", "report"]
DIRS = {"audit": "00_audit", "profile": "01_profile", "scale": "02_scale", "structural": "03_structural",
        "item": "04_item", "synthesis": "05_synthesis", "logs": "logs"}
# figures live INSIDE the folder of the stage they belong to (e.g. 02_scale/figures/), not in one flat folder,
# so that every level of the analysis is self-contained


class Pipeline:
    def __init__(self, cfg, out_dir: str | Path, registry: Registry | None = None):
        self.cfg = cfg
        self.out = Path(out_dir)
        self.registry = registry or load_registry()
        self.study: Study | None = None
        self.res: dict[str, pd.DataFrame] = {}
        self.harm: dict | None = None
        self.timings: dict[str, float] = {}
        self.notes: list[str] = []
        for d in DIRS.values():
            (self.out / d).mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ helpers
    def rng(self, stage: str) -> np.random.Generator:
        return np.random.default_rng([int(self.cfg.run.seed), STAGES.index(stage) if stage in STAGES else 99])

    def _w(self, stage: str, name: str, df: pd.DataFrame, index: bool = False):
        path = self.out / DIRS[stage] / f"{name}.csv"
        df.to_csv(path, index=index)
        self.res[name] = df

    def _figdir(self, stage: str) -> Path:
        d = self.out / DIRS[stage] / "figures"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def _get(self, name: str) -> pd.DataFrame | None:
        if name in self.res:
            return self.res[name]
        for d in DIRS.values():
            p = self.out / d / f"{name}.csv"
            if p.exists():
                self.res[name] = pd.read_csv(p)
                return self.res[name]
        return None

    def load(self) -> Study:
        if self.study is None:
            d = self.cfg.data
            self.study = load_study(d.data_dir, self.registry, d.instruments, d.id_column, d.legacy_file_names, d.missing)
            log.info("loaded %d respondents, %d scales, %d instruments", self.study.n, len(self.study.scale_ids),
                     len(self.study.instrument_ids))
        return self.study

    def _harmonized(self):
        if self.harm is None:
            hh, ht, hlog, hbad = prep.harmonize(self.load())
            self.harm = {"items_h": hh, "items_t": ht, "log": hlog, "bad": hbad}
        return self.harm

    def _write_harmonized_data(self, hm):
        """Save the harmonised paired data (same format as the input files) to ``00_audit/harmonized_data``.

        These are the data on which the structural models (CFA / ESEM) and - with ``item.data_version: harmonized``,
        the default and the setting of the original study - the item-level analysis are computed. Scale scores and the
        profile-level analysis always use the ORIGINAL (raw) input data.
        """
        st = self.load()
        d = self.out / DIRS["audit"] / "harmonized_data"
        d.mkdir(parents=True, exist_ok=True)
        for sid in st.scale_ids:
            k = st.registry.scales[sid].n_items
            df = pd.DataFrame({"ID": st.ids})
            for i in range(k):
                df[f"it{i + 1}_h"] = hm["items_h"][sid][:, i].astype(int)
            for i in range(k):
                df[f"it{i + 1}_t"] = hm["items_t"][sid][:, i].astype(int)
            df.to_csv(d / f"{sid}_paired_wide.csv", index=False)

    # ------------------------------------------------------------------ stages
    def stage_audit(self):
        st = self.load()
        self._w("audit", "data_audit", st.audit)
        reg = st.registry
        rows = []
        for sid in st.scale_ids:
            sc = reg.scales[sid]
            h, t = st.H[sid].to_numpy(), st.T[sid].to_numpy()
            rows.append({"scale_id": sid, "scale_name": sc.name_en, "instrument_id": sc.instrument_id, "n_items": sc.n_items,
                         "scoring": sc.scoring, "redundant": sc.redundant, "n": st.n,
                         "min_human": h.min(), "max_human": h.max(), "min_twin": t.min(), "max_twin": t.max(),
                         "zero_variance_human": bool(h.std() == 0), "zero_variance_twin": bool(t.std() == 0)})
        self._w("audit", "scale_score_qc", pd.DataFrame(rows))
        hm = self._harmonized()
        self._write_harmonized_data(hm)
        self._w("audit", "harmonization_log", hm["log"])
        self._w("audit", "harmonization_unusable", hm["bad"])
        self._w("audit", "harmonization_status", prep.harmonization_status(st, hm["log"], hm["bad"]))
        if len(st.instrument_ids) < 4:
            self.notes.append(f"Only {len(st.instrument_ids)} instruments are loaded: cluster (per-instrument) bootstrap CIs, "
                              "equal-instrument-weight statistics and leave-one-instrument-out analyses are "
                              "degenerate or skipped; results are illustrative and not comparable with a full study.")

    def stage_profile(self):
        st = self.load()
        cfg = self.cfg
        if len(st.nonredundant_ids) < 2:
            self.notes.append("profile stage skipped: fewer than 2 non-redundant scales")
            return
        res = profile_mod.identification(st, cfg, self.rng("profile"))
        for k, v in res.items():
            self._w("profile", k, v)
        cal = profile_mod.calibration(st, cfg, self.rng("profile"))
        for k, v in cal.items():
            self._w("profile", k, v)
        mm = profile_mod.mahalanobis_and_mean_profile(st, cfg, self.rng("profile"))
        for k, v in mm.items():
            self._w("profile", k, v)
        figures.profile_ccc_hist(cal["profile_calibration_per_person"], self._figdir("profile") / "profile_ccc_own_vs_cross", cfg.figures.dpi)

    def stage_structural(self):
        st, cfg = self.load(), self.cfg
        if not cfg.structural.enabled:
            return
        ok, msg = S.r_available(cfg.structural.rscript)
        if not ok:
            self.notes.append(f"structural stage skipped: {msg}")
            log.warning("structural stage skipped: %s", msg)
            return
        wd = self.out / DIRS["structural"] / "r_work"
        logs = self.out / DIRS["logs"]
        inp = S.prepare_inputs(st, cfg, wd)
        self.harm = {"items_h": inp["items_h"], "items_t": inp["items_t"], "log": inp["harmonization_log"], "bad": inp["harmonization_unusable"]}
        if len(inp["item_qc"]):
            self._w("structural", "item_ordinal_suitability_qc", inp["item_qc"])
        d3 = self.out / DIRS["structural"]

        cfa_status = None
        if cfg.structural.cfa and len(inp["manifest"]):
            res_path = wd / "cfa_results" / "invariance_results.csv"
            if not (cfg.structural.resume and res_path.exists()):
                S.run_r("paired_cfa_invariance.R", [wd / "cfa_manifest.csv", wd / "cfa_results", cfg.thresholds.cfa_delta_cfi_min,
                                                    cfg.thresholds.cfa_delta_rmsea_max], logs / "cfa.log", cfg.structural.rscript)
            results = pd.read_csv(res_path)
            self._w("structural", "invariance_results", results)
            for k, v in S.parse_cfa(results).items():
                self._w("structural", k, v)
            cfa_status = S.cfa_scale_status(results, cfg)
            self._w("structural", "cfa_scale_status", cfa_status)
            self._w("structural", "cfa_step_diagnostics", S.cfa_step_diagnostics(results, cfg))
            self._w("structural", "cfa_status_by_instrument", S.cfa_by_instrument(st, cfa_status))
            if cfg.structural.order_check_permutations > 0:
                oc = S.cfa_order_sensitivity(results, inp["manifest"], wd, cfg.structural.order_check_permutations,
                                             int(cfg.run.seed), cfg, logs)
                self._w("structural", "cfa_order_sensitivity", oc)
                if len(oc) and oc["order_sensitive"].any():
                    self.notes.append("CFA convergence is order-sensitive for: " + ", ".join(oc.loc[oc.order_sensitive, "scale_id"]) +
                                      " (their structural status depends on the row order of the data; interpret with caution).")

        esem_ids = [e["instrument_id"] for e in inp["esem_instruments"]]
        tucker = gate = None
        if cfg.structural.esem and esem_ids:
            diag_path = wd / "esem_model_diagnostics.csv"
            if not (cfg.structural.resume and diag_path.exists()):
                S.run_r("esem_fit.R", [wd], logs / "esem.log", cfg.structural.rscript)
            diag = pd.read_csv(diag_path)
            for name in ("esem_model_diagnostics", "esem_fit_by_instrument", "esem_loadings_long", "esem_factor_correlations", "esem_thresholds"):
                p = wd / f"{name}.csv"
                if p.exists() and p.stat().st_size > 2:
                    self._w("structural", name, pd.read_csv(p))
            pts = S.tucker_point_estimates(wd, diag, inp["esem_instruments"])
            tucker = pts["tucker_by_factor"]
            self._w("structural", "esem_aligned_loadings", pts["esem_aligned_loadings"])
            self._w("structural", "tucker_summary", pts["tucker_summary"])
            gate = S.admissibility(diag, esem_ids)
            self._w("structural", "esem_admissibility", gate)
            if cfg.structural.tucker_bootstrap and gate["paired_esem_admissible"].any():
                # bootstrap only where a bootstrap-supported conclusion is possible (both ESEMs admissible on the full sample)
                boot_ids = gate.loc[gate.paired_esem_admissible, "instrument_id"].tolist() if cfg.structural.bootstrap_admissible_only else esem_ids
                bwd = wd / "tucker_bootstrap"
                bwd.mkdir(exist_ok=True)
                for sub in ("item_battery", "target_matrices"):
                    (bwd / sub).mkdir(exist_ok=True)
                    for f in (wd / sub).glob("*"):
                        if int(f.name.split("_")[0]) in boot_ids:
                            (bwd / sub / f.name).write_bytes(f.read_bytes())
                pd.DataFrame([e for e in inp["esem_instruments"] if e["instrument_id"] in boot_ids]).to_csv(bwd / "esem_registry.csv", index=False)
                rep_path = bwd / "tucker_bootstrap_replicates.csv"
                if not (rep_path.exists() and cfg.structural.resume and pd.read_csv(rep_path)["rep"].max() >= cfg.structural.n_boot_tucker):
                    S.run_r("tucker_bootstrap.R", [bwd, cfg.structural.n_boot_tucker, cfg.run.seed], logs / "tucker_bootstrap.log", cfg.structural.rscript)
                boot = pd.read_csv(rep_path)
                tucker, mean_ci = S.merge_tucker_bootstrap(tucker, boot, cfg.thresholds.bootstrap_valid_min)
                self._w("structural", "tucker_by_factor", tucker)
                self._w("structural", "tucker_mean_phi_ci", mean_ci)
            else:
                self._w("structural", "tucker_by_factor", tucker)
        table = S.scale_structural_table(st, cfg, cfa_status, tucker, gate, esem_ids, self._get("cfa_order_sensitivity"))
        self._w("structural", "structural_scale_table", table)
        figures.structural_phi_forest(table, cfg.thresholds.phi_congruence, self._figdir("structural") / "tucker_phi_admissibility", cfg.figures.dpi)
        inst = S.instrument_master(st, self._get("cfa_status_by_instrument"), self._get("esem_fit_by_instrument"),
                                   self._get("tucker_summary"), self._get("tucker_mean_phi_ci"), gate, self._get("tucker_by_factor"))
        self._w("structural", "structural_instrument_master", inst)

    def stage_scale(self):
        st, cfg = self.load(), self.cfg
        rng = self.rng("scale")
        core = scale_mod.core_table(st, cfg, rng)
        self._w("scale", "scale_core_table", core)
        per_repeat, per_scale = scale_mod.recalibration(st, cfg, rng)
        self._w("scale", "scale_recalibration_per_repeat", per_repeat)
        self._w("scale", "scale_recalibration_summary", per_scale)
        rel = self._get("cfa_reliability")
        disatt = scale_mod.disattenuation(st, cfg, rng, core, rel)
        self._w("scale", "scale_disattenuation", disatt)
        master = scale_mod.build_master(core, per_scale, disatt, st)
        self._w("scale", "scale_level_master", master)
        self._w("scale", "scale_level_summary", scale_mod.summary_tables(master))
        figures.scale_fidelity_scatter(master, self._figdir("scale") / "scale_r_vs_cb", cfg.figures.dpi)

    def stage_item(self):
        """Item-level analysis.

        DATA VERSION: ``item.data_version`` (default ``harmonized``, as in the original study) decides whether the items are the
        harmonised ones (rare categories merged identically for Human and Twin when their category sets differ) or the raw
        input answers. Every item table carries the column ``data_version``.
        """
        st, cfg = self.load(), self.cfg
        rng = self.rng("item")
        if cfg.item.data_version == "harmonized":
            hm = self._harmonized()
            panel = item_mod.build_panel(st, hm["items_h"], hm["items_t"])
        else:
            panel = item_mod.build_panel(st)
        ag, bias, audit = item_mod.agreement_and_bias(panel, cfg, rng)
        dv = cfg.item.data_version
        tag = lambda df: df.assign(data_version=dv)
        ag, bias, audit = tag(ag), tag(bias), tag(audit)
        self._w("item", "item_agreement_error", ag)
        self._w("item", "item_bias_compression", bias)
        self._w("item", "item_level_data_audit", audit)
        self._w("item", "item_range_codebook", tag(item_mod.range_codebook(panel)))
        decomp, contr = item_mod.variance_decomposition(panel, cfg, rng)
        decomp, contr = tag(decomp), tag(contr)
        self._w("item", "item_variance_decomposition", decomp)
        self._w("item", "item_variance_decomposition_contrasts", contr)
        master = ag.merge(bias.drop(columns=[c for c in bias.columns if c in ag.columns and c != "item_id"]), on="item_id")  # keeps data_version
        master = master.merge(audit[["item_id", "twin_top2_categories_share", "twin_near_degenerate_flag"]], on="item_id")
        load = self._get("esem_aligned_loadings")
        gate = self._get("esem_admissibility")
        if load is not None and len(load) and gate is not None:
            adm = gate.loc[gate.paired_esem_admissible, "instrument_id"].tolist()
            disc, assoc = item_mod.loading_discrepancy(ag, bias, load, adm, cfg, rng)
            self._w("item", "item_loading_discrepancy", disc)
            if len(assoc):
                self._w("item", "item_loading_associations", assoc.assign(per_instrument_rho=assoc["per_instrument_rho"].astype(str)))
            master = master.merge(disc.drop(columns=["instrument_id"]), on="item_id", how="left")
        else:
            self.notes.append("item-level loading discrepancy skipped: no ESEM results (run the structural stage first)")
        self._w("item", "item_level_master", master)
        self._w("item", "item_level_summary", item_mod.summarize(master))
        rk = item_mod.reverse_keyed_analysis(ag, bias, cfg, rng)
        if len(rk):
            self._w("item", "item_reverse_keyed_analysis", rk)
        figures.item_fidelity_map(master, self._figdir("item") / "item_fidelity_map", cfg.figures.dpi)
        figures.variance_decomposition_bars(decomp, self._figdir("item") / "item_variance_decomposition", cfg.figures.dpi)

    def stage_synthesis(self):
        st, cfg = self.load(), self.cfg
        rng = self.rng("synthesis")
        sm, structural, im = self._get("scale_level_master"), self._get("structural_scale_table"), self._get("item_level_master")
        if sm is None:
            raise RuntimeError("synthesis needs the scale stage (scale_level_master.csv)")
        master = synth.build_master(st, sm, structural, im)
        self._w("synthesis", "cross_level_master", master)
        out = synth.associations(master, cfg, rng)
        for k, v in out.items():
            self._w("synthesis", "cross_level_" + k, v, index=(k == "correlation_matrix_exploratory"))
        if "structural_evidence_level" in master and master["structural_evidence_level"].notna().any():
            typ, summ, rules = synth.typology(master, cfg)
            self._w("synthesis", "scale_typology", typ)
            self._w("synthesis", "scale_typology_summary", summ)
            with open(self.out / DIRS["synthesis"] / "scale_typology_rules.yaml", "w", encoding="utf-8") as fh:
                yaml.safe_dump(rules, fh, allow_unicode=True, sort_keys=False)
            figures.typology_map(typ, synth.type_labels(cfg), self._figdir("synthesis") / "scale_typology_map", cfg.figures.dpi)
        else:
            self.notes.append("typology skipped: no structural results (run the structural stage; requires R + lavaan)")
        figures.association_heatmap(out["associations"], self._figdir("synthesis") / "cross_level_associations_heatmap", cfg.figures.dpi)
        prof = {k: self._get(k) for k in ("identification_summary", "profile_calibration_summary")}
        prof = {k: v for k, v in prof.items() if v is not None} or None
        rm = synth.results_matrix(st, prof, sm, structural, self._get("item_level_summary"), self._get("item_variance_decomposition_contrasts"))
        self._w("synthesis", "cross_level_results_matrix", rm)

    def stage_report(self):
        from .report import write_report
        write_report(self)

    # ------------------------------------------------------------------ driver
    def run(self, stages: list[str] | None = None):
        cfg = self.cfg
        order = [s for s in STAGES if stages is None or s in stages]
        for s in order:
            if s in ("profile", "scale", "item", "structural", "synthesis") and not cfg[s].get("enabled", True):
                log.info("stage %s disabled in the configuration", s)
                continue
            t = time.time()
            log.info("=== stage: %s ===", s)
            getattr(self, f"stage_{s}")()
            self.timings[s] = round(time.time() - t, 1)
        if not cfg.output.keep_intermediate:
            self._prune_intermediate()
        self._write_manifest(order)

    def _prune_intermediate(self):
        """Remove intermediate copies of the input data (harmonised paired files, item batteries) from the run folder."""
        import shutil
        w = self.out / DIRS["structural"] / "r_work"
        for rel in ("cfa_input", "item_battery", "order_check", "tucker_bootstrap/item_battery"):
            shutil.rmtree(w / rel, ignore_errors=True)
        # make the run folder relocatable / publishable: no absolute machine paths in the manifests and logs
        root = str(self.out.resolve())
        for f in [w / "cfa_manifest.csv", *(self.out / "logs").glob("*.log")]:
            if f.exists():
                txt = f.read_text(encoding="utf-8", errors="replace")
                if root in txt:
                    f.write_text(txt.replace(root, "<run_dir>"), encoding="utf-8")

    def _write_manifest(self, order):
        import scipy, sklearn, statsmodels
        st = self.study
        mpath = self.out / "run_manifest.json"
        old = json.loads(mpath.read_text(encoding="utf-8")) if mpath.exists() else {}
        # a stage that was only resumed (results reused) keeps the wall-clock time of its original run
        timings = {**old.get("timings_sec", {}), **{k: max(v, old.get("timings_sec", {}).get(k, 0)) for k, v in self.timings.items()}}
        stages_run = [x for x in STAGES if x in set(old.get("stages_run", [])) | set(order)]
        info = {"twinval_version": __version__, "python": platform.python_version(),
                "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
                "statsmodels": statsmodels.__version__, "scikit_learn": sklearn.__version__,
                "config_hash": config_hash(self.cfg), "seed": self.cfg.run.seed, "stages_run": stages_run,
                "timings_sec": timings, "notes": self.notes,
                "data_versions": {"profile_and_scale_level": "raw (input files)", "structural_models": "harmonized",
                                  "item_level": self.cfg.item.data_version},
                "n_respondents": st.n if st else None, "n_scales": len(st.scale_ids) if st else None,
                "instruments": [self.registry.instruments[i].name_en for i in st.instrument_ids] if st else None}
        ok, _ = S.r_available(self.cfg.structural.rscript)
        if ok:
            import subprocess
            p = subprocess.run([self.cfg.structural.rscript, "-e", "cat(as.character(getRversion()), as.character(packageVersion('lavaan')), as.character(packageVersion('semTools')))"], capture_output=True, text=True)
            info["R_lavaan_semTools"] = p.stdout.strip()
        with open(self.out / "run_manifest.json", "w", encoding="utf-8") as fh:
            json.dump(info, fh, indent=2, default=str)
        with open(self.out / "config_used.yaml", "w", encoding="utf-8") as fh:
            yaml.safe_dump(json.loads(json.dumps(self.cfg, default=dict)), fh, sort_keys=False, allow_unicode=True)
