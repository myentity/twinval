"""Markdown summary of a run (``report.md``)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def _fmt(v, d=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{d}f}"


def write_report(p) -> None:
    """Assemble ``report.md`` from the tables that exist in the output folder."""
    st = p.load()
    L = ["# twinval run report", "",
         f"* respondents (complete pairs): **{st.n}**", f"* scales loaded: **{len(st.scale_ids)}** "
         f"({len(st.nonredundant_ids)} non-redundant)", f"* instruments: **{', '.join(p.registry.instruments[i].name_en for i in st.instrument_ids)}**", ""]
    dv = p.cfg.item.data_version
    L += ["## Data versions", "",
          "* profile and scale level: **raw** input answers (scale scores are computed from the original items)",
          "* structural models (CFA, ESEM): **harmonised** categories (`00_audit/harmonized_data/`, log in `harmonization_log.csv`)",
          f"* item level: **{dv}** (`item.data_version`; `raw` = primary analysis of the study, `harmonized` = sensitivity analysis)", ""]
    if p.notes:
        L += ["## Notes and warnings", ""] + [f"* {n}" for n in p.notes] + [""]

    rm = p._get("cross_level_results_matrix")
    if rm is not None and len(rm):
        L += ["## Multilevel results matrix", "", "Different units of analysis are shown side by side and are never averaged into one index.", "",
              "| level | unit | component | metric | value | 95% CI | detail |", "|---|---|---|---|---|---|---|"]
        for r in rm.itertuples():
            ci = "" if pd.isna(r.ci_lo) else f"[{_fmt(r.ci_lo)}, {_fmt(r.ci_hi)}]"
            L.append(f"| {r.level} | {r.unit} | {r.component} | {r.metric} | {_fmt(r.value)} | {ci} | {r.detail if isinstance(r.detail, str) else ''} |")
        L.append("")

    summ = p._get("scale_typology_summary")
    if summ is not None and len(summ):
        L += ["## Scale typology (rule-based)", "", "| type | n | % | example scales |", "|---|---|---|---|"]
        for r in summ.itertuples():
            k = f"{r.final_type}. " if len(str(r.final_type)) == 1 else ""
            L.append(f"| {k}{r.label} | {r.n} | {r.pct} | {r.example_scales} |")
        L.append("")

    assoc = p._get("cross_level_associations")
    if assoc is not None and len(assoc):
        L += ["## Associations between levels (Spearman)", "", "| family | pair | rho | 95% CI (cluster bootstrap) | p (naive) | q (BH) | n | clusters |", "|---|---|---|---|---|---|---|---|"]
        for r in assoc.itertuples():
            ci = "n/a" if pd.isna(r.ci_lo) else f"[{_fmt(r.ci_lo)}, {_fmt(r.ci_hi)}]"
            L.append(f"| {r.family} | {r.pair} | {_fmt(r.rho)} | {ci} | {_fmt(r.p_naive, 4)} | {_fmt(getattr(r, 'p_BH', np.nan), 4)} | {r.n} | {r.n_clusters} |")
        L.append("")
    L += ["## Output folders", "",
          "Each stage folder holds its own `figures/` subfolder (png + pdf) - figures are never collected in one flat folder.",
          "", "* `00_audit` data checks and harmonisation log", "* `01_profile` respondent-level results, `figures/` profile CCC histogram",
          "* `02_scale` scale-level results, `figures/` r vs C_b scatter",
          "* `03_structural` CFA / ESEM / Tucker / admissibility, `figures/` Tucker's phi forest plot",
          "* `04_item` item-level results, `figures/` item fidelity map, variance decomposition bars",
          "* `05_synthesis` cross-level master, associations, typology, `figures/` typology map, association heat map",
          "* `run_manifest.json`, `config_used.yaml` reproducibility record"]
    (p.out / "report.md").write_text("\n".join(L), encoding="utf-8")
