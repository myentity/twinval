# Simulation study (Monte Carlo)

Independent replications of the twinval classification rules over a grid of twin processes. It complements the synthetic
worked example in `examples/synthetic`, which is a single constructed data set.

* `DESIGN.md` — ADEMP design (aims, data-generating mechanisms, estimands, methods, performance measures) and the
  population values.
* `conditions.py` — the 53 conditions and the scenario builder. Data are generated with `twinval.synth.simulate`, the same
  generator as `twinval simulate`, with fixed item parameters via `item_seed`.
* `mc_core.py` — the twinval decisions for one data set (same code path as `twinval run`).
* `population.py` — the three bases of comparison: generator properties with the implied Wu status, population score values (N = 200,000) and the large-sample reference classification (N = 20,000).
* `population/` — the population values used as estimands.
* `run_replications.py` — independent replications (one row per scale and replication in `replications/<condition>.csv`;
  resumable; `--order-check` runs the row-order subset into `replications_order_check/`).
* `summarize.py` — performance measures (`summary/`): bias, RMSE and CI coverage; estimation and step-decision rates;
  agreement with the expected and the reference classification with Monte Carlo SEs; sensitivity / specificity; D/F
  with cause; harmonisation.

Requirements: twinval installed (`pip install -e .`) and the R backend (`twinval check-r`).

```bash
cd simulation_study
python population.py                                   # about 3 minutes on 7 parallel processes
python run_replications.py --reps 300 --workers 10     # main grid (hours; resumable)
python run_replications.py --reps 100 --workers 10 --order-check
python summarize.py
```

Working data go to `work_*/`; they are not versioned and can be deleted.

Note on flags in `replications/`: the `cfa_flags` / `structural_flags` columns were written before the
`identification_scheme_sensitive` flag was split into `identification_scheme_sensitive` (W2 and M1 both evaluated and in
conflict) and `identification_route_dependent` (W2 not interpreted after a rejected threshold step, M1 supported). In these
files `identification_scheme_sensitive` therefore covers both cases. Flags never change a type, and no performance
measure in `summary/` uses them; the step decisions (`W1`, `W2`, `W3`, `M1`) are recorded and allow the current flags to
be recomputed.
