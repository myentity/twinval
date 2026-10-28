# Synthetic demonstration data

The real example (`examples/data`, two instruments) cannot show every kind of result the framework can produce: with two
instruments the cluster bootstrap and leave-one-instrument-out analyses are degenerate, and only some types of the
typology occur. `twinval simulate` therefore generates **fictional instruments with a twin whose behaviour is built in**.
The synthetic example (`examples/synthetic`) has 10 instruments, 29 scales, 500 respondents, and produces **all eight
outcomes of the typology** (A, B, C, D, E, F, *CFA misfit*, *model not estimable*), a cluster bootstrap over 10 instruments,
leave-one-instrument-out results, harmonisation of categories, reverse-keyed items and a redundant total score.
The data carry no personal information; the instruments do not exist.

```bash
# regenerate the data (already committed in examples/synthetic)
twinval simulate --scenario examples/synthetic/scenario.yaml --out examples/synthetic

# analyse them with the custom registry that `simulate` wrote
twinval run --data examples/synthetic/data --registry examples/synthetic/registry.yaml \
            --config configs/synthetic.yaml --out examples/synthetic/results

# designed versus recovered scale types
twinval compare-truth --truth examples/synthetic/scenario_truth.csv --results examples/synthetic/results
```

## Generative model

Human: instrument *m* has one latent factor per scale, factors correlate with `factor_correlation` (0.3).
Item response `x* = λ·η + √(1-λ²)·e` is cut into ordered categories by fixed item thresholds (a unimodal category
distribution, jittered per item). Twin: `η_T = μ + σ·(√(1-h²)·(ρ·η_H + √(1-ρ²)·ζ) + h·G)`, items are generated from `η_T`
with the *same* thresholds, so a shift of the twin is a real shift of its answers. Data are written in the framework's
input format (reverse-keyed items already reverse-scored). Every instrument has its own random stream
(`[seed, instrument id, instrument seed]`), so editing one instrument never changes another.

| parameter (per scale, under `twin:`) | effect | what it changes in the results |
|---|---|---|
| `rho` | latent correlation twin-human | individuation: Pearson r, AUC, identification |
| `mu`, `sigma` | twin mean shift (SD units), twin factor scale | calibration: C_b, bias, SD ratio (type **B**) |
| `halo` | share of an instrument-wide factor shared by all scales of the instrument | cross-scale contamination in ESEM |
| `loading` (scalar or list), `human.loading` | twin item loadings (default: the human ones) | unequal loadings break metric invariance -> **configural fit** (type **E**) |
| `reverse_loading_scale` | shrinks the twin loadings of reverse-keyed items | reverse-keyed item analysis; metric invariance |
| `item_shift` | moves single twin items up/down | scalar invariance |
| `resid_corr` | the twin reproduces the human item-specific residuals | item-level agreement beyond the trait |
| `duplicate: [[i, j], ...]`, `duplicate_corr` | twin answers item *j* like item *i* | with three identical answers: CFA and/or twin ESEM not estimable |
| `copy_items` | twin reproduces the human answer verbatim (leakage) | perfect item agreement |
| `near_constant`, `near_constant_keep` | twin gives (almost) everybody the same answer to an item | near-degenerate items (item level) |
| `never_lowest`, `never_highest` | twin never uses an extreme category | category harmonisation |

Instrument level: `response: [min, max]`, `key`, `name`, `total: {name, scoring}` (a redundant total score over all items),
`seed`, `item_seed` (fixes item thresholds and default human loadings across data sets that differ only in `seed`; used by the
Monte Carlo study in [`simulation_study/`](../simulation_study)). File level: `seed`, `n_respondents`, `factor_correlation`, `id_prefix`.

## The scenario

`intended` is the type the design aims at; it is a *design*, not a guarantee (`scenario_truth.csv`: column
`intended_type`). `twinval compare-truth` compares it with the recovered type (see
`examples/synthetic/results/expected_vs_recovered.csv`).

| instrument | design | intended types |
|---|---|---|
| 901 faithful twin | high individuation, calibrated, same loadings | A A A |
| 902 miscalibrated twin | good individuation, twin mean shifted / variance compressed | B B B |
| 903 weakly individuating twin | same structure, twin tracks the person weakly (rho 0.25-0.35) | C C C |
| 904 unequal twin loadings | good and calibrated, loadings alternate 0.3/0.9 (no metric invariance) | E E A |
| 905 repeated twin answers (congruent) | items 1-3 answered identically: scale CFA fails, ESEM admissible and congruent | E A A |
| 906 mixed instrument with total | reverse items handled worse (E), shifted twin never using the lowest category (B), weak twin (C), redundant total | E B C (+ total) |
| 907 short scales | 3-4 item scales | A C A |
| 908 single-scale instrument | items 1-4 answered identically: CFA fails and ESEM is not applicable | model not estimable |
| 909 repeated twin answers (discrepant) | items 1-3 identical, items 4-5 ignored: CFA fails, ESEM loadings differ (phi < 0.85) | D F (+ control) |
| 910 repeated twin answers (inadmissible) | items 1-3 identical: CFA and twin ESEM not estimable, the gate closes for the instrument | model not estimable, A, A |

## Result of the published run

`examples/synthetic/results` (settings: `configs/synthetic.yaml`). The typology of the 29
scales: **A 9, B 4, C 5, D 1, E 6, F 1, CFA misfit 1, model not estimable 2** (the redundant total score `906_Total` is typed
with the others; its baseline CFA misfits). The designed type is recovered for **26 of 28** non-redundant scales
(`expected_vs_recovered.csv`); the exceptions are two control scales, both designed A and recovered E:

* `909_Plain` - the Wu loading step is rejected (delta CFI -.009, delta RMSEA .017): borderline;
* `905_Plain2` - the Wu threshold step is rejected (delta RMSEA .021 > .015), so the CFA conclusion is configural fit; the
  Millsap/Liu sequence supports metric invariance, and the scale carries the flag `identification_route_dependent`.

The row-order check refits the five candidate scales on 3 random row orders; none is
order-sensitive. The baseline of `909_DupF` (type F) is nonconverged in some orders and not identified in
others, i.e. not estimable in every fit, so its structural status does not change. Cross-level associations have
cluster-bootstrap intervals over 10 instruments, and `01_profile/identification_leave_one_instrument_out.csv` has 10 rows.

## Read this before interpreting the results

* **Borderline numerics.** The paired ordinal CFA rarely fails on synthetic data. Types D, F, *model not estimable* and the
  congruence route to E need a scale-level CFA that is **not estimable** (here: does not converge), which happens only in
  extreme configurations (three or more identical twin answers) and then depends on the random draw. The seeds of
  instruments 905, 908, 909 and 910 (`seed:` in `scenario.yaml`) were chosen so that the intended outcome occurs; the
  same designs with other seeds give other, equally valid results (`E`, `C`, ...). This mirrors the fragility of the
  real study (see "Order sensitivity" in the README). The committed data files are authoritative; regenerating them on
  another platform can change floating-point details and hence such borderline outcomes.
* Control scales can end up in another type by chance (in the shipped example the control scales `909_Plain` and `905_Plain2`
  are type E because an invariance step is borderline). `expected_vs_recovered.csv` reports such mismatches openly.
* The two identification schemes can also differ numerically: with `structural.primary_branch: millsap` the scale `909_DupD`
  becomes E instead of D, because its Millsap baseline converges while the statistically equivalent Wu configural model does not.
* With ten instruments the cluster bootstrap and leave-one-instrument-out run, but ten clusters are still few:
  the intervals are wide and are meant to illustrate the procedure.
