# Monte Carlo evaluation of the twinval classification rules — design (ADEMP)

This study responds to reviewer comment 4. The synthetic worked example (`examples/synthetic`) shows how the procedure
behaves on one constructed data set. It does not estimate how often the procedure reaches a given classification when
data are regenerated; this study estimates that with independent replications over a grid of conditions.

The design follows Morris, White & Crowther (2019, *Statistics in Medicine* 38:2074–2102): **A**ims, **D**ata-generating
mechanisms, **E**stimands, **M**ethods, **P**erformance measures. The final grid was fixed before any replication was run
(see *Design history*).

## A — Aims

1. Estimate how well the procedure **detects the violations built into the generator**: unequal loadings, shifted item
   intercepts, local dependence, too few items, miscalibration, weak tracking.
2. Estimate how stably it **reproduces its own large-sample classification** (the reference) at realistic sample sizes.
3. Describe the estimation process: convergence, established identification, admissible and near-boundary solutions,
   and the support, rejection and non-evaluability of every invariance step.
4. Estimate bias, RMSE and CI coverage of the scale-level statistics that feed the I and C axes.
5. Estimate how often a substantive structural conclusion is possible at all.
6. Estimate the false-alarm rate of the row-order check.

## D — Data-generating mechanisms

The generator is `twinval.synth.simulate`, the same generator as `twinval simulate`. Its only change is the optional
`item_seed` key, which keeps item parameters fixed across replications. The shipped worked example is reproduced
byte-for-byte.

**One data set** consists of:
- one instrument with three scales of k items each (k = 5 unless stated);
- 5 response categories;
- human factors correlated at .30;
- the same twin process on all three scales.

**Fixed item parameters:**
- human loadings .78, .70, .62, .74, .66, or the first k of them;
- thresholds of a unimodal category distribution, jittered per item and fixed by `item_seed = 20261027`.

Twin factor: η_T = μ + σ(ρ η_H + √(1−ρ²) ζ). Twin items use the human thresholds.

**Factors:**

| Factor | Levels |
|---|---|
| N (human–twin pairs) | 250 / 500 / 1000 |
| ρ — tracking (latent) | .30 / .50 / .60 / .625 / .65 / .70 / .80 / .90 |
| μ — twin mean shift (bias), latent SD | 0 / .5 / .6 / 1 |
| σ — twin factor scale (variance compression) | 1 / .6 |
| Loading difference | none / twin loadings .45, .90, .40, .90, .45 |
| Intercept shift | none / twin items 1 and 5 shifted by +.5 and −.5 latent SD |
| Rare categories | none / twin never uses the highest category (triggers harmonisation) |
| Duplicated twin answers | none / twin items 2 and 3 copy item 1 (r = .98 before categorisation) |
| Items per scale | 5 / 3 / 2 |

**Fractional design: 53 conditions over 21 twin processes.** The base condition is ρ = .9, μ = 0, σ = 1, k = 5, with no
violation.

| Block | Conditions |
|---|---|
| I_tracking | ρ ∈ {.3, .5, .7, .8, .9} × N ∈ {250, 500, 1000} — 15 |
| I_threshold | ρ ∈ {.60, .625, .65} × N — 9; observed population r ≈ .48 / .50 / .52, i.e. at the I cut-off |
| C_calibration | (μ, σ) ∈ {(.5, 1), (.5, .6), (.6, .6), (1, 1)} × N — 12; population C_b ≈ .91, .82, .77, .71 |
| S_loadings, S_thresholds | loading difference × N, intercept shift × N — 6 |
| R_rare | never-highest with μ ∈ {0, .5}, N = 500 — 2 |
| X_combined | ρ = .5 + loading difference; σ = .6 + loading difference, N = 500 — 2 |
| P_duplicates | duplicated twin answers, N = 500 — 1 |
| K_short | k ∈ {3, 2} × N — 6 |

**Replications and seeds.**
- R = 300 independent replications per condition.
- Replication seed = 10⁶ + 1000 · condition index + replication number.
- Item parameters stay fixed.

**Scope.** Results are conditional on one set of item parameters and on one architecture: three scales per instrument,
5 categories. Two- and three-item scales are covered by the K_short block, but the full variety of the empirical battery
(different lengths and category numbers within one instrument, near-saturated models with correlated residuals) is not.

## E — Estimands: three bases of comparison

1. **Generator properties (known by construction).** For each process we record whether loadings, threshold spacing and
   item intercepts are equal between human and twin, whether twin items are locally dependent, the number of items and
   whether the twin omits a category. The expected status of the **primary Wu/Estabrook sequence** follows from these
   properties by the mapping below. Millsap/Liu plays no part in it.

   The Wu/Estabrook steps are fitted with theta parameterisation and `ID.cat = Wu.Estabrook.2016`:
   - **W0 (configural):** thresholds free on both sides; item intercepts fixed at 0 and residual variances at 1;
     residual covariance of the same item human–twin free.
   - **W1 (thresholds):** thresholds equal; the twin intercepts and residual variances are freed. W1 therefore tests only
     the *spacing* of the thresholds: with 5 categories, 4 equalities minus 2 freed parameters gives 2 constraints per item.
   - **W2 (loadings):** W1 + equal loadings.
   - **W3 (scalar):** W2 + equal intercepts.

   | Generator | Expected Wu decisions | Expected S |
   |---|---|---|
   | no violation | W1, W2, W3 supported | scalar invariance |
   | intercept shift (items 1, 5) | the location shift of x* is absorbed by the freed twin intercept at W1; W1, W2 hold, W3 fails | metric invariance |
   | twin loadings differ | x* keeps unit variance, so threshold spacing is unchanged; W1 holds, W2 fails | configural fit |
   | 2-item scales | the paired W0 has negative df | model not estimable (identification not established) |
   | duplicated twin answers | local dependence; the one-factor twin model is misspecified; estimability depends on the data | no expectation |

   Exact equality in the generator is not the same as support under the Δ-criteria (ΔCFI ≥ −.010, ΔRMSEA ≤ .015,
   ΔSRMR ≤ .030 / .010). A small real difference can stay below the practical cut-offs without the procedure being wrong.
   The violations were therefore sized so that they exceed the cut-offs with nearly unlimited data, which basis 3 checks.

2. **Population score values.** Pearson r, C_b, CCC and the scale-level AUC of the scale scores are computed on an
   independent sample of N = 200,000 (AUC from 2 million random human–twin pairs). These are the targets for bias, RMSE
   and coverage. Population axes: I_pop = (r ≥ .50), because with r > 0 the population AUC exceeds .50; C_pop = (C_b ≥ .80).

3. **Large-sample reference classification.** The rules (rules_version 2) are applied to one sample of N = 20,000. This
   is the procedure's own verdict with nearly unlimited data, **not an independent truth**.

**Two outcome families follow:**
- **detection** — agreement with the expected status / type from bases 1 + 2;
- **stability** — agreement with the reference status / type from bases 2 + 3.

For the I axis they coincide: the axis is a cut-off on a population quantity.

**Results of step 1 (`population.py`).** Mean of the three scales per process; full results are in
`population/population_truth.csv`.

| Twin process | Block | r | C_b | AUC | S expected (generator, Wu) | S reference (N = 20,000) | Type expected | Type reference |
|---|---|---|---|---|---|---|---|---|
| `rho0.90_mu0.5_sg0.6` | C_calibration | 0.624 | 0.818 | 0.643 | scalar_invariance | scalar_invariance | A | A |
| `rho0.90_mu0.5_sg1.0` | C_calibration | 0.720 | 0.909 | 0.676 | scalar_invariance | scalar_invariance | A | A |
| `rho0.90_mu0.6_sg0.6` | C_calibration | 0.623 | 0.772 | 0.635 | scalar_invariance | scalar_invariance | B | B |
| `rho0.90_mu1.0_sg1.0` | C_calibration | 0.712 | 0.709 | 0.623 | scalar_invariance | scalar_invariance | B | B |
| `rho0.60_mu0.0_sg1.0` | I_threshold | 0.480 | 1.000 | 0.611 | scalar_invariance | scalar_invariance | C | C |
| `rho0.625_mu0.0_sg1.0` | I_threshold | 0.500 | 1.000 | 0.618 | scalar_invariance | scalar_invariance | A | A/C |
| `rho0.65_mu0.0_sg1.0` | I_threshold | 0.520 | 1.000 | 0.624 | scalar_invariance | scalar_invariance | A | A |
| `rho0.30_mu0.0_sg1.0` | I_tracking | 0.239 | 1.000 | 0.550 | scalar_invariance | scalar_invariance | C | C |
| `rho0.50_mu0.0_sg1.0` | I_tracking | 0.400 | 1.000 | 0.589 | scalar_invariance | scalar_invariance | C | C |
| `rho0.70_mu0.0_sg1.0` | I_tracking | 0.561 | 1.000 | 0.637 | scalar_invariance | scalar_invariance | A | A |
| `rho0.80_mu0.0_sg1.0` | I_tracking | 0.642 | 1.000 | 0.665 | scalar_invariance | scalar_invariance | A | A |
| `rho0.90_mu0.0_sg1.0` | I_tracking | 0.723 | 1.000 | 0.699 | scalar_invariance | scalar_invariance | A | A |
| `rho0.90_mu0.0_sg1.0_k2` | K_short | 0.604 | 1.000 | 0.651 | model_not_estimable | model_not_estimable | model_not_estimable | model_not_estimable |
| `rho0.90_mu0.0_sg1.0_k3` | K_short | 0.640 | 1.000 | 0.664 | scalar_invariance | scalar_invariance | A | A |
| `rho0.90_mu0.0_sg1.0_dup` | P_duplicates | 0.689 | 0.994 | 0.680 | — (no expectation) | model_not_estimable | — | model_not_estimable |
| `rho0.90_mu0.0_sg1.0_rare` | R_rare | 0.706 | 0.969 | 0.692 | scalar_invariance | scalar_invariance | A | A |
| `rho0.90_mu0.5_sg1.0_rare` | R_rare | 0.685 | 0.933 | 0.687 | scalar_invariance | scalar_invariance | A | A |
| `rho0.90_mu0.0_sg1.0_dl` | S_loadings | 0.696 | 0.996 | 0.689 | configural_fit | configural_fit | E | E |
| `rho0.90_mu0.0_sg1.0_sh` | S_thresholds | 0.722 | 1.000 | 0.699 | metric_invariance | metric_invariance | A | A |
| `rho0.50_mu0.0_sg1.0_dl` | X_combined | 0.385 | 0.997 | 0.587 | configural_fit | configural_fit | C | C |
| `rho0.90_mu0.0_sg0.6_dl` | X_combined | 0.585 | 0.934 | 0.650 | configural_fit | configural_fit | E | E |

- **Expected vs reference S.** They agree for all 20 processes with an expectation. For duplicated twin answers the
  reference verdict is "not estimable": the baseline CFA did not converge and the paired ESEM was inadmissible.
- **Reference types over the 53 conditions:**
  - A: 26;
  - C: 10;
  - B: 6;
  - E: 4;
  - not estimable: 4;
  - at the cut-off: 3.
- **At the cut-off (ρ = .625).** The population r of the three scales is .4997–.5003, so the population I axis is not
  defined there. These conditions are reported as the probability of classifying I as adequate, which is expected to be
  near .50.

**Types D and F.** In this grid the procedure's ability to recognise true D and F cases is not evaluated. No twin process
has a D or F expectation: D needs a CFA without conclusion together with an admissible ESEM discrepancy; F is a residual
category. If D or F is recovered, it counts as a mismatch in the overall agreement. Its frequency is reported separately,
and so is its cause:
- the invariance decision (a step rejected or supported against the expectation);
- loss of CFA admissibility;
- ESEM discrepancy (φ < .85 with both ESEM solutions admissible);
- a combination outside the rules.

## M — Methods evaluated

In every replication the procedure runs exactly as in `twinval run` (`mc_core.py`):

- **Scale level:** Pearson r, scale-level AUC, CCC and C_b, each with a person-level bootstrap CI (B = 2,000; the study
  used 5,000). I axis = r ≥ .50 and lower AUC CI > .50; C axis = C_b ≥ .80.
- **Structural level:**
  - category harmonisation (its frequency and the number of categories kept are recorded per replication, because with
    rare categories different N can lead to different recodings);
  - paired ordinal CFA, both sequences;
  - ESEM with target rotation + Procrustes/Tucker point estimate;
  - admissibility gate;
  - rules_version 2 status;
  - final type.
- **Not repeated:** the Tucker bootstrap and the profile level, neither of which enters the type. As a result, the
  `bootstrap_unstable` flag is not evaluated.
- **Row-order check:** 100 replications in each of 7 conditions (base, both S blocks, both R blocks, P, K_short k = 2;
  N = 500), 3 random row orders each.

## P — Performance measures

**Non-estimable replications stay in every main denominator.** Accuracy conditional on an estimable model is reported
separately and labelled as such.

| Component | Measures |
|---|---|
| r, CCC, C_b, AUC | bias and RMSE against the population values; coverage of the 95% bootstrap CIs (for C_b only where the population C_b < .99: C_b ≤ 1 by construction, so at C_b = 1 no percentile interval of the estimate can contain the population value) |
| Structural estimation | rates of convergence, established identification, admissible and near-boundary solutions (baseline and every step); misfitting baseline; inadmissible paired ESEM |
| Invariance | per step (W1, W2, W3; Millsap M1): rates of supported / rejected / not evaluable / not interpreted |
| I and C axes | rate of correct decisions per axis; for the I_threshold block, P(I adequate) as a function of r |
| Structural status | expected vs recovered and reference vs recovered S: rate of agreement and transition matrix |
| Final type | expected → recovered and reference → recovered transition matrices; rate of agreement with MCSE; per type, sensitivity and specificity over the grid by N; D and F reported separately with their cause |
| Availability of a conclusion | share of scales with a substantive structural status (scalar, metric, configural fit, congruence, discrepant) |
| Row-order check | share of scales flagged order-sensitive when the reference baseline is admissible (false alarms) |
| Harmonisation (R block) | rate of harmonisation and number of categories kept, by N |

**Monte Carlo uncertainty.** The three scales of a replication share one data set. Each proportion is therefore first
averaged within the replication, and its empirical MCSE is SD(replication means) / √R.

With R = 300 and a proportion of .5 under independence the MCSE is √(.25/300) = .029, so a 95% Monte Carlo interval is
about ±5.7 percentage points. Near 0 or 1 it is much narrower, e.g. ±2.5 points at .95.

This precision distinguishes large differences between conditions. It does not reliably distinguish, for example, 92%
from 95% recovery, and the text will not draw such distinctions. If a comparison of that size becomes necessary, R will
be increased for the conditions concerned, fixed in advance from the required MCSE.

## Design history (all changes made before any replication was run)

| Version | Change | Reason |
|---|---|---|
| v0 | 34 conditions; calibration (μ, σ) ∈ {(.5, 1), (0, .6), (.5, .6)} | initial grid |
| v1 | calibration (μ, σ) ∈ {(.5, 1), (.5, .6), (.6, .6), (1, 1)}; duplicated-answers condition added | population values of v0: all calibration processes had C_b ≥ .82, so no process had a true type B |
| v2 (final) | three conditions at the I cut-off (ρ = .60 / .625 / .65 → r ≈ .48 / .50 / .52) × N; short-scale block (k = 3, 2) × N; generator properties separated from the large-sample reference; explicit Wu mapping; extended performance measures | review of the protocol |

Only population values (basis 2) and the reference (basis 3) were used to revise the grid; no replication results.

## Computing budget

Measured in the pilot run (10 replications in each of the 53 conditions, 10 parallel workers, no failed replication):
- 520 replications took 18 minutes;
- typical time per replication under load is 14–22 s;
- 2-item scales take about 45 s (all models are attempted although the baseline is not identified).

| Part | Approximate time |
|---|---|
| Main grid: 53 × 300 replications | ≈ 9–10 h on 10 parallel workers |
| Row-order subset: 7 × 100 replications, 4 CFA fits each | ≈ 1.5–2 h |

Runs are resumable: each replication writes one row per scale, and finished replications are skipped. The pilot
replications are replications 1–10 of the main series: same seeds, same code. No design decision was taken after the
pilot.
