# Methods reference

Formulas and conventions implemented in `twinval`. Notation: `H`, `T` are human and twin scores; `N` respondents;
`z = (X - mean_H) / SD_H` is the **human-based standardisation** applied to *both* sides of a scale, so systematic twin bias
stays visible (twin scores are never re-centred on their own mean, except in the explicitly "distinctive" profile variant).

## Scale set and weights
* All scales are analysed at level 2. Profile, ESEM and item level use the **non-redundant** scales (facets nested in a domain and
  total scores are excluded, registry flag `redundant`).
* Profile statistics give every instrument the same weight: a scale of instrument `m` with `k_m` scales gets weight `1/(k_m M)`.

## Level 1 - profile (`profile.py`)
* **Identification.** `D_ij = (1/M) sum_m sqrt( (1/k_m) sum_{s in m} (zH_is - zT_js)^2 )` (rows: human i, columns: twin j). Rank of the own
  twin, Top-1/5/10, per-person AUC (ties count 1/2, own cell excluded), person-level bootstrap CIs (joint resampling of respondents),
  label-permutation test `p = (1 + #{null >= obs}) / (B + 1)`. *Full* and *distinctive* (per-scale sample means of humans and twins removed) profiles.
  Robustness table: Pearson similarity on all/non-redundant scales, equal-weight Fisher-z Pearson, and the final RMS distance. Leave-one-instrument-out.
* **Calibration** in the same specification: Lin's CCC of a person's human and twin profile (own) and averaged over all other twins (cross), weighted
  covariance with weights `1/(k_m M)`; MAE and RMSE of the profile (full) and after removing the mean twin bias of each scale (= distinctive);
  mean absolute bias and RMS bias across scales; bias per scale with bootstrap CI.
* **Mahalanobis** distance with a Ledoit-Wolf shrunk covariance of the human z-profiles, and a comparison of the *mean* human and twin profiles
  (Procrustes disparity, Euclidean distance, pattern r; raw units, sum- and mean-scored scales differ in unit).

## Level 2 - scale (`scale.py`)
* Pearson r, Spearman rho (BH-adjusted p across scales); scale-level AUC from `|zH_i - zT_j|` with permutation p; CCC and its decomposition
  `CCC = r * C_b`, `C_b = 2 / (v + 1/v + u^2)`, `v = SD_T/SD_H`, `u = (mean_T - mean_H)/sqrt(SD_H SD_T)` (computed directly, not as CCC/r).
  All CIs from one set of person-level bootstrap resamples. Descriptive two-component classification by r and C_b (`Good reconstruction`, `Ordering preserved / inadequate
  calibration`, `Limited individuation`, `Combined limitations`).
* **CV recalibration.** Repeated (10 x) 5-fold CV; location-scale `b = SD_H,train/SD_T,train` (corrects location and scale only) versus OLS `b = cov/var(T_train)`
  (also shrinks by r). Out-of-fold CCC/C_b/MAE/RMSE (MAE/RMSE also in units of the human SD). Uncertainty: nested person-level bootstrap in which **all copies of
  a resampled respondent are kept in one fold** (otherwise duplicated respondents leak between training and test folds).
* **Disattenuation** `r_dis = r / sqrt(rel_H rel_T)`: Cronbach's alpha (k >= 3, primary), Spearman-Brown for two-item scales, composite reliability from the
  configural CFA (sensitivity). Values above 1 are reported as they are (never clipped); bootstrap CIs recompute r and reliabilities in every replicate.

## Level 3 - structure (`structural.py`, `r/`)
* **Harmonisation** (`prep.harmonize`): if the observed category sets of Human and Twin differ for an item, adjacent categories are merged with one shared recoding
  until every category occurs on both sides. Items that collapse to one category make the scale unusable for structural models.
* **Paired ordinal CFA** (`r/paired_cfa_invariance.R`, lavaan/semTools, theta parameterisation, WLSMV): Human_i and Twin_i are two measurements of one person
  (two correlated factors, free residual covariance of the same item on both sides).
  Scheme W (Wu & Estabrook 2016): W0 configural -> W1 thresholds -> W2 loadings -> W3 intercepts (scalar) -> W4 residuals.
  Scheme M (Millsap & Tein 2004; Liu et al. 2017; marker item `it1`): M0 baseline -> M1 loadings -> M2 thresholds -> M3 unique variances.
  **The two schemes rest on different identification constraints, are not nested and are never merged into one ladder.**
  Besides fit indices the R script writes, for every model, the diagnostics used below: scaled chi-square and df, CFI/TLI,
  RMSEA with its 90% CI, SRMR, smallest eigenvalue of the parameter covariance matrix (duplicated columns of equality-constrained
  parameters removed), of the residual (theta) and latent (psi) covariance matrices, largest |standardised loading|,
  largest |residual correlation| and |latent correlation|, and lavaan's `post.check`. Also extracted: CR, AVE, latent correlation,
  Tucker's phi of CFA loadings.
* **Structural decision rules** (`structural.py`: `model_status`, `baseline_fit`,
  `step_decision`, `cfa_scale_status`, `structural_evidence_level`). Three properties are kept apart:
  1. *Estimability* of every fitted model: `nonconverged`; `identification_not_established` (df < 0, or the parameter covariance
     matrix cannot be computed or is not positive definite beyond `cfa_pd_tolerance`); `inadmissible` (Heywood case, residual or
     latent covariance matrix not positive definite, |residual or latent correlation| > 1 + `cfa_bound_tolerance`); `admissible`.
     Estimates close to the boundary are **warnings, not exclusions**: flags `near_boundary_loading` (|standardised loading| >= .95),
     `near_boundary_residual_correlation`, `near_boundary_latent_correlation` (|r| >= .95).
  2. *Absolute fit* of an admissible baseline (W0): supported if SRMR <= `cfa_srmr_fit_max` (.08; an operational criterion for
     this framework, not a universal cut-off). CFI/TLI and the RMSEA CI are diagnostics: SRMR acceptable but CFI < .90 or lower RMSEA
     CI bound > .10 gives the flag `conflicting_fit_diagnostics`. A just-identified baseline (df = 0) is not blocked: its fit is
     `not_testable` (flag `baseline_fit_not_testable`) and the later, over-identified steps can still be tested.
  3. *Invariance*, step by step along the primary sequence (Wu). A step is `supported` if delta CFI >= -.010, delta RMSEA <= .015
     **and** delta SRMR <= .030 (steps adding loading constraints) or <= .010 (threshold / intercept steps). It is *not evaluable*
     (never "rejected") if either model of the comparison is not admissible, if the step adds no testable constraints (df difference
     <= 0) or if an index is missing; a step after a step that is not supported is not interpreted. Note: under the theta
     parameterisation as implemented here (W1 equates the thresholds and frees the twin intercepts and residual variances), the
     threshold step gave delta SRMR = 0 in every comparison of the study, so the SRMR condition did not affect any W1 decision.
  The **CFA conclusion** is `scalar`, `metric` or `configural_fit` only from an admissible baseline whose fit is supported (or not
  testable with a supported later step). A threshold step that is supported without a supported loading step gives
  `configural_fit` plus the flag `threshold_supported` (not a separate level). Invariance steps supported under a misfitting baseline
  give no conclusion (flag `invariance_supported_under_baseline_misfit`); an unidentified baseline gives no conclusion (flag
  `baseline_not_identified`). The Millsap/Liu sequence is computed as a sensitivity result and never replaces a Wu step; if its
  loading decision was evaluated and disagrees with the evaluated Wu loading decision (W2), the scale is flagged
  `identification_scheme_sensitive` (a direct conflict of the two tests). If the Wu sequence gives a CFA conclusion but does
  not interpret its loading step because the threshold step was rejected, while the Millsap/Liu loading step is supported,
  the flag is `identification_route_dependent`: the available conclusion depends on the sequence, but the tests do not
  contradict each other. A missing conclusion because of baseline misfit or estimability triggers neither flag.
* **Order sensitivity check** (`structural.order_check_permutations`): refits the CFA of every candidate scale (baseline not admissible,
  or first step of the primary sequence not evaluable) on random row orders and flags scales whose baseline admissibility or CFA
  conclusion is not identical across fits; the observed baseline statuses are listed in `baseline_statuses`.
* **ESEM** (`r/esem_fit.R`): `lavaan::efa()` with target rotation (the theoretical simple structure), ordinal WLSMV, separately for Human and Twin, same target and same
  factors on both sides. **Admissible** = converged, no negative residual variance, no |factor correlation| > 1.0001.
* **Tucker's phi**: Twin loadings are orthogonally Procrustes-rotated onto Human loadings, phi per factor; Hungarian matching is *verified* against the theoretical order.
  Bootstrap (`r/tucker_bootstrap.R`): every replicate re-fits both ESEMs; inadmissible replicates are *not* replaced; a CI is "stable" if >= 80% of replicates are valid.
  The CI of the instrument mean phi uses the per-replicate mean over factors.
* **Admissibility gate**: Tucker's phi is *interpreted* (`tucker_phi_primary`) only if BOTH ESEM solutions of the instrument are admissible; all point estimates are kept as
  `tucker_phi_exploratory` and never used for typology.
* **Structural status S** of a scale: the CFA conclusion if there is one (`scalar_invariance` > `metric_invariance` >
  `configural_fit`); otherwise, if both ESEM solutions of the instrument are admissible, Tucker's phi decides (`congruence_supported`
  if phi >= .85, else `structurally_discrepant`); otherwise an analysis status without structural evidence: `cfa_misfit` (baseline fit
  not supported), `cfa_fit_not_testable`, `model_not_estimable` (CFA and/or ESEM fitted but not estimable), `not_tested`.
  Flags that qualify but do not change S: `cfa_esem_conflict` (CFA conclusion while an interpretable phi < .85), `bootstrap_unstable`
  (ESEM-based S whose Tucker CI rests on < 80% valid replicates), `order_sensitive`, and the CFA flags above.
  `structural.structural_priority: esem` (sensitivity) lets an admissible paired ESEM decide S before the CFA.
* **Validation trap** (instrument level, `structural_instrument_master`): high phi (mean >= .95, min >= .90) from admissible, bootstrap-stable ESEM together with scalar invariance
  not established. Not assigned when configural CFA is unsupported or ESEM is inadmissible.

## Level 4 - item (`item.py`)
* Items of non-redundant scales, response range = **theoretical** range from the registry. Agreement/error (exact and adjacent agreement, MAE/RMSE in units of the range, Spearman rho,
  CCC, quadratic weighted kappa), bias (Wilcoxon signed-rank, BH), compression (SD ratio, normalised entropy, floor/ceiling, category coverage, calibration slope), person-level bootstrap CIs.
  SD ratio and entropy describe the compression of the distribution, not its form; the share of answers in the extreme categories of the range is reported for both sides (`extreme_share_human/_twin`).
* **Variance decomposition** (ANOVA-style crossed): `Y_ij = mu + P_i + I_j + (PI)_ij` on items rescaled to [0, 1]; with one observation per cell the interaction cannot be separated from
  measurement error and is called *person x item / residual*. Human, twin and error (twin - human); paired bootstrap contrast of the percentages; instrument-weighted contrast.
  Absolute variances are reported next to the shares (a share can grow because other components shrink). Similar person shares do not mean that the twin keeps the same
  persons in the same positions: `person_effect_r` is the correlation of the human and twin person effects (person means over the instrument's items).
* **Loading discrepancy**: `|delta loading|` on the target factor and the length of the difference vector of aligned loadings; associations with item fidelity only for instruments with admissible paired ESEM.
* **Reverse-keyed vs direct** items (cluster bootstrap over instruments).

## Level 5 - synthesis (`synthesis.py`)
* Master table with one row per scale (scale-level, structural, item aggregates by scale: medians).
* Spearman associations in pre-defined families: **A** within the scale level, **B** scale x item, **C2** exploratory Tucker (all point estimates), **C1** gated Tucker (descriptive).
  Cluster bootstrap over instruments; bootstrap p-value `2 (b+1)/(n+1)` (never exactly 0), BH within family; with fewer than `min_clusters_for_cluster_bootstrap` instruments CI and p are `NaN`.
* Typology (rules in the README) with the column `evidence_quality`, a short qualifier built from the structural flags
  (identification-sensitive, bootstrap-unstable, conflicting fit diagnostics, CFA-ESEM conflict, order-sensitive, near-boundary
  estimates) that never changes the type. Results matrix: one row per main result of every level. Different units are never
  averaged into one index.
* **Sensitivity analyses of the structural rules** are configuration switches, re-applied to the same model output:
  `structural.primary_branch: millsap`, `structural.use_delta_srmr: false`, `structural.exclude_near_boundary: true`,
  `thresholds.cfa_srmr_fit_max: 0.05`, `structural.structural_priority: esem`.

## Harmonised vs raw item data
In the example data 27 of 75 items have merged categories. Medians of item Spearman rho and CCC are identical in both versions; the largest per-item differences are 0.029 (exact agreement), 0.021 (normalised MAE), 0.015 (Spearman rho), 0.025 (CCC), 0.084 (SD ratio) and 0.107 (entropy difference). The median SD ratio is 0.825 on raw and 0.851 on harmonised categories.
