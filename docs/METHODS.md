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
  All CIs from one set of person-level bootstrap resamples. Descriptive two-component classification by r and C_b (`Good reconstruction`, ...).
* **CV recalibration.** Repeated (10 x) 5-fold CV; location-scale `b = SD_H,train/SD_T,train` (corrects location and scale only) versus OLS `b = cov/var(T_train)`
  (also shrinks by r). Out-of-fold CCC/C_b/MAE/RMSE (MAE/RMSE also in units of the human SD). Uncertainty: nested person-level bootstrap in which **all copies of
  a resampled respondent are kept in one fold** (otherwise duplicated respondents leak between training and test folds).
* **Disattenuation** `r_dis = r / sqrt(rel_H rel_T)`: Cronbach's alpha (k >= 3, primary), Spearman-Brown for two-item scales, composite reliability from the
  configural CFA (sensitivity). Values above 1 are reported as they are (never clipped); bootstrap CIs recompute r and reliabilities in every replicate.

## Level 3 - structure (`structural.py`, `r/`)
* **Harmonisation** (`prep.harmonize`): if the observed category sets of Human and Twin differ for an item, adjacent categories are merged with one shared recoding
  until every category occurs on both sides. Items that collapse to one category make the scale unusable for structural models.
* **Paired ordinal CFA** (`r/paired_cfa_invariance.R`, lavaan/semTools, theta parameterisation, WLSMV): Human_i and Twin_i are two measurements of one person.
  Scheme W (Wu & Estabrook 2016): configural -> thresholds -> loadings -> intercepts -> residuals, higher levels only if the previous one is supported.
  Scheme M (Millsap & Tein 2004; Liu et al. 2017): baseline -> loadings -> thresholds -> unique variances. **The two schemes are not nested and are never merged into
  one ladder**: threshold and scalar verdicts come from W, the unconditional metric (loading) verdict from M. "Supported" means delta CFI >= -0.010 and delta RMSEA <= 0.015 (configurable).
  Scalar invariance is interpretable only if thresholds are supported; otherwise the status is `not_tested`, never "not supported". Also extracted: CR, AVE, latent correlation, Tucker's phi of CFA loadings.
* **Order sensitivity check** (`structural.order_check_permutations`): refits the CFA of every non-clean scale (configural model not converged/inadmissible, or threshold step not testable) on random row orders and flags scales whose configural admissibility is not identical across fits.
* **ESEM** (`r/esem_fit.R`): `lavaan::efa()` with target rotation (the theoretical simple structure), ordinal WLSMV, separately for Human and Twin, same target and same
  factors on both sides. **Admissible** = converged, no negative residual variance, no |factor correlation| > 1.0001.
* **Tucker's phi**: Twin loadings are orthogonally Procrustes-rotated onto Human loadings, phi per factor; Hungarian matching is *verified* against the theoretical order.
  Bootstrap (`r/tucker_bootstrap.R`): every replicate re-fits both ESEMs; inadmissible replicates are *not* replaced; a CI is "stable" if >= 80% of replicates are valid.
  The CI of the instrument mean phi uses the per-replicate mean over factors.
* **Admissibility gate**: Tucker's phi is *interpreted* (`tucker_phi_primary`) only if BOTH ESEM solutions of the instrument are admissible; all point estimates are kept as
  `tucker_phi_exploratory` and never used for typology.
* **Structural status** of a scale: see the typology section of the README.
* **Validation trap** (instrument level, `structural_instrument_master`): high phi (mean >= .95, min >= .90) from admissible, bootstrap-stable ESEM together with scalar invariance
  not established. Not assigned when configural CFA is unsupported or ESEM is inadmissible.

## Level 4 - item (`item.py`)
* Items of non-redundant scales, response range = **theoretical** range from the registry. Agreement/error (exact and adjacent agreement, MAE/RMSE in units of the range, Spearman rho,
  CCC, quadratic weighted kappa), bias (Wilcoxon signed-rank, BH), compression (SD ratio, normalised entropy, floor/ceiling, category coverage, calibration slope), person-level bootstrap CIs.
* **Variance decomposition** (ANOVA-style crossed): `Y_ij = mu + P_i + I_j + (PI)_ij` on items rescaled to [0, 1]; with one observation per cell the interaction cannot be separated from
  measurement error and is called *person x item / residual*. Human, twin and error (twin - human); paired bootstrap contrast of the percentages; instrument-weighted contrast.
* **Loading discrepancy**: `|delta loading|` on the target factor and the length of the difference vector of aligned loadings; associations with item fidelity only for instruments with admissible paired ESEM.
* **Reverse-keyed vs direct** items (cluster bootstrap over instruments).

## Level 5 - synthesis (`synthesis.py`)
* Master table with one row per scale (scale-level, structural, item aggregates by scale: medians).
* Spearman associations in pre-defined families: **A** within the scale level, **B** scale x item, **C2** exploratory Tucker (all point estimates), **C1** gated Tucker (descriptive).
  Cluster bootstrap over instruments; bootstrap p-value `2 (b+1)/(n+1)` (never exactly 0), BH within family; with fewer than `min_clusters_for_cluster_bootstrap` instruments CI and p are `NaN`.
* Typology and results matrix (README). Different units are never averaged into one index.

## Harmonised vs raw item data
In the example data 27 of 75 items have merged categories. Medians of item Spearman rho and CCC are identical in both versions; the largest per-item differences are 0.029 (exact agreement), 0.021 (normalised MAE), 0.015 (Spearman rho), 0.025 (CCC), 0.084 (SD ratio) and 0.107 (entropy difference). The median SD ratio is 0.825 on raw and 0.851 on harmonised categories.
