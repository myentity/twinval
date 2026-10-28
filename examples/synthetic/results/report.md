# twinval run report

* respondents (complete pairs): **500**
* scales loaded: **29** (28 non-redundant)
* instruments: **Synthetic - Faithful twin, Synthetic - Miscalibrated twin, Synthetic - Weakly individuating twin, Synthetic - Unequal twin loadings, Synthetic - Repeated twin answers (congruent), Synthetic - Mixed instrument with a total score, Synthetic - Short scales, Synthetic - Single-scale instrument, Synthetic - Repeated twin answers (discrepant), Synthetic - Repeated twin answers (inadmissible)**

## Data versions

* profile and scale level: **raw** input answers (scale scores are computed from the original items)
* structural models (CFA, ESEM): **harmonised** categories (`00_audit/harmonized_data/`, log in `harmonization_log.csv`)
* item level: **harmonized** (`item.data_version`; `harmonized` = setting of the original study)

## Multilevel results matrix

Different units of analysis are shown side by side and are never averaged into one index.

| level | unit | component | metric | value | 95% CI | detail |
|---|---|---|---|---|---|---|
| profile | respondent | individuation | Own-twin identification AUC (distinctive) | 0.982 | [0.977, 0.986] |  |
| profile | respondent | individuation | Top-1 identification (%) | 43.200 | [38.800, 47.400] |  |
| profile | respondent | calibration/agreement | Profile CCC, own (full) | 0.511 | [0.499, 0.524] | same specification as the identification model (non-redundant scales, equal instrument weights) |
| profile | respondent | calibration/agreement | Profile CCC, own (distinctive / bias-corrected) | 0.527 | [0.515, 0.540] | same specification as the identification model (non-redundant scales, equal instrument weights) |
| scale | scale | individuation | Pearson r, Fisher-z mean | 0.583 |  | raw arithmetic mean = 0.561; median = 0.633 [IQR 0.545, 0.689] |
| scale | scale | individuation | Scale-level AUC, mean | 0.639 |  | median = 0.652 |
| scale | scale | individuation | Disattenuated r (alpha-based, k >= 3), median | 0.827 |  | IQR [0.698, 0.872] |
| scale | scale | calibration/agreement | CCC, mean | 0.533 |  | median = 0.589 |
| scale | scale | calibration/agreement | C_b, mean | 0.956 |  | median = 0.996 |
| scale | scale | calibration/agreement | Delta CCC after LS recalibration, mean | 0.026 |  | median = -0.000 |
| structure | factor | structural fidelity | Tucker phi, GATED (both ESEM admissible), median | 0.986 |  | descriptive; bootstrap-stable only where >= 80% replicates are admissible |
| structure | factor | structural fidelity | Tucker phi, EXPLORATORY (incl. inadmissible ESEM), median | 0.986 |  | NOT used for typology or main conclusions |
| structure | scale | structural fidelity | Metric invariance (Wu chain) | 14.000 |  | yes: 14, no: 6, not_tested: 9 |
| structure | scale | structural fidelity | Scalar invariance (Wu chain) | 14.000 |  | yes: 14, no: 0, not_tested: 15 |
| item | item | individuation | Item-level Spearman rho, median | 0.310 |  | mean = 0.298; IQR [0.211, 0.396] |
| item | item | calibration/agreement | Item-level CCC, median | 0.309 |  | mean = 0.297; IQR [0.212, 0.394] |
| item | item | calibration/agreement | Item-level SD ratio, median | 0.994 |  | mean = 0.977; IQR [0.958, 1.019] |
| item | instrument | mechanism | Delta %var item (twin - human), pp | 0.669 | [0.177, 1.152] | ANOVA-style decomposition, paired bootstrap contrast |

## Scale typology (rule-based)

| type | n | % | example scales |
|---|---|---|---|
| A. Score fidelity with measurement-invariance support | 9 | 31.0 | Alpha; Beta; Gamma |
| B. Individuation without adequate calibration | 4 | 13.8 | Shifted; Compressed; Lowered |
| C. Limited individuation with structural evidence | 5 | 17.2 | Generic1; Generic2; Generic3 |
| D. Score fidelity with structural discrepancy | 1 | 3.4 | DupD |
| E. Score fidelity with limited structural evidence | 6 | 20.7 | Uneven1; Uneven2; Duplicated |
| F. Unclassified residual | 1 | 3.4 | DupF |
| Structural evidence not established (baseline CFA misfit) | 1 | 3.4 | Total |
| Structural model not estimable (nonconverged, not identified or inadmissible) | 2 | 6.9 | Solo; Duplicated |

## Associations between levels (Spearman)

| family | pair | rho | 95% CI (cluster bootstrap) | p (naive) | q (BH) | n | clusters |
|---|---|---|---|---|---|---|---|
| A | AUC vs Pearson r | 0.878 | [0.586, 0.989] | 0.0000 | 0.0012 | 29 | 10 |
| A | AUC vs CCC | 0.970 | [0.863, 0.991] | 0.0000 | 0.0012 | 29 | 10 |
| A | Pearson r vs CCC | 0.854 | [0.593, 1.000] | 0.0000 | 0.0012 | 29 | 10 |
| A | Disattenuated r vs raw r | 0.857 | [0.730, 0.965] | 0.0000 | 0.0012 | 29 | 10 |
| A | C_b vs change of CCC after OLS recalibration | -0.394 | [-0.774, 0.097] | 0.0344 | 0.1209 | 29 | 10 |
| B | Median item CCC vs scale CCC | 0.821 | [0.448, 0.941] | 0.0000 | 0.0165 | 28 | 10 |
| B | Median item rho vs scale r | 0.747 | [0.361, 0.917] | 0.0000 | 0.0165 | 28 | 10 |
| B | Median normalised item MAE vs scale CCC | -0.571 | [-0.857, -0.045] | 0.0015 | 0.0962 | 28 | 10 |
| B | Median |item bias| vs C_b | -0.558 | [-0.848, 0.012] | 0.0021 | 0.1297 | 28 | 10 |
| B | Median item SD ratio vs C_b | 0.522 | [0.077, 0.799] | 0.0044 | 0.0550 | 28 | 10 |
| B | Share of compressed items vs scale AUC | 0.167 | [-0.266, 0.560] | 0.3958 | 0.4276 | 28 | 10 |
| B | Share of compressed items vs scale r | 0.297 | [-0.107, 0.607] | 0.1249 | 0.2565 | 28 | 10 |
| B | Share of compressed items vs scale CCC | 0.097 | [-0.389, 0.503] | 0.6219 | 0.5977 | 28 | 10 |
| B | Median entropy difference vs scale AUC | -0.180 | [-0.562, 0.191] | 0.3607 | 0.4276 | 28 | 10 |
| B | Median entropy difference vs scale r | -0.196 | [-0.548, 0.204] | 0.3176 | 0.4276 | 28 | 10 |
| B | Median entropy difference vs scale CCC | -0.105 | [-0.468, 0.350] | 0.5946 | 0.5977 | 28 | 10 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs scale r | -0.212 | [-0.667, 0.284] | 0.2888 | 0.7156 | 27 | 9 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs scale AUC | -0.136 | [-0.662, 0.351] | 0.4983 | 0.7156 | 27 | 9 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs scale CCC | -0.132 | [-0.652, 0.399] | 0.5101 | 0.7156 | 27 | 9 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs item loading discrepancy | -0.738 | [-0.884, -0.558] | 0.0000 | 0.0050 | 27 | 9 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs item CCC | -0.093 | [-0.518, 0.287] | 0.6431 | 0.7156 | 27 | 9 |

## Output folders

Each stage folder holds its own `figures/` subfolder (png + pdf) - figures are never collected in one flat folder.

* `00_audit` data checks and harmonisation log
* `01_profile` respondent-level results, `figures/` profile CCC histogram
* `02_scale` scale-level results, `figures/` r vs C_b scatter
* `03_structural` CFA / ESEM / Tucker / admissibility, `figures/` Tucker's phi forest plot
* `04_item` item-level results, `figures/` item fidelity map, variance decomposition bars
* `05_synthesis` cross-level master, associations, typology, `figures/` typology map, association heat map
* `run_manifest.json`, `config_used.yaml` reproducibility record