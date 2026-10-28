# twinval run report

* respondents (complete pairs): **484**
* scales loaded: **24** (8 non-redundant)
* instruments: **Big Five Inventory-2 (BFI-2), Spectrum of Psychological Health**

## Data versions

* profile and scale level: **raw** input answers (scale scores are computed from the original items)
* structural models (CFA, ESEM): **harmonised** categories (`00_audit/harmonized_data/`, log in `harmonization_log.csv`)
* item level: **harmonized** (`item.data_version`; `harmonized` = setting of the original study)

## Notes and warnings

* Only 2 instruments are loaded: cluster (per-instrument) bootstrap CIs, equal-instrument-weight statistics and leave-one-instrument-out analyses are degenerate or skipped; results are illustrative and not comparable with a full study.

## Multilevel results matrix

Different units of analysis are shown side by side and are never averaged into one index.

| level | unit | component | metric | value | 95% CI | detail |
|---|---|---|---|---|---|---|
| profile | respondent | individuation | Own-twin identification AUC (distinctive) | 0.813 | [0.795, 0.831] |  |
| profile | respondent | individuation | Top-1 identification (%) | 2.686 | [1.240, 4.132] |  |
| profile | respondent | calibration/agreement | Profile CCC, own (full) | 0.349 | [0.323, 0.376] | same specification as the identification model (non-redundant scales, equal instrument weights) |
| profile | respondent | calibration/agreement | Profile CCC, own (distinctive / bias-corrected) | 0.379 | [0.354, 0.405] | same specification as the identification model (non-redundant scales, equal instrument weights) |
| scale | scale | individuation | Pearson r, Fisher-z mean | 0.528 |  | raw arithmetic mean = 0.524; median = 0.536 [IQR 0.470, 0.581] |
| scale | scale | individuation | Scale-level AUC, mean | 0.619 |  | median = 0.618 |
| scale | scale | individuation | Disattenuated r (alpha-based, k >= 3), median | 0.623 |  | IQR [0.576, 0.683] |
| scale | scale | calibration/agreement | CCC, mean | 0.469 |  | median = 0.484 |
| scale | scale | calibration/agreement | C_b, mean | 0.896 |  | median = 0.941 |
| scale | scale | calibration/agreement | Delta CCC after LS recalibration, mean | 0.052 |  | median = +0.028 |
| structure | factor | structural fidelity | Tucker phi, GATED (both ESEM admissible), median | 0.939 |  | descriptive; bootstrap-stable only where >= 80% replicates are admissible |
| structure | factor | structural fidelity | Tucker phi, EXPLORATORY (incl. inadmissible ESEM), median | 0.939 |  | NOT used for typology or main conclusions |
| structure | scale | structural fidelity | Metric invariance (Wu chain) | 6.000 |  | yes: 6, no: 1, not_tested: 17 |
| structure | scale | structural fidelity | Scalar invariance (Wu chain) | 2.000 |  | yes: 2, no: 4, not_tested: 18 |
| item | item | individuation | Item-level Spearman rho, median | 0.373 |  | mean = 0.365; IQR [0.295, 0.421] |
| item | item | calibration/agreement | Item-level CCC, median | 0.337 |  | mean = 0.338; IQR [0.274, 0.405] |
| item | item | calibration/agreement | Item-level SD ratio, median | 0.851 |  | mean = 0.847; IQR [0.748, 0.934] |
| item | instrument | mechanism | Delta %var item (twin - human), pp | 16.931 | [15.338, 18.633] | ANOVA-style decomposition, paired bootstrap contrast |

## Scale typology (rule-based)

| type | n | % | example scales |
|---|---|---|---|
| A. Score fidelity with measurement-invariance support | 3 | 12.5 | Conscientiousness; Productiveness; Negative Emotionality |
| B. Individuation without adequate calibration | 3 | 12.5 | Agreeableness; Compassion; Psychological Well-Being |
| C. Limited individuation with structural evidence | 4 | 16.7 | Respectfulness; Responsibility; Emotional Volatility |
| D. Score fidelity with structural discrepancy | 1 | 4.2 | Emotional Well-Being |
| E. Score fidelity with limited structural evidence | 6 | 25.0 | Extraversion; Persistence; Anxiety |
| Structural evidence not established (baseline CFA misfit) | 1 | 4.2 | Overall Score |
| Structural model not estimable (nonconverged, not identified or inadmissible) | 6 | 25.0 | Sociability; Energy Level; Trust |

## Associations between levels (Spearman)

| family | pair | rho | 95% CI (cluster bootstrap) | p (naive) | q (BH) | n | clusters |
|---|---|---|---|---|---|---|---|
| A | AUC vs Pearson r | 0.841 | n/a | 0.0000 | n/a | 24 | 2 |
| A | AUC vs CCC | 0.890 | n/a | 0.0000 | n/a | 24 | 2 |
| A | Pearson r vs CCC | 0.790 | n/a | 0.0000 | n/a | 24 | 2 |
| A | Disattenuated r vs raw r | 0.876 | n/a | 0.0000 | n/a | 24 | 2 |
| A | C_b vs change of CCC after OLS recalibration | -0.790 | n/a | 0.0000 | n/a | 24 | 2 |
| B | Median item CCC vs scale CCC | 0.643 | n/a | 0.0856 | n/a | 8 | 2 |
| B | Median item rho vs scale r | 0.619 | n/a | 0.1017 | n/a | 8 | 2 |
| B | Median normalised item MAE vs scale CCC | 0.071 | n/a | 0.8665 | n/a | 8 | 2 |
| B | Median |item bias| vs C_b | -0.857 | n/a | 0.0065 | n/a | 8 | 2 |
| B | Median item SD ratio vs C_b | 0.571 | n/a | 0.1390 | n/a | 8 | 2 |
| B | Share of compressed items vs scale AUC | 0.062 | n/a | 0.8834 | n/a | 8 | 2 |
| B | Share of compressed items vs scale r | 0.171 | n/a | 0.6848 | n/a | 8 | 2 |
| B | Share of compressed items vs scale CCC | -0.296 | n/a | 0.4763 | n/a | 8 | 2 |
| B | Median entropy difference vs scale AUC | 0.214 | n/a | 0.6103 | n/a | 8 | 2 |
| B | Median entropy difference vs scale r | 0.095 | n/a | 0.8225 | n/a | 8 | 2 |
| B | Median entropy difference vs scale CCC | 0.429 | n/a | 0.2894 | n/a | 8 | 2 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs scale r | 0.476 | n/a | 0.2329 | n/a | 8 | 2 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs scale AUC | 0.357 | n/a | 0.3851 | n/a | 8 | 2 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs scale CCC | 0.595 | n/a | 0.1195 | n/a | 8 | 2 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs item loading discrepancy | -0.429 | n/a | 0.2894 | n/a | 8 | 2 |
| C2_exploratory | [EXPLORATORY] Tucker phi vs item CCC | 0.048 | n/a | 0.9108 | n/a | 8 | 2 |

## Output folders

Each stage folder holds its own `figures/` subfolder (png + pdf) - figures are never collected in one flat folder.

* `00_audit` data checks and harmonisation log
* `01_profile` respondent-level results, `figures/` profile CCC histogram
* `02_scale` scale-level results, `figures/` r vs C_b scatter
* `03_structural` CFA / ESEM / Tucker / admissibility, `figures/` Tucker's phi forest plot
* `04_item` item-level results, `figures/` item fidelity map, variance decomposition bars
* `05_synthesis` cross-level master, associations, typology, `figures/` typology map, association heat map
* `run_manifest.json`, `config_used.yaml` reproducibility record