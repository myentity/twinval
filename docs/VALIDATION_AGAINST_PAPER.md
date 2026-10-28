# Validation of the code against the published study

For the two bundled instruments (BFI-2 and the Mental Health Continuum, N = 484) the framework reproduces the values
of the original study. The comparison is part of the automatic test suite (`tests/test_reference_paper.py`; the reference values
are in `tests/reference/paper_*.csv`).

| Level | Result | Agreement |
|---|---|---|
| Scale | Pearson r, CCC, C_b per scale (24 scales) | to 1e-3 (published scores were rounded to 3 decimals) |
| Scale | scale-level AUC per scale | to 0.01 (rounded scores produce ties) |
| Scale | delta CCC after LS / OLS recalibration | to 5e-3 (different bootstrap streams do not enter the point estimates) |
| Scale | disattenuated r (alpha) | to 1e-3 |
| Scale | descriptive r x C_b class | identical |
| Item | agreement, MAE, rho, CCC, SD ratio, bias, entropy difference for 75 items (original categories) | to 1e-6 |
| Item | variance decomposition (human / twin / error) | to 1e-6 |
| Structure | ESEM admissibility, Tucker's phi for 8 factors | both ESEM admissible; phi to 1e-5 |
| Structure | CFA invariance decisions (thresholds, loadings, scalar) for 24 scales | identical |
| Synthesis | structural status and final type of all 24 scales | identical (types A-E all occur) |

## Other notes

Profile-level results need several instruments, so they are not part of the two-instrument tests. Note that a few scales of the original
study have small, nearly saturated CFA models that converge or fail depending only on the **row order** of the data (see *Caveats* in the
README and the option `structural.order_check_permutations`); for the example instruments this concerns only the BFI-2 facet
*Creative Imagination*, whose baseline is not estimable in every row order (nonconverged or not identified), so that its
structural status does not change.

## Not ported

Exploratory analyses of the earlier stage of the project that are not part of the article's results (dynamic-fit-index simulations,
"over-coherence" regressions, forward-regression reports) are not included.
