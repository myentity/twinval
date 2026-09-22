REPORTING-REVISED paired/dependent ordinal CFA.
Conventional status flags use DeltaCFI >= -0.01 AND DeltaRMSEA <= 0.015. Raw deltas and LRT p-values are also saved; DG-DMI should be run separately if required.
PRIMARY THRESHOLD RESULT: Wu & Estabrook (2016), theta parameterization, W0 configural -> W1 thresholds.
Wu W1->W2 loadings are CONDITIONAL on threshold invariance and are diagnostic only; do NOT report their rate as the unconditional metric/loading-invariance rate.
PRIMARY LOADING RESULT: Millsap & Tein (2004)/Liu et al. (2017), theta + marker identification, M0 baseline -> M1 loadings.
Millsap/Liu thresholds are tested only after loading invariance is supported; unique/residual variances only after loading+threshold invariance is supported.
The script intentionally skips more restrictive models when the preceding invariance level is rejected, inadmissible, nonconvergent, or not testable.
Millsap/Liu currently uses it1 as the marker in both Human and Twin factors.
Category labels must match exactly between Human and Twin for each repeated item; mismatches are stopped and reported.
Do NOT merge Wu-only and Millsap-only invariance decisions into one hybrid percentage. Report Wu threshold invariance and Millsap/Liu loading invariance separately.
invariance_reporting_table.csv is the manuscript-oriented table.
invariance_summary_counts.csv gives transparent denominators: only supported/rejected steps count as interpretable.
CR/AVE/Tucker/Fh-Ft are extracted separately from each scheme's least-restricted model; cross_scheme_common_admissible.csv can be used to assess robustness of these latent metrics across identification schemes.
Do not rerun DG-DMI on old models; if needed, rerun it using the corrected fitted-model specification.
