# Output tables

Column lists were generated from a run on `examples/data`. Every file is a CSV in the run folder; `*_per_person` tables contain the respondent `id`.

Conventions: `*_ci_lo/_ci_hi` are 95% percentile bootstrap limits; `q_BH` are Benjamini-Hochberg adjusted; `p_below_resolution` marks bootstrap p-values at the floor `2/(n+1)`; statuses `not_tested` mean *unknown*, never *not equivalent*.

## `00_audit/data_audit.csv`

which scale files were loaded, rows per file, number of complete respondents.

Columns: `scale_id`, `status`, `file`, `n_rows`, `n_complete_respondents`

## `00_audit/scale_score_qc.csv`

per-scale score ranges and zero-variance flags.

Columns: `scale_id`, `scale_name`, `instrument_id`, `n_items`, `scoring`, `redundant`, `n`, `min_human`, `max_human`, `min_twin`, `max_twin`, `zero_variance_human`, `zero_variance_twin`

## `00_audit/harmonization_log.csv`

one row per item whose categories were merged (before/after levels, merge map).

Columns: `scale_id`, `item_index`, `human_levels_before`, `twin_levels_before`, `n_cat_after`, `merge_map`

## `00_audit/harmonization_status.csv`

per scale: clean | harmonized | unusable.

Columns: `scale_id`, `category_status`

## `00_audit/harmonization_unusable.csv`

items that collapse to a single category (scale excluded from structural models).

Columns: `scale_id`, `item_index`, `human_levels`, `twin_levels`, `reason`

## `00_audit/harmonized_data/<scale_id>_paired_wide.csv`

the harmonised paired data, one file per scale, same format as the input (`ID, it1_h..itK_h, it1_t..itK_t`), rows sorted by `ID`. Items whose Human and Twin category sets differ have rare categories merged identically on both sides (merged categories are re-numbered 1..k'; the map is in `harmonization_log.csv`); all other items are unchanged. These are the data of the structural models and, with `item.data_version: harmonized` (default), of the item level. Scale and profile levels use the raw input files.

## `01_profile/identification_summary.csv`

Top-1/5/10, AUC, ranks, own vs cross distance, d_z with bootstrap CIs; rows `full` and `distinctive`.

Columns: `spec`, `n`, `top1_pct`, `top1_ci_lo`, `top1_ci_hi`, `top5_pct`, `top5_ci_lo`, `top5_ci_hi`, `top10_pct`, `top10_ci_lo`, `top10_ci_hi`, `auc`, `auc_ci_lo`, `auc_ci_hi`, `mean_rank`, `mean_rank_ci_lo`, `mean_rank_ci_hi`, `own_mean`, `own_mean_ci_lo`, `own_mean_ci_hi`, `cross_mean`, `cross_mean_ci_lo`, `cross_mean_ci_hi`, `median_rank`, `rank_iqr_lo`, `rank_iqr_hi`, `diff_mean`, `diff_ci_lo`, `diff_ci_hi`, `dz`, `dz_ci_lo`, `dz_ci_hi`, `paired_t`, `paired_t_p`, `n_instruments`, `n_scales`

## `01_profile/identification_permutation.csv`

label-permutation tests of Top-1, Top-5 and AUC.

Columns: `spec`, `n_permutations`, `obs_top1`, `obs_top5`, `obs_auc`, `null_top1_mean`, `null_auc_mean`, `p_top1`, `p_top5`, `p_auc`

## `01_profile/identification_per_person.csv`

per respondent: own/cross distance, rank of the own twin, AUC, percentile (full and distinctive).

Columns: `id`, `own_distance_full`, `cross_distance_full`, `rank_full`, `auc_full`, `percentile_full`, `own_distance_distinctive`, `cross_distance_distinctive`, `rank_distinctive`, `auc_distinctive`, `percentile_distinctive`

## `01_profile/identification_robustness.csv`

identification under alternative specifications.

Columns: `spec`, `n_scales`, `n_instruments`, `metric`, `centred`, `top1_pct`, `top5_pct`, `auc`, `dz`

## `01_profile/profile_calibration_summary.csv`

profile CCC (own/cross, full/distinctive), MAE/RMSE (full/bias-corrected), mean absolute bias, RMS bias.

Columns: `metric`, `statistic`, `value`, `ci_lo`, `ci_hi`, `median`, `iqr_lo`, `iqr_hi`, `n_scales`, `n_instruments`, `weighting`, `standardization`

## `01_profile/profile_calibration_per_person.csv`

per respondent: own/cross profile CCC, MAE, RMSE.

Columns: `id`, `own_ccc_full`, `cross_ccc_full_mean`, `own_ccc_distinctive`, `cross_ccc_distinctive_mean`, `mae_full`, `rmse_full`, `mae_bias_corrected`, `rmse_bias_corrected`

## `01_profile/profile_bias_by_scale.csv`

mean twin-minus-human bias per non-redundant scale in human SD units with bootstrap CI.

Columns: `scale_id`, `instrument_id`, `weight`, `bias_twin_minus_human_sd`, `bias_ci_lo`, `bias_ci_hi`, `scale_name`

## `01_profile/profile_mahalanobis_mean_profile.csv`

Mahalanobis own vs cross, mean-profile similarity.

Columns: `n_scales`, `ledoit_wolf_shrinkage`, `mahalanobis_own_mean`, `mahalanobis_own_ci_lo`, `mahalanobis_own_ci_hi`, `mahalanobis_cross_mean`, `mahalanobis_dz`, `mahalanobis_auc`, `mean_profile_procrustes_disparity`, `mean_profile_euclidean`, `mean_profile_pearson_r`, `mean_profile_cosine`, `note`

## `02_scale/scale_core_table.csv`

one row per scale: r, rho, AUC, CCC, C_b, u, v with CIs, permutation and BH results, descriptive r x C_b class.

Columns: `scale_id`, `instrument_id`, `n`, `n_items`, `is_nonredundant`, `pearson_r`, `pearson_r_ci_lo`, `pearson_r_ci_hi`, `pearson_p`, `spearman_rho`, `spearman_p`, `auc`, `auc_ci_lo`, `auc_ci_hi`, `auc_perm_p`, `auc_perm_null_mean`, `auc_perm_null_sd`, `d_own_mean`, `d_cross_mean`, `auc_cohens_dz`, `tie_share`, `ccc`, `ccc_ci_lo`, `ccc_ci_hi`, `cb`, `cb_ci_lo`, `cb_ci_hi`, `u`, `v`, `n_valid_bootstrap`, `pearson_q_bh`, `pearson_q_bh_sig`, `auc_perm_q_bh`, `auc_perm_q_bh_sig`, `typology_r_cb`

## `02_scale/scale_recalibration_per_repeat.csv`

CV recalibration metrics for every scale and repeat.

Columns: `scale_id`, `repeat`, `ccc_raw`, `cb_raw`, `r_raw`, `mae_raw`, `rmse_raw`, `mae_raw_std`, `rmse_raw_std`, `ccc_ls`, `cb_ls`, `r_ls`, `mae_ls`, `rmse_ls`, `mae_ls_std`, `rmse_ls_std`, `ccc_ols`, `cb_ols`, `r_ols`, `mae_ols`, `rmse_ols`, `mae_ols_std`, `rmse_ols_std`, `delta_ccc_ls`, `delta_ccc_ols`

## `02_scale/scale_recalibration_summary.csv`

delta CCC after LS/OLS recalibration (mean over repeats, nested-bootstrap CI), out-of-fold C_b, MAE/RMSE in human SD.

Columns: `scale_id`, `delta_ccc_ls`, `delta_ccc_ls_sd_repeats`, `delta_ccc_ls_n_positive_repeats`, `delta_ccc_ls_boot_mean`, `delta_ccc_ls_boot_ci_lo`, `delta_ccc_ls_boot_ci_hi`, `delta_ccc_ls_ci_excludes_zero`, `ccc_ls`, `cb_ls`, `mae_ls_std`, `rmse_ls_std`, `delta_ccc_ols`, `delta_ccc_ols_sd_repeats`, `delta_ccc_ols_n_positive_repeats`, `delta_ccc_ols_boot_mean`, `delta_ccc_ols_boot_ci_lo`, `delta_ccc_ols_boot_ci_hi`, `delta_ccc_ols_ci_excludes_zero`, `ccc_ols`, `cb_ols`, `mae_ols_std`, `rmse_ols_std`, `cb_raw`, `mae_raw_std`, `rmse_raw_std`

## `02_scale/scale_disattenuation.csv`

reliabilities (alpha / Spearman-Brown / CR) and disattenuated r with CI.

Columns: `scale_id`, `n_items`, `pearson_r`, `reliability_type`, `rel_human`, `rel_twin`, `reliability_ceiling`, `r_dis`, `r_dis_ci_lo`, `r_dis_ci_hi`, `undefined_reason`, `r_dis_alpha`, `r_dis_sb`, `flag_abs_r_dis_gt_1`, `cr_human`, `cr_twin`, `r_dis_cr`

## `02_scale/scale_level_master.csv`

merge of the three scale-level tables (one row per scale).

Columns: `scale_id`, `scale_name`, `instrument_name`, `instrument_id`, `n`, `n_items`, `is_nonredundant`, `pearson_r`, `pearson_r_ci_lo`, `pearson_r_ci_hi`, `pearson_p`, `spearman_rho`, `spearman_p`, `auc`, `auc_ci_lo`, `auc_ci_hi`, `auc_perm_p`, `auc_perm_null_mean`, `auc_perm_null_sd`, `d_own_mean`, `d_cross_mean`, `auc_cohens_dz`, `tie_share`, `ccc`, `ccc_ci_lo`, `ccc_ci_hi`, `cb`, `cb_ci_lo`, `cb_ci_hi`, `u`, `v`, `n_valid_bootstrap`, `pearson_q_bh`, `pearson_q_bh_sig`, `auc_perm_q_bh`, `auc_perm_q_bh_sig`, `typology_r_cb`, `delta_ccc_ls`, `delta_ccc_ls_sd_repeats`, `delta_ccc_ls_n_positive_repeats`, `delta_ccc_ls_boot_mean`, `delta_ccc_ls_boot_ci_lo`, `delta_ccc_ls_boot_ci_hi`, `delta_ccc_ls_ci_excludes_zero`, `ccc_ls`, `cb_ls`, `mae_ls_std`, `rmse_ls_std`, `delta_ccc_ols`, `delta_ccc_ols_sd_repeats`, `delta_ccc_ols_n_positive_repeats`, `delta_ccc_ols_boot_mean`, `delta_ccc_ols_boot_ci_lo`, `delta_ccc_ols_boot_ci_hi`, `delta_ccc_ols_ci_excludes_zero`, `ccc_ols`, `cb_ols`, `mae_ols_std`, `rmse_ols_std`, `cb_raw`, `mae_raw_std`, `rmse_raw_std`, `reliability_type`, `rel_human`, `rel_twin`, `reliability_ceiling`, `r_dis`, `r_dis_alpha`, `r_dis_sb`, `r_dis_ci_lo`, `r_dis_ci_hi`, `flag_abs_r_dis_gt_1`, `cr_human`, `cr_twin`, `r_dis_cr`

## `02_scale/scale_level_summary.csv`

mean, median, IQR, Fisher-z mean, instrument-weighted and non-redundant means of the main scale metrics.

Columns: `metric`, `n_scales`, `mean`, `median`, `iqr_lo`, `iqr_hi`, `min`, `max`, `fisher_z_mean`, `instrument_equal_weight_mean`, `nonredundant_mean`

## `03_structural/invariance_results.csv`

raw output of the paired CFA script: fit, deltas, statuses of every invariance step.

Columns: `scale_id`, `k_items`, `n_pairs_manifest`, `n_rows_data`, `category_mismatch_items`, `wu_configural_id`, `millsap_marker_item`, `wu_configural_converged`, `wu_configural_cfi`, `wu_configural_rmsea`, `wu_configural_srmr`, `wu_configural_df`, `wu_configural_heywood`, `wu_thresholds_converged`, `wu_thresholds_cfi`, `wu_thresholds_rmsea`, `wu_thresholds_srmr`, `wu_thresholds_df`, `wu_thresholds_heywood`, `wu_conditional_loading_converged`, `wu_conditional_loading_cfi`, `wu_conditional_loading_rmsea`, `wu_conditional_loading_srmr`, `wu_conditional_loading_df`, `wu_conditional_loading_heywood`, `wu_conditional_scalar_converged`, `wu_conditional_scalar_cfi`, `wu_conditional_scalar_rmsea`, `wu_conditional_scalar_srmr`, `wu_conditional_scalar_df`, `wu_conditional_scalar_heywood`, `wu_conditional_strict_converged`, `wu_conditional_strict_cfi`, `wu_conditional_strict_rmsea`, `wu_conditional_strict_srmr`, `wu_conditional_strict_df`, `wu_conditional_strict_heywood`, `ms_baseline_converged`, `ms_baseline_cfi`, `ms_baseline_rmsea`, `ms_baseline_srmr`, `ms_baseline_df`, `ms_baseline_heywood`, `ms_loading_converged`, `ms_loading_cfi`, `ms_loading_rmsea`, `ms_loading_srmr`, `ms_loading_df`, `ms_loading_heywood`, `ms_threshold_after_loading_converged`, `ms_threshold_after_loading_cfi`, `ms_threshold_after_loading_rmsea`, `ms_threshold_after_loading_srmr`, `ms_threshold_after_loading_df`, `ms_threshold_after_loading_heywood`, `ms_unique_after_threshold_converged`, `ms_unique_after_threshold_cfi`, `ms_unique_after_threshold_rmsea`, `ms_unique_after_threshold_srmr`, `ms_unique_after_threshold_df`, `ms_unique_after_threshold_heywood`, `wu_threshold_status`, `wu_threshold_delta_cfi`, `wu_threshold_delta_rmsea`, `wu_threshold_delta_srmr`, `wu_threshold_lrt_p`, `wu_threshold_df_diff`, `wu_conditional_loading_status`, `wu_conditional_loading_delta_cfi`, `wu_conditional_loading_delta_rmsea`, `wu_conditional_loading_delta_srmr`, `wu_conditional_loading_lrt_p`, `wu_conditional_loading_df_diff`, `wu_conditional_scalar_status`, `wu_conditional_scalar_delta_cfi`, `wu_conditional_scalar_delta_rmsea`, `wu_conditional_scalar_delta_srmr`, `wu_conditional_scalar_lrt_p`, `wu_conditional_scalar_df_diff`, `wu_conditional_strict_status`, `wu_conditional_strict_delta_cfi`, `wu_conditional_strict_delta_rmsea`, `wu_conditional_strict_delta_srmr`, `wu_conditional_strict_lrt_p`, `wu_conditional_strict_df_diff`, `wu_latent_cor_raw`, `wu_latent_cor_oriented`, `wu_tucker_phi`, `wu_cr_h`, `wu_cr_t`, `wu_ave_h`, `wu_ave_t`, `wu_negative_loadings_h`, `wu_negative_loadings_t`, `ms_loading_status`, `ms_loading_delta_cfi`, `ms_loading_delta_rmsea`, `ms_loading_delta_srmr`, `ms_loading_lrt_p`, `ms_loading_df_diff`, `ms_threshold_after_loading_status`, `ms_threshold_after_loading_delta_cfi`, `ms_threshold_after_loading_delta_rmsea`, `ms_threshold_after_loading_delta_srmr`, `ms_threshold_after_loading_lrt_p`, `ms_threshold_after_loading_df_diff`, `ms_unique_after_threshold_status`, `ms_unique_after_threshold_delta_cfi`, `ms_unique_after_threshold_delta_rmsea`, `ms_unique_after_threshold_delta_srmr`, `ms_unique_after_threshold_lrt_p`, `ms_unique_after_threshold_df_diff`, `ms_latent_cor_raw`, `ms_latent_cor_oriented`, `ms_tucker_phi`, `ms_cr_h`, `ms_cr_t`, `ms_ave_h`, `ms_ave_t`, `manifest_instrument`, `manifest_scale_name`, `manifest_method_count`, `manifest_category_status`

Model diagnostics, for every model prefix (`wu_configural`, `wu_thresholds`, `wu_conditional_loading`, `wu_conditional_scalar`, `wu_conditional_strict`, `ms_baseline`, `ms_loading`, `ms_threshold_after_loading`, `ms_unique_after_threshold`): `<prefix>_chisq`, `_npar`, `_tli`, `_rmsea_ci_lo`, `_rmsea_ci_hi`, `_post_check`, `_vcov_min_eigen`, `_theta_min_eigen`, `_psi_min_eigen`, `_max_sample_cor_h`, `_max_sample_cor_t`, `_max_std_loading`, `_min_std_resid_var`, `_max_abs_resid_cor`, `_abs_factor_cor`.

## `03_structural/invariance_decisions.csv`

statuses and fit deltas of the invariance tests per scale.

Columns: `scale_id`, `wu_threshold_status`, `wu_threshold_delta_cfi`, `wu_threshold_delta_rmsea`, `wu_threshold_delta_srmr`, `wu_threshold_lrt_p`, `wu_conditional_loading_status`, `wu_conditional_loading_delta_cfi`, `wu_conditional_loading_delta_rmsea`, `wu_conditional_scalar_status`, `wu_conditional_scalar_delta_cfi`, `wu_conditional_scalar_delta_rmsea`, `wu_conditional_strict_status`, `ms_loading_status`, `ms_loading_delta_cfi`, `ms_loading_delta_rmsea`, `ms_loading_delta_srmr`, `ms_loading_lrt_p`, `ms_threshold_after_loading_status`, `ms_threshold_after_loading_delta_cfi`, `ms_unique_after_threshold_status`

## `03_structural/paired_cfa_invariance_fit.csv`

fit indices of every model of every scale (long format).

Columns: `scale_id`, `scheme`, `step`, `converged`, `cfi`, `rmsea`, `srmr`, `df`, `heywood`

## `03_structural/cfa_scale_status.csv`

per scale: estimability and fit of the baseline, every step decision of the primary (Wu) sequence, the CFA conclusion, the Millsap/Liu sensitivity result and the CFA flags. `metric_invariance` / `scalar_invariance` follow the CFA conclusion (a step supported under a misfitting baseline is only a flag).

Columns: `scale_id`, `primary_branch`, `baseline_status`, `baseline_fit`, `threshold_decision`, `metric_decision`, `scalar_decision`, `cfa_conclusion`, `cfa_no_conclusion_reason`, `configural_cfa_admissible`, `metric_invariance`, `scalar_invariance`, `metric_test_available`, `scalar_test_available`, `sensitivity_millsap_baseline_status`, `sensitivity_millsap_metric_decision`, `sensitivity_millsap_scalar_decision`, `sensitivity_millsap_cfa_conclusion`, `cfa_flags`

## `03_structural/cfa_step_diagnostics.csv`

the basis of every CFA decision, one row per scale and fitted model of both sequences (Wu/Estabrook primary, Millsap/Liu sensitivity): estimability status and near-boundary flags with the numerical diagnostics, absolute fit, and for comparisons the change in df / CFI / RMSEA / SRMR with the step decision.

Columns: `scale_id`, `sequence`, `step`, `model`, `compared_with`, `status`, `near_boundary_flags`, `df`, `chisq`, `cfi`, `tli`, `rmsea`, `rmsea_ci_lo`, `rmsea_ci_hi`, `srmr`, `converged`, `heywood`, `post_check`, `vcov_min_eigen`, `theta_min_eigen`, `psi_min_eigen`, `max_std_loading`, `min_std_resid_var`, `max_abs_resid_cor`, `abs_factor_cor`, `decision`, `delta_df`, `delta_cfi`, `delta_rmsea`, `delta_srmr`

## `03_structural/cfa_status_by_instrument.csv`

instrument-level CFA synthesis. Every verdict follows the CFA conclusion of the scales - `configural_status` is `supported` if every scale has a CFA conclusion, otherwise the common reason (`misfit`, `not_estimable`, `fit_not_testable`) or `mixed`; `metric_status` / `scalar_status` are `supported` if every scale reaches that level, `not_supported` if a scale with a conclusion rejects the step, otherwise `not_established`; the counts include only scales with a CFA conclusion.

Columns: `instrument_id`, `instrument_name`, `n_scales`, `primary_branch`, `configural_status`, `n_scales_configural_ok`, `n_scales_cfa_misfit`, `n_scales_fit_not_testable`, `n_scales_not_estimable`, `metric_status`, `scalar_status`, `n_scales_metric_supported`, `n_scales_metric_rejected`, `n_scales_scalar_supported`, `n_scales_scalar_rejected`

## `03_structural/cfa_reliability.csv`

CR, AVE, Tucker's phi and latent correlation from the configural CFA.

Columns: `scale_id`, `cr_h`, `cr_t`, `ave_h`, `ave_t`, `wu_tucker_phi`, `wu_latent_cor_oriented`, `wu_negative_loadings_h`, `wu_negative_loadings_t`

## `03_structural/esem_model_diagnostics.csv`

status of every ESEM fit (ok, inadmissible_heywood, nonconverged, ...).

Columns: `instrument_id`, `instrument_name`, `side`, `status`, `note`, `elapsed_sec`

## `03_structural/esem_admissibility.csv`

per instrument: human / twin / paired ESEM admissible.

Columns: `instrument_id`, `human_esem_admissible`, `twin_esem_admissible`, `paired_esem_admissible`

## `03_structural/tucker_by_factor.csv`

Tucker's phi per factor (Procrustes-aligned), bootstrap CI and share of valid replicates.

Columns: `instrument_id`, `instrument_name`, `factor`, `matched_twin_factor`, `phi_full_column`, `phi_target_items_only`, `status_human`, `status_twin`, `category`, `n_valid_bootstrap`, `n_total_bootstrap`, `frac_valid_bootstrap`, `phi_boot_mean`, `phi_ci_lo`, `phi_ci_hi`, `ci_stable`

## `03_structural/tucker_summary.csv`

instrument-level phi mean/min and matching checks.

Columns: `instrument_id`, `instrument_name`, `n_factors`, `status_human`, `status_twin`, `matched_theoretical_order`, `procrustes_is_reflection`, `phi_mean`, `phi_median`, `phi_min`, `n_phi_ge_095`, `n_phi_090_095`, `n_phi_085_090`, `n_phi_lt_085`

## `03_structural/tucker_mean_phi_ci.csv`

bootstrap CI of the instrument mean phi.

Columns: `instrument_id`, `n_valid_reps_for_mean_phi`, `phi_mean_ci_lo`, `phi_mean_ci_hi`

## `03_structural/esem_aligned_loadings.csv`

Human and Procrustes-aligned Twin loadings (long).

Columns: `instrument_id`, `instrument_name`, `factor`, `item`, `loading_human`, `loading_twin_aligned`, `is_target`

## `03_structural/structural_scale_table.csv`

one row per scale: CFA conclusion and baseline statuses, ESEM admissibility, primary/exploratory phi, structural_evidence_level (S) and `structural_flags`.

Columns: `scale_id`, `instrument_id`, `configural_cfa_admissible`, `metric_invariance`, `scalar_invariance`, `metric_test_available`, `scalar_test_available`, `esem_applicable`, `human_esem_admissible`, `twin_esem_admissible`, `paired_esem_admissible`, `tucker_point_available`, `tucker_interpretable`, `tucker_bootstrap_stable`, `tucker_phi_exploratory`, `tucker_phi_primary`, `tucker_phi_ci_lo`, `tucker_phi_ci_hi`, `tucker_frac_valid_bootstrap`, `structural_evidence_level`, `structural_flags`, `baseline_status`, `baseline_fit`, `threshold_decision`, `metric_decision`, `scalar_decision`, `cfa_conclusion`, `cfa_no_conclusion_reason`, `cfa_order_sensitive`

## `03_structural/structural_instrument_master.csv`

instrument-level structural synthesis and validation-trap flag (assessed only when `configural_status` is `supported`).

Columns: `instrument_id`, `instrument_name`, `esem_applicable`, `paired_esem_admissible`, `phi_mean`, `phi_min`, `phi_ci_stable`, `configural_status`, `metric_status`, `scalar_status`, `validation_trap_strict`, `n_scales_metric_supported`, `n_scales_metric_rejected`, `n_valid_reps_for_mean_phi`, `phi_mean_ci_lo`, `phi_mean_ci_hi`

## `03_structural/cfa_order_sensitivity.csv`

refits of every candidate CFA model (baseline not admissible or first step not evaluable) on random row orders (`structural.order_check_permutations`); `order_sensitive` flags scales whose baseline admissibility or CFA conclusion depends on row order (copied as `cfa_order_sensitive` to `structural_scale_table.csv`).

Columns: `scale_id`, `n_fits`, `share_configural_admissible`, `n_distinct_verdicts`, `baseline_statuses`, `order_sensitive`

## `03_structural/esem_factor_correlations.csv`

latent factor correlations of the Human and Twin ESEM solutions.

Columns: `instrument_id`, `instrument_name`, `side`, `factor1`, `factor2`, `correlation`, `se`, `pvalue`

## `03_structural/esem_fit_by_instrument.csv`

fit indices of the ESEM models (Human and Twin) per instrument.

Columns: `instrument_id`, `instrument_name`, `side`, `status`, `n_factors`, `n_items`, `chisq`, `df`, `pvalue`, `cfi`, `tli`, `rmsea`, `rmsea_ci_lo`, `rmsea_ci_hi`, `srmr`

## `03_structural/esem_loadings_long.csv`

ESEM standardised loadings, long format (item x factor x side).

Columns: `instrument_id`, `instrument_name`, `side`, `factor`, `item`, `loading`, `se`, `pvalue`, `is_target`

## `03_structural/esem_thresholds.csv`

item thresholds of the ordinal ESEM models.

Columns: `instrument_id`, `instrument_name`, `side`, `item`, `threshold`, `value`, `se`

## `03_structural/item_ordinal_suitability_qc.csv`

category counts and sparsity checks behind the ordinal (WLSMV) models.

Columns: `instrument_id`, `item`, `side`, `n_categories_observed`, `min_category_freq`, `modal_category_share`, `empty_categories_within_range`, `near_zero_variance`

## `04_item/item_level_master.csv`

merge of agreement, bias/compression, audit and loading discrepancy (one row per item).

Columns: `item_id`, `n_pairs`, `exact_agreement`, `adjacent_agreement`, `mae_raw`, `mae_range`, `rmse_raw`, `rmse_range`, `spearman_rho`, `pearson_r`, `ccc`, `weighted_kappa_quadratic`, `exact_agreement_ci_lo`, `exact_agreement_ci_hi`, `mae_range_ci_lo`, `mae_range_ci_hi`, `spearman_rho_ci_lo`, `spearman_rho_ci_hi`, `ccc_ci_lo`, `ccc_ci_hi`, `instrument_id`, `instrument_name`, `scale_id`, `scale_name`, `reverse_keyed`, `item_min`, `item_max`, `n_categories`, `data_version`, `mean_human`, `mean_twin`, `bias_raw`, `bias_range`, `bias_p_wilcoxon`, `bias_q_bh`, `sd_human`, `sd_twin`, `sd_ratio`, `log_sd_ratio`, `var_ratio`, `entropy_human_norm`, `entropy_twin_norm`, `entropy_difference`, `n_categories_used_human`, `n_categories_used_twin`, `modal_category_share_twin`, `floor_share_twin`, `ceiling_share_twin`, `calibration_slope`, `calibration_intercept`, `floor_share_human`, `ceiling_share_human`, `extreme_share_human`, `extreme_share_twin`, `extreme_share_difference`, `category_coverage_ratio`, `bias_range_ci_lo`, `bias_range_ci_hi`, `sd_ratio_ci_lo`, `sd_ratio_ci_hi`, `entropy_difference_ci_lo`, `entropy_difference_ci_hi`, `twin_top2_categories_share`, `twin_near_degenerate_flag`, `abs_delta_loading_target`, `loading_vector_distance`, `n_factors`, `esem_admissible_both_sides`

## `04_item/item_agreement_error.csv`

exact/adjacent agreement, MAE, RMSE, rho, CCC, kappa per item with CIs.

Columns: `item_id`, `n_pairs`, `exact_agreement`, `adjacent_agreement`, `mae_raw`, `mae_range`, `rmse_raw`, `rmse_range`, `spearman_rho`, `pearson_r`, `ccc`, `weighted_kappa_quadratic`, `exact_agreement_ci_lo`, `exact_agreement_ci_hi`, `mae_range_ci_lo`, `mae_range_ci_hi`, `spearman_rho_ci_lo`, `spearman_rho_ci_hi`, `ccc_ci_lo`, `ccc_ci_hi`, `instrument_id`, `instrument_name`, `scale_id`, `scale_name`, `reverse_keyed`, `item_min`, `item_max`, `n_categories`, `data_version`

## `04_item/item_bias_compression.csv`

bias, SD ratio, entropy, floor/ceiling, calibration slope per item.

Columns: `item_id`, `mean_human`, `mean_twin`, `bias_raw`, `bias_range`, `bias_p_wilcoxon`, `bias_q_bh`, `sd_human`, `sd_twin`, `sd_ratio`, `log_sd_ratio`, `var_ratio`, `entropy_human_norm`, `entropy_twin_norm`, `entropy_difference`, `n_categories_used_human`, `n_categories_used_twin`, `modal_category_share_twin`, `floor_share_twin`, `ceiling_share_twin`, `calibration_slope`, `calibration_intercept`, `floor_share_human`, `ceiling_share_human`, `extreme_share_human`, `extreme_share_twin`, `extreme_share_difference`, `category_coverage_ratio`, `bias_range_ci_lo`, `bias_range_ci_hi`, `sd_ratio_ci_lo`, `sd_ratio_ci_hi`, `entropy_difference_ci_lo`, `entropy_difference_ci_hi`, `instrument_id`, `instrument_name`, `scale_id`, `scale_name`, `reverse_keyed`, `item_min`, `item_max`, `n_categories`, `data_version`

## `04_item/item_range_codebook.csv`

theoretical vs observed response range per item.

Columns: `item_id`, `instrument_id`, `instrument_name`, `scale_id`, `reverse_keyed`, `theoretical_min`, `theoretical_max`, `theoretical_n_categories`, `observed_min`, `observed_max`, `observed_range_narrower_than_theoretical`, `data_version`

## `04_item/item_level_data_audit.csv`

zero variance and near-degenerate twin flags.

Columns: `item_id`, `n_categories_human`, `n_categories_twin`, `zero_variance_human`, `zero_variance_twin`, `twin_top2_categories_share`, `twin_near_degenerate_flag`, `instrument_id`, `instrument_name`, `scale_id`, `scale_name`, `reverse_keyed`, `item_min`, `item_max`, `n_categories`, `data_version`

## `04_item/item_variance_decomposition.csv`

percentages of person / item / person x item-residual variance (human, twin, error).

Columns: `instrument_id`, `instrument_name`, `side`, `n_persons`, `n_items`, `var_person`, `var_item`, `var_person_item_residual`, `var_total`, `pct_person`, `pct_item`, `pct_person_item_residual`, `data_version`

## `04_item/item_variance_decomposition_contrasts.csv`

twin minus human differences in percentage points with paired bootstrap CI.

Columns: `instrument_id`, `instrument_name`, `n_valid_reps`, `person_effect_r`, `delta_item_point`, `delta_item_ci_lo`, `delta_item_ci_hi`, `delta_resid_point`, `delta_resid_ci_lo`, `delta_resid_ci_hi`, `data_version`

## `04_item/item_loading_discrepancy.csv`

|delta loading| and loading-vector distance per item.

Columns: `instrument_id`, `item_id`, `abs_delta_loading_target`, `loading_vector_distance`, `n_factors`, `esem_admissible_both_sides`

## `04_item/item_loading_associations.csv`

Spearman association of loading discrepancy with item fidelity (admissible instruments only).

Columns: `metric`, `rho`, `p_naive`, `n`, `n_clusters`, `ci_lo`, `ci_hi`, `p_boot`, `p_below_resolution`, `n_valid_boot`, `per_instrument_rho`, `q_bh`

## `04_item/item_reverse_keyed_analysis.csv`

reverse-keyed vs direct items.

Columns: `metric`, `n_reverse`, `n_direct`, `mean_reverse`, `mean_direct`, `diff_reverse_minus_direct`, `cluster_bootstrap_ci_lo`, `cluster_bootstrap_ci_hi`, `mannwhitney_p`, `q_bh`

## `04_item/item_level_summary.csv`

item-weighted and instrument-weighted summaries.

Columns: `metric`, `n_items`, `item_weighted_mean`, `item_weighted_median`, `item_weighted_iqr_lo`, `item_weighted_iqr_hi`, `item_weighted_min`, `item_weighted_max`, `n_instruments`, `instrument_weighted_mean`, `instrument_weighted_median`

## `05_synthesis/cross_level_master.csv`

one row per scale: all levels joined.

Columns: `scale_id`, `scale_name`, `instrument_id`, `instrument_name`, `n_items_scale`, `is_nonredundant`, `r`, `pearson_r_ci_lo`, `pearson_r_ci_hi`, `pearson_q_bh_sig`, `AUC`, `AUC_ci_low`, `AUC_ci_high`, `CCC`, `Cb`, `u`, `v`, `r_dis_alpha`, `delta_ccc_ls`, `delta_ccc_ols`, `typology_r_cb`, `configural_cfa_admissible`, `metric_invariance`, `scalar_invariance`, `metric_test_available`, `scalar_test_available`, `esem_applicable`, `human_esem_admissible`, `twin_esem_admissible`, `paired_esem_admissible`, `tucker_point_available`, `tucker_interpretable`, `tucker_bootstrap_stable`, `tucker_phi_exploratory`, `tucker_phi_primary`, `tucker_phi_ci_lo`, `tucker_phi_ci_hi`, `tucker_frac_valid_bootstrap`, `structural_evidence_level`, `structural_flags`, `baseline_status`, `baseline_fit`, `threshold_decision`, `metric_decision`, `scalar_decision`, `cfa_conclusion`, `cfa_no_conclusion_reason`, `cfa_order_sensitive`, `n_items_in_item_level`, `item_median_spearman_rho`, `item_median_ccc`, `item_median_mae_range`, `item_median_abs_bias_range`, `item_median_sd_ratio`, `item_pct_sd_ratio_lt1`, `item_median_entropy_difference`, `item_pct_near_degenerate_twin`, `item_median_loading_vector_distance`, `item_median_abs_delta_loading`

## `05_synthesis/cross_level_associations.csv`

families A, B (and C2 exploratory) with rho, cluster-bootstrap CI, p, BH q.

Columns: `family`, `pair`, `x`, `y`, `rho`, `p_naive`, `n`, `n_clusters`, `ci_lo`, `ci_hi`, `p_boot`, `p_below_resolution`, `n_valid_boot`, `p_BH`

## `05_synthesis/cross_level_associations_gated_descriptive.csv`

family C1: Tucker's phi from admissible ESEM vs accuracy (descriptive).

Columns: `family`, `pair`, `x`, `y`, `rho`, `p_naive`, `n`, `n_clusters`

## `05_synthesis/cross_level_gated_scale_detail.csv`

the factors entering C1.

Columns: `scale_id`, `scale_name`, `instrument_name`, `tucker_phi_primary`, `tucker_phi_ci_lo`, `tucker_phi_ci_hi`, `tucker_bootstrap_stable`, `r`, `AUC`, `CCC`

## `05_synthesis/cross_level_invariance_group_comparison.csv`

scale accuracy by metric / scalar invariance status.

Columns: `status`, `metric`, `n_yes`, `n_no`, `median_yes`, `median_no`

## `05_synthesis/cross_level_correlation_matrix_exploratory.csv`

descriptive Spearman matrix, exploratory only.

Columns: `Unnamed: 0`, `r`, `AUC`, `CCC`, `Cb`, `r_dis_alpha`, `delta_ccc_ls`, `delta_ccc_ols`, `tucker_phi_exploratory`, `tucker_phi_primary`, `item_median_spearman_rho`, `item_median_ccc`, `item_median_mae_range`, `item_median_abs_bias_range`, `item_median_sd_ratio`, `item_median_entropy_difference`, `item_median_loading_vector_distance`

## `05_synthesis/scale_typology.csv`

one row per scale: axes I, C, S, final type, reason.

Columns: `scale_id`, `scale_name`, `instrument_name`, `r`, `AUC`, `AUC_ci_low`, `CCC`, `Cb`, `structural_evidence_level`, `paired_esem_admissible`, `tucker_phi_primary`, `tucker_phi_exploratory`, `metric_invariance`, `scalar_invariance`, `baseline_status`, `baseline_fit`, `threshold_decision`, `metric_decision`, `scalar_decision`, `axis_I_detectable`, `axis_I_practical`, `axis_C`, `axis_S`, `final_type`, `final_type_label`, `classification_reason`, `structural_flags`, `evidence_quality`

## `05_synthesis/scale_typology_summary.csv`

counts and percentages per type.

Columns: `final_type`, `n`, `pct`, `label`, `example_scales`

## `05_synthesis/cross_level_results_matrix.csv`

navigation table with the main result of every level.

Columns: `level`, `unit`, `n`, `component`, `metric`, `value`, `ci_lo`, `ci_hi`, `detail`, `status`

## Figures

Every stage that produces a figure writes it to a `figures/` subfolder **inside its own result folder** (e.g.
`02_scale/figures/`) - figures are never collected in one flat top-level folder. Each is saved as both `.png` (for
quick viewing, e.g. directly on GitHub) and `.pdf` (vector, for publication).

| file | folder | what it shows |
|---|---|---|
| `profile_ccc_own_vs_cross` | `01_profile/figures/` | histogram of profile CCC (human-based z, equal instrument weights): own twin vs a random other person's twin |
| `scale_r_vs_cb` | `02_scale/figures/` | individuation (Pearson r) against calibration (C_b) for every scale, coloured by instrument, with the typology thresholds (`r_practical`, `cb_adequate`) marked |
| `tucker_phi_admissibility` | `03_structural/figures/` | Tucker's phi with its bootstrap CI for every scale with an admissible paired ESEM (gated); colour = structural route (congruence supported / discrepant / other), marker = bootstrap-stable or not, with the congruence threshold (`phi_congruence`) marked |
| `item_fidelity_map` | `04_item/figures/` | item-level Spearman rho against SD ratio, point size = normalised MAE, coloured by instrument |
| `item_variance_decomposition` | `04_item/figures/` | ANOVA-style variance decomposition (person / item / person x item + residual) as stacked bars, Human vs Twin vs difference, per instrument |
| `scale_typology_map` | `05_synthesis/figures/` | scales in the (AUC, CCC) and (r, CCC) planes, colour = final type, marker = structural status |
| `cross_level_associations_heatmap` | `05_synthesis/figures/` | Spearman rho between item-level/structural indicators and scale-level accuracy, with BH-significant cells starred |
