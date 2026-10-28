# ============================================================================
# Corrected paired/dependent ordinal CFA measurement invariance
# Human_i <-> Twin_i
#
# PRIMARY THRESHOLD TEST: Wu & Estabrook (2016) identification
#   W0 configural
#   W1 threshold invariance
#   W2 threshold + loading invariance  [ONLY fitted/interpreted if W1 is supported]
#   W3 + intercept invariance          [ONLY if W2 is supported]
#   W4 + residual-variance invariance  [ONLY if W3 is supported]
#
# PRIMARY LOADING TEST: Millsap & Tein (2004) / Liu et al. (2017) longitudinal scheme
#   M0 baseline (includes identification threshold constraints required by scheme)
#   M1 loading invariance              [THIS is the reported loading-invariance test]
#   M2 loading + threshold invariance  [ONLY fitted/interpreted if M1 is supported]
#   M3 + unique/residual variance invariance [ONLY if M2 is supported]
#
# REPORTING PRINCIPLE:
# - threshold invariance is reported from Wu W0 -> W1;
# - loading invariance is reported separately from Millsap/Liu M0 -> M1;
# - Wu W1 -> W2 is CONDITIONAL on already-supported thresholds and is NOT used as
#   the unconditional "metric invariance rate";
# - higher-order models are not interpreted when the preceding invariance level fails.
#
# IMPORTANT DIFFERENCE FROM THE OLD SCRIPT:
# - categorical identification is handled by semTools rather than manual labels;
# - theta parameterization is used;
# - for Wu-Estabrook, thresholds are tested BEFORE loadings;
# - category labels (not only the number of categories) must match Human/Twin;
# - no silent use of min(number_of_categories_h, number_of_categories_t).
#
# References:
# Wu, H., & Estabrook, R. (2016). Psychometrika, 81, 1014-1045.
# Millsap, R. E., & Tein, J.-Y. (2004). Multivariate Behavioral Research, 39, 479-515.
# Liu, Y., et al. (2017). Psychological Methods, 22, 486-506.
# semTools::measEq.syntax documentation.
#
# Usage (called by `twinval` - see src/twinval/structural/cfa.py):
#   Rscript paired_cfa_invariance.R manifest.csv out_dir [delta_cfi_min] [delta_rmsea_max]
#
# Manifest must contain columns:
#   scale_id, file, n_items, n_pairs
# Each scale CSV must contain:
#   it1_h ... itK_h, it1_t ... itK_t
# ============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(semTools)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- if (length(args) >= 1) args[[1]] else "cfa_manifest.csv"
out_dir <- if (length(args) >= 2) args[[2]] else "cfa_results"

# Conventional change-in-fit criteria used ONLY to create transparent status flags.
# Keep the raw fit indices/deltas in the output and replace these cutoffs if the
# preregistered/target-journal criteria differ. DG-DMI, if used, should be run
# separately on the corrected model specifications.
DELTA_CFI_MIN <- if (length(args) >= 3) as.numeric(args[[3]]) else -0.010   # support if delta CFI >= -0.010
DELTA_RMSEA_MAX <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.015  # support if delta RMSEA <= +0.015

# ------------------------------
# Output directories
# ------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (scheme in c("wu", "millsap")) {
  for (sub in c("model_syntax", "model_summaries", "parameter_estimates", "thresholds")) {
    dir.create(file.path(out_dir, scheme, sub), showWarnings = FALSE, recursive = TRUE)
  }
}

# Save package/session versions for reproducibility
sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
required_manifest_cols <- c("scale_id", "file", "n_items", "n_pairs")
miss_cols <- setdiff(required_manifest_cols, names(manifest))
if (length(miss_cols) > 0) {
  stop("Manifest is missing required columns: ", paste(miss_cols, collapse = ", "))
}

# ============================================================================
# Helpers
# ============================================================================

# Stable ordering of observed response codes. Numeric-like labels are ordered
# numerically; otherwise lexical ordering is used.
observed_levels <- function(x) {
  vals <- unique(as.character(x[!is.na(x)]))
  if (length(vals) == 0) return(character(0))
  num <- suppressWarnings(as.numeric(vals))
  if (all(!is.na(num))) vals[order(num)] else sort(vals)
}

# Harmonize Human/Twin item categories STRICTLY.
# We do not silently collapse categories because equality of threshold t_j only
# has a meaningful interpretation when the same observed cut points are present.
prepare_ordered_pairs <- function(d, h_items, t_items) {
  diagnostics <- vector("list", length(h_items))
  mismatch <- integer(0)
  too_few <- integer(0)

  for (i in seq_along(h_items)) {
    lh <- observed_levels(d[[h_items[i]]])
    lt <- observed_levels(d[[t_items[i]]])

    diagnostics[[i]] <- data.frame(
      item_index = i,
      human_levels = paste(lh, collapse = "|"),
      twin_levels = paste(lt, collapse = "|"),
      ncat_human = length(lh),
      ncat_twin = length(lt),
      labels_identical = identical(lh, lt),
      stringsAsFactors = FALSE
    )

    if (length(lh) < 2 || length(lt) < 2) too_few <- c(too_few, i)
    if (!identical(lh, lt)) mismatch <- c(mismatch, i)

    if (identical(lh, lt) && length(lh) >= 2) {
      # Force exactly the same ordered level coding on both measurements.
      d[[h_items[i]]] <- ordered(as.character(d[[h_items[i]]]), levels = lh)
      d[[t_items[i]]] <- ordered(as.character(d[[t_items[i]]]), levels = lh)
    }
  }

  list(
    data = d,
    diagnostics = do.call(rbind, diagnostics),
    mismatch = mismatch,
    too_few = too_few
  )
}

build_base_model <- function(k, marker_index = NULL) {
  h_items <- paste0("it", seq_len(k), "_h")
  t_items <- paste0("it", seq_len(k), "_t")

  # For Millsap ID.fac="marker", the first indicator in each factor is the marker.
  # Reordering is only used when marker_index is explicitly supplied.
  ord <- seq_len(k)
  if (!is.null(marker_index)) {
    ord <- c(marker_index, setdiff(ord, marker_index))
  }

  paste0(
    "Fh =~ ", paste(h_items[ord], collapse = " + "), "\n",
    "Ft =~ ", paste(t_items[ord], collapse = " + ")
  )
}

make_long_names <- function(k) {
  h_items <- paste0("it", seq_len(k), "_h")
  t_items <- paste0("it", seq_len(k), "_t")
  long_ind <- setNames(
    lapply(seq_len(k), function(i) c(h_items[i], t_items[i])),
    paste0("item", seq_len(k))
  )
  list(
    longFacNames = list(F = c("Fh", "Ft")),
    longIndNames = long_ind
  )
}

fit_generated_syntax <- function(syntax_obj, d, all_items) {
  model_txt <- as.character(syntax_obj)
  tryCatch(
    cfa(
      model_txt,
      data = d,
      ordered = all_items,
      estimator = "WLSMV",
      parameterization = "theta"
    ),
    error = function(e) e
  )
}

fit_ok <- function(fit) {
  inherits(fit, "lavaan") && isTRUE(tryCatch(lavInspect(fit, "converged"), error = function(e) FALSE))
}

get_fit <- function(fit) {
  if (!fit_ok(fit)) return(NULL)
  tryCatch(
    fitMeasures(fit, c(
      "chisq.scaled", "df.scaled", "pvalue.scaled",
      "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"
    )),
    error = function(e) NULL
  )
}

has_heywood <- function(fit) {
  if (!fit_ok(fit)) return(NA)
  pe <- tryCatch(parameterEstimates(fit), error = function(e) NULL)
  if (is.null(pe)) return(NA)
  vr <- pe[pe$op == "~~" & pe$lhs == pe$rhs, , drop = FALSE]
  any(vr$est < 0, na.rm = TRUE)
}

# Save the full summary, full parameter table, and thresholds.
save_fit_details <- function(fit, scale_id, scheme, model_name) {
  safe_id <- gsub("[^[:alnum:]_.-]", "_", scale_id)

  summ_path <- file.path(out_dir, scheme, "model_summaries",
                         paste0(safe_id, "_", model_name, ".txt"))
  sink(summ_path)
  if (inherits(fit, "error")) {
    cat("ERROR:\n", conditionMessage(fit), "\n")
  } else {
    print(summary(fit, fit.measures = TRUE, standardized = TRUE))
  }
  sink()

  if (fit_ok(fit)) {
    pe <- parameterEstimates(fit, standardized = TRUE)
    write.csv(
      pe,
      file.path(out_dir, scheme, "parameter_estimates",
                paste0(safe_id, "_", model_name, ".csv")),
      row.names = FALSE
    )
    thr <- pe[pe$op == "|", , drop = FALSE]
    write.csv(
      thr,
      file.path(out_dir, scheme, "thresholds",
                paste0(safe_id, "_", model_name, "_thresholds.csv")),
      row.names = FALSE
    )
  }
}

save_syntax <- function(syntax_obj, scale_id, scheme, model_name) {
  safe_id <- gsub("[^[:alnum:]_.-]", "_", scale_id)
  writeLines(
    as.character(syntax_obj),
    con = file.path(out_dir, scheme, "model_syntax",
                    paste0(safe_id, "_", model_name, ".lav"))
  )
}

# Robust WLSMV nested-model comparison.
compare_nested <- function(fit_less, fit_more) {
  out <- list(p = NA_real_, chisq_diff = NA_real_, df_diff = NA_real_)
  if (!fit_ok(fit_less) || !fit_ok(fit_more)) return(out)
  tab <- tryCatch(lavTestLRT(fit_less, fit_more), error = function(e) NULL)
  if (is.null(tab) || nrow(tab) < 2) return(out)

  pcol <- grep("Pr\\(>Chisq\\)", names(tab), value = TRUE)
  if (length(pcol) > 0) out$p <- as.numeric(tab[[pcol[1]]][2])

  # Depending on lavaan version, difference columns can vary. Compute df diff
  # from model df as a stable descriptive fallback.
  fm1 <- get_fit(fit_less)
  fm2 <- get_fit(fit_more)
  if (!is.null(fm1) && !is.null(fm2)) {
    out$df_diff <- as.numeric(fm2[["df.scaled"]] - fm1[["df.scaled"]])
  }
  # If a Chisq diff column is printed by lavTestLRT, capture it.
  cands <- grep("Chisq diff|Chisq Diff|Chisq.diff", names(tab), value = TRUE)
  if (length(cands) > 0) out$chisq_diff <- as.numeric(tab[[cands[1]]][2])
  out
}

fit_delta <- function(fit_less, fit_more) {
  a <- get_fit(fit_less)
  b <- get_fit(fit_more)
  if (is.null(a) || is.null(b)) {
    return(c(delta_cfi = NA_real_, delta_rmsea = NA_real_, delta_srmr = NA_real_))
  }
  c(
    delta_cfi = as.numeric(b[["cfi.scaled"]] - a[["cfi.scaled"]]),
    delta_rmsea = as.numeric(b[["rmsea.scaled"]] - a[["rmsea.scaled"]]),
    delta_srmr = as.numeric(b[["srmr"]] - a[["srmr"]])
  )
}

# Decide whether an invariance step is supported using conventional change-in-fit
# cutoffs.  The function deliberately distinguishes REJECTED from INADMISSIBLE
# and NOT_TESTABLE.  A constrained Heywood solution is NOT counted as evidence
# against invariance; it is treated as an inadmissible model.
decide_step <- function(fit_less, fit_more,
                        delta_cfi_min = DELTA_CFI_MIN,
                        delta_rmsea_max = DELTA_RMSEA_MAX) {
  if (!fit_ok(fit_less) || !fit_ok(fit_more)) return("not_testable_nonconvergence")
  if (isTRUE(has_heywood(fit_less)) || isTRUE(has_heywood(fit_more))) return("inadmissible_heywood")

  d <- fit_delta(fit_less, fit_more)
  lrt <- compare_nested(fit_less, fit_more)

  # If the equality constraints add no positive df, the step cannot provide an
  # invariance test (common for very short/saturated scales).
  if (!is.finite(lrt$df_diff) || lrt$df_diff <= 0) return("not_testable_df")
  if (!is.finite(d[["delta_cfi"]]) || !is.finite(d[["delta_rmsea"]])) return("not_testable_fit_indices")

  if (d[["delta_cfi"]] >= delta_cfi_min && d[["delta_rmsea"]] <= delta_rmsea_max) {
    "supported"
  } else {
    "rejected"
  }
}

# Create a placeholder error object so downstream save/summary functions can
# handle models intentionally skipped by the hierarchical testing rule.
skipped_model <- function(reason) {
  structure(list(message = paste0("SKIPPED: ", reason)), class = c("simpleError", "error", "condition"))
}

extract_factor_metrics <- function(fit) {
  empty <- list(
    latent_cor_raw = NA_real_, latent_cor_oriented = NA_real_,
    tucker_phi = NA_real_,
    cr_h = NA_real_, cr_t = NA_real_, ave_h = NA_real_, ave_t = NA_real_,
    n_negative_loadings_h = NA_integer_, n_negative_loadings_t = NA_integer_
  )
  if (!fit_ok(fit)) return(empty)

  ss <- tryCatch(standardizedSolution(fit), error = function(e) NULL)
  if (is.null(ss)) return(empty)

  lh_df <- ss[ss$op == "=~" & ss$lhs == "Fh", , drop = FALSE]
  lt_df <- ss[ss$op == "=~" & ss$lhs == "Ft", , drop = FALSE]
  if (nrow(lh_df) == 0 || nrow(lt_df) == 0) return(empty)

  # Match loadings by item index rather than by position in syntax.
  idx_h <- sub("^it([0-9]+)_h$", "\\1", lh_df$rhs)
  idx_t <- sub("^it([0-9]+)_t$", "\\1", lt_df$rhs)
  common <- intersect(idx_h, idx_t)
  lh <- sapply(common, function(z) lh_df$est.std[match(z, idx_h)])
  lt <- sapply(common, function(z) lt_df$est.std[match(z, idx_t)])

  # Orient each loading vector so its average loading is positive. Factor sign
  # is arbitrary; this prevents a pure sign flip from making Tucker phi negative.
  sh <- ifelse(mean(lh, na.rm = TRUE) < 0, -1, 1)
  st <- ifelse(mean(lt, na.rm = TRUE) < 0, -1, 1)
  lh_o <- sh * lh
  lt_o <- st * lt

  denom <- sqrt(sum(lh_o^2, na.rm = TRUE) * sum(lt_o^2, na.rm = TRUE))
  phi <- ifelse(is.finite(denom) && denom > 0,
                sum(lh_o * lt_o, na.rm = TRUE) / denom,
                NA_real_)

  # Latent factor covariance/correlation from standardized solution.
  lr <- ss[ss$op == "~~" & ss$lhs == "Fh" & ss$rhs == "Ft", , drop = FALSE]
  if (nrow(lr) == 0) lr <- ss[ss$op == "~~" & ss$lhs == "Ft" & ss$rhs == "Fh", , drop = FALSE]
  latent_raw <- if (nrow(lr) > 0) lr$est.std[1] else NA_real_
  latent_oriented <- ifelse(is.na(latent_raw), NA_real_, latent_raw * sh * st)

  # CR/AVE from standardized latent-response loadings and residual variances.
  # For ordered indicators under theta, these are metrics of the underlying
  # continuous response variables. Cross-time residual covariances do not enter
  # within-source CR/AVE denominators.
  resid <- ss[ss$op == "~~" & ss$lhs == ss$rhs, , drop = FALSE]
  theta_h <- sapply(common, function(z) {
    nm <- paste0("it", z, "_h")
    rr <- resid[resid$lhs == nm, , drop = FALSE]
    if (nrow(rr) > 0) rr$est.std[1] else NA_real_
  })
  theta_t <- sapply(common, function(z) {
    nm <- paste0("it", z, "_t")
    rr <- resid[resid$lhs == nm, , drop = FALSE]
    if (nrow(rr) > 0) rr$est.std[1] else NA_real_
  })

  calc_cr <- function(lam, theta) {
    ok <- is.finite(lam) & is.finite(theta)
    if (sum(ok) < 2) return(NA_real_)
    (sum(lam[ok])^2) / (sum(lam[ok])^2 + sum(theta[ok]))
  }
  calc_ave <- function(lam, theta) {
    ok <- is.finite(lam) & is.finite(theta)
    if (sum(ok) < 2) return(NA_real_)
    sum(lam[ok]^2) / (sum(lam[ok]^2) + sum(theta[ok]))
  }

  list(
    latent_cor_raw = latent_raw,
    latent_cor_oriented = latent_oriented,
    tucker_phi = phi,
    cr_h = calc_cr(lh_o, theta_h),
    cr_t = calc_cr(lt_o, theta_t),
    ave_h = calc_ave(lh_o, theta_h),
    ave_t = calc_ave(lt_o, theta_t),
    n_negative_loadings_h = sum(lh_o < 0, na.rm = TRUE),
    n_negative_loadings_t = sum(lt_o < 0, na.rm = TRUE)
  )
}

fit_summary_row <- function(fit, prefix) {
  fm <- get_fit(fit)
  vals <- c(
    converged = fit_ok(fit),
    cfi = if (is.null(fm)) NA_real_ else fm[["cfi.scaled"]],
    rmsea = if (is.null(fm)) NA_real_ else fm[["rmsea.scaled"]],
    srmr = if (is.null(fm)) NA_real_ else fm[["srmr"]],
    df = if (is.null(fm)) NA_real_ else fm[["df.scaled"]],
    heywood = has_heywood(fit),
    fit_diagnostics(fit)
  )
  names(vals) <- paste0(prefix, "_", names(vals))
  vals
}

# Model-level diagnostics used to separate estimability, absolute fit and invariance.
# Reported only; they do not enter any invariance decision of this script.
#   chisq/npar/tli/rmsea CI : absolute fit (scaled versions, as the decision indices)
#   post_check              : lavaan check of admissible variances / positive-definite latent matrices
#   vcov_min_eigen          : smallest eigenvalue of the parameter covariance matrix (identification)
#   theta/psi_min_eigen     : smallest eigenvalue of the residual / latent covariance matrix (< 0: not PD)
#   max_sample_cor_h/_t     : largest sample polychoric correlation among Human / among Twin items
#   max_std_loading, min_std_resid_var, max_abs_resid_cor : boundary diagnostics (std.all)
#   abs_factor_cor          : |cor(Fh, Ft)|
fit_diagnostics <- function(fit) {
  out <- c(chisq = NA_real_, npar = NA_real_, tli = NA_real_, rmsea_ci_lo = NA_real_, rmsea_ci_hi = NA_real_,
           post_check = NA, vcov_min_eigen = NA_real_, theta_min_eigen = NA_real_, psi_min_eigen = NA_real_,
           max_sample_cor_h = NA_real_, max_sample_cor_t = NA_real_, max_std_loading = NA_real_,
           min_std_resid_var = NA_real_, max_abs_resid_cor = NA_real_, abs_factor_cor = NA_real_)
  if (!fit_ok(fit)) return(out)
  one <- function(name) tryCatch(as.numeric(fitMeasures(fit, name)), error = function(e) NA_real_)
  out[["chisq"]] <- one("chisq.scaled")
  out[["npar"]] <- one("npar")
  out[["tli"]] <- one("tli.scaled")
  out[["rmsea_ci_lo"]] <- one("rmsea.ci.lower.scaled")
  out[["rmsea_ci_hi"]] <- one("rmsea.ci.upper.scaled")
  out[["post_check"]] <- tryCatch(isTRUE(suppressWarnings(lavInspect(fit, "post.check"))), error = function(e) NA)
  # equality-constrained parameters appear as identical columns of vcov; keep one copy of each
  out[["vcov_min_eigen"]] <- tryCatch({
    v <- suppressWarnings(lavInspect(fit, "vcov"))
    keep <- !duplicated(t(round(v, 10)))
    min(eigen(v[keep, keep, drop = FALSE], symmetric = TRUE, only.values = TRUE)$values)
  }, error = function(e) NA_real_)
  est <- tryCatch(lavInspect(fit, "est"), error = function(e) NULL)
  min_eig <- function(m) tryCatch(min(eigen(m, symmetric = TRUE, only.values = TRUE)$values),
                                  error = function(e) NA_real_)
  if (!is.null(est)) {
    out[["theta_min_eigen"]] <- min_eig(est$theta)
    out[["psi_min_eigen"]] <- min_eig(est$psi)
  }
  # largest off-diagonal sample polychoric correlation within the Human and within the Twin items
  sc <- tryCatch(lavInspect(fit, "sampstat")$cov, error = function(e) NULL)
  if (!is.null(sc)) {
    side_max <- function(suffix) {
      ix <- grep(paste0(suffix, "$"), rownames(sc))
      if (length(ix) < 2) return(NA_real_)
      m <- abs(sc[ix, ix]); diag(m) <- NA
      max(m, na.rm = TRUE)
    }
    out[["max_sample_cor_h"]] <- side_max("_h")
    out[["max_sample_cor_t"]] <- side_max("_t")
  }
  ss <- tryCatch(standardizedSolution(fit), error = function(e) NULL)
  if (!is.null(ss)) {
    items <- ss$op == "~~" & grepl("^it", ss$lhs) & grepl("^it", ss$rhs)
    lam <- ss$est.std[ss$op == "=~"]
    rv <- ss$est.std[items & ss$lhs == ss$rhs]
    rc <- ss$est.std[items & ss$lhs != ss$rhs]
    fc <- ss$est.std[ss$op == "~~" & ss$lhs %in% c("Fh", "Ft") & ss$rhs %in% c("Fh", "Ft") & ss$lhs != ss$rhs]
    if (length(lam)) out[["max_std_loading"]] <- max(abs(lam), na.rm = TRUE)
    if (length(rv)) out[["min_std_resid_var"]] <- min(rv, na.rm = TRUE)
    if (length(rc)) out[["max_abs_resid_cor"]] <- max(abs(rc), na.rm = TRUE)
    if (length(fc)) out[["abs_factor_cor"]] <- abs(fc[1])
  }
  out
}

# ============================================================================
# Fit one scale: Wu threshold test + Millsap/Liu loading test
# ============================================================================
fit_one_scale <- function(scale_id, data_path, k, n_pairs) {
  if (k < 2) stop("At least 2 indicators are required. scale_id=", scale_id)
  d <- read.csv(data_path, stringsAsFactors = FALSE)

  h_items <- paste0("it", seq_len(k), "_h")
  t_items <- paste0("it", seq_len(k), "_t")
  all_items <- c(h_items, t_items)
  miss <- setdiff(all_items, names(d))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  prep <- prepare_ordered_pairs(d, h_items, t_items)
  d <- prep$data

  safe_id <- gsub("[^[:alnum:]_.-]", "_", scale_id)
  write.csv(prep$diagnostics,
            file.path(out_dir, paste0(safe_id, "_category_diagnostics.csv")),
            row.names = FALSE)

  if (length(prep$too_few) > 0) {
    stop("Item(s) with <2 observed categories: ", paste(prep$too_few, collapse = ";"))
  }
  if (length(prep$mismatch) > 0) {
    stop(
      "Human/Twin observed category LABELS differ for item(s): ",
      paste(prep$mismatch, collapse = ";"),
      ". Harmonize/collapse categories symmetrically before invariance testing."
    )
  }

  long_names <- make_long_names(k)
  base_model <- build_base_model(k)

  # --------------------------------------------------------------------------
  # PRIMARY THRESHOLD TEST: Wu & Estabrook 2016
  # --------------------------------------------------------------------------
  # IMPORTANT: the unconditional Wu result reported in the manuscript is W0->W1
  # (threshold invariance). W1->W2 is a CONDITIONAL loading test and is only
  # estimated if threshold invariance is first supported.
  wu_specs <- list(
    configural = character(0),
    thresholds = c("thresholds"),
    conditional_loading = c("thresholds", "loadings"),
    conditional_scalar = c("thresholds", "loadings", "intercepts"),
    conditional_strict = c("thresholds", "loadings", "intercepts", "residuals")
  )

  wu_syn <- list()
  wu_fit <- list()

  # W0 configural
  wu_syn$configural <- measEq.syntax(
    configural.model = base_model,
    data = d,
    ordered = all_items,
    parameterization = "theta",
    ID.fac = "std.lv",
    ID.cat = "Wu.Estabrook.2016",
    longFacNames = long_names$longFacNames,
    longIndNames = long_names$longIndNames,
    long.equal = wu_specs$configural,
    auto = "all",
    return.fit = FALSE
  )
  save_syntax(wu_syn$configural, scale_id, "wu", "configural")
  wu_fit$configural <- fit_generated_syntax(wu_syn$configural, d, all_items)
  save_fit_details(wu_fit$configural, scale_id, "wu", "configural")

  # W1 threshold invariance -- PRIMARY Wu test
  wu_syn$thresholds <- measEq.syntax(
    configural.model = base_model,
    data = d,
    ordered = all_items,
    parameterization = "theta",
    ID.fac = "std.lv",
    ID.cat = "Wu.Estabrook.2016",
    longFacNames = long_names$longFacNames,
    longIndNames = long_names$longIndNames,
    long.equal = wu_specs$thresholds,
    auto = "all",
    return.fit = FALSE
  )
  save_syntax(wu_syn$thresholds, scale_id, "wu", "thresholds")
  wu_fit$thresholds <- fit_generated_syntax(wu_syn$thresholds, d, all_items)
  save_fit_details(wu_fit$thresholds, scale_id, "wu", "thresholds")

  wu_threshold_status <- decide_step(wu_fit$configural, wu_fit$thresholds)

  # Higher Wu levels are hierarchical/conditional.  They are NOT used to
  # estimate the unconditional loading-invariance rate in the paper.
  if (identical(wu_threshold_status, "supported")) {
    wu_syn$conditional_loading <- measEq.syntax(
      configural.model = base_model,
      data = d,
      ordered = all_items,
      parameterization = "theta",
      ID.fac = "std.lv",
      ID.cat = "Wu.Estabrook.2016",
      longFacNames = long_names$longFacNames,
      longIndNames = long_names$longIndNames,
      long.equal = wu_specs$conditional_loading,
      auto = "all",
      return.fit = FALSE
    )
    save_syntax(wu_syn$conditional_loading, scale_id, "wu", "conditional_loading")
    wu_fit$conditional_loading <- fit_generated_syntax(wu_syn$conditional_loading, d, all_items)
    save_fit_details(wu_fit$conditional_loading, scale_id, "wu", "conditional_loading")
    wu_cond_loading_status <- decide_step(wu_fit$thresholds, wu_fit$conditional_loading)
  } else {
    wu_fit$conditional_loading <- skipped_model(paste0("threshold invariance status = ", wu_threshold_status))
    wu_cond_loading_status <- "not_interpreted_threshold_not_supported"
  }

  if (identical(wu_cond_loading_status, "supported")) {
    wu_syn$conditional_scalar <- measEq.syntax(
      configural.model = base_model,
      data = d,
      ordered = all_items,
      parameterization = "theta",
      ID.fac = "std.lv",
      ID.cat = "Wu.Estabrook.2016",
      longFacNames = long_names$longFacNames,
      longIndNames = long_names$longIndNames,
      long.equal = wu_specs$conditional_scalar,
      auto = "all",
      return.fit = FALSE
    )
    save_syntax(wu_syn$conditional_scalar, scale_id, "wu", "conditional_scalar")
    wu_fit$conditional_scalar <- fit_generated_syntax(wu_syn$conditional_scalar, d, all_items)
    save_fit_details(wu_fit$conditional_scalar, scale_id, "wu", "conditional_scalar")
    wu_cond_scalar_status <- decide_step(wu_fit$conditional_loading, wu_fit$conditional_scalar)
  } else {
    wu_fit$conditional_scalar <- skipped_model(paste0("conditional loading status = ", wu_cond_loading_status))
    wu_cond_scalar_status <- "not_interpreted_prior_level_not_supported"
  }

  if (identical(wu_cond_scalar_status, "supported")) {
    wu_syn$conditional_strict <- measEq.syntax(
      configural.model = base_model,
      data = d,
      ordered = all_items,
      parameterization = "theta",
      ID.fac = "std.lv",
      ID.cat = "Wu.Estabrook.2016",
      longFacNames = long_names$longFacNames,
      longIndNames = long_names$longIndNames,
      long.equal = wu_specs$conditional_strict,
      auto = "all",
      return.fit = FALSE
    )
    save_syntax(wu_syn$conditional_strict, scale_id, "wu", "conditional_strict")
    wu_fit$conditional_strict <- fit_generated_syntax(wu_syn$conditional_strict, d, all_items)
    save_fit_details(wu_fit$conditional_strict, scale_id, "wu", "conditional_strict")
    wu_cond_strict_status <- decide_step(wu_fit$conditional_scalar, wu_fit$conditional_strict)
  } else {
    wu_fit$conditional_strict <- skipped_model(paste0("conditional scalar status = ", wu_cond_scalar_status))
    wu_cond_strict_status <- "not_interpreted_prior_level_not_supported"
  }

  wu_d1 <- fit_delta(wu_fit$configural, wu_fit$thresholds)
  wu_lrt1 <- compare_nested(wu_fit$configural, wu_fit$thresholds)
  wu_d2 <- fit_delta(wu_fit$thresholds, wu_fit$conditional_loading)
  wu_lrt2 <- compare_nested(wu_fit$thresholds, wu_fit$conditional_loading)
  wu_d3 <- fit_delta(wu_fit$conditional_loading, wu_fit$conditional_scalar)
  wu_lrt3 <- compare_nested(wu_fit$conditional_loading, wu_fit$conditional_scalar)
  wu_d4 <- fit_delta(wu_fit$conditional_scalar, wu_fit$conditional_strict)
  wu_lrt4 <- compare_nested(wu_fit$conditional_scalar, wu_fit$conditional_strict)
  wu_fm <- extract_factor_metrics(wu_fit$configural)

  # --------------------------------------------------------------------------
  # PRIMARY LOADING TEST: Millsap & Tein 2004 / Liu et al. 2017
  # --------------------------------------------------------------------------
  # The manuscript's unconditional loading-invariance result comes from M0->M1.
  # Subsequent threshold and unique-variance tests are hierarchical and are only
  # interpreted after the preceding Millsap/Liu level is supported.
  millsap_base <- build_base_model(k, marker_index = 1L)
  millsap_specs <- list(
    baseline = character(0),
    loading = c("loadings"),
    threshold_after_loading = c("loadings", "thresholds"),
    unique_after_threshold = c("loadings", "thresholds", "residuals")
  )

  ms_syn <- list()
  ms_fit <- list()

  # M0 baseline
  ms_syn$baseline <- measEq.syntax(
    configural.model = millsap_base,
    data = d,
    ordered = all_items,
    parameterization = "theta",
    ID.fac = "marker",
    ID.cat = "Millsap.Tein.2004",
    ID.thr = c(1L, 2L),
    longFacNames = long_names$longFacNames,
    longIndNames = long_names$longIndNames,
    long.equal = millsap_specs$baseline,
    auto = "all",
    return.fit = FALSE
  )
  save_syntax(ms_syn$baseline, scale_id, "millsap", "baseline")
  ms_fit$baseline <- fit_generated_syntax(ms_syn$baseline, d, all_items)
  save_fit_details(ms_fit$baseline, scale_id, "millsap", "baseline")

  # M1 loading invariance -- PRIMARY loading test
  ms_syn$loading <- measEq.syntax(
    configural.model = millsap_base,
    data = d,
    ordered = all_items,
    parameterization = "theta",
    ID.fac = "marker",
    ID.cat = "Millsap.Tein.2004",
    ID.thr = c(1L, 2L),
    longFacNames = long_names$longFacNames,
    longIndNames = long_names$longIndNames,
    long.equal = millsap_specs$loading,
    auto = "all",
    return.fit = FALSE
  )
  save_syntax(ms_syn$loading, scale_id, "millsap", "loading")
  ms_fit$loading <- fit_generated_syntax(ms_syn$loading, d, all_items)
  save_fit_details(ms_fit$loading, scale_id, "millsap", "loading")

  ms_loading_status <- decide_step(ms_fit$baseline, ms_fit$loading)

  if (identical(ms_loading_status, "supported")) {
    ms_syn$threshold_after_loading <- measEq.syntax(
      configural.model = millsap_base,
      data = d,
      ordered = all_items,
      parameterization = "theta",
      ID.fac = "marker",
      ID.cat = "Millsap.Tein.2004",
      ID.thr = c(1L, 2L),
      longFacNames = long_names$longFacNames,
      longIndNames = long_names$longIndNames,
      long.equal = millsap_specs$threshold_after_loading,
      auto = "all",
      return.fit = FALSE
    )
    save_syntax(ms_syn$threshold_after_loading, scale_id, "millsap", "threshold_after_loading")
    ms_fit$threshold_after_loading <- fit_generated_syntax(ms_syn$threshold_after_loading, d, all_items)
    save_fit_details(ms_fit$threshold_after_loading, scale_id, "millsap", "threshold_after_loading")
    ms_threshold_status <- decide_step(ms_fit$loading, ms_fit$threshold_after_loading)
  } else {
    ms_fit$threshold_after_loading <- skipped_model(paste0("loading invariance status = ", ms_loading_status))
    ms_threshold_status <- "not_interpreted_loading_not_supported"
  }

  if (identical(ms_threshold_status, "supported")) {
    ms_syn$unique_after_threshold <- measEq.syntax(
      configural.model = millsap_base,
      data = d,
      ordered = all_items,
      parameterization = "theta",
      ID.fac = "marker",
      ID.cat = "Millsap.Tein.2004",
      ID.thr = c(1L, 2L),
      longFacNames = long_names$longFacNames,
      longIndNames = long_names$longIndNames,
      long.equal = millsap_specs$unique_after_threshold,
      auto = "all",
      return.fit = FALSE
    )
    save_syntax(ms_syn$unique_after_threshold, scale_id, "millsap", "unique_after_threshold")
    ms_fit$unique_after_threshold <- fit_generated_syntax(ms_syn$unique_after_threshold, d, all_items)
    save_fit_details(ms_fit$unique_after_threshold, scale_id, "millsap", "unique_after_threshold")
    ms_unique_status <- decide_step(ms_fit$threshold_after_loading, ms_fit$unique_after_threshold)
  } else {
    ms_fit$unique_after_threshold <- skipped_model(paste0("threshold-after-loading status = ", ms_threshold_status))
    ms_unique_status <- "not_interpreted_prior_level_not_supported"
  }

  ms_d1 <- fit_delta(ms_fit$baseline, ms_fit$loading)
  ms_lrt1 <- compare_nested(ms_fit$baseline, ms_fit$loading)
  ms_d2 <- fit_delta(ms_fit$loading, ms_fit$threshold_after_loading)
  ms_lrt2 <- compare_nested(ms_fit$loading, ms_fit$threshold_after_loading)
  ms_d3 <- fit_delta(ms_fit$threshold_after_loading, ms_fit$unique_after_threshold)
  ms_lrt3 <- compare_nested(ms_fit$threshold_after_loading, ms_fit$unique_after_threshold)
  ms_fm <- extract_factor_metrics(ms_fit$baseline)

  # --------------------------------------------------------------------------
  # One-row summary
  # --------------------------------------------------------------------------
  out <- c(
    scale_id = scale_id,
    k_items = k,
    n_pairs_manifest = n_pairs,
    n_rows_data = nrow(d),
    category_mismatch_items = paste(prep$mismatch, collapse = ";"),
    wu_configural_id = "Wu.Estabrook.2016",
    millsap_marker_item = "it1"
  )

  # Fit summaries
  for (nm in names(wu_fit)) out <- c(out, fit_summary_row(wu_fit[[nm]], paste0("wu_", nm)))
  for (nm in names(ms_fit)) out <- c(out, fit_summary_row(ms_fit[[nm]], paste0("ms_", nm)))

  # Wu deltas, LRT, and hierarchical reporting statuses.
  # Only wu_threshold_status is an unconditional primary invariance result.
  out <- c(out,
    wu_threshold_status = wu_threshold_status,
    wu_threshold_delta_cfi = wu_d1[["delta_cfi"]],
    wu_threshold_delta_rmsea = wu_d1[["delta_rmsea"]],
    wu_threshold_delta_srmr = wu_d1[["delta_srmr"]],
    wu_threshold_lrt_p = wu_lrt1$p,
    wu_threshold_df_diff = wu_lrt1$df_diff,

    wu_conditional_loading_status = wu_cond_loading_status,
    wu_conditional_loading_delta_cfi = wu_d2[["delta_cfi"]],
    wu_conditional_loading_delta_rmsea = wu_d2[["delta_rmsea"]],
    wu_conditional_loading_delta_srmr = wu_d2[["delta_srmr"]],
    wu_conditional_loading_lrt_p = wu_lrt2$p,
    wu_conditional_loading_df_diff = wu_lrt2$df_diff,

    wu_conditional_scalar_status = wu_cond_scalar_status,
    wu_conditional_scalar_delta_cfi = wu_d3[["delta_cfi"]],
    wu_conditional_scalar_delta_rmsea = wu_d3[["delta_rmsea"]],
    wu_conditional_scalar_delta_srmr = wu_d3[["delta_srmr"]],
    wu_conditional_scalar_lrt_p = wu_lrt3$p,
    wu_conditional_scalar_df_diff = wu_lrt3$df_diff,

    wu_conditional_strict_status = wu_cond_strict_status,
    wu_conditional_strict_delta_cfi = wu_d4[["delta_cfi"]],
    wu_conditional_strict_delta_rmsea = wu_d4[["delta_rmsea"]],
    wu_conditional_strict_delta_srmr = wu_d4[["delta_srmr"]],
    wu_conditional_strict_lrt_p = wu_lrt4$p,
    wu_conditional_strict_df_diff = wu_lrt4$df_diff,

    wu_latent_cor_raw = wu_fm$latent_cor_raw,
    wu_latent_cor_oriented = wu_fm$latent_cor_oriented,
    wu_tucker_phi = wu_fm$tucker_phi,
    wu_cr_h = wu_fm$cr_h,
    wu_cr_t = wu_fm$cr_t,
    wu_ave_h = wu_fm$ave_h,
    wu_ave_t = wu_fm$ave_t,
    wu_negative_loadings_h = wu_fm$n_negative_loadings_h,
    wu_negative_loadings_t = wu_fm$n_negative_loadings_t
  )

  # Millsap/Liu deltas, LRT, and hierarchical reporting statuses.
  # ms_loading_status (M0->M1) is the manuscript's primary loading-invariance result.
  out <- c(out,
    ms_loading_status = ms_loading_status,
    ms_loading_delta_cfi = ms_d1[["delta_cfi"]],
    ms_loading_delta_rmsea = ms_d1[["delta_rmsea"]],
    ms_loading_delta_srmr = ms_d1[["delta_srmr"]],
    ms_loading_lrt_p = ms_lrt1$p,
    ms_loading_df_diff = ms_lrt1$df_diff,

    ms_threshold_after_loading_status = ms_threshold_status,
    ms_threshold_after_loading_delta_cfi = ms_d2[["delta_cfi"]],
    ms_threshold_after_loading_delta_rmsea = ms_d2[["delta_rmsea"]],
    ms_threshold_after_loading_delta_srmr = ms_d2[["delta_srmr"]],
    ms_threshold_after_loading_lrt_p = ms_lrt2$p,
    ms_threshold_after_loading_df_diff = ms_lrt2$df_diff,

    ms_unique_after_threshold_status = ms_unique_status,
    ms_unique_after_threshold_delta_cfi = ms_d3[["delta_cfi"]],
    ms_unique_after_threshold_delta_rmsea = ms_d3[["delta_rmsea"]],
    ms_unique_after_threshold_delta_srmr = ms_d3[["delta_srmr"]],
    ms_unique_after_threshold_lrt_p = ms_lrt3$p,
    ms_unique_after_threshold_df_diff = ms_lrt3$df_diff,

    ms_latent_cor_raw = ms_fm$latent_cor_raw,
    ms_latent_cor_oriented = ms_fm$latent_cor_oriented,
    ms_tucker_phi = ms_fm$tucker_phi,
    ms_cr_h = ms_fm$cr_h,
    ms_cr_t = ms_fm$cr_t,
    ms_ave_h = ms_fm$ave_h,
    ms_ave_t = ms_fm$ave_t
  )

  as.data.frame(as.list(out), stringsAsFactors = FALSE)
}

# ============================================================================
# Main loop
# ============================================================================
results <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  cat(sprintf("[%d/%d] %s (k=%d, n=%d)\n",
              i, nrow(manifest), row$scale_id, row$n_items, row$n_pairs))

  results[[i]] <- tryCatch(
    fit_one_scale(
      scale_id = row$scale_id,
      data_path = row$file,
      k = as.integer(row$n_items),
      n_pairs = as.integer(row$n_pairs)
    ),
    error = function(e) {
      data.frame(
        scale_id = row$scale_id,
        k_items = row$n_items,
        n_pairs_manifest = row$n_pairs,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
}

# rbind with non-identical columns (errors vs successes)
all_names <- unique(unlist(lapply(results, names)))
results_aligned <- lapply(results, function(x) {
  missing <- setdiff(all_names, names(x))
  for (m in missing) x[[m]] <- NA
  x[, all_names, drop = FALSE]
})
final <- do.call(rbind, results_aligned)
# The row constructor mixes text and numeric fields; recover native column types.
final[] <- lapply(final, function(x) type.convert(x, as.is = TRUE))

# Preserve any optional manifest metadata (for example clean/harmonized status)
# without requiring a particular column name.  Prefix to avoid collisions.
extra_manifest_cols <- setdiff(names(manifest), required_manifest_cols)
if (length(extra_manifest_cols) > 0) {
  idx <- match(final$scale_id, manifest$scale_id)
  for (nm in extra_manifest_cols) {
    final[[paste0("manifest_", nm)]] <- manifest[[nm]][idx]
  }
}

write.csv(final, file.path(out_dir, "invariance_results.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Manuscript-oriented RQ2 table: threshold and loading invariance are reported
# from DIFFERENT identification schemes, rather than merged into a hybrid rate.
# ---------------------------------------------------------------------------
report_cols <- c(
  "scale_id", "k_items", "n_pairs_manifest",
  "wu_configural_converged", "wu_configural_heywood",
  "wu_thresholds_converged", "wu_thresholds_heywood",
  "wu_threshold_status", "wu_threshold_delta_cfi", "wu_threshold_delta_rmsea", "wu_threshold_delta_srmr",
  "ms_baseline_converged", "ms_baseline_heywood",
  "ms_loading_converged", "ms_loading_heywood",
  "ms_loading_status", "ms_loading_delta_cfi", "ms_loading_delta_rmsea", "ms_loading_delta_srmr",
  "ms_threshold_after_loading_status", "ms_threshold_after_loading_delta_cfi", "ms_threshold_after_loading_delta_rmsea",
  "wu_latent_cor_oriented", "ms_latent_cor_oriented",
  "wu_tucker_phi", "ms_tucker_phi",
  "wu_cr_h", "wu_cr_t", "wu_ave_h", "wu_ave_t",
  grep("^manifest_", names(final), value = TRUE)
)
report_cols <- report_cols[report_cols %in% names(final)]
rq2_table <- final[, report_cols, drop = FALSE]
write.csv(rq2_table, file.path(out_dir, "invariance_reporting_table.csv"), row.names = FALSE)

# Compact counts. Denominator = only steps that were actually interpretable
# (supported or rejected), not nonconverged/inadmissible/not-testable models.
summarize_status <- function(x, label, role) {
  x <- as.character(x)
  data.frame(
    analysis = label,
    role = role,
    n_total_rows = length(x),
    n_interpretable = sum(x %in% c("supported", "rejected"), na.rm = TRUE),
    n_supported = sum(x == "supported", na.rm = TRUE),
    n_rejected = sum(x == "rejected", na.rm = TRUE),
    n_inadmissible_heywood = sum(x == "inadmissible_heywood", na.rm = TRUE),
    n_not_testable_nonconvergence = sum(x == "not_testable_nonconvergence", na.rm = TRUE),
    n_not_testable_df = sum(x == "not_testable_df", na.rm = TRUE),
    n_not_testable_fit_indices = sum(x == "not_testable_fit_indices", na.rm = TRUE),
    n_not_interpreted = sum(grepl("^not_interpreted", x), na.rm = TRUE),
    supported_pct_of_interpretable = ifelse(
      sum(x %in% c("supported", "rejected"), na.rm = TRUE) > 0,
      100 * sum(x == "supported", na.rm = TRUE) / sum(x %in% c("supported", "rejected"), na.rm = TRUE),
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
if ("wu_threshold_status" %in% names(final)) {
  summary_rows[[length(summary_rows) + 1L]] <- summarize_status(
    final$wu_threshold_status,
    "Wu W0 -> W1: threshold invariance",
    "PRIMARY threshold result"
  )
}
if ("ms_loading_status" %in% names(final)) {
  summary_rows[[length(summary_rows) + 1L]] <- summarize_status(
    final$ms_loading_status,
    "Millsap/Liu M0 -> M1: loading invariance",
    "PRIMARY loading result"
  )
}
if ("ms_threshold_after_loading_status" %in% names(final)) {
  summary_rows[[length(summary_rows) + 1L]] <- summarize_status(
    final$ms_threshold_after_loading_status,
    "Millsap/Liu M1 -> M2: thresholds after loading invariance",
    "CONDITIONAL: only after loading supported"
  )
}
if ("ms_unique_after_threshold_status" %in% names(final)) {
  summary_rows[[length(summary_rows) + 1L]] <- summarize_status(
    final$ms_unique_after_threshold_status,
    "Millsap/Liu M2 -> M3: unique/residual variances",
    "CONDITIONAL: only after loading+threshold supported"
  )
}
if ("wu_conditional_loading_status" %in% names(final)) {
  summary_rows[[length(summary_rows) + 1L]] <- summarize_status(
    final$wu_conditional_loading_status,
    "Wu W1 -> W2: loadings conditional on threshold invariance",
    "DIAGNOSTIC only; do NOT use as unconditional metric rate"
  )
}

if (length(summary_rows) > 0) {
  rq2_summary <- do.call(rbind, summary_rows)
  write.csv(rq2_summary, file.path(out_dir, "invariance_summary_counts.csv"), row.names = FALSE)
}

# Cross-scheme common-admissible subset for robustness of latent metrics.  This
# does NOT combine invariance decisions across schemes.
if (all(c("wu_configural_converged", "wu_configural_heywood",
          "ms_baseline_converged", "ms_baseline_heywood") %in% names(final))) {
  common_ok <- final$wu_configural_converged %in% TRUE & !(final$wu_configural_heywood %in% TRUE) &
               final$ms_baseline_converged %in% TRUE & !(final$ms_baseline_heywood %in% TRUE)
  cross <- final[common_ok, c(
    "scale_id", "k_items",
    "wu_latent_cor_oriented", "ms_latent_cor_oriented",
    "wu_tucker_phi", "ms_tucker_phi",
    "wu_cr_h", "ms_cr_h", "wu_cr_t", "ms_cr_t",
    "wu_ave_h", "ms_ave_h", "wu_ave_t", "ms_ave_t"
  ), drop = FALSE]
  if (nrow(cross) > 0) {
    cross$delta_latent_cor_wu_minus_ms <- cross$wu_latent_cor_oriented - cross$ms_latent_cor_oriented
    cross$delta_tucker_wu_minus_ms <- cross$wu_tucker_phi - cross$ms_tucker_phi
    write.csv(cross, file.path(out_dir, "cross_scheme_common_admissible.csv"), row.names = FALSE)
  }
}

# A compact README so that output columns are not misinterpreted.
readme <- c(
  "REPORTING-REVISED paired/dependent ordinal CFA.",
  paste0("Conventional status flags use DeltaCFI >= ", DELTA_CFI_MIN,
         " AND DeltaRMSEA <= ", DELTA_RMSEA_MAX,
         ". Raw deltas and LRT p-values are also saved; DG-DMI should be run separately if required."),
  "PRIMARY THRESHOLD RESULT: Wu & Estabrook (2016), theta parameterization, W0 configural -> W1 thresholds.",
  "Wu W1->W2 loadings are CONDITIONAL on threshold invariance and are diagnostic only; do NOT report their rate as the unconditional metric/loading-invariance rate.",
  "PRIMARY LOADING RESULT: Millsap & Tein (2004)/Liu et al. (2017), theta + marker identification, M0 baseline -> M1 loadings.",
  "Millsap/Liu thresholds are tested only after loading invariance is supported; unique/residual variances only after loading+threshold invariance is supported.",
  "The script intentionally skips more restrictive models when the preceding invariance level is rejected, inadmissible, nonconvergent, or not testable.",
  "Millsap/Liu currently uses it1 as the marker in both Human and Twin factors.",
  "Category labels must match exactly between Human and Twin for each repeated item; mismatches are stopped and reported.",
  "Do NOT merge Wu-only and Millsap-only invariance decisions into one hybrid percentage. Report Wu threshold invariance and Millsap/Liu loading invariance separately.",
  "invariance_reporting_table.csv is the manuscript-oriented table.",
  "invariance_summary_counts.csv gives transparent denominators: only supported/rejected steps count as interpretable.",
  "CR/AVE/Tucker/Fh-Ft are extracted separately from each scheme's least-restricted model; cross_scheme_common_admissible.csv can be used to assess robustness of these latent metrics across identification schemes.",
  "Do not rerun DG-DMI on old models; if needed, rerun it using the corrected fitted-model specification."
)
writeLines(readme, file.path(out_dir, "README.txt"))

cat("\nDone. Main outputs:\n")
cat("  ", file.path(out_dir, "invariance_results.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "invariance_reporting_table.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "invariance_summary_counts.csv"), "\n", sep = "")
cat("See also model syntax, summaries, parameter estimates, thresholds, category diagnostics, and README in:\n  ", out_dir, "\n", sep = "")
