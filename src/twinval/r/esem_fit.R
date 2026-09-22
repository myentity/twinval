## ESEM (exploratory structural equation model) with TARGET rotation, ordinal indicators, WLSMV.
##
## Fitted separately for the Human and the Twin responses of every multidimensional instrument, with the SAME
## theoretical target and the SAME list of factors on both sides (rotation and number of factors are not tuned
## towards a better Human-Twin similarity).
##
## Usage: Rscript esem_fit.R <workdir>
## Input  <workdir>/esem_registry.csv               columns: instrument_id, instrument_name
##        <workdir>/item_battery/<id>_item_battery.csv   columns: id, <item>_h, <item>_t
##        <workdir>/target_matrices/<id>_target_matrix.csv  items x factors (0/1), first column = item
## Output <workdir>/esem_model_diagnostics.csv, esem_fit_by_instrument.csv, esem_loadings_long.csv,
##        esem_thresholds.csv, esem_factor_correlations.csv, esem_loadings_matrices/<id>_<side>_loadings.csv

suppressPackageStartupMessages(library(lavaan))

args <- commandArgs(trailingOnly = TRUE)
wd <- if (length(args) >= 1) args[[1]] else "."
reg <- read.csv(file.path(wd, "esem_registry.csv"), stringsAsFactors = FALSE)
dir.create(file.path(wd, "esem_loadings_matrices"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(wd, "esem_summaries"), showWarnings = FALSE, recursive = TRUE)

fit_one_side <- function(battery, target_mat, items, nfactors) {
  d <- battery[items]
  for (col in items) d[[col]] <- ordered(d[[col]])
  warn_msgs <- character(0)
  res <- withCallingHandlers(
    tryCatch({
      fit <- efa(data = d, nfactors = nfactors, rotation = list("target", target = target_mat),
                 ordered = items, estimator = "WLSMV", output = "lavaan")
      if (!isTRUE(lavInspect(fit, "converged"))) {
        list(status = "nonconverged", fit = fit, note = "converged=FALSE")
      } else {
        theta <- lavInspect(fit, "theta")
        psi <- lavInspect(fit, "cov.lv")
        heywood_resid <- any(diag(theta) < 0, na.rm = TRUE)
        heywood_corr <- any(abs(psi[lower.tri(psi)]) > 1.0001, na.rm = TRUE)
        status <- if (heywood_resid || heywood_corr) "inadmissible_heywood" else "ok"
        list(status = status, fit = fit,
             note = paste0("heywood_resid=", heywood_resid, "; heywood_corr=", heywood_corr))
      }
    }, error = function(e) list(status = "error", fit = NULL, note = conditionMessage(e))),
    warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  if (length(warn_msgs) > 0) res$note <- paste0(res$note, " | warnings: ", paste(unique(warn_msgs), collapse = "; "))
  res
}

fit_rows <- list(); load_rows_all <- list(); fc_rows <- list(); thr_rows <- list(); diag_rows <- list()

for (i in seq_len(nrow(reg))) {
  tid <- reg$instrument_id[i]; tname <- reg$instrument_name[i]
  cat(sprintf("\n=== [%d/%d] %s (instrument_id=%d) ===\n", i, nrow(reg), tname, tid))
  battery <- read.csv(file.path(wd, "item_battery", sprintf("%d_item_battery.csv", tid)), stringsAsFactors = FALSE, check.names = FALSE)
  target_df <- read.csv(file.path(wd, "target_matrices", sprintf("%d_target_matrix.csv", tid)), row.names = 1, check.names = FALSE)
  factor_names <- colnames(target_df); nfac <- length(factor_names)
  h_items <- grep("_h$", names(battery), value = TRUE)
  t_items <- grep("_t$", names(battery), value = TRUE)
  stopifnot(length(h_items) == length(t_items), length(h_items) == nrow(target_df))
  tgt_h <- as.matrix(target_df[sub("_h$", "", h_items), , drop = FALSE]); rownames(tgt_h) <- h_items
  tgt_t <- as.matrix(target_df[sub("_t$", "", t_items), , drop = FALSE]); rownames(tgt_t) <- t_items

  for (side_info in list(list(side = "human", items = h_items, target = tgt_h),
                          list(side = "twin", items = t_items, target = tgt_t))) {
    side <- side_info$side; items <- side_info$items; target_mat <- side_info$target
    t0 <- Sys.time()
    res <- fit_one_side(battery, target_mat, items, nfac)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  %s: status=%s, %.1fs, %s\n", side, res$status, elapsed, substr(res$note, 1, 120)))
    diag_rows[[length(diag_rows) + 1]] <- data.frame(instrument_id = tid, instrument_name = tname, side = side,
      status = res$status, note = res$note, elapsed_sec = elapsed, stringsAsFactors = FALSE)
    if (is.null(res$fit)) next
    fit <- res$fit
    fm <- tryCatch(fitMeasures(fit, c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled", "tli.scaled",
                                      "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr")),
                   error = function(e) rep(NA, 9))
    fit_rows[[length(fit_rows) + 1]] <- data.frame(instrument_id = tid, instrument_name = tname, side = side,
      status = res$status, n_factors = nfac, n_items = length(items), chisq = fm[1], df = fm[2], pvalue = fm[3],
      cfi = fm[4], tli = fm[5], rmsea = fm[6], rmsea_ci_lo = fm[7], rmsea_ci_hi = fm[8], srmr = fm[9],
      stringsAsFactors = FALSE)
    ss <- tryCatch(standardizedSolution(fit), error = function(e) NULL)
    if (is.null(ss)) next
    lr <- ss[ss$op == "=~", c("lhs", "rhs", "est.std", "se", "pvalue")]
    fac_map <- setNames(factor_names, paste0("f", seq_along(factor_names)))
    lr$factor_name <- fac_map[lr$lhs]
    lr$is_target <- mapply(function(it, fac) target_mat[it, fac], lr$rhs, lr$factor_name)
    load_rows_all[[length(load_rows_all) + 1]] <- data.frame(instrument_id = tid, instrument_name = tname, side = side,
      factor = lr$factor_name, item = lr$rhs, loading = lr$est.std, se = lr$se, pvalue = lr$pvalue,
      is_target = lr$is_target, stringsAsFactors = FALSE)
    fc <- ss[ss$op == "~~" & ss$lhs != ss$rhs & grepl("^f[0-9]+$", ss$lhs) & grepl("^f[0-9]+$", ss$rhs), ]
    if (nrow(fc) > 0) fc_rows[[length(fc_rows) + 1]] <- data.frame(instrument_id = tid, instrument_name = tname,
      side = side, factor1 = fac_map[fc$lhs], factor2 = fac_map[fc$rhs], correlation = fc$est.std, se = fc$se,
      pvalue = fc$pvalue, stringsAsFactors = FALSE)
    thr <- ss[ss$op == "|", c("lhs", "rhs", "est.std", "se")]
    if (nrow(thr) > 0) thr_rows[[length(thr_rows) + 1]] <- data.frame(instrument_id = tid, instrument_name = tname,
      side = side, item = thr$lhs, threshold = thr$rhs, value = thr$est.std, se = thr$se, stringsAsFactors = FALSE)
    suffix <- paste0("_", substr(side, 1, 1), "$")
    lm <- matrix(NA_real_, nrow = length(items), ncol = nfac, dimnames = list(sub(suffix, "", items), factor_names))
    for (r in seq_len(nrow(lr))) lm[sub(suffix, "", lr$rhs[r]), lr$factor_name[r]] <- lr$est.std[r]
    write.csv(lm, file.path(wd, "esem_loadings_matrices", sprintf("%d_%s_loadings.csv", tid, side)))
    tryCatch({ sink(file.path(wd, "esem_summaries", sprintf("%d_%s_summary.txt", tid, side)))
               print(summary(fit, fit.measures = TRUE, standardized = TRUE)) },
             error = function(e) cat("summary() failed:", conditionMessage(e), "\n"), finally = sink())
  }
}
w <- function(x, f) write.csv(if (length(x)) do.call(rbind, x) else data.frame(), file.path(wd, f), row.names = FALSE)
w(diag_rows, "esem_model_diagnostics.csv"); w(fit_rows, "esem_fit_by_instrument.csv")
w(load_rows_all, "esem_loadings_long.csv"); w(fc_rows, "esem_factor_correlations.csv"); w(thr_rows, "esem_thresholds.csv")
cat("\nDone.\n")
