## Person-level bootstrap of Tucker's congruence coefficient (phi) between Human and Twin ESEM solutions.
##
## In every replicate the N respondent pairs are resampled with replacement; for each instrument the Human and the
## Twin ESEM are re-estimated (same target rotation and specification as the point estimate), the Twin loadings are
## orthogonally Procrustes-rotated onto the Human loadings and phi is computed for every factor. Replicates in which
## either side is inadmissible (non-convergence, Heywood case) are NOT replaced by the point estimate: they are
## recorded as invalid and the share of valid replicates is reported (CIs with < 80% valid replicates are unstable).
##
## Resumable and reproducible: replicate b uses set.seed(seed + b), so a run can be stopped and continued.
##
## Usage: Rscript tucker_bootstrap.R <workdir> <B> <seed>
## Input  <workdir>/esem_registry.csv, item_battery/, target_matrices/ (as for esem_fit.R)
## Output <workdir>/tucker_bootstrap_replicates.csv, tucker_bootstrap_progress_log.csv

suppressMessages(library(lavaan))

args <- commandArgs(trailingOnly = TRUE)
wd <- if (length(args) >= 1) args[[1]] else "."
B <- if (length(args) >= 2) as.integer(args[[2]]) else 500L
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 20260915L

reg <- read.csv(file.path(wd, "esem_registry.csv"), stringsAsFactors = FALSE)
n_instr <- nrow(reg)
batteries <- list(); targets <- list(); h_items <- list(); t_items <- list()
for (i in seq_len(n_instr)) {
  key <- as.character(reg$instrument_id[i])
  batteries[[key]] <- read.csv(file.path(wd, "item_battery", sprintf("%s_item_battery.csv", key)), stringsAsFactors = FALSE, check.names = FALSE)
  targets[[key]] <- as.matrix(read.csv(file.path(wd, "target_matrices", sprintf("%s_target_matrix.csv", key)), row.names = 1, check.names = FALSE))
  h_items[[key]] <- grep("_h$", names(batteries[[key]]), value = TRUE)
  t_items[[key]] <- grep("_t$", names(batteries[[key]]), value = TRUE)
}
N <- nrow(batteries[[1]])

tucker_phi <- function(a, b) {
  den <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (is.na(den) || den <= 0) return(NA_real_)
  sum(a * b) / den
}

fit_side <- function(bb, items, target_mat) {
  d <- bb[items]
  for (col in items) d[[col]] <- ordered(d[[col]])
  nf <- ncol(target_mat)
  withCallingHandlers(
    tryCatch({
      fit <- efa(data = d, nfactors = nf, rotation = list("target", target = target_mat),
                 ordered = items, estimator = "WLSMV", output = "lavaan")
      if (!isTRUE(lavInspect(fit, "converged"))) return(list(status = "nonconverged", loadings = NULL))
      theta <- lavInspect(fit, "theta"); psi <- lavInspect(fit, "cov.lv")
      if (any(diag(theta) < 0, na.rm = TRUE) || any(abs(psi[lower.tri(psi)]) > 1.0001, na.rm = TRUE))
        return(list(status = "inadmissible_heywood", loadings = NULL))
      ss <- tryCatch(standardizedSolution(fit), error = function(e) NULL)
      if (is.null(ss)) return(list(status = "extract_error", loadings = NULL))
      lr <- ss[ss$op == "=~", c("lhs", "rhs", "est.std")]
      fac_map <- setNames(colnames(target_mat), paste0("f", seq_len(nf)))
      sfx <- "_[ht]$"
      lm <- matrix(NA_real_, nrow = length(items), ncol = nf, dimnames = list(sub(sfx, "", items), colnames(target_mat)))
      for (r in seq_len(nrow(lr))) lm[sub(sfx, "", lr$rhs[r]), fac_map[lr$lhs[r]]] <- lr$est.std[r]
      list(status = "ok", loadings = lm)
    }, error = function(e) list(status = "error", loadings = NULL)),
    warning = function(w) invokeRestart("muffleWarning"))
}

out_path <- file.path(wd, "tucker_bootstrap_replicates.csv")
prog_path <- file.path(wd, "tucker_bootstrap_progress_log.csv")
start <- 1L
if (file.exists(out_path)) {
  prev <- read.csv(out_path, stringsAsFactors = FALSE)
  if (nrow(prev)) start <- max(prev$rep) + 1L
}
progress <- if (file.exists(prog_path) && start > 1L) read.csv(prog_path, stringsAsFactors = FALSE) else data.frame()
cat(sprintf("Tucker bootstrap: %d instruments, B=%d, starting at replicate %d\n", n_instr, B, start))
t_start <- Sys.time()

for (b in seq(start, B)) {
  if (b > B) break
  set.seed(seed + b)
  idx <- sample.int(N, N, replace = TRUE)
  rows <- list()
  for (i in seq_len(n_instr)) {
    tid <- reg$instrument_id[i]; tname <- reg$instrument_name[i]; key <- as.character(tid)
    bb <- batteries[[key]][idx, ]
    tm <- targets[[key]]
    th <- tm[sub("_h$", "", h_items[[key]]), , drop = FALSE]; rownames(th) <- h_items[[key]]
    tt <- tm[sub("_t$", "", t_items[[key]]), , drop = FALSE]; rownames(tt) <- t_items[[key]]
    rh <- fit_side(bb, h_items[[key]], th); rt <- fit_side(bb, t_items[[key]], tt)
    if (rh$status == "ok" && rt$status == "ok") {
      H <- rh$loadings; T <- rt$loadings
      sv <- svd(t(T) %*% H); R <- sv$u %*% t(sv$v)
      Ta <- T %*% R
      for (f in seq_len(ncol(H))) rows[[length(rows) + 1]] <- data.frame(instrument_id = tid, instrument_name = tname,
        factor = colnames(H)[f], rep = b, phi = tucker_phi(H[, f], Ta[, f]), valid = TRUE)
    } else {
      for (f in seq_len(ncol(tm))) rows[[length(rows) + 1]] <- data.frame(instrument_id = tid, instrument_name = tname,
        factor = colnames(tm)[f], rep = b, phi = NA_real_, valid = FALSE)
    }
    progress <- rbind(progress, data.frame(rep = b, instrument_id = tid, instrument_name = tname,
                                            status_human = rh$status, status_twin = rt$status))
  }
  df <- do.call(rbind, rows)
  write.table(df, out_path, sep = ",", row.names = FALSE, col.names = !file.exists(out_path), append = file.exists(out_path))
  el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  cat(sprintf("[%s] replicate %d/%d done (%.1f min elapsed, ETA ~%.0f min)\n", format(Sys.time(), "%H:%M:%S"), b, B,
              el, el / (b - start + 1) * (B - b)))
  if (b %% 5 == 0 || b == B) write.csv(progress, prog_path, row.names = FALSE)
}
write.csv(progress, prog_path, row.names = FALSE)
cat("Bootstrap finished.\n")
