# Mutate the fix and run tests on the mutant. Sites are for the final source. Each mutant is a source copy
# outside the worktree, installed into its own library outside the shared ones.
#   Rscript dev/impsamp-mutants.R make <mutant>
#   Rscript dev/impsamp-mutants.R test <mutant> <regexp for test names or "">
WT <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp"
SRC <- "C:/Users/adf44/source/r/wt-impsamp-mutsrc"
LIB <- "C:/Users/adf44/source/r/wt-impsamp-mutlib"
args <- commandArgs(TRUE)
cmd <- args[[1]]
mut <- args[[2]]
muts <- list(
  # fix round 2. A parameter with no column steps by its own size (F1)
  own_size = list(
    c("  h <- rep(imp_cross_step, length(at))\n",
      "  h <- imp_cross_step * pmax(1, abs(unname(at)))\n")),
  # F2: the precision read at the unit vector
  t_one = list(
    c("      if (is.finite(v) && v > 0) t0 / sqrt(v) else 1\n", "      1\n")),
  # F3: the rounds start where the first optimizer stops
  no_escape = list(
    c("  if (!is.null(lap_esc)) lap_opt <- lap_esc\n", "")),
  # F4: a settled count is named whatever its worst group reads
  low_ess_named = list(
    c("    if (next_try$settled && min(next_try$ess$ess) >= ess_floor) {\n",
      "    if (next_try$settled) {\n")),
  # F6: no response unit in the step
  no_units = list(
    c("    u <- if (key %in% names(units)) units[[key]] else 1\n", "    u <- 1\n")),
  # the curvature's partial derivative: the mode does not move with the step
  partial_rate = list(
    c("      up[rnd] <- up[rnd] - hsolve(d) * h\n", ""),
    c("      up[rnd] <- up[rnd] + du\n      dn[rnd] <- dn[rnd] - du\n", "")),
  # the old step: absolute, in the coefficient's own units
  abs_step = list(
    c("  steps <- imp_outer_steps(frame, full[fix])\n",
      "  steps <- imp_cross_step * pmax(1, abs(unname(full[fix])))\n")),
  # the column's scale for beta only, not for a dispersion's column
  beta_only = list(
    c("    if (is.null(scale[[cp]]) || is.null(X) || !ncol(X)) next\n",
      "    if (is.null(scale[[cp]]) || is.null(X) || !ncol(X) || cp != \"beta\") next\n")),
  no_rate = list(
    c("    hscale <- exp(as.vector(hlik_rate %*% dth))\n", "    hscale <- 1\n")),
  no_cross = list(
    c("    lik <- as.vector(cross_d %*% dth)\n",
      "    lik <- as.vector(0 * (cross_d %*% dth))\n")),
  frozen = list(
    c("    list(u = u, logdet = logdet, prior = tm)\n",
      paste0("    list(u = lapply(seq_len(q), function(j) as.vector(\n",
             "      plan[[\"u\"]][plan[[\"idx\"]][j, ], , drop = FALSE])),\n",
             "      logdet = plan[[\"logdet_l\"]], prior = tm)\n"))),
  no_jacobian = list(
    c("    list(u = u, logdet = logdet, prior = tm)\n",
      "    list(u = u, logdet = plan[[\"logdet_l\"]], prior = tm)\n")),
  half_jacobian = list(
    c("    list(u = u, logdet = logdet, prior = tm)\n",
      "    list(u = u, logdet = 0.5 * (logdet + plan[[\"logdet_l\"]]), prior = tm)\n")),
  map_recycle = list(
    c("      out[[cp]][keep] <- unname(par[pn == cp])[fm[keep]]\n",
      "      out[[cp]][] <- unname(par[pn == cp])\n")),
  no_prior_mode = list(
    c("        s <- s + (pjk - prec_a[j, k]) * uhat_g[[k]]\n",
      "        s <- s + 0 * (pjk - prec_a[j, k]) * uhat_g[[k]]\n")),
  no_refusal = list(
    c("  if (it$rose) {\n    imp_refuse_draws(it, n_draw, iterate, lay, frame[[\"n_obs\"]])\n  }\n", "")),
  one_round = list(
    c("      if (moved < imp_move_tol) break\n", "      break\n")),
  grad_always = list(
    c("      if (is.finite(grad) && grad < control$grad_tol) break\n", "      break\n")),
  no_cut = list(
    c("    hlik[1L, 1L, ] <- pmax(hd - prec[1L, 1L], 0)\n",
      "    hlik[1L, 1L, ] <- abs(hd - prec[1L, 1L])\n")),
  # a capped run counts as settled
  capped_settles = list(
    c("rose = rose, capped = capped, settled = !rose && !capped)",
      "rose = rose, capped = capped, settled = !rose)")),
  # the prior left out of the objective, and out of the pin so that it runs
  no_prior = list(
    c("         if (is.null(neg_log_prior)) val else val + neg_log_prior(pars)\n",
      "         val\n"))
)
if (cmd == "make") {
  m <- muts[[mut]]
  dst <- file.path(SRC, mut)
  lib <- file.path(LIB, mut)
  unlink(dst, recursive = TRUE)
  dir.create(file.path(dst, "R"), recursive = TRUE)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  for (f in c("DESCRIPTION", "NAMESPACE", "LICENSE", ".Rbuildignore")) {
    if (file.exists(file.path(WT, f))) file.copy(file.path(WT, f), dst)
  }
  file.copy(list.files(file.path(WT, "R"), full.names = TRUE),
            file.path(dst, "R"))
  for (d in c("src", "inst", "data", "man")) {
    if (dir.exists(file.path(WT, d))) {
      file.copy(file.path(WT, d), dst, recursive = TRUE)
    }
  }
  p <- file.path(dst, "R", "importance.R")
  s <- paste0(paste(readLines(p), collapse = "\n"), "\n")
  for (e in m) {
    if (!grepl(e[[1L]], s, fixed = TRUE)) stop("no site: ", substr(e[[1L]], 1, 50))
    s <- sub(e[[1L]], e[[2L]], s, fixed = TRUE)
  }
  writeBin(charToRaw(s), p)
  st <- system2(file.path(R.home("bin"), "R.exe"),
                c("CMD", "INSTALL", paste0("--library=", shQuote(lib)),
                  "--no-multiarch", "--no-docs", "--no-test-load", shQuote(dst)),
                stdout = TRUE, stderr = TRUE)
  cat(tail(st, 2), sep = "\n")
} else {
  f <- file.path(WT, "tests/testthat/test-importance.R")
  Sys.setenv(NOT_CRAN = "true")
  .libPaths(c(file.path(LIB, mut), "C:/Users/adf44/source/r/wt-impsamp-lib",
              "C:/Users/adf44/source/r/rellib-r16",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages({
    library(frmtmb)
    library(testthat)
  })
  cat("lib:", find.package("frmtmb"), "\n")
  stopifnot(grepl("mutlib", find.package("frmtmb"), fixed = TRUE))
  desc <- if (length(args) >= 3L && nzchar(args[[3L]])) args[[3L]]
  res <- as.data.frame(testthat::test_file(
    f, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent",
    stop_on_failure = FALSE,
    desc = desc))
  bad <- res[res$failed > 0 | res$error, c("test", "failed", "error")]
  cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
              mut, sum(res$passed), sum(res$failed), sum(res$error),
              sum(res$skipped), sum(res$warning)))
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  FAILED: %s (fail %d, error %s)\n", bad$test[i], bad$failed[i],
                bad$error[i]))
  }
}
