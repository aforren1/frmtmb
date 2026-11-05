# Item 3: mutate the fix and see whether the lane's tests fail. Each mutant is a
# copy of the source under rev-impsamp-src with one edit, installed into its own
# directory of rev-impsamp-lib. Never touches the lane's tree or library.
# Usage: Rscript dev/review-impsamp-mutants.R make <mutant>
#        Rscript dev/review-impsamp-mutants.R test <mutant> <file.R>
WT <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp"
SRC <- "C:/Users/adf44/source/r/rev-impsamp-src"
LIB <- "C:/Users/adf44/source/r/rev-impsamp-lib"
args <- commandArgs(TRUE)
cmd <- args[[1]]
mut <- args[[2]]
muts <- list(
  # the partial derivative of the likelihood's curvature, not the total
  partial_rate = list(
    from = "    up[rnd] <- up[rnd] + du\n    dn[rnd] <- dn[rnd] - du\n",
    to = ""),
  # no scale at all: the curvature frozen, the mode and prior still follow
  no_rate = list(
    from = "    hscale <- exp(as.vector(hlik_rate %*% dth))\n",
    to = "    hscale <- 1\n"),
  # the mode does not take the likelihood's slope
  no_cross = list(
    from = "    lik <- as.vector(cross_d %*% dth)\n",
    to = "    lik <- as.vector(0 * (cross_d %*% dth))\n"),
  # the whole proposal frozen where it was anchored: the old estimator
  frozen = list(
    from = "    list(u = u, logdet = logdet, prior = tm)\n",
    to = paste0(
      "    list(u = lapply(seq_len(q), function(j) as.vector(\n",
      "      plan[[\"u\"]][plan[[\"idx\"]][j, ], , drop = FALSE])),\n",
      "      logdet = plan[[\"logdet_l\"]], prior = tm)\n")),
  # the Jacobian frozen while the draws move: a wrong density off the anchor
  no_jacobian = list(
    from = "    list(u = u, logdet = logdet, prior = tm)\n",
    to = "    list(u = u, logdet = plan[[\"logdet_l\"]], prior = tm)\n"),
  # half the Jacobian: wrong by a constant factor of the log-determinant move
  half_jacobian = list(
    from = "    list(u = u, logdet = logdet, prior = tm)\n",
    to = paste0("    list(u = u, logdet = 0.5 * (logdet + ",
                "plan[[\"logdet_l\"]]), prior = tm)\n")),
  # the old imp_par_list: recycle over mapped entries
  map_recycle = list(
    from = "      out[[cp]][keep] <- unname(par[pn == cp])[fm[keep]]\n",
    to = "      out[[cp]][] <- unname(par[pn == cp])\n"),
  # the prior's move left out of the mode
  no_prior_mode = list(
    from = "        s <- s + (pjk - prec_a[j, k]) * uhat_g[[k]]\n",
    to = "        s <- s + 0 * (pjk - prec_a[j, k]) * uhat_g[[k]]\n"),
  # the refusal never fires
  no_refusal = list(
    from = "  if (it$rose) {\n    imp_refuse_draws(it, n_draw, iterate, lay, frame[[\"n_obs\"]])\n  }\n",
    to = ""),
  # one round, whatever the fresh gradient says: "convergence by construction"
  one_round = list(
    from = "      if (moved < imp_move_tol) break\n",
    to = "      break\n"),
  # the stationarity criterion always true
  grad_always = list(
    from = "      if (is.finite(grad) && grad < control$grad_tol) break\n",
    to = "      break\n"),
  # NOT a mutation of the test power: the proposed FIX for the fixed step.
  # A step scale-free in the coordinate, from the diagonal of the joint Hessian
  fix_h = list(
    from = "    h <- imp_cross_step * max(1, abs(full[[fix[[j]]]]))
",
    to = paste0(
      "    h <- imp_cross_step * max(1, abs(full[[fix[[j]]]]))
",
      "    if (is.finite(djoint[[j]]) && djoint[[j]] > 0) {
",
      "      h <- min(h, imp_cross_step / sqrt(djoint[[j]]))
",
      "    }
")),
  # no cut of a negative likelihood curvature on q = 1
  no_cut = list(
    from = "    hlik[1L, 1L, ] <- pmax(hd - prec[1L, 1L], 0)\n",
    to = "    hlik[1L, 1L, ] <- abs(hd - prec[1L, 1L])\n")
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
  s <- paste(readLines(p), collapse = "\n")
  s <- paste0(s, "\n")
  if (lengths(gregexpr(m$from, s, fixed = TRUE)) != 1L ||
      !grepl(m$from, s, fixed = TRUE)) {
    stop("mutation site not found once: ", mut)
  }
  s <- sub(m$from, m$to, s, fixed = TRUE)
  writeBin(charToRaw(s), p)
  st <- system2(file.path(R.home("bin"), "R.exe"),
                c("CMD", "INSTALL", paste0("--library=", shQuote(lib)),
                  "--no-multiarch", "--no-docs", "--no-test-load",
                  shQuote(dst)),
                stdout = TRUE, stderr = TRUE)
  cat(tail(st, 3), sep = "\n")
} else {
  f <- file.path(WT, "tests/testthat", args[[3]])
  Sys.setenv(NOT_CRAN = "true")
  .libPaths(c(file.path(LIB, mut), "C:/Users/adf44/source/r/wt-impsamp-lib",
              "C:/Users/adf44/source/r/rellib-r16",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages({
    library(frmtmb)
    library(testthat)
  })
  cat("lib:", find.package("frmtmb"), "\n")
  res <- as.data.frame(test_file(f, package = "frmtmb",
                                 env = test_env("frmtmb"), reporter = "silent",
                                 stop_on_failure = FALSE))
  bad <- res[res$failed > 0 | res$error, c("test", "failed", "error")]
  cat(sprintf("RESULT %s %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
              mut, basename(f), sum(res$passed), sum(res$failed),
              sum(res$error), sum(res$skipped), sum(res$warning)))
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  FAILED: %s (fail %d, error %s)\n", bad$test[i], bad$failed[i],
                bad$error[i]))
  }
}
