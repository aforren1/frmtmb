# Reviewer, item 2.4: re-run a SUBSET of the lane's 200 replicates on
# the lane's own seeds, and check the numbers land where the tables say.
#
# Three seeds, chosen so that each answers something different:
#   20260910  the first replicate, which the new test-lca.R block uses:
#             the identity numbers.
#   20260970  the one of the eight where nlminb returns 0 and ?lca
#             claims diagnose() prints "No convergence problems
#             detected": the claim that a silent wrong answer is silent.
#   20261013  the one of the eight where the JITTER recipe found nothing
#             better in ten refits while the DOC recipe found the
#             optimum on all ten.
#
# The doc recipe is run exactly as ?lca prints it, `$` reads and all,
# because the help page's text is the deliverable of a Phase 2 row.
#
#   Rscript dev/rev-latent-lca.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/latent-lca-sim.R")
suppressMessages(loadNamespace("poLCA"))
hr <- function(s) cat("\n======== ", s, " ========\n", sep = "")

K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
SEEDS <- c(20260910L, 20260970L, 20261013L)

for (seed in SEEDS) {
  hr(paste("seed", seed))
  s <- lca_sim(seed = seed)
  dp <- data.frame(x1 = s$dd$x1, x2 = s$dd$x2)
  for (j in seq_len(ncol(s$Y))) dp[[paste0("I", j)]] <- s$Y[, j]
  ff <- stats::as.formula(paste0("cbind(",
    paste(paste0("I", seq_len(ncol(s$Y))), collapse = ", "),
    ") ~ x1 + x2"))

  t0 <- Sys.time()
  fit <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd, control = tight)))
  sfrm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ll <- as.numeric(logLik(fit))
  d1 <- frmtmb::diagnose(fit, quiet = TRUE)
  msg <- utils::capture.output(frmtmb::diagnose(fit))

  set.seed(700000L + (seed %% 1000L))
  t1 <- Sys.time()
  pl <- poLCA::poLCA(ff, dp, nclass = K, nrep = 10L, verbose = FALSE,
                     maxiter = 20000, tol = 1e-12)
  spl <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  cat("frm  logLik:", format(ll, digits = 14), sprintf(" (%.1f s)\n", sfrm))
  cat("poLCA nrep=10:", format(pl$llik, digits = 14),
      sprintf(" (%.1f s)\n", spl))
  cat("gap          :", format(ll - pl$llik, digits = 8), "\n")
  cat("|d|/|ll|     :", format(abs(ll - pl$llik) / abs(ll), digits = 3), "\n")
  cat("nlminb code  :", d1$convergence, "  pdHess:", isTRUE(d1$pdHess),
      "  bad_se:", length(d1$bad_se), "  flat:", length(d1$flat), "\n")
  cat("max_grad/|ll|:", format(d1$max_grad / abs(ll), digits = 4), "\n")
  cat("diagnose(), verbatim:\n")
  cat(paste0("  | ", msg, collapse = "\n"), "\n")
  cat("prints 'No convergence problems detected':",
      any(grepl("No convergence problems detected", msg)), "\n")

  ## ---- the recipe ?lca prints, run exactly as printed --------------
  cat("\n-- the ?lca recipe, verbatim --\n")
  p0 <- fit$frame$par_template[fit$frame$extra_names]
  cat("  extra_names           :",
      paste(fit$frame$extra_names, collapse = " "), "\n")
  cat("  p0 components         :", paste(names(p0), collapse = " "), "\n")
  cat("  p0 lengths            :", paste(vapply(p0, length, 1L),
                                         collapse = " "), "\n")
  cat("  p0 is the ITEM logits, not the gating coefficients:",
      !any(c("beta", "betad") %in% names(p0)), "\n")
  set.seed(500000L + (seed %% 1000L))
  t2 <- Sys.time()
  refits <- replicate(10, simplify = FALSE, {
    st <- lapply(p0, function(v) v + stats::rnorm(length(v)))
    try(suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
          control = tight, start = st))), silent = TRUE)
  })
  sdoc <- as.numeric(difftime(Sys.time(), t2, units = "secs"))
  lls <- vapply(refits, function(f) {
    if (inherits(f, "try-error")) NA_real_ else as.numeric(logLik(f))
  }, 1)
  cat("  ten refit logLiks     :\n")
  cat(paste0("    ", format(lls, digits = 12), collapse = "\n"), "\n")
  best <- max(c(ll, lls), na.rm = TRUE)
  cat("  best over the ten     :", format(best, digits = 14), "\n")
  cat("  gap to poLCA(nrep=10) :", format(best - pl$llik, digits = 4), "\n")
  cat("  reached poLCA's optimum:",
      abs(best - pl$llik) < 1e-6 * abs(pl$llik), "\n")
  cat("  refits better than cold:", sum(lls > ll + 1e-9 * abs(ll),
                                        na.rm = TRUE), "of 10\n")
  cat("  errors                 :", sum(is.na(lls)), "\n")
  cat("  seconds for ten refits :", format(sdoc, digits = 4),
      "  (?lca says 'about 5 seconds')\n")

  flush(stdout())
}
