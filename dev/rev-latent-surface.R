# Reviewer, item 2.4: is it the SURFACE, or is it this implementation's
# start?
#
# The lane's claim rests on one number: poLCA's own single-start EM
# "lands low on 0 to 4 of 10 seeds on the same data". Two of the eight
# read ZERO of ten. On those two, every one of poLCA's random starts
# found the global optimum while frmtmb's deterministic start did not,
# which is a weaker statement than the one on the help page. This
# re-runs the single-start arm on the lane's own seeds for the two ends
# of that range.
#
#   Rscript dev/rev-latent-surface.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/latent-lca-sim.R")
suppressMessages(loadNamespace("poLCA"))

K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
# 20261010 the lane recorded as 0 of 10, 20260999 as 4 of 10
for (seed in c(20261010L, 20260999L)) {
  cat("\n======== seed", seed, "========\n")
  s <- lca_sim(seed = seed)
  dp <- data.frame(x1 = s$dd$x1, x2 = s$dd$x2)
  for (j in seq_len(ncol(s$Y))) dp[[paste0("I", j)]] <- s$Y[, j]
  ff <- stats::as.formula(paste0("cbind(",
    paste(paste0("I", seq_len(ncol(s$Y))), collapse = ", "),
    ") ~ x1 + x2"))
  fit <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd, control = tight)))
  ll <- as.numeric(logLik(fit))
  set.seed(700000L + (seed %% 1000L))
  pl10 <- poLCA::poLCA(ff, dp, nclass = K, nrep = 10L, verbose = FALSE,
                       maxiter = 20000, tol = 1e-12)
  cat("frm cold      :", format(ll, digits = 12), "\n")
  cat("poLCA nrep=10 :", format(pl10$llik, digits = 12), "\n")
  lo <- numeric(10)
  for (i in 1:10) {
    set.seed(800000L + i)
    p1 <- try(poLCA::poLCA(ff, dp, nclass = K, nrep = 1L, verbose = FALSE,
                           maxiter = 20000, tol = 1e-12), silent = TRUE)
    lo[i] <- if (inherits(p1, "try-error")) NA_real_ else p1$llik
  }
  low <- sum(lo < pl10$llik - 1e-6 * abs(pl10$llik), na.rm = TRUE)
  cat("poLCA nrep=1, ten seeds:\n")
  cat(paste0("  ", format(lo, digits = 12), collapse = "\n"), "\n")
  cat("landed low  :", low, "of 10\n")
  cat("landed at frm's own local optimum (within 1 unit):",
      sum(abs(lo - ll) < 1, na.rm = TRUE), "of 10\n")
  cat("worst single start:", format(min(lo, na.rm = TRUE), digits = 12),
      "  frm's cold start:", format(ll, digits = 12), "\n")
  cat("is frm's cold start worse than EVERY poLCA single start?",
      ll < min(lo, na.rm = TRUE) - 1e-6, "\n")
}
