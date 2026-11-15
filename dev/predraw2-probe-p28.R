# Lane predraw2: mixture(poisson(), poisson()) on one-component data, sweep
# generator seeds 1 and 28: the check does not fire (median excess 0.04,
# -0.04), yet predict()'s row means are far from the truth (mean absolute
# error 17.8 and 342 against 0.27 and 0.14 for fitted()). Where do the
# draws go? The excess at all 1000 draws of predict() (seed 1), and the
# largest draws of the response.
# Usage: Rscript dev/predraw2-probe-p28.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
for (seed in c(1, 28)) {
  d <- gen("mixture_poisson", "boundary", seed)
  fit <- suppressWarnings(frm(bf(y ~ x), family = mixture(poisson(),
                                                          poisson()),
                              data = d))
  V <- vcov(fit, full = TRUE)
  cat("\n== seed", seed, "\n")
  print(round(cbind(est = fit$opt$par, sd = sqrt(diag(V))), 3))
  set.seed(1)
  dr <- suppressMessages(ns$predict_par_drawer(fit, TRUE, 1000L))
  D <- vapply(1:1000, function(s) ns$fit_outer_vector(dr(s)),
              numeric(length(fit$opt$par)))
  f0 <- fit$obj$fn(fit$opt$par)
  drop <- apply(D, 2, function(v) fit$obj$fn(v) - f0)
  Q <- solve(V)
  m <- apply(D - fit$opt$par, 2, function(z) sum(z * (Q %*% z)))
  e <- drop - m / 2
  cat("excess quantiles 0.5 0.9 0.99 0.999 max:",
      format(stats::quantile(e, c(0.5, 0.9, 0.99, 0.999, 1)), digits = 3),
      "\n")
  set.seed(1)
  p <- suppressWarnings(predict(fit, ndraws = 1000, summary = FALSE))
  big <- order(-apply(p, 1, max))[1:5]
  cat("largest draws of the response, by replicate:\n")
  for (s in big) {
    cat(sprintf("  replicate %4d max y %.4g excess %.3g params %s\n", s,
                max(p[s, ]), e[s], paste(format(D[, s], digits = 3),
                                         collapse = " ")))
  }
}
