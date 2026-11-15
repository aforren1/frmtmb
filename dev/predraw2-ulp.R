# Lane predraw2: the new tests' assertions under rounding. x perturbed by
# (1 + 2e-15 k), k = 0..20, on the fixtures of
# tests/testthat/test-predict-draw-check.R: the mixture's verdict, excess
# and the ratio of its median drop at predict()'s draws to that at the
# normal law's (asserted below 1/2), and the poisson GLM's excess
# (asserted below 1/10 of the threshold).
# Usage: Rscript dev/predraw2-ulp.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
for (k in 0:20) {
  set.seed(101)
  x <- stats::rnorm(300)
  invisible(stats::runif(300))
  y <- stats::rpois(300, exp(1 + 0.3 * x))
  d <- data.frame(x = x * (1 + 2e-15 * k), y = y)
  fit <- frm(bf(y ~ x), family = mixture(poisson(), poisson()), data = d)
  set.seed(1)
  dr <- suppressWarnings(suppressMessages(ns$predict_par_drawer(fit, TRUE,
                                                                1000L)))
  f0 <- fit$obj$fn(fit$opt$par)
  drops <- vapply(1:1000, function(s) {
    fit$obj$fn(ns$fit_outer_vector(dr(s))) - f0
  }, 0)
  V <- vcov(fit, full = TRUE)
  set.seed(2)
  Z <- matrix(stats::rnorm(nrow(V) * 1000), nrow(V), 1000)
  D <- fit$opt$par + crossprod(chol(V), Z)
  f0 <- fit$obj$fn(fit$opt$par)
  raw <- apply(D, 2L, function(v) fit$obj$fn(v) - f0)
  chk <- fit$cache[["draw_check"]]
  set.seed(3)
  dg <- data.frame(x = stats::rnorm(100))
  dg$y <- stats::rpois(100, exp(0.5 + 0.4 * dg$x))
  dg$x <- dg$x * (1 + 2e-15 * k)
  g <- frm(bf(y ~ x), family = poisson(), data = dg)
  invisible(predict(g, ndraws = 5))
  cat(sprintf(paste("k %2d mixture fire %s excess %.4f n %d drop ratio %.4f",
                    "| glm excess %.5f\n"), k, chk$fire, chk$excess, chk$n,
              stats::median(drops) / stats::median(raw),
              g$cache[["draw_check"]]$excess))
}
