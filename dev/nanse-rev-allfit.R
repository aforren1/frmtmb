# Reviewer: frm_allfit() verdicts on ordinary models. Is grad_tol (a
# gradient tolerance, 1e-3 by default) a sound yardstick for "converged
# elsewhere"? Prints each table and, per refit below the best, its own
# convergence verdict and its distance from the best estimates.
#   Rscript dev/nanse-rev-allfit.R
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(1)
n <- 2000
dg <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(1:50, n, TRUE)))
dg$y <- rpois(n, exp(0.2 + 0.3 * dg$x + rnorm(50, 0, 0.5)[dg$g]))
set.seed(11)
d2 <- data.frame(f = factor(sample(1:40, 3000, TRUE)), x = rnorm(3000),
                 g = factor(sample(1:30, 3000, TRUE)))
d2$y <- rpois(3000, exp(0.3 + 0.2 * d2$x + rnorm(40, 0, 0.3)[d2$f] +
                          rnorm(30, 0, 0.5)[d2$g]))
sleep <- lme4::sleepstudy
models <- list(
  pois_glmm_2000 = function() frm(y ~ x + z + (1 | g), data = dg,
                                  family = poisson()),
  glmm_f40 = function() frm(y ~ x + f + (1 | g), data = d2,
                            family = poisson()),
  lmm_slope = function() frm(Reaction ~ Days + (Days | Subject),
                             data = sleep),
  negbin = function() frm(count ~ zAge + zBase * Trt + (1 | patient),
                          data = brms::epilepsy, family = negbinomial())
)
for (nm in names(models)) {
  fit <- suppressWarnings(models[[nm]]())
  for (sfm in c(TRUE, FALSE)) {
    cat("\n==", nm, "start_from_mle =", sfm, "\n")
    a <- frm_allfit(fit, start_from_mle = sfm)
    print(a)
    ll <- vapply(a$fits, function(f) if (is.null(f)) NA_real_ else
      as.numeric(logLik(f)), 0)
    best <- max(c(ll, a$original$logLik), na.rm = TRUE)
    for (k in which(ll < best - a$tol)) {
      f <- a$fits[[k]]
      chk <- suppressWarnings(ns$check_convergence(f, f$control))
      cat(sprintf("  %s: %.3g below; own check: %s; max |par - best| %.3g\n",
                  names(a$fits)[k], best - ll[k],
                  if (length(chk$warnings)) substr(chk$warnings[1], 1, 70)
                  else "passes",
                  max(abs(f$opt$par - fit$opt$par))))
    }
  }
}
