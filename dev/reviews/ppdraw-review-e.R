# Review E: brms 2.23.0 posterior_predict against frmtmb predict() on
# (1) the 3269 fit (held zi: is the interval too narrow on the identified
# side?) and (2) review B's seed 3 (zi interior at -7.0, not held).
# Usage: Rscript ppdraw-review-e.R <lib>
args <- commandArgs(TRUE)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(brms); library(frmtmb)})
cat("brms", as.character(packageVersion("brms")), " frmtmb",
    find.package("frmtmb"), "\n")
summ <- function(label, p) {
  lo <- apply(p, 2, stats::quantile, 0.05)
  hi <- apply(p, 2, stats::quantile, 0.95)
  lo1 <- apply(p, 2, stats::quantile, 0.01)
  hi1 <- apply(p, 2, stats::quantile, 0.99)
  cat(sprintf("%-28s mean %.4f  P(y = 0) %.4f  mean 90%% width %.4f  mean 98%% width %.4f\n",
              label, mean(p), mean(p == 0), mean(hi - lo), mean(hi1 - lo1)))
}
run <- function(label, d, fbrms, ffrm, fam_b, fam_f) {
  cat("\n==", label, "\n")
  b <- brms::brm(fbrms, family = fam_b, data = d, chains = 4, iter = 2000,
                 seed = 11, refresh = 0)
  zi <- brms::as_draws_df(b)$zi
  cat("brms zi quantiles (0.5, 0.9, 0.99):",
      format(stats::quantile(zi, c(0.5, 0.9, 0.99)), digits = 4), "\n")
  set.seed(1)
  pb <- brms::posterior_predict(b)
  summ("brms posterior_predict", pb)
  cat(sprintf("brms posterior_epred mean %.4f\n",
              mean(brms::posterior_epred(b))))
  fit <- suppressWarnings(frmtmb::frm(ffrm, family = fam_f, data = d))
  set.seed(1)
  pf <- withCallingHandlers(
    predict(fit, ndraws = 4000, summary = FALSE),
    message = function(m) {
      cat("  M:", substr(conditionMessage(m), 1, 90), "\n")
      invokeRestart("muffleMessage")
    })
  summ("frmtmb predict (lane)", pf)
  set.seed(1)
  p0 <- predict(fit, ndraws = 4000, summary = FALSE, propagate_error = FALSE)
  summ("frmtmb propagate = FALSE", p0)
  cat(sprintf("frmtmb fitted mean %.4f\n", mean(fitted(fit)[, "Estimate"])))
}
set.seed(7)
dz <- data.frame(x = rnorm(200), N = sample(5:15, 200, TRUE))
dz$y <- rbinom(200, dz$N, plogis(0.3 + 0.5 * dz$x))
run("3269 zib, zi held", dz, brms::bf(y | trials(N) ~ x),
    frmtmb::bf(y | trials(N) ~ x), brms::zero_inflated_binomial(),
    frmtmb::zero_inflated_binomial())
set.seed(3)
dp <- data.frame(x = stats::rnorm(200))
dp$y <- stats::rpois(200, exp(1 + 0.3 * dp$x))
run("zip null data seed 3, zi -7.0 not held", dp, brms::bf(y ~ x),
    frmtmb::bf(y ~ x), brms::zero_inflated_poisson(),
    frmtmb::zero_inflated_poisson())
