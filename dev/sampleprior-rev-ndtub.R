# Reviewer of lane sampleprior: the ndt default with brms's upper bound
# given as ub = min(y), so Stan's transform holds log(ndt) below
# log(min(y)) and the NaN wall is never reached. Same data and seed as
# dev/sampleprior-rev-shiftdiv.R.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261116)
n <- 60
d <- data.frame(x = rnorm(n))
d$y <- 0.4 + exp(rnorm(n, -0.5 + 0.2 * d$x, 0.4))
my <- min(d$y)
fit <- frm(y ~ x, data = d, family = frmtmb::brmsfamily("shifted_lognormal"))
pr <- set_prior(sprintf("uniform(0, %.17g)", my), class = "ndt", ub = my)
r <- tryCatch({
  ds <- suppressMessages(suppressWarnings(frm_sample(
    fit, prior = pr, chains = 4, iter = 2000, refresh = 0, seed = 9)))
  v <- ds$draws[, "ndt"]
  sprintf("ndt q05/50/95 %s, max/min(y) %.6f, divergences %d of %d",
          paste(sprintf("%.4f", quantile(v, c(.05, .5, .95))), collapse = " "),
          max(v) / my, sum(rstan::get_divergent_iterations(ds$stanfit)),
          nrow(ds$draws))
}, error = function(e) paste("ERROR", conditionMessage(e)))
cat("fit route, uniform(0, min_Y) with ub = min_Y:", r, "\n")
