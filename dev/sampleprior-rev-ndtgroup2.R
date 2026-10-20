# Reviewer of lane sampleprior, fix round 1: the posterior of the ndt
# fraction under ndt_group() with the defaults, against the same defaults
# with a uniform(0, 1) on the fraction. Data seed 20261120, sampler seed 4.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample); library(frmtmb.eam)})
set.seed(20261120)
n <- 120
d <- frmtmb.eam::ddm_simulate(n, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
d$subj <- factor(rep(1:4, length.out = n))
d$x <- rnorm(n)
f <- rt | dec(upper) + ndt_group(subj) ~ x
fit <- frm(f, data = d, family = frmtmb.eam::wiener())
cat("ML ndt fraction", frm_linpred(fit, dpar = "ndt", type = "response")[1],
    " true ndt / fastest per subject", 0.3 / tapply(d$rt, d$subj, min), "\n")
for (lab in c("defaults", "uniform(0, 1) on ndt")) {
  pr <- if (lab == "defaults") NULL else
    set_prior("uniform(0, 1)", class = "ndt")
  ds <- suppressMessages(suppressWarnings(frm_sample(
    fit, prior = pr, chains = 2, iter = 600, refresh = 0, seed = 4)))
  v <- ds$draws[, "ndt"]
  cat(sprintf("%-22s ndt fraction q05/50/95 %s max %.4f\n", lab,
              paste(sprintf("%.4f", quantile(v, c(.05, .5, .95))), collapse = " "),
              max(v)))
}
