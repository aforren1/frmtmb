# The draws' column names of a model whose location predictor has no
# column (beta is empty), on the base and the lane builds: an ordinal
# y ~ 1 with a modeled disc. Seed 20261005, sampler seed 3.
# Usage: Rscript dev/ordmix-emptybeta.R [base|lane]
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "base"
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(20261005)
n <- 200
d <- data.frame(z = rnorm(n))
u <- stats::rlogis(n) / exp(0.4 * d$z)
d$y <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
fit <- frm(bf(y ~ 1, disc ~ 0 + z), family = cumulative(), data = d,
           prior = set_prior("student_t(3, 0, 2.5)", class = "Intercept"))
cat("fit variables:", variables(fit), "\n")
cat("labels:", frmtmb:::brms_par_labels(fit), "\n")
ds <- suppressMessages(suppressWarnings(
  frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3)))
cat("draws variables:", variables(ds), "\n")
