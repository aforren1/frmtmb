# Reviewer of lane sampleprior: a uniform() prior on a slope, sampled,
# when the inits fall inside its support, and with lb/ub given.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
d <- D$wide
fit <- frm(y ~ x, data = d, family = gaussian())
try1 <- function(lab, pr, init_jitter = 0.25) {
  r <- tryCatch({
    ds <- suppressMessages(suppressWarnings(frm_sample(
      fit, prior = pr, chains = 2, iter = 1000, refresh = 0, seed = 5,
      init_jitter = init_jitter)))
    v <- ds$draws[, "b_x"]
    sprintf("OK range %.4f %.4f divergences %d", min(v), max(v),
            sum(rstan::get_divergent_iterations(ds$stanfit)))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 60)))
  cat(lab, ":", r, "\n")
}
try1("uniform(-10, 10)", set_prior("uniform(-10, 10)", class = "b", coef = "x"))
try1("uniform(2.2, 2.6) jitter 0", set_prior("uniform(2.2, 2.6)", class = "b",
                                          coef = "x"), 0)
try1("uniform(2.2, 2.6) lb ub", set_prior("uniform(2.2, 2.6)", class = "b",
                                       coef = "x", lb = 2.2, ub = 2.6))
