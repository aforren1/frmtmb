# Reviewer of lane sampleprior: a uniform() prior on a slope, sampled in a
# fresh session (no tmbstan capture), formula and fit routes.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
d <- D$wide
for (route in c("formula", "fit")) {
  r <- tryCatch({
    m <- if (route == "formula") y ~ x else frm(y ~ x, data = d, family = gaussian())
    ds <- suppressMessages(suppressWarnings(frm_sample(
      m, data = if (route == "formula") d, family = if (route == "formula") gaussian(),
      prior = set_prior("uniform(2.2, 2.6)", class = "b", coef = "x"),
      chains = 2, iter = 1000, refresh = 0, seed = 5)))
    v <- ds$draws[, "b_x"]
    sprintf("range %.4f %.4f divergences %d", min(v), max(v),
            sum(rstan::get_divergent_iterations(ds$stanfit)))
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat(route, ":", r, "\n")
}
# the same with normal(0, 10), to see whether the error is the uniform's
r <- tryCatch({
  ds <- suppressMessages(suppressWarnings(frm_sample(
    y ~ x, data = d, family = gaussian(),
    prior = set_prior("normal(0, 10)", class = "b", coef = "x"),
    chains = 2, iter = 1000, refresh = 0, seed = 5)))
  "OK"
}, error = function(e) paste("ERROR", conditionMessage(e)))
cat("formula normal(0,10):", r, "\n")
