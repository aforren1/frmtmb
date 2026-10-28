# The prior draw names of the ported suite's draws fixtures 1 and 2
# (brms's example fits), and of an ar() model, on the lane's build.
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat); library(frmtmb); library(frmtmb.sample)
})
td <- "extensions/frmtmb.sample/tests/testthat"
source(file.path(td, "helper-brms-suite.R"))
source(file.path(td, "helper-brms-suite-draws.R"))
for (k in 1:2) {
  ds <- brms_fixture_draws(k)
  cat("fixture", k, "prior_draws:", names(prior_draws(ds)), "\n")
  cat("  ndraws", ndraws(ds), "warmup rows", NROW(ds$warmup), "\n")
}
want <- c("Intercept", "b", paste0("simo_moExp1[", 1:4, "]"), "bsp",
          "bs", "sds_sAge", "b_sigma", "Intercept_sigma", "nu",
          "sd_visit", "cor_visit")
cat("fixture 1 equals the suite's prior_names:",
    identical(names(prior_draws(brms_fixture_draws(1))), want), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), h = factor(rep(1:5, each = 12)),
                t = rep(1:12, 5))
d$y <- 1 + 0.5 * d$x + rnorm(60)
ds <- suppressWarnings(suppressMessages(frm_sample(
  bf(y ~ x + ar(time = t, gr = h)), data = d, family = gaussian(),
  prior = set_prior("normal(0, 0.5)", class = "ar"), sample_prior = "yes",
  chains = 1, iter = 200, seed = 1, refresh = 0)))
cat("ar model:", names(prior_draws(ds)), "\n")
cat("ar range:", range(prior_draws(ds)$ar), "\n")
