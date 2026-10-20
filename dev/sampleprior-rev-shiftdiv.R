# Reviewer of lane sampleprior: divergences of frm_sample() on a shifted
# lognormal with ndt near min(y), the wall the likelihood and the uniform
# default both put at min(y). Data seed 20261116, sampler seed 9.
lane <- !identical(Sys.getenv("REV_ARM"), "r9")
.libPaths(c(if (lane) "C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb.sample"), "\n")
set.seed(20261116)
n <- 60
d <- data.frame(x = rnorm(n))
d$y <- 0.4 + exp(rnorm(n, -0.5 + 0.2 * d$x, 0.4))
for (route in c("formula", "fit")) {
  m <- if (route == "formula") y ~ x else
    frm(y ~ x, data = d, family = frmtmb::brmsfamily("shifted_lognormal"))
  ds <- suppressMessages(suppressWarnings(frm_sample(
    m, data = if (route == "formula") d,
    family = if (route == "formula") frmtmb::brmsfamily("shifted_lognormal"),
    chains = 4, iter = 2000, refresh = 0, seed = 9)))
  v <- ds$draws[, "ndt"]
  cat(sprintf("%s: ndt q05/50/95 %s of min(y) %.4f, divergences %d of %d\n",
              route, paste(sprintf("%.4f", quantile(v, c(.05, .5, .95))),
                           collapse = " "), min(d$y),
              sum(rstan::get_divergent_iterations(ds$stanfit)), nrow(ds$draws)))
}
