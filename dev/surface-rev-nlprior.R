# Reviewer: nlpar coef prior after update() changes an nlpar formula
source("dev/surface-rev-env.R")
arm <- rev_env(commandArgs(TRUE)[1])
suppressPackageStartupMessages(library(frmtmb))
set.seed(11); n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$yn <- (2 + 0.3 * d$z) * exp(0.3 * d$x) + rnorm(n, 0, 0.3)
pn <- c(set_prior("normal(2, 1)", nlpar = "a", coef = "Intercept"),
        set_prior("normal(0, 1)", nlpar = "a", coef = "z"),
        set_prior("normal(0, 1)", nlpar = "b"))
fn <- rev_show("nl fit a ~ 1 + z",
  frm(bf(yn ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
      family = gaussian(), data = d, prior = pn))
if (!inherits(fn, "rev_err")) {
  u <- rev_show("update(a ~ 1)",
    update(fn, bf(yn ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE)))
  if (!inherits(u, "rev_err")) print(prior_summary(u))
}
