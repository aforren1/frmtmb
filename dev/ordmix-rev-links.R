# Reviewer of lane ordmix: brms 2.23.0's Links line for ordinal mixtures
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
for (f in list(bf(y ~ x, family = mixture(cumulative("probit"), cumulative())),
               bf(y ~ x, theta1 ~ z, family = mixture(cumulative(), sratio())),
               bf(y ~ x, family = mixture(gaussian(), gaussian())))) {
  bt <- suppressMessages(brms:::validate_formula(f, data = data.frame(y = sample(1:4, 50, TRUE), x = rnorm(50), z = rnorm(50))))
  cat(brms:::summarise_links(bt, wsp = 0), "\n")
}
