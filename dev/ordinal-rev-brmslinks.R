# Reviewer: what brms 2.23.0's summary prints on the Links line of an
# ordinal model, without fitting (brms:::summarise_links on the formula).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
d <- data.frame(y = rep(1:4, 10), x = rnorm(40))
for (f in list(cumulative(), acat(), cumulative(link_disc = "softplus"))) {
  bf1 <- brms:::validate_formula(bf(y ~ x), data = d, family = f)
  cat(f$family, ": ", brms:::summarise_links(bf1, wsp = 0), "\n")
  bf2 <- brms:::validate_formula(bf(y ~ x, disc ~ 0 + x), data = d, family = f)
  cat(f$family, "with disc ~ x: ", brms:::summarise_links(bf2, wsp = 0), "\n")
}
print(body(brms:::summarise_links))
