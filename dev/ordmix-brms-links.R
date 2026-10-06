# brms 2.23.0's Links line for ordinal mixtures. Output:
# dev/ordmix-log-brms-links.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
sl <- function(f, fam) cat(deparse1(f$formula), "|",
                            brms:::summarise_links(bf(f, family = fam)), "\n")
sl(bf(y ~ x), mixture(cumulative, cumulative))
sl(bf(y ~ x), mixture(cumulative("probit"), sratio))
sl(bf(y ~ x, disc1 ~ z, theta1 ~ z), mixture(cumulative, cumulative))
sl(bf(y ~ x, hu1 ~ z), mixture(hurdle_cumulative, hurdle_cumulative))
sl(bf(y ~ x), cumulative())
sl(bf(y ~ x), mixture(gaussian, gaussian))
