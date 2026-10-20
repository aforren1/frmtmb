# Reviewer of lane sampleprior: why a uniform() prior on a slope stops
# tmbstan: the bounds frm_sample() derives from the resolved prior.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
ns <- asNamespace("frmtmb.sample")
fit <- frm(y ~ x, data = D$wide, family = gaussian())
for (pr in list(set_prior("uniform(2.2, 2.6)", class = "b", coef = "x"),
                set_prior("normal(0, 10)", class = "b", coef = "x"))) {
  rp <- suppressMessages(ns$sample_resolve_priors(fit, pr, base = NULL,
                                                  defaults = TRUE))
  cat("prior", unclass(pr)[[1]]$prior, "\n")
  str(rp$ri$lower); str(rp$ri$upper)
  b <- ns$resolve_bounds(fit, unlist(as.list(rp$ri$lower)),
                         unlist(as.list(rp$ri$upper)))
  str(b)
}
