# Lane sampleprior: where a mixture component's intercept default comes
# from, in brms 2.23.0 and in frmtmb.sample (data seed 7).
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(7)
d <- data.frame(y = rpois(80, 20), x = rnorm(80))
cat("median log y", median(log(d$y)), "\n")
bp <- brms::default_prior(brms::bf(y ~ x), data = d,
                          family = brms::mixture("poisson", "poisson"))
print(bp[bp$class == "Intercept", c("prior", "class", "dpar")])
bt <- brms::brmsterms(brms::bf(y ~ x), family = brms::mixture("poisson",
                                                               "poisson"))
cat("brms mu1 link:", bt$dpars$mu1$family$link, " family link:",
    bt$family$link, "\n")
fp <- frmtmb::default_prior(bf(y ~ x), data = d,
                            family = mixture(poisson(), poisson()),
                            route = "sample")
print(fp[fp$class == "Intercept", c("prior", "class", "dpar")])
uf <- frm(bf(y ~ x), data = d, family = mixture(poisson(), poisson()),
          dry_run = "objective")
lp <- uf$frame$linpreds[[1]]
cat("frmtmb mu1 link:", lp$link$name, " type:", uf$spec$responses[[1]]$family$type,
    "\n")
str(frmtmb.sample:::default_prior_scale(uf, uf$spec$responses[[1]],
                                         link = lp$link$name))
