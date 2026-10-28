# Lane reslot, fix round 1 (review finding 3): frm() under the sampling
# table of ar(cov = TRUE), against ML, the uniform and the run's prior
# object. The data construction is the review's (dev/review-reslot-p4.R
# section iii): seed 13, 6 series of 8, AR(1) 0.6.
# Usage: Rscript dev/reslot-boxfit.R before|after
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(13)
G <- 6L; n <- 8L
da <- data.frame(g = factor(rep(seq_len(G), each = n)), t = rep(seq_len(n), G))
da$y <- 1 + unlist(lapply(seq_len(G), function(j)
  as.numeric(stats::arima.sim(list(ar = 0.6), n))))
fa <- bf(y ~ 1 + ar(t, g, cov = TRUE))
tab <- default_prior(fa, data = da, route = "sample")
print(as.data.frame(tab)[, c("prior", "class", "lb", "ub")], row.names = FALSE)
f_ml <- q(frm(fa, data = da))
f_tab <- q(frm(fa, data = da, prior = tab))
rp <- q(frmtmb.sample:::sample_resolve_priors(f_ml, NULL))
f_ps <- q(frm(fa, data = da, prior = rp$effective))
f_uni <- q(frm(fa, data = da,
               prior = set_prior("uniform(-1, 1)", class = "ar", lb = -1,
                                 ub = 1) +
                 set_prior(tab$prior[tab$class == "sigma"], class = "sigma") +
                 set_prior(tab$prior[tab$class == "Intercept"],
                           class = "Intercept")))
est <- function(f) unlist(f$estimates[c("beta", "betad", "thetaac")])
cat("  ML                  :", format(est(f_ml), digits = 8), "\n")
cat("  prior = table       :", format(est(f_tab), digits = 8), "\n")
cat("  prior = run's prior :", format(est(f_ps), digits = 8), "\n")
cat("  uniform(-1, 1) on ar:", format(est(f_uni), digits = 8), "\n")
