# Lane wt-resmooth, nits round. The frmtmb-to-mgcv relative gap on the
# t2() block of test-smooth-population.R, so the tolerance there can be
# the file's own 1e-5 or carry a measured reason for being looser.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
rel_gap <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  max(abs(a - b)) / max(abs(c(a, b)))
}
set.seed(11)
n <- 240L
d <- data.frame(t = runif(n), g = factor(sample(letters[1:8], n, TRUE)))
d$y <- sin(3 * d$t) + rnorm(8, 0, 0.7)[d$g] * d$t + rnorm(n, 0, 0.3)
fit <- suppressWarnings(frm(bf(y ~ t2(t, g, bs = c("cr", "re"), k = 5)),
                            family = gaussian(), data = d))
gm <- suppressWarnings(mgcv::gam(y ~ t2(t, g, bs = c("cr", "re"), k = 5),
                                 data = d, method = "ML"))
cat("t2 frmtmb vs mgcv relative gap:",
    sprintf("%.3e", rel_gap(frm_linpred(fit, re_formula = NA),
                            predict(gm))), "\n")
cat("t2 NA vs NULL relative gap:",
    sprintf("%.3e", rel_gap(frm_linpred(fit, re_formula = NA),
                            frm_linpred(fit))), "\n")
cat("logLik frmtmb", sprintf("%.8f", as.numeric(logLik(fit))),
    "| mgcv ML", sprintf("%.8f", -gm$gcv.ubre), "\n")
