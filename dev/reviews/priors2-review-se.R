# Review of lane priors2, item 6: standard errors and the reported fit
# after the open chart of a zero-density end, against the same
# posterior boxed just inside the support.
# Rscript priors2-review-se.R
.libPaths(c("C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(seq_len(n / 12), each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(n / 12)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
one <- function(lab, pr, f = yneg ~ z) {
  w <- character(0)
  fit <- withCallingHandlers(frm(f, data = d, prior = pr),
                             warning = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
  s <- fixef(fit)
  cat(sprintf("%-34s code %d obj %.12f z %.10g se %.10g warn %d", lab,
              fit$opt$convergence, fit$opt$objective, s["z", "Estimate"],
              s["z", "Est.Error"], length(w)), "\n")
  if (length(w)) cat("   ", substr(w, 1, 140), sep = "\n   ")
  invisible(fit)
}
for (dn in c("lognormal(0, 1)", "gamma(2, 2)", "beta(2, 2)",
             "inv_gamma(3, 1)")) {
  ub <- if (startsWith(dn, "beta")) 1 else NA
  one(paste(dn, "lb=0"), set_prior(dn, class = "b", coef = "z", lb = 0,
                                   ub = ub))
  one(paste(dn, "lb=1e-8"), set_prior(dn, class = "b", coef = "z",
                                      lb = 1e-8, ub = ub))
}
# a random-effect model, the open chart on b beside the sd
one("glmm lognormal lb=0", set_prior("lognormal(0, 1)", class = "b",
                                     coef = "z", lb = 0),
    yneg ~ z + (1 | g))
one("glmm lognormal lb=1e-8", set_prior("lognormal(0, 1)", class = "b",
                                        coef = "z", lb = 1e-8),
    yneg ~ z + (1 | g))
# verbose: which run is reported
invisible(capture.output(
  fv <- frm(yneg ~ z, data = d, prior = set_prior("beta(2, 2)", class = "b",
                                                  coef = "z", lb = 0, ub = 1),
            verbose = TRUE),
  type = "message") -> vb)
cat(grep("open chart", vb, value = TRUE), sep = "\n")
# the chart on a two-sided box where only one end is at zero density
# (beta(1, 2): the density is finite at 0 and zero at 1)
one("beta(2,1) lb=0 ub=1 (zero at 0 only)",
    set_prior("beta(2, 1)", class = "b", coef = "z", lb = 0, ub = 1))
one("beta(2,1) lb=1e-8 ub=1",
    set_prior("beta(2, 1)", class = "b", coef = "z", lb = 1e-8, ub = 1))
