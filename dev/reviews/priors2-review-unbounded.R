# Review of lane priors2, item 2: does check_map_unbounded() refuse a
# fit that has a maximum? A prior infinite at lb = 0 on a coefficient
# whose likelihood goes to 0 at that end.
# Rscript priors2-review-unbounded.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
set.seed(11)
n <- 80
d <- data.frame(x = runif(n, 0.5, 2))
d$y <- rpois(n, 1.5 * d$x)
show <- function(lab, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR ", conditionMessage(r), "\n", sep = "")
  } else {
    cat(lab, ": code ", r$opt$convergence, " objective ",
        format(r$opt$objective, digits = 15), " par ",
        paste(format(r$opt$par, digits = 10), collapse = " "), "\n",
        sep = "")
  }
  invisible(r)
}
# A: poisson identity link, y ~ 0 + x, mu = b x. As b -> 0 the
# log-likelihood falls like sum(y) log b, the prior rises like
# -0.5 log b: the posterior goes to 0 at the end, a maximum inside.
pA <- set_prior("gamma(0.5, 1)", class = "b", coef = "x", lb = 0)
show("A gamma(0.5,1) lb=0, poisson(identity) y ~ 0 + x",
     frm(y ~ 0 + x, data = d, family = poisson(link = "identity"),
         prior = pA))
# the same posterior with the bound moved off the end (no refusal path)
show("A' gamma(0.5,1) lb=1e-10",
     frm(y ~ 0 + x, data = d, family = poisson(link = "identity"),
         prior = set_prior("gamma(0.5, 1)", class = "b", coef = "x",
                           lb = 1e-10)))
show("A'' no prior (ML)",
     frm(y ~ 0 + x, data = d, family = poisson(link = "identity")))
# the log posterior along b, by hand, to show the interior maximum
b <- c(1e-12, 1e-6, 1e-3, 0.1, 0.5, 1, 1.5, 2, 3)
lp <- sapply(b, function(bb) {
  sum(dpois(d$y, bb * d$x, log = TRUE)) + dgamma(bb, 0.5, 1, log = TRUE)
})
cat("A log posterior by hand at b =", b, "\n  ", format(lp, digits = 8),
    "\n")
# B: beta(0.5, 2) on a proportion-like slope with the likelihood 0 at 0
pB <- set_prior("beta(0.5, 2)", class = "b", coef = "x", lb = 0, ub = 1)
d$y2 <- rpois(n, 0.6 * d$x)
show("B beta(0.5,2) lb=0 ub=1, poisson(identity) y2 ~ 0 + x",
     frm(y2 ~ 0 + x, data = d, family = poisson(link = "identity"),
         prior = pB))
show("B' beta(0.5,2) lb=1e-10 ub=1",
     frm(y2 ~ 0 + x, data = d, family = poisson(link = "identity"),
         prior = set_prior("beta(0.5, 2)", class = "b", coef = "x",
                           lb = 1e-10, ub = 1)))
# C: the lane's own refusal case, a gaussian slope whose likelihood is
# finite at 0 (the refusal is right there)
set.seed(1)
dc <- data.frame(z = rnorm(120))
dc$y <- 1 + 0.3 * dc$z + rnorm(120, sd = 0.3)
show("C gamma(0.5,1) lb=0 gaussian y ~ z (likelihood finite at 0)",
     frm(y ~ z, data = dc, prior = set_prior("gamma(0.5, 1)", class = "b",
                                             coef = "z", lb = 0)))
