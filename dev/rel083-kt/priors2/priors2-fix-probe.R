# Lane priors2, fix round: the review's F1 cases, the 5005 pushed cases,
# a random-effect sd under gamma(0.5, 1), and a caller's start on a
# zero-density end (famend's case), on both builds.
# PRIORS2_LIB unset = lane library; PRIORS2_LIB=none = rellib-r21 only.
# Rscript dev/priors2-fix-probe.R > dev/priors2-out/fix-probe-<arm>.txt
lib <- Sys.getenv("PRIORS2_LIB", "C:/Users/adf44/source/r/rellib-r22")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
show <- function(lab, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", substr(conditionMessage(w), 1, 110), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR ", substr(conditionMessage(r), 1, 200), "\n", sep = "")
  } else {
    cat(lab, ": code ", r$opt$convergence, " objective ",
        format(r$opt$objective, digits = 15), " par ",
        paste(format(r$opt$par, digits = 10), collapse = " "), "\n",
        sep = "")
  }
  invisible(r)
}
# the review's data (dev/reviews/priors2-review-unbounded.R)
set.seed(11)
n <- 80
d <- data.frame(x = runif(n, 0.5, 2))
d$y <- rpois(n, 1.5 * d$x)
d$y2 <- rpois(n, 0.6 * d$x)
show("A gamma(0.5,1) lb=0 poisson(identity)",
     frm(y ~ 0 + x, data = d, family = poisson(link = "identity"),
         prior = set_prior("gamma(0.5, 1)", class = "b", coef = "x",
                           lb = 0)))
show("A' gamma(0.5,1) lb=1e-10",
     frm(y ~ 0 + x, data = d, family = poisson(link = "identity"),
         prior = set_prior("gamma(0.5, 1)", class = "b", coef = "x",
                           lb = 1e-10)))
show("B beta(0.5,2) lb=0 ub=1 poisson(identity)",
     frm(y2 ~ 0 + x, data = d, family = poisson(link = "identity"),
         prior = set_prior("beta(0.5, 2)", class = "b", coef = "x", lb = 0,
                           ub = 1)))
show("B' beta(0.5,2) lb=1e-10 ub=1",
     frm(y2 ~ 0 + x, data = d, family = poisson(link = "identity"),
         prior = set_prior("beta(0.5, 2)", class = "b", coef = "x",
                           lb = 1e-10, ub = 1)))
set.seed(1)
dc <- data.frame(z = rnorm(120))
dc$y <- 1 + 0.3 * dc$z + rnorm(120, sd = 0.3)
show("C gamma(0.5,1) lb=0 gaussian (review's)",
     frm(y ~ z, data = dc, prior = set_prior("gamma(0.5, 1)", class = "b",
                                             coef = "z", lb = 0)))
# the lane's own data (dev/priors2-bounds.R)
set.seed(1)
n <- 120
e <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
e$y <- 1 + e$z * 0.3 + rnorm(10)[e$g] * 0.3 + rnorm(n, sd = 0.3)
e$yneg <- e$y - 0.8 * e$z
for (f in c("y", "yneg")) {
  show(paste("D gamma(0.5,1) lb=0 gaussian", f),
       frm(stats::as.formula(paste(f, "~ z")), data = e,
           prior = set_prior("gamma(0.5, 1)", class = "b", coef = "z",
                             lb = 0)))
}
for (dn in c("lognormal(0, 1)", "beta(2, 2)", "gamma(2, 2)")) {
  ub <- if (startsWith(dn, "beta")) 1 else NA
  show(paste("E", dn, "lb=0 yneg (pushed)"),
       frm(yneg ~ z, data = e, prior = set_prior(dn, class = "b",
                                                 coef = "z", lb = 0,
                                                 ub = ub)))
  show(paste("E'", dn, "lb=1e-8 yneg"),
       frm(yneg ~ z, data = e, prior = set_prior(dn, class = "b",
                                                 coef = "z", lb = 1e-8,
                                                 ub = ub)))
}
# famend's case: a caller's start on the zero-density end
show("F beta(2,2) lb=0 ub=1 start z = 0",
     frm(yneg ~ z, data = e,
         prior = set_prior("beta(2, 2)", class = "b", coef = "z", lb = 0,
                           ub = 1), start = list(beta = c(1, 0))))
# a random-effect sd under gamma(0.5, 1) with lb = 0
show("G sd gamma(0.5,1) lb=0",
     frm(y ~ z + (1 | g), data = e,
         prior = set_prior("gamma(0.5, 1)", class = "sd", lb = 0)))
e$y0 <- 1 + e$z * 0.3 + rnorm(n, sd = 0.3)
show("G' sd gamma(0.5,1) lb=0, no group effect in the data",
     frm(y0 ~ z + (1 | g), data = e,
         prior = set_prior("gamma(0.5, 1)", class = "sd", lb = 0)))
# fix round 2: the short review's case H (dev/reviews/priors2-short-f1.R),
# an exposure per row of 1e5 to 1e9 that puts the maximum b* =
# (S - 0.5) / (sum x + 1) within the slack of lb = 0
set.seed(5)
hx <- runif(80, 0.5, 2)
for (scale in c(1e5, 1e7, 1e9)) {
  for (S in c(1L, 3L)) {
    h <- data.frame(x = hx * scale, y = 0L)
    h$y[seq_len(S)] <- 1L
    cat(sprintf("H b* = %.10g\n", (S - 0.5) / (sum(h$x) + 1)))
    show(sprintf("H gamma(0.5,1) lb=0 scale %g S %d", scale, S),
         frm(y ~ 0 + x, data = h, family = poisson(link = "identity"),
             prior = set_prior("gamma(0.5, 1)", class = "b", coef = "x",
                               lb = 0)))
    show(sprintf("H beta(0.5,2) lb=0 ub=1 scale %g S %d", scale, S),
         frm(y ~ 0 + x, data = h, family = poisson(link = "identity"),
             prior = set_prior("beta(0.5, 2)", class = "b", coef = "x",
                               lb = 0, ub = 1)))
  }
}
# the end label: beta(0.5, 2) is infinite at lb = 0, zero at ub = 1
show("D' beta(0.5,2) lb=0 ub=1 gaussian yneg",
     frm(yneg ~ z, data = e, prior = set_prior("beta(0.5, 2)", class = "b",
                                               coef = "z", lb = 0, ub = 1)))
# a NaN lower bound above the support (short review, item 3)
set.seed(1)
dz <- data.frame(z = rnorm(120))
invisible(rnorm(10))
invisible(rnorm(120))
dz$cnt <- rpois(120, 2) * rbinom(120, 1, 0.5)
for (lb in c(2, -1)) {
  show(sprintf("I zi lb = %g", lb),
       frm(cnt ~ z, data = dz, family = zero_inflated_poisson(),
           prior = set_prior("", class = "zi", lb = lb)))
}
# fix round 3: an interior fit beside an infinite-density end at which
# the likelihood is finite (the short review's re-check, P3)
set.seed(1)
d3 <- data.frame(x = rnorm(150), z = rnorm(150))
d3$y <- 1 + 0.4 * d3$x - 0.3 * d3$z + rnorm(150, sd = 0.7)
for (pr in c("gamma(0.5, 1)", "lognormal(0, 1)", "gamma(2, 2)")) {
  show(sprintf("J %s lb=0 on x, y ~ x + z", pr),
       frm(y ~ x + z, data = d3, prior = set_prior(pr, class = "b",
                                                   coef = "x", lb = 0)))
}
