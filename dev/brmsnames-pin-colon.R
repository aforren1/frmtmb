## Punch round 1, BLOCKER: the colon pin of tests/testthat/test-brms-names.R
## ("hypothesis() reads x:fe as brms does"), run as a script on either
## build, so the old code's behavioral failure is on record: it RETURNS
## A NUMBER, b_x, for b_x:fe. The estimate is read from whichever object
## the build returns (0.58.0's data frame, or brms's list).
##   Rscript dev/brmsnames-pin-colon.R base|lane      data seed 8
arm <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/brmsnames-libs.R")
brmsnames_libs(arm)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
set.seed(8)
n <- 400
d <- data.frame(x = stats::rnorm(n),
                f = factor(sample(c("a", "e"), n, TRUE)))
d$y <- 1 + 0.3 * d$x + 0.8 * (d$f == "e") + 1.0 * d$x * (d$f == "e") +
  stats::rnorm(n)
fit <- q(frm(bf(y ~ x * f) + gaussian(), data = d))
fe <- fixef(fit)$mu
est_of <- function(h) {
  if (inherits(h, "brmshypothesis")) h$hypothesis$Estimate else h$estimate
}
got <- tryCatch(est_of(q(hypothesis(fit, "x:fe > 0"))),
                error = function(e) paste("ERROR:", conditionMessage(e)))
want <- unname(fe[["x:fe"]])
cat("arm", arm, "frmtmb from", dirname(system.file(package = "frmtmb")), "\n")
cat(sprintf("|b_x - b_x:fe| = %.4f (< 1, so R's `:` gives one number)\n",
            abs(fe[["x"]] - fe[["x:fe"]])))
cat("b_x =", format(fe[["x"]], digits = 7), " b_x:fe =",
    format(want, digits = 7), "\n")
cat("hypothesis(fit, \"x:fe > 0\") estimate:", format(got, digits = 7), "\n")
ok <- is.numeric(got) && abs(got - want) <= 64 * .Machine$double.eps * abs(want)
cat("pin:", if (ok) "PASS" else "FAIL", "\n")
