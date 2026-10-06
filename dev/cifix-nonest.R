# A rank-deficient fit: does frm_curve() read the rows core cannot
# estimate (frm_lp_basis()$nonest) the way frm_linpred() reports them?
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3)
n <- 200
d <- data.frame(x1 = stats::rnorm(n), z = stats::runif(n))
d$x2 <- 2 * d$x1
d$y <- 1 + d$x1 + sin(2 * pi * d$z) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ x1 + x2 + s(z, k = 6)), data = d)
nd <- data.frame(z = seq(0.05, 0.95, length.out = 7), x1 = 0.5,
                 x2 = c(1, 1, 1, 2, 1, 1, 1))
lb <- frm_lp_basis(fit, newdata = nd)
cat("nonest:", lb$nonest, "\n")
p <- suppressWarnings(frm_linpred(fit, newdata = nd, se.fit = TRUE))
cat("linpred se:", format(p$se.fit, digits = 4), "\n")
r <- tryCatch(suppressWarnings(frm_curve(fit, newdata = nd,
                                         simultaneous = FALSE)),
              error = function(e) conditionMessage(e))
if (is.character(r)) cat("frm_curve refused:", substr(r, 1, 160), "\n") else
  cat("frm_curve se:", format(r$.se, digits = 4), "\n")
