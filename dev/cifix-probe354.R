# The numbers behind test-difference.R's "a gp() difference at ONE
# position cancels the kriging residual": the spread of the difference's
# standard error and estimate, in ulps, and the largest |eta| subtracted.
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(21)
n <- 90
d <- data.frame(x = sort(stats::runif(n, 0, 10)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.5, 0) + stats::rnorm(n, 0, 0.3)
gx <- d$x[-1] - diff(d$x) / 2
fit <- frm(bf(y ~ fac + gp(x)), family = stats::gaussian(), data = d)
A <- data.frame(x = gx, fac = factor("A", levels = levels(d$fac)))
B <- data.frame(x = gx, fac = factor("B", levels = levels(d$fac)))
dif <- frm_curve(fit, newdata = A, contrast = B, simultaneous = FALSE)
ulp <- function(v) 2^(floor(log2(abs(v))) - 52)
se <- dif$.se
cat("se spread:", format(diff(range(se)), digits = 4), "=",
    diff(range(se)) / ulp(se[1]), "ulps of", format(se[1], digits = 17),
    "\n")
cat("rows off row 1:", sum(se != se[1]), "\n")
est <- dif$.estimate
eta <- frm_linpred(fit, newdata = rbind(A, B))
cat("estimate spread:", format(diff(range(est)), digits = 4),
    " mean:", format(mean(est), digits = 6),
    " max |eta|:", format(max(abs(eta)), digits = 6),
    " spread / (eps max|eta|):",
    format(diff(range(est)) / (.Machine$double.eps * max(abs(eta))),
           digits = 4), "\n")
