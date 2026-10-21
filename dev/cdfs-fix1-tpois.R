# Fix round 1: the truncated poisson of test-trunc-postfit.R, whose mean
# residual moved from below 1e-8 to 1.3e-8. Is it the optimizer's stop
# or the gradient? Compare the objective and gradient with the exact
# truncated poisson likelihood.
#   Rscript dev/cdfs-fix1-tpois.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
y <- stats::rpois(6000, 3)
y <- y[y >= 2 & y <= 6][1:400]
fit <- frm(bf(y | trunc(lb = 2, ub = 6) ~ 1) + poisson(),
           data = data.frame(y = y))
b <- fit$opt$par
nll <- function(b) {
  mu <- exp(b)
  -sum(stats::dpois(y, mu, log = TRUE) -
         log(stats::ppois(6, mu) - stats::ppois(1, mu)))
}
ex <- stats::optimize(nll, c(0, 2), tol = 1e-14)
cat(sprintf("frmtmb b %.15f  exact argmin %.15f  diff %.2e\n", b, ex$minimum,
            b - ex$minimum))
cat(sprintf("objective frmtmb %.12f exact at b %.12f\n", fit$obj$fn(b), nll(b)))
h <- 1e-5
cat(sprintf("gradient AD %.6e  FD of exact %.6e\n", fit$obj$gr(b),
            (nll(b + h) - nll(b - h)) / (2 * h)))
cat(sprintf("mean residual %.3e  (sd(y) %.3f)\n",
            mean(residuals(fit)[, "Estimate"]), sd(y)))
cat("opt message:", fit$opt$message, "\n")
cat("END\n")
