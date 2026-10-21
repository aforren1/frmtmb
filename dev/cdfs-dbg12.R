# test-trunc-postfit.R's truncated poisson: the score at the stop, base
# build against this lane's, in one process each (argument: library)
a <- commandArgs(TRUE)
.libPaths(c(if (length(a)) a[1], "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
y <- stats::rpois(6000, 3)
y <- y[y >= 2 & y <= 6][1:400]
fit <- frm(bf(y | trunc(lb = 2, ub = 6) ~ 1) + poisson(), data = data.frame(y = y))
cat(dirname(find.package("frmtmb")), "\n")
cat(sprintf("b %.15g  grad %.3e  mean resid %.3e  iter %d  msg %s\n", fit$opt$par[1],
            fit$obj$gr(fit$opt$par)[1], mean(residuals(fit)[, "Estimate"]),
            fit$opt$iterations, fit$opt$message))
# exact root of the score by uniroot on the objective's own gradient
r <- uniroot(function(b) fit$obj$gr(b), fit$opt$par + c(-0.01, 0.01), tol = 1e-14)
cat(sprintf("root %.15g\n", r$root))
# the truncated mean at the root against the sample mean: the score
# equation of a one-parameter exponential family
mu <- exp(r$root)
k <- 2:6
w <- dpois(k, mu)
cat(sprintf("tmean(root) - mean(y) = %.3e\n", sum(k * w) / sum(w) - mean(y)))
