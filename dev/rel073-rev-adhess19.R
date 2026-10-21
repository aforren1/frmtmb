# Release reviewer, 0.73.0: the adhess item filed at 0.72.0,
# "`y | trunc(lb = 1) ~ x` with poisson: the gradient at the optimum is
# not finite and obj$fn() at last.par.best is Inf (seed 19)". The case
# of dev/adhess-rev-classes.R (whose common file left with its
# worktree), with the checks that item names, and an independent fit of
# the zero-truncated Poisson on stats::dpois / stats::ppois.
#   Rscript dev/rel073-rev-adhess19.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
set.seed(19)
x <- rnorm(400, 300, 60)
y <- rpois(400, exp(0.2 + 0.004 * (x - 300)))
d <- data.frame(x, y)[y >= 1, ]
w <- character()
fit <- withCallingHandlers(frm(y | trunc(lb = 1) ~ x, family = poisson(),
                               data = d),
                           warning = function(e) {
                             w <<- c(w, conditionMessage(e))
                             invokeRestart("muffleWarning")
                           })
cat("rows", nrow(d), "| convergence", fit$opt$convergence,
    fit$opt$message, "| warnings", length(w), "\n")
if (length(w)) cat("  ", substr(w, 1, 160), sep = "\n   ")
cat("\n")
lpb <- fit$obj$env$last.par.best
cat("obj$fn(last.par.best):", format(fit$obj$fn(lpb), digits = 12), "\n")
g <- fit$obj$gr(fit$opt$par)
cat("gradient at opt$par:", format(g, digits = 6), "| finite",
    all(is.finite(g)), "\n")
cat("obj$fn(opt$par):", format(fit$obj$fn(fit$opt$par), digits = 12),
    "| logLik", format(as.numeric(logLik(fit)), digits = 12), "\n")
# independent: the zero-truncated poisson on stats, by optim
nll <- function(b) {
  lam <- exp(b[1] + b[2] * d$x)
  -sum(dpois(d$y, lam, log = TRUE) - ppois(0, lam, lower.tail = FALSE,
                                             log.p = TRUE))
}
o <- optim(c(-1, 0.004), nll, method = "BFGS",
           control = list(reltol = 1e-14, maxit = 1000, parscale = c(1, 0.001)))
cat("stats reference: logLik", format(-o$value, digits = 12), "| b",
    format(o$par, digits = 8), "\n")
cat("frmtmb b:", format(fit$estimates$beta, digits = 8), "\n")
se <- sqrt(diag(vcov(fit)))
H <- optimHess(o$par, nll)
cat("SE frmtmb:", format(se, digits = 6), "| SE from optimHess of the stats",
    "reference:", format(sqrt(diag(solve(H))), digits = 6), "\n")
cat("END\n")
# the observed information on the eta scale, per row, by a central
# second difference in one dimension (step 1e-4 on eta), then X'hX: no
# step in the slope, whose x sits near 300
X <- cbind(1, d$x)
lr <- function(eta) {
  lam <- exp(eta)
  dpois(d$y, lam, log = TRUE) - ppois(0, lam, lower.tail = FALSE, log.p = TRUE)
}
e0 <- drop(X %*% o$par); h <- 1e-4
hh <- -(lr(e0 + h) - 2 * lr(e0) + lr(e0 - h)) / h^2
I <- crossprod(X * hh, X)
cat("SE from X'hX at the stats optimum:", format(sqrt(diag(solve(I))),
                                                  digits = 6), "\n")
