# Reviewer, punch round 1: test-nl-rtmb-scope.R's `a * squash(b * x)`
# fit, where the projection rule in fixes' nl_flat_message() now warns
# on the trial merge. Is the a, b direction really flat?
#   Rscript dev/nanse-rev2-squash.R merge|release
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  release = "C:/Users/adf44/source/r/rellib-r6",
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
set.seed(3)
d <- data.frame(x = stats::runif(200, -2, 2))
d$y <- 1.4 * stats::pnorm(0.8 * d$x - 0.3) + stats::rnorm(200, 0, 0.15)
squash <- function(u) u / (1 + abs(u))
w <- character()
fit <- withCallingHandlers(
  frm(bf(y ~ a * squash(b * x), a ~ 1, b ~ 1, nl = TRUE),
      family = gaussian(), data = d, start = list(beta = c(1, 1))),
  warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
p <- fit$opt$par
cat("code", fit$opt$convergence, "par", signif(p, 6), "logLik",
    format(logLik(fit), digits = 10), "\n")
for (x in w) cat("warn:", substr(x, 1, 120), "\n")
H <- fit$obj$he(p)
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
cat("exact unit-diag eigenvalues:", signif(e$values, 4), "\n")
V <- suppressWarnings(vcov(fit, full = TRUE))
cat("reported SEs:", signif(sqrt(diag(V)), 4), "\n")
# the product a * b along the fit's own path: refit with b fixed larger
# shows whether the likelihood moves
f <- function(bb) {
  -fit$obj$fn(c(p[1] * p[2] / bb, bb, p[3]))
}
cat("logLik with a * b held, b = p_b * (1, 2, 0.5):",
    signif(sapply(c(1, 2, 0.5), function(s) f(p[2] * s)), 10), "\n")
