# Lane nlstop2, fix round F1: the explin test's new body (any warning, and
# nl_limit_message() non-NULL on the fit) on ulp draws of y,
# y * (1 + 2e-15 k), k = 0 to 11, as the review's
# dev/reviews/nlstop2-review-explintest.R draws them.
# Usage: NLSTOP2_ARM=before|after Rscript dev/nlstop2-explintest.R
#   (or under an OpenBLAS copy's Rscript, with NLSTOP2_ARM set)
source("dev/nlstop2-common.R")
cat("blas:", R.home(), "\n")
ns <- asNamespace("frmtmb")
set.seed(2)
d <- data.frame(x = runif(1000, 0, 3))
d$y <- 2 * d$x + rnorm(1000, 0, 0.5)
pass <- 0
for (k in 0:11) {
  dk <- d
  dk$y <- dk$y * (1 + 2e-15 * k)
  r <- tt(frm(bf(y ~ a * (exp(b * x) - 1), a + b ~ 1, nl = TRUE),
              data = dk, start = list(beta = c(2, 0.5))))
  warned <- length(attr(r, "wc")) > 0
  lim <- !is.null(ns$nl_limit_message(r))
  ok <- warned && lim
  pass <- pass + ok
  cat(sprintf("ROW k %d code %s warn %s limit_message %s test-body %s\n", k,
              r$opt$convergence, wsum(r), lim, if (ok) "PASS" else "FAIL"))
}
cat("explin test body: pass", pass, "of 12\n")
cat("DONE\n")
