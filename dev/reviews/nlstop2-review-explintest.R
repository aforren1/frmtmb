# Review of lane nlstop2, F1: the explin test's body on ulp draws of y.
# Rscript dev/reviews/nlstop2-review-explintest.R
source("dev/reviews/nlstop2-review-common.R")
library(testthat)
set.seed(2)
d <- data.frame(x = runif(1000, 0, 3))
d$y <- 2 * d$x + rnorm(1000, 0, 0.5)
for (k in 0:5) {
  dk <- d
  dk$y <- dk$y * (1 + 2e-15 * k)
  r <- tt(frm(bf(y ~ a * (exp(b * x) - 1), a + b ~ 1, nl = TRUE),
              data = dk, start = list(beta = c(2, 0.5))))
  code <- r$v$opt$convergence
  res <- tryCatch({
    expect_true(length(r$w) > 0)
    if (identical(as.integer(code), 0L)) {
      expect_true("frmtmb_nl_limit" %in% r$w)
    }
    "PASS"
  }, expectation_failure = function(e) "FAIL")
  cat("ROW k", k, "code", code, "warn", paste(unique(r$w), collapse = "+"),
      "test-body", res, "\n")
}
cat("DONE\n")
