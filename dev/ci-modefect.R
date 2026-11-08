# The absent case of test-mo-coding.R's rewritten line 236: with new
# data coded from THEIR OWN minimum (the defect bb221eb3 fixed), both new
# assertions must fail; with the real coder both must pass.
#
#   Rscript dev/ci-modefect.R <lib1;lib2;...>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[[1L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$k0 <- d$k - 1L
d$k2 <- d$k + 1L
d$kn <- d$k - 3L
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
f1 <- frm(y ~ mo(k), data = d)
p1 <- fitted(f1, newdata = data.frame(k = 1:4))
own_min <- function(fit, lp, mi, newdata) {
  v <- eval(mi$expr, newdata, environment())
  as.integer(v - min(v))
}
outcome <- function(expr) {
  tryCatch({
    expr
    "pass"
  }, expectation_failure = function(e) "FAIL", error = function(e) {
    paste("ERROR:", conditionMessage(e))
  })
}
checks <- function() {
  a <- outcome(expect_identical(
    fitted(f1, newdata = data.frame(k = c(2L, 2:4)))[2:4, ], p1[2:4, ]))
  p3 <- fitted(f1, newdata = data.frame(k = 2:4))
  r <- max(abs(p3 - p1[2:4, ]) / abs(p1[2:4, ]))
  b <- outcome(expect_lt(r, 64 * .Machine$double.eps))
  cat(sprintf("  identity on 4 rows: %s; 3-row bound: %s (ratio %.3g)\n",
              a, b, r))
}
cat("real coder:\n")
checks()
cat("coder from the new data's own minimum:\n")
with_mocked_bindings(checks(), mo_codes = own_min, .package = "frmtmb")
