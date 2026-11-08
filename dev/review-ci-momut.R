# Reviewer mutants for test-mo-coding.R line 236: which defects the OLD
# line (3-row identity) catches and the NEW pair (4-row identity plus
# 3-row 64-eps bound) does not, and the reverse.
#
#   Rscript dev/review-ci-momut.R <lib1;lib2;...>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[[1L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), "\n")
# the test file's own data (mo_code_data(), seed 5)
src <- parse("tests/testthat/test-mo-coding.R")
env <- new.env(parent = asNamespace("frmtmb"))
for (e in src) {
  if (is.call(e) && identical(e[[1L]], as.name("<-"))) eval(e, env)
}
d <- env$mo_code_data()
f1 <- frm(y ~ mo(k), data = d)
p1 <- fitted(f1, newdata = data.frame(k = 1:4))
real_codes <- frmtmb:::mo_codes
real_vals <- frmtmb:::mo_col_values
eps <- .Machine$double.eps
outcome <- function(expr) {
  tryCatch({ expr; "pass" }, expectation_failure = function(e) "FAIL",
           error = function(e) paste("ERROR:", conditionMessage(e)))
}
checks <- function() {
  old <- outcome(expect_identical(fitted(f1, newdata = data.frame(k = 2:4)),
                                  p1[2:4, ]))
  a <- outcome(expect_identical(
    fitted(f1, newdata = data.frame(k = c(2L, 2:4)))[2:4, ], p1[2:4, ]))
  p3 <- fitted(f1, newdata = data.frame(k = 2:4))
  r <- max(abs(p3 - p1[2:4, ]) / abs(p1[2:4, ]))
  b <- outcome(expect_lt(r, 64 * eps))
  sprintf("old 3-row identity: %-4s | new 4-row identity: %-4s | new 3-row bound: %-4s (ratio %.3g = %.1f eps)",
          old, a, b, r, r / eps)
}
muts <- list(
  none = list(),
  own_min = list(mo_codes = function(fit, lp, mi, newdata) {
    v <- eval(mi$expr, newdata, environment()); as.integer(v - min(v))
  }),
  # codes by row position: right on 1:4 and on c(2, 2:4) rows 2:4
  by_position = list(mo_codes = function(fit, lp, mi, newdata) {
    cd <- real_codes(fit, lp, mi, newdata); seq_along(cd) - 1L
  }),
  # codes from the new data's distinct values (factor coding)
  by_levels = list(mo_codes = function(fit, lp, mi, newdata) {
    cd <- real_codes(fit, lp, mi, newdata); as.integer(factor(cd)) - 1L
  }),
  # a whole category off, only when new data lack the fitted minimum
  shift_if_no_min = list(mo_codes = function(fit, lp, mi, newdata) {
    cd <- real_codes(fit, lp, mi, newdata)
    if (min(cd) > 0L) pmax(cd - 1L, 0L) else cd
  }),
  # 8 ulps on the term's value when new data lack the fitted minimum
  ulp8_if_no_min = list(mo_col_values = function(fit, mi, cd) {
    v <- real_vals(fit, mi, cd)
    if (min(cd) > 0L) v * (1 + 8 * eps) else v
  }),
  # 8 ulps on the term's value on a 3-row prediction only: the rounding
  # the lane's change is built to tolerate
  ulp8_if_3rows = list(mo_col_values = function(fit, mi, cd) {
    v <- real_vals(fit, mi, cd)
    if (length(cd) == 3L) v * (1 + 8 * eps) else v
  }),
  # 1e-12 relative on a 3-row prediction only: above the bound
  rel1e12_if_3rows = list(mo_col_values = function(fit, mi, cd) {
    v <- real_vals(fit, mi, cd)
    if (length(cd) == 3L) v * (1 + 1e-12) else v
  })
)
for (nm in names(muts)) {
  m <- muts[[nm]]
  res <- if (!length(m)) checks() else
    do.call(with_mocked_bindings, c(list(quote(checks())), m,
                                    list(.package = "frmtmb")))
  cat(sprintf("%-17s %s\n", nm, res))
}
