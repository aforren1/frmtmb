# Lane adhess: the margins of test-se-exact-hessian.R's bounds. Runs the
# test file with expect_lt() wrapped so that each bound prints its
# value, its limit and their ratio.
#   ADHESS_ARM=lane|base Rscript dev/adhess-tolmargin.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
suppressMessages(library(testthat))
tf <- "C:/Users/adf44/source/r/frmtmb-wt-adhess/tests/testthat/test-se-exact-hessian.R"
env <- new.env(parent = testthat::test_env("frmtmb"))
env$expect_lt <- function(object, expected, ...) {
  lab <- paste(deparse(substitute(object)), collapse = "")
  cat(sprintf("  %-22s value %.3e limit %.3e ratio %.3e\n", lab, object,
              expected, object / expected))
  testthat::expect_lt(object, expected, ...)
}
env$test_that <- function(desc, code) {
  cat(desc, "\n")
  r <- try(eval(substitute(code), new.env(parent = env)), silent = TRUE)
  if (inherits(r, "try-error")) cat("  ERROR", conditionMessage(attr(r, "condition")), "\n")
}
sys.source(tf, envir = env)
