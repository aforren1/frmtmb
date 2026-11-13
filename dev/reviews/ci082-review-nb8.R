# Review of lane ci082: the code-8 branch of the new "ignore" assertion
# of test-se-limit.R, run on data where nlminb stops with code 8 on the
# reference BLAS (x scaled by 1 + 2e-15 * k, k = 3 and 4,
# dev/ci082-out/ref-nbk.txt). Prints the warning in full.
#   Rscript dev/reviews/ci082-review-nb8.R
.libPaths(c("C:/Users/adf44/source/r/wt-ci082-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
limit_conds <- function(expr) {
  w <- character()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}
for (k in c(0, 3, 4)) {
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.4 * d$x))
  d$x <- d$x * (1 + 2e-15 * k)
  rw <- limit_conds(frm(y ~ x, family = negbinomial(), data = d))
  ri <- limit_conds(frm(y ~ x, family = negbinomial(), data = d,
                        control = frmtmb_control(check_se = "ignore")))
  cat(sprintf("k %d code %s (%s)\n", k, ri$value$opt$convergence,
              ri$value$opt$message))
  cat("  W:", ri$warnings, sep = "\n  ")
  cat(sprintf("  par equal to the warning fit's: %s\n",
              identical(ri$value$opt$par, rw$value$opt$par)))
  # the lane's assertion, verbatim
  res <- tryCatch({
    expect_length(ri$messages, 0L)
    if (isTRUE(ri$value$opt$convergence == 0)) {
      expect_length(ri$warnings, 0L)
    } else {
      expect_length(ri$warnings, 1L)
      expect_match(ri$warnings, "Optimizer did not report convergence: ",
                   fixed = TRUE)
    }
    "pass"
  }, expectation_failure = function(e) paste("FAIL", conditionMessage(e)))
  cat(sprintf("RESULT nb8 k %d lane-assertion %s\n", k, res))
}
