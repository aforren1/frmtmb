# Lane ci082, fix round: the "ignore" assertion of test-se-limit.R, as
# the file now has it, on data where nlminb stops with code 8 on the
# reference BLAS (x scaled by 1 + 2e-15 * k, k = 3 and 4) and at code 0
# (k = 0). A second arm runs it on a fit forced to stop short (nlminb's
# iteration limit), which it must FAIL. After dev/reviews/ci082-review-nb8.R.
#   Rscript dev/ci082-nb8.R
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
assertion <- function(ri, r) {
  tryCatch({
    expect_length(ri$messages, 0L)
    if (isTRUE(ri$value$opt$convergence == 0)) {
      expect_length(ri$warnings, 0L)
    } else {
      expect_length(ri$warnings, 1L)
      expect_match(ri$warnings, "Optimizer did not report convergence: ",
                   fixed = TRUE)
      expect_match(ri$warnings, "false convergence (8)", fixed = TRUE)
      expect_match(ri$warnings, "so the estimates are at the maximum",
                   fixed = TRUE)
      expect_equal(unname(ri$value$opt$par[1:2] / r$value$opt$par[1:2]),
                   c(1, 1), tolerance = 1e-6)
    }
    "pass"
  }, expectation_failure = function(e) {
    paste("FAIL", gsub("\n", " | ", conditionMessage(e)))
  })
}
for (arm in c("normal", "short")) {
  for (k in c(0, 3, 4)) {
    set.seed(2)
    d <- data.frame(x = stats::rnorm(200))
    d$y <- stats::rpois(200, exp(1 + 0.4 * d$x))
    d$x <- d$x * (1 + 2e-15 * k)
    ctl <- if (arm == "short") {
      list(optCtrl = list(iter.max = 3, eval.max = 6), restarts = 0)
    } else list()
    r <- limit_conds(frm(y ~ x, family = negbinomial(), data = d,
                         control = do.call(frmtmb_control, ctl)))
    ri <- limit_conds(frm(y ~ x, family = negbinomial(), data = d,
                          control = do.call(frmtmb_control,
                                            c(ctl, check_se = "ignore"))))
    cat(sprintf("%s k %d code %s (%s)\n", arm, k, ri$value$opt$convergence,
                ri$value$opt$message))
    for (w in ri$warnings) cat("  W:", substr(w, 1, 300), "\n")
    cat(sprintf("RESULT nb8 %s k %d assertion %s\n", arm, k,
                assertion(ri, r)))
  }
}
