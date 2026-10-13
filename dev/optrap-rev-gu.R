# Reviewer of lane optrap: the convergence check's trip-wire read in the
# full par_units (response's share included) instead of grad_units(),
# patched in this process only. Runs a test file, or the student
# sigma ~ x case at y * 1e-6 over seeds 1 to 20 (dev/optrap-rev-scale.R's
# generator).
#   Rscript dev/optrap-rev-gu.R <lib> <patch: full|none> <test file|student>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
if (a[2] == "full") {
  utils::assignInNamespace("grad_units", function(fit) {
    fit$par_units %||% 1
  }, "frmtmb")
}
cat("patch", a[2], "lib", dirname(find.package("frmtmb")), "\n")
if (a[3] == "student") {
  for (sc in c(1e-6, 1e-3, 1, 1e3, 1e6)) {
    nw <- 0L
    for (s in 1:20) {
      set.seed(s)
      d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
      b <- rnorm(20)[d$g]
      d$y <- (1 + 0.5 * d$x + 0.5 * b + exp(0.3 * d$x) * rt(100, 5)) * sc
      w <- character()
      f <- withCallingHandlers(
        suppressMessages(frm(bf(y ~ x + (1 | g), sigma ~ x),
                             family = student(), data = d)),
        warning = function(x) {
          w <<- c(w, conditionMessage(x))
          invokeRestart("muffleWarning")
        })
      nw <- nw + any(grepl("maximum absolute gradient|did not report", w))
    }
    cat("student sigma ~ x scale", sc, ": fits with a convergence warning",
        nw, "of 20\n")
  }
} else {
  res <- test_file(a[3], package = "frmtmb",
                   env = testthat::test_env("frmtmb"), reporter = "silent")
  r <- as.data.frame(res)
  for (t in res) for (x in t$results) {
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat("  FAIL [", t$test, "]", substr(conditionMessage(x), 1, 200), "\n")
    }
  }
  cat(sprintf("RESULT pass=%d fail=%d error=%d skip=%d warn=%d\n",
              sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
              sum(r$warning)))
}
