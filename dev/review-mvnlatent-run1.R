# Reviewer's runner: one test file, one process, lane library first.
#   Rscript dev/review-mvnlatent-run1.R <package> <test file>
local({
  .libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
              "C:/Users/adf44/source/r/rellib-r11",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(NOT_CRAN = "true")
  suppressMessages(library(testthat))
  a <- commandArgs(trailingOnly = TRUE)
  p <- a[1]
  f <- a[2]
  suppressMessages(library(p, character.only = TRUE))
  cat("lib: ", p, " ", format(utils::packageVersion(p)), " ",
      dirname(find.package(p)), "; frmtmb ",
      format(utils::packageVersion("frmtmb")), " ",
      dirname(find.package("frmtmb")), "\n", sep = "")
  for (q in c("frmtmb.latent", "frmtmb.sample")) {
    cat("lib: ", q, " ", dirname(find.package(q)), "\n", sep = "")
  }
  res <- test_file(f, package = p, env = testthat::test_env(p),
                   reporter = "silent")
  r <- as.data.frame(res)
  for (t in res) for (x in t$results) {
    k <- if (inherits(x, "expectation_failure")) "FAIL" else
      if (inherits(x, "expectation_error")) "ERROR" else
        if (inherits(x, "expectation_skip")) "SKIP" else
          if (inherits(x, "expectation_warning")) "WARN" else NA
    if (!is.na(k)) cat("DETAIL ", k, " [", t$test, "] ",
                       gsub("\n", " | ", conditionMessage(x)), "\n", sep = "")
  }
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), " tests=", nrow(r), "\n", sep = "")
})
