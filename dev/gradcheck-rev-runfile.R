## Reviewer: run ONE test file in ONE process against a chosen core.
## usage: Rscript gradcheck-rev-runfile.R <core-lib> <test-file> [pkg]
## Counts come from the testthat results object, so an aborted file
## reports its error rather than a clean zero, and every failure or
## error message is printed (no ten-failure cap).
a <- commandArgs(TRUE)
CORE <- a[1]
FILE <- a[2]
PKG <- if (length(a) >= 3) a[3] else NA_character_
.libPaths(c(CORE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
# the scale tier is the one that reads diagnose()$max_grad on every row,
# so it is exactly what this review has to compare between the builds
Sys.setenv(FRMTMB_SCALE_TESTS = "true")
cat("CORE  :", CORE, "\n")
cat("FILE  :", FILE, "\n")
cat("frmtmb:", find.package("frmtmb"), "\n")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
tenv <- NULL
if (!is.na(PKG)) {
  cat("ext   :", find.package(PKG), "\n")
  suppressMessages(library(PKG, character.only = TRUE))
  # the suites call the package's OWN internals, which only resolve when
  # the test environment's parent is its namespace (what test_check does)
  tenv <- new.env(parent = asNamespace(PKG))
}
if (is.null(tenv)) tenv <- new.env(parent = asNamespace("frmtmb"))
res <- try(testthat::test_file(FILE, reporter = "silent", env = tenv),
           silent = TRUE)
if (inherits(res, "try-error")) {
  cat("RESULT", basename(FILE), "ABORTED\n")
  cat(as.character(res), "\n")
} else {
  df <- as.data.frame(res)
  p <- sum(df$passed); f <- sum(df$failed)
  e <- sum(df$error); s <- sum(df$skipped)
  cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d\n",
              basename(FILE), p, f, e, s))
  for (tb in res) {
    for (ex in tb$results) {
      if (inherits(ex, "expectation_failure") ||
          inherits(ex, "expectation_error")) {
        cat("### ", class(ex)[1], " in: ", tb$test, "\n", sep = "")
        cat(conditionMessage(ex), "\n")
      } else if (inherits(ex, "expectation_skip")) {
        cat("--- skip in: ", tb$test, ": ",
            conditionMessage(ex), "\n", sep = "")
      }
    }
  }
}
