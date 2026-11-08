# Review lane osa: run one test file against the library named by
# REV_OSA_LIB (default wt-osa-lib) in front of rellib-r18.
#   Rscript dev/review-osa-test1.R <package> <test file>
local({
  lib1 <- Sys.getenv("REV_OSA_LIB", "C:/Users/adf44/source/r/wt-osa-lib")
  .libPaths(c(if (nzchar(lib1) && lib1 != "none") lib1,
              "C:/Users/adf44/source/r/rellib-r18",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(testthat))
  a <- commandArgs(trailingOnly = TRUE)
  p <- a[1]
  f <- a[2]
  suppressMessages(library(p, character.only = TRUE))
  cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
      dirname(find.package(p)), "; frmtmb ",
      as.character(utils::packageVersion("frmtmb")), " ",
      dirname(find.package("frmtmb")), "\n", sep = "")
  res <- NULL
  r <- tryCatch({
    res <- test_file(f, package = p, env = testthat::test_env(p),
                     reporter = "silent")
    as.data.frame(res)
  }, error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })
  if (!is.null(res)) {
    for (t in res) {
      for (x in t$results) {
        k <- if (inherits(x, "expectation_failure")) "FAIL" else
          if (inherits(x, "expectation_error")) "ERROR" else
            if (inherits(x, "expectation_skip")) "SKIP" else
              if (inherits(x, "expectation_warning")) "WARN" else NA
        if (!is.na(k)) {
          cat("DETAIL ", k, " [", t$test, "] ",
              substr(gsub("\n", " | ", conditionMessage(x)), 1, 400), "\n",
              sep = "")
        }
      }
    }
  }
  if (!is.null(r)) {
    cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
        sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
        " warn=", sum(r$warning), "\n", sep = "")
  }
})
