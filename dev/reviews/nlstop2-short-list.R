# Short check of lane nlstop2's fix round: run one test file and list
# every expectation it records, with its test and source line, so two
# BLAS builds can be compared expectation by expectation.
# Usage: Rscript nlstop2-short-list.R <test file>; NLSTOP2_SHORT_ARM=before
# drops the lane library (rellib-r21 alone).
local({
  base <- "C:/Users/adf44/source/r/rellib-r21"
  user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
  lane <- "C:/Users/adf44/source/r/wt-nlstop2-lib"
  arm <- Sys.getenv("NLSTOP2_SHORT_ARM", "after")
  .libPaths(c(if (arm != "before") lane, base, user))
  suppressMessages(library(testthat))
  suppressMessages(library(frmtmb))
  f <- commandArgs(trailingOnly = TRUE)[1]
  cat("blas: ", R.home(), " | ", extSoftVersion()[["BLAS"]], " | ",
      La_library(), " | coretype ", Sys.getenv("OPENBLAS_CORETYPE"), "\n",
      sep = "")
  cat("lib: frmtmb ", as.character(utils::packageVersion("frmtmb")), " ",
      dirname(find.package("frmtmb")), " arm ", arm, "\n", sep = "")
  res <- test_file(f, package = "frmtmb",
                   env = testthat::test_env("frmtmb"), reporter = "silent")
  for (t in res) {
    for (x in t$results) {
      k <- class(x)[1]
      ln <- if (!is.null(x$srcref)) x$srcref[1] else NA
      cat("EXP ", k, " line ", ln, " [", t$test, "] ",
          substr(gsub("\n", " | ", conditionMessage(x)), 1, 160), "\n",
          sep = "")
    }
  }
  r <- as.data.frame(res)
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
})
