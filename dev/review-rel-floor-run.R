# Release review 0.80.0: one test file of an extension, its 0.80.0 build
# (a copy of rellib-r19's) loaded over rellib-r18 (frmtmb 0.79.1), every
# gate on. Tests whether an extension that kept its floor runs on the
# previous core.
#
#   Rscript dev/review-rel-floor-run.R <copy lib> <package> <test file>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
local({
  suppressMessages(library(testthat))
  suppressMessages(library(a[2], character.only = TRUE))
  cat("lib:", a[2], format(packageVersion(a[2])), dirname(find.package(a[2])),
      "; frmtmb", format(packageVersion("frmtmb")),
      dirname(find.package("frmtmb")), "\n")
  r <- as.data.frame(testthat::test_file(a[3], package = a[2],
                                         env = testthat::test_env(a[2]),
                                         reporter = "silent"))
  cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
              basename(a[3]), sum(r$passed), sum(r$failed),
              sum(r$error), sum(r$skipped), sum(r$warning)))
})
