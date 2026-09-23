# Lane wt-correct: one test file, one R process, against either arm.
# Usage: Rscript dev/correct-runtest.R <base|lane> <package> <test-file>
# The base arm drops the lane library, so it reads rellib-r3 (0.61.0).
a <- commandArgs(trailingOnly = TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
if (identical(a[1], "base")) .libPaths(.libPaths()[-1L])
suppressMessages(library(testthat))
suppressMessages(library(a[2], character.only = TRUE))
cat("arm", a[1], a[2], format(packageVersion(a[2])), "from",
    dirname(find.package(a[2])), "\n")
r <- as.data.frame(test_file(a[3], package = a[2],
                             env = testthat::test_env(a[2]),
                             reporter = ProgressReporter$new(
                               show_praise = FALSE)))
cat("RESULT ", basename(a[3]), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    "\n", sep = "")
