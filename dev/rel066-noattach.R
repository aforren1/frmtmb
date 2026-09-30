# Run one test file with or without attaching the package, to check the
# lane-rules claim that an unattached runner fails some files.
.libPaths(c("C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
f <- a[1]; attach_pkg <- identical(a[2], "attach")
suppressMessages(library(testthat))
if (attach_pkg) suppressMessages(library(frmtmb))
r <- as.data.frame(test_file(f, package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent"))
cat("NOATTACH ", basename(f), " attach=", attach_pkg, " pass=",
    sum(r$passed), " fail=", sum(r$failed), " err=", sum(r$error), "\n",
    sep = "")
