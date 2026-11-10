# Release review of 0.81.0: one test file in one process.
#   Rscript dev/review-rel081-tf.R <lib> <pkg> <test file path>
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]; pkg <- a[2]; f <- a[3]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", dirname(find.package(pkg)), pkg, format(packageVersion(pkg)),
    "frmtmb", dirname(find.package("frmtmb")),
    format(packageVersion("frmtmb")), "\n")
r <- as.data.frame(testthat::test_file(f, package = pkg,
                                       env = testthat::test_env(pkg),
                                       reporter = "silent"))
cat("RESULT", basename(f), "pass", sum(r$passed), "fail", sum(r$failed),
    "err", sum(r$error), "skip", sum(r$skipped), "warn", sum(r$warning), "\n")
bad <- r[r$failed > 0 | r$error | r$warning > 0, "test"]
if (length(bad)) cat("BAD:", paste(bad, collapse = " | "), "\n")
