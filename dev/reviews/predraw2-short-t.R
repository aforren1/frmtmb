# Short check of predraw2's fix round, item 5: one test file on the lane
# library, reference BLAS. Usage: Rscript predraw2-short-t.R <test file>
local({
.libPaths(c("C:/Users/adf44/source/r/wt-predraw2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
f <- commandArgs(trailingOnly = TRUE)[1]
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "| La:", La_library(), "\n")
r <- as.data.frame(testthat::test_file(f, package = "frmtmb",
                                       env = testthat::test_env("frmtmb"),
                                       reporter = "silent"))
cat(sprintf("%s: tests %d pass %d fail %d err %d skip %d warn %d\n",
            basename(f), nrow(r), sum(r$passed), sum(r$failed),
            sum(r$error), sum(r$skipped), sum(r$warning)))
bad <- r[r$failed > 0 | r$error | r$warning > 0, "test"]
if (length(bad)) cat("not clean:", paste(bad, collapse = " | "), "\n")
})
