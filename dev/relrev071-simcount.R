# Release reviewer of 0.71.0: an independent count of the statistical
# assertions of tests/testthat/test-simulate-density.R. Every expect_lt()
# call is recorded with the expression of its bound; a statistical
# assertion is one whose bound is made from SIM_ALPHA or SIM_Z (the
# file's alpha), whatever its label says.
#   Rscript dev/relrev071-simcount.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n")
rec <- character(0)
trace("expect_lt", where = asNamespace("testthat"), print = FALSE,
      tracer = quote(rec <<- c(rec, paste(deparse(substitute(expected)),
                                          collapse = " "))))
f <- "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat/test-simulate-density.R"
res <- as.data.frame(test_file(f, package = "frmtmb",
                               env = test_env("frmtmb"), reporter = "silent"))
untrace("expect_lt", where = asNamespace("testthat"))
cat("file: tests", nrow(res), "fail", sum(res$failed), "error",
    sum(res$error), "skip", sum(res$skipped), "\n")
stat <- grepl("SIM_ALPHA|SIM_Z", rec)
cat("expect_lt calls", length(rec), "| bound from the file's alpha",
    sum(stat), "\n")
print(table(rec))
