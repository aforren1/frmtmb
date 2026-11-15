# Lane famend: run one test file of one package, one R process per file.
#   Rscript dev/famend-runtest.R <lib or "base"> <package> <test file>
# Prints the lib line, one line per failing or erroring expectation, and
# the RESULT line with pass, fail, error, skip and warning counts.
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
Sys.setenv(NOT_CRAN = "true")
p <- a[2]
f <- a[3]
suppressPackageStartupMessages({
  library(testthat)
  library(p, character.only = TRUE)
})
cat("lib:", find.package(p), as.character(packageVersion(p)),
    "| core:", find.package("frmtmb"), "\n")
t0 <- proc.time()[[3]]
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
for (i in seq_len(nrow(df))) {
  if (df$failed[i] > 0 || isTRUE(df$error[i])) {
    cat("FAIL/ERROR in:", df$test[i], "\n")
  }
}
for (r in res) {
  for (e in r$results) {
    if (inherits(e, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      cat("  ", class(e)[1], ":", substr(gsub("\\s+", " ",
                                              conditionMessage(e)), 1, 300),
          "\n")
    }
  }
}
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d time=%.0fs\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning),
            proc.time()[[3]] - t0))
