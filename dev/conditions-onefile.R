# Lane wt-conditions: run ONE test file in ONE process on a library arm,
# counting failures AND errors (dev/lane-rules.md).
#   Rscript dev/conditions-onefile.R lane|base <file>
av <- commandArgs(trailingOnly = TRUE)
arm <- av[1L]; file <- av[2L]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             mustWork = FALSE))
suppressMessages(library(testthat))
pkg <- regmatches(file, regexpr("frmtmb[.][a-z]+", file))
if (!length(pkg)) pkg <- "frmtmb"
suppressMessages(library(pkg, character.only = TRUE))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")
res <- test_file(file, reporter = "silent", package = pkg)
df <- as.data.frame(res)
cat(sprintf("BLOCKS %d PASS %d FAIL %d ERROR %d SKIP %d  %s\n", nrow(df),
            sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped),
            file))
for (i in seq_len(nrow(df))) {
  if (df$failed[i] > 0 || df$error[i]) {
    cat("  failing: ", df$test[i], "\n", sep = "")
    for (r in res[[i]]$results) {
      if (inherits(r, c("expectation_failure", "expectation_error")))
        cat("      ", substr(gsub("[\r\n]+", " ", conditionMessage(r)),
                             1, 500), "\n")
    }
  }
}
