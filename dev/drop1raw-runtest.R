# Run one test file of one package against the lane or the base build.
# Args: arm (base | lane), package, test file path.
# Prints the library line, every failure and error with its test, and
# one RESULT line with the counts (the silent reporter caps nothing).
a <- commandArgs(TRUE)
arm <- a[1]; pkg <- a[2]; f <- a[3]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "| pkg lib:",
    dirname(find.package(pkg)), "\n")
r <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                         reporter = "silent", stop_on_failure = FALSE)
for (t in r) {
  for (e in t$results) {
    if (inherits(e, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      cat(sprintf("[%s] %s: %s\n", class(e)[1], t$test,
                  substr(gsub("\n", " ", conditionMessage(e)), 1, 300)))
    }
  }
}
df <- as.data.frame(r)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d tests=%d\n",
            basename(f), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning), nrow(df)))
