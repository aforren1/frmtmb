# Lane slopetrig: run ONE test file of one package against a library
# arm and print a RESULT line with pass, fail, error, skip and warning
# counts.
# Usage: Rscript slopetrig-testfile.R <base|lane> <pkg> <test file path>
a <- commandArgs(TRUE)
arm <- a[1]
pkg <- a[2]
f <- a[3]
stopifnot(arm %in% c("base", "lane"))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("arm", arm, pkg, as.character(packageVersion(pkg)),
    "lib:", dirname(find.package(pkg)), "frmtmb lib:",
    dirname(find.package("frmtmb")), "LAPACK:", La_library(), "\n")
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
for (t in res) for (e in t$results) {
  if (inherits(e, c("expectation_failure", "expectation_error",
                    "expectation_warning"))) {
    cat("--", class(e)[1], "in", t$test, ":\n",
        substr(conditionMessage(e), 1, 1500), "\n")
  }
}
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
