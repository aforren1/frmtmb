# Reviewer of lane optstart: run one test file against a library and
# print every failed, errored, skipped or warned expectation.
#   Rscript dev/optstart-rev-runtest.R <lib> <package> <file>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
p <- a[2]
suppressMessages(library(p, character.only = TRUE))
suppressMessages(library(testthat))
cat("lib:", dirname(find.package(p)), " frmtmb from",
    dirname(find.package("frmtmb")), "\n")
res <- testthat::test_file(a[3], package = p, env = testthat::test_env(p),
                           reporter = "silent")
df <- as.data.frame(res)
for (r in res) {
  for (e in r$results) {
    if (!inherits(e, "expectation_success")) {
      cat("--", class(e)[1], "|", r$test, "\n   ",
          substr(gsub("\n", " ", conditionMessage(e)), 1, 400), "\n")
    }
  }
}
cat(sprintf("RESULT %s pass %d fail %d err %d skip %d warn %d\n",
            basename(a[3]), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
