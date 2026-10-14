# Reviewer of lane postfix: run one test file against the lane library.
# Usage: Rscript dev/postfix-rev-testfile.R <pkg> <file>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-postfix-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
p <- a[1]
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", find.package(p), "\n")
r <- as.data.frame(testthat::test_file(a[2], package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent"))
cat(sprintf("RESULT pass=%d fail=%d err=%d skip=%d warn=%d\n",
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
bad <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
if (nrow(bad)) print(bad)
