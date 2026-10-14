# Lane illse: run one test file of one package against the lane library
# (or, with ILLSE_ARM=base, against rellib-r8 alone), printing the
# library it loaded and the counts.
#   Rscript dev/illse-run1.R <package> <test file path>
args <- commandArgs(TRUE)
pkg <- args[1]
f <- args[2]
arm <- Sys.getenv("ILLSE_ARM", "lane")
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-illse-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", dirname(find.package(pkg)), "| frmtmb from",
    dirname(find.package("frmtmb")), "\n")
res <- as.data.frame(testthat::test_file(
  f, package = pkg, env = testthat::test_env(pkg), reporter = "summary",
  stop_on_failure = FALSE))
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(res$passed), sum(res$failed),
            sum(res$error), sum(res$skipped), sum(res$warning)))
