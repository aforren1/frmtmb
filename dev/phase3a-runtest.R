# One test file per process. Usage:
#   Rscript phase3a-runtest.R <pkg source dir> <file> <arm> <out tsv>
# <arm> is "base" (the released build in rellib-r3) or "lane" (this
# lane's private library first). Appends one row per test_that() block
# to <out tsv> and prints a FILE line with the totals.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
Sys.setenv(PHASE3A_ARM = args[3L])
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
Sys.setenv(NOT_CRAN = "true")
pkg_dir <- args[1L]
pkg <- read.dcf(file.path(pkg_dir, "DESCRIPTION"), fields = "Package")[1L]
setwd(pkg_dir)
suppressMessages({
  library(testthat)
  library(frmtmb)
  library(pkg, character.only = TRUE)
})
phase3a_where("frmtmb"); phase3a_where(pkg)
# test_check() evaluates each file with the package NAMESPACE as parent,
# so a test that reaches an internal function works; a bare test_file()
# does not, and would report the runner's error as the code's.
# max_reports = Inf because the summary reporter otherwise stops
# listing at ten failures.
res <- testthat::test_file(
  file.path("tests/testthat", args[2L]),
  reporter = testthat::SummaryReporter$new(max_reports = Inf),
  env = testthat::test_env(pkg))
df <- as.data.frame(res)
rows <- data.frame(arm = args[3L], file = args[2L], test = df$test,
                   nb = df$nb, passed = df$passed, failed = df$failed,
                   error = df$error, skipped = df$skipped,
                   warning = df$warning)
utils::write.table(rows, args[4L], sep = "\t", quote = FALSE,
                   row.names = FALSE, append = file.exists(args[4L]),
                   col.names = !file.exists(args[4L]))
cat("\nFILE", args[2L], "ARM", args[3L], "TESTS", nrow(df),
    "PASS", sum(df$passed), "FAIL", sum(df$failed),
    "ERROR", sum(df$error), "SKIP", sum(df$skipped),
    "WARN", sum(df$warning), "\n")
