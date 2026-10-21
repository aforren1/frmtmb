# Run the fuzz specs that came up REAL-NEW in the cdfs gated run on a
# given library, to tell whether they are new with the lane or older.
#   Rscript dev/cdfs-fuzz-base.R <lib> <seed> [<seed> ...]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10", .libPaths()))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
testthat::source_test_helpers("tests/testthat", env = globalenv())
plan <- fuzz_plan(seed = 20260901L, size = 300L)
sub <- plan[plan$seed %in% as.integer(a[-1]), , drop = FALSE]
cat("specs:", nrow(sub), "\n")
res <- fuzz_run(sub)
sm <- fuzz_summary(res)
for (f in sm$triaged) cat(f$class, "|", fuzz_format(f), "\n")
cat("END\n")
