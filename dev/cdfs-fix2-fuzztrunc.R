# Fix round 2, N3: the fuzz plan's trunc() specs on a given library,
# with every triaged finding printed, now that the stale `trunc-postfit`
# entry no longer absorbs them.
#   Rscript dev/cdfs-fix2-fuzztrunc.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10", .libPaths()))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
invisible(capture.output(testthat::source_test_helpers("tests/testthat",
                                                       env = globalenv())))
plan <- fuzz_plan(seed = 20260901L, size = 300L)
sub <- plan[plan$aterm == "trunc", , drop = FALSE]
cat("trunc specs:", nrow(sub), "\n")
res <- fuzz_run(sub)
sm <- fuzz_summary(res)
cls <- vapply(sm$triaged, function(f) f$class, "")
print(table(cls))
for (f in sm$triaged[cls != "pass"]) cat(f$class, "|", fuzz_format(f), "\n")
cat("END\n")
