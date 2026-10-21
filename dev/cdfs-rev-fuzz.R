# Reviewer, fix round 2: the gated fuzz plan (300 specs, seed 20260901)
# on the lane build, every triaged finding printed with its class
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_FUZZ = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", find.package("frmtmb"), "\n")
testthat::source_test_helpers("C:/Users/adf44/source/r/frmtmb-wt-cdfs/tests/testthat", env = globalenv())
plan <- fuzz_plan(seed = 20260901L, size = 300L)
res <- fuzz_run(plan)
sm <- fuzz_summary(res)
cls <- vapply(sm$triaged, function(f) f$class, "")
print(table(cls))
for (f in sm$triaged) cat(f$class, if (!is.null(f$known_id)) paste0("[", f$known_id, "]"), "|", fuzz_format(f), "\n")
cat("END\n")
