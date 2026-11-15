# Review of lane selimit2, item 6: does test-gp-noncentered.R's block at
# line 481 (and the file's other blocks) reach any function the lane
# changed? Counts calls under trace() while the whole file runs.
.libPaths(c("C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
fns <- c("se_limit_joint_dirs", "se_limit_joint_dir", "se_limit_separating",
         "se_limit_ranges_of", "se_limit_ranges", "jc_limit_extrap",
         "se_boundary_stop", "joint_cov_white_inner", "se_lost_why",
         "jc_nan_lost", "se_edge_limit", "se_tier3", "joint_cov_repair",
         "sdr_unwhiten")
cnt <- new.env()
for (f in fns) {
  assign(f, 0L, cnt)
  local({
    nm <- f
    suppressMessages(trace(nm, where = ns, print = FALSE,
                           tracer = bquote(assign(.(nm), get(.(nm), envir = .(cnt)) + 1L,
                                                  envir = .(cnt)))))
  })
}
f <- "C:/Users/adf44/source/r/frmtmb-wt-selimit2/tests/testthat/test-gp-noncentered.R"
res <- testthat::test_file(f, package = "frmtmb",
                           env = testthat::test_env("frmtmb"),
                           reporter = "silent")
df <- as.data.frame(res)
cat("tests", nrow(df), "failed", sum(df$failed), "error", sum(df$error), "\n")
for (f in fns) cat(sprintf("%-22s %d\n", f, get(f, envir = cnt)))
