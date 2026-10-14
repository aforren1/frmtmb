# Run ONE frmtmb.eam test file, against the worktree source (load_all)
# or against rellib-r8's installed, unmodified eam 0.11.3.
# Usage: Rscript dev/famuni-wiener-runtest.R worktree|rellib <test file>
# Prints testthat's summary reporter, then one RESULT line with counts.
args <- commandArgs(TRUE)
mode <- args[1]
f <- args[2]
stopifnot(mode %in% c("worktree", "rellib"), file.exists(f))
Sys.setenv(NOT_CRAN = "true")
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
eam <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam"
if (mode == "worktree") {
  suppressMessages(pkgload::load_all(eam, quiet = TRUE))
  env <- NULL
} else {
  suppressMessages(library(frmtmb.eam))
  env <- testthat::test_env("frmtmb.eam")
}
cat("mode:", mode, " frmtmb", format(packageVersion("frmtmb")), "lib:",
    find.package("frmtmb"), "\n")
cat("frmtmb.eam from:", if (mode == "worktree") eam else
  find.package("frmtmb.eam"), "\n")
r <- if (is.null(env)) {
  testthat::test_file(f, reporter = "summary", stop_on_failure = FALSE)
} else {
  testthat::test_file(f, reporter = "summary", package = "frmtmb.eam",
                      env = env, stop_on_failure = FALSE)
}
df <- as.data.frame(r)
cat(sprintf("RESULT %s %s tests=%d expectations=%d failed=%d error=%d skipped=%d warning=%d\n",
            mode, basename(f), nrow(df), sum(df$nb), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
