# Runs tests/testthat/test-drmtmb-agreement.R against the installed
# frmtmb, one file per process as dev/lane-rules.md requires. The first
# argument chooses the library set: "with" puts drmtmb-lib first,
# "without" leaves it off so the file must skip. An optional second
# argument names another test file, for the mutation runs.
args <- commandArgs(TRUE)
libs <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(args[1], "with")) {
  libs <- c("C:/Users/adf44/source/r/drmtmb-lib", libs)
}
.libPaths(libs)
cat("drmTMB visible:", nzchar(system.file(package = "drmTMB")),
    " NOT_CRAN:", Sys.getenv("NOT_CRAN"), "\n")
library(testthat)
library(frmtmb)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-drmtmb/tests/testthat"
tf <- if (length(args) > 1) {
  args[2]
} else file.path(wt, "test-drmtmb-agreement.R")
t0 <- proc.time()[["elapsed"]]
r <- testthat::test_file(tf, reporter = "summary", package = "frmtmb",
                         load_package = "installed")
df <- as.data.frame(r)
cat("\nelapsed s:", round(proc.time()[["elapsed"]] - t0, 1), "\n")
cat("tests:", nrow(df), " expectations:", sum(df$nb), " failed:",
    sum(df$failed), " error:", sum(df$error), " skipped:", sum(df$skipped),
    " warnings:", sum(df$warning), "\n")
print(df[, c("test", "nb", "failed", "skipped", "error", "warning")])
