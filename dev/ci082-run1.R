# Lane ci082: one test file in one process.
#   Rscript dev/ci082-run1.R <lib|base> <pkg> <test file path>
# Prints the lib line, then one RESULT line: pass, fail, error, skip,
# warn, and the elapsed time.
args <- commandArgs(TRUE)
arm <- args[1]
p <- args[2]
f <- args[3]
base <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-ci082-lib", base)
          else base)
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
cat("pkg lib:", dirname(system.file(package = p)), "\n")
cat("R:", R.home(), " CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
t0 <- proc.time()[["elapsed"]]
rep <- if (identical(Sys.getenv("CI082_REPORTER"), "location")) {
  testthat::LocationReporter$new()
} else testthat::ProgressReporter$new(show_praise = FALSE)
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = rep)
df <- as.data.frame(res)
el <- proc.time()[["elapsed"]] - t0
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d secs=%.1f\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning), el))
