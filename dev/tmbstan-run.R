# lane tmbstan: run ONE test file in this process and print a count
# that includes errors. testthat's default summary caps failures at ten
# and a runner that sums `failed` alone prints a clean line for a file
# that aborted, so every field is printed here.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
args <- commandArgs(trailingOnly = TRUE)
f <- args[1L]
if (length(args) > 1L && nzchar(args[2L])) {
  Sys.setenv(FRMTMB_BRMS_FIT_TESTS = args[2L])
}
# poison the build detector's memo before anything loads a test file,
# which is how the BROKEN case is constructed at file scope
poison <- length(args) > 2L && identical(args[3L], "broken")

library(testthat)
library(frmtmb)
library(frmtmb.sample)
setwd("C:/Users/adf44/source/r/frmtmb-wt-tmbstan/extensions/frmtmb.sample")
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath(
  "../../dev/stan-cache", mustWork = FALSE))

if (poison) {
  env <- environment(frmtmb.sample:::tmbstan_build_broken)
  assign("cached", TRUE, envir = env)
  cat("### detector poisoned: tmbstan_build_broken() ==",
      frmtmb.sample:::tmbstan_build_broken(), "\n")
} else {
  cat("### detector as installed: tmbstan_build_broken() ==",
      frmtmb.sample:::tmbstan_build_broken(), "\n")
}

t0 <- proc.time()[["elapsed"]]
# package = names the namespace whose internals the blocks see.
# Without it test_file() parents the test environment on
# globalenv() and every unqualified call to an INTERNAL (here
# stan_cores() and loo_matrix()) errors, which is a property of
# the runner and not of the suite: test_check() does pass it.
r <- as.data.frame(test_file(file.path("tests/testthat", f),
                             package = "frmtmb.sample",
                             reporter = "silent"))
el <- proc.time()[["elapsed"]] - t0
cat("\n### FILE:", f, " arm:", if (poison) "BROKEN" else "clean", "\n")
cat("### blocks:", nrow(r), "\n")
for (k in c("failed", "error", "skipped", "warning", "passed")) {
  v <- if (k %in% names(r)) r[[k]] else NA
  cat("###  ", k, ": ",
      if (is.logical(v)) sum(v, na.rm = TRUE) else sum(v, na.rm = TRUE),
      "\n", sep = "")
}
cat("### elapsed_s:", round(el, 1), "\n")
if (any(r$failed > 0) || any(r$error)) {
  cat("\n### non-clean blocks:\n")
  bad <- r[r$failed > 0 | r$error, c("file", "test", "failed", "error")]
  print(bad, row.names = FALSE)
}
if (sum(r$skipped) > 0) {
  cat("\n### skipped blocks:", sum(r$skipped), "\n")
}
