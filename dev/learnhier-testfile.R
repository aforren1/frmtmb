# Lane `learnhier`: run ONE test file in ONE process, with the lane's
# library stack and the gates the file needs.
#
#   Rscript dev/learnhier-testfile.R <file> [gates]
#
# `gates` is a comma-separated list from: stan (the identity tier),
# scale (the Phase 0 measurement tier).
source("dev/learnhier-env.R")
a <- commandArgs(trailingOnly = TRUE)
f <- a[[1L]]
gates <- if (length(a) >= 2L) strsplit(a[[2L]], ",")[[1L]] else character(0)
Sys.setenv(NOT_CRAN = "true")
if ("stan" %in% gates) Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
if ("scale" %in% gates) Sys.setenv(FRMTMB_SCALE_TESTS = "true")
# A third argument is a path the identity tier appends its residuals to,
# so the published table is that file rather than a transcription of it.
if (length(a) >= 3L) Sys.setenv(FRMTMB_LEARN_LP_LOG = normalizePath(
  a[[3L]], winslash = "/", mustWork = FALSE))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(frmtmb.learn))
setwd("extensions/frmtmb.learn/tests/testthat")
# IN THE PACKAGE'S NAMESPACE. `R CMD check` runs a suite through
# test_check(), which evaluates tests where the package's INTERNAL
# functions are visible. A bare test_file() after library() does not,
# and this lane's first full-suite run reported 5 ERRORs in
# test-counterfactual.R reading `ln_aterms`, `ln_family()` and
# `ln_counterfactual_of()`, none of which are exported. That was the
# runner, not the package: a runner that cannot see what the real one
# sees manufactures failures and, worse, could hide real ones behind
# them.
t0 <- Sys.time()
res <- test_file(f, reporter = "summary",
                 env = new.env(parent = asNamespace("frmtmb.learn")))
df <- as.data.frame(res)
cat("\nFILE ", f, "\n", sep = "")
cat("PASS=", sum(df$passed), " FAIL=", sum(df$failed),
    " ERROR=", sum(df$error), " SKIP=", sum(df$skipped),
    " WARN=", sum(df$warning),
    " secs=", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    "\n", sep = "")
# A runner that sums `failed` and not `error` prints a clean line for a
# file that aborted halfway; this one exits non-zero on either.
if (sum(df$failed) > 0 || sum(df$error) > 0) quit(status = 1L)
