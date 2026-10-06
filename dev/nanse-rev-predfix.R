# Reviewer: which fit of test-predfix.R warns twice (convergence, then
# the SE warning)? Traces se_check() and frm_warning() with the test name.
#   Rscript dev/nanse-rev-predfix.R <merge|lane> <test file>
a <- commandArgs(trailingOnly = TRUE)
libs <- switch(a[1],
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"),
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
suppressMessages(trace("se_check", where = ns, print = FALSE,
  tracer = quote(cat("  [se_check] code", fit$opt$convergence,
                     "explained", format(fit$cache$se_explained),
                     "autoscaled", !is.null(fit$par_units),
                     "| call", substr(paste(deparse(fit$call), collapse = ""),
                                      1, 160), "\n"))))
suppressMessages(trace("frm_warning", where = ns, print = FALSE,
  tracer = quote(cat("  [warn]", substr(paste0(..., collapse = ""), 1, 120),
                     "\n"))))
suppressMessages(trace("check_convergence", where = ns, print = FALSE,
  tracer = quote(cat("  [check_convergence] code", fit$opt$convergence,
                     "\n"))))
setwd(dirname(a[2]))
r <- test_file(a[2], package = "frmtmb", env = test_env("frmtmb"),
               reporter = LocationReporter$new())
