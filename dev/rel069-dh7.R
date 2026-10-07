# Ported row data-helpers:7 under the OpenBLAS 0.3.32 scan: which
# warning does fitted(fixture 1, newdata) give, and is it the build or
# the BLAS? Prints every warning of the fixture fit and of the call.
#   Rscript dev/rel069-dh7.R <lib dir>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "| R", R.home(), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat"
for (h in list.files(root, "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
withr::local_seed(3L)
w <- character(); m <- character()
fit <- withCallingHandlers(env$brms_fixture(1),
  warning = function(x) { w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning") },
  message = function(x) { m <<- c(m, conditionMessage(x))
    invokeRestart("muffleMessage") })
cat("fixture: code", fit$opt$convergence, "objective",
    sprintf("%.10f", fit$opt$objective), "\n")
cat("fixture warnings:", length(w), "\n"); for (x in w) cat("  W:", substr(x, 1, 400), "\n")
for (x in m) cat("  M:", substr(x, 1, 300), "\n")
fit$data$fac <- factor(sample(1:3, nrow(fit$data), TRUE))
nd <- fit$data[1:5, ]
w <- character()
invisible(withCallingHandlers(env$brms_shim_validate_newdata(nd, fit),
  warning = function(x) { w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning") }))
cat("fitted(newdata) warnings:", length(w), "\n")
for (x in w) cat("  W:", substr(x, 1, 600), "\n")
# what the fit lost, and what the fixture's own handler recorded
cat("fixture recorded warnings:\n")
print(attr(fit, "brms_fixture_warnings") %||% env$brms_port_state$fixture_warnings)
cat("spec:\n"); print(fit$call$formula %||% fit$formula)
sdr <- suppressWarnings(frmtmb:::sdr_of(fit))
cat("se_lost:\n"); print(sdr$se_lost)
cat("se_explained:", format(fit$cache$se_explained), "\n")
print(summary(fit))
