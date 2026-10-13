# rel070 seflags: the B1 fixture (brms_fixture(1), seed 3, as in
# dev/seflags-rev3-pred3x.R) on one build: code, logLik, every fit
# warning, lost set, optrap's unit and polish; on a build with units,
# again with predictor_units() switched off in this process only.
#   Rscript rel070-seflags-fx1.R <lib> [<lib2>]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(testthat)})
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
ns <- asNamespace("frmtmb")
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "| R", R.home(), "\n")
root <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/tests/testthat"
show <- function(tag) {
  env <- new.env(parent = ns)
  for (h in list.files(root, "^helper-.*[.]R$", full.names = TRUE)) {
    sys.source(h, envir = env)
  }
  withr::local_seed(3L)
  fit <- suppressMessages(env$brms_fixture(1))
  l <- suppressWarnings(ns$sdr_of(fit)$se_lost)
  u <- if (exists("predictor_units", ns)) ns$predictor_units(fit$frame)
  cat(sprintf("%s | %s | code %d | logLik %.6f | lost %s | units %s | polish %s\n",
              tag, deparse(fit$call$formula)[1], fit$opt$convergence,
              as.numeric(logLik(fit)),
              if (length(l)) paste(names(l), l, sep = ":", collapse = " ") else "-",
              if (length(u)) paste(names(u), u) else "-",
              if (is.null(fit$opt$polish)) "-" else
                paste(signif(fit$opt$polish, 3), collapse = " ")))
  for (w in attr(fit, "brmsport_warnings")) {
    cat("   W:", substr(gsub("\n", " ", w), 1, 300), "\n")
  }
}
show("as installed")
if (exists("predictor_units", ns)) {
  unlockBinding("predictor_units", ns)
  assign("predictor_units", function(frame) numeric(0), envir = ns)
  show("units off   ")
}
