# Fix round 1: the trial merge's core, with this tree's tests and gpgap's two files, and the fresh row.
.libPaths(c("C:/Users/adf44/source/r/rev-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true", R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
td <- "C:/Users/adf44/source/r/rev-ncgp-src/mtree/tests/testthat"
for (f in c("test-gp-start.R", "test-gp-noncentered.R", "test-gp-by.R", "test-brms-lik-gp-exact.R")) {
  p <- file.path(td, f)
  if (!file.exists(p)) { src <- file.path("C:/Users/adf44/source/r/frmtmb-wt-gpgap/tests/testthat", f); if (file.exists(src)) file.copy(src, p) }
  r <- as.data.frame(test_file(p, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE))
  cat(sprintf("RMM %s tests %d pass %d fail %d err %d skip %d\n", f, nrow(r), sum(r$nb - r$failed - r$warning - r$skipped), sum(r$failed), sum(r$error), sum(r$skipped)))
  print(r[r$failed > 0 | r$error, "test"])
}
cat("DONE\n")
