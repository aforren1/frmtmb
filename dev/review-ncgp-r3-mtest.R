# Short check: the trial merge, one test file per argument. And a gp beside mi() through the reader that retapes, on the merged core.
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/rev-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true", R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
td <- "C:/Users/adf44/source/r/rev-ncgp-src/mtree3/tests/testthat"
for (f in args) {
  r <- as.data.frame(test_file(file.path(td, f), package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE))
  cat(sprintf("RM3 %s tests %d pass %d fail %d err %d skip %d\n", f, nrow(r), sum(r$nb - r$failed - r$warning - r$skipped), sum(r$failed), sum(r$error), sum(r$skipped)))
  bad <- r[r$failed > 0 | r$error, "test"]
  if (length(bad)) print(bad)
}
set.seed(58001); x <- sort(runif(60, 0, 10)); z <- rnorm(60); zm <- z; zm[seq(5, 60, 12)] <- NA
d <- data.frame(x, zm, y = 2 * sin(3 * x) + x + 0.3 * z + 0.1 * rnorm(60))
f <- suppressWarnings(frm(bf(y ~ mi(zm) + gp(x)) + bf(zm | mi() ~ 1) + set_rescor(FALSE), d))
cat(sprintf("RM3 mi + gp ll %.4f conv %d; compress after %d matmul %s\n", logLik(f), f$opt$convergence, TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress, RTMB::TapeConfig()$matmul))
print(tail(ls(asNamespace("frmtmb"), pattern = "^mi_"), 3))
cat("DONE\n")
