# Lane loglik, fix round 1, item 4: the test that sar() uses its slot,
# seen to fail on a mutant. sar() is registered again without its
# jacobian slot (the same call register_sar() makes, less that one
# argument), and test-sar.R runs against it. Output:
# dev/loglik-mutant-sarslot.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb.spatial)})
ns <- asNamespace("frmtmb.spatial")
frmtmb_register_autocor(
  "sar", args = c("M", "type"), parse = ns$sar_parse, build = ns$sar_build,
  whiten = ns$sar_whiten, unwhiten = ns$sar_unwhiten, npar = 1L,
  start = ns$sar_start, natural = ns$sar_natural, sigma_varies = FALSE,
  bounds = ns$sar_bounds, usage = "sar(M, type = \"lag\")",
  compat_note = ns$sar_compat_note, replace = TRUE,
  package = "frmtmb.spatial")
cat("slot after the mutant registration:",
    is.null(frmtmb:::ac_entry("sar")$jacobian), "(TRUE = no slot)\n")
r <- as.data.frame(test_file(
  "C:/Users/adf44/source/r/frmtmb-wt-loglik/extensions/frmtmb.spatial/tests/testthat/test-sar.R",
  package = "frmtmb.spatial", env = test_env("frmtmb.spatial"),
  reporter = "silent", stop_on_failure = FALSE))
cat("RESULT pass=", sum(r$passed), " fail=", sum(r$failed), " err=",
    sum(r$error), "\n", sep = "")
cat("FAILING:", r$test[r$failed > 0 | r$error], sep = "\n  ")
