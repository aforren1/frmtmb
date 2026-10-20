# Lane illse2, fix round 1: test-se-reparam-reml.R with the REML
# Jacobian term removed (install 17's beta_whiten_nll()), so that the
# new logLik tests are seen to fail on the defect they pin.
#   Rscript dev/illse2-fx-nojac.R
.libPaths(c("C:/Users/adf44/source/r/wt-illse2-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
orig <- ns$beta_whiten_nll
assignInNamespace("beta_whiten_nll",
                  function(nll, bw, jacobian = FALSE) orig(nll, bw, FALSE),
                  "frmtmb")
res <- as.data.frame(test_file(
  "C:/Users/adf44/source/r/frmtmb-wt-illse2/tests/testthat/test-se-reparam-reml.R",
  package = "frmtmb", env = test_env("frmtmb"), reporter = "summary",
  stop_on_failure = FALSE))
print(res[, c("test", "passed", "failed")], row.names = FALSE)
