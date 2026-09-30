# The four new rows of test-brms-likelihood.R alone, each reported by
# name, with the measured constant and gradient printed.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/stan-cache")
options(frmtmb.brms_lp_report = TRUE)
library(testthat)
library(frmtmb)
f <- "C:/Users/adf44/source/r/frmtmb-wt-formula2/tests/testthat/test-brms-likelihood.R"
for (d in c("row 17c: sigma1 = \"sigma2\" equates the components' sigma",
            "row 24: cmc = FALSE on a population- and a group-level term",
            "row 25: y ~ . expands against the data as in brms",
            "row 26: a list of families, one per response")) {
  cat("ROW", d, "\n")
  r <- as.data.frame(test_file(f, package = "frmtmb", desc = d,
                               env = test_env("frmtmb"),
                               reporter = "silent"))
  cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d\n", d,
              sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped)))
}
