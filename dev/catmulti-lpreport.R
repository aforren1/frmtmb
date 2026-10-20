# test-refcat-brms.R with brms_lp_check()'s report switched on, so the
# measured constant and gradient of each check are printed.
# Output: dev/catmulti-out/lpreport.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib", "C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-catmulti/dev/stan-cache")
options(frmtmb.brms_lp_report = TRUE)
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
r <- as.data.frame(test_file("tests/testthat/test-refcat-brms.R", package = "frmtmb",
                             env = test_env("frmtmb"), reporter = "silent"))
print(r[, c("test", "passed", "failed", "error")])
