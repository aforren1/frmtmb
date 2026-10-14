# The measured constant and gradient of every brms_lp_check() row of
# tests/testthat/test-simplex-brms.R. Output: dev/famsimplex-lpcheck.log
local({
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/stan-cache",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", find.package("frmtmb"), "\n")
options(frmtmb.brms_lp_report = TRUE)
f <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/tests/testthat/test-simplex-brms.R"
setwd(dirname(f))
r <- as.data.frame(test_file(f, package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent"))
print(r[, c("test", "nb", "failed", "error")])
})
