# The dirichlet_multinomial tests with lgamma_shift_diff() replaced by
# the plain lgamma difference brms writes: the arm the first build was.
#   Rscript dev/nb2dm-naive-arm.R > dev/nb2dm-out/naive-arm.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-nb2dm-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
utils::assignInNamespace("lgamma_shift_diff",
                         function(a, s) lgamma(a + s) - lgamma(a), "frmtmb")
r <- as.data.frame(test_file("tests/testthat/test-dirichlet-multinomial.R",
                             package = "frmtmb", env = test_env("frmtmb"),
                             reporter = "silent", stop_on_failure = FALSE))
print(r[r$failed > 0 | r$error, c("test", "failed", "error")])
