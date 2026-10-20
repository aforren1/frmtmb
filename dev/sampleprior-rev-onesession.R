# Reviewer of lane sampleprior: frmtmb.sample's whole suite in ONE
# session, alphabetical, gated, against the lane library.
#   Rscript dev/sampleprior-rev-onesession.R
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/stan-cache")
suppressPackageStartupMessages(library(frmtmb.sample))
cat("lib:", find.package("frmtmb.sample"), " frmtmb:",
    find.package("frmtmb"), "\n")
d <- "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/extensions/frmtmb.sample/tests/testthat"
t0 <- Sys.time()
r <- testthat::test_dir(d, package = "frmtmb.sample", load_package = "none",
                        stop_on_failure = FALSE, reporter = "summary")
x <- as.data.frame(r)
cat(sprintf("ONE SESSION: files %d tests %d pass %d fail %d error %d skip %d warn %d (%.0f s)\n",
            length(unique(x$file)), nrow(x), sum(x$passed), sum(x$failed),
            sum(x$error), sum(x$skipped), sum(x$warning),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
bad <- x[x$failed > 0 | x$error | x$warning > 0, c("file", "test")]
if (nrow(bad)) print(bad)
