# Lane sampleprior: frmtmb.sample's whole suite in ONE session, in
# alphabetical order, as covr and R CMD check run it (lane-env.md, "Test
# isolation"). Every extension is on the library path.
#   Rscript dev/sampleprior-onesession.R > dev/sampleprior-log/onesession.txt
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
r <- testthat::test_dir(
  "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/extensions/frmtmb.sample/tests/testthat",
  package = "frmtmb.sample", load_package = "none", stop_on_failure = FALSE,
  reporter = testthat::ProgressReporter$new(show_praise = FALSE))
d <- as.data.frame(r)
cat(sprintf("ONE SESSION: files %d tests %d pass %d fail %d error %d skip %d warn %d\n",
            length(unique(d$file)), nrow(d), sum(d$passed), sum(d$failed),
            sum(d$error), sum(d$skipped), sum(d$warning)))
bad <- d[d$failed > 0 | d$error, c("file", "test")]
if (nrow(bad)) print(bad, row.names = FALSE)
