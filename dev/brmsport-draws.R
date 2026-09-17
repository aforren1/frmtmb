# Draw the frmtmb.sample stand-ins for brms's example fits 1, 2, 3 and 5,
# the four the sample half of tests.brmsfit-methods.R reads, and report
# the counts brms's suite asserts (25 draws, 1 chain) and the time.
#
#   Rscript dev/brmsport-draws.R > dev/brmsport-log/draws.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(identical(as.character(utils::packageVersion("StanHeaders")),
                    "2.32.10"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.sample)
})
h <- new.env()
sys.source("tests/testthat/helper-brms-suite.R", envir = h)
sys.source(file.path("extensions/frmtmb.sample/tests/testthat",
                     "helper-brms-suite-draws.R"), envir = h)
for (k in c(1, 2, 3, 5)) {
  t0 <- proc.time()[["elapsed"]]
  ds <- tryCatch(h$brms_fixture_draws(k), error = function(e) e)
  if (inherits(ds, "error")) {
    cat(k, "ERROR", conditionMessage(ds), "\n")
    next
  }
  cat(sprintf("fixture %d: %s, ndraws %d, nchains %d, %.1f s\n", k,
              class(ds)[1], ndraws(ds), nchains(ds),
              proc.time()[["elapsed"]] - t0))
}
