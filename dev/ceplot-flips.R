# Lane ceplot: the ported brms-suite rows whose recorded verdict the
# lane build contradicts ("now HOLDS" or "no longer holds"), read from
# the gated methods files run with the summary reporter.
#   Rscript dev/ceplot-flips.R <lane|base>
a <- commandArgs(trailingOnly = TRUE)
arm <- if (length(a)) a[1] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib", libs)
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/stan-cache")
suppressMessages(library(testthat))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot"
for (p in c("frmtmb", "frmtmb.sample")) {
  suppressMessages(library(p, character.only = TRUE))
  f <- if (p == "frmtmb") {
    file.path(wt, "tests/testthat/test-brms-suite-methods.R")
  } else {
    file.path(wt, "extensions/frmtmb.sample/tests/testthat",
              "test-brms-suite-methods.R")
  }
  old <- setwd(dirname(f))
  r <- test_file(f, package = p, env = test_env(p),
                 reporter = ListReporter$new())
  setwd(old)
  for (t in r) {
    for (res in t$results) {
      if (inherits(res, "expectation_failure") ||
          inherits(res, "expectation_error")) {
        cat(arm, p, ":", sub("\n.*", "", conditionMessage(res)), "\n")
      }
    }
  }
}
