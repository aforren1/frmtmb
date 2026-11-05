# Item 9: the lane's new test files on rellib-r16: which assertions fail, and do they
# fail on behavior or on a missing symbol.
#   Rscript dev/review-predread-seen.R <package> <file>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-predread/dev/stan-cache",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(testthat))
p <- a[1]
suppressMessages(library(p, character.only = TRUE))
cat("lib:", p, dirname(find.package(p)), "\n")
setwd(if (p == "frmtmb") "C:/Users/adf44/source/r/frmtmb-wt-predread/tests/testthat" else
  "C:/Users/adf44/source/r/frmtmb-wt-predread/extensions/frmtmb.sample/tests/testthat")
r <- testthat::test_file(basename(a[2]), package = p, env = testthat::test_env(p),
                         reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(r)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n", basename(a[2]),
            sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped), sum(d$warning)))
for (t in r) {
  for (x in t$results) {
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat(sprintf("%s | %s | %s\n", class(x)[1], t$test,
                  substr(gsub("[\r\n]+", " ", conditionMessage(x)), 1, 170)))
    }
  }
}
