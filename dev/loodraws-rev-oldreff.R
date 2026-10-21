# The gated same-draws test of test-loo-predict.R must fail with the
# frmtmb.sample 0.20.0 r_eff rule put back (r_eff = 1 on a subset of
# the draws). The old loo_r_eff() is mocked in, the file is run, and
# the failing blocks and their messages are printed.
#   Rscript dev/loodraws-rev-oldreff.R (copy of the lane's dev/loodraws-oldreff-check.R on the reviewer build)
#   -> dev/loodraws-rev-out/fix2-oldreff.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
suppressMessages({library(testthat); library(frmtmb)
  library(frmtmb.sample)})
old_loo_r_eff <- function(ll, x = NULL, chain_id = attr(ll, "chain_id"),
                          multiple = !is.null(x) &&
                            frmtmb.sample:::draws_multiple(x)) {
  # 0.20.0: 1 on combined imputations, NULL (so r_eff = 1) on a subset
  if (multiple) return(rep(1, ncol(ll)))
  if (is.null(chain_id)) return(rep(1, ncol(ll)))
  loo::relative_eff(exp(ll), chain_id = chain_id)
}
testthat::local_mocked_bindings(loo_r_eff = old_loo_r_eff,
                                .package = "frmtmb.sample",
                                .env = globalenv())
f <- "C:/Users/adf44/source/r/frmtmb-wt-loodraws/extensions/frmtmb.sample/tests/testthat/test-loo-predict.R"
res <- test_file(f, package = "frmtmb.sample",
                 env = testthat::test_env("frmtmb.sample"),
                 reporter = "silent", stop_on_failure = FALSE)
for (t in res) {
  for (x in t$results) {
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat("FAIL [", t$test, "] ",
          substr(gsub("\n", " | ", conditionMessage(x)), 1, 300), "\n",
          sep = "")
    }
  }
}
r <- as.data.frame(res)
cat("RESULT pass=", sum(r$passed), " fail=", sum(r$failed), " err=",
    sum(r$error), " skip=", sum(r$skipped), "\n", sep = "")
