# tests/testthat/test-brms-lik-cdfs.R with brms_lp_check()'s report on,
# so each case prints "LPCHECK const <brms lp - logLik> grad <max |grad|>
# ours <logLik>", labelled here in the order the file runs them.
# Compiles come from FRMTMB_STAN_CACHE. Output: stdout (redirect to
# dev/cdfs-out/brms-lik-report.txt)
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
cat("CASE order: student skew_normal exgaussian shifted_lognormal",
    "von_mises Gamma inverse.gaussian Beta | asym_laplace |",
    "negbinomial geometric binomial beta_binomial poisson |",
    "zero_inflated_poisson zero_inflated_negbinomial",
    "zero_inflated_binomial zero_inflated_beta_binomial hurdle_poisson",
    "hurdle_negbinomial hurdle_gamma hurdle_lognormal zero_inflated_beta |",
    "trunc: negbinomial student von_mises hurdle_gamma",
    "zero_inflated_poisson poisson+cens poisson-lb0\n")
options(frmtmb.brms_lp_report = TRUE)
r <- as.data.frame(test_file("tests/testthat/test-brms-lik-cdfs.R",
                             package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent"))
cat("RESULT test-brms-lik-cdfs.R pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped), "\n",
    sep = "")
