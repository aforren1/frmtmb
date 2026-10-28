# Reviewer: the lane's frmtmb.sample on frmtmb 0.76.0 (rellib-r14 core),
# which has no prior_only field: does "only" sample the posterior
# silently, and does test-prior-draws.R catch it?
#   Rscript dev/review-priordraws-floor.R
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-rv-lib/floor",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-stancache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "|",
    find.package("frmtmb.sample"), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
d$y <- 1 + 0.5 * d$x + rnorm(60)
ws <- character(0)
do <- withCallingHandlers(frm_sample(
  bf(y ~ x + (1 | g)), data = d, family = gaussian(),
  prior = set_prior("normal(0, 1)", class = "b"), sample_prior = "only",
  chains = 2, iter = 1000, seed = 5, refresh = 0),
  warning = function(w) {ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")},
  message = function(m) invokeRestart("muffleMessage"))
bx <- as.data.frame(do)$b_x
cat("only on 0.76.0 core: mean(b_x)", mean(bx), " sd(b_x)", sd(bx),
    " (prior normal(0,1)); warnings:", length(ws), "\n")
library(testthat)
r <- testthat::test_file(
  "C:/Users/adf44/source/r/frmtmb-wt-priordraws/extensions/frmtmb.sample/tests/testthat/test-prior-draws.R",
  package = "frmtmb.sample", env = testthat::test_env("frmtmb.sample"),
  reporter = "summary")
df <- as.data.frame(r)
cat(sprintf("RESULT floor test-prior-draws.R pass=%d fail=%d error=%d\n",
            sum(df$passed), sum(df$failed), sum(df$error)))
