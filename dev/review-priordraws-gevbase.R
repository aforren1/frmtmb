# Reviewer: does frm_sample() from a formula with gen_extreme_value()
# fail on rellib-r14 too (pre-existing)?
.libPaths(c("C:/Users/adf44/source/r/rellib-r14", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb.sample"), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60)); d$ye <- 1 + 0.5 * d$x + rnorm(60) + rexp(60)
r <- tryCatch({suppressWarnings(suppressMessages(frm_sample(bf(ye ~ x), data = d,
  family = gen_extreme_value(), prior = set_prior("normal(0, 1)", class = "b"),
  chains = 2, iter = 300, seed = 13, refresh = 0))); "OK"}, error = function(e) conditionMessage(e))
cat("base gev formula route:", r, "\n")
