# The release review's short check of fix 2: bayes_R2() and loo_R2() with
# resp = character(0), against posterior_epred() on the same draws.
# Data seed 3, sampler seed 3, 1 chain of 300; rellib-r16 and rellib-r15.
#   Rscript dev/rel078-fix3-emptyresp.R <lib> > dev/rel078-log/fix3-emptyresp-<lib>.txt
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb.sample"), format(packageVersion("frmtmb.sample")), "\n")
set.seed(3)
d <- data.frame(x = rnorm(80)); d$y <- d$x + rnorm(80); d$y2 <- -d$x + rnorm(80)
ds <- suppressWarnings(suppressMessages(frm_sample(
  mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE), data = d, family = gaussian(),
  chains = 1, iter = 300, seed = 3, refresh = 0)))
tr <- function(e) tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
cat("bayes_R2(resp = character(0)):", tr(paste(rownames(bayes_R2(ds, resp = character(0))), collapse = " ")), "\n")
cat("loo_R2(resp = character(0)):", tr(paste(rownames(suppressWarnings(loo_R2(ds, resp = character(0)))), collapse = " ")), "\n")
cat("posterior_epred(resp = character(0)) responses:",
    tr(paste(dimnames(posterior_epred(ds, resp = character(0)))[[3]], collapse = " ")), "\n")
