# Reviewer round 2, lane dblbar: is frm_sample(seed = ) reproducible
# within one process? Usage: Rscript ... base|lane. Data seed 20261013.
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
a <- as.matrix(frm_sample(fit, chains = 1, iter = 300, seed = 3,
                          refresh = 0))
b <- as.matrix(frm_sample(fit, chains = 1, iter = 300, seed = 3,
                          refresh = 0))
set.seed(1)
c1 <- as.matrix(frm_sample(fit, chains = 1, iter = 300, seed = 3,
                           refresh = 0))
set.seed(1)
c2 <- as.matrix(frm_sample(fit, chains = 1, iter = 300, seed = 3,
                           refresh = 0))
cat(arm, "same seed identical:", identical(a, b),
    " with set.seed before each:", identical(c1, c2), "\n")
