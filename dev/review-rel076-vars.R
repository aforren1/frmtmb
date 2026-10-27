# Review: does 0.75.0 (rellib-r13) report thetaac_1 among the draws of
# an ar(cov = TRUE) and a sar() model, as 0.76.0 does?
#   Rscript dev/review-rel076-vars.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", format(packageVersion("frmtmb")), "\n")
set.seed(1)
d <- expand.grid(t = 1:5, g = 1:10); d$y <- rnorm(50)
ds <- suppressWarnings(suppressMessages(frm_sample(y ~ ar(t, g, cov = TRUE),
  data = d, chains = 1, iter = 200, seed = 1, refresh = 0)))
cat("ar variables:", variables(ds), "\n")
