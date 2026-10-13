# variables() of a mo() fit on rellib-r7 (no simo_ names there).
#   Rscript dev/moparity-variables.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r7", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1); d <- data.frame(x = sample(0:3, 100, TRUE)); d$y <- d$x + rnorm(100)
f <- frm(bf(y ~ mo(x)), family = gaussian(), data = d)
cat("rellib-r7 variables():", variables(f), "\n")
