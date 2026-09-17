.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.spline))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_curve
### Title: A fitted curve on a grid, with pointwise and simultaneous bands
### Aliases: frm_curve

### ** Examples

set.seed(1)
dd <- data.frame(x = sort(runif(200)))
dd$y <- 2 * sin(pi * dd$x) + rnorm(200, 0, 0.4)
fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
                   family = stats::gaussian(), data = dd)
cv <- frm_curve(fit, newdata = data.frame(x = seq(0, 1, length.out = 25)),
                nsim = 2000)
head(cv[, c("x", ".estimate", ".se", ".lower_ci", ".lower_sim")])



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
