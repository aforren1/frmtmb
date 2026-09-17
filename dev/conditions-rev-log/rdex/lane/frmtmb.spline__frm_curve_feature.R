.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.spline))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_curve_feature
### Title: Features of a fitted curve: a peak, a trough, or a level
###   crossing
### Aliases: frm_curve_feature

### ** Examples

set.seed(1)
dd <- data.frame(x = sort(runif(300)))
dd$y <- 2 * sin(pi * dd$x) + rnorm(300, 0, 0.3)
fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
                   family = stats::gaussian(), data = dd)
# 2 sin(pi x) peaks at x = 0.5
frm_curve_feature(fit, var = "x", type = "maximum",
                  newdata = data.frame(x = seq(0.05, 0.95, length.out = 41)))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
