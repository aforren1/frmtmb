.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.spline))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_curve_deriv
### Title: Derivatives of a fitted curve
### Aliases: frm_curve_deriv

### ** Examples

set.seed(1)
dd <- data.frame(x = sort(runif(200)))
dd$y <- 2 * sin(pi * dd$x) + rnorm(200, 0, 0.4)
fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
                   family = stats::gaussian(), data = dd)
g <- data.frame(x = seq(0.05, 0.95, length.out = 19))
d1 <- frm_curve_deriv(fit, var = "x", order = 1, newdata = g,
                      simultaneous = FALSE)
# the derivative of 2 sin(pi x) is 2 pi cos(pi x)
head(cbind(g, fitted = d1$.estimate, truth = 2 * pi * cos(pi * g$x)))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
