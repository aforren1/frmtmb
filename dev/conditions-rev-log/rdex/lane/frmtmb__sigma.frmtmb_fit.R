.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: sigma.frmtmb_fit
### Title: Residual standard deviation
### Aliases: sigma.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# the residual SD, on the response scale
sigma(fit)
# which is what the standardized residuals divide by
max(abs(residuals(fit, type = "pearson") -
          residuals(fit) / sigma(fit)))

# a poisson fit has no dispersion parameter, so sigma() is 1
dd$cnt <- rpois(100, exp(0.5 + 0.3 * dd$x))
sigma(frm(bf(cnt ~ x) + poisson(), data = dd))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
