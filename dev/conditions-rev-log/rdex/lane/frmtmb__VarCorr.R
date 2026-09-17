.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: VarCorr.frmtmb_fit
### Title: Extract random-effect standard deviations and correlations
### Aliases: VarCorr.frmtmb_fit VarCorr

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)

vc <- VarCorr(fit)
names(vc)                       # "g" and "residual__", as in brms
vc$g$sd                         # standard deviations
vc$g$cor["Intercept", "Estimate", "x"]
vc$residual__$sd



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
