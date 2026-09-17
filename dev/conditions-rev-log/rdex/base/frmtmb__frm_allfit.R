.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_allfit
### Title: Refit a model with every available optimizer
### Aliases: frm_allfit

### ** Examples

set.seed(6)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rpois(60, exp(0.3 + 0.4 * dd$x + rnorm(6, 0, 0.4)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
frm_allfit(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
