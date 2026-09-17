.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: bayes_R2
### Title: Bayesian R-squared
### Aliases: bayes_R2 bayes_R2.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)
try(bayes_R2(fit))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
