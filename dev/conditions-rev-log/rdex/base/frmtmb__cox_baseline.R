.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: cox_baseline
### Title: The fitted baseline-hazard simplex of a 'cox()' fit.
### Aliases: cox_baseline

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(200))
dd$time <- rexp(200, exp(-0.5 + 0.7 * dd$x))
fit <- frm(bf(time ~ x), family = cox(), data = dd)
cox_baseline(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
