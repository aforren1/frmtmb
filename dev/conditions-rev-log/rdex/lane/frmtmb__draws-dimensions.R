.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: draws-dimensions
### Title: Size of a draws object
### Aliases: draws-dimensions ndraws nchains niterations nvariables
###   ndraws.frmtmb_fit nchains.frmtmb_fit niterations.frmtmb_fit
###   nvariables.frmtmb_fit

### ** Examples

dd <- data.frame(y = rnorm(40), x = rnorm(40))
fits <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))
try(ndraws(fits))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
