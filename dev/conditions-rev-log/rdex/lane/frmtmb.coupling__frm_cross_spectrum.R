.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.coupling))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_cross_spectrum
### Title: The cross-spectrum of a pair of signals, as rows a model can
###   read
### Aliases: frm_cross_spectrum

### ** Examples

set.seed(1)
src <- rnorm(2048)
a <- src + rnorm(2048)
b <- 0.8 * src + rnorm(2048)
xs <- frm_cross_spectrum(a, b, sfreq = 256, segments = 8)
head(xs)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
