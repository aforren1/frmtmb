.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb_priorlist-columns
### Title: Column access on a prior specification
### Aliases: frmtmb_priorlist-columns $.frmtmb_priorlist
###   $<-.frmtmb_priorlist

### ** Examples

pr <- set_prior("normal(0, 2)", class = c("b", "sd"))
pr$class
pr$prior[2] <- "exponential(1)"
pr



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
