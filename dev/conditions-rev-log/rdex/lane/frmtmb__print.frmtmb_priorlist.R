.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: print.frmtmb_priorlist
### Title: Print a prior specification
### Aliases: print.frmtmb_priorlist

### ** Examples

set_prior("normal(0,1)", coef = "x")
set_prior("cauchy(0,1)", class = "sd", group = "g")
set_prior("normal(0, 2)", class = c("b", "sd"))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
