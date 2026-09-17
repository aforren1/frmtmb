.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ndt_bound_of
### Title: The non-decision-time bound a family already carries
### Aliases: ndt_bound_of

### ** Examples

# a family built outside frm() has no settled bound yet
ndt_bound_of(wiener())
ndt_bound_of(wiener(max_ndt = 0.4))[["ub"]]



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
