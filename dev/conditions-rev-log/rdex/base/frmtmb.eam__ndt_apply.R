.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ndt_apply
### Title: The non-decision time on the response's own scale
### Aliases: ndt_apply

### ** Examples

# no grouping: `ndt` is already a time
ndt_apply(list(ndt = 0.2))
# grouped: a fraction of each row's own bound
ndt_apply(list(ndt = 0.5), list(ndt_floor = c(0.4, 0.6)))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
