.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ndt_bound
### Title: The bound a non-decision time is measured against
### Aliases: ndt_bound

### ** Examples

rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
g <- c("a", "a", "a", "b", "b")
# `[[` and not `$`: the record has five names and `bd$s` would
# partial-match `sizes`
ndt_bound(rt)[["ub"]]
ndt_bound(rt, list(ndt_group = ndt_bound_key(g)))[["floors"]]



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
