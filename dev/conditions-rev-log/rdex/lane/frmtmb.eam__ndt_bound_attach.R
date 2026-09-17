.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ndt_bound_attach
### Title: Put a non-decision-time bound on a family
### Aliases: ndt_bound_attach

### ** Examples

# what a family_finalize() written against this seam looks like
finalize <- function(fam, y, aterms) {
  if (!is.null(ndt_bound_of(fam))) return(fam)
  ndt_bound_attach(fam, ndt_bound(y, aterms, what = "my_family"))
}



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
