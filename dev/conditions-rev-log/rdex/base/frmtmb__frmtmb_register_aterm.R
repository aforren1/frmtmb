.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb_register_aterm
### Title: Add an addition term from another package
### Aliases: frmtmb_register_aterm

### ** Examples

# a decision indicator, given as a factor and delivered as 0/1
## Not run: 
##D frmtmb_register_aterm("dec", arity = 1, coerce = function(x) {
##D   as.integer(factor(x)) - 1L
##D })
## End(Not run)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
