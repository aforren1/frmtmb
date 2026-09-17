.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_compat
### Title: Query the feature compatibility registry
### Aliases: frm_compat

### ** Examples

# one pair
frm_compat("rescor", "cens()")

# everything known about truncation
frm_compat("trunc()", status = c("refused", "broken"))

# the pairs to avoid
frm_compat(status = "broken")[, 1:5]



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
