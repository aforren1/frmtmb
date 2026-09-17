.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.ode))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_ode_failures
### Title: Groups whose ODE solve failed in the last 'frm_ode()' call
### Aliases: frm_ode_failures

### ** Examples

# NULL until a solve fails
frm_ode_failures()



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
