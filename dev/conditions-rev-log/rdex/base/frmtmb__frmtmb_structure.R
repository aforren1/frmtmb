.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb_structure
### Title: Declare a non-rowwise likelihood to the core
### Aliases: frmtmb_structure

### ** Examples

# a structure whose likelihood is the rowwise one, written out: the
# smallest thing the protocol accepts
st <- frmtmb_structure(
  loglik = function(y, dpars, aterms, weights, block, extra) {
    sum(weights * RTMB::dnorm(y, dpars$mu, dpars$sigma, log = TRUE))
  },
  supports = list(conditional_effects = TRUE)
)
st



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
