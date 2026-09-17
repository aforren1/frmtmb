.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: as_draws
### Title: Convert to a posterior draws object
### Aliases: as_draws as_draws_matrix as_draws_array as_draws_df
###   as_draws_list as_draws_rvars as_draws.frmtmb_fit
###   as_draws_matrix.frmtmb_fit as_draws_array.frmtmb_fit
###   as_draws_df.frmtmb_fit as_draws_list.frmtmb_fit
###   as_draws_rvars.frmtmb_fit as_draws.frmtmb_multiple

### ** Examples

# a maximum-likelihood fit carries no draws, and says so by name
dd <- data.frame(y = rnorm(40), x = rnorm(40))
try(as_draws_df(frm(bf(y ~ x) + gaussian(), data = dd)))

# frm_multiple() pools estimates rather than carrying draws, so it
# answers with the reason rather than a matrix
fits <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))
try(as_draws(fits))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
