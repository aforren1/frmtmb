.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: posterior_summary
### Title: Summaries and intervals of draws
### Aliases: posterior_summary posterior_summary.default
###   posterior_summary.frmtmb_fit posterior_summary.frmtmb_multiple

### ** Examples

# any matrix of draws: rows are draws, columns are variables
m <- cbind(a = rnorm(500), b = rnorm(500, 2))
posterior_summary(m)
posterior_summary(m, robust = TRUE)

# a maximum-likelihood fit has no draws to summarize
dd <- data.frame(y = rnorm(40), x = rnorm(40))
try(posterior_summary(frm(bf(y ~ x) + gaussian(), data = dd)))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
