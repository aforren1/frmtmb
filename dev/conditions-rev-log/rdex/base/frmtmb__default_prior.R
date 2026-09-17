.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: default_prior
### Title: Default priors: the slots a prior can target
### Aliases: default_prior get_prior

### ** Examples

dd <- data.frame(y = rnorm(60), x = rnorm(60),
                 g = factor(rep(1:6, 10)))
# what frm() applies: flat, whatever else is loaded. For what
# frm_sample() applies, see "Which route the defaults describe"
default_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
# the same table under brms's older name
get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
