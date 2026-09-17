.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: pp_check
### Title: Predictive check against simulated responses
### Aliases: pp_check pp_check.frmtmb_fit

### ** Examples

if (requireNamespace("bayesplot", quietly = TRUE)) {
  set.seed(1)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

  # the observed density against draws from the fit
  pp_check(fit, ndraws = 20)

  # any bayesplot ppc_* check, named by its suffix. A statistic the
  # model was not fitted to is the informative one: here, the share
  # of zeros, which is how zero inflation shows up.
  pp_check(fit, type = "stat", stat = function(y) mean(y == 0),
           ndraws = 50)
}



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
