.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: confint.frmtmb_fit
### Title: Confidence intervals for frmtmb fits
### Aliases: confint.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# Wald intervals for every parameter, covariance ones included
confint(fit)

# the likelihood profile does not assume a quadratic log-likelihood,
# so it is the one to trust for a variance component
confint(fit, parm = "theta_1", method = "profile")

# confint_varcorr() puts the same information on the SD scale
confint_varcorr(fit)
## No test: 
# a parametric bootstrap, the lme4 confint(method = "boot") analog
confint(fit, parm = "x", method = "boot", nsim = 50, seed = 1)
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
