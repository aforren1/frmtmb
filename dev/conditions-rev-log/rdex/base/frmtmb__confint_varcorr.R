.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: confint_varcorr
### Title: Natural-scale confidence intervals for covariance parameters
### Aliases: confint_varcorr

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)

# one row per SD and per correlation, on the scale they are read on
confint_varcorr(fit)

# confint() reports the same parameters on their internal scale, so
# the bounds there are log-SDs and Fisher-z correlations
confint(fit)[grep("^theta", rownames(confint(fit))), ]

# a fit with no random effects has no covariance parameters
confint_varcorr(frm(bf(y ~ x) + gaussian(), data = dd))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
