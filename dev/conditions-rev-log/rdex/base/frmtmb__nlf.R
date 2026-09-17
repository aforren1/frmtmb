.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: nlf
### Title: Add a nonlinear parameter formula to a model formula
### Aliases: nlf

### ** Examples

# the composed spelling of a nonlinear model
bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1)

# a nonlinear sigma with a linear mu
bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1)

# bodies chain: cc feeds a, a feeds mu
bf(y ~ a, nl = TRUE) + nlf(a ~ cc * x) + nlf(cc ~ exp(b)) + lf(b ~ 1)

# a variance function of the fitted mean: sd = exp(ls) * |mu|^th
bf(y ~ x) + nlf(sigma ~ ls + th * log(abs(mu))) + lf(ls ~ 1, th ~ 1)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
