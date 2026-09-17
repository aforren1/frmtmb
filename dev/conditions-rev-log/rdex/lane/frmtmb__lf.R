.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: lf
### Title: Add parameter formulas to a model formula
### Aliases: lf

### ** Examples

# the two spellings are the same model
bf(y ~ x) + lf(sigma ~ z)
bf(y ~ x, sigma ~ z)

# nonlinear parameter formulas can arrive the same way
bf(y ~ a * exp(-b * x), a ~ 1, nl = TRUE) + lf(b ~ 1 + (1 | g))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
