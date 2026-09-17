.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: bf
### Title: Set up a model formula
### Aliases: bf

### ** Examples

# brms-style model formulas: attach a family with `+`
bf(y ~ x + (1 | g)) + gaussian()
# distributional parameters get their own formulas or constants
bf(y ~ x, sigma ~ x)
bf(y ~ x, shape = 2) + Gamma()
# nonlinear models declare parameter formulas and nl = TRUE
bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1 + (1 | g), nl = TRUE)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
