.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: validate_prior
### Title: Check a prior against a model
### Aliases: validate_prior

### ** Examples

dd <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60),
                 g = factor(rep(1:6, 10)))
validate_prior(prior(normal(0, 10), class = b) +
                 prior(cauchy(0, 2), class = sd),
               y ~ x + z + (1 | g), data = dd)
# a prior on a coefficient the model does not have is refused
try(validate_prior(prior(normal(0, 1), coef = w), y ~ x, data = dd))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
