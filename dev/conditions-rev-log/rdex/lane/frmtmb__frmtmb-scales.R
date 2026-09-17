.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-scales
### Title: Which scale each method reports
### Aliases: frmtmb-scales

### ** Examples

set.seed(2026)
dd <- data.frame(x = rnorm(200))
dd$y <- exp(rnorm(200, 8 + 0.4 * dd$x, 0.4))
fit <- frm(bf(y ~ x) + lognormal(), data = dd)

# the default is the linear predictor, not the outcome
head(predict(fit))
head(fitted(fit))

# and the two are related by the family's own mean
head(exp(predict(fit) + sigma(fit)^2 / 2) - fitted(fit))

# sigma is printed on its log link and back-transformed by sigma()
summary(fit)$coefficients$sigma[1, 1]
sigma(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
