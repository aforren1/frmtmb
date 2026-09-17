.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: residuals.frmtmb_fit
### Title: Residuals from a frmtmb fit
### Aliases: residuals.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

# raw and variance-standardized residuals
head(residuals(fit))
head(residuals(fit, type = "pearson"))
# the usual overdispersion check for a poisson fit
sum(residuals(fit, type = "pearson")^2) / df.residual(fit)

# one-step-ahead quantile residuals are standard normal under a
# correctly specified model, whatever the family
r <- residuals(fit, type = "osa")
qqnorm(r); qqline(r)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
