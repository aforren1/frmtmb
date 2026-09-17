.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: vcov.frmtmb_fit
### Title: Covariance matrix of the fixed-effect estimates
### Aliases: vcov.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# standard errors of the fixed effects
sqrt(diag(vcov(fit)))
# the covariance parameters join the block on their internal scale
rownames(vcov(fit, full = TRUE))

# the matrix is what a delta-method calculation needs
V <- vcov(fit)
a <- c(1, 2)                       # prediction at x = 2, no group
sqrt(drop(t(a) %*% V[1:2, 1:2] %*% a))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
