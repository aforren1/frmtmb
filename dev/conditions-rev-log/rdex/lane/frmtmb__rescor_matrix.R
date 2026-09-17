.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: rescor_matrix
### Title: Estimated residual correlation matrix (rescor fits), else NULL
### Aliases: rescor_matrix

### ** Examples

set.seed(2)
n <- 80
dd <- data.frame(x = rnorm(n))
# two responses that share a residual disturbance
e <- rnorm(n)
dd$y1 <- 1 + 0.5 * dd$x + e + rnorm(n, 0, 0.5)
dd$y2 <- 2 - 0.3 * dd$x + e + rnorm(n, 0, 0.5)
fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE) + gaussian(),
           data = dd)
rescor_matrix(fit)

# a fit without rescor has no residual correlation to report
fit0 <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + gaussian(), data = dd)
rescor_matrix(fit0)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
