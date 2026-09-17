.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_multiple
### Title: Fit a model across multiply imputed datasets
### Aliases: frm_multiple

### ** Examples

set.seed(8)
n <- 80
x <- rnorm(n)
y <- rnorm(n, 1 + 0.5 * x, 1)
x[sample(n, 15)] <- NA
imps <- lapply(1:3, function(i) {
  xi <- x
  xi[is.na(xi)] <- sample(x[!is.na(x)], sum(is.na(xi)), TRUE)
  data.frame(y = y, x = xi)
})
frm_multiple(bf(y ~ x) + gaussian(), data = imps)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
