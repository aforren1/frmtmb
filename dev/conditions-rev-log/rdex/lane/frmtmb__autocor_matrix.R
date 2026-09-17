.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: autocor_matrix
### Title: Estimated within-group residual correlation matrix
### Aliases: autocor_matrix

### ** Examples

set.seed(1)
d <- expand.grid(week = 1:5, subj = factor(1:25))
e <- as.vector(apply(matrix(rnorm(125), 5, 25), 2, function(z) {
  as.vector(stats::filter(z, 0.6, "recursive"))
}))
d$x <- rnorm(125)
d$y <- 1 + 0.5 * d$x + e
fit <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
           data = d)
autocor_matrix(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
