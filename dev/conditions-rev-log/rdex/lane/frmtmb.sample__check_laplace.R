.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.sample))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: check_laplace
### Title: Check the Laplace/Wald approximation against NUTS
### Aliases: check_laplace

### ** Examples

## No test: 
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
# a binary GLMM with small clusters: the regime where the Laplace
# approximation and Wald intervals are least reliable
set.seed(4)
dd <- data.frame(x = rnorm(120), g = factor(rep(1:30, 4)))
dd$y <- rbinom(120, 1,
               plogis(0.3 + 0.5 * dd$x + rnorm(30, 0, 1)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + frmtmb::bernoulli(), data = dd)

cl <- check_laplace(fit, chains = 1, iter = 500, refresh = 0)
cl
# |z_shift| well above 0 or sd_ratio far from 1 marks the parameters
# whose Wald interval to replace with a profile or bootstrap one
cl[abs(cl$z_shift) > 0.3 | cl$sd_ratio > 1.3, ]
}
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
