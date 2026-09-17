.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.sample))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: posterior_epred
### Title: Expected-value and predictive draws from sampled parameters
### Aliases: posterior_epred posterior_epred.frmtmb_draws posterior_linpred
###   posterior_linpred.frmtmb_draws posterior_predict
###   posterior_predict.frmtmb_draws

### ** Examples

## No test: 
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rpois(80, exp(0.3 + 0.4 * dd$x + rnorm(8, 0, 0.5)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

nd <- data.frame(x = c(-1, 0, 1),
                 g = factor(1, levels = levels(dd$g)))

# the expected response per draw: uncertainty in the mean
ep <- posterior_epred(ds, newdata = nd)
apply(ep, 2, quantile, c(0.025, 0.5, 0.975))

# the predictive distribution adds the family's own noise, so its
# intervals are wider
pp <- posterior_predict(ds, newdata = nd)
apply(pp, 2, quantile, c(0.025, 0.5, 0.975))

# the linear predictor itself, on the link scale by default
head(posterior_linpred(ds, newdata = nd, ndraws = 5))
}
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
