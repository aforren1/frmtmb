.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.sample))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: log_lik
### Title: Pointwise log-likelihood of posterior draws
### Aliases: log_lik log_lik.frmtmb_draws

### ** Examples

## No test: 
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

  ll <- log_lik(ds)
  dim(ll)
  # the column means are the per-observation expected log-densities
  head(colMeans(ll))
}
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
