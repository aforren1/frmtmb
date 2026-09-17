.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.sample))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-loo-refusals
### Title: Refusals for the refit-based and marginal-likelihood brmsfit
###   methods
### Aliases: frmtmb-loo-refusals loo_moment_match
###   loo_moment_match.frmtmb_draws loo_subsample
###   loo_subsample.frmtmb_draws reloo reloo.frmtmb_draws kfold
###   kfold.frmtmb_draws bridge_sampler bridge_sampler.frmtmb_draws
###   bayes_factor bayes_factor.frmtmb_draws post_prob
###   post_prob.frmtmb_draws

### ** Examples

## No test: 
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(1)
  dd <- data.frame(x = rnorm(40))
  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
  # each refusal names its reason and the replacement
  try(reloo(ds))
  try(bayes_factor(ds, ds))
}
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
