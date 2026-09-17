.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.sample))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: sample-as_draws
### Title: Convert draws to a posterior draws object
### Aliases: sample-as_draws as_draws.frmtmb_draws
###   as.data.frame.frmtmb_draws as.array.frmtmb_draws
###   as_draws_matrix.frmtmb_draws as_draws_array.frmtmb_draws
###   as_draws_df.frmtmb_draws as_draws_list.frmtmb_draws
###   as_draws_rvars.frmtmb_draws as.mcmc as.mcmc.frmtmb_draws
###   as.matrix.frmtmb_draws

### ** Examples

## No test: 
if (requireNamespace("posterior", quietly = TRUE) &&
    requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(9)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

  # hands the draws to the posterior package, keeping the frmtmb
  # parameter names
  dm <- as_draws(ds)
  posterior::summarise_draws(dm)
  # which is what variables() lists
  head(variables(ds))
}
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
