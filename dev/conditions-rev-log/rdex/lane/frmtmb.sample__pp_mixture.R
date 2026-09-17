.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.sample))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: pp_mixture
### Title: Posterior mixture-component probabilities
### Aliases: pp_mixture pp_mixture.frmtmb_draws

### ** Examples

## No test: 
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(4)
  dd <- data.frame(y = c(rnorm(60, -2), rnorm(60, 3)))
  fit <- frm(bf(y ~ 1), family = frmtmb::mixture(gaussian(), gaussian()),
             data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
  head(pp_mixture(ds)[, "Estimate", ])
}
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
