.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.sample))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-draws-refusals
### Title: Methods a ported brms script may call that frmtmb does not have
### Aliases: frmtmb-draws-refusals stancode stancode.frmtmb_draws standata
###   standata.frmtmb_draws expose_functions.frmtmb_draws plot.frmtmb_draws
###   update.frmtmb_draws restructure restructure.frmtmb_draws
###   posterior_samples posterior_samples.frmtmb_draws nsamples
###   nsamples.frmtmb_draws parnames parnames.frmtmb_draws

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
  try(stancode(ds))
  try(nsamples(ds))
}
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
