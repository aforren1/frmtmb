.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: response_mean
### Title: The two fitting options every mixture-type family refuses, in
###   its own name.
### Aliases: response_mean mixture_multimodal_refusals mixture_posterior
###   as_frmtmb_family eval_dpars single_response fit_extras dpar_linpred
###   structure_supports_all frmtmb-extension-api frame_block_of

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(50))
dd$y <- rnorm(50, 1 + 2 * dd$x, 0.5)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)

rspec <- single_response(fit, "my_family_probs()")
dp <- eval_dpars(fit)[[rspec$resp_name]]
str(dp)
head(response_mean(rspec$family, dp, list()))
fit_extras(fit)                       # NULL: no extra parameters
head(dpar_linpred(fit$frame, fit$estimates, rspec$resp_name, "mu"))
as_frmtmb_family(gaussian)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
