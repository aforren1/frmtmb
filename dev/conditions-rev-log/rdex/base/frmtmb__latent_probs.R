.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: latent_probs
### Title: Posterior probabilities of the latent states
### Aliases: latent_probs

### ** Examples

set.seed(1)
dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)))
fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
head(latent_probs(fit))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
