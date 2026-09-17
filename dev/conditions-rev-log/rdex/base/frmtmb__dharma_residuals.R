.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: dharma_residuals
### Title: DHARMa residual diagnostics
### Aliases: dharma_residuals

### ** Examples

if (requireNamespace("DHARMa", quietly = TRUE)) {
  set.seed(1)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

  # scaled quantile residuals: uniform under a correct model, which
  # is what Pearson residuals cannot give for a discrete family
  res <- dharma_residuals(fit, nsim = 100, seed = 1)
  plot(res)
  DHARMa::testUniformity(res, plot = FALSE)
  DHARMa::testDispersion(res, plot = FALSE)
}



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
