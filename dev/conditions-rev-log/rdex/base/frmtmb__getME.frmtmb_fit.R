.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: getME.frmtmb_fit
### Title: Extract components of a fit, lme4 style
### Aliases: getME.frmtmb_fit

### ** Examples

if (requireNamespace("lme4", quietly = TRUE)) {
  set.seed(1)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  # the designs, for downstream code written against merMod objects
  dim(lme4::getME(fit, "X"))
  dim(lme4::getME(fit, "Zt"))

  # a vector of names returns a named list
  str(lme4::getME(fit, c("n_rtrms", "n_rfacs", "sigma")))

  # the conditional modes in coefficient space, aligned with Z
  head(lme4::getME(fit, "b"))
  # note: "lower" is all -Inf here, because the internal covariance
  # parameterization is unbounded. Use diagnose() to spot a singular
  # fit, not theta == lower.
  lme4::getME(fit, "lower")
}



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
