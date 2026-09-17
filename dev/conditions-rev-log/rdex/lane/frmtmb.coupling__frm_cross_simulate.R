.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.coupling))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_cross_simulate
### Title: Draw whole cross-spectral matrices from a fitted coupling model
### Aliases: frm_cross_simulate

### ** Examples

set.seed(6)
src <- rnorm(2048)
xs <- frm_cross_spectrum(src + rnorm(2048), 0.9 * src + rnorm(2048),
                         segments = 16)
xs <- xs[xs$freq < 0.1, ]
fit <- frmtmb::frm(
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
             pow2 ~ 1, coh ~ 1, phase ~ 1),
  family = cross_wishart(), data = xs)
str(frm_cross_simulate(fit, nsim = 2, seed = 1)[[1]])



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
