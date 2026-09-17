.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.coupling))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: cross_wishart
### Title: The complex Wishart family for a pair of signals
### Aliases: cross_wishart

### ** Examples

set.seed(4)
src <- rnorm(4096)
xs <- frm_cross_spectrum(src + rnorm(4096), 0.9 * src + rnorm(4096),
                         segments = 16)
xs <- xs[xs$freq < 0.1, ]
fit <- frmtmb::frm(
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
             pow2 ~ 1, coh ~ 1, phase ~ 1),
  family = cross_wishart(), data = xs)
frm_coherence(fit)[1, ]



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
