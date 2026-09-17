.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_series_draw
### Title: Draw a time series from a fitted spectrum
### Aliases: frm_series_draw

### ** Examples

set.seed(1)
y <- as.numeric(arima.sim(list(ar = 0.6), 512))
pg <- frm_periodogram(y, fs = 128)
fit <- frm(bf(pgram ~ s(freq, k = 8)), family = whittle(), data = pg)

# a new series with the same estimated spectrum, not the same series
z <- frm_series_draw(fit, nsim = 2, seed = 1)
dim(z)
stats::frequency(z)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
