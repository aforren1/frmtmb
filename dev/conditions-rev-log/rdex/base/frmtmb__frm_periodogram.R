.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_periodogram
### Title: Periodogram for a Whittle-likelihood fit
### Aliases: frm_periodogram

### ** Examples

set.seed(1)
y <- as.numeric(arima.sim(list(ar = 0.7), 512))
pg <- frm_periodogram(y, fs = 256)
str(pg)

# the ordinates average to the variance, in per-Hz units
c(mean(pg$pgram), var(y) / 256)

# one call for many series: columns are series
pg2 <- frm_periodogram(cbind(a = y, b = rev(y)), fs = 256)
table(pg2$series)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
