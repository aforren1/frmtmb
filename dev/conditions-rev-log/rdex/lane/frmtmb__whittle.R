.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: whittle
### Title: Whittle likelihood for a periodogram response
### Aliases: whittle

### ** Examples

set.seed(1)
y <- as.numeric(arima.sim(list(ar = 0.6), 512))
pg <- frm_periodogram(y)

# a nonparametric log spectrum: the smoothing parameter is estimated
# by the same Laplace marginal likelihood as any other smooth
fit <- frm(bf(pgram ~ s(freq, k = 8)), family = whittle(), data = pg)
fixef(fit)

# an averaged periodogram needs the shape it was averaged with
pg4 <- frm_periodogram(y, segments = 4)
frm(bf(pgram ~ log(freq)), family = whittle(tapers = 4), data = pg4)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
