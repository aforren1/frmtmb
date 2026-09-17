.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

f <- function(x) warning("careful ", x)
withCallingHandlers(f(1), warning = function(w) { cat("got:", conditionMessage(w), deparse(conditionCall(w)), "\n"); invokeRestart("muffleWarning") })
