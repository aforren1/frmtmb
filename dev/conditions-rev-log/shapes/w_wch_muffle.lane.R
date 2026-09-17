.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

f <- function(x) frm_warning("careful ", x)
withCallingHandlers(f(1), warning = function(w) { cat("got:", conditionMessage(w), deparse(conditionCall(w)), "\n"); invokeRestart("muffleWarning") })
