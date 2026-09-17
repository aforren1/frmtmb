.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

f <- function(x) frm_warning("careful ", x)
v <- suppressWarnings(f(1))
print(v)
