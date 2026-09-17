.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

f <- function(i) frm_warning("w", i)
g <- function() for (i in 1:12) f(i)
g()
print(warnings())
