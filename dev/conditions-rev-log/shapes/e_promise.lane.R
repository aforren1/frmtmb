.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

g <- function(a) a
f <- function() g(frm_stop("p"))
f()
