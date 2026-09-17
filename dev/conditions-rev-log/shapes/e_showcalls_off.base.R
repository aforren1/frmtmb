.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

options(showErrorCalls = FALSE)
h <- function() stop("deep")
g <- function() h()
g()
