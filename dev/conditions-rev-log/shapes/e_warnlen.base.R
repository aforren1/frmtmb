.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

options(warning.length = 100)
f <- function() stop(strrep("x", 300))
f()
