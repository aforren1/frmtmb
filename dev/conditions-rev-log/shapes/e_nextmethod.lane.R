.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

gen <- function(x) UseMethod("gen")
gen.a <- function(x) NextMethod()
gen.default <- function(x) frm_stop("s3n")
gen(structure(1, class = "a"))
