.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

setGeneric("gen4", function(x) standardGeneric("gen4"))
setMethod("gen4", "numeric", function(x) x)
h <- function() frm_stop("in arg")
gen4(h())
