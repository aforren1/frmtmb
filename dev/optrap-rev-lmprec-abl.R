# Reviewer of lane optrap: dev/optrap-rev-lmprec.R on the lane's build
# with the coefficient units taken out (predictor_par_units() NULL), to
# tell whether they cause the precision change.
#   Rscript dev/optrap-rev-lmprec-abl.R
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
utils::assignInNamespace("predictor_par_units", function(...) NULL, "frmtmb")
a <- c("C:/Users/adf44/source/r/wt-optrap-lib", "lane-nounits")
commandArgs <- function(trailingOnly = FALSE) a
source("dev/optrap-rev-lmprec.R")
