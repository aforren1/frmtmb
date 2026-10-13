# Reviewer of lane optrap, re-check: run another script with the MAD
# patch (dev/optrap-rev2-mad.R) applied first, on the lane's build.
#   Rscript dev/optrap-rev2-with-mad.R <script> [its args ...]
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/optrap-rev2-mad.R")
args0 <- commandArgs(TRUE)
a <- args0[-1]
commandArgs <- function(trailingOnly = FALSE) a
source(args0[1])
