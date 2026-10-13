# Reviewer of lane optrap, re-check: dev/optrap-rev2-b1.R without the least-squares cap
# patch (dev/optrap-rev2-nocap.R) applied first.
#   Rscript dev/optrap-rev2-b1-mad.R <case> [seeds]
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/optrap-rev2-nocap.R")
args0 <- commandArgs(TRUE)
a <- c("C:/Users/adf44/source/r/wt-optrap-lib", args0[1], "lanenocap", args0[-1])
commandArgs <- function(trailingOnly = FALSE) a
source("dev/optrap-rev2-b1.R")
