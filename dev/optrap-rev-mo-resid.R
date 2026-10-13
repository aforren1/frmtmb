# Reviewer of lane optrap: dev/optrap-mo-study.R with dev/optrap-rev-resid.R
# applied first (the unit from the residual sd).
#   Rscript dev/optrap-rev-mo-resid.R optrap <seeds> <out.tsv>
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/optrap-rev-resid.R")
source("dev/optrap-mo-study.R")
