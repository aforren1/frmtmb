# Source the lane's two R files into the global environment for fast
# iteration, with the packages frm_ode() needs attached.
.libPaths(c("C:/Users/adf44/source/r/lincmt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB)
  library(frmtmb)
})
LINROOT <- "C:/Users/adf44/source/r/frmtmb-wt-lincmt/extensions/frmtmb.ode"
source(file.path(LINROOT, "R", "ode.R"))
source(file.path(LINROOT, "R", "lincmt.R"))
