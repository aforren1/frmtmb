# Punch round 1, m2: dev/gpby-brms-epred2.R (frmtmb's conditional law of
# the field at new positions against brms's, draw by draw, on brms's own
# draws, and the spread of the sampled field) with frmtmb's nugget set to
# the value given. The fit's estimates are replaced by each brms draw, so
# only the frame and the nugget come from frmtmb.
# Usage: Rscript dev/gpby-p1-nugget-brms.R <nugget>
nug <- as.numeric(commandArgs(TRUE)[1])
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
utils::assignInNamespace("gp_nugget", nug, "frmtmb")
cat("gp_nugget", asNamespace("frmtmb")$gp_nugget, "\n")
setwd("C:/Users/adf44/source/r/frmtmb-wt-gpby")
source("dev/gpby-brms-epred2.R", echo = FALSE)
