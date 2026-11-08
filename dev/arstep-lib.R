# Lane arstep: the review's helpers (dev/review-ncgp-lib.R) pointed at this
# round's libraries. "lane" is wt-arstep-lib over rellib-r18; "base" is
# rellib-r18 alone (frmtmb 0.79.1 at 3cd2a0c8).
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-ncgp-lib.R")
REV_LANE <- "C:/Users/adf44/source/r/wt-arstep-lib"
REV_BASE <- "C:/Users/adf44/source/r/rellib-r18"
REV_WT <- "C:/Users/adf44/source/r/frmtmb-wt-arstep"
REV_OUT <- file.path(REV_WT, "dev/arstep-out")
dir.create(REV_OUT, showWarnings = FALSE)
