# Short check of lane quadseam2's fix round: dev/quadseam2-hmm.R (item
# 5222, hmm and mixture seeds 1 to 3 started at alpha = 0) on the lane
# library. Git Bash mangles a path passed through cmd, so the library is
# given here instead.
#   Rscript quadseam2-short-hmm.R
commandArgs <- function(trailingOnly = FALSE) {
  "C:/Users/adf44/source/r/wt-quadseam2-lib"
}
source("C:/Users/adf44/source/r/frmtmb-wt-quadseam2/dev/quadseam2-hmm.R")
