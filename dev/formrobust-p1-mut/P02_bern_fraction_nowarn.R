# m2's warning removed: two values inside (0, 1) are coded silently.
source("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut/helper.R", local = TRUE)
mut("extract_y", "all(uy > 0 & uy < 1)", "FALSE")
