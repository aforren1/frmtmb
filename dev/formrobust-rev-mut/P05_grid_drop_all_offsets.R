# The grid route drops the selected predictor's offset too.
source("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut/helper.R", local = TRUE)
mut("emm_drop_offsets", "setdiff(names(lps), keep)", "names(lps)")
