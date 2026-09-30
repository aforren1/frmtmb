# m1's refusal removed: refit() takes any bernoulli value.
source("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut/helper.R", local = TRUE)
mut("refit.frmtmb_fit", "!is.null(lv) && any(bad)", "FALSE")
