source("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut/helper.R", local = TRUE)
old <- get("bernoulli_levels", asNamespace("frmtmb")); put("bernoulli_levels", function(y) rev(old(y)))
