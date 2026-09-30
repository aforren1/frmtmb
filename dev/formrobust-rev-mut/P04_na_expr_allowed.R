# m5's refusal removed: an NA from an expression reaches the fit.
source("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut/helper.R", local = TRUE)
mut("assemble_frame", "if (anyNA(v)) {", "if (FALSE) {")
