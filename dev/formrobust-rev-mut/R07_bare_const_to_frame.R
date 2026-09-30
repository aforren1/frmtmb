source("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut/helper.R", local = TRUE)
mut("assemble_frame", "if (is.name(a)) {\n for (v in expr_frame_vars(a, data, resp$formula_env)) {\n add_part(as.name(v))\n }", "if (is.name(a)) {\n add_part(a)")
