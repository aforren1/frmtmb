# Lane ceplot punch 1: the deparsed text of two installed functions, to
# write mutant patterns that match.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat(deparse(frmtmb:::ce_model_vars), sep = "\n")
cat(deparse(frmtmb.sample:::ps_select), sep = "\n")
