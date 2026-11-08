# Review of lane latperf: the deparsed text of a mutant site.
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
x <- deparse(frmtmb:::jcs_rowvar_inner, control = "useSource")
cat(grep("e1|e2", x, value = TRUE), sep = "\n")
