.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib","C:/Users/adf44/source/r/rellib-r4","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
x <- deparse(frmtmb:::ce_plan_part); writeLines(grep("label|relabel", x, value = TRUE))
writeLines(grep("pick", deparse(frmtmb:::fitted_old_levels), value = TRUE))