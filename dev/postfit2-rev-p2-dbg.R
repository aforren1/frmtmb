.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib","C:/Users/adf44/source/r/rellib-r3","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
tx <- deparse(get("ce_plan_part", asNamespace("frmtmb")))
print(grep("unique[(]vals|match[(]vals|new:", tx, value = TRUE))
