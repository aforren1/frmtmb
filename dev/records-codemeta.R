LIB <- "C:/Users/adf44/source/r/wt-records-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-records"
cat("codemetar", format(packageVersion("codemetar")), "\n")
codemetar::write_codemeta(WT)
cat("DONE\n")
