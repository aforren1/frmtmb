LIB <- "C:/Users/adf44/source/r/wt-records-lib"
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-records"
cat("roxygen2", format(packageVersion("roxygen2")), "\n")
roxygen2::roxygenise(WT)
cat("DONE roxygenise\n")
