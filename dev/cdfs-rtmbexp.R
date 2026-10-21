.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
x <- getNamespaceExports("RTMB"); print(sort(grep("^[pdq]", x, value=TRUE)))
print(packageVersion("RTMB")); print(packageVersion("RTMBdist"))
y <- getNamespaceExports("RTMBdist"); print(sort(grep("^p", y, value=TRUE)))
