.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
x <- getNamespaceExports("RTMB")
print(sort(grep("logspace|CondExp|pmax|pmin|ifelse|max", x, value = TRUE)))
