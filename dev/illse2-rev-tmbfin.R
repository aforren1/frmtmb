# Re-check of lane illse2, item 6: TMB's ADFun finalizer and the lazy
# load of TMB:::isNullPointer.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("TMB")
print(packageVersion("TMB")); print(find.package("TMB"))
print(get("registerFinalizer", ns))
print(get("isNullPointer", ns))
print(packageVersion("RTMB"))
cat(grep("registerFinalizer|isNullPointer|finaliz", deparse(get("MakeADFun", ns)), value = TRUE), sep = "\n")
print(grep("inaliz|NullPointer", ls(asNamespace("RTMB"), all.names = TRUE), value = TRUE))
print(grep("TMB", c(packageDescription("RTMB")$Imports, packageDescription("RTMB")$Depends), value = TRUE))
