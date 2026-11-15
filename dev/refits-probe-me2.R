.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
print(marginaleffects:::get_modeldata)
ns <- asNamespace("marginaleffects")
cat("\nexported generics:\n")
print(grep("^(get_|set_|sanitize)", getNamespaceExports(ns), value = TRUE))
cat("\nmerge_original_data:\n")
print(marginaleffects:::merge_original_data)
