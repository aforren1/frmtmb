.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (f in c("get_vcov.default", "sanitize_vcov", "get_vcov")) { cat("=====", f, "\n"); if (exists(f, ns)) print(get(f, ns)) }
cat(grep("vcov|sanitize_model_specific|supported|re_formula", ls(ns), value = TRUE), sep = "\n")
