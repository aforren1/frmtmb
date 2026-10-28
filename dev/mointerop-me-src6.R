.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (f in c("get_jacobian", "get_se_delta", "get_jacobian_numderiv", "get_jacobian_fdforward")) { cat("=====", f, "\n"); if (exists(f, ns)) print(get(f, ns)) }
cat(grep("jacobian|set_coef", ls(ns), value = TRUE), sep = "\n")
