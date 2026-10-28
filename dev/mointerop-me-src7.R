.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (f in c("get_delta_jacobian", "get_jacobian_fd", "std_error_from_jacobian")) { cat("=====", f, "\n"); if (exists(f, ns)) print(get(f, ns)) }
