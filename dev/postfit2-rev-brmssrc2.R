.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("brms")
for (n in c("combine_prefix", "check_prefix", "conditional_smooths.brmsfit", "conditional_smooths.mvbrmsterms", "posterior_smooths.btl", "fill_newdata")) {
  cat("=====", n, "\n"); if (exists(n, envir = ns)) print(get(n, envir = ns)) else cat("absent\n")
}
