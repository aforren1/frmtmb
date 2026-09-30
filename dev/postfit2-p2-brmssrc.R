.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("brms")
hits <- Filter(function(n) { f <- get(n, ns); is.function(f) && any(grepl("get_new_rdraws|new_levels", deparse(f))) }, ls(ns, all.names = TRUE))
print(hits)
for (fn in c("prepare_predictions_re", "prepare_predictions_re_global", "prepare_Z_mm", "validate_newdata")) {
  if (exists(fn, ns)) { cat("=====", fn, "\n"); print(get(fn, ns)) }
}
