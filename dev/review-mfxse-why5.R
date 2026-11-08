# Review of lane mfxse: marginaleffects 1.0.0's get_predict_error() and its callers.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
print(get("get_predict_error", ns))
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns); if (!is.function(f)) next
  if (grepl("get_predict_error[(]", paste(deparse(f), collapse = "\n"))) cat("caller:", nm, "\n")
}
