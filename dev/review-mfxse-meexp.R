# Review of lane mfxse: which marginaleffects 1.0.0 extension hooks are
# exported, and who sets the analytic Jacobian's model-matrix attribute.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ex <- getNamespaceExports("marginaleffects")
cat("marginaleffects", format(packageVersion("marginaleffects")), "\n")
print(sort(grep("^get_|^set_|jacobian|model_matrix", ex, value = TRUE)))
ns <- asNamespace("marginaleffects")
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns); if (!is.function(f)) next
  txt <- paste(deparse(f), collapse = "\n")
  if (grepl("marginaleffects_model_matrix|get_model_matrix", txt)) cat("uses mm:", nm, "\n")
}
print(methods("get_jacobian_analytic"))
print(methods("get_model_matrix"))
