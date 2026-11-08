# Review of lane mfxse: marginaleffects 1.0.0 functions that catch warnings.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns); if (!is.function(f)) next
  t <- deparse(f)
  i <- grep("warning = function|suppressWarnings[(]", t)
  if (length(i) && any(grepl("get_predict|FUN|do_call|inner|func[(]", t))) cat(nm, "\n")
}
cat("== callers of myTryCatch or withCallingHandlers(warning):\n")
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns); if (!is.function(f)) next
  t <- paste(deparse(f), collapse = "\n")
  if (grepl("myTryCatch[(]|muffleWarning", t)) cat(nm, "\n")
}
