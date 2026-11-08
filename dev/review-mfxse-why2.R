# Review of lane mfxse: which marginaleffects 1.0.0 function swallows the
# warning get_predict() raises inside the comparison plan.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (nm in c("comparison_plan_predict", "prediction_plan_predict", "myTryCatch")) {
  if (exists(nm, ns)) { cat("=====", nm, "\n"); print(get(nm, ns)) }
}
cat("\n== functions with suppressWarnings or myTryCatch, plan/jacobian/comparison/slopes:\n")
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns); if (!is.function(f)) next
  if (!grepl("plan|jacobian|compar|slope|get_se|std_error", nm)) next
  t <- deparse(f)
  i <- grep("suppressWarnings|myTryCatch", t)
  if (length(i)) cat(nm, ":", trimws(t[i]), sep = "\n  ")
}
