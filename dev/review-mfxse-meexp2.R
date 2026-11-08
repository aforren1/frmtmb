# Review of lane mfxse: the analytic-Jacobian route of marginaleffects
# 1.0.0: get_model_matrix (exported), its attribute, and the analytic
# generic's methods.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (nm in c("get_model_matrix", "get_model_matrix.default",
             "add_model_matrix_attribute", "model_matrix_retry",
             "get_jacobian_analytic.default", "jacobian_analytic_scale")) {
  if (exists(nm, ns)) { cat("=====", nm, "\n"); print(get(nm, ns)) }
}
cat("analytic methods:", grep("^get_jacobian_analytic", ls(ns), value = TRUE), "\n")
