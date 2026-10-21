.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (g in c("loo_R2","loo_epred","loo_linpred","loo_predict","loo_predictive_interval")) {
  f <- getExportedValue("rstantools", g); cat("rstantools::", g, ": "); print(args(f))
}
f <- loo::loo_model_weights; print(args(f)); print(args(loo:::loo_model_weights.default))
print(args(loo::E_loo)); print(args(loo:::E_loo.default)); print(args(loo:::E_loo.matrix))
print(args(loo:::loo.matrix)); print(args(loo:::psis.default))
ns <- asNamespace("brms")
for (g in c("loo_R2","loo_epred","loo_linpred","loo_predict","loo_predictive_interval","loo_model_weights")) {
  f <- getExportedValue("brms", g); cat("brms export", g, "env:", environmentName(environment(f)), "\n")
}
