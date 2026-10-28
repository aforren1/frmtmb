.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (f in c("get_vcov.glmmTMB", "get_vcov.brmsfit", "align_vcov_to_coef", "sanitize_model_specific", "sanitize_model_specific.default", "sanitize_model_specific.glmmTMB", "sanitize_model_specific.merMod","sanitize_model_specific.brmsfit", "sanitize_vcov_request", "get_vcov_label")) { cat("=====", f, "\n"); if (exists(f, ns)) print(get(f, ns)) }
