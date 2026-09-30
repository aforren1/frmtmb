.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
show <- function(nm) {
  cat("\n######## ", nm, "\n")
  print(get(nm, envir = ns))
}
for (nm in c("data_response.brmsframe", "get_rate_denom",
             "multiply_dpar_rate_denom", "log_lik_poisson",
             "log_lik_negbinomial", "log_lik_geometric",
             "posterior_epred_poisson", "posterior_epred_negbinomial",
             "posterior_epred_trunc_poisson",
             "posterior_predict_poisson", "posterior_predict_negbinomial",
             "stan_log_lik_poisson", "stan_log_lik_negbinomial",
             "stan_log_lik_multiply_rate_denom", "use_glm_primitive",
             "prepare_predictions_data", "terms_ad", "frame_resp",
             "has_subset", "subset_data", "data_sp", "predictor_sp",
             "prepare_predictions_sp", "rename_Ymi", "stan_sp",
             "prepare_predictions.brmsframe", "prepare_conditions",
             "na_omit", "valid_ad_terms", "ad_families")) {
  if (exists(nm, envir = ns)) show(nm) else cat("\n## missing:", nm, "\n")
}
