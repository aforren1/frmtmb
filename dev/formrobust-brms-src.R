# Print brms 2.23.0 functions this lane reads, to
# dev/formrobust-brms-src/<name>.R
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("brms")
cat("brms", format(packageVersion("brms")), "\n")
out <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-brms-src"
dir.create(out, showWarnings = FALSE)
fns <- c("data_response", "data_response.brmsterms", "data_response.mvbrmsterms",
         "get_ad_values", "get_ad_vars", "get_ad_expr", "eval_rhs",
         "resp_weights", "resp_rate", "resp_trials", "resp_se", "resp_cens",
         "resp_trunc", "resp_subset", "resp_index", "resp_mi", "resp_thres",
         "resp_cat", "resp_dec", "resp_vint", "resp_vreal",
         "validate_data", "subset_data", "brmsterms.brmsformula",
         "terms_ad", "all_vars", "get_all_vars", "validate_newdata",
         "prepare_predictions.brmsfit", "prepare_predictions_ac",
         "prepare_predictions.bframel", "get_ac_params",
         "predictive_error.brmsfit", "posterior_predict.brmsfit",
         "posterior_predict_ac", "predictor.bprepl", "predictor_ac",
         "data_ac", "check_brm_ac", "validate_resp", "check_ac",
         "arma", "ar", "ma", "cosy", "autocor", "autocor.brmsfit",
         "`+.bform`", "plus_brmsformula", "validate_formula",
         "validate_formula.default", "validate_formula.brmsformula",
         "update.brmsfit", "update.brmsformula", "update_re_terms",
         "bf", "brmsformula", "prepare_family", "check_cor_by",
         "data_fe", "fixef_predictors", "rm_intercept",
         "validate_fixef", "frame_fe", "get_model_matrix",
         "check_reserved_vars", "rename_intercept",
         "standata.default", "make_standata", "standata",
         "drop_unused_factor_levels", "validate_data2",
         "extract_offset", "get_offset", "offset_data",
         "conditional_effects.brmsfit", "prepare_conditions",
         "make_conditions", "prepare_marg_data", "get_int_vars",
         "get_all_effects", "emmeans.brmsfit", "emm_basis.brmsfit",
         "recover_data.brmsfit", "prepare_predictions_offset",
         "prepare_predictions_fe", "is_like_factor", "as_factor",
         "validate_y", "exclude_pars", "is.cor_arma", "cor_arma",
         "as_cor_brms", "brms_to_cor", "autocor", "ac_term_data",
         "get_ac_vars", "get_ac_obj", "tidy_acef", "data_gr_global",
         "resp_bernoulli", "family_bounds", "check_regex",
         "add_ac_terms", "terms_ac")
for (f in fns) {
  f0 <- gsub("`", "", f)
  if (!exists(f0, envir = ns, inherits = FALSE)) {
    cat("MISSING", f0, "\n"); next
  }
  obj <- get(f0, envir = ns)
  writeLines(c(paste0(f0, " <- "), deparse(obj, control = "useSource")),
             file.path(out, paste0(gsub("[^A-Za-z0-9._]", "_", f0), ".R")))
}
cat(paste(ls(ns, pattern = "^resp_|^data_response|offset|^ac_|arma|autocor|intercept", all.names = TRUE), collapse = "\n"))
