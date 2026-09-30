.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
show <- function(f) {
  cat("\n==== ", f, " ====\n", sep = "")
  if (exists(f, ns, inherits = FALSE)) print(get(f, ns)) else cat("(none)\n")
}
for (f in c(".brmsfamily", "fix_intercepts", "stan_center_X",
            "prepare_predictions_thres", "posterior_epred_ordinal",
            "rename_thres", "subset_thres", "has_extra_cat",
            "prepare_predictions.bframel", "dpar_family.default",
            "stan_hurdle_ordinal_lpmf", "has_cs", "check_cs",
            "needs_ordered_cs", "has_thres_minus_eta",
            "has_eta_minus_thres", "log_lik_cumulative", "log_lik_sratio",
            "log_lik_acat", "dcumulative", "dsratio", "dacat",
            "inv_link_cumulative", "inv_link_acat",
            "posterior_predict_ordinal", "pordinal",
            "check_fdpars", "summarise_families", "print.brmsfit",
            "change_effects.btl", "change_fe", "prior_fe",
            "prior_predictor.bframel", "prior_predictor.default",
            "def_scale_prior.brmsterms", "validate_family",
            "check_family", "family_info.brmsfamily",
            "family_info.mixfamily", "prior_Intercept")) show(f)
