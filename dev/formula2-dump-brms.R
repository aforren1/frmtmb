.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("brms")
out <- "dev/formula2-brms-src"
dir.create(out, showWarnings = FALSE)
fns <- c("brmsformula", "validate_formula", "validate_formula.default",
         "validate_formula.brmsformula", "validate_formula.mvbrmsformula",
         "brmsterms", "brmsterms.brmsformula", "brmsterms.mvbrmsformula",
         "lf", "set_nl", "validate_par_formula", "lhs", "dpar_class",
         "prepare_auxformula", "brmsterms.default", "terms_ad",
         "validate_family", "is_nlpar", "check_fixed_dpars",
         "get_cmc", "prepare_mixture", "validate_data", "expand_dot_formula",
         "update.brmsfit", "update.brmsformula", "update_re_terms",
         "mvbrmsformula", "mvbf", "plus_brmsformula", "+.bform",
         "bf_to_mvbf", "validate_formula.default",
         "stan_log_lik_mixture", "data_fe", "frame_fe", "frame_basis",
         "prepare_design_matrix", "get_model_matrix", "fixef_predictors",
         "cs_predictors", "rm_int_fe", "check_fdpars", "mixture",
         "dpar_family", "get_fixed_dpars", "is_equal_dpar")
for (f in fns) {
  if (exists(f, envir = ns, inherits = FALSE)) {
    writeLines(deparse(get(f, envir = ns)), file.path(out, paste0(f, ".R")))
  } else {
    cat("missing:", f, "\n")
  }
}
# find every function whose body mentions cmc, or equating
all <- ls(ns, all.names = TRUE)
hit <- function(pat) {
  Filter(function(f) {
    o <- get(f, envir = ns)
    is.function(o) && any(grepl(pat, deparse(o)))
  }, all)
}
cat("cmc:", hit("cmc"), "\n")
cat("expand_dot:", hit("expand_dot"), "\n")
cat("equate:", hit("[Ee]quat"), "\n")
cat("pforms:", head(hit("pfix"), 40), "\n")
