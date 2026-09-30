.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("brms")
out <- "dev/formula2-brms-src"
fns <- c("check_fdpars", "expand_dot_formula", "terms_fe", "terms_re",
         "no_cmc", "validate_terms", "lf", "update.brmsfit",
         "update.brmsformula", "get_model_matrix", "frame_fe",
         "data_fe", "validate_family", "validate_formula.default",
         "validate_formula.mvbrmsformula", "brmsterms.mvbrmsformula",
         "mvbrmsformula", "plus_brmsformula", "+.bform", "terms_lf",
         "data_re", "frame_re", "get_re_terms", "stan_fe", "stan_mixture",
         "stan_dpar_defs", "fixef_predictors", "is_equal_formula",
         "combine_formulas", "stan_log_lik_mixture", "stan_predictor",
         "stan_dpars", "validate_resp_formula")
for (f in fns) {
  if (exists(f, envir = ns, inherits = FALSE)) {
    writeLines(deparse(get(f, envir = ns)), file.path(out, paste0(f, ".R")))
  } else {
    cat("missing:", f, "\n")
  }
}
all <- ls(ns, all.names = TRUE)
hit <- function(pat) {
  Filter(function(f) {
    o <- get(f, envir = ns)
    is.function(o) && any(grepl(pat, deparse(o)))
  }, all)
}
cat("pfix:", hit("pfix|fdpars"), "\n")
cat("mvbf family list:", hit("is.list\\(family\\)|family\\[\\[i\\]\\]"), "\n")
