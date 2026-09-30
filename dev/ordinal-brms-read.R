.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
ns <- asNamespace("brms")
nm <- ls(ns, all.names = TRUE)
cat("== names with thres ==\n"); print(grep("thres", nm, value = TRUE))
cat("== names with disc ==\n"); print(grep("disc", nm, value = TRUE))
cat("== names with ordinal/ordered ==\n")
print(grep("ordinal|ordered", nm, value = TRUE))
show <- function(f) {
  cat("\n==== ", f, " ====\n", sep = "")
  print(get(f, ns))
}
for (f in c("cumulative", "sratio", "acat", "hurdle_cumulative",
            ".family_cumulative", "check_thres",
            "stan_thres", "has_thres", "has_ordered_thres",
            "has_equidistant_thres", "has_sum_to_zero_thres",
            "get_thres", "prior_thres", "stan_ordinal_lpmf",
            "use_ordered_builtin", "stan_log_lik_ordinal",
            "stan_log_lik_cumulative", "data_thres",
            "extract_thres_names", "brmsterms.brmsformula",
            "valid_dpars.default", "prior_predictor.btl",
            "def_dpar_prior", "family_bounds.brmsterms",
            "dpar_class", "get_ad_values", "fixed_dpars",
            "stan_predictor.brmsterms", "stan_fe")) {
  if (exists(f, ns, inherits = FALSE)) show(f) else cat("\n(no ", f, ")\n")
}
