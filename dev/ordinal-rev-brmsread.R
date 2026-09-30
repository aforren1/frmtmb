# Reviewer, lane ordinal: print brms 2.23.0 functions the review reads.
# Output: dev/ordinal-rev-log-brmsread.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
cat("brms", format(packageVersion("brms")), "\n")
ns <- asNamespace("brms")
for (f in c("stan_thres", "prior_thres", "stan_ordinal_lpmf",
            "inv_link_acat", "inv_link", "dcumulative", "dacat", "dsratio",
            "dcratio", "posterior_epred_ordinal", "def_dpar_prior",
            "has_ordered_thres", "stan_center_X", "has_thres_groups",
            "has_sum_to_zero_thres", "has_equidistant_thres",
            "prepare_family", "validate_family", "is_equal",
            "log_lik_ordinal", "posterior_epred_cumulative")) {
  if (exists(f, envir = ns, inherits = FALSE)) {
    cat("\n########", f, "\n")
    print(get(f, envir = ns))
  } else cat("\n######## (missing)", f, "\n")
}
