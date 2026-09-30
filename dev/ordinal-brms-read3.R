.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
for (f in c("posterior_linpred.brmsfit", "posterior_linpred.mvbrmsprep",
            "posterior_linpred.brmsprep", "fitted.brmsfit",
            "get_dpar", "predictor.bprepl", "posterior_epred_hurdle_cumulative",
            "posterior_predict_hurdle_cumulative", "stan_log_lik_hurdle_cumulative",
            "summary.brmsfit")) {
  cat("\n==== ", f, " ====\n", sep = "")
  if (exists(f, ns, inherits = FALSE)) print(get(f, ns)) else cat("(none)\n")
}
