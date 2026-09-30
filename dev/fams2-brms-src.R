# Dump brms 2.23.0 source for the families this lane adds.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(brms)
cat("brms", format(packageVersion("brms")), "\n")
ns <- asNamespace("brms")
show <- function(nm) {
  cat("\n=====", nm, "=====\n")
  if (exists(nm, ns)) print(get(nm, ns)) else cat("(absent)\n")
}
for (nm in c("xbeta", "zero_inflated_beta_binomial", "hurdle_cumulative",
             "beta_binomial", "sratio", "cumulative", "cse", "cs",
             "posterior_epred_xbeta", "posterior_predict_xbeta",
             "log_lik_xbeta", "posterior_epred_zero_inflated_beta_binomial",
             "posterior_predict_zero_inflated_beta_binomial",
             "log_lik_zero_inflated_beta_binomial",
             "posterior_epred_beta_binomial", "posterior_predict_beta_binomial",
             "log_lik_beta_binomial",
             "posterior_epred_hurdle_cumulative",
             "posterior_predict_hurdle_cumulative",
             "log_lik_hurdle_cumulative", "posterior_epred_ordinal",
             "posterior_epred_cumulative",
             "dhurdle_cumulative", "dcumulative", "rxbeta", "dxbeta",
             "posterior_predict_ordinal", "log_lik_ordinal",
             "log_lik_cumulative", "posterior_predict_cumulative",
             "dbeta_binomial", "rbeta_binomial")) show(nm)
fams <- ls(ns, pattern = "^\\.family_")
cat("\nfamily info fns:\n"); print(fams)
for (nm in c(".family_xbeta", ".family_zero_inflated_beta_binomial",
             ".family_hurdle_cumulative", ".family_beta_binomial",
             ".family_cumulative", ".family_sratio")) show(nm)
fdir <- system.file("chunks", package = "brms")
for (f in c("fun_xbeta.stan", "fun_zero_inflated_beta_binomial.stan",
            "fun_hurdle_cumulative.stan", "fun_beta_binomial.stan",
            "fun_cumulative_logit.stan")) {
  cat("\n=====", f, "=====\n")
  p <- file.path(fdir, f)
  if (file.exists(p)) cat(readLines(p), sep = "\n") else cat("(absent)\n")
}
cat("\nchunk files matching:\n")
print(grep("xbeta|beta_binomial|hurdle_cum|cumulative", list.files(fdir),
           value = TRUE))
