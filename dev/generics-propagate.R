# Does core's fix reach frmtmb.sample's RE-EXPORTS?
#
# frmtmb.sample does `importFrom(frmtmb, as_draws, ...)` and exports
# the result. With the round-1 binding SWAP that re-export was a COPY
# taken at frmtmb.sample's load, so a later rebinding never reached it.
# With an active binding the active-ness propagates through
# `importIntoEnv()`, so the re-exported binding is active too and
# resolves on every access.
#
# The instrument is the method-table lookup, not an error string, for
# the reason given at the top of dev/generics-check.R.
#
#   Rscript dev/generics-propagate.R <LIB> <mode>
#     S  library(brms); library(frmtmb.sample)
#     T  library(frmtmb.sample); library(brms)
#     U  library(frmtmb.sample) with brms merely LOADED, which is what
#        happens the moment any other attached package Imports brms.
#        Review found this order missing from round 1.
av <- commandArgs(trailingOnly = TRUE)
.libPaths(c(av[1], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mode <- av[2]
if (mode == "S") {
  suppressMessages(library(brms))
  suppressMessages(library(frmtmb.sample))
} else if (mode == "T") {
  suppressMessages(library(frmtmb.sample))
  suppressMessages(library(brms))
} else {
  suppressMessages(library(frmtmb.sample))
  suppressMessages(loadNamespace("brms"))
}

reexported <- c("as_draws", "as_draws_array", "as_draws_df",
                "as_draws_list", "as_draws_matrix", "as_draws_rvars",
                "bayes_R2", "conditional_effects", "fixef", "hypothesis",
                "loo", "LOO", "loo_compare", "nchains", "ndraws",
                "ngrps", "niterations", "nvariables",
                "posterior_summary", "pp_check", "prior_summary",
                "ranef", "VarCorr", "variables", "waic", "WAIC")
own_to_sample <- c("log_lik", "posterior_epred", "posterior_linpred",
                   "posterior_predict", "posterior_interval",
                   "predictive_interval", "predictive_error", "kfold",
                   "psis", "reloo", "bridge_sampler", "bayes_factor",
                   "post_prob", "nuts_params", "log_posterior", "rhat",
                   "neff_ratio", "as.mcmc", "mcmc_plot", "stancode",
                   "standata", "restructure", "nsamples", "parnames",
                   "posterior_samples", "pp_mixture", "loo_moment_match",
                   "loo_subsample")

reach <- function(gen, k) {
  f <- tryCatch(get(gen), error = function(e) NULL)
  if (!is.function(f)) return(FALSE)
  e <- environment(f)
  if (is.null(e)) e <- baseenv()
  if (!exists(".__S3MethodsTable__.", envir = e, inherits = FALSE)) {
    return(FALSE)
  }
  tb <- get(".__S3MethodsTable__.", envir = e, inherits = FALSE)
  exists(paste0(gen, ".", k), envir = tb, inherits = FALSE)
}
lost <- function(gens) gens[!vapply(gens, reach, NA, "brmsfit")]

a <- lost(reexported)
b <- lost(own_to_sample)
cat(sprintf("MODE %s\n", mode))
cat(sprintf("re-exported from core      lost %d of %d\n",
            length(a), length(reexported)))
# which of the lost ones are answered by SOMEBODY'S default rather
# than by nothing: a count that only looks for "no method at all"
# will not count these, which is the method difference that makes
# the review's 22 and this script's 23 disagree by one
dflt <- a[vapply(a, reach, NA, "default")]
cat(sprintf("  of which answered by a .default  %d  %s\n",
            length(dflt), paste(dflt, collapse = ",")))
cat(sprintf("defined by frmtmb.sample   lost %d of %d\n",
            length(b), length(own_to_sample)))
act <- function(nm, e) tryCatch(bindingIsActive(nm, e),
                                error = function(x) NA)
cat(sprintf("as_draws binding active in package:frmtmb.sample  %s\n",
            act("as_draws", as.environment("package:frmtmb.sample"))))
cat(sprintf("as_draws resolves to       %s\n",
            environmentName(topenv(environment(get("as_draws"))))))
cat("DONE\n")
