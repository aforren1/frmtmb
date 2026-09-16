## brms method bodies for the ten generics, so each positional slot's
## MEANING is read off brms rather than guessed.
.libPaths(c("C:/Users/adf44/source/r/brmsmatch-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
loadNamespace("brms")
for (nm in c("mcmc_plot.brmsfit", "posterior_predict.brmsfit",
             "predictive_error.brmsfit", "psis.brmsfit",
             "pp_mixture.brmsfit", "as.matrix.brmsfit",
             "as_draws_array.brmsfit", "log_lik.brmsfit")) {
  cat("\n======== brms:::", nm, "\n", sep = "")
  print(get(nm, envir = asNamespace("brms")))
}

## The two that were left alone: brms's own methods delegate to
## bayesplot on the stanfit, so matching bayesplot there IS matching
## brms. Printed rather than asserted from memory.
for (nm in c("nuts_params.brmsfit", "log_posterior.brmsfit")) {
  cat("\n======== brms:::", nm, "\n", sep = "")
  print(get(nm, envir = asNamespace("brms")))
}
