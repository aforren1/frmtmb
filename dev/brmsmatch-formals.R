## Prints brms's method formals for the ten generics named in item 2.5f,
## so the positional slots are read off the real package rather than
## recalled. Run: Rscript dev/brmsmatch-formals.R
.libPaths(c("C:/Users/adf44/source/r/brmsmatch-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

show_fm <- function(pkg, nm) {
  ns <- asNamespace(pkg)
  if (!exists(nm, envir = ns, inherits = FALSE)) {
    cat(pkg, "::", nm, " MISSING\n", sep = "")
    return(invisible(NULL))
  }
  f <- get(nm, envir = ns)
  cat(pkg, "::", nm, "\n  ", sep = "")
  cat(paste(names(formals(f)), collapse = " "), "\n")
}

loadNamespace("brms")
cat("brms ", as.character(packageVersion("brms")), "\n", sep = "")
for (nm in c("as.mcmc.brmsfit", "log_lik.brmsfit", "mcmc_plot.brmsfit",
             "posterior_epred.brmsfit", "posterior_interval.brmsfit",
             "posterior_linpred.brmsfit", "posterior_predict.brmsfit",
             "pp_mixture.brmsfit", "predictive_error.brmsfit",
             "psis.brmsfit", "rhat.brmsfit", "neff_ratio.brmsfit",
             "nuts_params.brmsfit", "log_posterior.brmsfit",
             "predictive_interval.brmsfit", "pp_check.brmsfit")) {
  show_fm("brms", nm)
}

cat("\n-- generics (owners) --\n")
for (p in c("coda", "rstantools", "brms", "loo", "bayesplot",
            "posterior")) {
  loadNamespace(p)
}
show_fm("coda", "as.mcmc")
show_fm("rstantools", "log_lik")
show_fm("brms", "mcmc_plot")
show_fm("rstantools", "posterior_epred")
show_fm("rstantools", "posterior_interval")
show_fm("rstantools", "posterior_linpred")
show_fm("rstantools", "posterior_predict")
show_fm("brms", "pp_mixture")
show_fm("rstantools", "predictive_error")
show_fm("loo", "psis")
show_fm("posterior", "rhat")
show_fm("bayesplot", "rhat")
show_fm("bayesplot", "neff_ratio")

cat("\n-- bodies --\n")
print(brms:::rhat.brmsfit)
print(brms:::neff_ratio.brmsfit)
print(brms:::as.mcmc.brmsfit)
print(brms:::posterior_interval.brmsfit)
