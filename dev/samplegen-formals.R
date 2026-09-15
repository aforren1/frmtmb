# The full formals, defaults included, of every owner's generic for
# the 28 names, beside brms's method for the same name, and the
# versions of the owners. brms is the tiebreaker where two disagree.
#
#   Rscript dev/samplegen-formals.R
RR <- "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(RR, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
tab <- list(
  as.mcmc = "coda", bayes_factor = "bridgesampling",
  bridge_sampler = "bridgesampling", post_prob = "bridgesampling",
  kfold = "loo", loo_moment_match = "loo", loo_subsample = "loo",
  psis = "loo", log_lik = "rstantools", nsamples = "rstantools",
  posterior_epred = "rstantools", posterior_interval = "rstantools",
  posterior_linpred = "rstantools", posterior_predict = "rstantools",
  predictive_error = "rstantools", predictive_interval = "rstantools",
  log_posterior = "bayesplot", neff_ratio = "bayesplot",
  nuts_params = "bayesplot", rhat = c("posterior", "bayesplot"),
  mcmc_plot = "brms", parnames = "brms",
  posterior_samples = c("brms", "gratia"), pp_mixture = "brms",
  reloo = "brms", restructure = "brms", stancode = "brms",
  standata = "brms")
suppressMessages(loadNamespace("brms"))
suppressMessages(loadNamespace("frmtmb.sample"))
fm <- function(f) {
  a <- formals(f)
  paste(vapply(names(a), function(n) {
    d <- paste(deparse(a[[n]]), collapse = "")
    if (!nzchar(d)) n else paste0(n, " = ", d)
  }, ""), collapse = ", ")
}
for (g in names(tab)) {
  cat(sprintf("\n%s\n", g))
  cat(sprintf("  %-15s (%s)\n", "frmtmb.sample",
              fm(get(g, envir = asNamespace("frmtmb.sample")))))
  m <- get0(paste0(g, ".frmtmb_draws"),
            envir = asNamespace("frmtmb.sample"))
  cat(sprintf("  %-15s (%s)\n", "  its method", fm(m)))
  for (p in tab[[g]]) {
    suppressMessages(loadNamespace(p))
    cat(sprintf("  %-15s (%s)  body: %s\n", p, fm(getExportedValue(p, g)),
                paste(deparse(body(getExportedValue(p, g))), collapse = " ")))
  }
  bm <- getS3method(g, "brmsfit", optional = TRUE,
                    envir = asNamespace("brms"))
  if (is.null(bm)) bm <- get0(paste0(g, ".brmsfit"), envir = asNamespace("brms"))
  cat(sprintf("  %-15s (%s)\n", "brmsfit method",
              if (is.null(bm)) "(none)" else fm(bm)))
}
cat("\n")
for (p in c("coda", "bridgesampling", "loo", "rstantools", "bayesplot",
            "posterior", "brms", "gratia")) {
  cat(sprintf("%-15s %s\n", p, as.character(packageVersion(p))))
}
