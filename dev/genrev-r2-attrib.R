# genrev round 2: is "no owner loaded, frmtmb_fit methods not frmtmb's:
# SWAP 10/26, ACTIVE 0/26" a win of the ACTIVE BINDING, or of the ten
# refusal methods the same round added?  Take the SWAP build, register
# ten refusals in-session exactly where the ACTIVE build's NAMESPACE
# puts them (frmtmb's table and posterior's), and recount with the
# lane's own list and rule: the method that would run must be a
# frmtmb_fit method owned by frmtmb.
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]; patch <- identical(a[2], "patch")
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
ten <- c("as_draws", "as_draws_df", "as_draws_array", "as_draws_list",
         "as_draws_matrix", "as_draws_rvars", "ndraws", "nchains",
         "niterations", "nvariables")
if (patch) for (g in ten) {
  f <- eval(bquote(function(x, ...) stop(.(g), "() refused")))
  environment(f) <- ns
  registerS3method(g, "frmtmb_fit", f, envir = ns)
}
own <- c("fixef", "ranef", "VarCorr", "ngrps", "prior_summary", "loo",
         "waic", "bayes_R2", "LOO", "WAIC", "hypothesis",
         "conditional_effects", "expose_functions", "variables", ten,
         "pp_check", "refit")
bad <- character()
for (g in own) {
  f <- get(g, envir = globalenv())
  tb <- get(".__S3MethodsTable__.", envir = environment(f))
  nm <- paste0(g, ".frmtmb_fit")
  ok <- exists(nm, envir = tb, inherits = FALSE) &&
    environmentName(topenv(environment(get(nm, envir = tb)))) == "frmtmb"
  if (!ok) bad <- c(bad, g)
}
cat(sprintf("LIB %s patch %s owners loaded %s\n", basename(LIB), patch,
            paste(intersect(c("posterior", "loo", "brms", "bayesplot",
                              "rstantools", "lme4"), loadedNamespaces()),
                  collapse = ",")))
cat(sprintf("FITMETHOD_NOT_OURS %d of %d  %s\n", length(bad), length(own),
            paste(bad, collapse = ",")))
cat("GENREVDONE\n")
