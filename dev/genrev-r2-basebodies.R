# Which borrowed generics carried work in their body at the BASE commit?
# Any such work that moved into a method must reach every method,
# including frmtmb.sample's.
e <- new.env()
lazyLoad(file.path("C:/Users/adf44/source/r/rellib-r3", "frmtmb", "R", "frmtmb"), envir = e)
gens <- c("as_draws", "as_draws_array", "as_draws_df", "as_draws_list",
  "as_draws_matrix", "as_draws_rvars", "nchains", "ndraws", "niterations",
  "nvariables", "variables", "loo", "loo_compare", "waic", "bayes_R2",
  "prior_summary", "pp_check", "conditional_effects", "expose_functions",
  "hypothesis", "posterior_summary", "LOO", "WAIC", "ngrps", "fixef",
  "ranef", "VarCorr", "refit")
for (g in gens) {
  if (!exists(g, envir = e, inherits = FALSE)) { cat(g, ": absent\n"); next }
  b <- body(get(g, envir = e))
  if (is.call(b) && identical(b[[1L]], as.name("{")) && length(b) == 2L) b <- b[[2L]]
  bare <- is.call(b) && identical(b[[1L]], as.name("UseMethod"))
  if (!bare) cat(sprintf("%-20s WORK: %s\n", g, paste(deparse(b), collapse = " ")))
}
cat("GENREVDONE\n")
