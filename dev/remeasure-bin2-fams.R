# Lane remeasure, bin 2P groups G3 and G4: which brms family names of
# tests.log_lik.R and tests.posterior_predict.R frmtmb can answer, and
# under which dpar names. Read before writing dev/brmsport-bin2-adapters.R.
#
#   Rscript dev/remeasure-bin2-fams.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "\n")
fams <- c("gaussian", "student", "lognormal", "shifted_lognormal",
          "exgaussian", "skew_normal", "asym_laplace",
          "zero_inflated_asym_laplace", "binomial", "beta_binomial",
          "discrete_weibull", "poisson", "negbinomial", "negbinomial2",
          "geometric", "com_poisson", "exponential", "gamma", "Gamma",
          "frechet", "inverse.gaussian", "weibull", "gen_extreme_value",
          "bernoulli", "Beta", "beta", "von_mises", "hurdle_poisson",
          "hurdle_negbinomial", "hurdle_gamma", "zero_inflated_poisson",
          "zero_inflated_negbinomial", "zero_inflated_binomial",
          "zero_inflated_beta_binomial", "zero_inflated_beta",
          "zero_one_inflated_beta", "xbeta", "wiener")
ns <- asNamespace("frmtmb")
for (f in fams) {
  ex <- f %in% getNamespaceExports("frmtmb")
  fun <- if (exists(f, envir = ns)) get(f, envir = ns) else NULL
  if (!is.function(fun)) {
    cat(sprintf("%-28s exported=%s  NOT A FUNCTION IN frmtmb\n", f, ex))
    next
  }
  env <- environment(fun)
  owner <- if (is.null(env)) "?" else environmentName(env)
  fam <- tryCatch(fun(), error = function(e) e)
  if (inherits(fam, "error")) {
    cat(sprintf("%-28s exported=%s owner=%s ERROR %s\n", f, ex, owner,
                conditionMessage(fam)))
    next
  }
  cat(sprintf("%-28s exported=%s owner=%s class=%s family=%s\n", f, ex,
              owner, paste(class(fam), collapse = "/"),
              paste(fam$family, collapse = "/")))
  cat(sprintf("   dpars: %s  links: %s  lpdf=%s sim=%s lcdf=%s\n",
              paste(fam$dpars, collapse = ","),
              paste(vapply(fam$links, function(l) {
                if (is.character(l)) l else
                  paste(l$name %||% class(l)[1], collapse = "")
              }, ""), collapse = ","),
              is.function(fam$lpdf), is.function(fam$sim),
              is.function(fam$lcdf)))
  if (is.function(fam$lpdf)) {
    cat("   lpdf formals:", names(formals(fam$lpdf)), " sim formals:",
        if (is.function(fam$sim)) names(formals(fam$sim)), "\n")
  }
}
cat("\nbrmsfamily('frechet'):\n")
print(tryCatch(brmsfamily("frechet")$dpars, error = conditionMessage))
cat("stats::Gamma masked? ", identical(Gamma, stats::Gamma), "\n")
