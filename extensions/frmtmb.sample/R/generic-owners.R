# Why this file exists.
#
# This package defines 28 S3 generics of its own, and every one of
# those names is owned by another package: rstantools, loo,
# bridgesampling, bayesplot, posterior, coda, gratia or brms. A generic
# defined here and exported is a rival, not an alias. `UseMethod()`
# reads the method table of the namespace where the generic IT REACHED
# was defined, so after `library(brms); library(frmtmb.sample)` a
# `brmsfit` handed to `log_lik()` went to this package's table, which
# has no method for it, while brms's method sat registered in
# rstantools' table the whole time. Measured on the base commit with
# dev/samplegen-check.R: 28 of 28 lost in that order, 28 of 28 with
# brms merely LOADED, and a gratia user lost `posterior_samples()` for
# every class gratia has a method for.
#
# The `S3method(owner::generic, frmtmb_draws)` directives this package
# already carried did not prevent any of it. They put this package's
# method in the owner's table, which helps the owner's generic, and do
# nothing about the rival one here.
#
# The repair is frmtmb's own, called with this package's table:
# `frmtmb::frm_install_generics()` turns each name into an active
# binding that returns the owner's generic while the owner is loaded,
# and this package's generic while it is not. One implementation, two
# callers; see R/generic-owners.R in frmtmb for why it is an active
# binding and not a swap. It has to run from THIS package's `.onLoad`,
# the one moment this namespace is unsealed.
#
# What that asks of the code here, each one guarded by a test in
# tests/testthat/test-generic-collision.R:
#
#  * Every generic below stays a bare `UseMethod()`. While its owner is
#    loaded it is not the generic that runs, so work put in its body
#    would silently stop happening.
#  * Its formals are the owner's. Where two owners disagree, brms
#    decides: `rhat` takes posterior's `x`, which brms imports, and not
#    bayesplot's `object`; `posterior_samples` takes brms's `x, pars`
#    and not gratia's `model`.
#  * Every `frmtmb_draws` method is registered in this package's table
#    AND in each owner's, or it is unreachable while that owner is
#    loaded.

#' The generic names this package defines that another package owns,
#' with the owner to defer to.
#'
#' Audited with `parseNamespaceFile()` over every installed package and
#' then by loading each candidate and reading the environment of its
#' exported function, because a package that exports a name may only
#' have imported it: brms exports 20 of these 28 and defines 8
#' (dev/samplegen-audit.R). The owner is the package that DEFINES the
#' generic, because that is the table brms's own methods are in.
#'
#' Two names have rival definers whose generics are different closures.
#' `rhat` is defined by posterior and by bayesplot, and brms imports
#' posterior's; `posterior_samples` is defined by brms and by gratia.
#' The binding takes whichever the user would have reached without this
#' package: the search path first, then this order, which puts brms's
#' choice first. bbmle exports a `parnames` that is not a generic, so it
#' is not an owner and the binding would refuse it anyway.
#'
#' @noRd
sample_generic_owners <- list(
  as.mcmc = "coda",
  bayes_factor = "bridgesampling",
  bridge_sampler = "bridgesampling",
  post_prob = "bridgesampling",
  kfold = "loo",
  loo_moment_match = "loo",
  loo_subsample = "loo",
  psis = "loo",
  nsamples = "rstantools",
  posterior_epred = "rstantools",
  posterior_interval = "rstantools",
  posterior_linpred = "rstantools",
  posterior_predict = "rstantools",
  predictive_error = "rstantools",
  predictive_interval = "rstantools",
  log_posterior = "bayesplot",
  neff_ratio = "bayesplot",
  nuts_params = "bayesplot",
  rhat = c("posterior", "bayesplot"),
  mcmc_plot = "brms",
  parnames = "brms",
  posterior_samples = c("brms", "gratia"),
  pp_mixture = "brms",
  reloo = "brms",
  restructure = "brms",
  stancode = "brms",
  standata = "brms"
)
