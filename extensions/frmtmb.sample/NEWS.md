# frmtmb.sample 0.1.0

First release, extracted from frmtmb 0.46.0.

* NUTS sampling for frmtmb models through tmbstan: `frm_sample()` on
  a fit or a formula, with brms 2.23's default priors on both routes,
  priors stacking most-explicit-first, and non-centered sampling of
  qualifying random-effect blocks. Draws carry frmtmb parameter
  names.
* The full posterior method surface on core's generics and on the
  brms, loo, posterior, bayesplot, rstantools, coda and
  bridgesampling generics: `log_lik()`, `loo()`, `waic()`,
  `posterior_epred()`, `posterior_predict()`, `hypothesis()`,
  `conditional_effects()` on draws, `bayes_R2()`, `mcmc_plot()`, and
  the rest.
* `check_laplace()` measures the Laplace and Wald approximations
  against NUTS on the same objective, centered and with no added
  priors.
* At load, the package registers its sampling default priors with
  core's `get_prior()` and its feature rows with the compatibility
  matrix, so both stay truthful whether or not this package is
  installed.
* Refuses a tmbstan built against StanHeaders 2.39 or later, whose
  code generator leaves one log-density overload unpatched so every
  chain silently samples a standard normal instead of the model.
* Parallel chains work on every platform, Windows included: tmbstan
  rebuilds the tape on each socket worker from the serialized
  objective closure, and frmtmb's generated closures are
  self-contained.
