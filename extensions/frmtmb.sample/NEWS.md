# frmtmb.sample (development version)

* The defaults this package registers with `frmtmb::get_prior()` now
  answer `route = "sample"` only. Loading this package used to change
  the table `get_prior()` returned for every caller, an author asking
  about `frm()` included; it no longer does. Ask for
  `get_prior(..., route = "sample")` to see what `frm_sample()`
  applies, and with this package unloaded that route is refused rather
  than reported as flat. `frm_sample()` itself is unchanged: the same
  defaults under the same call-over-fit-over-defaults precedence.

* The frmtmb floor rises to 0.51.0. The registration in `.onLoad()`
  now passes `frmtmb_register_compat(expects = )`, which frmtmb gained
  in 0.51.0, and an older frmtmb stops the load outright with
  `unused argument (expects = c("hmm", "lca"))`. The old floor of
  0.46.0 permitted an installation that could not load, so this is a
  hard requirement rather than a preference.

* `frmtmb.latent` joins `Suggests`. `hmm()` and `lca()` left frmtmb
  for that package, so the tests here that sample one of them qualify
  the call and skip when it is absent. The two compatibility rules this
  package registers for those pairs are declared in
  `frmtmb_register_compat(expects = )`, which says they name a feature
  another package owns. Before frmtmb checked, those rules were dropped
  without a word in every session that did not also load
  `frmtmb.latent`; now `frm_compat()` reports them as unresolved until
  it is loaded.

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
