# Changelog

## frmtmb.sample (development version)

- [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  refuses parallel chains on Windows when the model’s family comes from
  a namespace built by
  [`pkgload::load_all()`](https://pkgload.r-lib.org/reference/load_all.html),
  naming the package and the two remedies. A worker process cannot load
  such a namespace, and R silently substitutes the global environment
  for it, so the family failed at its first internal call and rstan
  reported only that the chains returned no draws.

- The defaults this package registers with
  [`frmtmb::get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.html)
  now answer `route = "sample"` only. Loading this package used to
  change the table
  [`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.html)
  returned for every caller, an author asking about
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
  included; it no longer does. Ask for
  `get_prior(..., route = "sample")` to see what
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  applies, and with this package unloaded that route is refused rather
  than reported as flat.
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  itself is unchanged: the same defaults under the same
  call-over-fit-over-defaults precedence.

- The frmtmb floor rises to 0.51.0. The registration in `.onLoad()` now
  passes `frmtmb_register_compat(expects = )`, which frmtmb gained in
  0.51.0, and an older frmtmb stops the load outright with
  `unused argument (expects = c("hmm", "lca"))`. The old floor of 0.46.0
  permitted an installation that could not load, so this is a hard
  requirement rather than a preference.

- `frmtmb.latent` joins `Suggests`. `hmm()` and `lca()` left frmtmb for
  that package, so the tests here that sample one of them qualify the
  call and skip when it is absent. The two compatibility rules this
  package registers for those pairs are declared in
  `frmtmb_register_compat(expects = )`, which says they name a feature
  another package owns. Before frmtmb checked, those rules were dropped
  without a word in every session that did not also load
  `frmtmb.latent`; now
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html)
  reports them as unresolved until it is loaded.

## frmtmb.sample 0.1.0

First release, extracted from frmtmb 0.46.0.

- NUTS sampling for frmtmb models through tmbstan:
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  on a fit or a formula, with brms 2.23’s default priors on both routes,
  priors stacking most-explicit-first, and non-centered sampling of
  qualifying random-effect blocks. Draws carry frmtmb parameter names.
- The full posterior method surface on core’s generics and on the brms,
  loo, posterior, bayesplot, rstantools, coda and bridgesampling
  generics:
  [`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md),
  [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html),
  [`waic()`](https://aforren1.github.io/frmtmb/reference/loo.html),
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html),
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html)
  on draws,
  [`bayes_R2()`](https://aforren1.github.io/frmtmb/reference/bayes_R2.html),
  [`mcmc_plot()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
  and the rest.
- [`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md)
  measures the Laplace and Wald approximations against NUTS on the same
  objective, centered and with no added priors.
- At load, the package registers its sampling default priors with core’s
  [`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.html)
  and its feature rows with the compatibility matrix, so both stay
  truthful whether or not this package is installed.
- Refuses a tmbstan built against StanHeaders 2.39 or later, whose code
  generator leaves one log-density overload unpatched so every chain
  silently samples a standard normal instead of the model.
- Parallel chains work on every platform, Windows included: tmbstan
  rebuilds the tape on each socket worker from the serialized objective
  closure, and frmtmb’s generated closures are self-contained.
