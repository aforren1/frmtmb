# Changelog

## frmtmb.sample (development version)

- `frmtmb.latent` joins `Suggests`. `hmm()` and `lca()` left frmtmb for
  that package, so the tests here that sample one of them qualify the
  call and skip when it is absent. The compatibility rules this package
  registers for those two pairs are unchanged and need no guard: a rule
  naming a feature the matrix does not have contributes nothing, so
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html)
  is correct whether or not the other package is loaded.

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
