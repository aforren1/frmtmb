# frmtmb.sample (development version)

Needs the frmtmb of the same release, whose `sigma` under the covariance
form of `ar()`, `ma()` and `arma()` is the innovation sd.

## Breaking changes

* BREAKING: under the covariance form of `ar()`, `ma()` and `arma()`,
  the draws of `sigma` are of the innovation sd, as brms's are. Before,
  they were of the marginal sd. On an AR(1) with `ar` near 0.92, the
  posterior mean of `sigma` was 1.154 against brms's 0.432. It is now
  0.433. `summary()`, `posterior_summary()`, `variables()` and the draws
  read the same parameter as brms's.
* The default `student_t(3, 0, s)` on `sigma` and a prior of your own
  are densities on the innovation sd itself. A `lb` or `ub` on `sigma`
  is taken, as in brms. Under `sigma ~ x`, the default and a prior of
  your own on the intercept of `sigma` are read at the log innovation
  sd, as brms's `Intercept_sigma`. With `normal(-1, 0.5)` there, the
  posterior mean of `ar` was 7.95 Monte Carlo standard errors from
  brms's and is now 0.55.
* We compared the draws with brms's own (4 chains of 4000 on each
  side): `ar()`, `ma()` and `arma()` under the defaults, `ar()` under
  user priors, and `ar()` under `sigma ~ x`. Every posterior mean is
  within 1.8 Monte Carlo standard errors of brms's. Every sd ratio is
  within 0.969 to 1.037. The bulk effective sample size of `sigma` is
  2.0, 0.94 and 1.2 times what it was under `ar()`, `ma()` and
  `arma()`.
* The message of the default priors no longer marks `sigma` with "on
  brms's innovation sd": it is brms's sigma.

