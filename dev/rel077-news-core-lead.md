The main changes:

- BREAKING: under the covariance form of `ar()`, `ma()` and `arma()`,
  `sigma` is the innovation sd, as in brms 2.23.0. A fit saved by
  0.76.0 or earlier with such a term is refused where an answer reads
  `sigma`.
- BREAKING: each random-effect parameter and the residual correlation
  of a multivariate model have one prior slot, as in brms. The
  named-list prior spelling is refused. `set_prior("")` is brms's
  `(flat)` on its slot.
- `mo()` takes a variable with one or two values, as brms does.
- marginaleffects reads `vcov =` and brms's prediction arguments, and
  gives finite standard errors on a fit with an empty `mo()` category.
- `build_objective()` reads a `prior_only` field, for
  `sample_prior = "only"` in frmtmb.sample. The sampling-API function
  `autocor_innov_sd()` is removed.
- Update the extensions together with frmtmb. frmtmb.sample 0.24.0
  needs frmtmb 0.77.0.
