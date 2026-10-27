# frmtmb.spatial 0.3.0

Needs frmtmb 0.76.0, which takes `frmtmb_register_autocor(jacobian =)`
and accepts a bound at an end of a residual correlation parameter's
interval.

## New features

* `sar()` gives its whitening map's Jacobian, `(I - rho W) / sigma`, as
  a sparse matrix. So frmtmb.sample's `log_lik()` returns brms 2.23.0's
  leave-one-out conditional density for both types, gaussian and
  student. It equals brms's `log_lik()` at the same parameters to
  6.3e-16 relative. One draw on a 40 x 40 lattice (1600 rows) takes 1.9
  ms; brms's dense algorithm takes 4.1 s for the same draw.

## Bug fixes

* The `default_prior()` table of a `sar()` model passes back as
  `prior =`. Its `lagsar` or `errorsar` row lists the ends of the
  interval as `lb` and `ub`, as brms's does. `frm()` refused them as
  outside the interval. The fix is in frmtmb, which now accepts a bound
  at an end of the interval as no bound. Under `frm_sample()`, a
  bounds-only prior on `lagsar` keeps the flat default density; see
  frmtmb.sample's NEWS.

# frmtmb.spatial 0.2.0

Needs frmtmb 0.75.0, which exports `frmtmb_register_autocor()`, the
seam through which this package registers `sar()`.

## New features

* `sar(M, type = )`, brms's spatial simultaneous autoregressive term, in
  its lag and error forms. As in brms 2.23.0, it is for `gaussian()` and
  `student()` with a constant `sigma` and `nu`. It is registered through
  `frmtmb::frmtmb_register_autocor()`. brms leaves out `-N/2 log(2 pi)`
  (`-N/2 log(pi)` for student). Less that constant, the log-likelihood
  equals brms's Stan program to 1.4e-14 (gaussian, at frmtmb's
  estimates) and 4e-15 (student, at set parameter values with `nu = 5`).
  brms fits of 64 rows agree with the maximum likelihood fits to within
  0.44 standard errors in every parameter. `fitted()` reports brms's
  `posterior_epred()`, `(I - rho W)^-1 mu` for the lag form. `summary()`
  and `confint_varcorr()` report `lagsar` or `errorsar`, and
  `set_prior()` takes them as classes: a density on `rho`, carried to
  the internal parameter.

* `fitted()`, `predict()` and `simulate()` take `newdata` with one row
  per row of `M`, read as the locations in their order, as brms reads
  it. When the fit dropped rows for missing values, `newdata` covers
  them too, and `I - rho M` is then over every row of `M`. A
  counterfactual such as `x + 1` everywhere then gives the lag model's
  spillover, `(I - rho W)^-1 mu_new`. That equals brms's
  `posterior_epred()` on the same `newdata` at the same parameter
  values. Other `newdata` is refused, as brms refuses it without
  `newdata2`.

* `rho` is estimated on `(1 / min eigenvalue, 1 / max eigenvalue)` of
  `M`, where `I - rho M` has a positive determinant. brms bounds it by
  `(min eigenvalue, max eigenvalue)`. For a row-standardized `M` that
  cuts off negative values below the smallest eigenvalue. So frmtmb can
  estimate `rho` beyond brms's bound, and the two are then different
  models. An example is a row-standardized 8 x 8 queen lattice
  (eigenvalues -0.5006 to 1) with data drawn at `rho = -0.8`. frmtmb's
  estimate is -1.103, with a logLik 3.74 above the best inside brms's
  interval. brms's posterior median is -0.443, with a fifth of its draws
  within 0.02 of its bound. For a matrix that is not row-standardized,
  brms's interval reaches past `1 / max eigenvalue`, where its
  log-determinant is not a number and Stan rejects the proposal. So
  brms's support there is frmtmb's interval, at the cost of rejected
  proposals.

* A matrix whose smallest or largest eigenvalue repeats for all rows but
  one, such as a complete graph, gives a warning. At that end of the
  interval `I - rho M` keeps one direction only. With an intercept the
  likelihood has no maximum there, so `rho` runs to the end. The model
  is not refused, because without an intercept it has its maximum.

* A row-standardized `M` of an irregular graph, spdep's default, has
  real eigenvalues that `eigen()` returns with a rounding imaginary
  part. brms 2.23.0 passes them to Stan as they come and its fit fails
  ("variable does not exist; ... eigenMsar"); frmtmb takes the real part
  when the imaginary part is below 1e-8 of the largest modulus, and
  refuses a matrix with complex eigenvalues.

* spdep's `nb` and `listw` objects are read as brms reads them
  (`nb2mat()` row-standardized, `listw2mat()`) without needing spdep.

# frmtmb.spatial 0.1.0

The first release. Load it with `library(frmtmb.spatial)`, which
attaches frmtmb as well. Needs frmtmb 0.74.0, which exports
`frmtmb_register_covstruct()`.

## New features

* New package. Two fields moved here from frmtmb.
  `car(M, gr = g, type = )` is brms's conditional autoregressive field
  in its four types (`"escar"`, `"esicar"`, `"icar"`, `"bym2"`).
  `spde(fem, gr = node)` is a Matern field over a finite-element mesh.
  Both are registered through `frmtmb::frmtmb_register_covstruct()`, and
  the package reads no frmtmb internal. Fits are unchanged. On the moved
  tests' models, these are identical to frmtmb 0.73.0's: the
  log-likelihood, every estimate and standard error, `fitted()` with and
  without an unseen location, `ranef()` conditional standard deviations,
  `confint_varcorr()` and `simulate()`. `frmtmb.sample::frm_sample()`
  gives the same draws from the same seed. `frm_compat()` gives `car`
  and `spde` the statuses frmtmb gave them.
