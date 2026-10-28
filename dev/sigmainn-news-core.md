# frmtmb (development version)

## Breaking changes

* BREAKING: under the covariance form of `ar()`, `ma()` and `arma()`,
  `sigma` is the innovation sd, as in brms 2.23.0. Before, it was the
  marginal sd of a residual. The model was the same, but the parameter
  was not brms's. On an AR(1) with `ar` 0.93, `sigma(fit)` was 1.078
  and is now 0.406, against a brms posterior mean of 0.432. For an
  AR(1), the marginal sd is `sigma / sqrt(1 - ar^2)`. At any order, it
  is `sigma` times the square root of the process variance at unit
  innovation variance. The logLik and the other estimates do not
  change.
* So `sigma()`, `VarCorr()`, `summary()`, `confint()`, `vcov()`, the
  standard errors and `fitted(dpar = "sigma")` report the innovation
  sd. Under `sigma ~ x` the coefficients of `sigma` are those of the log
  innovation sd, as brms's `Intercept_sigma` and `b_sigma` are. Pearson
  residuals, `simulate()` and the pointwise log-density read the
  marginal sd, as before.
* A prior on class `sigma` is a density on the innovation sd itself,
  and its `lb` and `ub` bound it, as in brms. Before, the density was
  read at the marginal sd times a factor of `ar` and `ma`. A `lb` above
  0 or a `ub` was refused. Under `sigma ~ x`, a prior on the intercept
  of `sigma` is now read at the log innovation sd, as brms reads it.
  Before, it was read at the log marginal sd.
* The standard error of `log(sigma)` is the one of the innovation sd:
  0.109 against 0.213 for the marginal sd on the AR(1) above. Its
  correlation with the `ar` parameter falls from 0.88 to -0.34.
* `ar(p = 2, cov = TRUE)` and the other higher orders, which brms does
  not fit, read `sigma` as the innovation sd too.
* A fit saved by frmtmb 0.76.0 or earlier with one of these terms is
  refused by the methods that report `sigma` or read the model:
  `summary()`, `sigma()`, `VarCorr()`, `fixef()`, `vcov()`, `confint()`,
  `predict()`, `residuals()` and `simulate()`. `update()` fits it again.
* The sampling-API function `autocor_innov_sd()` is removed: no prior
  needs it now.

