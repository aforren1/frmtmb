### brms-24. `com_poisson`: the R-side pmf is not the pmf the Stan program fits

- Seen: 2.23.0. Confirmed 2026-10-13 (source read and a fit; lane
  remeasure, R6).
- Repro: brms's Stan function `com_poisson_log_lpmf(y | log_mu, nu)`
  (`chunks/fun_com_poisson.stan`) is `y * log_mu - nu * lgamma(y + 1)
  - log_Z_com_poisson(log_mu, nu)`, with series terms
  `(k - 1) * log_mu - nu * lgamma(k)`: the pmf is proportional to
  `mu^y / (y!)^nu`, so the fitted `mu` is Shmueli's rate lambda.
  brms's R side, `brms:::dcom_poisson(x, mu, shape)` (used by
  `log_lik()`, `posterior_predict()` and `loo()`), is
  `shape * (x * log(mu) - lgamma(x + 1)) - log_Z`, with
  `brms:::log_Z_com_poisson()` summing `shape * (k * log_mu - lfac)`:
  proportional to `(mu^x / x!)^shape`, so there `mu = lambda^(1/shape)`.
  At `mu = 3, shape = 2` the two pmfs have different means: 1.453549
  on the Stan side, 2.737078 on the R side (exact sums over 0:200,
  `dev/remeasure-compois-means.R`).
- Expected: one pmf. Actual: post-processing of a `com_poisson` fit
  evaluates and simulates a different pmf from the one sampled, for
  every `shape != 1`; `shape = 1` agrees (Poisson). The R-side
  normalizer is also approximate in its asymptotic branch (mass
  1.00014 at `mu = 3, shape = 2`).
- Workaround: none needed in frmtmb, which calls neither. Its own `mu`
  is the mean (`RTMB::dcompois2()`), a third parameterization, and
  `?brmsfamily` and `?frmtmb-families` say so.
- Record: `dev/remeasure-findings.md` R6; `dev/reviews/2026-10-13-remeasure.md`
  section C, `dev/remeasure-rev-compois.R` (brms fit: 2 chains of 1000,
  seed 20261013).
- Status: not reported.

### brms-25. The simplex response check is a mean over rows

- Seen: 2.23.0, 2026-10-13 (lane famsimplex, brms-A).
- Repro: a 300-row `dirichlet()` data set with row 1 off the simplex
  by 1e-7 (`dev/famsimplex-brms-tol.R`, `dev/famsimplex-brms-tol2.R`).
  `brms:::data_response()` tests `is_equal(rowSums(Y), 1)`, which is
  `all.equal()`: the mean absolute departure, here 3.3e-10, against
  1.5e-8.
- Expected: refuse the row, "Response values in simplex models must
  sum to 1." Actual: `standata()` accepts it. `dirichlet()` then fails
  inside Stan with an unrelated message ("Exception: no more scalars
  to read") and `brm()` returns a fit with no draws;
  `logistic_normal()` samples, scoring a row that is not a
  composition.
- Workaround: frmtmb tests each row at brms's 1.5e-8 and refuses such a
  row by number. In one case this refuses data brms fits: with one row
  of 300 off by 1e-7, brms's `logistic_normal()` samples.
- Record: `dev/famsimplex-findings.md`, section 6, brms-A.
- Status: not reported.

### brms-26. A categorical category's own formula stops `stancode()`

- Seen: 2.23.0, 2026-10-13 (lane famsimplex, brms-B;
  `dev/famsimplex-brms-first.R`).
- Repro: `stancode(bf(f ~ x, mub ~ z), data = d, family =
  categorical())`, `f` a three-level factor.
- Expected: Stan code with `X_mub` on `z`, as brms writes for the same
  formula on `dirichlet()`. Actual: "Error in use_glm_primitive(tmp) &&
  all.equal(tmp$dpars$mu$fe, bterms$dpars[[1]]$fe) : invalid 'y' type
  in 'x && y'" from `use_glm_primitive_categorical()`: `all.equal()`
  returns a string when the two designs differ, and `&&` refuses it.
- Workaround: frmtmb fits it.
- Record: `dev/famsimplex-findings.md`, section 6, brms-B.
- Status: not reported.

### brms-27. `posterior_epred()` for `wiener()` is right only at `bs = 1`

- Seen: 2.23.0, 2026-10-13 (lane famuni; `dev/famuni-wiener-findings.md`,
  section 5).
- Repro: `brms:::posterior_epred_wiener()` is
  `ndt - bias / mu + bs / mu * (exp(-2 mu bias) - 1) /
  (exp(-2 mu bs) - 1)`. On a 48-cell grid against the exact
  unconditional mean decision time (frmtmb.eam's conditional means
  weighted by the boundary probabilities) its relative error is at
  most 4.07e-16 at `bs = 1` and up to 0.637 at `bs != 1` (`bs = 2.5`,
  `mu = 1.5`, `bias = 0.7`). With the absolute starting point
  `z = bias * bs` in place of `bias` the formula agrees to 1.04e-15. A
  20000-draw RWiener Monte Carlo agrees with the exact value on four
  cells (0.44867 against 0.44902, 0.79461 against 0.79217). The
  formula is also `NaN` at `mu = 0`.
- Expected: the mean of the Wiener first-passage time. Actual: a mean
  for `bs = 1` with the starting point read as relative, so every fit
  with a boundary separation other than 1 gets a wrong expected
  response time.
- Workaround: frmtmb.eam's `fitted()` is conditional on the observed
  boundary and does not use this formula; `?wiener` documents the
  difference.
- Record: `dev/famuni-wiener-findings.md`, sections 5 and 9.
- Status: not reported.
