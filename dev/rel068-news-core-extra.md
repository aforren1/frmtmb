## Breaking changes

* **`conditional_effects(method = "predict")` reports the median of
  its simulated responses as `estimate__`**, as brms reports the median
  of its predictive draws. It reported the expected response, which
  `method = "epred"` still shows. The band is unchanged, and the
  estimate is the 50% point of the same draws: a whole count for a
  discrete family when `ndraws` is odd, inside the bounds of a
  truncated one. At a Poisson fit's x = -1, 0, 1 it is 1, 2, 4, where
  it was 1.46, 2.52, 4.35, and brms's median at the same parameters is
  1, 2, 4 (`dev/fixes-ce-pred.R`, `test-ce-predict-median.R`).
  frmtmb.sample's draws route reported the median already. The user's
  decision of 2026-10-06.

## New features

* **`cs()` on `cumulative()`**, as brms 2.23.0 fits it, with brms's
  warning that category-specific effects for the family are
  experimental. It was refused. Row `i` reads the thresholds
  `tau_k - cs_ik`, brms's `Intercept - transpose(mucs[n])`, and the log
  density equals brms's compiled program to at most 0.9 ulp at the
  optimum and at three perturbed points on six shapes: the logit and
  the probit, a predictor beside `cs()`, two `cs()` terms, a factor in
  `cs()`, `disc ~ 0 + z`, equidistant thresholds and a mixture with
  `sratio()` (`dev/rel068-cs-lpcheck.R`). The offsets can make a row's
  thresholds cross: such a row's density and `fitted()` probabilities
  are `NaN`, as brms's density is (brms's `posterior_epred()` returns
  the negative difference), and `simulate()` gives `NA` there.
  `predict()` gives `NA` proportions for a row that crosses at the
  estimates; a row that crosses only in some of its simulated
  parameter draws reports the proportions over the other draws, and
  the call warns with the draws dropped per row, here and on
  `hurdle_cumulative()` (the user's decision of 2026-10-06).
  `hurdle_cumulative()` and a `cumulative()` component of an ordinal
  mixture take `cs()` the same way and warn as brms does, once per
  component. A family that takes no `cs()` is refused in brms's words,
  "Category specific effects are not supported for this family". The
  user's decision of 2026-10-06.

* **Ordinal mixtures and the per-threshold features meet.**
  `default_prior()` lists the per-threshold class `"Intercept"` rows of
  each component under `dpar = "mu<k>"`, and of the vector an
  `order = "mu"` mixture shares once with no `dpar`, as brms lists
  them; `set_prior(class = "Intercept", coef = "2", dpar = "mu1")` puts
  a density on that threshold alone. `confint()` and
  `vcov(full = TRUE)` name a mixture's thresholds by what they are,
  `mu1_Intercept[1]` and `log(mu1_Intercept[2] - mu1_Intercept[1])`,
  and the shared vector once, as `Intercept[1]`. The compatibility
  table's cells for `disc`, `"equidistant"` and `"sum_to_zero"` with
  `mixture` say what the release does: equidistant is refused under
  `order = "mu"`, as in brms.

## Performance

* **The kriging draw on a high-rank grid is faster.** Past a quarter
  of the unseen positions, `gp_krig_factor()` takes the whole
  conditional covariance and LAPACK's pivoted Cholesky instead of its
  column loop: 0.020 s against 0.071 s at rank 400 of 400 positions and
  0.161 s against 0.293 s at rank 499 of 1000, with the law unchanged
  (`dev/rel068-krig-rank.R`). The low-rank regime the loop was built
  for is untouched.

## Bug fixes

* `?refit` says that `refit()` does not repeat the warning `frm()`
  gives when a nonlinear model's likelihood is flat at the optimum: the
  design and its flat direction are the original fit's.

* `vignette("brms-migration")` was stale on `hurdle_cumulative()` (it
  takes the threshold structures, `thres(gr = )` and `cs()`), on
  `acat()`'s links, on `disc`, and on `order = "mu"` for ordinal
  mixtures, and now says what `cs()` on the cumulative families does.
