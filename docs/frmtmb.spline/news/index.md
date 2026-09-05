# Changelog

## frmtmb.spline 0.1.0

First release. Two things a spline needs that a fit does not give you:
inference about the CURVE rather than about its coefficients, and a
survival family whose parameter is a spline.

### Curve inference

- [`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
  evaluates a fitted linear predictor on a grid and returns it with two
  intervals. The pointwise one is the usual delta-method interval. The
  simultaneous one covers the whole curve at once, by the max-deviation
  simulation of Ruppert, Wand and Carroll (2003, ch. 6): draw the
  curve’s own deviation process from its joint covariance, standardize
  by the pointwise standard error, take the largest absolute value over
  the grid, and use that distribution’s quantile in place of 1.96. This
  is the construction `gratia::confint(type = "simultaneous")` uses, and
  this package’s simulation reproduces gratia’s critical value on the
  same mgcv fit to inside two Monte Carlo standard errors at 10000,
  100000 and 500000 draws.
- The critical value is reported WITH its Monte Carlo standard error,
  which is the quantile standard error `sqrt(p(1-p)/n) / f(q)`. A
  critical value without one invites a comparison that simulation noise
  alone would fail: two packages agreeing to 0.007 at nsim = 10000 have
  agreed to half a standard error, and that is the honest way to say it.
- The comparison itself needed care. gratia standardizes a SMOOTH-ONLY
  deviation by `smooth_estimates()`’s `.se`, which is the FULL linear
  predictor’s standard error, intercept column included; the two differ
  by up to 13.5 percent on the test model and the critical value that
  comes back is 8 percent smaller than the self-standardized one. Both
  bands are exact, because the same divisor calibrates and scales, but
  the critical VALUES are not comparable unless the divisor is. The
  package’s own bands standardize by their own standard error.
- [`frm_curve_deriv()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_deriv.md)
  gives the first or second derivative of the same curve with
  delta-method standard errors and the same two intervals. The
  derivative is taken of the DESIGN, not of the fitted values, so the
  estimate and its standard error describe one function.
- The step size is not gratia’s. Swept against the exact derivative of a
  known function, a central difference bottoms out near `eps = 1e-6` at
  4.0e-10 for the first derivative, and near `eps = 1e-4` at 2.2e-07 for
  the second, because a second difference divides by `eps^2` and turns
  cancellation into the dominant error. gratia’s fixed `eps = 1e-7` is
  14 times off the optimum at order 1, which is negligible, and eight
  orders of magnitude off at order 2, which is not: on the test fit
  gratia’s second derivative differs from mgcv’s own same-stencil answer
  by 1.76, 52 percent at the grid edge, and this package’s agrees with
  it to 8.6e-05. So `eps` is scaled by the covariate’s range and set per
  order.
- [`frm_curve_feature()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_feature.md)
  locates a peak, a trough or a level crossing and gives its POSITION a
  standard error, by the implicit-function delta method:
  `var(t*) = var(f'(t*)) / f''(t*)^2` for a stationary point and
  `var(f(t*)) / f'(t*)^2` for a crossing. These are the numbers a
  movement paper reports, the time of peak velocity and movement onset,
  and they are usually reported without an interval because nothing
  hands one over.
- At a stationary point the curve’s own value has a simplification worth
  knowing and the function reports it: the derivative of `f(t*)` with
  respect to the coefficients loses its `f'(t*) dt*/dc` term because
  `f'(t*)` is zero there, so the standard error of the PEAK HEIGHT is
  the ordinary pointwise standard error at `t*`, not inflated by the
  uncertainty in where the peak is.
- Every root the grid brackets is returned, so a curve with two peaks
  gives two rows and a curve with none gives a zero-row answer rather
  than an error. “This curve does not peak in this window” is an answer.

### The covariance these three needed, and the core seam

- All three need the same object, the covariance of the whole grid
  prediction, and frmtmb exports no route to it. A penalized smooth’s
  wiggly part is a random-effect block even when the smooth is a
  population term, so the curve covariance needs the joint covariance of
  the fixed AND random coefficients. `vcov(full = TRUE)` excludes `b`
  under both of its branches, by construction: its documented invariant
  is that its row names are
  [`confint()`](https://rdrr.io/r/stats/confint.html)’s.
  `predict(se.fit = TRUE)` forms the grid covariance internally and
  returns its diagonal.
- So the package rebuilds it, out of two facts. The linear predictor is
  LINEAR in the coefficients, so the difference between a prediction and
  the same prediction with one coefficient raised by one is that
  coefficient’s design column, exactly and with nothing to tune. And the
  joint covariance is the inverse of the fit’s own joint precision.
- Neither piece came from an exported function, so neither is trusted.
  Every call recomputes `sqrt(diag(Sigma))` and compares it with
  `predict(se.fit = TRUE)`, and refuses when they disagree by more than
  `tol`. Measured agreement on the package’s own models is 3e-15 to
  9e-15 relative, which is machine precision, and it is reported by
  [`print()`](https://rdrr.io/r/base/print.html) rather than hidden.
- **The covariance is core’s own, read from its cache rather than
  recomputed.** `predict(se.fit = TRUE)`, which this package calls
  anyway for the check, memoizes the inverted joint precision on the
  fit; the curve then subsets it. Recomputing it instead, which the
  first draft did, paid for a second `sdreport()` per call, took a Schur
  complement over the coefficients the curve does not touch on a
  densified copy of a sparse matrix, and bypassed
  `autoscale_sdreport()`. Measured at 8006 random coefficients: 114 s
  and 2.1 GB before, 8.2 s and 1.2 GB after, with the covariance subset
  itself at 0.00 s. An autoscaled fit now works rather than being
  refused, at 3.7e-16 against `predict(se.fit = TRUE)`.
- **The cost is the joint-precision solve, not the call count.**
  Measured at `re.form = NA` on a 20-point grid: design rebuild 0.01 s
  against `predict(se.fit = TRUE)` 0.29 s at 8 random coefficients, 0.07
  s against 0.98 s at 2006, and 0.28 s against 6.87 s at 8006. The term
  the call count measures is a tenth of the cost at every size. Size a
  job from the number of coefficients in the FIT, not from the grid.
- The design rebuild is one
  [`predict()`](https://rdrr.io/r/stats/predict.html) call per
  contributing coefficient plus one probe per block of 24 that
  contributes nothing. Under the default `re.form = NA` a per-subject
  grouping block contributes nothing and is skipped in blocks rather
  than one coefficient at a time: measured, the vignette’s 20-subject
  factor-smooth model costs 32 calls against 110 random coefficients,
  and the count does not grow with the number of LEVELS. It is not the
  ideal 13 either, for two stated reasons: a chunk that straddles the
  boundary between a live block and a dead one is expanded whole, and a
  component short enough to skip the probe is expanded whether or not it
  contributes.
- `dev/spline-seam-proposal.md` says what core would have to export to
  make the rebuild unnecessary, and what a penalized coefficient block
  handed to a NONLINEAR body would need beyond that.

### The Royston-Parmar family

- [`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
  writes the log cumulative hazard as a natural cubic spline in log
  time, which is the flexible parametric survival model of Royston and
  Parmar (2002). `scale = "hazard"` is the proportional hazards version,
  `"odds"` the proportional odds version and `"normal"` the probit one;
  all three are here, because the second and third are one `switch` arm
  each once the first is written.
- It is parameterized exactly as `flexsurv::flexsurvspline()`
  parameterizes it, and that is testable rather than claimed. Taking
  flexsurv’s fitted coefficients on its own `bc` data (686 rows, 299
  events) and evaluating frmtmb’s objective at them reproduces
  flexsurv’s log likelihood to between 1.4e-16 and 1.6e-15 relative
  across all three scales and 0, 1 and 3 interior knots. That is an
  identity, not an agreement. The two optima then match as well, and
  frmtmb’s is never the worse of the two.
- `mu` is `gamma0` and the rest are `gamma1`, `gamma2`, and so on, all
  with identity links. A formula on `mu` is proportional hazards; a
  formula on any other coefficient is a TIME-VARYING effect, which is
  what flexsurv spells `anc =`. Each is an ordinary distributional
  parameter, so `s()`, a random effect and a prior all reach them.
- Knots go at equally spaced quantiles of the log UNCENSORED times, as
  Royston and Parmar place them and as flexsurv defaults. The quantiles
  need the response, so they are taken at frame assembly through
  `family_finalize()` and the family the fit carries has them baked in;
  `knots =` and `bknots =` pin them at construction instead. `df` counts
  interior knots plus one, so flexsurv’s `k` is `df - 1`.
- The family declares an `lcdf` as well as an `lpdf`, so `cens()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) both reach the
  likelihood, and right, left and interval censoring all work. Both are
  registered CONDITIONAL rather than working, for the reason below. Note
  the coding: frmtmb’s `cens()` reads 0 as an observed event and 1 as
  right censored, which is the OPPOSITE of a `Surv()` status column.
- **[`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md)
  is the check this family cannot do without, and it refuses by
  default.** Two things in
  [`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
  are floors rather than answers and both are silent in the fitted
  object. This function recomputes them at the fitted parameters and
  stops, naming the row count, the maximum, the threshold, the reason
  and the remedy.
  [`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
  and its two companions call it, so the documented way to inspect this
  family will not draw a curve off a fit whose likelihood is a floor
  artifact.
- The first floor is the one that made the family unshippable without
  this check. frmtmb forms a right-censored contribution as
  `log(1 - F(y))` on the PROBABILITY scale, and core offers a family no
  complementary log-CDF slot to hand back `log S` directly, so the
  scored `log S` carries absolute error about `.Machine$double.eps / S`
  whatever the family does internally: exact to 1.3e-13 at `-log S` of
  10, wrong by 1.7e-04 at 30, floored at -35.127363 past 36. Past 30 the
  term is FLAT, its gradient exactly zero, so the optimizer prices such
  a row at a constant and fits the rest as if it were free. Measured: a
  600-subject fit with one subject censored far beyond every event time
  converges without a warning and puts its treatment coefficient tens of
  percent out, on data flexsurv declines to fit at all. How far the
  reported log likelihood is from the model’s depends on the data, not
  on the family: a floored row contributes -35.127363 instead of its own
  `-log S`, so the shortfall is about `-log S - 35` per row, and two
  runs of that design differing only in seed give 2.4e+03 and 2.166e+04.
  The threshold is one quantity for all three scales, because the error
  depends only on `S`: `-log S` is `H` on `"hazard"`,
  `log(1 + exp(eta))` on `"odds"` and `-log(Phi(-eta))` on `"normal"`,
  and 19.2 is where `eps / S` passes 1e-8.
- The refusal is POST-FIT and cannot be otherwise.
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) reads
  `object$opt$objective` directly and the family protocol has no hook
  that runs when a fit finishes, so nothing in this package can make
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) or
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) refuse on their own, and
  the optimizer may have walked through or stopped inside the flat
  region before the check is called. The real fix is a core one, an
  `lccdf` slot, and it is five sites rather than one: see
  `dev/spline-seam-proposal.md`.
- Monotonicity of the cumulative hazard is NOT enforced, and flexsurv
  does not enforce it either. What this family does instead is refuse to
  answer `NaN` for it: where the spline’s derivative in log time goes
  non-positive, the log density is a large finite negative number, since
  `NaN` stops the optimizer and one `NaN` component poisons a mixture’s
  log-sum-exp. That floor is NOT inert when it is used: it turns
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) and
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) into a pseudo-likelihood,
  measured at 6 floored rows and 3952 units on a cure-fraction dataset.
  [`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md)
  counts those rows too and refuses on them. The floor is a smooth
  positive part rather than a branch, because RTMB refuses a comparison
  on an AD type outright and has no `CondExp`. It changes the log
  density by 2.4e-15 at a derivative of 1.
- That floor had to be written twice, and the measurement is what caught
  it. The textbook spelling, `0.5 * (u + sqrt(u^2 + eps2))`, is right on
  paper and wrong in double precision on the side it exists for: at
  `u = -35.75` the square root rounds to `|u|`, the sum cancels to
  exactly zero and [`log()`](https://rdrr.io/r/base/Log.html) returns
  `-Inf` after all. Measured, 647 of 686 rows reached `-Inf` that way.
  Recovering the small branch through `(u + s)(s - u) = eps2` and ADDING
  it back does not help either, because adding 1.4e-16 to 71.5 loses it.
  The form that works SELECTS between the two cancellation-free branches
  with [`sign()`](https://rdrr.io/r/base/sign.html), which is available
  on an advector and whose zero derivative is the correct one because
  the branches agree at `u = 0`.
- The spline basis is written branch-free, with the truncated power
  spelled `0.5 * (e + abs(e))`, so that it tapes if the response is ever
  promoted to a parameter. It is a function of the response, which is
  data, so nothing normally needs that; it costs one line to have it
  anyway.

### Deliberate omissions

- No `post$mean_fn`, so
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  `predict(type = "response")` are refused. The mean of a Royston-Parmar
  survival time has no closed form, and core’s `cox()` refuses for the
  same reason. Read the fitted log cumulative hazard with
  [`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
  instead.
- No exact basis derivative in
  [`frm_curve_deriv()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_deriv.md).
  The design is rebuilt through
  [`predict()`](https://rdrr.io/r/stats/predict.html), which evaluates a
  basis and never differentiates one, so an exact derivative would mean
  reading mgcv smooth objects out of the fitted frame. That is a deeper
  reach into core than anything else here makes, for an error already at
  the tenth significant figure.
- Derivatives of order 3 and higher are refused rather than offered: a
  third central difference divides by `eps^3` and there is no step size
  at which it is accurate.
- No accurate log survival past `-log S` of about 19.2, and none at all
  past 36. See the
  [`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md)
  entries above: the squeeze that keeps the term finite rather than
  `-Inf` cannot also keep it accurate, because the number core asks a
  family for is `F` and the complement of a double near 1 is not
  representable. flexsurv computes `log S` directly and has no such
  region. The package’s answer is refusal by name rather than a silent
  floor, and the real fix is a core `lccdf` slot.
