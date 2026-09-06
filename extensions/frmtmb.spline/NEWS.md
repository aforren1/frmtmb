# frmtmb.spline 0.3.0

A curve past the fitted knot span warns once under the function that
drew it, the feature search refuses, and `frm_curve_deriv()` no longer
dies on its default simultaneous band. Requires frmtmb 0.53.0 for the
classed span warning.

* The three curve functions say when the grid leaves a `ps()` knot
  span. They are the doors a user actually draws a curve through, and
  they were the ones that said nothing useful: `frm_curve()` and
  `frm_curve_deriv()` read the curve through `frm_lp_basis()`, which
  did not raise it at all, while `frm_curve_feature()` reached it
  through `predict()`, which did, ELEVEN times on one call. The first
  two now warn once, with the span in the message and under the name of
  the function that was called. The third refuses instead, below.

* The span surfacing needs the frmtmb this round ships, not the 0.52.0
  the `Depends:` floor names today. It is driven by the
  `frmtmb_ps_span_warning` condition class, which is new in core this
  round; against a released frmtmb 0.52.0 the class is never raised, so
  `frm_curve()` and `frm_curve_deriv()` fall silent again and
  `frm_curve_feature()` stops refusing, with nothing to say so. The
  floor rises with core's version at consolidation.

* FIX: `frm_curve_deriv()` past a `ps()` knot span works on its own
  DEFAULT arguments. It warned and then died in `quantile.default()`
  with "missing values and NaN's not allowed", naming neither the span
  nor the function. Past the outer knot the derivative design is exactly
  zero, so those rows carry a standard error of exactly zero, and the
  max-deviation simulation standardized their (exactly zero) deviation
  by it. Such a point is covered with probability one and cannot be the
  argmax, so it now leaves the maximization; the band is computed over
  the rows that carry uncertainty and the zero rows get a zero-width
  one. A grid on which EVERY row is past the outer knot refuses by name
  instead. Dropping nothing is bit-identical to the old arithmetic, so
  an ordinary grid is unaffected.

* FIX: a grid laid exactly on the knot span is inside it. The span was
  checked on the widened difference stencil, which reaches a millionth
  of the grid's range past both ends, so a grid ending on a knot drew a
  warning from `frm_curve_deriv()` and a REFUSAL from
  `frm_curve_feature()`, and narrowing the grid to the span, which is
  what the refusal tells the user to do, reproduced the refusal. Both
  now re-ask on the grid that was passed, and only when the widened call
  reported something, so an ordinary call pays nothing. The
  `frm_curve_deriv()` warning counts grid rows rather than stencil rows
  as a result.

* BEHAVIOR CHANGE: `frm_curve_feature()` REFUSES a search whose bracket
  leaves the span, where it used to warn after the fact. A band drawn
  past the span is visible on the page; a peak or a crossing located
  past it leaves the function as a number with a standard error beside
  it and nothing to say which curve it came off, because the decaying
  partial sum has peaks and crossings of its own. The bracket is
  checked at the grid scan, before any root is refined, and again on
  the difference stencil at the located roots. Narrow `newdata` to the
  span, or refit with a larger `pad =`.

* BEHAVIOR CHANGE: `rp_floored()` returns `n_censored_deep` where it
  returned `n_censored_floored`. Nothing about the count changed; the
  name did. Since 0.2.0 the censored term is scored exactly through
  `lccdf` and no row past the threshold is floored, so a field that
  says "floored" reports the opposite of what the slot bought. On the
  review's 600-subject design the count is 1, at `-log S` of 55.73, and
  that row's contribution is -55.7302136 exactly where the old
  probability-scale form gave `-Inf` and the floor gave -35.127363.

* The package documentation no longer describes seams core has since
  supplied. `?frmtmb.spline-package` said `fit$cache$Vjoint` was an
  internal this package reaches into, which has not been true since
  0.2.0, and named "an `lccdf` slot, or a post-fit family hook" as
  fixes core owed, both of which landed and both of which this package
  uses. What is still missing is named for what it is: a per-row
  log-likelihood slot, which `loo()` and `waic()` need, and a per-group
  one, which `frm(importance =)` corrects with. The
  `vignette("royston-parmar")` paragraphs that predated the same seams
  are corrected with them.

# frmtmb.spline 0.2.0

The three things this package had to work around are seams in frmtmb
0.52.0 now, so it stops working around them. Requires frmtmb
(>= 0.52.0).

* The one internal this package reached into is gone. `frm_curve()` and
  its two companions read `fit$cache$Vjoint` directly, which was
  written by an `@noRd` function and documented nowhere; they now call
  the exported `frmtmb::frm_lp_basis()`. `?frm_curve`'s section on it
  changes from "The one internal this reaches into" to "The route to
  the covariance".
* The design rebuild is gone with it. `frm_curve()` used to rebuild the
  grid design by unit perturbation, one `predict()` call per
  contributing coefficient plus one probe per block of 24 that
  contributed nothing; `frm_lp_basis()` returns the design core already
  had. The `predict()` call count no longer depends on the number of
  coefficients at all, and the linearity probe that guarded the
  perturbation is no longer needed.
* `frm_curve()` now works on a NONLINEAR (`nl = TRUE`) linear predictor,
  which it used to refuse. `A` is a Jacobian there rather than a design,
  and `predict(se.fit = TRUE)` is refused for a nonlinear predictor, so
  there is no second route to check against: `cov_rel_error` is `NA` and
  `print()` says the check did not run rather than reporting a passing
  one that never did.
* `frm_curve()` now works on an `rr()` block at `re.form = NULL`, which
  it used to refuse through the covariance check. A reduced-rank block's
  loadings live in `theta`; the perturbation could not see the
  derivative with respect to them and the assembled standard errors came
  out 27 percent away from `predict(se.fit = TRUE)`'s.
  `frm_lp_basis()` carries the loading columns, and the two now agree
  exactly.
* An exact `gp()` term works too. Its kriging variance is not
  coefficient uncertainty, and `frm_lp_basis()` returns it separately as
  `extra_var` rather than folding it into `A V A'`.
* `royston_parmar()` declares frmtmb's new `lccdf` slot, so a
  right-censored row is scored from `log S` directly, in closed form on
  all three scales: `-exp(eta)` on `"hazard"`, `-log1p(exp(eta))` on
  `"odds"` and `pnorm(eta, lower.tail = FALSE, log.p = TRUE)` on
  `"normal"`. The floor at -35.127363 is gone, and with it the region
  past `-log S = 30` where the term was flat and its gradient exactly
  zero. `cens()` moves from `conditional` to `works` in the
  compatibility table.
* `rp_floored()`'s censored count is a DIAGNOSTIC rather than a
  refusal. It still reports censored rows whose fitted `-log S` passes
  19.2, because such a row is one the data barely constrain, but it no
  longer stops: the arithmetic is exact there now. The MONOTONICITY
  floor is unchanged and still refuses, because a non-positive
  `d(eta)/d(log t)` means no hazard exists and the reported likelihood
  is a pseudo-likelihood rather than a density.
* `royston_parmar()` declares `post$fit_check`, frmtmb's new fit-end
  hook, so a fit with a non-monotone row warns as it is returned rather
  than only when someone calls `rp_floored()`.

# frmtmb.spline 0.1.0

First release. Two things a spline needs that a fit does not give you:
inference about the CURVE rather than about its coefficients, and a
survival family whose parameter is a spline.

## Curve inference

* `frm_curve()` evaluates a fitted linear predictor on a grid and
  returns it with two intervals. The pointwise one is the usual
  delta-method interval. The simultaneous one covers the whole curve at
  once, by the max-deviation simulation of Ruppert, Wand and Carroll
  (2003, ch. 6): draw the curve's own deviation process from its joint
  covariance, standardize by the pointwise standard error, take the
  largest absolute value over the grid, and use that distribution's
  quantile in place of 1.96. This is the construction
  `gratia::confint(type = "simultaneous")` uses, and this package's
  simulation reproduces gratia's critical value on the same mgcv fit to
  inside two Monte Carlo standard errors at 10000, 100000 and 500000
  draws.
* The critical value is reported WITH its Monte Carlo standard error,
  which is the quantile standard error `sqrt(p(1-p)/n) / f(q)`. A
  critical value without one invites a comparison that simulation noise
  alone would fail: two packages agreeing to 0.007 at nsim = 10000 have
  agreed to half a standard error, and that is the honest way to say it.
* The comparison itself needed care. gratia standardizes a SMOOTH-ONLY
  deviation by `smooth_estimates()`'s `.se`, which is the FULL linear
  predictor's standard error, intercept column included; the two differ
  by up to 13.5 percent on the test model and the critical value that
  comes back is 8 percent smaller than the self-standardized one. Both
  bands are exact, because the same divisor calibrates and scales, but
  the critical VALUES are not comparable unless the divisor is. The
  package's own bands standardize by their own standard error.
* `frm_curve_deriv()` gives the first or second derivative of the same
  curve with delta-method standard errors and the same two intervals.
  The derivative is taken of the DESIGN, not of the fitted values, so
  the estimate and its standard error describe one function.
* The step size is not gratia's. Swept against the exact derivative of
  a known function, a central difference bottoms out near `eps = 1e-6`
  at 4.0e-10 for the first derivative, and near `eps = 1e-4` at 2.2e-07
  for the second, because a second difference divides by `eps^2` and
  turns cancellation into the dominant error. gratia's fixed `eps =
  1e-7` is 14 times off the optimum at order 1, which is negligible,
  and eight orders of magnitude off at order 2, which is not: on the
  test fit gratia's second derivative differs from mgcv's own
  same-stencil answer by 1.76, 52 percent at the grid edge, and this
  package's agrees with it to 8.6e-05. So `eps` is scaled by the
  covariate's range and set per order.
* `frm_curve_feature()` locates a peak, a trough or a level crossing
  and gives its POSITION a standard error, by the implicit-function
  delta method: `var(t*) = var(f'(t*)) / f''(t*)^2` for a stationary
  point and `var(f(t*)) / f'(t*)^2` for a crossing. These are the
  numbers a movement paper reports, the time of peak velocity and
  movement onset, and they are usually reported without an interval
  because nothing hands one over.
* At a stationary point the curve's own value has a simplification
  worth knowing and the function reports it: the derivative of `f(t*)`
  with respect to the coefficients loses its `f'(t*) dt*/dc` term
  because `f'(t*)` is zero there, so the standard error of the PEAK
  HEIGHT is the ordinary pointwise standard error at `t*`, not inflated
  by the uncertainty in where the peak is.
* Every root the grid brackets is returned, so a curve with two peaks
  gives two rows and a curve with none gives a zero-row answer rather
  than an error. "This curve does not peak in this window" is an
  answer.

## The covariance these three needed, and the core seam

* All three need the same object, the covariance of the whole grid
  prediction, and frmtmb exports no route to it. A penalized smooth's
  wiggly part is a random-effect block even when the smooth is a
  population term, so the curve covariance needs the joint covariance of
  the fixed AND random coefficients. `vcov(full = TRUE)` excludes `b`
  under both of its branches, by construction: its documented invariant
  is that its row names are `confint()`'s. `predict(se.fit = TRUE)`
  forms the grid covariance internally and returns its diagonal.
* So the package rebuilds it, out of two facts. The linear predictor is
  LINEAR in the coefficients, so the difference between a prediction and
  the same prediction with one coefficient raised by one is that
  coefficient's design column, exactly and with nothing to tune. And the
  joint covariance is the inverse of the fit's own joint precision.
* Neither piece came from an exported function, so neither is trusted.
  Every call recomputes `sqrt(diag(Sigma))` and compares it with
  `predict(se.fit = TRUE)`, and refuses when they disagree by more than
  `tol`. Measured agreement on the package's own models is 3e-15 to
  9e-15 relative, which is machine precision, and it is reported by
  `print()` rather than hidden.
* **The covariance is core's own, read from its cache rather than
  recomputed.** `predict(se.fit = TRUE)`, which this package calls anyway
  for the check, memoizes the inverted joint precision on the fit; the
  curve then subsets it. Recomputing it instead, which the first draft
  did, paid for a second `sdreport()` per call, took a Schur complement
  over the coefficients the curve does not touch on a densified copy of
  a sparse matrix, and bypassed `autoscale_sdreport()`. Measured at 8006
  random coefficients: 114 s and 2.1 GB before, 8.2 s and 1.2 GB after,
  with the covariance subset itself at 0.00 s. An autoscaled fit now
  works rather than being refused, at 3.7e-16 against
  `predict(se.fit = TRUE)`.
* **The cost is the joint-precision solve, not the call count.**
  Measured at `re.form = NA` on a 20-point grid: design rebuild 0.01 s
  against `predict(se.fit = TRUE)` 0.29 s at 8 random coefficients,
  0.07 s against 0.98 s at 2006, and 0.28 s against 6.87 s at 8006. The
  term the call count measures is a tenth of the cost at every size.
  Size a job from the number of coefficients in the FIT, not from the
  grid.
* The design rebuild is one `predict()` call per contributing coefficient
  plus one probe per block of 24 that contributes nothing. Under the default
  `re.form = NA` a per-subject grouping block contributes nothing and is
  skipped in blocks rather than one coefficient at a time: measured, the
  vignette's 20-subject factor-smooth model costs 32 calls against 110
  random coefficients, and the count does not grow with the number of
  LEVELS. It is not the ideal 13 either, for two stated reasons: a chunk
  that straddles the boundary between a live block and a dead one is
  expanded whole, and a component short enough to skip the probe is
  expanded whether or not it contributes.
* `dev/spline-seam-proposal.md` says what core would have to export to
  make the rebuild unnecessary, and what a penalized coefficient block
  handed to a NONLINEAR body would need beyond that.

## The Royston-Parmar family

* `royston_parmar()` writes the log cumulative hazard as a natural cubic
  spline in log time, which is the flexible parametric survival model of
  Royston and Parmar (2002). `scale = "hazard"` is the proportional
  hazards version, `"odds"` the proportional odds version and
  `"normal"` the probit one; all three are here, because the second and
  third are one `switch` arm each once the first is written.
* It is parameterized exactly as `flexsurv::flexsurvspline()`
  parameterizes it, and that is testable rather than claimed. Taking
  flexsurv's fitted coefficients on its own `bc` data (686 rows, 299
  events) and evaluating frmtmb's objective at them reproduces
  flexsurv's log likelihood to between 1.4e-16 and 1.6e-15 relative
  across all three scales and 0, 1 and 3 interior knots. That is an
  identity, not an agreement. The two optima then match as well, and
  frmtmb's is never the worse of the two.
* `mu` is `gamma0` and the rest are `gamma1`, `gamma2`, and so on, all
  with identity links. A formula on `mu` is proportional hazards; a
  formula on any other coefficient is a TIME-VARYING effect, which is
  what flexsurv spells `anc =`. Each is an ordinary distributional
  parameter, so `s()`, a random effect and a prior all reach them.
* Knots go at equally spaced quantiles of the log UNCENSORED times, as
  Royston and Parmar place them and as flexsurv defaults. The quantiles
  need the response, so they are taken at frame assembly through
  `family_finalize()` and the family the fit carries has them baked in;
  `knots =` and `bknots =` pin them at construction instead. `df` counts
  interior knots plus one, so flexsurv's `k` is `df - 1`.
* The family declares an `lcdf` as well as an `lpdf`, so `cens()` and
  `trunc()` both reach the likelihood, and right, left and interval
  censoring all work. Both are registered CONDITIONAL rather than
  working, for the reason below. Note the coding: frmtmb's `cens()`
  reads 0 as an observed event and 1 as right censored, which is the
  OPPOSITE of a `Surv()` status column.
* **`rp_floored()` is the check this family cannot do without, and it
  refuses by default.** Two things in `royston_parmar()` are floors
  rather than answers and both are silent in the fitted object. This
  function recomputes them at the fitted parameters and stops, naming the
  row count, the maximum, the threshold, the reason and the remedy.
  `frm_curve()` and its two companions call it, so the documented way to
  inspect this family will not draw a curve off a fit whose likelihood
  is a floor artifact.
* The first floor is the one that made the family unshippable without
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
  reported log likelihood is from the model's depends on the data, not
  on the family: a floored row contributes -35.127363 instead of its own
  `-log S`, so the shortfall is about `-log S - 35` per row, and two
  runs of that design differing only in seed give 2.4e+03 and
  2.166e+04. The threshold is
  one quantity for all three scales, because the error depends only on
  `S`: `-log S` is `H` on `"hazard"`, `log(1 + exp(eta))` on `"odds"`
  and `-log(Phi(-eta))` on `"normal"`, and 19.2 is where `eps / S`
  passes 1e-8.
* The refusal is POST-FIT and cannot be otherwise. `logLik()` reads
  `object$opt$objective` directly and the family protocol has no hook
  that runs when a fit finishes, so nothing in this package can make
  `logLik()` or `AIC()` refuse on their own, and the optimizer may have
  walked through or stopped inside the flat region before the check is
  called. The real fix is a core one, an `lccdf` slot, and it is five
  sites rather than one: see `dev/spline-seam-proposal.md`.
* Monotonicity of the cumulative hazard is NOT enforced, and flexsurv
  does not enforce it either. What this family does instead is refuse to
  answer `NaN` for it: where the spline's derivative in log time goes
  non-positive, the log density is a large finite negative number, since
  `NaN` stops the optimizer and one `NaN` component poisons a mixture's
  log-sum-exp. That floor is NOT inert when it is used: it turns
  `logLik()` and `AIC()` into a pseudo-likelihood, measured at 6 floored
  rows and 3952 units on a cure-fraction dataset. `rp_floored()` counts
  those rows too and refuses on them. The floor is a smooth positive part rather than a branch,
  because RTMB refuses a comparison on an AD type outright and has no
  `CondExp`. It changes the log density by 2.4e-15 at a derivative of 1.
* That floor had to be written twice, and the measurement is what caught
  it. The textbook spelling, `0.5 * (u + sqrt(u^2 + eps2))`, is right on
  paper and wrong in double precision on the side it exists for: at
  `u = -35.75` the square root rounds to `|u|`, the sum cancels to
  exactly zero and `log()` returns `-Inf` after all. Measured, 647 of
  686 rows reached `-Inf` that way. Recovering the small branch through
  `(u + s)(s - u) = eps2` and ADDING it back does not help either,
  because adding 1.4e-16 to 71.5 loses it. The form that works SELECTS
  between the two cancellation-free branches with `sign()`, which is
  available on an advector and whose zero derivative is the correct one
  because the branches agree at `u = 0`.
* The spline basis is written branch-free, with the truncated power
  spelled `0.5 * (e + abs(e))`, so that it tapes if the response is ever
  promoted to a parameter. It is a function of the response, which is
  data, so nothing normally needs that; it costs one line to have it
  anyway.

## Deliberate omissions

* No `post$mean_fn`, so `fitted()` and `predict(type = "response")` are
  refused. The mean of a Royston-Parmar survival time has no closed
  form, and core's `cox()` refuses for the same reason. Read the fitted
  log cumulative hazard with `frm_curve()` instead.
* No exact basis derivative in `frm_curve_deriv()`. The design is
  rebuilt through `predict()`, which evaluates a basis and never
  differentiates one, so an exact derivative would mean reading mgcv
  smooth objects out of the fitted frame. That is a deeper reach into
  core than anything else here makes, for an error already at the tenth
  significant figure.
* Derivatives of order 3 and higher are refused rather than offered: a
  third central difference divides by `eps^3` and there is no step size
  at which it is accurate.
* No accurate log survival past `-log S` of about 19.2, and none at all
  past 36. See the `rp_floored()` entries above: the squeeze that keeps
  the term finite rather than `-Inf` cannot also keep it accurate,
  because the number core asks a family for is `F` and the complement of
  a double near 1 is not representable. flexsurv computes `log S`
  directly and has no such region. The package's answer is refusal by
  name rather than a silent floor, and the real fix is a core `lccdf`
  slot.
