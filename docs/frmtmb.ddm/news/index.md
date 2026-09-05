# Changelog

## frmtmb.ddm 0.2.0

Three families where there was one: Ratcliff’s full diffusion model as
an extension of
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md),
the generalized drift-diffusion model
[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md),
and the linear ballistic accumulator
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/lba.md),
plus the queued defects. The three share one floor idiom and one
compatibility table.

### The full diffusion model

- `wiener(variability = )` adds across-trial variability to the family
  rather than forking it. Naming any of `"sv"` (drift rate), `"sz"`
  (start point) and `"st"` (non-decision time) turns that one into an
  ordinary distributional parameter, with its own link and its own
  formula.
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
  with no arguments is the model it always was: on a grid of 11520
  parameter combinations its log density is byte-for-byte identical to
  the previous release everywhere the density is defined, and its tests
  are unchanged. The one difference is below the non-decision time,
  where it now returns `-Inf` instead of `NaN`, which is the deliberate
  fix below.
- The likelihood is the analytic Wiener density averaged over those
  distributions, and the three are done three different ways. The drift
  integral is Gaussian against an exponential-quadratic and is evaluated
  in CLOSED FORM: it agrees with adaptive quadrature of the same thing
  to better than 1e-13 relative and takes no nodes, so estimating `sv`
  is FREE relative to the plain density. Measured at the fit level over
  three seeds, `sv` costs 0.97 to 1.07 times plain, and at the density
  level both are 0.0002 s per call. The start-point and
  non-decision-time integrals are uniform and use fixed-node
  Gauss-Legendre quadrature, and those do cost: 0.0056 s per call for
  `sz` at 7 nodes and 0.0170 s for `st` at 21.
- Node counts are the `nodes` argument and the defaults are measured:
  `sz` reaches machine precision at 7 nodes, `st` reaches 1e-9 at 21.
  They differ because the non-decision-time range is cut by the response
  time on a fast trial, and the integrand turns on sharply at the cut.
  Node positions and counts are fixed when the family object is built,
  because an automatic-differentiation tape cannot record a branch on a
  parameter; a parameter only rescales the interval they map onto.
- Variability parameters at zero reproduce the plain Wiener density to
  better than 1e-13 in the log density, which is floating-point rounding
  on a differently associated sum rather than a quadrature error.
- Recovery on simulated data, 8 replicates of 1500 trials with a normal
  drift rate and a uniform non-decision time: every parameter’s Monte
  Carlo mean is within 2.2 of its own standard errors of the value it
  was generated from.
- Across-trial variability is NOT frmtmb’s `quadrature = TRUE`. That
  marginalizes random effects by Gauss-Kronrod, is wired to the
  random-effect coefficient vector by name, and refuses a model with no
  random-effect block. This integral shares nothing between trials, has
  no level to estimate, and exists in models with no grouping factor at
  all, so it lives inside the density.
- [`ddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/ddm_simulate.md)
  takes `sv`, `sz` and `st`, drawing each trial’s own parameters and
  then running the ordinary process, so the simulator states the model
  independently of the density.
- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) follow the
  variability rather than ignoring it. Conditioning on the boundary a
  row ended at reweights which per-trial parameters that row could have
  had, so the fitted mean is a ratio of two quadratures and the
  simulator accepts a drawn drift rate and start point with the boundary
  probability it implies. The plain closed form would not do: at an
  unbiased start point it returns the same mean for both boundaries,
  where 40000 simulated trials put them 0.06 s apart. Both are checked
  against simulated data, which knows nothing about how either is
  computed.

### The generalized drift-diffusion model

- [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)
  is the generalized drift-diffusion family of Shinn, Lam and Murray
  (2020): a drift that may depend on the accumulator’s own level and on
  a covariate, boundaries that may collapse within a trial, and a
  starting distribution that may be a point or an interval. There is no
  closed-form first-passage density, so every likelihood evaluation
  solves the Fokker-Planck equation forward in time and reads the
  probability flux through each boundary. It needs no change to core
  frmtmb.
- The components are chosen by argument and are extensible.
  [`gddm_drift_constant()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md),
  [`gddm_drift_coherence()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md)
  and
  [`gddm_drift_leak()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md)
  are summed to make a drift;
  [`gddm_bound_constant()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md),
  [`gddm_bound_exponential()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md)
  and
  [`gddm_bound_linear()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md)
  give the boundary;
  [`gddm_start_point()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-start.md)
  and
  [`gddm_start_uniform()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-start.md)
  give the start.
  [`gddm_drift_term()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_drift_term.md),
  [`gddm_bound_term()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_bound_term.md)
  and
  [`gddm_start_term()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_start_term.md)
  are the documented seams for writing more. Every free quantity is a
  dpar that takes a formula, as in
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md).
- `bs` is the boundary SEPARATION and `bias` the relative start point,
  both as in
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md),
  so estimates are directly comparable between the analytic family and
  the generalized one.
- The substitution `y = x / B(t)` pins the moving boundaries at plus and
  minus one, which keeps the grid fixed while the boundary collapses and
  is what makes the likelihood differentiable: nothing on the taped path
  branches on a parameter. With the walls stationary the scheme is
  Crank-Nicolson, which a solver that chases a moving bound cannot use.
- The likelihood is an ordinary rowwise family, not a
  `frmtmb_structure()`. frmtmb calls `lpdf` once per objective
  evaluation with full-length vectors, and the condition a trial belongs
  to is data, so a rowwise density does one solve per condition, which
  is all a structure would have bought, and it keeps
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) rather than
  defaulting them to refused. Measured: 6 solves for 2400 trials over 6
  conditions.
- [`gddm_control()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_control.md)
  carries the grid and what is done with the answer: `dt`, `ny`,
  `t_max`, `max_ndt`, `renormalize` and `tridiagonal`. Renormalizing the
  defective density is on by default and should stay on: the discretized
  solve loses mass in a parameter-dependent way, so a likelihood that
  does not divide it out rewards fast absorption. On data simulated from
  the model, turning it off more than doubles the fitted leak and
  shrinks the boundary separation by a fifth.
- `tridiagonal` picks how the solve inside each step reaches the tape.
  `"recorded"`, the default, builds a large tape that runs in compiled
  code; `"atomic"` collapses the solve into one node with a hand-written
  adjoint, building about twelve times faster and evaluating about
  twelve times slower. Both give the same derivative to machine
  precision.
- The published coherence nonlinearity has no derivative at zero
  coherence, which a motion design normally contains. The coherence is
  data, so the zero condition is resolved once when the tape is built
  and never reaches it; the gradient in the exponent is finite, and is
  exactly zero there, because the drift is zero whatever the exponent
  is. Signed coherences are supported for stimulus coding.
- The family admits exactly two responses, because one accumulator
  between two absorbing boundaries has two walls. A decision indicator
  with more than two levels is refused at frame assembly, naming how
  many levels the data has and pointing at
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/lba.md),
  which fits the racing accumulators that more than two alternatives
  need, rather than being folded into one of the two.
- The boundary is read from `dec()` when it is there and from `vint()`
  otherwise, so the spelling brms uses works on this family as it does
  on
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md).
  `vint()` numbers its values positionally, so the condition index is
  the first `vint()` value alongside `dec()` and the second inside
  `vint(upper, cond)`; the family reads whichever it is. Neither can be
  declared through `required_aterms`, which names the terms a density
  needs ALL of, so both refusals are written out.
- [`gddm_conditions()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_conditions.md)
  builds the condition index.
  [`gddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_simulate.md)
  draws choices and response times from the model’s own solved density.
- Validated in the package’s own suite: against this package’s Wiener
  density with a constant drift and fixed bounds, to better than 0.01 in
  the log density at the shipped grid over decision times from 0.2 s on,
  degrading on a coarser grid and improving on a finer one; automatic
  gradients against numDeriv at parameter points including a collapsing
  boundary; parameter recovery with Monte Carlo standard errors; the
  size of the renormalization bias; and the zero-coherence gradient.
- The density is floored before it is logged, as the Wiener density in
  this package already was. Where the solved density at a trial’s own
  response time underflows, the log density is a large finite negative
  number instead of `NaN`: the optimizer gets a value it can use, and a
  mixture’s log-sum-exp is not poisoned by one component. The floored
  row is flat, so its gradient is exactly zero; `-Inf` would not do,
  because `-Inf` differentiates to `NaN`.
- `gddm_floored(fit)` is where that goes to be read. It returns the
  number of rows answered by the floor rather than by the solver at the
  fitted parameters, with the row indices in the `"rows"` attribute.
  Zero is the ordinary case and means the grid represented every
  observation. A few rows means a few trials sit within a few time steps
  of the fitted non-decision time, where a fixed grid cannot resolve a
  density climbing through orders of magnitude. Many rows means the fit
  is not to be trusted: shrink `dt`, or add a lapse component. The count
  replaces what used to surface as repeated optimizer warnings about
  `NaN` function evaluations, which a user could not act on and which
  masked real warnings in a test run.
- Known limits, measured and stated rather than hidden: the density at
  decision times of only a few time steps is far larger than the truth,
  because an implicit scheme spreads a little mass everywhere at once
  where the true density is exponentially small. A lapse component
  (`gddm(lapse = "uniform")`) floors it in the model rather than in the
  arithmetic. Cost scales with the number of conditions, so a design
  with many distinct parameter settings is where this becomes painful.
- [`vignette("gddm")`](https://aforren1.github.io/frmtmb/frmtmb.ddm/articles/gddm.md)
  fits one of the paper’s models end to end and says plainly where this
  is slower than
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
  and why someone would pay that.

### The linear ballistic accumulator

- `lba(n)` is a race between `n` accumulators, each rising in a straight
  line from a uniform start point on `(0, A)` to a common threshold at a
  normally distributed rate. It is the family for choices with more than
  two alternatives: a diffusion between two absorbing boundaries admits
  exactly two responses, and no reparameterization of
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
  reaches a third. The likelihood is closed form for any `n`.
- Each accumulator’s drift mean is its own distributional parameter
  (`v1`, …, `vn`, identity link) and all of them are primary, so the
  main formula reaches every drift and each gets its own coefficients.
  Give one its own formula and a covariate moves that alternative alone,
  which is the capability neither two-choice family can offer.
- The remaining parameters are `A` (start-point range, log), `k`
  (threshold above the start-point range, log) and `ndt` (non-decision
  time, a logit scaled onto `(0, max_ndt)` as in
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)).
  The threshold is `A + k` rather than a parameter of its own, so
  `b > A` is structural: a threshold inside the start-point range, where
  trials would begin already finished, is not a state the optimizer can
  reach.
- The drift standard deviation is fixed, because the model is identified
  only up to a common rescaling of `A`, the threshold, the drifts and
  the drift standard deviation. It is the family argument `sd_v` rather
  than a hidden default, and is carried on the family object.
- Drift rates are truncated at zero by default, following `rtdists` and
  its `posdrift = TRUE`. Every accumulator then arrives eventually and
  the choice probabilities sum to one. `posdrift = FALSE` gives the
  untruncated convention, under which the response distribution is
  defective; the two are different models, not a rescaling of each
  other, which matters when comparing against another package.
- The single-accumulator density agrees with
  [`rtdists::dlba_norm()`](https://rdrr.io/pkg/rtdists/man/single-LBA.html)
  to better than 1e-11 relative wherever `rtdists` is itself accurate,
  and the race with
  [`rtdists::n1PDF()`](https://rdrr.io/pkg/rtdists/man/LBA-race.html)
  for two, three and four accumulators. In the fast tail the two diverge
  by design and this one is the better: `rtdists` writes the normal
  difference `Phi(g) - Phi(h)` as a subtraction of two lower tails,
  which returns exactly zero once `pnorm(h)` saturates. Written in log
  space instead, as `exp(la) * -expm1(lb - la)` on the two upper-tail
  logs, nothing saturates and no comparison is needed, which matters
  because RTMB refuses comparison on AD types. Against a 200-bit Rmpfr
  reference on 1144 points, the subtractive form returns exactly zero on
  68 of the 80 tail rows and is already 2.7e-3 wrong on 8 rows outside
  the tail; the log-space form holds 5.0e-14 in the bulk and 3.6e-4 in
  the tail. Four rows still exceed 1e-6, so it is an improvement, not a
  proof.
- That mattered beyond accuracy. The subtractive form errs in the value
  only, and a tape differentiates the function that was written, so
  value and gradient described different surfaces near the
  non-decision-time bound. It did not move the point estimate, which the
  likelihood keeps out of that region, but Hessian-based standard errors
  were up to 16.5 percent off before the change.
- The survival function is written out directly rather than as
  `1 - plba_norm()`, which returns exactly zero, and so a
  log-contribution of `-Inf`, while the true survival is still around
  1e-19; the direct form agrees with a quadrature of the density down to
  survivals of 1e-23. Under `posdrift = FALSE` the never-arriving mass
  it adds back is spelled `Phi(-v/s)` rather than `1 - Phi(v/s)`, which
  is the same saturation trap one line further on.
- The choice reaches the density through `vint()`, as the Wiener
  family’s boundary indicator does, and is declared with
  `required_aterms`, so omitting it is refused by name. A choice outside
  `1..n`, a non-positive response time and a non-decision-time bound
  above the fastest response are each refused with their own message.
- The non-decision-time bound is derived from the response through
  `family_finalize()`, so the family object carries the link it will
  actually use rather than an environment filled in later.
- [`lba_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/lba_simulate.md)
  draws choice and time jointly from the generative process.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) on a fitted
  object redraws times conditional on each row’s observed choice, since
  the choice is data.
- Deliberate omissions: no `lcdf`, so `cens()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) are refused; no
  `post$mean_fn`, because the mean of the race has no closed form and a
  quadrature per row was not worth writing for
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

### Defects fixed

- **The density returns `-Inf` below the non-decision time, not `NaN`.**
  A likelihood of zero is the right answer there; `NaN` is not an answer
  and propagates through every component of a mixture. The normalized
  time is held at the smallest positive double, which is bit-for-bit
  inert everywhere the density was already right.
- **`wiener(max_ndt = )` above the smallest response time is usable in a
  mixture.** The refusal is correct for the family on its own and wrong
  inside a mixture, where the contaminant component is exactly what
  covers the trials the diffusion cannot produce, so
  `allow_unreachable = TRUE` lifts it and the refusal now names it. A
  trial below the non-decision time then gets a log density that
  exponentiates to zero AND differentiates to zero, which a true `-Inf`
  does not: `-Inf` yields a `NaN` gradient and stops the fit.

### What frmtmb 0.49.0 let this package delete

- **`dec()` is the spelling now.** frmtmb gained
  `frmtmb_register_aterm()`, this package registers `dec` when it loads,
  and `rt | dec(response) ~ x` works and takes a factor, a character
  vector or a logical the way brms does. `vint()` carries the same thing
  as a 0/1 integer and is unchanged.
- **The environment the link closures read is gone.** The bound on the
  non-decision time is a property of the response, and the family object
  is built before `frm()` has any; this package used to have `valid_y()`
  write the bound into an environment, which worked only for as long as
  an undocumented slot order held. `family_finalize()` is the documented
  slot for it and the family now derives itself from the data there.
- One hand-rolled check remains, and is not `required_aterms`’s fault:
  that argument names the terms a density needs ALL of, and this family
  needs EITHER `dec()` or `vint()`. See `dev-findings.md`.

## frmtmb.ddm 0.1.0

First release. A Wiener first-passage time family for two-choice
response times, written entirely against frmtmb’s exported extension
API.

- [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
  is the drift-diffusion family, with brms’s dpar names and links: `mu`
  (drift rate, identity), `bs` (boundary separation, log), `ndt`
  (non-decision time, bounded) and `bias` (relative start point, logit).
  A model written for
  [`brms::wiener()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  reads the same here, and the parameterization is pinned against brms
  in the test suite.
- The density is the Navarro and Fuss (2009) pair of series. Both are
  evaluated at a fixed truncation and their logs blended with a logistic
  weight in the normalized time, because an AD tape cannot choose
  between them on a parameter. It agrees with
  [`RWiener::dwiener()`](https://rdrr.io/pkg/RWiener/man/wienerdist.html)
  to better than 1e-11 relative over normalized times from below 1e-3 to
  50, where a fixed truncation of the small-time series alone is wrong
  by tens of percent past about 8.
- The decision indicator reaches the density through `vint()`, because
  frmtmb’s addition terms are a closed set and brms’s `dec()` is not one
  of them. Omitting it is refused with a message naming both spellings,
  rather than silently producing a log likelihood over no rows.
- `ndt` uses a logit scaled onto `(0, max_ndt)`, so the support
  constraint is structural rather than something the optimizer has to
  discover. `max_ndt` defaults to the smallest response time in the
  data, found at frame assembly.
- `post$mean_fn` gives the mean response time in closed form,
  conditional on the boundary the row ended at, with the zero-drift
  limits handled explicitly. `sim` draws conditionally by inverse
  transform through
  [`RWiener::qwiener()`](https://rdrr.io/pkg/RWiener/man/wienerdist.html),
  with a discretized forward simulation as the fallback when RWiener is
  absent.
- [`ddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/ddm_simulate.md)
  draws response times and boundary choices jointly, for building
  example and test data sets.
- The package registers its compatibility rows with the core at load
  through
  [`frmtmb::frmtmb_register_compat()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_compat.html),
  so `frm_compat("wiener")` states what was exercised and what was not.
- Deliberate omissions: no `lcdf`, so `cens()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) are refused; no
  variance function and no unit deviance, so
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) answers
  `type = "response"` only.
- `dev-findings.md` in the package source records what building this
  from outside frmtmb cost, as an acceptance test of the extension API.
