# frmtmb.eam (development version)

* **`gddm()` now refuses a parameter that varies inside a condition,
  instead of ignoring it.** One solve of the Fokker-Planck equation
  serves a whole condition and every parameter is read at that
  condition's FIRST ROW. `?gddm` has always required every row of a
  condition to share every parameter value, and nothing enforced it, so
  a model that broke the contract fitted with no error and no warning
  and the varying term reached nothing. Measured on 120 rows in two
  conditions, at a FITTED parameter vector, adding 100 to a covariate
  on 118 of the 120 rows:

  | model | objective as drawn | with the covariate moved | identical |
  |---|---|---|---|
  | `mu ~ x` | 18.962737109172043 | 18.962737109172043 | yes |
  | `ndt ~ x` | 18.810414088214987 | 18.810414088214987 | yes |
  | `bs ~ x` | 18.948550650674157 | 18.948550650674157 | yes |
  | `bias ~ x` | 17.746051487984953 | 17.746051487984953 | yes |

  Bitwise identical, on every one. A random effect whose grouping
  crosses the index fails the same way from the other side: with four
  subjects across two conditions, the two subjects that never sit on a
  condition's first row keep deviations of exactly 0, and the fit
  reports a between-subject standard deviation of 7.6e-11 rather than
  refusing.

  The refusal names the parameter, the variable and a condition it
  varies inside, and the remedy is to name that variable in the index
  as well, which `gddm_conditions()` does. The cost of the extra
  resolution is one solve per condition, and it is why this is a
  refusal rather than a silent refinement: at 120 rows on a coarse
  grid, one solve per row against one per condition is 41 times the
  tape build and about 52 times the evaluation.

  The check compares each parameter's model-frame columns against its
  condition's first row, so it is sufficient rather than exact, in the
  direction that refuses. It carries a tolerance,
  `1e-8 * max(|a|, |b|) + 1e-11 * max|column|`, because `poly()`
  returns rows 8.7e-14 apart, relative to the column maximum, for
  bitwise identical inputs. The band under that tolerance is the guard's
  limit and it is stated rather than hidden: the check misses a
  within-condition difference `s` on a column whose largest entry is
  `M` once `M / s` passes about 1e11, which takes a sentinel value a
  thousand billion times the data to reach.
  On 21 correct designs covering drift by coherence, boundary
  by block, start point by cue, non-decision time by subject, subject
  random intercepts and slopes, a spline, a monotonic ordered predictor,
  an offset, a collapsing boundary, a lapse rate and the coherence drift
  nonlinearity, it fires 0 times; on the 26 ways of dropping one
  variable from those same indices it fires 26 times and names the
  dropped variable every time. One refusal names every offending
  parameter and variable and hands back the `gddm_conditions()` call
  that fixes them.

  A random effect is where this bites hardest, and `frm_compat("gddm")`
  now says so: `(1 | g)` works when `g` is in the condition index and is
  refused when it is not. Before the refusal, four subjects across two
  conditions left the two not on a condition's first row with deviations
  of exactly 0, and the fit reported a between-subject standard
  deviation of 7.6e-11.

* **`gddm_simulate()` refuses a parameter that varies between trials
  sharing a coherence, instead of drawing every trial from the first
  one.** It solves once per distinct value of `coh` and read every
  parameter at that value's first trial, which is the same defect one
  layer over. Measured on 400 trials with
  `mu = c(rep(-2.5, 200), rep(2.5, 200))` and `coh = 0`: both halves
  were drawn from `mu = -2.5`, giving upper-boundary rates of 0.010 and
  0.000, where the same values with `coh` separating them give 0.010 and
  0.990. No error and no warning. Give `coh` a distinct value per
  parameter setting, or call the function once per setting. A scalar
  parameter, and a length-`n` one that is constant within each `coh`,
  are unaffected.

# frmtmb.eam 0.7.0

The non-decision time can now be bounded PER GROUP, which fixes the
defect 0.6.0 disclosed. A model that does not ask for it is unchanged.

* **`ndt_group()` is a new addition term: the grouping the bound is
  taken per.** Write `rt | dec(r) + ndt_group(subject) ~ ...` and each
  row's non-decision time is bounded by its own subject's fastest
  response instead of by the whole data set's. That is what a random
  effect on `ndt` needs, and 0.6.0 disclosed what happens without it.
  On the Phase 0 design of `dev/extension-gaps-plan.md`, 30 subjects by
  400 trials, the same data and the same model:

  | | 0.6.0, one global bound | with `ndt_group(s)` |
  |---|---|---|
  | convergence code | 1 | 0 |
  | maximum absolute gradient | 1.25e11 | 9.9e-04 |
  | positive definite Hessian | no | yes |
  | standard errors that are `NaN` | 7 of 7 | 0 of 7 |
  | log-likelihood | -7148.81 | -7003.01 |
  | population `ndt`, truth 0.25 | 0.2261 | 0.2469 (se 0.0020) |
  | condition effect, truth 0.9 | 0.891, no interval | 0.911 (0.858, 0.964) |
  | whole `frm()` call | 269 s | 108 s |

  145.8 log-likelihood units better at the same parameter count, and
  faster. At the optimum all 30 subjects sit below their own fastest
  response, the tightest by 27.3 ms.

  **The argument that settles it is what happens as data accumulates.**
  Same design, same truths, at 100, 200 and 400 trials per subject, the
  per-subject root mean squared error of the fitted non-decision times:

  | trials | with `ndt_group(s)` | one global bound |
  |---|---|---|
  | 100 | 20.91 ms | 27.23 ms |
  | 200 | 14.08 ms | 20.51 ms |
  | 400 | **7.67 ms** | **30.37 ms** |

  The per-group bound converges on the truth and the global bound does
  not, because more data lowers the global minimum and tightens the
  ceiling on every subject at once.

  Give the grouping as a factor, a character vector, a logical or
  integer codes: it is keyed on the group's LABEL, so subsetting,
  `droplevels()`, `relevel()` and a prediction grid you build yourself
  all pair each row with the same bound the fit used. A group the fit
  never saw is refused rather than given the global bound.

* **Under `ndt_group()`, and only there, `ndt` is a FRACTION of the
  row's own bound rather than a time.** A per-row bound cannot live in
  a link, so with a grouping the `ndt` link is a plain logit and the
  density multiplies. What that changes, for a model that uses the new
  term: `predict(fit, dpar = "ndt", type = "response")` returns a
  number in `(0, 1)`, a `prior(class = "ndt")` is a density on that
  fraction, and a `bf(ndt = )` constant is a fraction. The new
  `ndt_time()` returns the non-decision time in the units of the
  response for either parameterization and is the call to reach for.
  `st` follows `ndt`: with a grouping it is a fraction of twice the
  row's bound.

  **Without `ndt_group()` nothing moved.** The bound is one number and
  stays in the link exactly as before, so `ndt` is still a time,
  `predict(dpar = "ndt", type = "response")` still reports seconds, a
  ported `prior(normal(0.30, 0.01), class = "ndt")` still means 300 ms,
  and `bf(ndt = 0.2)` is still 0.2 s and is still refused when the
  family has no bound to measure it against.

* `ndt_time(fit, newdata = )` reports the non-decision time in the
  units of the response, for `wiener()`, `lba()`, `rdm()` and
  `wiener_gng()`. The bound each row was measured against is on the
  fitted family, at `family(fit)$ndt_bound`, with the trial count per
  group beside it. The bound it uses is a property of the data the model
  was FITTED to, so a prediction on new rows uses the bound the fit
  used.

* `max_ndt` still means one absolute upper bound applied to every row,
  and a `max_ndt` above the fastest response is still refused outside a
  mixture. Combining it with `ndt_group()` is refused: the two set the
  same bound to different things. An `ndt_group()` no family reads is
  refused too, which is what a grouping inside a `mixture()` would be,
  because a mixture never finalizes its components.

* **On upgrading a model that does not use `ndt_group()`**: nothing
  changes. The bound stays in the link, where 0.6.0 put it, so this is
  not an equivalent parameterization but the same one. The objective
  and its gradient are BITWISE identical to 0.6.0's, 0 ulp over 50
  fixed-parameter probes taken at the starting values and off the
  optimum, across `wiener()` plain and with `max_ndt`, with `st` and
  with `sv`, `sz` and `st` together, `rdm()` plain and censored,
  `lba()`, `wiener_gng()` plain and with `st`, and `gddm()`. Refitting
  fourteen models and reading every reachable quantity off each gives
  119 of 120 identical, the one difference being the wording of an
  error message.

* **`sd(ndt)` reads as recovered and the floors do most of the work.**
  The table above is a real improvement and this is the caveat that
  belongs beside it: because a grouped model estimates `ndt` as a
  fraction of each group's own floor, the fitted per-subject times vary
  with the floors even when the variance component is zero. At the
  design above an estimator with NO random effect on `ndt` returns
  `sd(ndt)` = 0.02748 against a truth of 0.02629, where the full model
  returns 0.02543, and at 100 trials per subject the two are identical
  in every digit. The variance component is not empty at 400 trials, it
  buys 8.44 log-likelihood units and cuts the per-subject RMSE from
  11.90 ms to 7.67 ms, but `sd(ndt)` is the wrong statistic to read
  that off. Read the per-subject error and the log-likelihood instead.

* `gddm()` does NOT take `ndt_group()`, and keeps the single scaled
  logit, so its `ndt` is a time. The exclusion is on SCOPE: its solver
  reads every parameter at the first row of each condition and `?gddm`
  already requires every row sharing a condition to share every
  parameter value, so a valid model whose non-decision time varies by
  subject already carries a condition per subject and a per-condition
  bound would reach the density exactly as `ndt` does. What it would
  cost is a Fokker-Planck solve per subject. `frm_compat("gddm")` says
  so.

# frmtmb.eam 0.6.0

A random effect on the non-decision time is broken, and with
`variability` set it is broken SILENTLY. Read the first bullet before
fitting a hierarchical DDM.

* `wiener()` bounds the non-decision time by `min(rt)`, the GLOBAL
  fastest response in the data, through a scaled logit. With a random
  effect on `ndt` this is the wrong constraint: the information about
  a subject's non-decision time is that subject's own fastest
  response, and a subject whose true `ndt` is above the global minimum
  cannot be represented at any value of the random effect. At 30
  subjects by 400 trials with a between-subject `ndt` spread of 26 ms
  on a mean of 250 ms, 20 of the 30 subjects are in that position
  while NONE is inconsistent with its own data. Without `variability`
  the fit does not converge and says so: maximum absolute gradient
  1.3e11, Hessian not positive definite, all seven standard errors
  `NaN`. With `variability = "sv"` it CONVERGES: code 0, maximum
  gradient 7.5e-05, positive definite Hessian, no bad standard errors,
  and `diagnose()` reports nothing. The population non-decision time
  then comes back pinned at the bound, 0.2236 against a truth of 0.25
  and 0.42 of a standard error below `min(rt)`, with a delta-method
  standard error of 7.2e-06 on it, because the scaled logit's
  derivative vanishes where the estimate has been pushed. Do not put a
  random effect on `ndt` until the per-subject bound lands; the same
  defect reaches `rlddm()` in frmtmb.learn, which takes its diffusion
  parameterization from this package. Fitting the same data with the
  bound raised above every subject's truth recovers everything and
  finds a log likelihood 121.4 units higher at the same parameter
  count.

* `valid_y` warns when the fastest response exceeds 20 seconds and
  names milliseconds as the likely cause, across all five families.
  Zero false alarms over 192 designs a two-choice task produces. Where
  it can fire on a correct model the help reports the rate with the
  fastest response beside it, because the rate is not a function of
  the median: at about 50 seconds it spans 0.185 to 0.935 depending on
  drift.

* **A units guard.** Every default in this package reads the response as
  a time in SECONDS: the starting values, the bound the `ndt` link is
  scaled onto, `gddm_control(dt = 0.01)` and the window `t_max` takes
  from the data. Nothing in any of the five likelihoods refused
  milliseconds. The fit converged and reported a boundary separation
  three orders of magnitude out, which is the silent wrong answer this
  project ranks first. `valid_y` now warns when the fastest response in
  the data is above 20, names milliseconds as the likely cause and says
  to divide by 1000. All five families raise it, in their own name.

  Where 20 comes from. It is a ceiling on the fastest response in the
  WHOLE data set rather than on any one trial, so one slow trial does
  not reach it and cannot. Over 81 cells at the standard published range
  (drift 0.5 to 3, boundary 0.8 to 2.5, non-decision time 0.15 to 0.6,
  at 200, 400 and 12000 trials) the fastest response ran from 0.157 to
  0.801 seconds: 0 of 81 false alarms, and 81 of 81 caught when the same
  data is read as milliseconds. In the other direction the miss rate is
  0 of 81 at every threshold up to 100 and 14 of 81 at 200, so 20 sits a
  factor of five inside the band where nothing is missed.

  **Where it CAN fire on a correct model**, which the first sweep hid.
  That sweep held the trial count at 200 and up, where the fastest
  response is pinned just above the non-decision time and the guard
  reduces to "is the non-decision time above 20 seconds". A SHORT
  session of SLOW decisions is the exposed case, because the minimum of
  a few draws sits far above the floor. Measured with `ddm_simulate()`
  at a drift of 0.18 and a non-decision time of 1.5 seconds, 200
  replicates a cell. The FASTEST response is in the table beside the
  median, because the fastest is what the guard reads:

  | median | 20 trials | 60 trials | 200 trials |
  | --- | --- | --- | --- |
  | 28 s | 0.000 (9.2 s) | 0.000 (7.3 s) | 0.000 (6.0 s) |
  | 39 s | 0.045 (13.4 s) | 0.000 (10.3 s) | 0.000 (8.5 s) |
  | 45 s | 0.155 (15.7 s) | 0.000 (12.2 s) | 0.000 (10.1 s) |
  | 61 s | 0.690 (22.7 s) | 0.285 (18.2 s) | 0.000 (15.1 s) |
  | 100 s | 1.000 (41.7 s) | 0.995 (34.8 s) | 0.995 (29.1 s) |

  A 20-trial session with a median near 40 seconds is a deliberation,
  insight or matrix-reasoning design, not a stress test, and its
  non-decision time is 1 to 2 seconds rather than 20.

  **The median does not determine the rate**, which is why the fastest
  response is in the table. Holding the median near 50 seconds and
  reaching it three ways at 20 trials: drift 0.10 and boundary 18.7
  gives a median of 57.0 s, a fastest of 15.7 s and a rate of 0.185;
  drift 0.18 and boundary 21.5 gives 49.4 s, 17.5 s and 0.260; drift
  0.35 and boundary 37.6 gives 52.4 s, 26.9 s and **0.935**. A task
  that is slow because the boundary is far and the evidence is strong
  has a tight response time distribution, so the fastest of twenty
  trials sits close to the median; a task that is slow because the
  evidence is weak has a long right tail and a fast minimum. At one
  median the rate spans 0.185 to 0.935, and it tracks the fastest
  response throughout.

  So the honest one-line cost is the one the guard implements: it fires
  when the FASTEST response passes 20 seconds, and how far the fastest
  sits below the median depends on the spread as much as on the
  center.

  Nothing a standard two-choice task produces reaches the ceiling. Over
  192 cells at boundary separations of 2 to 5, non-decision times of 0.5
  to 5 seconds, drifts of 0.3 to 2 and 12 to 80 trials, the fastest
  response ran from 0.599 to 6.605 seconds and there were 0 false
  alarms.

  A warning and not a refusal for exactly that reason, and it carries
  the class `frmtmb_eam_units_warning`, so a genuinely slow design can
  silence this one condition without also hiding the convergence
  warnings beside it. It is raised once per FIT: counted over eleven
  entry points, `frm()` raises it once and `update()` once, because
  `update()` reassembles the frame, and the nine post-fit methods raise
  it zero times.

  **What it misses**, stated because it is measurable: a millisecond
  record with fast-guess contamination in it. At 12000 rows with 5
  percent of trials drawn uniformly over the observed range, only 2 of
  30 were caught, because one contaminant below 20 ms hides the whole
  data set.

  That is not softened by the record being broken in seconds too. On one
  such record the SECONDS reading raises no warning of any kind, reports
  "No convergence problems detected" with a maximum gradient of 1.2e-05
  and a positive definite Hessian, and estimates a boundary separation
  of 2.06 against a truth of 1.4, which is 47 percent out. So it is a
  second SILENT wrong answer rather than a visible failure, and the
  right closure for it is the contaminant mixture of item 3.5 rather
  than a wider units guard. A median-based arm would catch the
  millisecond half, because a millisecond median sits near 600 and
  contamination cannot move a median; it was not shipped because it was
  not measured for false alarms, and an unmeasured second heuristic
  inside a guard is what this project's rules forbid.

# frmtmb.eam 0.5.1

Requires frmtmb 0.55.0, for the hazard-container lint that now
runs in this package's own check.

# frmtmb.eam 0.5.0

`wiener()` refuses a boundary given twice, one export for a
sibling package, and the compatibility table says refused where it
used to say untested. Requires frmtmb 0.54.0 for the exclusivity
declaration.

* `wiener()` declares `dec()` and `vint1` mutually exclusive, so a model
  supplying both is refused by name instead of fitted with the second
  column unread. Measured before the change on 120 rows:
  `rt | dec(u) ~ 1`, `rt | vint(u) ~ 1` and
  `rt | dec(u) + vint(1 - u) ~ 1` all gave a log-likelihood of
  -76.0486443897369, the third with the two columns CONTRADICTING each
  other. `ddm_indicator()` reads `dec` and falls back to `vint1`, and
  the declaration is that precedence written down. Needs frmtmb's new
  `frmtmb_family(exclusive_aterms =)`.

  `gddm()` deliberately declares no such set and is why the rule is
  opt-in: there `dec()` and `vint1` are two data, the boundary and the
  condition index, not two spellings of one.

* Every family x addition-term cell in `frm_compat()` is now decided.
  The refusals each family has already declared in `accepts_aterms` are
  derived through frmtmb's `compat_aterm_rules()` rather than written
  out again, so `trials()`, and `vreal()` where the family does not read
  it, stop reading `untested`. Rows written by hand keep their own
  notes: the derivation defers to any pair already refused.

  The allow-lists themselves move to one `ddm_accepts` list that the
  five constructors and the compatibility rows both read, so the table
  cannot promise a term frame assembly refuses.

* `wiener_lpdf()` is now exported: the Wiener first-passage log density
  with `wiener()`'s parameterization and `wiener()`'s tape safety.
  `frmtmb.learn::rlddm()` is a delta learning rule whose value
  difference drives the drift rate of this density, and until now the
  only route to it was `frmtmb.eam:::ddm_lpdf_both()`, which is a
  promise nobody made. This is that promise, made deliberately and kept
  to one function: the series truncations, the blend between them, the
  across-trial variability integrals and the CDF all stay internal.
  Nothing about the package's own behavior changes.

# frmtmb.eam 0.4.0

`wiener_gng()` gains across-trial variability with the go branch
identical to `wiener()`'s; `rdm()` and `wiener_gng()` take `cens()`;
every family declares the addition terms it reads. Requires frmtmb
0.53.0 for the allow-list.

* `wiener_gng()` gains `variability =`, and takes Ratcliff's `sv`, `sz`
  and `st` under the same names, links and argument [wiener()] takes
  them. The go branch is [wiener()]'s averaged density and is not merely
  equivalent to it: evaluated at the same parameters, the two families'
  log densities are BIT-IDENTICAL on every row, for every combination of
  the three. Making that true is why the family now has a builder that
  is handed the non-decision-time bound and the unreachable-row margin,
  as `wiener()` has; at 0.3.0 it floored the decision time at a flat
  1e-12 where `wiener()` used `1e-9 * min(y)`, and no identity survives
  two different margins.

  The no-go branch is the part that had to be written, because it
  averages a DISTRIBUTION FUNCTION rather than a density. `sz` goes
  through Gauss-Legendre nodes, `st` shifts the deadline rather than
  cutting a range, and `sv` needs a Gauss-Hermite quadrature of its own:
  the drift enters the density only as an exponential-quadratic, which a
  normal average integrates in closed form, and it enters the
  distribution function through the eigenvalues as well, where nothing
  does.

  Verified against a 200-bit `Rmpfr` integration of the distribution
  function, against 40000-trial simulation from the process (the
  observed no-go rate is within 1.7 standard errors at every setting,
  and a ten-cell chi-square on the go response times reaches 16.0
  against a 0.999 critical value of 27.9), and against EMC2, whose
  `DDMGNG` carries the same three: the no-go probability matches
  `1 - EMC2:::pDDM()` to 3.9e-13 or better and the go density matches
  `EMC2:::dDDM()` to 1.6e-15. EMC2's `st0` is a uniform on
  `[t0, t0 + st0]` where this family's `st` is centred on `ndt`,
  following brms; the comparison applies that shift and `?wiener_gng`
  records it.

* `wiener_gng(nogo_nodes =)` is a second node-count argument, because
  the two branches are two integrands. The density's `st` count is 21
  because its range is cut at the response time and its integrand turns
  on sharply at the cut; the probability has no cut and saturates at 7.
  That would be a curiosity if the two cost the same, and they do not:
  the probability's rule is a three-dimensional product, so a count
  carried over from the density is multiplied by every other count in
  the grid. With the density's own counts a 500-row three-variability
  model exhausted memory outright.

* **`sv` is weakly identified in a go/no-go design, and `?wiener_gng`
  now says so.** On 3000 simulated trials with a true `sv` of 0.6 this
  family returns 0.001 while `wiener()`, on the same generative
  parameters and seeing both boundaries, returns 0.612. It is not an
  optimizer failure: the fit reaches a HIGHER log likelihood at `sv`
  near zero by trading it against the drift and the boundary
  separation. Read a small fitted `sv` here as "the data did not pin
  it".

* A start point that `sz` pushes past a boundary is now the boundary
  case it is, rather than a `NaN`. `wiener()` documents that a wide `sz`
  at a biased start can leave the boundaries and calls the density there
  a barrier; that is true of the density, which stays finite, and it was
  NOT true of the no-go probability, which took `log1p(-w)` of a
  negative above one and returned a probability ABOVE one below zero. A
  `NaN` is not a barrier, it is the end of the tape, and one node of one
  row took the whole fit with it.

  The clamped value is the correct limit for THE BRANCH IT IS APPLIED
  TO. It does not repair the pair, and the documentation no longer says
  it does: the go branch is left unclamped so that the bit-identity with
  `wiener()` survives, so past a boundary the two branches average
  different start-point distributions and the go mass plus the no-go
  probability falls short of one - measured, 0.970639 at `sz` = 0.5 with
  `bias` = 0.85, and 0.843370 at `sz` = 0.9 with `bias` = 0.90, so up to
  16 percent of the mass. It is tolerable because the region is strictly
  downhill: profiled at `bias` = 0.85, the best point inside the
  boundaries beats the best point outside by 189 log units and the
  surface is monotone across the crossing, so the clamp is a barrier an
  optimizer walks away from rather than a corner it can be pulled into.

  Inside the boundaries, where the model is defined, the two branches
  still sum to one - to 8.9e-16 in the plain family, and otherwise to
  whatever the quadrature gives, which is 3.9e-14 at `sv` = 0.6 and
  8.1e-05 at `sv` = 2.0 with default nodes. `?wiener_gng` carries both
  tables.

* `rdm()` declares an `lccdf` and an `lcdf`, so `cens()` and `trunc()`
  both work through core 0.52.0's slot. Neither needed new algebra: a
  race is unfinished exactly when every accumulator is, so the log
  survivor is the sum of the same per-accumulator survivals the density
  already forms for the losers of an observed trial, over all `n`
  instead of `n - 1`. Checked against the likelihood written out by
  hand at the fitted parameters: 2.2e-15 relative on a right-censored
  data set, 5.0e-16 with all four censoring codes, and 9.4e-16 on a
  left-truncated fit. A censored row still needs a `vint()` winner,
  which the likelihood does not read, and does not read EXACTLY:
  moving every censored row's winner changes the log likelihood by
  zero.

* `wiener_gng()` declares an `lccdf` and, deliberately, no `lcdf`. Right
  censoring is the same statement the family already makes, so it is
  exact rather than close: the same rows scored as no-go trials at the
  deadline and as trials right-censored at the deadline give log
  likelihoods that differ by no bits at all. Left censoring, interval
  censoring and `trunc()` are refused by name, because this likelihood
  is a defective density plus a point mass and a window normalizer on
  the response scale would renormalize the density while saying nothing
  about the mass.

* `wiener()` and `gddm()` declare the decision indicator through core's
  any-of `required_aterms` instead of checking for it by hand. Both read
  the boundary from `dec()` or from `vint1`, which core 0.51.0 spells
  `list(c("dec", "vint1"))`, so frame assembly refuses a model that
  supplies neither BEFORE the frame is built rather than after. The
  refusal a user sees changes: it names the term values `dec` and
  `vint1` rather than the spellings `dec(decision)` and `vint(upper)`,
  and writes `rt | dec(<column>) ~ ...` as the example. Four pinned
  expectations changed with it. **This retires the 0.2.0 note below
  that "one hand-rolled check remains, and is not `required_aterms`'s
  fault"** - core grew the seam, and the check is gone. `gddm()`'s
  CONDITION index is still checked by hand, and still cannot be
  declared: which slot carries it moves with the boundary's spelling,
  so the requirement is a disjunction of conjunctions and no
  declaration says that.

* `ddm_cdf_ks`, the half-width of the no-go distribution function's
  image sum, is 4 rather than 12. The blend gives the small-time route a
  non-zero weight only below `u = 0.197`, where `tanh` has not yet
  saturated, and at that `u` the `|j| = 2` term is already 9e-23.
  Measured over 840 rows spanning `t` in 0.05 to 15, every truncation
  from 2 to 12 gives bit-identical values AND bit-identical gradients.
  It is a pure cost change and it is worth 2.3x on a go/no-go
  variability fit, which evaluates that function on a three-dimensional
  node grid. `ddm_cdf_kl` is not reducible the same way and is
  unchanged.

* `wiener_gng_simulate()` gains `sv`, `sz` and `st`, and `simulate()` on
  a fitted go/no-go model follows the variability, by drawing each
  trial's parameters before it runs the process. The rejection is on the
  OUTCOME, so it reweights all three at once where `wiener()` has to
  reweight the drift by the boundary probability it implies.

* Every family declares `frmtmb_family(accepts_aterms = )`, the core's
  new addition-term allow-list, and the hand-written `dec()` refusal
  the two race families shared is deleted. `lba()` and `rdm()` take
  `vint()` and `weights()`, `wiener()` adds `dec()`, `wiener_gng()`
  takes `dec()` and `vreal()`, and `gddm()` takes all four. A term
  outside a family's list is refused by name at frame assembly.

* BEHAVIOR CHANGE: `rdm()` refuses `vreal()`. The compatibility table
  recorded that pair as working on the ground that a model supplying
  one fitted and gave the same answer as one that did not, which is
  the defect rather than the feature: the column travelled into the
  fit and changed nothing. The row now reads refused.

* The `dec()` refusal for `lba()` and `rdm()` is now the core's
  generic sentence, which names the term and lists what the family
  does take, in place of the family-specific one. Same refusal, same
  families, one seam instead of a check per package.

# frmtmb.eam 0.3.0

* RENAMED from frmtmb.ddm. The package holds the linear ballistic
  accumulator, the racing diffusion model and the go/no-go diffusion
  beside the Wiener and generalized drift-diffusion families, so its name
  now says evidence accumulation models. Family constructors and the
  `ddm_` helper names are unchanged; `library(frmtmb.eam)` replaces
  `library(frmtmb.ddm)`, and the site is at frmtmb/frmtmb.eam/ with the
  old address redirecting. EMC2 leaves Suggests: no shipped code uses it.

* New `tests/testthat/test-gddm-reference.R`: `gddm()`'s generalized
  components are now checked against PyDDM, the reference
  implementation of Shinn, Lam and Murray (2020), rather than only
  against this package. Nothing outside the solver checked them
  before: the analytic Wiener comparison covers a constant drift and
  fixed bounds and no more, the gradient tests show the derivative
  matches the value it differentiates, and recovery is circular
  because `gddm_simulate()` draws from the solver's own density on
  purpose. Ten cases, covering leak at both signs, exponential and
  linear collapse, the coherence nonlinearity at two coherences,
  start-point variability and a lapse, are frozen into a fixture under
  `tests/testthat/fixtures/`; no Python runs at test time, and
  `dev/gddm-pyddm-reference.py` in the source repository regenerates
  it. Each case compares the density at both boundaries, the two
  boundary masses and the log-likelihood of a small fixed dataset,
  and the test also requires
  the disagreement to shrink as the grid is refined, so the reference
  is what the solver converges to and not merely something it lands
  near. One case is also checked against an Euler-Maruyama simulation
  of the equation itself, so that two grid solvers cannot be wrong
  together. The tolerances are measured, not guessed, and
  `vignette("gddm")` records what the reference covers, what it does
  not, and the one case where the two disagree.

* `frm_simulate()` now works for `gddm()`. The family installs its
  density and its simulator in `family_finalize()`, so before frame
  assembly the family object carries neither, and frmtmb's
  `frm_simulate()` read the family as written and refused the model
  for having no simulator. `simulate()` on a fitted `gddm()` was
  unaffected, because a fit carries the finalized family. The fix is
  in frmtmb rather than here: nothing in this package's wiring was
  wrong, and `wiener()` and `lba()`, which pass `sim` in the
  constructor, were never affected.
* New `tests/testthat/test-simulate-density.R`: the agreement tier for
  `wiener()`, `gddm()` and `lba()`, matching the one frmtmb holds its
  built-in families to. The other tests here check the density against
  RWiener and the exported `ddm_simulate()` against its generative
  process; neither of those is the seam frmtmb uses, which is the
  family's `sim` slot. The check is shaped for a choice-RT model: the
  density is defective, so the reference is that density renormalized
  on the row's own boundary, and the two boundaries' masses are
  asserted to sum to one.

## Two more families

* New `rdm()`, the racing diffusion model of Tillman, Van Zandt and
  Logan (2020): `n` independent Wiener accumulators with their own
  positive drifts, a start point uniform on `(0, A)` and a common
  threshold `A + k`, so each finishing time is an inverse Gaussian
  averaged over the start point. It is `lba()`'s geometry with the
  ballistic assumption removed, and it takes `lba()`'s spelling for
  that reason: a model moves between the two by changing one word. The
  drifts have a LOG link rather than an identity one, which is the one
  place the two part company, because a racing-diffusion drift is the
  rate itself and an accumulator with a rate of zero never finishes.
  `rdm_simulate()` draws from the generative process.
* New `wiener_gng()`, the go/no-go diffusion of Gomez, Ratcliff and
  Perea (2007), which EMC2 calls `DDMGNG`: a two-boundary diffusion
  where only the upper boundary produces an observable response. A go
  trial contributes the ordinary upper-boundary density; a no-go trial
  contributes the probability of no upper crossing before a deadline.
  The deadline goes on the family when every trial shares it and
  through `vreal()` when it does not, and `dec()` says which trials
  produced a response. `wiener_gng_simulate()` draws from the
  generative process. Across-trial variability is deliberately not
  offered; `?wiener_gng` says why.
* The Wiener defective distribution function, which this package did
  not have. `wiener()` declares no `lcdf` and its compatibility table
  says so, but the go/no-go family's no-go branch needs one, so it is
  written in `R/wiener-cdf.R` as two series blended in `log(u)` the way
  the density's two are. The blend is centred at `u = 0.02` rather than
  the density's 0.35, and that is the substantive choice: the
  small-time route reaches the no-go probability as `1 - F_upper`,
  which cancels exactly where a no-go trial is surprising, and the
  large-time route computes it directly and never subtracts. Handing
  over as early as the large-time route is accurate holds 2.2e-12
  relative against a 260-bit reference over 1200 points, with one row
  worse than 1e-12 and none worse than 1e-9.

## Measured against EMC2

* Both families agree with EMC2 where EMC2 has digits, and the
  comparison lives in `dev/rdm-gng-emc2-reference.R` rather than in the
  suite, because every EMC2 function that computes either likelihood is
  internal. Composed exactly as `EMC2:::log_likelihood_ddmgng` does,
  the go/no-go log likelihood agrees to 6.1e-16 per go row and 1.3e-15
  per no-go row. Composed as `EMC2:::log_likelihood_race` does, the
  racing-diffusion log likelihood agrees to 1.1e-12 over 200 rows at
  two, three and four accumulators.
* Where they disagree, the disagreement is adjudicated rather than
  asserted. EMC2 writes a race loser's survival as `1 - pWald(...)`,
  which returns EXACTLY ZERO on 35 of 315 grid rows, the first where
  the true survival is near 1e-13; `rdm()` writes it out directly and
  returns no zeros anywhere. Against `statmod`, which is party to
  neither, the single-accumulator density is 4.6e-13 from the truth for
  this package and 2.1e-04 for EMC2 on the rows where they differ.
  Similarly `WienR`, which is EMC2's own distribution function, is 4.4
  percent wrong at a no-go probability of 1.19e-13 at its default
  precision and at every setting down to 1e-12, and 0.033 percent
  wrong at its tightest; EMC2 calls it with `precision = 0.005`,
  looser than any of those, so EMC2 sits at the 4.4 percent end. This
  family's large-time route is 3.0e-15 there.
* The suite itself uses only EXPORTED references: `statmod`'s inverse
  Gaussian averaged over the start point for `rdm()`, and `WienR` and
  `RWiener` for the go/no-go distribution function. New in `Suggests`:
  `EMC2`, `statmod` and `WienR`.

## Fixed

* `fitted()`, `predict(type = "response")` and
  `residuals(type = "response")` now refuse on both new families
  instead of returning a drift rate. A family that declares no
  `post$mean_fn` gets frmtmb's fallback, "the first primary dpar on the
  response scale is the mean", and for these two that number is a
  drift: `rdm(3)` returned 3.51 for data whose response times average
  0.36, and `wiener_gng()` returned a constant 1.05 for data whose go
  response times average 0.6, with nothing to signal it. Both now
  declare a mean that stops with a reason. For `wiener_gng()` the
  refusal is also the right answer on the merits: a go/no-go trial
  produces a pair, and the no-go rows have no response time to average.
* `rdm()` refuses `dec()` by name rather than ignoring it. Measured on
  `lba()`, the sibling race family, a `dec()` term supplied alongside
  `vint()` is silently dropped and the model fits; a model ported over
  from `wiener()` would otherwise fit while quietly meaning something
  else.

## Also

* `tests/testthat/test-simulate-density.R` covers the two new families
  on the same terms as the other three: each family's `sim` slot
  against that family's own log density, and the defective masses
  asserted to sum to one. For `wiener_gng()` the go branch's support
  ends at the deadline rather than running to infinity, and a no-go
  row's draw is its deadline.
* The compatibility rows for both families were RUN rather than
  reasoned about, and several first guesses were wrong: `REML` and
  `quadrature` both work where they had been written down as refused or
  untested, `mixture()` refuses `rdm()` outright because its components
  need a dpar called `mu`, and neither family makes its addition terms
  mandatory on newdata for a link-scale prediction.
* `vignette("ddm")` gains a section per family, each with a fitted
  example, a recovery check and the comparison against EMC2.

# frmtmb.eam 0.2.0

Three families where there was one: Ratcliff's full diffusion model
as an extension of `wiener()`, the generalized drift-diffusion model
`gddm()`, and the linear ballistic accumulator `lba()`, plus the
queued defects. The three share one floor idiom and one
compatibility table.

## The full diffusion model
* `wiener(variability = )` adds across-trial variability to the family
  rather than forking it. Naming any of `"sv"` (drift rate), `"sz"`
  (start point) and `"st"` (non-decision time) turns that one into an
  ordinary distributional parameter, with its own link and its own
  formula. `wiener()` with no arguments is the model it always was: on
  a grid of 11520 parameter combinations its log density is
  byte-for-byte identical to the previous release everywhere the
  density is defined, and its tests are unchanged. The one difference
  is below the non-decision time, where it now returns `-Inf` instead
  of `NaN`, which is the deliberate fix below.
* The likelihood is the analytic Wiener density averaged over those
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
* Node counts are the `nodes` argument and the defaults are measured:
  `sz` reaches machine precision at 7 nodes, `st` reaches 1e-9 at 21.
  They differ because the non-decision-time range is cut by the response
  time on a fast trial, and the integrand turns on sharply at the cut.
  Node positions and counts are fixed when the family object is built,
  because an automatic-differentiation tape cannot record a branch on a
  parameter; a parameter only rescales the interval they map onto.
* Variability parameters at zero reproduce the plain Wiener density to
  better than 1e-13 in the log density, which is floating-point rounding
  on a differently associated sum rather than a quadrature error.
* Recovery on simulated data, 8 replicates of 1500 trials with a normal
  drift rate and a uniform non-decision time: every parameter's Monte
  Carlo mean is within 2.2 of its own standard errors of the value it
  was generated from.
* Across-trial variability is NOT frmtmb's `quadrature = TRUE`. That
  marginalizes random effects by Gauss-Kronrod, is wired to the
  random-effect coefficient vector by name, and refuses a model with no
  random-effect block. This integral shares nothing between trials, has
  no level to estimate, and exists in models with no grouping factor at
  all, so it lives inside the density.
* `ddm_simulate()` takes `sv`, `sz` and `st`, drawing each trial's own
  parameters and then running the ordinary process, so the simulator
  states the model independently of the density.
* `fitted()` and `simulate()` follow the variability rather than
  ignoring it. Conditioning on the boundary a row ended at reweights
  which per-trial parameters that row could have had, so the fitted mean
  is a ratio of two quadratures and the simulator accepts a drawn drift
  rate and start point with the boundary probability it implies. The
  plain closed form would not do: at an unbiased start point it returns
  the same mean for both boundaries, where 40000 simulated trials put
  them 0.06 s apart. Both are checked against simulated data, which
  knows nothing about how either is computed.


## The generalized drift-diffusion model

* `gddm()` is the generalized drift-diffusion family of Shinn, Lam and
  Murray (2020): a drift that may depend on the accumulator's own level
  and on a covariate, boundaries that may collapse within a trial, and a
  starting distribution that may be a point or an interval. There is no
  closed-form first-passage density, so every likelihood evaluation
  solves the Fokker-Planck equation forward in time and reads the
  probability flux through each boundary. It needs no change to core
  frmtmb.
* The components are chosen by argument and are extensible.
  `gddm_drift_constant()`, `gddm_drift_coherence()` and
  `gddm_drift_leak()` are summed to make a drift; `gddm_bound_constant()`,
  `gddm_bound_exponential()` and `gddm_bound_linear()` give the
  boundary; `gddm_start_point()` and `gddm_start_uniform()` give the
  start. `gddm_drift_term()`, `gddm_bound_term()` and
  `gddm_start_term()` are the documented seams for writing more. Every
  free quantity is a dpar that takes a formula, as in `wiener()`.
* `bs` is the boundary SEPARATION and `bias` the relative start point,
  both as in `wiener()`, so estimates are directly comparable between
  the analytic family and the generalized one.
* The substitution `y = x / B(t)` pins the moving boundaries at plus and
  minus one, which keeps the grid fixed while the boundary collapses and
  is what makes the likelihood differentiable: nothing on the taped path
  branches on a parameter. With the walls stationary the scheme is
  Crank-Nicolson, which a solver that chases a moving bound cannot use.
* The likelihood is an ordinary rowwise family, not a
  `frmtmb_structure()`. frmtmb calls `lpdf` once per objective
  evaluation with full-length vectors, and the condition a trial belongs
  to is data, so a rowwise density does one solve per condition, which
  is all a structure would have bought, and it keeps `fitted()`,
  `predict()`, `simulate()` and `residuals()` rather than defaulting
  them to refused. Measured: 6 solves for 2400 trials over 6 conditions.
* `gddm_control()` carries the grid and what is done with the answer:
  `dt`, `ny`, `t_max`, `max_ndt`, `renormalize` and `tridiagonal`.
  Renormalizing the defective density is on by default and should stay
  on: the discretized solve loses mass in a parameter-dependent way, so
  a likelihood that does not divide it out rewards fast absorption. On
  data simulated from the model, turning it off more than doubles the
  fitted leak and shrinks the boundary separation by a fifth.
* `tridiagonal` picks how the solve inside each step reaches the tape.
  `"recorded"`, the default, builds a large tape that runs in compiled
  code; `"atomic"` collapses the solve into one node with a hand-written
  adjoint, building about twelve times faster and evaluating about
  twelve times slower. Both give the same derivative to machine
  precision.
* The published coherence nonlinearity has no derivative at zero
  coherence, which a motion design normally contains. The coherence is
  data, so the zero condition is resolved once when the tape is built
  and never reaches it; the gradient in the exponent is finite, and is
  exactly zero there, because the drift is zero whatever the exponent
  is. Signed coherences are supported for stimulus coding.
* The family admits exactly two responses, because one accumulator
  between two absorbing boundaries has two walls. A decision indicator
  with more than two levels is refused at frame assembly, naming how
  many levels the data has and pointing at `lba()`, which fits the
  racing accumulators that more than two alternatives need, rather than
  being folded into one of the two.
* The boundary is read from `dec()` when it is there and from `vint()`
  otherwise, so the spelling brms uses works on this family as it does
  on `wiener()`. `vint()` numbers its values positionally, so the
  condition index is the first `vint()` value alongside `dec()` and the
  second inside `vint(upper, cond)`; the family reads whichever it is.
  Neither can be declared through `required_aterms`, which names the
  terms a density needs ALL of, so both refusals are written out.
* `gddm_conditions()` builds the condition index. `gddm_simulate()`
  draws choices and response times from the model's own solved density.
* Validated in the package's own suite: against this package's Wiener
  density with a constant drift and fixed bounds, to better than 0.01 in
  the log density at the shipped grid over decision times from 0.2 s on,
  degrading on a coarser grid and improving on a finer one; automatic
  gradients against numDeriv at parameter points including a collapsing
  boundary; parameter recovery with Monte Carlo standard errors; the
  size of the renormalization bias; and the zero-coherence gradient.
* The density is floored before it is logged, as the Wiener density in
  this package already was. Where the solved density at a trial's own
  response time underflows, the log density is a large finite negative
  number instead of `NaN`: the optimizer gets a value it can use, and a
  mixture's log-sum-exp is not poisoned by one component. The floored
  row is flat, so its gradient is exactly zero; `-Inf` would not do,
  because `-Inf` differentiates to `NaN`.
* `gddm_floored(fit)` is where that goes to be read. It returns the
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
* Known limits, measured and stated rather than hidden: the density at
  decision times of only a few time steps is far larger than the truth,
  because an implicit scheme spreads a little mass everywhere at once
  where the true density is exponentially small. A lapse component
  (`gddm(lapse = "uniform")`) floors it in the model rather than in the
  arithmetic. Cost scales with the number of conditions, so a design
  with many distinct parameter settings is where this becomes painful.
* `vignette("gddm")` fits one of the paper's models end to end and says
  plainly where this is slower than `wiener()` and why someone would pay
  that.
* That vignette now carries four figures, under a `tinyplot` gate: the
  three shipped boundaries over a trial, the fit against the data it was
  fitted to as defective cumulative distributions, the solved density
  against the analytic Wiener density with the log-density error on two
  grids beside it, and the leading edge of the fitted density at three
  time steps, which is what `gddm_floored()` counts. They are drawn from
  the vignette's own simulated data and its one fit, and together they
  cost about two seconds of a knit that the fit dominates.

## The linear ballistic accumulator

* `lba(n)` is a race between `n` accumulators, each rising in a straight
  line from a uniform start point on `(0, A)` to a common threshold at a
  normally distributed rate. It is the family for choices with more than
  two alternatives: a diffusion between two absorbing boundaries admits
  exactly two responses, and no reparameterization of `wiener()` reaches
  a third. The likelihood is closed form for any `n`.
* Each accumulator's drift mean is its own distributional parameter
  (`v1`, ..., `vn`, identity link) and all of them are primary, so the
  main formula reaches every drift and each gets its own coefficients.
  Give one its own formula and a covariate moves that alternative alone,
  which is the capability neither two-choice family can offer.
* The remaining parameters are `A` (start-point range, log), `k`
  (threshold above the start-point range, log) and `ndt` (non-decision
  time, a logit scaled onto `(0, max_ndt)` as in `wiener()`). The
  threshold is `A + k` rather than a parameter of its own, so `b > A` is
  structural: a threshold inside the start-point range, where trials
  would begin already finished, is not a state the optimizer can reach.
* The drift standard deviation is fixed, because the model is identified
  only up to a common rescaling of `A`, the threshold, the drifts and
  the drift standard deviation. It is the family argument `sd_v` rather
  than a hidden default, and is carried on the family object.
* Drift rates are truncated at zero by default, following `rtdists` and
  its `posdrift = TRUE`. Every accumulator then arrives eventually and
  the choice probabilities sum to one. `posdrift = FALSE` gives the
  untruncated convention, under which the response distribution is
  defective; the two are different models, not a rescaling of each
  other, which matters when comparing against another package.
* The single-accumulator density agrees with `rtdists::dlba_norm()` to
  better than 1e-11 relative wherever `rtdists` is itself accurate, and
  the race with `rtdists::n1PDF()` for two, three and four
  accumulators. In the fast tail the two diverge by design and this one
  is the better: `rtdists` writes the normal difference
  `Phi(g) - Phi(h)` as a subtraction of two lower tails, which returns
  exactly zero once `pnorm(h)` saturates. Written in log space instead,
  as `exp(la) * -expm1(lb - la)` on the two upper-tail logs, nothing
  saturates and no comparison is needed, which matters because RTMB
  refuses comparison on AD types. Against a 200-bit Rmpfr reference on
  1144 points, the subtractive form returns exactly zero on 68 of the
  80 tail rows and is already 2.7e-3 wrong on 8 rows outside the tail;
  the log-space form holds 5.0e-14 in the bulk and 3.6e-4 in the tail.
  Four rows still exceed 1e-6, so it is an improvement, not a proof.
* That mattered beyond accuracy. The subtractive form errs in the value
  only, and a tape differentiates the function that was written, so
  value and gradient described different surfaces near the
  non-decision-time bound. It did not move the point estimate, which
  the likelihood keeps out of that region, but Hessian-based standard
  errors were up to 16.5 percent off before the change.
* The survival function is written out directly rather than as
  `1 - plba_norm()`, which returns exactly zero, and so a
  log-contribution of `-Inf`, while the true survival is still around
  1e-19; the direct form agrees with a quadrature of the density down
  to survivals of 1e-23. Under `posdrift = FALSE` the never-arriving
  mass it adds back is spelled `Phi(-v/s)` rather than `1 - Phi(v/s)`,
  which is the same saturation trap one line further on.
* The choice reaches the density through `vint()`, as the Wiener
  family's boundary indicator does, and is declared with
  `required_aterms`, so omitting it is refused by name. A choice outside
  `1..n`, a non-positive response time and a non-decision-time bound
  above the fastest response are each refused with their own message.
* The non-decision-time bound is derived from the response through
  `family_finalize()`, so the family object carries the link it will
  actually use rather than an environment filled in later.
* `lba_simulate()` draws choice and time jointly from the generative
  process. `simulate()` on a fitted object redraws times conditional on
  each row's observed choice, since the choice is data.
* Deliberate omissions: no `lcdf`, so `cens()` and `trunc()` are
  refused; no `post$mean_fn`, because the mean of the race has no closed
  form and a quadrature per row was not worth writing for `fitted()`.

## Defects fixed

* **The density returns `-Inf` below the non-decision time, not `NaN`.**
  A likelihood of zero is the right answer there; `NaN` is not an answer
  and propagates through every component of a mixture. The normalized
  time is held at the smallest positive double, which is bit-for-bit
  inert everywhere the density was already right.
* **`wiener(max_ndt = )` above the smallest response time is usable in a
  mixture.** The refusal is correct for the family on its own and wrong
  inside a mixture, where the contaminant component is exactly what
  covers the trials the diffusion cannot produce, so
  `allow_unreachable = TRUE` lifts it and the refusal now names it. A
  trial below the non-decision time then gets a log density that
  exponentiates to zero AND differentiates to zero, which a true `-Inf`
  does not: `-Inf` yields a `NaN` gradient and stops the fit.

## What frmtmb 0.49.0 let this package delete

* **`dec()` is the spelling now.** frmtmb gained
  `frmtmb_register_aterm()`, this package registers `dec` when it loads,
  and `rt | dec(response) ~ x` works and takes a factor, a character
  vector or a logical the way brms does. `vint()` carries the same thing
  as a 0/1 integer and is unchanged.
* **The environment the link closures read is gone.** The bound on the
  non-decision time is a property of the response, and the family object
  is built before `frm()` has any; this package used to have `valid_y()`
  write the bound into an environment, which worked only for as long as
  an undocumented slot order held. `family_finalize()` is the documented
  slot for it and the family now derives itself from the data there.
* One hand-rolled check remains, and is not `required_aterms`'s fault:
  that argument names the terms a density needs ALL of, and this family
  needs EITHER `dec()` or `vint()`. See `dev-findings.md`.

# frmtmb.eam 0.1.0

First release. A Wiener first-passage time family for two-choice
response times, written entirely against frmtmb's exported extension
API.

* `wiener()` is the drift-diffusion family, with brms's dpar names and
  links: `mu` (drift rate, identity), `bs` (boundary separation, log),
  `ndt` (non-decision time, bounded) and `bias` (relative start point,
  logit). A model written for `brms::wiener()` reads the same here, and
  the parameterization is pinned against brms in the test suite.
* The density is the Navarro and Fuss (2009) pair of series. Both are
  evaluated at a fixed truncation and their logs blended with a
  logistic weight in the normalized time, because an AD tape cannot
  choose between them on a parameter. It agrees with
  `RWiener::dwiener()` to better than 1e-11 relative over normalized
  times from below 1e-3 to 50, where a fixed truncation of the
  small-time series alone is wrong by tens of percent past about 8.
* The decision indicator reaches the density through `vint()`, because
  frmtmb's addition terms are a closed set and brms's `dec()` is not one
  of them. Omitting it is refused with a message naming both spellings,
  rather than silently producing a log likelihood over no rows.
* `ndt` uses a logit scaled onto `(0, max_ndt)`, so the support
  constraint is structural rather than something the optimizer has to
  discover. `max_ndt` defaults to the smallest response time in the
  data, found at frame assembly.
* `post$mean_fn` gives the mean response time in closed form,
  conditional on the boundary the row ended at, with the zero-drift
  limits handled explicitly. `sim` draws conditionally by inverse
  transform through `RWiener::qwiener()`, with a discretized forward
  simulation as the fallback when RWiener is absent.
* `ddm_simulate()` draws response times and boundary choices jointly,
  for building example and test data sets.
* The package registers its compatibility rows with the core at load
  through `frmtmb::frmtmb_register_compat()`, so `frm_compat("wiener")`
  states what was exercised and what was not.
* Deliberate omissions: no `lcdf`, so `cens()` and `trunc()` are
  refused; no variance function and no unit deviance, so
  `residuals()` answers `type = "response"` only.
* `dev-findings.md` in the package source records what building this
  from outside frmtmb cost, as an acceptance test of the extension API.

