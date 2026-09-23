# skew_normal() stops at alpha = 0, and what fixed it

Lane `wt-skewinit`, 2026-09-22, off main 40a1611 at frmtmb 0.61.0.
Base for every comparison is the read-only reference library
`C:/Users/adf44/source/r/rellib-r3` (frmtmb 0.61.0). The lane build is
`C:/Users/adf44/source/r/skewinit-lib`. Every script named below is in
`dev/`, every log in `dev/skewinit-log/`.

## The answer

`alpha = 0` is a stationary point of the skew-normal likelihood at
every sample, and the information is singular there, so a fit that
reaches it converges with code 0 and a log-likelihood up to 35 units
below the maximum. Two changes remove it.

1. **`alpha` and `sigma` now start from the RESIDUAL.** An
   `init_dpars` function that declares a third argument is handed the
   response with its `mu` predictor removed by least squares, and
   `skew_normal()` declares one. This alone takes the two 40-seed
   acceptance streams from 28 and 33 stalls to 0 and 0.
2. **A fit that lands on the point anyway is refit from +2 and -2.**
   Declared by the family as `post$stationary`, executed by
   `escape_stationary()` in `R/fit.R`. It alone, with the OLD start,
   also fixes 80 of 80.

There is no `R/estimate.R` in this tree; estimation lives in `R/fit.R`.
No second optimizer driver was written: `escape_stationary()` calls
`optimize_obj()`, the single entry point every fit mode already uses,
through its existing `start_par` argument, which is the same argument
`fit_recovery_starts()` uses to restart a failed fit. `optimize_obj()`
already applies bounds, autoscale units, `control$restarts` and the
restart-from-best recovery, so the escape inherits all four.

Both are kept. The start does the work and the refit is the guard,
because the start's success is not robust: it depends on `sigma`'s
start as well as `alpha`'s (see "Why both", below).

Nothing is printed. The reasoning is in "What the user is told".

## Reproduction on the base build

`dev/skewinit-repro.R`, log `skewinit-log/repro-base-stream-stated.txt`.
Seed 1 of the stream as `dev/drmtmb-skewnormal.R` now writes it
(n = 200, `xs <- -abs(rnorm(n)) * 3`, `y <- xs + (abs(rnorm(n)) -
sqrt(2/pi)) * 1.5 + rnorm(n, 0, 0.3)`):

| | value |
|---|---|
| skewness of the raw response | -0.424213543522 |
| skewness of `lm(y ~ xs)` residuals | 0.913191947002 |
| default logLik | -284.026142739 |
| default alpha | -0.000836862957676 |
| convergence code, message | 0, "relative convergence (4)" |
| max abs gradient at the stall | 1.40473334642e-05 |
| conditions raised by the fit | 0 |
| standard error of alpha | 5.46 |
| logLik from `start = list(betad = c(0, 2))` | -259.338935404 |
| gap | 24.6872073344 |
| `sn::selm(y ~ xs, family = "SN")` logLik | -259.338935404 |
| `sn::selm` alpha | 7.50645961055 |
| sn against the started fit | 1.855e-10 |

The two skewnesses have opposite signs, which is the whole defect: the
start rule read the response's.

Over both streams on the base build (`dev/skewinit-gaprange.R`, log
`skewinit-log/gaprange.txt`), reproducing what `wt-drmtmb` recorded to
the digit:

| stream | stalls | gaps below `sn::selm` |
|---|---|---|
| stated | 28 of 40 | 6.406930 to 35.137329 |
| dead draw | 33 of 40 | 6.107917 to 32.692143 |

`sn` 2.1.3 came from the user library; nothing was installed for it.

## The calibration of the detection rule

`dev/skewinit-calibrate.R` fits 240 models on the BASE build across six
arms of 40 seeds: the two acceptance streams, alpha truly 0 at n = 200
and n = 50, and alpha 1 and 0.5 at n = 200. For each it records the
fitted alpha, its standard error, the objective at alpha displaced by
0.25, 0.5, 1 and 2, `sn::selm`'s optimum, and the best of a +2/-2
refit. Log: `skewinit-log/calib.csv`, summaries
`calib-summary-1.txt` and `threshold.txt`.

### What does NOT work

**The standard error does not separate.** At a stall it is 0.48 to
35.7; on genuinely symmetric data 0.49 to 48.8. That is expected rather
than unlucky: the information is singular at `alpha = 0` whether or not
the sample is skewed, so a large standard error is what BOTH look like.

**The objective probe does not separate either.** Displacing alpha by
0.25 from the optimum and keeping every other parameter changes the
objective by 0.0024 to 0.1487 on rows a refit improves and by up to
0.0026 on rows it does not. The ranges overlap, so the probe cannot
gate the refit. Recorded because it is the obvious idea.

### What does work

**Where the fit landed.** Across THESE 240 fits, `|alpha|` never fell
in the open interval (0.006442, 0.091507): the largest `|alpha|` at a
stationary point was 0.006305, and the smallest genuine optimum was
0.006442 (itself a case where `sn` agrees the optimum is there).

**That band is a property of these six arms and NOT of
`skew_normal()`.** The reviewer of punch round 1 ran 270 fits over nine
designs this lane did not use (n from 50 to 5000, covariate scales
1e-3 to 1e3, true alpha 0, 0.1 and 0.2, t and uniform residuals) and
found the fitted `|alpha|` inside the band on 10 of them. So the band
is not empty in general and the threshold cannot be justified by it
alone. The reviewer could not turn any of the 10 into a wrong answer,
and the reason is structural: a fit inside the window is refitted from
+2 and -2 and the BETTER optimum is kept, so being inside the window
costs two optimizer runs and can never lower the likelihood. The
consequence of the band being non-empty is cost, not correctness.

**The rule, as implemented.** Refit from +2 and from -2 when the fitted
`alpha` PREDICTOR is within **0.05** of 0 at every row, on the link
scale. The justification that does carry beyond these arms is the
interpretive one: a skew normal with `alpha = 0.05` has skewness
2.7e-5, so a fit inside the window is reporting a distribution that is
numerically symmetric, and checking that claim from both sides costs
two fits. What a fit inside the band costs is exactly that: two extra
optimizer runs, about 1.3x the wall clock of the fit (measured below),
and no risk to the answer.

It evaluates the objective zero times when it does not fire, by
construction: `stationary_escapes()` touches `opt$par`, the design and
the family, and nothing else.

### False negatives

Measured by forcing the OLD start on both acceptance streams, which is
what the escape has to survive on its own
(`skewinit-log/accept-escape-only.txt`, 80 fits):

| stream | stalls | detected | rescued to `sn`'s optimum |
|---|---|---|---|
| stated | 28 of 40 | 28 | 28 |
| dead draw | 33 of 40 | 33 | 33 |

0 of 61 missed. Worst remaining shortfall against `sn` over the 80:
7.474e-10.

Over the whole 240-fit calibration set, the rule misses 1 of the 157
fits that are not at the optimum. That one is `alpha = -0.551280`, a
second local optimum 0.0238 below `sn`'s, not the stationary point. It
is filed in `dev/test-backlog.md`; a window wide enough to catch it
would refit on genuine optima.

### False positives

A false positive costs two optimizer runs and nothing else: the refit
keeps the better optimum, so it can never make a fit worse.

On the FIXED build, where the residual start is in place
(`dev/skewinit-falsealarm.R`, `skewinit-log/falsealarm-fixed.txt`,
40 seeds per arm):

| arm | fires | gains > 0.01 | gains > 1 | fires and gains nothing |
|---|---|---|---|---|
| alpha = 0, n = 200 | 1 | 1 | 0 | 0 |
| alpha = 0, n = 50 | 2 | 1 | 0 | 1 |
| alpha = 1, n = 200 | 2 | 2 | 1 | 0 |
| alpha = 0.5, n = 200 | 5 | 4 | 1 | 1 |
| half-normal residual | 0 | 0 | 0 | 0 |

10 of 200 fits fire; 2 of the 200 fire and gain nothing. On genuinely
symmetric data the rate is 3 of 80 firing and 1 of 80 wasted.

On the BASE start, which is the worst case for the detector, the
symmetric arms fired on 20 of 40 (n = 200) and 18 of 40 (n = 50), and
6 of those 38 gained nothing. So the residual start is what keeps the
refit rare; the threshold is unchanged between the two measurements.

One thing worth saying plainly, because it surprises: on genuinely
symmetric data the refit often does find a better optimum, at a large
`|alpha|` with a large standard error. That is the skew normal's known
behavior, not a defect of the refit. The old code reported `alpha = 0`
there, and `alpha = 0` was not the maximum.

## Acceptance

`dev/skewinit-accept.R`, log `skewinit-log/accept-fixed.txt`. Both
40-seed streams, 80 fits, against `sn::selm` on the same data. A fit
counts as reaching the optimum when its logLik is within 1e-6 of
`sn`'s, an absolute tolerance chosen because 0.01 logLik units is the
scale at which a likelihood difference means anything.

| stream | fits | reached `sn` | worst shortfall | escape fired |
|---|---|---|---|---|
| stated | 40 | 40 | 1.884e-10 | 0 |
| dead draw | 40 | 40 | 1.037e-09 | 0 |

"stated" is the stream `dev/drmtmb-skewnormal.R` now writes; "dead
draw" is the first run's, with one extra `rnorm(n)` before the
covariate. Both are acceptance.

**Worst remaining gap over all 80: 1.03653974293e-09**, relative to the
log-likelihood 3.92e-12. No fit exceeded `sn` by more than 1.14e-13.

The escape fired on 0 of the 80, because the residual start already
avoids the stationary point on all of them.

## Why both pieces are kept

`dev/skewinit-starts.R` measures four start rules on the BASE build,
each passing BOTH `betad` entries explicitly, so the `sigma` start
behind every number is stated. 40 seeds per stream, 80 fits per arm.
Log: `skewinit-log/starts-base.txt`.

| arm | alpha from | sigma from | stalls, stated | stalls, dead draw |
|---|---|---|---|---|
| `raw_sdy` (= the 0.61.0 default) | raw skew | log sd(y) | 28 | 33 |
| `res_sdy` | residual skew | log sd(y) | 0 | 1 (seed 34) |
| `res_sdr` (= what ships) | residual skew | log sd(resid) | 0 | 0 |
| `raw_sdr` | raw skew | log sd(resid) | 31 | 28 |

Three things follow.

`res_sdy`'s one failure is seed 34 of the dead-draw stream, at a gap of
25.307078. Its residual skew is 1.1048193 and the start rule gives
2.5524097, the correct sign and a long way from zero, reproducing the
wt-drmtmb reviewer's 2.5524 (`dev/skewinit-claims.R`). A start that is
right about alpha can still fall into the point, which is why the
refit stays.

The sigma start was wrong on its own terms. `sigma` in this
parameterization is the CONDITIONAL standard deviation
(`RTMBdist::dskewnorm2` takes mean, sd, alpha), so `sd(y)` overstates
it by whatever the mean structure explains. Seed 1: `log(sd(y))` is
0.662949 and the fitted `sigma_(Intercept)` is -0.014760, while
`log(sd(resid))` is 0.003698.

Fixing sigma alone makes things WORSE (`raw_sdr`, 31 stalls against
28). The two starts interact, so neither is a fix by itself.

Every other family still starts `sigma` from `sd(y)`. That is filed in
`dev/test-backlog.md` rather than fixed here: those families' optimizers
recover without help, and changing it moves the path of nearly every
fit in the package.

## What the user is told: nothing

Silence, plus a record on the fit and a `verbose = TRUE` line. The
argument, against the alternatives:

- **A message.** It would fire on roughly one skew-normal fit in twenty
  under the shipped start, and on about half of all symmetric-data fits
  if the start ever regresses. Nothing else in frmtmb narrates an
  optimizer restart: `control$restarts` is silent, and
  `optimizer_from_best()` speaks only under `verbose`. A message here
  would be telling the user that the optimizer did its job.
- **A `diagnose()` line.** Measured rather than assumed
  (`dev/skewinit-claims.R`, log `skewinit-log/claims.txt`). On a
  symmetric-data fit where `alpha` is weakly identified but finite
  (seed 3, `alpha` 0.2477, standard error 5.09) `diagnose()` reports
  convergence 0, `pdHess = TRUE` and a smallest covariance eigenvalue
  of 0.0045: it says nothing, and the only signal the user gets is the
  standard error itself. When `alpha` runs away instead (seed 6,
  `alpha` 2.7e6) `diagnose()` DOES speak: convergence 1, "false
  convergence (8)", `pdHess = FALSE`, smallest eigenvalue -0.378, and
  the fit warns as well.

  So the claim that `diagnose()` already covers this is FALSE in the
  finite case, and I am not adding a line anyway. The reason is that
  what would be reported is weak identification of `alpha`, which is a
  property of the skew normal that exists with or without this change
  and is not about the refit. Filed in `dev/test-backlog.md` as its own
  question.
- **What `verbose = TRUE` prints** (`dev/skewinit-verbose.R`, log
  `skewinit-log/verbose.txt`), on a fit forced onto the point:

      frmtmb: optimize [0.00s]: objective 284.02614
      frmtmb: stationary escape [0.00s]: 2 restarts, gain 24.7
      frmtmb: done [0.10s]: objective 259.33894, ...

  The same fit from the shipped start prints no escape line at all.
- **What is recorded.** `fit$opt$stationary_escape` is
  `c(starts =, gain =)` whenever the escape RAN, gain or no gain, so a
  fit that paid for two restarts and kept its own optimum is
  measurable. `verbose = TRUE` prints it as a stage beside the other
  optimizer stages. brms has no analogue, so there is no parity answer
  to defer to.

The NEWS entry states plainly that a `skew_normal()` fit with a `mu`
covariate can now report a different logLik, different estimates and
different standard errors, and that where it differs the new value is
the higher likelihood.

## What did not move

`dev/skewinit-nomove.R` fits 27 jobs on both builds and
`dev/skewinit-nomove-cmp.R` compares logLik, the whole optimizer
parameter vector and every standard error with `identical()`. Log:
`skewinit-log/nomove-compare.txt`.

**20 of 20 jobs that are not `skew_normal` are bitwise identical in all
three**, covering gaussian (ML, random effects, REML, `profile = TRUE`,
`sigma ~ x`), student, Gamma(log), lognormal, weibull,
inverse.gaussian(log), exgaussian, asym_laplace, poisson (ML, random
effects, `profile = TRUE`), negbinomial, bernoulli and Beta.

**Both `skew_normal` jobs without a `mu` covariate are bitwise
identical too**, including one with a random effect. That is by
construction: `mu_residuals()` returns the response untouched when the
`mu` design is an intercept alone, so the initializers compute the same
floating-point numbers they computed from `y` before. Centring the
response instead moved `sd()` in its last bit on 5 of 5000 random
samples (`dev/skewinit-sdulp.R`, log `skewinit-log/sdulp.txt`), and a
start that moves in its last bit moves every iterate after it.

**The seven `skew_normal` jobs WITH a `mu` covariate all moved**, which
is the change. Four of them have a finite alpha and the two builds
agree on the optimum:

| job | relative d logLik | max relative d par | max relative d se |
|---|---|---|---|
| `skew_fin` (ML) | 3.24e-15 | 2.18e-07 | 1.59e-07 |
| `skew_fin_reml` (REML = TRUE) | 7.61e-16 | 6.31e-08 | 4.75e-08 |
| `skew_fin_prof` (profile = TRUE) | 3.24e-15 | 4.52e-07 | 8.08e-08 |
| `skew_fin_sig` (`sigma ~ x`) | 9.87e-15 | 1.94e-06 | 1.08e-06 |

The other three (`skew_cov*`) were built with a half-normal residual,
so their MLE for alpha is unbounded: both builds stop on a flat ridge
at alpha between 3.7e5 and 5.2e7, and they differ by up to 1.24e-06
relative in logLik. `sn::selm` caps alpha at 183.4 and reports a logLik
1.48 BELOW both, so neither build is short of a finite optimum there;
the ridge is where the difference lives. Probed in
`dev/skewinit-prof.R`, log `skewinit-log/profile-probe.txt`.

## Punch round 1: what the review found and what it changed

Two blockers and four majors, all reproduced before anything was
touched (`dev/skewinit-punch1.R`, logs `punch1-before.txt` and
`punch1-after.txt`).

**B1, a dpar design with no intercept got NEITHER half of the fix.**
Both halves keyed on a column literally named `"(Intercept)"` and
skipped the design when there was none, so `alpha ~ 0 + g` started
every coefficient at zero, which IS the stationary point, and was
never escaped. Measured on this lane's own design, n = 300 seed 7:
`alpha ~ 0 + g` gave -394.1388237 with the fitted alphas at 9.5e-15,
convergence 0 and no warning, against -340.6074122 for `alpha ~ g` on
the same column space, a loss of **53.53 units**.

The fix reads and writes the PREDICTOR instead of an intercept
coefficient. `predictor_at()` returns the intercept answer exactly
when there is an intercept, so those models do not move at all, and
otherwise the least-squares solution of `X b = value`, which puts the
value in every cell of `~ 0 + g`. Detection likewise compares the
fitted predictor against `at` at every row. After the fix
`alpha ~ 0 + g` gives -340.6075687, the same optimum as `alpha ~ g`.

A design whose columns cannot represent a constant at all
(`alpha ~ 0 + x` with a centred x) has no least-squares answer but is
still stallable, so `predictor_away()` moves along the design's
largest column, scaled so the predictor's root mean square is the
value asked for. Checked that this design is not left short: from 13
starts between -20 and +20 it reaches -392.5188098 every time, which
is its own optimum, and the default reaches it too.

**B2, an atomic `post$stationary` killed the fit.** The dpar was
indexed out of the field before the field's own shape was checked, so
`post$stationary = "alpha"`, `c(0, 1)`, `c(alpha = 0)`, `TRUE` and
`NA` all raised "subscript out of bounds" and took the fit with them.
The six malformed cases already tested all wrapped the declaration in
a list, so none reached it. The shape check moved outside the `[[`,
and the five atomic shapes are now in the test.

**M1, the empty band.** Restated above as a property of these six
arms rather than of `skew_normal()`, with what a fit inside it costs.

**M2, `offset()` was not removed from the residual.** `offset()` is
part of the mu predictor and is not in `X`, so entering the acceptance
covariate as `offset(o)` left its skew in the residual and the start
failed exactly as the raw response did: the escape fired on 6 of 6
seeds with gains of 14.42 to 34.56. `mu_residuals()` now subtracts it,
and the escape fires on 0 of 6. Blast radius, measured: `mu_residuals()`
is reached only by an initializer that declares a third argument, and
`skew_normal()` is the only core family that does, so the six
non-skew_normal offset models in the blast battery are **bitwise
identical**; the two skew_normal offset models agree to 1.1e-13 and
2.3e-14 relative.

**M3, a new convergence warning.** Answered with numbers rather than
a change, in its own section below.

**M4, stale base numbers.** Re-measured twice, and the figure below is
the punch-round-2 one, on a clean process against the file as it
ships: `rellib-r3` gives **46 pass, 16 fail, 4 error**
(`skewinit-log/test-new-on-base.txt`), and the lane build **70 pass,
0 fail, 0 error** (`test-new-on-lane.txt`). The earlier
`test-new-on-fixed.txt` citation is dropped. The file grew in both
rounds, so these supersede round 1's 32/13/3 rather than contradict
it. The round-1 two-pass discrepancy against the reviewer (34/13/4
against this runner's 32/13/3) is **closed**: it was the reviewer's
runner not setting the working directory before `test_file()`. With
that aligned both runners measure base **46 / 16 / 4** and lane
**70 / 0 / 0**. The numbers agree; neither instrument was wrong about
the code.

### A regression the fix introduced, and how it was caught

`predictor_at()` was first written with an `is.matrix(X)` gate. A
`frmtmb_control(sparse_x = TRUE)` design is a **Matrix**, not a base
matrix, so that gate returned NULL and `make_start()` then placed no
initializer at all: every sparse model silently lost its starting
values. `tests/testthat/test-sparsex.R` caught it at **49 pass, 15
fail** in the tier rerun; nothing in this lane's own battery would
have. The helpers now judge a design by `dim()` and convert only where
`base::qr()` needs a dense matrix, and `test-sparsex.R` is back to
**64 pass, 0 fail**. A sparse case is pinned in the lane's own test
file so it is not caught only by someone else's.

The lesson is the tier's, not the fix's: a start change reaches every
family and every design class, and the only instrument wide enough to
see that is the whole suite.

### Punch round 2: the placement was not scoped to alpha (BLOCKER)

The round-1 fix replaced "skip a design with no intercept" with
`predictor_at()` for **every** linear predictor of **every** family,
not just for a dpar that declares a stationary point. On an
ill-conditioned design that silently loses the optimum: the reviewer's
`ypois ~ 0 + xt` with xt at 1e-6, n = 250, seed 23, went from base's
-5189.03446731 (matching `stats::glm()`) to **-7156.06063609**, a loss
of **1967.03 units**, at convergence 0 with no warning, because the
placement starts the coefficient at 1.6e6 and nlminb reports
X-convergence before the fit moves. Right on the predictor scale,
ruinous on the parameter scale.

**Fix**: the placement is now confined to a dpar whose family declares
`post$stationary`; every other dpar keeps the previous rule, an
intercept or nothing. The escape itself was already correctly scoped.

Proved with a new cross-family battery, `dev/skewinit-noint.R` and
`dev/skewinit-noint-cmp.R` (log `noint-compare.txt`), 19 no-intercept
jobs plus a four-point poisson scale sweep against `glm()`:

- **16 of 16 non-skew_normal jobs bitwise identical to base**:
  gaussian, gaussian `sigma ~ 0 + g`, student, student `nu ~ 0 + g`,
  Gamma, Gamma `shape ~ 0 + g`, exgaussian, lognormal, poisson,
  negbinomial, bernoulli, Beta, Beta `phi ~ 0 + g`, zero-inflated
  poisson and its `zi ~ 0 + g`, weibull.
- the poisson sweep is bitwise identical to base at **every** scale
  (1, 1e-2, 1e-4, 1e-6), and matches `glm()` exactly at 1, 1e-2 and
  1e-4.
- the three skew_normal jobs move, and `dev/skewinit-noint-why.R`
  shows WHY from the starts themselves: `sigma ~ 0 + g` is
  `0 0 0 0` and `mu ~ 0 + g` is `0 0 0 0` on **both** builds, so a
  non-declaring dpar is still skipped; only `alpha ~ 0 + g` is placed,
  at 2.31835 per cell against base's zeros.

The round-1 blast battery also improves from 7 to **15 of 20** bitwise
identical, and the one job that had been LOWER (`stu_sig_noint`,
-1.87e-08) is identical again.

**What I did NOT take.** Several other families' no-intercept starts
got better under the unscoped version (gaussian +263.96, Gamma +624.94,
exgaussian +16.85 on ordinary cell-means models). That is a real
finding and a separate change with its own risk, exactly as the poisson
case shows; it is filed in `dev/test-backlog.md` with the
counterexample and with what a safe version would need.

One claim of mine was WITHDRAWN in round 3 rather than quietly edited:
the backlog said "Autoscale does not rescue it" of the pre-existing
no-intercept poisson defect. It is false. `frmtmb_control(autoscale =
TRUE)` closes the 62.578284 gap to **0.000000** at every scale
including 1e-6 and 1e-8, verified here in `dev/skewinit-punch3.R`. The
entry now carries the table and says the open question is why the
DEFAULT does not engage autoscale there, not why the fit is wrong. A
confident wrong diagnosis in a backlog costs more than an open
question, because it tells the next reader to stop looking.

Round 3 also filed the residual hazard inside the declaring dpar:
`alpha ~ 0 + xs` at scale 1e-6 starts the coefficient at -6.769e+05
and stops 5.886 units short of a predictor-scale sweep with
"X-convergence (3)" and no warning, while still being 9.65 to 18.91
units BETTER than base on the same design, and better on 48 of 48
paired fits. It is the blocker's mechanism surviving inside the one
dpar that opted in, and a cell-means test cannot see it.

Two records corrected in round 2: `from = c(Inf, -2)` now drops
the non-finite entry AT THE GUARD, so it yields one start by decision
rather than by a downstream rejection, and the Rd says so. The base
test-file counts are re-measured below on a clean process.

### The blast radius of the two start changes

`dev/skewinit-blast.R`, 20 jobs, both builds, log
`blast-compare.txt`. Seven jobs are bitwise identical, including every
non-skew_normal `offset()` model and both no-offset intercept
controls. Of the rest:

- `sn_alpha_noint` is the B1 fix: base -339.552804089, lane
  -329.346656684, **+10.206147**.
- every other no-intercept job agrees to at most **4.7e-11** relative
  in log-likelihood.
- `stu_sig_noint` is the one job where the lane is lower, by 1.87e-08
  absolute and 4.7e-11 relative. Checked rather than waved past: the
  only parameter that moved is student's `nu`, 18.94687 to 18.49573,
  which has no finite standard error in that fit. It is a flat
  direction, not a worse optimum.
- **0 new warnings and 0 new nonzero convergence codes** across all 20.

## M3: the convergence codes, measured

Over the 80 acceptance fits the lane raises **0** nonzero codes, **0**
warnings and **0** messages, and so does base
(`dev/skewinit-convrate.R`). Acceptance now records the optimizer's
code, the condition counts and the gradient, not just the
log-likelihood.

Over the 200-fit spread the lane raises nonzero codes on 37 against
base's 7 (`dev/skewinit-convrate200.R`). That looked bad until the 30
extra ones were looked at (`dev/skewinit-convwhy.R`):

- **all 30 are runaway-alpha fits**, and on **all 30 the lane's
  log-likelihood is HIGHER than base's**, by a median of **38.62**
  units. Base looked clean there because it was sitting on the
  stationary point at `|alpha|` near 0.001.
- the association is exact on both builds: `|alpha| > 1e3` against a
  nonzero code is 37 of 37 and 163 of 163 on the lane, 7 of 7 and 193
  of 193 on base.

So the bulk of the increase is the optimizer correctly reporting an
unbounded parameter on fits base never reached. That leaves the
reviewer's own case, which is different: finite alpha 5.79,
log-likelihood equal to `sn::selm` to 2.3e-10, max abs gradient
4.81e-04 against `grad_tol` 1e-3, and a "false convergence (8)"
warning. Ruled out by measurement: autoscale (`par_units` is NULL, it
never engaged) and the escape (never fired); five extra restarts do
not clear it. The cause is which start the optimizer arrives from, and
the lane's start is the BETTER one here (sigma 0.392 against base's
2.388, true value 0.405; alpha +2.47 against base's -2.05, true value
4). I did not change the warning: suppressing a nonzero PORT code
whenever frmtmb's own gradient criterion is met would hide real
failures in every family and every fit mode, which is not this lane's
decision. Filed with the reproduction.

## The release tiers

### The run that completed, on the PRE-punch build

Counts generated by `dev/skewinit-tiers.R` from the logs, never typed.

<!-- generated by dev/skewinit-tiers.R, do not type these -->
| tier | files | pass | fail | error | skip | no RESULT line |
|---|---|---|---|---|---|---|
| ungated, 1 process per file | 270 | 15952 | 0 | 0 | 161 | 0 |
| gated, all env vars set | 38 | 3130 | 0 | 0 | 0 | 0 |
<!-- end generated -->

0 files with any failure or error in either tier, and 0 files produced
no RESULT line. The 15 ungated zero-count results are all
`test-brms-suite-*`, and all 15 appear in the gated log; checked, not
assumed. Four of those names occur twice because the same file exists
in core and in `frmtmb.sample`.

The ungated tier runs with four parallel workers
(`dev/skewinit-release/run-suite-worker.ps1`, a separate file; the
release scripts are unchanged), still one R process per file, because
the serial pass measured about two minutes per file with three other
lanes on the box. Its four slots report 68 + 68 + 67 + 67 = 270, the
file count in the tree.

### The rerun, on the SHIPPED build: IN FLIGHT, not complete

Punch round 1 changed `R/` code, so both tiers and `R CMD check` were
relaunched against the final build. **They had not finished when this
document was written, and no completed count for them is claimed
here.** State at that moment, from the same logs:

- ungated: 139 of 270 files, **0 with any failure, error or missing
  RESULT line**.
- gated: 14 of 38, same, 0.
- `R CMD check --as-cran`: still in `R CMD build`, which builds the
  vignettes.

Whoever consolidates reads the finished counts from
`dev/skewinit-log/suite-0.log` through `suite-3.log`,
`dev/skewinit-log/gated.log` and `dev/skewinit-log/check/check.log`,
and regenerates the table above with `dev/skewinit-tiers.R`. Do not
take the pre-punch table as the rerun's result.

One thing the rerun has ALREADY produced that the pre-punch run could
not: it is the instrument that caught the sparse regression described
above, at `test-sparsex.R` 49 pass / 15 fail, which every battery in
this lane had missed.

## The seventh library loss, 2026-09-22, and what its signature says

Recorded here because `dev/machine-library.md` asks for exactly this
evidence and asks for it BEFORE a restore overwrites it. Collected by
`dev/skewinit-loss-evidence.R`, log `skewinit-log/library-loss-7.txt`.

| | this loss |
|---|---|
| user-library directories emptied | 223 of 410 |
| window | 97 s, 17:05:29 to 17:07:06 |
| order | alphabetical from `abind` |
| `ZZZ-canary.txt` | SURVIVED, mtime 2026-09-09 |
| `00LOCK` directories | none, in any library |
| `rellib-r3`, outside `%LOCALAPPDATA%` | untouched, 8 of 8 |
| `pinlib`, outside `%LOCALAPPDATA%` | untouched |
| `SilentCleanup` last run | 18:26:12 the same day, result 0 |
| free space on C: | 38.06 GB of about 951 GB, about 4 percent |

This separates the two standing suspects, which the fourth and fifth
losses could not. It is NOT an interrupted install: an install hollows
the one directory it is unpacking into and leaves a `00LOCK`, and there
is no `00LOCK` anywhere while 223 directories went. It IS a file-level
sweep confined to `%LOCALAPPDATA%`, walking the tree alphabetically for
97 seconds, and the low-disk-space condition that triggers
`SilentCleanup` is currently TRUE at 4 percent free.

No R process of this lane was installing at 17:05. Its four ungated
workers and its gated run had all finished by 16:25, and the session
was cut by an account limit at about 17:20.

What it cost the lane: nothing measured before 16:25, because both
tiers had already finished and their logs are on disk.

### The restore moved the base build, so the key results were re-measured

`dev/machine-library.md` says to verify a restore with a FIT, not a
version string. The reference used here is this lane's own 27-job base
battery, measured against `rellib-r3` at 15:14, BEFORE the loss
(`dev/skewinit-verify-restore.R`, log `skewinit-log/verify-restore.txt`).

**21 of 27 base jobs came back bitwise identical; 6 moved**, by at most
1.71e-13 in log-likelihood (`gaussian_ml`, `gaussian_reml`,
`gaussian_prof`, `gaussian_sig`, `lognormal_ml`, `skew_icpt_re`). So
the restored toolchain is not bit-for-bit the one the earlier numbers
were taken on, and a base-versus-fixed pair straddling the restore
would not be trustworthy.

Both arms were therefore re-measured on the restored machine:

| result | pre-loss | post-restore |
|---|---|---|
| bitwise identical jobs | 20 of 27 | 20 of 27, same 7 moved |
| `skew_fin` relative d logLik | 3.237e-15 | 3.237e-15 |
| `skew_fin_sig` relative d logLik | 9.871e-15 | 9.871e-15 |
| acceptance, worst gap over 80 | 1.03653974293e-09 | 1.03653974293e-09 |
| acceptance, escape fired | 0 of 80 | 0 of 80 |
| `test-skew-normal-start.R` | 33 pass | 33 pass |

Every figure reproduced to the digit, so the conclusions stand. Logs:
`nomove-base-postrestore.rds`, `nomove-fixed-postrestore.rds`,
`nomove-compare-postrestore.txt`, `accept-fixed-postrestore.txt`.

Two ten-digit reference values for the NEXT restore, from this
battery's post-restore run:

    gaussian_re  (bf(ygau ~ x + (1 | g)), gaussian())  -412.711340227
    skew_icpt    (bf(ysn0 ~ 1, sigma ~ 1, alpha ~ 1))  -332.221783272

built by `dev/skewinit-nomove.R` with `set.seed(11)` and n = 300.

## Cost

**Extra optimizer runs.** 0 over the 80 acceptance fits. 20 over the
200-fit spread in `skewinit-log/falsealarm-fixed.txt`, which is 0.10
extra runs per fit, from 10 fits that each ran two.

Both of those, and the false-positive table above, were measured on the
pre-loss toolchain. They were not re-run after the restore. What WAS
re-run is the acceptance over the same 80 fits and the bitwise battery
over the same 27 jobs, and both reproduced to the digit, so the code
path these counts come from is unchanged on the current machine. Say so
rather than imply a fresh measurement.

**Wall clock.** `dev/skewinit-timing.R`, log `skewinit-log/timing.txt`.
Arms interleaved in one process, a fixed repeat count chosen once in a
warmup so no round pays for calibration, arm order rotated each round,
each block past 1.2 s, minimum of 9 rounds. On the fixed build, a fit
forced onto the stationary point takes **1.3034x** the time of the same
fit from the shipped start. The CONTROL, the same fit entered as a
second arm, reads 0.9101 rather than 1.0, so this instrument resolves a
ratio to about ten percent: call it 1.3 give or take 0.1.

It is 1.3x and not 3x because the escape reuses the taped objective:
`MakeADFun` runs once and only the `nlminb` calls repeat.

The arithmetic control read 0.00706 s per call in the fixed-build
process and 0.01154 s in the base-build process, a factor of 1.63 on
identical work, so nothing here compares absolute seconds across the
two builds. Only within-process ratios are reported.

**When it does not fire**, `stationary_escapes()` performs one
comparison per linear predictor whose family declares a stationary
point and evaluates the objective zero times. That is a statement about
the code, not a clock measurement: this box cannot resolve it.

## Other families with the same shape of defect

Audited three ways and found none:

- `sign(` across `R/families.R` and every extension: one hit inside an
  initializer, the `skew_normal` one this lane fixed. The other hits
  are `sign(y - E[Y])` in the deviance residual and a random sign in a
  simulator.
- `dpar_links_signed` (the link set a dpar that can be negative must
  use): one family, `skew_normal`. Every `alpha` in `frmtmb.learn` is a
  learning rate on a logit link, not a signed shape.
- "skew" across every extension's `R/`: no family.

Two things were found and FILED in `dev/test-backlog.md` rather than
fixed here: the `sd(y)` sigma start shared by many families, and a
second local optimum of `skew_normal` away from the stationary point.

## What I did not do

- **The escape does not run under `importance =`.** Its effective
  sample size and its proposal belong to the optimum
  `importance_fit()` returned, and replacing that optimum would leave
  both describing a point the fit no longer reports. It DOES run under
  `quadrature = TRUE`, on the Gauss-Kronrod tape, but NOT to the same
  tolerance, and the earlier 4.8e-07 was one design rather than a
  guarantee. On acceptance seed 1 with `(1 | g)` added
  (`dev/skewinit-punch1c.R`, log `skewinit-log/punch1c.txt`) the escape
  recovers 24.262 of the 24.272 units it is worth and lands
  **1.003e-02 BELOW** the unstalled fit, three orders above this
  lane's 1e-6 acceptance tolerance; the same model under Laplace lands
  9.5e-09 below. The cause is that `quad_fit()` calibrates the
  Gauss-Kronrod tape AT the Laplace optimum, which on a stalled fit is
  the stationary point, and the escape re-optimizes that frozen tape
  without recalibrating it. Keeping the escape is still far better than
  not having it, 24.262 against 0, but the honest guarantee under
  `quadrature = TRUE` is "recovers nearly all of the loss", not "to
  1e-6". Filed. Under `importance =` the residual start
  still applies, since that is in `make_start()`.
- **Autoscale.** `dev/skewinit-autoscale.R`, log
  `skewinit-log/autoscale.txt`. With `autoscale = TRUE` and the OLD
  start forced, the fit still reaches -259.338935404 at alpha 7.5065,
  and the outer escape does not fire: the autoscale pre-fit's own
  escape has already carried the warm start past the point. I could
  not construct a case where the OUTER escape fires with `par_units`
  non-NULL, so that combination is unexercised. The start vector it
  builds is in the same natural units as `opt$par`, which
  `optimize_obj()` converts through `par_units` exactly as it does for
  every other start handed to it, including `fit_recovery_starts()`'s.
- **No other family got the residual start**, including the `sigma`
  initializers that have the same weakness. Filed.
- **The second local optimum at `alpha = -0.55` is not fixed.** Filed.
- **The two test tiers ran BEFORE the library loss, and were not
  repeated after it.** Both had already completed (270 of 270 and 38 of
  38, both clean) by 16:25. They were rerun after punch round 1,
  because that round changed `R/` code; the counts in the tiers table
  are from the rerun. What was repeated post-restore is everything the
  conclusions rest on: the bitwise battery on both arms, the 80-fit
  acceptance, and `test-skew-normal-start.R`, all reproducing to the
  digit.
- **`sn::selm` is the oracle for the acceptance set only.** On the
  half-normal designs `sn` caps alpha and reports a lower likelihood
  than frmtmb, so it is not an oracle everywhere; that is stated where
  it matters above.
