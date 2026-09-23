# Lane wt-reunc: the group effects' uncertainty at a known level, and
# what `re_formula` keeps

Worktree `C:\Users\adf44\source\r\frmtmb-wt-reunc`, base `2d70b33`
(frmtmb 0.61.0, frmtmb.sample 0.9.0). Private library
`C:/Users/adf44/source/r/reunc-lib`, built by `dev/reunc-install.R`; the
round's read-only base build is `C:/Users/adf44/source/r/rellib-r3` and
nothing is installed into it. Every number below names the script that
produced it and the log under `dev/reunc-log/`.

## 0. Answer first

1. **`predict()` carries the group effects at a level the fit saw.**
   Each replicate draws the WHOLE vector of group effects from its
   conditional law given the data and that replicate's parameters, read
   from the fit's joint precision, so rows of one group share the
   group's draw and rows of different groups carry exactly the
   correlation the model implies. Measured against the analytic
   conditional variance and against `a_i' V a_j` (section 2).
2. **`fitted()` at a known level already carried it, on the scalar
   route.** The brief says both methods were conditional on the modes;
   that is right for `predict()` and WRONG for `fitted()`, whose
   delta-method standard error has read the joint covariance of
   `(beta, b)` since before this lane. What did NOT carry it is
   `fitted()` on an ordinal or categorical fit, which differenced the
   outer parameters only: 0.59 of the Monte Carlo value on a cumulative
   fit. The group effects join that difference now (section 2.4).
3. **A partial `re_formula` is honored**, with brms's rule, in
   `predict()`, `fitted()`, `frm_linpred()`, `frm_lp_basis()` and
   everything that reads them, and only the kept terms contribute their
   uncertainty. Measured against brms 2.23.0 on one data set, five
   formulas including two random-slope selections (section 4).
4. **`re_formula = ~(1 | nosuch)` is refused**, classed and naming the
   term. brms does NOT refuse it: it drops an unmatched term silently,
   which makes the answer the population-level one. This is a
   deliberate departure and it needs the user's eye (sections 5, 11).
5. **On draws (frmtmb.sample) nothing was added**: each draw carries
   its own sampled `b`, so the known-level uncertainty was already
   there. What was broken there is the same partial `re_formula`, and
   the core fix reaches it (section 6).
6. **Nothing else moved.** 330 recorded quantities are `identical()`
   between base and lane across five models and three estimation modes,
   including `vcov()`, `vcov(full = TRUE)`, `fixef()`, `summary()`,
   every `re_formula = NA` prediction, `predict()` at rows whose level
   the fit never saw, and the caller's RNG position after four kinds of
   call. The nine that moved are the conditional `predict()`
   (section 7).
7. **The tiers are green, and section 9.0 says which build each one
   covers.** Ungated suite 272 of 272 files and 16,027 assertions with
   no failure, and gated tier 38 of 38 files and 3,130 assertions, the
   release's own figure, with every gate variable set: both of those
   are the build BEFORE the `propagate_error` rename. What covers the
   shipped tree is `R CMD check --as-cran` on frmtmb and
   frmtmb.sample, whose `checking tests` runs both full suites with
   `NOT_CRAN` set (FAIL 0, PASS 10530 and 1791), plus the review's
   comparison of all six extension suites against base, which the
   check does not reach. The brms-port ledger's totals are unchanged.
   All on RTMB 2.0, against `reunc2-lib`.
8. **The machine lost its user library again at 17:06**, after every
   tier here had finished and with nothing of this lane running.
   Section 9a has the evidence `dev/machine-library.md` asks for.

## 1. The coverage claim, stated before the measurement

`predict()`'s interval at a grouping level the fit SAW is a prediction
interval for a NEW observation at that level. Its claim is coverage
AVERAGED over data sets and over the group effects: the truth draws
`u ~ N(0, tau^2)` afresh in every replicate, the fit is made, and the
new points land at levels chosen at random.

It does NOT claim coverage conditional on one group's realized effect.
The predictor shrinks a group toward the population, so for a group
whose effect is far from zero the interval sits too close to the
population and covers less than nominal, and for a group near zero it
covers more. This is a property of the quantity, not of the
implementation: the exact interval built from Henderson's prediction
error variance at the TRUE variance parameters behaves the same way,
which the tercile table in section 3 shows.

`fitted()`'s interval claims the same average coverage for the
CONDITIONAL MEAN at that level, `x'beta + u_g`, with no observation
noise.

The positive control is exact rather than approximate. With the
variance parameters known, the mixed model equations give the best
linear unbiased predictor and its prediction error variance
`a' C^-1 a`, and the prediction error is exactly normal with variance
`sigma^2 + a' C^-1 a` marginally over `u` and the noise. An interval
built from it covers 0.95 by construction, so an arm that reports
anything else says the harness is wrong rather than `predict()`.

## 2. Item 1: the construction

### 2.1 What is drawn, and why that and not a per-row variance

`sdreport()` gives the joint precision `Q` of every estimated
parameter, fixed and random. Partition it into the parameters a
replicate draws, `d` (the outer vector, plus `beta` under `REML = TRUE`
or `control(profile = TRUE)`, which is what `fit_draw_space()` already
assembled), and the rest `r`, which holds `b`:

    r | d  ~  N(r_hat - Q_rr^-1 Q_rd (d - d_hat),  Q_rr^-1).

For a linear mixed model that is Henderson's prediction error variance
of the BLUP, jointly with the fixed effects: `Q_rr` is the inner Hessian
and the off-diagonal block carries how the modes move with the
parameters, so the draw keeps the negative correlation between a group
effect and the intercept instead of adding an independent variance on
top of it.

**What the replicate's parameters do and do not change.** `Q` is read
ONCE, at the estimates. A replicate's own parameters enter through
`d - d_hat`, which shifts the conditional MEAN; the conditional
VARIANCE stays `Q_rr^-1` at the estimates rather than being
re-evaluated at each draw, which would mean a fresh joint precision per
replicate. So this is not "Henderson's PEV with the variance parameters
at their draw", which an earlier wording of this section and of
`predict_b_drawer()`'s comment claimed; the construction is the
linearization around the estimates, and the uncertainty in the variance
parameters themselves is the term section 3 measures as still missing
from every interval here.

The draw is of the whole vector, which is what makes
`predict(summary = FALSE)` right across rows: two rows of one group
load one draw of that group's effect. A per-row added variance would
get each row's summary right and the joint law wrong, which is the
shared-new-level-draw defect of the last round in reverse.

`predict_b_drawer()` (`R/predict-brms.R`) takes the sparse Cholesky of
`Q_rr` once per call and draws `P' L^-T z` per replicate, with the
conditional mean shift solved against the replicate's own parameter
deviation. Only the blocks `re_formula = NA` would drop are written
back: a population smooth, a `gp()` and an `hsgp()` curve stay at their
modes, because `re_formula = NA` keeps them too.

`propagate_error = FALSE` holds the group effects too, along with the
parameters. This was the other way round when the lane was built, on
the reasoning that an effect is not a parameter but an unobserved
random variable, and the user settled it on 2026-09-23 with the rename
(section 13): the argument names an AXIS, whether the error in the
estimates is propagated, and a fitted group effect is an estimate. So
`FALSE` leaves the observation noise alone.

### 2.2 The joint structure, measured

`dev/reunc-joint.R`, logs `dev/reunc-log/joint-ML.txt`, `joint-REML.txt`
and `joint-profile.txt`. Gaussian one-way random intercept, 6 groups of
4, so the conditional law at the estimates is analytic:
`Var(b_g | y, beta, sigma, tau) = sigma^2 tau^2 / (sigma^2 + n_g tau^2)`.
Four rows: two in group 1, one in group 2, one in group 3. 20,000
draws.

| mode | same-group cov | analytic | ratio | other-group cov | MC se |
|---|---|---|---|---|---|
| ML | 0.11383 | 0.11482 | 0.9914 | -0.00432 | 0.00506 |
| REML | 0.12499 | 0.12605 | 0.9916 | -0.00459 | 0.00536 |
| profile | 0.11383 | 0.11482 | 0.9914 | -0.00432 | 0.00506 |

**Remeasured after the 2026-09-23 rename, with a changed
construction.** Holding the parameters while drawing `b` used to be
`param_uncertainty = FALSE`, and `propagate_error = FALSE` now holds
`b` as well, so that call no longer isolates this law. The script
drives `predict_b_drawer()` directly instead and adds the noise
itself, which is what the old flag did internally. The numbers moved
within Monte Carlo error, because the noise now comes off a different
point in the stream: ML same-group cov 0.11146 before and 0.11383
after, against an analytic 0.11482 and an MC standard error of
0.00506. The construction is named here because a number that changed
for a reason should say what the reason was.

**What that table does and does not cover, now that it drives the
drawer.** It says the DRAWER implements the conditional law. It no
longer says, on its own, that `predict()` is wired to the drawer,
because `predict()` is no longer on its path. The end-to-end claim is
carried by the next table instead, which goes through public
`predict()` with no internals in the call, and by
`tests/testthat/test-predict-re-uncertainty.R`, which checks both
halves separately. Read the two tables as one argument: the first
fixes the law, the second fixes the wiring.

With the parameter draw on, the model-implied covariance of two rows is
`a_i' V a_j` from `frm_lp_basis()`, which is the delta method's own
matrix and therefore an outside number for the draws. This one is
`predict()` end to end, at its default `propagate_error = TRUE`:

| pair | draws | `a_i' V a_j` | z |
|---|---|---|---|
| rows 1, 2 (same group) | 0.09988 | 0.10778 | -1.44 |
| rows 1, 3 (different groups) | 0.01164 | 0.00907 | 0.47 |
| rows 3, 4 (different groups) | 0.02192 | 0.02039 | 0.24 |

Rows of different groups covary because they share the fixed effects,
and the amount is the one the joint covariance implies, not zero.

**A signal that dissolved, recorded because it is a result about the
count.** At 20,000 draws the plug-in variance of one row came out at
0.9794 of `Var(b_g | y) + sigma^2`, z = -2.1, on both designs, and the
same-group covariance at 0.9708. Two instruments settled it.
`dev/reunc-bdraw.R` measures the group-effect draw DIRECTLY against
`Q_bb^-1`, without the response on the path: 50,000 draws on three
seeds give diagonal ratios 0.9926 to 1.0170, against a Monte Carlo
standard error of 0.0063 per entry, so the draw itself is right.
`dev/reunc-plugvar.R` then repeated the response variance at 100,000
draws on three seeds: 0.9910, 1.0034, 1.0038 (z = -2.02, 0.77, 0.85),
and after the rename and the same construction change, 1.0038, 0.9981
and 0.9966 (z = 0.84, -0.43, -0.76).
The 0.9794 was an unlucky stream at too small a count, and the two runs
that produced it shared a seed. A third instrument, `dev/reunc-seedcorr.R`,
ruled out the mechanism I suspected first, that the group-effect draw
and the simulation draw take two streams seeded from one `sample.int()`
call and could be correlated: 200,000 pairs correlate at 0.0055, and
the sign is the wrong way for a shortfall.

The variance of one row exceeds `a' V a + sigma^2` by 0.03 to 0.045 on
all four rows. That is not a defect and it is not Monte Carlo: the
parameter draw perturbs `log sigma` too, so the draws' variance is
`sigma^2 E[exp(2 (log sigma draw - log sigma hat))]`, about
`sigma^2 (1 + 2 Var(log sigma))`, and `Var(log sigma)` is about 0.025
on 24 rows. `0.6 * 0.05 = 0.03`.

### 2.3 fitted(), and a correction to the brief

The brief says `predict()` AND `fitted()` are conditional on the modes
today. For `fitted()`'s scalar route that is wrong, and the base build
says so: on a `(1 + x | g) + (1 | h)` fit the base `Est.Error` at a
known level is 0.331315 against 0.2416937 at `re_formula = NA`
(`dev/reunc-explore.R`). `lp_delta_A()` puts the `Z` columns into `A`
and `get_joint_cov()` inverts the joint precision, so the `b` block is
in the quadratic form. `?fitted.frmtmb_fit` already said as much in its
ordinal section ("which the scalar route does"), while its own
description said "conditional on the modes"; the description is
corrected rather than the code.

So `fitted()`'s scalar route is UNCHANGED by this lane, bitwise
(section 7), and what it reports is the delta-method twin of
`predict()`'s draw: both read the same conditional law.

### 2.4 The ordinal and categorical route, which did not carry it

A category probability depends on the thresholds and the `cs()`
coefficients as well as on the predictor, so `fitted()` differences it
(`fit_fd_se()`). That route covered the OUTER parameters only. The
group effects of the kept blocks join the differenced vector now, with
their covariance taken jointly with the parameters from the joint
precision.

`dev/reunc-ordinal.R`, logs `dev/reunc-log/ordinal-base.txt` and
`ordinal-lane.txt`. Cumulative fit, 10 groups of 12, two prediction
rows. The reference is a brute-force Monte Carlo over the same joint
law, 4000 draws pushed through `fitted()`'s own estimate:

| row, category | base Est.Error | lane | Monte Carlo sd | base / MC | lane / MC |
|---|---|---|---|---|---|
| 1, P(Y = 1) | 0.05042 | 0.08445 | 0.08516 | 0.592 | 0.992 |
| 1, P(Y = 2) | 0.05429 | 0.05707 | 0.05708 | 0.951 | 1.000 |
| 1, P(Y = 3) | 0.07380 | 0.11451 | 0.10740 | 0.687 | 1.066 |
| 2, P(Y = 1) | 0.04513 | 0.06495 | 0.07029 | 0.642 | 0.924 |
| 2, P(Y = 2) | 0.05747 | 0.07131 | 0.06516 | 0.882 | 1.095 |
| 2, P(Y = 3) | 0.08444 | 0.12282 | 0.11578 | 0.729 | 1.061 |

The delta method is first order and a probability is curved, so the
right ratio is near 1 and not exactly 1. The base build also reported a
SMALLER standard error at a known level than at `re_formula = NA` for
`P(Y = 1)` (0.05042 against 0.06520), which is the visible symptom.

The cost is one function evaluation per group effect and direction:
`fitted()` on that fit went from 0.14 s to 0.39 s.

### 2.5 What the draw costs predict()

`dev/reunc-cost.R`, interleaved arms in one process, each block grown
past 1.2 s and the minimum of five rounds taken, because `proc.time()`
ticks at 10 ms here. 1,200 rows, 60 groups, 50 prediction rows, 200
draws. The `population` arm (`re_formula = NA`) is the control this
lane does not touch, and a fixed arithmetic loop is the second control.

| arm | base | lane | ratio |
|---|---|---|---|
| `predict()` conditional | 0.4867 s | 0.5567 s | 1.14 |
| `predict(re_formula = NA)` | 0.2920 s | 0.2940 s | 1.01 |
| arithmetic control | 0.0043 s | 0.0041 s | 0.95 |

So the group-effect draw costs about 14 percent of a conditional
`predict()` on this model, and the two controls say the machine
contributed about 1 percent of that. A fit also needs its joint
precision, which an ML `sdreport()` does not carry: 0.060 s on this fit
(minimum of five fresh fits), which is the same work `fitted()` already
does for its standard errors, and it is memoized per fit.

## 3. Coverage

`dev/reunc-coverage.R` on two designs, each run on the BASE build and
on the lane build with the same seeds, and summarized by
`dev/reunc-coverage-sum.R` into `dev/reunc-log/coverage-summary.txt`.
300 replicates of 10 fresh points is 3,000 predictions per arm; the
points of one replicate share a fit, so the standard error quoted is
the replicate-level one and z is against 0.95 with it. The binomial
interval is printed beside it and is too narrow for the same reason.

* `mixed12`: 60 rows, 12 groups of 5, `sigma = 1`, `tau = 0.7`, which
  is `dev/shapes-coverage.R`'s mixed design.
* `few4`: 40 rows, 4 groups of 10, the same variances. Few groups is
  where `tau` is estimated worst.

The controls first, because nothing below means anything without them.
`oracle_pred` covers 0.9523 and 0.9513 (z = +0.59, +0.34), so the
harness can report 0.95; `oracle_mean` covers 0.9480 and 0.9480 for the
conditional mean. `wald_obs`, `fitted()`'s interval asked to cover an
OBSERVATION, covers 0.5340 and 0.4433 (z = -46, -52), so a number near
0.95 elsewhere is not an artifact.

| design | arm | base | lane | base z | lane z |
|---|---|---|---|---|---|
| mixed12 | `predict()` ML | 0.9313 | **0.9407** | -3.67 | -1.96 |
| mixed12 | `predict()` REML | 0.9343 | **0.9450** | -3.14 | -1.06 |
| mixed12 | `predict()` profile | 0.9303 | **0.9407** | -3.84 | -1.99 |
| mixed12 | `predict()` plug-in, OLD flag | 0.9133 | 0.9317 | -6.60 | -3.68 |
| mixed12 | oracle, true variances | 0.9523 | 0.9523 | +0.59 | +0.59 |
| few4 | `predict()` ML | 0.9393 | 0.9383 | -2.24 | -2.45 |
| few4 | `predict()` REML | 0.9457 | 0.9420 | -0.96 | -1.75 |
| few4 | `predict()` profile | 0.9400 | 0.9383 | -2.10 | -2.43 |
| few4 | `predict()` plug-in, OLD flag | 0.9177 | 0.9243 | -5.96 | -4.87 |
| few4 | oracle, true variances | 0.9513 | 0.9513 | +0.34 | +0.34 |

The replicate-level standard error is 0.0045 to 0.0056 on every
`predict()` arm.

**The two plug-in rows measure a setting that no longer exists**, and
they are kept rather than deleted because they are the evidence for
the decision that replaced it. They are
`param_uncertainty = FALSE`, which held the parameters and still drew
the group effects. Since the 2026-09-23 rename (section 13),
`propagate_error = FALSE` holds both, so the same arm would now be an
interval carrying the observation noise alone, which under-covers by
construction and is not worth 300 replicates to confirm. Those rows
are therefore the OLD semantics, labelled as such, and NOT rerun. The
`predict()` rows above them are the shipped behavior and are
unaffected, because they are the `TRUE` setting, which the rename does
not change.

**mixed12 is the design the change was for, and it moves the right
way**: +0.0094 under ML and +0.0107 under REML, from z = -3.67 to
-1.96. The median width goes from 4.0254 to 4.1839.

**few4 does NOT improve, and that is not a failure of the
construction.** It is the design where holding the group effects at
their modes was wrong in TWO directions at once. With 4 groups the
intercept is badly determined, `Var(beta_0 hat)` is large, and the base
build drew `beta_0` while holding `b_g` fixed: that OVERSTATES the
spread of `x'beta + b_g`, because the two are strongly negatively
correlated and the group mean is much better determined than either
part. Base's 0.9393 is two errors cancelling. The lane draws the pair
jointly, which is the correct law, and the median width falls from
4.1176 to 4.0525 while the coverage stays at 0.9383.

**What is left, measured rather than asserted.** `dev/reunc-attrib.R`
runs the same designs, the same seeds and the same data, and builds the
EXACT interval with the fit's own estimated `sigma` and `tau` in place
of the true ones. No frmtmb prediction code is on its path.

| design | exact at true variances | exact at the ML estimates | at the REML estimates |
|---|---|---|---|
| mixed12 | 0.9523 (z +0.59) | 0.9403 (z -2.06) | 0.9437 (z -1.40) |
| few4 | 0.9513 (z +0.34) | 0.9337 (z -3.38) | 0.9380 (z -2.56) |

`predict()`'s own arms are 0.9407 and 0.9450 on mixed12 and 0.9383 and
0.9420 on few4. They sit on top of the exact interval computed at the
same estimates, within a replicate-level standard error, in all four
cells. So what remains after this lane is the uncertainty in the
variance parameters themselves, which no interval built at a point
estimate of `tau` can carry, and not the construction; and REML, whose
`tau` is less biased, covers better than ML in every arm and in both
designs.

**A shortfall this lane did NOT close, and did not introduce.**
`fitted()`'s interval covers the conditional mean at 0.9170 on mixed12
and 0.8977 on few4, identically on both builds (the scalar route does
not change here), against the oracle's 0.9480. It is the same cause:
the interval is a Wald interval at an estimated `tau`, and with 4 or 12
groups that estimate is poor. brms's posterior carries the uncertainty
in `tau` itself, which is what a frequentist interval at a point
estimate cannot. Filed in `dev/test-backlog.md`; closing it means a
wider interval than the delta method gives, not a different group-level
term.

**The claim is the average, not the group.** The tercile tables in
`dev/reunc-log/coverage-summary.txt` show the conditional behavior the
claim explicitly does not make. On mixed12, `fitted()`'s interval
covers 0.9590 for the third of points whose group effect is smallest in
absolute value and 0.8690 for the largest third; the ORACLE at the true
variances does the same, 0.9670 against 0.9180. Shrinkage, not an
implementation fault.

**What the last round's two arms measured** (`dev/shapes-coverage.R`,
condvar 0.9609; `dev/shapes-rev-cov.R`, joint_plus_condvar 0.9477).
Both measured this same average claim: both redraw `u` per replicate
and both put the fresh points at levels chosen at random. What differs
is the construction: both widened the modes-conditional interval by
`v_g = sigma^2 tau^2 / (sigma^2 + n_g tau^2)`, the conditional variance
of `u_g` GIVEN beta, on top of an interval that already carried
`x' Var(beta hat) x`. That sum is not the prediction error variance,
because the BLUP moves against the intercept. `dev/reunc-condvar.R`
computes exactly how much, at the true variances, so the arithmetic is
separated from the estimation:

```
script dev/reunc-condvar.R, seed 1, 4000 replicates of 10 points
shapes  sigma 1.0 tau 0.7 groups 12: known-parameter coverage of
  sigma^2 + x'Vx + v_g = 0.9547 (exact PEV interval 0.9500);
  (x'Vx + v_g) / PEV mean 1.311
review  sigma 0.8 tau 0.7 groups 12: 0.9567; ratio 1.416
few4    sigma 1.0 tau 0.7 groups  4: 0.9625; ratio 2.315
```

So the sum over-counts the variance by 31 to 42 percent on those two
designs and would over-cover by 0.005 to 0.007 even with the variance
parameters known, while estimating them pulls coverage down. The two
arms are on opposite sides of nominal for that reason, and neither was
a proposal: each was a test of an attribution, which they settled in
the direction this lane implemented.

## 4. Item 2: a partial re_formula, against brms

`R/re-formula.R` resolves the argument once per call and returns a VIEW
of the fit: the dropped terms are taken out of the prediction design,
in sample by zeroing those columns of each predictor's `Z` and on new
data by removing the component from its block (or masking the columns
it keeps). Every prediction path already reads the design from the
frame, so one construction covers the estimate, the delta-method
standard error, the new-level variance and `predict()`'s draws. The
resolution then hands the callers `re_formula = NULL` on the reduced
design.

The rule is brms's `check_re_formula()`: a formula term is kept when
the fit has a term with the same grouping factor whose columns include
the formula term's columns. `~ (1 | g)` on a `(1 + x | g)` fit
therefore keeps the intercept and drops the slope. brms's spellings are
read as brms reads them: an id (`(1 | p | g)`), `gr()`, `||` and a
nested group `a/b`, with a non-bar term ignored.

`dev/reunc-brms.R`, logs `dev/reunc-log/brms-lane.txt` and
`brms-base.txt`. One gaussian data set, 120 rows, fit
`y ~ x + (1 + x | g) + (1 | h)` with 15 `g` levels and 5 `h` levels;
brms 2.23.0, 2 chains of 3000. Two measurements per formula: brms
against the hand construction from its OWN draws (which says what brms
keeps), and the dropped contribution `fitted(NULL) - fitted(formula)`
of frmtmb against brms's posterior mean of the same difference.

| re_formula | brms vs its hand construction | dropped: brms sd | frmtmb sd | max abs diff | cor |
|---|---|---|---|---|---|
| `~(1 \| g)` | 0 | 0.7716 | 0.7542 | 0.1890 | 0.9983 |
| `~(0 + x \| g)` | 0 | 1.0576 | 1.0449 | 0.0594 | 0.9999 |
| `~(1 + x \| g)` | 8.9e-16 | 0.7671 | 0.7541 | 0.0241 | 0.9999 |
| `~(1 \| h)` | 0 | 0.7393 | 0.7304 | 0.2106 | 0.9983 |
| `~(1 \| g) + (1 \| h)` | 8.9e-16 | 0.0785 | 0.0336 | 0.2107 | 0.9709 |

The first column says brms keeps exactly the columns the rule predicts,
to machine precision. The rest is a comparison of two ESTIMATORS of the
same contribution, a posterior mean against a mode, so it is close and
not equal; on the base build the frmtmb column is 0.0000 in every row,
because the base keeps every term. Two of the five formulas select
within a random-slope term, which is what the brief asked for.

A dropped term's grouping column is no longer needed in `newdata`:
brms predicts `nd` without the `h` column under `~(1 | g)`; the base
build refuses it ("newdata has no column `h`") and the lane answers
0.2737, 1.5885, 1.3493 against brms's 0.2550, 1.5980, 1.3095. That also
removes the divergence lane `wt-adefects` filed with its D6 fix
(`dev/adefects-findings.md` section 11, item 7).

Two cases are REFUSED rather than guessed.

* A formula term that matches no term of the fit (section 5).
* A formula that keeps SOME terms on a fit that also has group-level
  content a formula cannot name: a factor-smooth term, `car()`,
  `spde()`. `NA` drops that content and `NULL` keeps it, and a partial
  formula has no way to say which. Naming every bar term is still
  `NULL`, so the smooth is kept there, which is what brms does under
  any `re_formula` and what frmtmb does under `NULL`.

`re_formula = ~1` and `~0` now mean the same thing, "no group-level
effects", as they do in brms. `~1` used to keep every term.

## 5. Item 3: a term the fit does not have

brms does not refuse it. `brms:::check_re_formula(~(1 | nosuch), f)`
returns `~1`, and `posterior_epred(bfit, re_formula = ~(1 | nosuch))`
is `identical()` to `re_formula = NA`
(`dev/reunc-log/brms-lane.txt`). The same happens to a term whose
grouping factor matches but whose columns do not: `~(1 + z | g)` on a
`(1 + x | g)` fit is dropped silently too.

frmtmb refuses, classed `frmtmb_error`, naming the term and listing the
fit's own terms. The reason is the one the whole item is about: a
misspelled grouping factor would otherwise change the answer with
nothing said, and here it would change it to a DIFFERENT wrong thing
from before (the population prediction rather than the full one).
`re_formula = NA` says "no group effects" in one unambiguous way.

This is a deliberate divergence from the tiebreaker and section 11
puts it in front of the user.

## 6. Item 4: draws (frmtmb.sample)

Verified rather than assumed. Each draws method evaluates one draw
through `frmtmb::frm_linpred()`, so the partial-`re_formula` defect was
there in full and the core fix reaches it with no change to
frmtmb.sample's own code. `dev/reunc-brms.R` section 3, 60 draws of the
same model, checking `posterior_epred()` draw by draw against the hand
construction from that draw's own `b`:

| re_formula | base max abs diff | lane |
|---|---|---|
| `~(1 \| g)` | 2.98 | 0 |
| `~(0 + x \| g)` | 3.86 | 0 |
| `~(1 + x \| g)` | 2.78 | 8.9e-16 |
| `~(1 \| h)` | 3.44 | 0 |
| `~(1 \| g) + (1 \| h)` | 1.75 | 8.9e-16 |

`fitted()` on draws is `colMeans(posterior_epred())` under every
formula, and `predict()` is the summary of `posterior_predict()` at the
same seed, in both builds.

The group effects' uncertainty at a known level needs nothing on
draws: every draw carries its own sampled `b`, which is the posterior
brms's interval carries. `dev/reunc-log/lane-test-re-formula-draws.R.txt`
pins it: the draws of `posterior_predict()` scatter around the kept
terms' `posterior_epred()` and not around the full one.

## 7. The bitwise-unchanged proof

`dev/reunc-bitwise.R` records 21 quantities per fit on five models
(gaussian with `(1 | g)`, poisson with `(1 | g)`, a distributional
`sigma ~ x` with `(1 | g)`, a model with NO group-level term, and a
population-smooth model) in each of ML, `REML = TRUE` and
`control(profile = TRUE)`: `vcov()`, `vcov(full = TRUE)`, `fixef()`,
`summary()$fixed`, `$spec_pars`, `$random`, `logLik()`, `frm_linpred()`
with and without `newdata` at `re_formula = NA` with standard errors,
`fitted()` at `NA` and at `NULL`, three `predict()` calls at `NA`
(summary, draws, plug-in), the caller's RNG position after them,
`predict()` at rows whose level the fit never saw, and `fitted()` there.
Punch round 1 added four more per mixed fit, all stream positions:
after a conditional `predict(newdata = )`, after
`predict(allow_new_levels = TRUE)` at unseen levels, after an in-sample
`predict(summary = FALSE)` and after `fitted()`. Round 1 recorded the
stream only after a POPULATION call, where the group-effect drawer
returns before taking anything, and a reviewer's own rebuild found 9
differing positions there (section 12, M2).

`dev/reunc-bitwise-cmp.R` compares the base and lane files with
`identical()`. Both sides were regenerated by the second session of
punch round 1, on `rellib-r3` and on `reunc2-lib`, RTMB 2.0, after the
last source edit of the round:

    identical 330, moved as intended 9, differs 0

The nine are `predict_null`, the conditional `predict()` on the three
mixed models in the three modes, which is the change. Everything else,
including `predict()` at an unseen level and the RNG stream position
after a `re_formula = NA` call, is bit for bit the base build's. The
group-effect draw takes its seeds from the caller's stream only when it
has something to draw, and each replicate draws from its own seed, so
the parameter draws and the simulation seeds are exactly what they were.

## 8. Every test seen failing first

Run against the base build `rellib-r3` with the lane's copy of the
release runner (`dev/reunc-run-tests.R`, `REUNC_LIB` and
`REUNC_REPORTER=check`), logs kept under `dev/reunc-log/`:

| file | base | lane |
|---|---|---|
| `tests/testthat/test-re-formula-partial.R` | pass 8, fail 25, error 3 | pass 39, fail 0 |
| `tests/testthat/test-predict-re-uncertainty.R` | pass 15, fail 12, error 3 | pass 47, fail 0 |
| `frmtmb.sample/test-re-formula-draws.R` | pass 10, fail 12 | pass 22, fail 0 |

The predict file grew across both punch rounds, from 24 assertions to
47, and its base counts moved with it (the lane's first session
recorded pass 11, fail 13). Of the four blocks added, three ERROR on
base rather than failing, because `re_governed_b()` and
`block_b_positionwise()` do not exist there, which `dev/lane-rules.md`
calls the weak form of seeing a test fail; the behavioral evidence for
those is the ulp and relative-error tables in section 12, each with
the defect reconstructed and measured. The M2 block PASSES on base,
which is also not evidence, and section 12 says why and what was seen
to fail instead.

`seen-failing-base-<file>.txt` holds the failures in full. The base
failures are behavioral, not "no such function": the partial formulas
return the full prediction (five cases), the unmatched term is accepted
(nine expectations over three functions), a dropped term's grouping
column is still demanded by `newdata` (an error), and on the
uncertainty file the same-group covariance of the draws is 17.8 and
18.6 standard errors from the analytic value under ML and REML, the
different-group covariance 14 to 16, and the ordinal ratio 0.598 rather
than near 1.

One existing file needed a one-line fix of its own, and it found a real
hazard: `test-review-v25.R` replaces the stored `sdreport`'s joint
precision with a singular matrix and expects the next `se.fit` to
degrade. `joint_precision()` memoized the matrix, so the replacement
was ignored and the test failed (`pass=58 fail=2`). The memo now holds
only the separate joint-precision `sdreport()` an ML fit needs, and the
stored one wins whenever it has a matrix; the file passes 60 of 60.

## 9. Suites and checks

The lane's copies of the release drivers change three things and
nothing else: the tree, the runner (`dev/reunc-run-tests.R`, which
puts the lane library first and is otherwise `dev/release/run-tests.R`),
and the log path. The release scripts themselves are untouched. They
take the library from `REUNC_LIB`, as the runner already did, so a
second session does not install into the first's library.

### 9.0 Which tier covers which build

**Not every tier below covers the shipped code, and this is where that
is stated rather than left to be inferred.** The `propagate_error`
rename (section 13) landed after the ungated suite and the gated tier
had run, so those two describe the build as it stood before it.

| tier | log written | covers |
|---|---|---|
| ungated suite, 272 files | 12:09:51 | PRE-rename |
| gated tier, 38 files | 12:51:42 | PRE-rename |
| `R CMD check --as-cran`, both packages | 15:25:10 | the shipped tree |

The last source edit is `R/predict-brms.R` at 14:05:41 and the last
test edit is `test-predict-re-uncertainty.R` at 14:21:27, both after
the first two logs and both before the check.

**Quantified, so the gap is a number and not an adjective.**
`suite.log` records `test-predict-re-uncertainty.R pass=47`; the final
tree gives 65; and 65 - 47 = 18 is exactly the rename block, the same
18 that takes the check's frmtmb total from 10512 to 10530. So the
"272 files, 16,027 passing" below is the PRE-rename total and the
post-rename total is **16,045**. The gated 3,130 is in the same
position: it names no file the rename touched, so the figure carries
over, but it was not rerun against the shipped tree.

**What does cover the shipped code**, and it is enough:

* **`R CMD check --as-cran` at 15:25:10**, after every edit. Its
  `checking tests` runs both packages' full suites in one process with
  `NOT_CRAN` set: frmtmb FAIL 0, PASS 10530 and frmtmb.sample FAIL 0,
  PASS 1791. That is the whole ungated tier for the two packages the
  lane changes.
* **The six extension suites**, which `R CMD check` on frmtmb and
  frmtmb.sample does NOT reach. The review ran all six against the
  base build and found them identical, which closes the one real gap
  in the reasoning above. `param_uncertainty` never appeared under
  `extensions/` in any of the seven packages, so none of them could
  name the removed argument, and the comparison confirms nothing
  reached them another way.
* **The bitwise battery**, rerun after the rename:
  `identical 330, moved as intended 9, differs 0`, including
  `predict_na_plug`, the one entry that passes the flag.

**No rerun of the two pre-rename tiers was done here**, deliberately:
the consolidation rerun covers the merged build anyway, and repeating
an hour of tiers to move one file's count by 18 buys less than saying
which build each number describes. The obligation this section
discharges is to say so plainly, because as written before it invited
the reader to think the 272-file suite had seen the rename.

**Every log in this document postdates the source it describes, and
that was checked rather than assumed.** The review found two logs
carrying numbers taken before the last source edit of punch round 1,
while section 12 opened by listing exactly that class of staleness for
the round before it. They were rerun and they reproduced. Punch round
2 then checked the whole set by mtime and found four more in the same
position, behind a roxygen-comment edit to `R/predict-brms.R`: the
positionwise enumeration, the car battery, the stream reconstruction
and both memory logs. The package was reinstalled and all five were
rerun. Every number reproduced exactly, including the memory table to
the decimal and the esicar ratios to four places.

The ungated suite and the gated tier ran on the build from before that
comment edit as well as before the rename, which section 9.0 tabulates.
The comment gap is backed by a measurement: after the reinstall, the
bitwise battery
returned `identical 330, moved as intended 9, differs 0` against the
same base file, and `test-predict-re-uncertainty.R` (47),
`test-ordinal.R` (109), `test-ordinal-fitted.R` (91),
`test-car-spde.R` (136), `test-re-formula-partial.R` (39),
`test-review-v25.R` (60) and `test-covstruct.R` (7) all returned the
counts they returned before it. The two builds are behaviorally the
same object.

**A shell trap that cost this round half an hour, and that
`dev/lane-rules.md` does not list.** On Windows, holding a log open
with `tail -f` blocks `Add-Content` on the same file. The driver runs
under `$ErrorActionPreference = "Continue"` around the R calls, for
the reason its own comment gives, so each failed write printed an
`IOException` to the driver's stdout and the loop carried on. Its
`$ran` counter increments on the RESULT line and not on the write, so
the run would have finished and reported `SUITE ran 272 of 272` over a
log holding 17 of them. This is the "a count is not a count" failure
in a new spelling: the COUNT was right and the EVIDENCE was gone. Do
not tail a driver's log while it runs. That run was discarded and
repeated.

**Ungated suite** (`dev/reunc-suite.ps1`, log `dev/reunc-log/suite.log`),
one R process per test file, all eight packages:

    SUITE ran 272 of 272

272 files, 16,027 passing assertions, 161 skips, 0 failures and 0
errors, summed from the log's own RESULT lines, with a RESULT line for
every file. **This is the pre-rename build** (section 9.0): the
post-rename total is 16,045, the difference being the 18 assertions
the rename added to one file.

272 against the release's 269 because this lane adds three files:
frmtmb 159, frmtmb.sample 31, eam 26, spline 14, learn 14,
coupling 9, latent 9, ode 10.

Where the assertion total comes from, because a total that moves for
an unnamed reason is exactly what this record exists to catch:

| run | assertions | what changed |
|---|---|---|
| 0.61.0 release | 15,919 | |
| lane, first session | 16,004 | the lane's three new files |
| punch round 1 | 16,014 | +10, its M1, M2 and factor-smooth guards |
| punch round 2 | 16,027 | +13, the esicar guard and 3 preconditions |
| the rename | 16,045 | +18, section 13; NOT rerun as a suite |

`test-predict-re-uncertainty.R` carries all of the movement: 24
assertions when the lane was built, 34 after punch round 1, 47 after
punch round 2 and 65 after the rename. The last step is the one the
suite log above predates, and the check at 15:25:10 is what covers it.

**Gated tier** (`dev/reunc-gated.ps1`, log `dev/reunc-log/gated.log`),
with every environment variable the release driver sets: `NOT_CRAN`,
`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_FUZZ`, `FRMTMB_STAN_CACHE`:

    GATED ran 38 of 38

38 files, 3,130 passing assertions, 0 failures, 0 errors, 0 skips,
summed from the log's own RESULT lines rather than typed. That is the
0.61.0 release figure exactly (`dev/round-handoff.md`: "Gated: 38 of 38
files, 3,130 assertions, 0 fail"), and it includes the four
brms-agreement files whose gate an earlier lane forgot to set. Rerun
by the second session on `reunc2-lib`, after the `re_row_support()`
fix, with the same total. **Pre-rename** (section 9.0): no gated file
names `propagate_error`, so the figure carries over, but it was not
rerun against the shipped tree.

**`R CMD check --as-cran`** (`dev/reunc-check.ps1`, log
`dev/reunc-log/check.log`), built WITH vignettes:

| package | run 1 | run 2 | run 3 | run 4, final |
|---|---|---|---|---|
| frmtmb | 2 NOTEs | 1 NOTE | 2 NOTEs | **1 NOTE** |
| frmtmb.sample | OK | OK | OK | **OK** |

One NOTE is the environmental one the base build has: the HTML
manual's "Skipping checking math rendering: package 'V8' unavailable".
The vignettes rebuild OK in both packages.

**`checking tests` is OK in both, and the counts are kept.** A passing
`checking tests` echoes no counts into the check log, so the only
record of what ran inside the check is `<pkg>.Rcheck/tests/*.Rout`,
and deleting the tree destroys it. `dev/reunc-check.ps1` copies those
files out to `dev/reunc-log/check-<pkg>-*.Rout` before anything is
deleted:

| package | inside the check |
|---|---|
| frmtmb | FAIL 0, WARN 16, SKIP 154, PASS 10530 |
| frmtmb.sample | FAIL 0, WARN 15, SKIP 11, PASS 1791 |

frmtmb's 10530 is 10512 plus the 18 assertions the `propagate_error`
rename added to `test-predict-re-uncertainty.R` (section 13). The
check sets `NOT_CRAN`, so this is the whole ungated tier for the two
packages the lane touches, run in one process inside the check.

**The second NOTE: a claim of mine that three runs refuted.** The
first session got an examples-timing NOTE, `residuals.frmtmb_fit` at
13.61 s of CPU against 25.20 s elapsed, and attributed it to load
because that CPU-to-elapsed ratio is the signature
`dev/lane-rules.md` describes. The second run did not reproduce it and
this record said so: "the second NOTE was load and this run is the
control that says so". **That was too strong and it is withdrawn
here**, on the page, because the third run reproduced it with a ratio
that does not fit the load story at all.

| run | examples section | `residuals.frmtmb_fit` | CPU / elapsed |
|---|---|---|---|
| 1 | NOTE | 13.61 CPU, 25.20 elapsed | 0.54 |
| 2 | OK, 42 s for the file | did not cross | |
| 3 | NOTE | 5.55 CPU, 5.58 elapsed | 0.99 |
| 4 | OK, 47 s for the file | did not cross | |

Four runs of the same tree, two NOTEs and two clean. Run 3 crossed the
5 second threshold on a quiet CPU, by 0.58 s. So the
example SITS ON the threshold and load moves it across and back; it is
not that the NOTE is purely an artifact. The distinction matters for
whoever ships this: the example is worth trimming, rather than
explaining away each time it appears. It is `residuals.frmtmb_fit`, a
poisson fit on the scalar standard-error route, which this lane does
not touch, and the route this lane changed is the ordinal and
categorical one, which no example uses.

**The ported brms tier and the ledger.** Recorded on the lane build
(`FRMTMB_PORT_ROOT` the worktree, `FRMTMB_PORT_LIB` the lane library)
with `dev/brmsport-record.sh`, then rebuilt with
`dev/brmsport-ledger.R`, which stopped at nothing. The totals do not
move:

| outcome | count |
|---|---|
| pass | 243 |
| defect | 71 |
| divergence | 35 |
| cannot transfer | 145 |
| bin 1 total | 494 |

`dev/brmsport-ledger.tsv`, `dev/brmsport-verdicts.tsv` and
`dev/brmsport-log/ledger-summary.md` differ from the committed copies in
LINE ENDINGS only, and `git diff` reports no content change in them. One
recorded file changes content,
`dev/brmsport-log/rec-frmtmb-methods.tsv`: four assertions
(`brmsfit-methods:345`, `:350`, `:764`, `:775`) record `pass` where the
committed record says `defect` or `pending 2.6d`. Those are the four
rows `dev/shapes-findings.md` section 13 already reported as holding on
the merged build, so the committed record predates the merge rather than
this lane changing them; the ledger's own verdicts and totals are
unchanged.

The recording above is the first session's, taken before punch round
1's `re_row_support()` fix, and it was NOT repeated. What was repeated
is the tier the ledger is derived from: the 15 `test-brms-suite-*.R`
files run inside the gated tier, which the second session reran on
`reunc2-lib` for a total of 3,130 assertions, the same figure to the
unit. A verdict that had moved would have moved that total. Said here
rather than left implied, because a record carried across a code
change needs its warrant stated.

### 9a. The user library was destroyed again on 2026-09-22 at 17:06

Recorded here because `dev/machine-library.md` asks the next session to
collect this evidence before a restore overwrites it.

**What this section used to claim, and why it is withdrawn.** It said
"AFTER every tier", and that "every number in this document was
measured before it", on the grounds that the last tier log was
`check.log` at 16:29:56. That was true when it was written and is not
true now. Punch round 1's second session reran the ungated suite, the
gated tier and `R CMD check` on 2026-09-23, after the loss and after
the restore another session ran that evening, on a rebuilt user
library and a fresh private library `reunc2-lib`. The tiers in section
9 are therefore on the FAR side of the loss. Withdrawn on the page
rather than quietly corrected, because the next reader would otherwise
use this sentence to date a number.

The three reruns agree with the pre-loss run wherever they can be
compared, which is the check `dev/machine-library.md` asks for: gated
3,130 assertions to the unit, and the suite's assertion total moving
only by the guards each punch round added, which section 9
tabulates. `R CMD check`
returned the same 2 NOTEs, both of them the ones the base build has;
section 9 has all three of its runs and what they settled.

What follows is the loss evidence itself.

- **When.** 217 package directories in
  `C:/Users/adf44/AppData/Local/R/win-library/4.6` went hollow in two
  minutes: 134 with an mtime of 2026-09-22 17:06 and 83 at 17:07, out
  of 410. One more at 18:33, which is probably my own failed load.
  223 were hollow when I measured at 18:32.
- **Outside `%LOCALAPPDATA%`: nothing.** `rellib-r3` (8 of 8),
  `reunc-lib` (2 of 2) and `pinlib` all have every `DESCRIPTION`. This
  is the most discriminating fact the machine notes ask for, and it
  points AWAY from the interrupted-install theory that the fifth loss
  supported.
- **No `00LOCK` directory anywhere in the library**, which an R
  installer would have left.
- **The canary survived.** `ZZZ-canary.txt` is still there, dated
  2026-09-09.
- **Free space was 38 GB of about 951 GB, under 5 percent**, which is
  the low-disk condition that triggers `SilentCleanup`. At the
  2026-09-09 loss the figure was 48.8 GB and was called "near five
  percent". `Get-ScheduledTaskInfo` for the DiskCleanup task returns
  nothing under this account, so the trigger is not confirmed, only the
  condition.
- **What was running.** No process of this lane: its last R exited at
  16:30. Another session's `R CMD check --as-cran` wrote
  `C:/Users/adf44/source/r/correct-check` at 17:02, four minutes
  before.
- **A restore was already running by 18:34** from another session: the
  hollow count fell from 223 to 32 in three minutes. I did not restore
  anything, because the user library is read-only to a lane.

The signature (two minutes, nothing outside `%LOCALAPPDATA%`, no
`00LOCK`, low disk) is the Disk Cleanup shape rather than the
interrupted-install shape, which is the opposite of the fifth loss and
the same as the fourth.

## 10. What I did NOT do, and why

* **`simulate(re_formula = )` is untouched.** Its switch is lme4's:
  `NA` redraws the group effects and anything else keeps them, so `~0`
  and a partial formula both simulate conditionally, silently. brms has
  no `simulate()`, so there is no tiebreaker for what a partial formula
  should MEAN there (redraw the dropped terms? drop them?), and the
  answer is a decision rather than a fix. Filed in
  `dev/test-backlog.md`. `dharma_residuals()` and `pp_check()` on a fit
  reach it.
* **`frm_linpred()`'s definition is unchanged.** It honors
  `re_formula` the way the other methods do, which is the item, and it
  still returns the linear predictor conditional on the modes with the
  delta-method standard error over the joint covariance. No case was
  found for changing what it returns.
* **A population smooth's and a `gp()` curve's coefficients are NOT
  drawn** by `predict()`, in either build. They live in `b` but
  `re_formula = NA` keeps them, so drawing them at `NULL` and not at
  `NA` would make the two settings differ by more than the group
  effects. brms draws them (they are parameters of its posterior), so a
  `predict()` interval on a smooth-only model is narrower here than
  brms's by the smooth's own uncertainty. That is a pre-existing
  divergence, unchanged and now stated; it is a separate decision from
  this one.
* **`re_formula` on a structured family** (a hidden-state or
  group-level-class likelihood) is still gated by `structure_gate()`:
  the gate asks whether the caller set the argument at all, before the
  resolution turns a partial formula into a reduced design.

## 11. What needs the user

1. **`re_formula = ~(1 | nosuch)` is refused here and silently dropped
   in brms** (section 5). The brief said to match brms unless brms does
   something else, in which case to report it; brms does something
   else. The refusal is what this package does with every other
   argument that would otherwise change an answer quietly. If the user
   wants brms's behavior instead, the one-line inverse is in
   `re_resolve()`: treat an unmatched term as contributing nothing.

   **The same decision covers a second case** the review found
   (section 12, minor 1): a nested group the fit does not have.
   `~(1 | a/b)` expands to `a` and `a:b`, and on a fit with `a` alone
   brms returns `~(1 | gr(a))`, keeping the half it recognizes and
   dropping the other half silently, while frmtmb refuses the whole
   formula. One rule decides both: refuse an unmatched term, or drop
   it the way brms does.
2. **A partial `re_formula` on a fit with a factor-smooth term is
   refused** (section 4). brms would keep the smooth, because a smooth
   is not a group-level term there; frmtmb's `re_formula = NA` drops
   it. Refusing is the only answer that is not a guess, but if the user
   prefers brms's reading the rule is one line.
3. **DECIDED 2026-09-23 by the user, and done: the argument is
   `propagate_error` and it holds the group effects too.** This item
   asked whether `param_uncertainty = FALSE` should keep drawing the
   group effects, on the argument that an effect is not a parameter.
   The answer is no, and the name was the reason the question arose.
   Section 13 has the rename, what it breaks and the measurement.

## 12. Punch round 1

**Two sessions did this round.** The first applied every item and was
stopped with its tiers part run, so a second session took it over,
built its own private library `C:/Users/adf44/source/r/reunc2-lib` and
remeasured. No code changed hands part written: the first session's
last install postdates its last source edit. What the second session
found stale was the EVIDENCE, and it says so here rather than leaving
it: the cost table in M1 predated the log it quotes, the bitwise
battery predated the last edit to `R/predict-brms.R`, and the suite,
the gated tier and `R CMD check` were part run or predated the fixes.
Everything in this section and in section 9 is now the second
session's own run on `reunc2-lib`, against the same read-only base
`rellib-r3`. The two sessions used the same toolchain, **RTMB 2.0**,
TMB 1.9.25, Matrix 1.7.6, brms 2.23.0, R 4.6.1. Sections 1 to 6 and
the coverage numbers in section 3 were measured by the first session,
also on RTMB 2.0; parts of its `R CMD check` ran on RTMB 1.9 and that
check is superseded by the rerun in section 9.

MERGEABLE with two majors and five minors, all of them costs or
records rather than wrong answers. Every item is fixed below, with what
it was measured at and what it is measured at now. The statistical core
reproduced on the reviewer's own simulations, including a `many40`
design (40 groups of 3) where the two errors of the `few4` design
cannot cancel: base 0.9290 ML and 0.9305 REML against this lane's
0.9515 and 0.9535, with a median width of 3.998 against 4.383 and the
exact oracle's 4.396.

### M1. The finite-difference route was O(levels)

`fitted()` on an ordinal or categorical fit differences the group
effects, and it differenced EVERY kept level: one pair of model
evaluations per level whether or not the rows being predicted load it.
The reviewer measured, at 6 observations per group, 20 groups 0.07 to
0.27 s, 100 groups 0.18 to 3.16 s and 200 groups 0.31 to 10.92 s,
growing about as `ng^1.6`, with a dense Jacobian of
`(n K) x (p + levels)` and a dense `(p + levels)^2` covariance behind
it.

Two things fix it, and they differ in kind.

**The bound.** `re_used_b()` asks the DESIGN which levels the predicted
rows load, from the same matrices `re_eta()` and `lp_delta_A()` read:
each predictor's `Z` in sample, the rebuilt `re_design_matrix()` and
the smooth parts on new data. Every other kept level has a Jacobian
column of exactly zero, so this is exact arithmetic: dropping a zero
column changes no sum.

**The batch.** In sample every level IS loaded, so the bound cannot
help there, and that is the case that was worst. Within one block a row
loads at most one level, so perturbing every level of the block at once
changes each row by exactly its own level's perturbation, and the
difference read off that row is the one a single-coefficient
perturbation would have produced. A block therefore costs one pair of
evaluations rather than one per level. The condition is CHECKED, not
assumed (`re_b_batches()`): a multi-membership term, where a row loads
two levels of one block, and an `rr` block, whose coefficients are a
function of several `b` entries, both fall back to one coefficient at a
time. Batching also lets the quadratic form be taken over the block
structure rather than over a dense Jacobian, which is what removes
that matrix. It does NOT remove the dense `(p + levels)^2` joint
covariance, which the quadratic form needs and which 0.61.0 already
builds for the scalar route; the measurement is below.

**What still costs every level, and why.** Three cases, all of them
named in the code rather than discovered by a timer.

* **In sample.** Every level is loaded, so the bound is the whole set
  by arithmetic and cannot help. The batch is what carries this case,
  and the `in` column below is what it costs.
* **A multi-membership term and an `rr` block.** A row of an `mm()`
  term loads two levels of one block, so a batched difference cannot
  be attributed to a level; an `rr` block's coefficients are a
  function of several `b` entries, so a loaded column does not name
  one `b` position. `re_b_batches()` checks both and returns `NULL`,
  and the route falls back to one coefficient at a time over the whole
  block. These are the small blocks in practice, and the fallback is
  the 0.61.0-era cost rather than a new one.
* **A design that cannot be rebuilt on `newdata`.** `re_used_b()`
  returns `NULL`, which means "no bound", and every kept level is
  perturbed as before. Correct but slow, which is the right way round.

`dev/reunc-fdcost.R`, logs `dev/reunc-log/fdcost-lane.txt` and
`fdcost-base.txt`, both rerun by the second session on `reunc2-lib`
and `rellib-r3`, RTMB 2.0. The answer first, on a 20-level fit where
the two sets differ (2 rows load 2 of 20 levels):

| comparison | identical | max relative | max ulp |
|---|---|---|---|
| 2 newdata rows, bounded and batched, against one at a time | FALSE | 2.05e-16 | 0.9 |
| in sample, batched, against one at a time | FALSE | 2.22e-16 | 1.0 |
| two blocks, batched, against one at a time | FALSE | 4.14e-16 | 1.9 |
| a random slope, 40 coefficients in 2 batches | FALSE | 3.50e-16 | 1.6 |
| multi-membership, which refuses to batch | TRUE | 0 | 0 |

The bound alone is bitwise; the batch adds the same terms in a
different ORDER, so the two agree to at most 2 ulps rather than
exactly. That is stated rather than hidden behind an `identical()` that
would have failed.

**The categorical half, which the table above does not reach.** Every
fit there is cumulative, and a cumulative fit has ONE linear
predictor. A categorical fit has one per non-reference category, so
`re_row_support()` sums several `Z` matrices and a block spans two
predictors, which is the case the per-position batch exists for.
`dev/reunc-fdcat.R`, log `dev/reunc-log/fdcat.txt`, 25 levels, 3
categories, against the same one-at-a-time reference:

| fit | case | identical | max ulp |
|---|---|---|---|
| `y ~ x + (1 \| g)` | 2 newdata rows | TRUE | 0.0 |
| `y ~ x + (1 \| g)` | in sample | FALSE | 1.1 |
| `+ muc ~ x + (1 \| h)` | in sample | FALSE | 1.3 |
| `+ muc ~ x + (1 \| h)` | 2 newdata rows | FALSE | 0.6 |
| `+ muc ~ x + (1 \| g)`, same factor twice | in sample | FALSE | 1.1 |

The first fit keeps 50 effects for 25 levels, because the block spans
both predictors, and it batches into 2, one per column position. The
2 predicted rows load 4 of the 50. Nothing here needed a fallback.

### M1, and a wrong answer the punch round's own review found

Widening the proof past the cumulative one-block fit found a defect
the batch had introduced, and it is recorded here in full because it
is the kind of thing a performance fix is supposed to be checked for
and the round-1 proof did not reach it.

**What was wrong.** `re_b_batches()` attributes a row's difference to
a level through `re_row_support()`, an indicator of which rows each
coefficient reaches. In sample that comes from each predictor's `Z`,
which carries every block. On NEW DATA the design is rebuilt in two
pieces, ordinary group parts and smooth parts, and `re_row_support()`
read only the first. A factor smooth, `s(x, g, bs = "fs")`, lives in
the second, and `re_governed_b()` keeps it, so every one of its
columns looked like a coefficient no row loads. A batch whose rows all
say "no owner" does not refuse; it writes a derivative of exactly
zero. So `fitted(newdata = )` zeroed the smooth's Jacobian columns in
`Est.Error`, silently.

**Zeroing a column is not the same as removing its variance, and the
error is NOT one-signed.** An earlier wording here said the
contribution was "dropped" and the standard error "understated by 68
percent". That is what happened on the design below, and it is not the
property of the failure. The variance is
`Jd' Vdd Jd + 2 Jd' Vdb Jb + Jb' Vbb Jb`, so setting a `Jb` column to
zero removes a positive third term AND a cross term of either sign,
and a cell whose cross term is negative comes out LARGER. Measured on
the review's own design, one cell was 0.1020 where the reference was
0.0777, and the same effect shows on the esicar defect below, where
the pre-fix ratio to the reference runs from 0.4170 to 1.0204 on one
fit (`dev/reunc-log/fdcar.txt`). "Understated" invites the reader to
treat the failure as conservative and it is not.

**What it cost, measured.** `dev/reunc-fdsmooth.R`, 6 levels, a
cumulative fit, 3 newdata rows, against the one-at-a-time reference.
Pre-fix log `dev/reunc-log/fdsmooth-prefix.txt`, post-fix
`fdsmooth.txt`:

| fit | case | before | after |
|---|---|---|---|
| `s(x, g, bs = "fs")` alone | in sample | 1.8 ulp | 1.8 ulp |
| `s(x, g, bs = "fs")` alone | 3 newdata rows | 0 ulp | 1.6 ulp |
| `+ (1 \| h)` | in sample | 2.5 ulp | 2.5 ulp |
| `+ (1 \| h)` | 3 newdata rows | **max relative 0.676** | 0.9 ulp |

**Why the smooth-only fit was right, which is the instructive part.**
With no ordinary group term there are no `re_parts` on newdata at all,
so `re_row_support()` returned `NULL`, batching was refused for the
whole call, and the route fell back to one coefficient at a time. It
was correct by accident. It took a fit with a factor smooth AND an
ordinary group term, so that the ordinary term gave the batch
something real to attribute while the smooth gave it nothing, to make
the defect visible. A proof built only on the simpler fit would have
passed.

**The fix.** `re_row_support()` adds the smooth parts on newdata, the
same two pieces `re_used_b()` already reads. The batch's own condition
then decides as it was meant to: a row that loads two columns of one
batch refuses and falls back. The smooth-only newdata case now batches
into 15 rather than refusing, and is still right to 1.6 ulp.

`tests/testthat/test-predict-re-uncertainty.R` pins it, on the
`smooth + (1 | h)` fit, against the one-at-a-time reference and with
the tolerance written as ulps of what the run produced.

`tests/testthat/test-predict-re-uncertainty.R`
pins all five, with the tolerance written as a count of ulps of the
values the run itself produced and never as a fixed number, and it
asserts that the two sets DIFFER on the newdata case, because a test
where the bound is the whole set would assert nothing.

Then the cost. `in` and `nd` are `fitted()` as shipped; `in_full` and
`nd_full` are the same two calls differencing every kept level one
coefficient at a time, which is what this lane did before the bound
and the batch; `nd_NA` differences no group effect at all, which is
what 0.61.0 costs.

**The primary figure is a COUNT, because the clock on this machine
measures load.** Punch round 2 remeasured this table while another
lane ran its suite, and every cell moved by three to five times,
including `nd_NA`, which is 0.61.0's own route and which this lane
does not touch. Across six runs `in_full` at 1000 levels read 95.82,
102.82, 190.10, 261.46, 296.30 and 323.52 seconds on code that did not
change between them, a spread of 3.4 to 1. `dev/lane-rules.md`
says to count something load-independent where the clock cannot
settle a claim, and here there is an exact one: the NUMBER OF MODEL
EVALUATIONS the route makes. `fit_fd_se()` calls its function once for
the point, twice per outer parameter, and then twice per BATCH, or
twice per COEFFICIENT when there is no batch.

| levels | in | nd | in_full | nd_full |
|---|---|---|---|---|
| 20 | **11** | **11** | 49 | 49 |
| 100 | **11** | **11** | 209 | 209 |
| 200 | **11** | **11** | 409 | 409 |
| 1000 | **11** | **11** | 2009 | 2009 |

This is the whole of M1 in one table, and it is exact rather than
measured: the counts came back identical on every run of the script,
while the clock below moved by a factor of three across the same
runs. The shipped route costs **11 evaluations at every level
count**, in sample and on newdata: 9 for the point and the 4 outer
parameters, which is what 0.61.0 already spent, plus 2 for the one
batch this fit's single block needs. The unbounded route costs
`2 * levels + 9`, which the four rows confirm at 49, 209, 409 and
2009. So the route is O(1) in the number of levels where it was
O(levels), the constant it adds over 0.61.0 is two evaluations, and
none of that can be moved by what else is running.

The clock agrees and is reported second, from one run, with `nd_NA` as
the in-process control. Cells are NOT comparable across runs:

| levels | rows | in | nd | in_full | nd_full | nd_NA |
|---|---|---|---|---|---|---|
| 20 | 120 | 0.03 | 0.03 | 0.11 | 0.09 | 0.01 |
| 100 | 600 | 0.09 | 0.03 | 1.63 | 0.44 | 0.01 |
| 200 | 1200 | 0.17 | 0.03 | 6.08 | 0.85 | 0.01 |
| 1000 | 6000 | 1.08 | 0.05 | 190.10 | 4.64 | 0.01 |

`in_full` at 1000 levels is ONE round, not the best of three, because
it is minutes; every other cell is the best of three.

Three instrument faults are recorded, because each produced a number
first. The first version of the timer took its argument as an
EXPRESSION, and an R promise evaluates once, so rounds 2 and 3 of
"best of 3" measured nothing and every cell read 0.00 s; it takes a
function now. The second version had no in-sample unbounded arm, so it
could not be compared with the review's own numbers, which are in
sample; `in_full` is that arm. The third is the clock itself, and the
evaluation count above is what replaced it.

**The memory, and a claim this section used to overstate.** The review
named two dense objects, and they have different owners. An earlier
wording here said batching "removes the memory", which is true of one
of them and not of the other. `dev/reunc-fdmem.R`, logs
`dev/reunc-log/fdmem-lane.txt` and `fdmem-base.txt`, R's own
`gc()` high-water mark in megabytes over the call:

| levels | in | in_full | nd | scalar control | base in | base scalar |
|---|---|---|---|---|---|---|
| 200 | 278.0 | 261.7 | 175.5 | 179.6 | 212.4 | 231.1 |
| 1000 | 231.9 | 732.9 | 184.8 | 332.0 | 257.2 | 323.5 |

* **The dense Jacobian, `(n K) x (p + levels)`, was this lane's and is
  gone.** It is the only cell that moves decisively: `in_full` at 1000
  levels peaks at 732.9 Mb against 231.9 Mb for the shipped route on
  the same fit in the same process.
* **The dense `(p + levels)^2` joint covariance stays, and it is not
  this lane's.** `get_joint_cov()` inverts the joint precision, and
  0.61.0 already does that for `fitted()`'s SCALAR route on any mixed
  fit, which is the control column: a gaussian fit of the same shape
  peaks at 323.5 Mb on the BASE build at 1000 levels. The quadratic
  form needs every pairwise covariance of the levels a row's block
  perturbs, so this one is inherent to the answer rather than to the
  implementation, and at 1000 levels it is 8 Mb of the total.

**The instrument is weak and says so.** A `gc()` high-water mark
includes the fit itself, and the 200-level `in` cell (278.0) is larger
than the 1000-level one (231.9), which cannot be a property of the
route. So the table supports one conclusion, that `in_full` at 1000
levels allocates about three times what the shipped route does, and it
does not support reading any other pair of cells as a difference.

Rerun on the fixed build in punch round 2, both columns, every cell
reproduced to the decimal except that same 200-level `in` cell, which
moved from 277.9 to 278.0. A tenth of a megabyte of jitter in the one
cell the paragraph above already calls unreliable is the instrument
agreeing with its own caveat.

### B1, punch round 2. The same defect for car(type = "esicar")

The factor-smooth fix above closed one way the batch's premise fails.
The review found a third, and it is the one that matters most, because
the guard was written as a list of two special cases rather than as
the property those cases are instances of.

**The property.** Perturbing a block's levels together and reading
each row's difference off that row is right only if `expand_b()`
carries `b` to the coefficient vector POSITIONWISE, so that
`cvec[c_idx[i]]` is a function of `b[b_idx[i]]` alone. `expand_b()`
has exactly three branches: `rr`, where a level's coefficients are its
latent factors through the loadings; esicar, where `car_center(b, a)`
is `b - m[comp]`, the connected component's mean removed; and a plain
copy for everything else.

**What the guard tested instead.** `length(c_idx) != length(b_idx)`.
That is a PROXY for the rr case, and it is not the property. It
catches `rr` below full rank, where the lengths differ, and it misses
`rr` AT full rank and it misses esicar, where the lengths agree and
the map is still not the identity. Under centering a batch moves each
row by its own step minus the mean of all the steps, and
`re_used_b()`'s bound fails the same way, because every `b` of a
component reaches every coefficient of that component, so bounding
drops Jacobian columns that are not zero.

**The fix is the property, not a third name.**
`block_b_positionwise(bk)` lives in `R/covstruct.R` beside
`block_is_esicar()` and `expand_b()`, which is where it can be kept
honest, and `frame_needs_expand()` now calls it rather than repeating
the same disjunction. Both `re_b_batches()` and `re_used_b()` ask it.
A block that is not positionwise does not batch at all, and its whole
block joins the bound.

**The general form, proved rather than asserted.**
`dev/reunc-positionwise.R`, log `dev/reunc-log/positionwise.txt`,
enumerates every covstruct the registry offers and every car type the
package accepts, and compares the guard's verdict against what
`expand_b()`'s own branch conditions say, read out of the function
body:

    24 covstructs, 4 car types, 27 rows, 0 disagreements
    not positionwise: car(esicar), rr

It also asserts that `expand_b()` has exactly three assignments to
`cvec[c_idx]`, so a fourth special case added later breaks this script
instead of breaking a standard error. `frame_needs_expand()` agrees on
all 27.

**The measurement.** `dev/reunc-fdcar.R`, log `dev/reunc-log/fdcar.txt`.
Cumulative fit on a 3x3 rook lattice, 20 rows per location, against
`fit_fd_se(b_idx = gov, b_batch = NULL)`, which is one coefficient at
a time. That reference is VALIDATED first rather than assumed: driving
the same function with `b_batch = list()`, which is the code path a
refusal takes, gives a bitwise equal answer.

| fit | in sample | 3 newdata rows | batches |
|---|---|---|---|
| icar | 1.1 ulp | 0.9 ulp | 1 |
| icar + `(1 \| h)` | 2.7e-16 rel | 2.7e-16 rel | 2 |
| escar | 1.8 ulp | 0.9 ulp | 1 |
| escar + `(1 \| h)` | 4.1e-16 rel | 3.2e-16 rel | 2 |
| esicar | **0, exactly** | **0, exactly** | refused |
| esicar + `(1 \| h)` | **0, exactly** | **0, exactly** | refused |
| rr, rank 2 of 2 | **0, exactly** | | refused |
| rr, rank 1 of 2 | **0, exactly** | | refused |

esicar and rr are exact because they no longer batch, so the shipped
route IS the reference. icar and escar still batch and are still
right, which is what makes this a property of the TYPE rather than of
`car()`.

**The defect reconstructed on the fixed build**, the same technique
M2 uses: `block_b_positionwise()` is replaced in the namespace by the
old length proxy, and the swap is shown to change the BATCH and not
the MODEL, by checking `fitted_point()` across it.

| fit | old batches | old bound | in sample | 3 newdata rows |
|---|---|---|---|---|
| esicar | 1 | 3 of 9 | 0.583 rel | 0.409 rel |
| esicar + `(1 \| h)` | 2 | 4 of 13 | 0.603 rel | 0.394 rel |

The review measured 0.614 and 0.269 on its own design and 0.385 and
0.173 with the extra term; same size, different lattice and seed. Its
Monte Carlo arbitration over the joint law put the reference at 0.826
to 1.090 of the truth and the shipped route at 0.704 to 0.961, so the
reference is the right target and not merely the other number.

**The error is not one-signed here either.** Old over reference in
sample runs from 0.4170 to 1.0204 on the plain esicar fit and 0.3967
to 1.0171 with `(1 | h)`: most cells too small, at least one too
large. The cross term `2 Jd' Vdb Jb` is why, and it is the same
reason the smooth wording above is corrected.

**Not a regression in absolute terms, and that is not the point.**
0.61.0 differences no group effect at all, so on the review's design
it gives 0.0848 where this lane gives 0.1009 and the truth is 0.1202:
the lane is closer even with the defect. What was wrong is the CLAIM.
Section 12's M1 says the batch gives "the same number, not an
approximation of it" and that "the condition is CHECKED, not
assumed", and the test file pins that claim. The check was one
predicate short of it.

`tests/testthat/test-predict-re-uncertainty.R` pins the whole of it:
that the block IS esicar and the predicate calls it non-positionwise,
that `re_b_batches()` refuses and `re_used_b()` returns the whole
block, that the answer is `identical()` to the one-at-a-time reference
in sample and on newdata, and that icar on the same data still batches
and is still right. The preconditions are asserted, so a future change
that made esicar refuse for some unrelated reason could not pass this
by falling back.

**What it costs.** An esicar or rr block now differences one
coefficient at a time. Those blocks are small in practice, a lattice
of locations rather than thousands of levels, and correctness is not
tradeable against the cost. `rr` at full rank batched correctly before
and does not now; that is accepted deliberately, because a guard
written as the general property is worth more than a guard that is
right for two enumerated cases and silent about the third.

### M2. predict() spent the caller's random numbers

`predict_b_drawer()` took one seed per replicate with `sample.int()`
BEFORE `predict_simulate()` captured the `.Random.seed` it restores on
exit, so a conditional `predict()` moved the caller's stream by
`ndraws` uniforms. The predictions were identical and the stream
position was not, so a `runif()` after `predict()` generated a
different data set on the two builds. The reviewer's independent
330-quantity rebuild found 18 differing: the 9 intended plus 9 stream
positions.

The capture now happens before the drawer takes its seeds, so a
conditional `predict()` leaves `.Random.seed` exactly where 0.61.0 left
it. The battery asks about it now: `dev/reunc-bitwise.R` records the
stream after a conditional `predict()` on newdata, after
`predict(allow_new_levels = TRUE)` at unseen levels, after an in-sample
`predict(summary = FALSE)` and after `fitted()`, on every mixed model
in every mode. Against base:

    identical 330, moved as intended 9, differs 0

330 against round 1's 294, the difference being the 36 new stream
checks, every one of them identical: 15 `stream_after` (a population
call, on all five models in all three modes) and, on the three mixed
models in the three modes, 9 each of `stream_after_cond`,
`stream_after_newlev`, `stream_after_insample` and
`stream_after_fitted`. No consumption remains, so NEWS needs no
exception for it.

**The defect was reconstructed and seen failing, on the fixed build.**
It cannot be seen on `rellib-r3`, which has no drawer at all, and the
build that had it is round 1 of this lane, which no longer exists. So
`dev/reunc-stream.R` (log `dev/reunc-log/stream.txt`) rebuilds the
round-1 ordering out of the shipped function: it takes
`predict_simulate()`'s body apart, moves the `.Random.seed` capture and
its `on.exit()` back to AFTER the drawer, and installs that in the
namespace, so the defect and the fix run in one process on one build.
Conditional `predict()` against `re_formula = NA`, 25 draws, gaussian
`(1 | g)` on 10 groups of 6:

| ordering | same stream position | next `runif()` |
|---|---|---|
| shipped | TRUE | 0.5671238950 both |
| round 1, reconstructed | FALSE | 0.5671238950 against 0.1416320675 |

and the predictions of the two orderings at one seed are
`identical()`, which is the reviewer's point exactly: the answer was
right and the caller's next random number was not.

`tests/testthat/test-predict-re-uncertainty.R` pins it, comparing the
stream position and the next uniform after a conditional `predict()`
against the same after `re_formula = NA`, in sample and at an unseen
level. That test PASSES on `rellib-r3`, which is not evidence of
anything and is said here for that reason: the base build has no
drawer to consume, so the test is a guard against reintroducing the
defect and the reconstruction above is what has been seen to fail.

### The five minors

Four of the five are claims about behavior, and the second session
rechecked all four on its own build rather than reading the first
session's note: `dev/reunc-minors.R`, log `dev/reunc-log/minors.txt`,
`reunc2-lib`, RTMB 2.0. Every line below that says what the build does
is from that log.

1. **A second undisclosed departure from brms**, the same root as the
   refusal in section 5 and recorded beside it in section 11. On a fit
   with `a` but not `a:b`, `brms:::check_re_formula(~(1 | a/b), ...)`
   returns `~(1 | gr(a))`: brms keeps the part it recognizes and drops
   the rest silently. Measured here too: `~(1 | a:b)` returns `~1`, and
   `~(1 | a) + (1 | a:b)` returns `~(1 | gr(a))`. frmtmb refuses the
   whole formula, because one expanded term matches nothing. The user
   decides this together with the `nosuch` refusal. Rechecked on this
   build: `~(1 | nosuch)`, `~(1 | g/h)` and `~(1 | g:h)` on a
   `(1 | g) + (1 | h)` fit are all three refused, classed and naming
   the term.
2. **`re_formula = ~1` means opposite things in `predict()` and
   `simulate()`.** Verified rather than taken on trust: at one seed,
   `simulate(re_formula = ~1)`, `~0` and `~ (1 | g)` are all
   `identical()` to `NULL` and only `NA` differs, while
   `predict(re_formula = ~1)` is `identical()` to `re_formula = NA`.
   The NEWS bullet names `simulate()` as the exception now, says the
   two readings differ, and points at `dharma_residuals()` and
   `pp_check()`, which reach it. Rechecked on this build on a
   `(1 | g) + (1 | h)` fit: `simulate()` at `~1`, `~0` and `~ (1 | g)`
   is `identical()` to `NULL` and only `NA` differs, while
   `predict(~1)` and `fitted(~1)` are `identical()` to `NA` and
   `predict(~1)` is NOT `identical()` to `NULL`.
3. **The b-perturbed fits shared the original's cache.** `fp <- fit`
   copies the list and not the environment, so `joint_precision()` on a
   perturbed fit returned the memo. Latent, because nothing on that
   path reads the cache today, and fixed the way `fit_set_outer()`
   already did it: a fresh `new.env(parent = emptyenv())` per
   perturbation. Rechecked on this build with the memo POPULATED
   first, which is the case the defect needs: after
   `joint_precision(fit)` the fit's cache is non-empty, and the
   perturbed copy's cache is a different environment and is empty, as
   is `fit_set_outer()`'s.
4. **Doc overreach, corrected in both places.** Section 2.1 and
   `predict_b_drawer()`'s comment said the draw is Henderson's PEV
   "with the variance parameters at their draw". It is not: `Q` is read
   once, at the estimates, and a replicate's parameters shift the
   conditional MEAN through `d - d_hat`, not the conditional variance.
   Both say so now, and section 2.1 connects it to the term section 3
   measures as still missing from every interval here.
5. **A quadrature fit is silent on the SCALAR route.** Filed with a
   reproduction: on `y ~ x + (1 | g)` with `bernoulli()` and
   `quadrature = TRUE`, `fitted()` reports `Est.Error` 0.09407 and
   0.08913 at a known level against 0.09340 and 0.09095 at
   `re_formula = NA`, so the group-effect term is missing and nothing
   is said; that fit's joint covariance carries `beta` and `theta` and
   no `b`. A Laplace fit of the same data reports 0.15115 against
   0.09158, which is the size of the term. Pre-existing, 0.61.0 does
   the same, and more visible now that the finite-difference route
   warns. Filed in `dev/test-backlog.md`. Rechecked on this build, to
   four more digits and on both rows: quadrature 0.09407 and 0.08913
   at a known level against 0.09340 and 0.09095 at `re_formula = NA`,
   Laplace 0.15115 and 0.14998 against 0.09158 and 0.08906, and
   `fitted()` on the quadrature fit raises no warning.

## 13. predict(propagate_error = ), a rename and a decision

The user's decision of 2026-09-23. `param_uncertainty` is now
`propagate_error`, default `TRUE`, and it switches off MORE than the
old flag did: `propagate_error = FALSE` holds the parameters AND the
group effects at their estimates. The old name is removed outright
rather than deprecated, under the break-compatibility rule: there are
no users yet, so an alias would only carry a wrong name forward.

**Why the name.** It states the AXIS. The old one named a subset of
the things being held, which is what made section 11 item 3 a question
at all: if the argument is about "parameters", then a group effect is
arguably not one, and the flag should leave it alone. If it is about
whether the error in the ESTIMATES is propagated, a fitted group
effect is plainly an estimate and the answer falls out. It does not
freeze the observation noise or an unseen level's fresh draw, because
neither is an estimate.

**Two arguments, two questions, neither able to state the other.**
This is the part the documentation now carries, because the pair is
easy to confuse.

* `re_formula` chooses WHICH terms are in the prediction. Only
  `re_formula = NA` can say "predict for an average group".
* `propagate_error` chooses whether the error in the estimates is
  carried into the interval. Only `propagate_error = FALSE` can say
  "include this group's own effect, but treat it as known", which is
  0.61.0's conditional-on-the-modes behavior.

Measured rather than asserted. `dev/reunc-propagate.R`, log
`dev/reunc-log/propagate.txt`: a gaussian `(1 | g)` fit, 8 groups of
5, predicted at the level whose fitted effect is largest in absolute
value (1.622), 4000 draws at one seed, all four corners.

| `re_formula` | `propagate_error` | estimate | interval width |
|---|---|---|---|
| `NULL` | `TRUE` | 2.4417 | 2.9297 |
| `NULL` | `FALSE` | 2.4356 | 2.5616 |
| `NA` | `TRUE` | 0.8093 | 3.1097 |
| `NA` | `FALSE` | 0.8138 | 2.5616 |

Read the two `FALSE` rows together: the same width, because both carry
the observation noise and nothing else, and centres apart by that
group's own fitted effect. Read the two `NULL` rows together: the same
centre, and a width that grows by the estimation error. The first is
`re_formula` moving the centre where `propagate_error` cannot; the
second is `propagate_error` moving the width where `re_formula`
cannot. The noise-only claim has its own check: the `FALSE` width is
0.9731 of `2 * 1.96 * sigma` at the fit's own `sigma`, the shortfall
being the Monte Carlo quantile at 4000 draws.

**How equal, in a document that spends pages on that question.** Four
printed decimals cannot carry either comparison, and this section used
to say "identical width" and "the fitted effect to four decimals". The
first overstated and the second understated. Both corners have a
mechanism, and it is the same one: with `propagate_error = FALSE` and
one seed the two calls simulate from the SAME stream, and their linear
predictors differ by the constant `b_g`. So the draws differ by
exactly that constant.

* The **centres** therefore differ by `b_g` as an identity, not as an
  average: 1.6217512329993795 against a fitted effect of
  1.6217512329993797, **0.62 ulp apart**, equal to 16 significant
  figures rather than to four decimals.
* The **widths** are a difference of two quantiles of the same set
  shifted by a constant, so they agree up to the rounding of that
  subtraction. Here they came out bitwise equal, `identical()` TRUE at
  2.5615970090632789, **0.00 ulp**. On the review's own design they
  landed 0.56 ulp apart. Both are the same statement and neither is an
  exact identity of construction, so the claim is "equal to the
  rounding of a quantile subtraction", not "identical".

`tests/testthat/test-predict-re-uncertainty.R` asserts the centre
identity draw by draw rather than on the mean, as a ratio to the fit's
own `sigma`.

brms needs no such argument, because a posterior draw carries its own
parameters and its own group effects, so there is nothing to switch
off. `?predict.frmtmb_fit` says so, in one sentence, so that a brms
user does not go looking for the missing argument.

**What it cost elsewhere, and the one thing that got harder.** There
is no public call left that draws `b` while holding the parameters.
That combination was how the lane isolated the analytic conditional
law in section 2.2, so `dev/reunc-joint.R`, `dev/reunc-plugvar.R` and
`tests/testthat/test-predict-re-uncertainty.R` now drive
`predict_b_drawer()` directly, which is what the old flag did
internally. Section 2.2 records the numbers before and after, and they
agree within Monte Carlo error.

**Two coverage scripts, renamed, now mean opposite things, and each
says so at its call site.** `dev/reunc-coverage.R`'s plug-in arm no
longer reproduces the rows it produced, because the flag it passes
switches off more than it used to; a comment there says the recorded
rows are the OLD semantics and that re-running produces a different
quantity rather than a check on them. `dev/shapes-coverage.R`'s
identical-looking arm DOES still reproduce, because that script
predates the group-effect draw: when it ran there was no `b` draw to
switch off, so the two readings of the flag coincided there. Both
facts are now comments in the scripts, not only in this document,
because the next person to run one of them will be reading the script.

**Two dev scripts were pinned to a stale library and this found them.**
`dev/reunc-bdraw.R` and `dev/reunc-plugvar.R` hardcoded the first
session's private library instead of taking one, so the rename showed
up there as "unused argument (TRUE)" against a build from a day
earlier rather than as a result. Both take the library as an argument
now, with the same default the other scripts use. Worth recording
because a pinned library is the quiet version of the stale-log problem
this round has already paid for twice: the script runs, it prints
numbers, and the numbers are from the wrong build.

**What was rerun.** The bitwise battery, which must not move for a
rename: `identical 330, moved as intended 9, differs 0`, and that
includes `predict_na_plug`, the one battery entry that passes the flag.
`test-predict-re-uncertainty.R` 65 of 65; `test-brms-shapes-punch2.R`
22, which is the other test that passes it; and
`test-brms-shapes.R` 72, `test-brms-shapes-punch.R` 55,
`test-predict-newdata.R` 12, `test-predict-lp-basis.R` 30,
`test-re-formula-partial.R` 39, `test-arg-refusal.R` 123,
`test-api-spellings.R` 36, `test-methods.R` 64,
`test-conditions.R` 150 and `test-backlog.R` 35, all passing.

`R CMD check --as-cran` was rerun too, at 15:25:10, after every edit:
frmtmb 1 NOTE, frmtmb.sample OK. Its `checking tests` runs the whole
ungated tier for both packages in one process with `NOT_CRAN` set, and
came back FAIL 0, PASS 10530 for frmtmb and FAIL 0, PASS 1791 for
frmtmb.sample. So the whole of both packages did run against the
shipped tree, and no separate suite or gated tier was launched. The
reasoning for not launching them: only two test files named the
argument and both were rerun; only two functions gained or changed a
parameter and no other caller exists; the default is unchanged, so a
file that never names the flag cannot see the change; and the bitwise
battery is the direct evidence that the shared path did not move.
`test-arg-refusal.R`, which is where a removed argument would surface,
passes 123. Section 9.0 records which tier covers which build, and
what the two pre-rename tiers do and do not stand for.

**The one gap in that reasoning, and who closed it.** `R CMD check` on
frmtmb and frmtmb.sample does not reach the other six extensions, so
"no extension names the argument" was a grep and not a run. The review
ran all six extension suites against the base build and found them
identical, which turns the grep into a measurement. `param_uncertainty`
appears nowhere under `extensions/` or `vignettes/` in any of the
seven packages.

**A packaging defect the rename surfaced.** testthat writes extracted
reproductions to `tests/testthat/_problems/` when a file fails under
the `check` reporter, and this lane's seen-failing runs against the
base build left 26 of them. The directory is in `.gitignore`, so it
never showed in `git status` and nothing flagged it. It was NOT in
`.Rbuildignore`, and a tarball built from the working tree carried all
27 entries, two of them calling `param_uncertainty`: a removed
argument, shipped in a source package, in files that would never run
because testthat does not recurse into subdirectories. `R CMD check`
passed over them for the same reason, so the check was not going to
find this.

Fixed by deleting the directory and adding `^tests/testthat/_problems$`
to `.Rbuildignore`. Verified in BOTH directions rather than by reading
the rule, the way the correctness lane verified its own `Rplots` rule:

| state | `_problems` entries in the tarball |
|---|---|
| before, 26 scratch files present | 27 |
| after, directory recreated with a canary file | 0 |

The second row is the one that matters. Deleting the files alone would
make the count zero without the rule doing anything, so the rule was
tested with the thing it excludes PRESENT, which is this project's
standing requirement for a guard. The 159 real test files still ship,
and `param_uncertainty` now appears nowhere in the tarball except
`NEWS.md`, where it is the historical entry for 0.61.0.

**Names rejected, recorded so nobody reopens it.** `plug_in`, jargon,
and its negative default reads backwards. `incl_uncertainty`, which
reads as "which things are included" rather than naming the axis.
`propagate_err`, because none of core's 191 documented argument names
contains a clipped word.
