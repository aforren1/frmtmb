# Lane slopescale: autoscale and random-slope columns

Round of 2026-11-13, base ec8178d3 (frmtmb 0.82.0), branch
`wt-r83-slopescale`, library `C:/Users/adf44/source/r/wt-slopescale-lib`.
Item 1018 of `dev/test-backlog.md` (silent-wrong). The "before" arm is
rellib-r21 (0.82.0) throughout. Logs are local and gitignored; every
number below is pasted from the log named beside it.

## Summary

- Autoscale planned a random-effect (Z) column only when a fixed-effect
  (X) column with the same values was planned too
  (`autoscale_plan_z()` walked the planned X columns and matched Z
  against them), and only in `us`, `diag`, `us_t` and `diag_t` blocks.
  So `y ~ 1 + (0 + x | g)` and `y ~ 1 + (1 + x | g)` were never
  rescaled, and neither was a slope in `cs()`, `ar1()`, `mm()` and the
  rest, with or without a fixed twin.
- Now every random-effect column with more than two distinct values is
  planned on its own, scaled by its root mean square, in every structure
  where that is an exact reparameterization. The fits of the triage go
  from 42 of 72 short (21 silent, up to 111.2 units) to 0 of 72; 13
  structure cases x 6 seeds go from 76 of 78 short to 0 of 78.
- The wrong boundary verdict ("sd_g__x is at zero" for a stall at sd 1)
  came from the SE check's upward probe, whose grid is on the response's
  scale and does not reach a slope's maximum. Fixed with a small hunk in
  `R/se-check.R`: a gain at the probe's far point counts as a gain.
  16 of 48 such verdicts before, 0 after; those fits now warn that they
  stopped short.
- A second form of the same problem is filed, not fixed: under
  `autoscale = FALSE` a slope that does move can leave its correlation
  run to the limit at a short stop, and the fit calls it a boundary.

## Fix round (review `dev/reviews/2026-11-13-slopescale.md`: FIX)

### F1 (regression): the default engaged on the root mean square

The first pass engaged the default when a Z column's root mean square
was below 0.05. For an uncentered column, x = c0 + s z, that value is
near c0, so `y ~ x + (1 + x | g)` with x = 5 + 0.03 z no longer engaged
where 0.82.0 did (through the fixed twin's sd). Now
`autoscale_block_cols()` also returns `spread`, the sd over the rows
where the column is nonzero, each plan entry carries it, and
`autoscale_decide()` engages on `spread < autoscale_slope_sd`. The scale
stays the root mean square.

- New test "an uncentered slope column of small sd engages the
  default": the reviewer's three cases (c0 5 seeds 2 and 4, c0 7 seed 4,
  s 0.03), against the fit on the centered, standardized column.
  On the first-pass build (`dev/slopescale-tests/autoscale-fix-before.log`):
  `c0 5 seed 2 ... 9.837156`, `c0 5 seed 4 ... 6.051766`,
  `c0 7 seed 4 ... 6.091467` short, and `fd$autoscaled` FALSE on all
  three (`RESULT test-autoscale.R pass=175 fail=7 err=0 skip=0 warn=1`;
  the 7th fail is the F2 assertion below).
  After (`dev/slopescale-tests/autoscale-fix-after.log`):
  `RESULT test-autoscale.R pass=182 fail=0 err=0 skip=0 warn=0`.
- The reviewer's 80-fit table rerun on the fixed build
  (`dev/reviews/slopescale-review-uncentered.R`, log
  `dev/slopescale-fix-uncentered.log`; "short" is against the best of
  base default, lane default and lane TRUE):

  ```
  review (first pass):
  fits: 80 ; short > 0.001 : base 58 (silent 34 ) ; lane 44 (silent 22 ) ; lane TRUE 0
  lane short where base was not: 3 ; base short where lane is not: 17 ; max lane short 41.1 ; max TRUE short 5.94e-09
  fix round:
  fits: 80 ; short > 0.001 : base 58 (silent 34 ) ; lane 0 (silent 0 ) ; lane TRUE 0
  lane short where base was not: 0 ; base short where lane is not: 58 ; max lane short 5.94e-09 ; max TRUE short 5.94e-09
  ```

  No fit is short on the lane default; 0.82.0 is short on 58 (34
  silent). Where 0.82.0 engaged (the twin form) it was short too,
  because it scaled the Z column by the sd (0.03), which makes an
  uncentered column of mean 5 about 167 wide.
- Ordinary fits, bit for bit against 0.82.0 in one process
  (`dev/reviews/slopescale-review-ordinary.R`, log
  `dev/slopescale-fix-ordinary.log`): `base deterministic: 21 of 21 ;
  lane identical to base (all fields): 21 of 21`.

### F1b (pre-existing, not changed): poorly centered slope columns

x = 2000 + 3 z and x = 5 + 0.1 z have sd above 0.05, so the default
does not engage, on 0.82.0 or now. Filed HIGH in the backlog with the
reviewer's numbers (10 of 10 short for the year, 7 silent, 3.3 to 41.1
units; 8 of 10 short for 5 + 0.1 z, 7 silent, 3.1 to 15.3;
`autoscale = TRUE` within 1.1e-8 on all 30). The trigger is a user
decision; this round does not change it.

### F2 (LOW): the default-mode note fired on fits that are not short

Dropped: under the default, `autoscale_note_unscaled()` names a slope
only when the X plan rescales its column, as 0.82.0 did. The test
"names every structure" now asserts that the default does NOT name
`homdiag(1 + x | g)` on a small column (it failed on the first-pass
build, as above). I did not look for a fit in that case that stops
short; the reviewer's 4 such designs all reach the maximum.

### Docs

- `autoscale_block_exact()`'s roxygen says it checks at one point and
  names the structure that passes falsely.
- The callers of the `se_sd_gain_up0()` hunk are named above.
- `frmtmb_control(autoscale =)` and NEWS: the default engages on the
  column's standard deviation; the note under the default is as in
  0.82.0.

### Reruns after the fix: see "Acceptance checks, fix round" below.

## What changed

### `R/autoscale.R` (owned by this lane)

- `autoscale_plan()`: the X part is unchanged. The Z part is planned by
  the new `autoscale_plan_z(frame)` over the blocks, independent of X.
  A linear predictor can now carry `z` entries and no `cols` (its `par`
  is then NULL); `autoscale_frame()`, `autoscale_map()` and
  `autoscale_units()` skip the X part when `cols` is empty.
- `autoscale_block_cols()`: per block position, the Z columns over the
  levels, the matching `b` entries, and the root mean square of the
  column over the rows where it is nonzero; `NA` where the column has 2
  or fewer distinct values (intercepts, dummies) or is a `mo()`
  component's.
  - Root mean square, not sd, as the SCALE: a random slope is not
    centered, so its effect on the predictor is sd times the column
    itself. For a centered column the two agree; for `x + 5` with
    spread 1e-4 the sd would make the column worse.
  - The sd over the same rows (`spread`) is what ENGAGES the default
    (fix round, F1 below). The first pass engaged on the root mean
    square, which was a regression.
- `autoscale_plan_z()`: a block is rescaled when its coefficients are
  its parameters positionwise (`block_b_positionwise()`), it is not a
  smooth, gp, hsgp, rr, latres or equalto block, and:
  - a structure with one log sd per column: each qualifying column by
    its own factor;
  - a structure with one log sd: all columns by one factor (geometric
    mean of the spreads), only when every column qualifies and all are
    in one linear predictor;
  - and the structure is one of `autoscale_exact_core` (us, diag, us_t,
    diag_t, cs, toep, hetar1, gr_cov, gr_prec; homdiag, ar1, homcs,
    homtoep, ou, exp, gau, mat), or, for a structure another package
    registers, its density passes `autoscale_block_exact()`.
  - A theta entry mapped to a constant keeps its block as it is; a
    non-centered block (`frame$ncp_blocks`) keeps its `b`.
- `autoscale_block_exact(bk, cpos)`: checks that
  `nll(b * cpos, theta + log cpos) == nll(b, theta) - n_levels *
  sum(log cpos)` at ONE generic theta and b, built without the RNG, to
  1e-8 relative. It is a spot check, not a proof: the review built a
  structure equivariant near its start and not near where the slope
  lands, which passes; a false pass costs only a worse start
  (`dev/reviews/slopescale-review-blockexact.log`, 6 of 6 seeds within
  1.3e-9). The roxygen now says so. The test file runs it on every structure of
  `autoscale_exact_core`, so a change of parameterization fails there.
- `b` and `theta` move by `bcols` (b index) rather than by the Z column
  index. The old code indexed `b` by the Z column, which is right only
  while no rr block sits before the block.
- `autoscale_note_unscaled(frame, plan, control)`: names the slopes that
  are not rescaled (one-sd structure over mixed columns, `mo()`, rr, a
  registered structure that fails the check). It is now called whenever
  `autoscale` is not FALSE, also where nothing is planned. Under the
  default it names a slope only when the X plan rescales its column, as
  0.82.0 did; under TRUE it names every one. (The first pass also named,
  under the default, a slope whose column was below 0.05; dropped in
  the fix round, F2.) New wording: "autoscale did not rescale the random
  slopes of ...".

### `R/fit.R` (shared; two hunks)

- The call site of `autoscale_note_unscaled()` in `fit_assembled()`
  (3 lines): called when `!objective_only && is.null(template)`, with
  `control`. The `cold_start` rule is not touched.
- The roxygen of `frmtmb_control(autoscale =)`: which columns and
  structures are rescaled, and the default's rule. `man/frmtmb_control.Rd`
  regenerated.

### `R/se-check.R` (selimit2 owns the file; one hunk of 9 lines)

- In `se_sd_gain_up0()`, after the far point (`se_sd_far`, +10 on the
  log scale) is evaluated: if it gains more than every grid point, it
  counts as the gain and as the best point (`at`). Nothing else in the
  function changed.
- `se_sd_gain_up()` has three callers, so the hunk reaches all three:
  `se_sd_at_zero()` and `se_boundary_names()` in R/se-check.R (the
  verdict), `sd_up_escape()` (R/fit.R:4432, which restarts from `at`
  when the sd is below 1% of its term's scale and the gain is above
  `grad_tol`), and `headroom_drop_rows()` (R/fit.R:6177, which adds the
  gain to the convergence headroom). The review saw no change through
  the last two: 260 singular fits and 21 ordinary fits identical to
  0.82.0 (`dev/reviews/slopescale-review-singular.log`,
  `-ordinary.log`).

### Tests: `tests/testthat/test-autoscale.R`

- "autoscale names every structure whose slope it leaves" rewritten: a
  registered `mydiag` now passes the density check and is rescaled; a
  registered `mybump` (variance `exp(2 theta) + 0.5`, which does not
  scale) is named; `homdiag(1 + x | g)` is named under TRUE and under
  the default; absent cases `(1 + x || g)`, `mydiag`, `homdiag(0 + x | g)`
  and `autoscale = FALSE`.
- New: "a random slope with no fixed effect on its column is rescaled"
  (gaussian and poisson, both forms, seed 2, spread 1e-4: logLik within
  1e-8 relative, sds, ranef/coef, fitted, predict draws seed for seed).
- New: "a random slope is rescaled in every covariance structure"
  (cs, ar1, homdiag, gr(cor = FALSE), mm, dpar, |ID|, nested; skip on
  CRAN).
- New: "each structure autoscale trusts has a density that scales"
  (all 17 core structures, the inverse for one-sd structures with
  unequal factors, and the registered pair, with `.Random.seed`
  unchanged).
- New: "a slope column of spread 1 leaves the fit as autoscale = FALSE"
  (guard absent: plan holds the slope, decide returns NULL, `opt$par`
  identical).
- New: "a slope stopped at its start is not called at zero".

## Evidence

### Step 1: reproduction on rellib-r21, and why the column was skipped

`dev/slopescale-repro.R` (the triage's `slope-unscaled.R` with the
library arm as an argument and no time cap; seeds 1 to 6; spreads 1e-3,
1e-4, 1e-5; forms `(0 + x | g)` and `(1 + x | g)`):

```
dev/slopescale-repro-gau-before.log:
short > 1e-6: 24 of 36 ; short > 1e-6 with no warning: 12 ; boundary message: 14 ; NA: 0 ; max short 12.67278
dev/slopescale-repro-pois-before.log:
short > 1e-6: 18 of 36 ; short > 1e-6 with no warning: 9 ; boundary message: 9 ; NA: 0 ; max short 111.2185
dev/slopescale-repro-gau-after.log:
short > 1e-6: 0 of 36 ; short > 1e-6 with no warning: 0 ; boundary message: 0 ; NA: 0 ; max short 2.861611e-09 ; max |short| 2.861611e-09
dev/slopescale-repro-pois-after.log:
short > 1e-6: 0 of 36 ; short > 1e-6 with no warning: 0 ; boundary message: 0 ; NA: 0 ; max short 8.727739e-10 ; max |short| 7.791043e-07
```

The before counts are the triage's (`slope-gau.log`, `slope-pois.log`)
exactly. The 7.79e-7 after is negative short: the spread-1 reference of
poisson `(0 + x | g)` seed 2 stops 7.79e-7 below the rescaled fit (lme4
and glmmTMB at spread 1 stop at the same point, below).

Cause: `autoscale_plan()` skipped every linear predictor without a
planned X column (`if (!length(cols)) next`), and `autoscale_plan_z()`
took a Z column only by matching its values against a planned X column.
`(0 + x | g)` has no X column but the intercept, so the plan was NULL
and the default never engaged.

### Step 2: every covariance structure

`dev/slopescale-structs.R`, 13 cases x seeds 1 to 6, gaussian, slope sd
0.6, spread 1e-4 against spread 1 on the same data. Per case, `short`
counts logLik(spread 1) - logLik(spread 1e-4) > 1e-6; `silent` those
with no warning; `boundary` a `frmtmb_boundary_fit` message on the 1e-4
fit:

```
before (dev/slopescale-structs-before.log)     after (dev/slopescale-structs-after.log)
ar1      n=6 short=6 silent=0 boundary=1 max=65.79   short=0 boundary=0 max=5.754e-10
cs       n=6 short=6 silent=0 boundary=0 max=45.76   short=0 boundary=0 max=2.507e-11
dbar     n=6 short=6 silent=0 boundary=0 max=50.05   short=0 boundary=0 max=1.057e-11
dpar     n=6 short=6 silent=0 boundary=0 max=29.77   short=0 boundary=0 max=4.32e-10
grcov    n=6 short=5 silent=1 boundary=2 max=51.57   short=0 boundary=0 max=5.605e-11
grnocor  n=6 short=6 silent=0 boundary=0 max=50.05   short=0 boundary=0 max=1.057e-11
homdiag  n=6 short=6 silent=0 boundary=0 max=50.05   short=0 boundary=0 max=1.057e-11
idlink   n=6 short=5 silent=0 boundary=0 max=32.99   short=0 boundary=0 max=6.191e-09
mm       n=6 short=6 silent=4 boundary=6 max=5.427   short=0 boundary=0 max=2.603e-11
nested   n=6 short=6 silent=3 boundary=3 max=50.16   short=0 boundary=6 max=8.151e-08
student  n=6 short=6 silent=2 boundary=2 max=34.72   short=0 boundary=0 max=1.439e-09
us       n=6 short=6 silent=3 boundary=3 max=33.42   short=0 boundary=0 max=3.857e-10
us0      n=6 short=6 silent=0 boundary=0 max=48.88   short=0 boundary=0 max=1.589e-10
before: short > 1e-6: 76 of 78 ; with no warning: 13 ; max short 65.78902
after:  short > 1e-6: 0 of 78 ; with no warning: 0 ; max short 8.150852e-08
```

Cases: us `(1 + x | g)`, us0 `(0 + x | g)`, dbar `(1 + x || g)`, cs
`cs(1 + x | g)`, ar1 `(1 | g) + ar1(0 + tf:x | g)`, homdiag
`(1 | g) + homdiag(0 + x | g)`, grnocor `gr(g, cor = FALSE)`, grcov
`gr(g, cov = A)`, student `gr(g, dist = "student")`, mm `mm(g, h)`,
nested `(1 + x | a/g)`, dpar `sigma ~ 1 + (0 + x | g)`, idlink
`(1 + x | p | g)` with `sigma ~ (1 | p | g)`.

Back-transformation, after, the largest over the 72 non-nested fits:
VarCorr sds (slope's times s) 3.61e-5 relative; ranef() and coef()
(slope's times s) 1.07e-5 of the largest effect; fitted values 1.33e-5
and predict() draws with one seed 3.46e-5 of sd(y). Standard errors
(rerun with them, `dev/slopescale-structs-after2.log`, same 78 fits):

```
short > 1e-6: 0 of 78 ; with no warning: 0 ; errors: 0 ; max short 8.150852e-08
non-nested max: vc 3.6e-05 re 1.07e-05 fit 1.33e-05 pred 3.46e-05 se_fe 2.29e-05 se_vc 0.000117 se_nan mismatches 0
```

`se_fe` is the largest relative difference of a fixed effect's
standard error, `se_vc` of a VarCorr sd's Est.Error (the slope's times
s); no sd has a NaN error on one scale and not the other. The reported fit is
the unscaled model warm-started from the pre-fit, so these are the
original fit's own accessors; no accessor reads the plan.

The nested case: the outer factor `a` has 5 levels and its correlation
is at its limit on both scales in all 6 seeds (the boundary message is
the spread-1 fit's too; seen by hand on seed 1). logLik agrees within
2.0e-7; on seed 3 the `a` sds differ by 0.30 relative and its effects
by 0.45 at equal logLik, a flat ridge of a singular fit, not a scaling
defect. The test compares only the nested logLik.

### Step 3: the boundary verdict

`dev/slopescale-verdict.R`, `autoscale = FALSE`, gaussian and poisson,
both forms, spreads 1e-4 and 1e-5, seeds 1 to 6:

```
dev/slopescale-verdict-before.log: fits: 48 ; 'sd_g__x is at zero': 16 ; with no warning: 27
dev/slopescale-verdict-after.log:  fits: 48 ; 'sd_g__x is at zero': 0 ; with no warning: 11
```

Why: `se_sd_at_zero()` asks `se_sd_gain_up()` for a gain above the
estimate on a grid of 0.01 to 1 times the term's scale (the residual
sd, 1.17 here). At spread 1e-5 the slope's maximum is near sd 5e4, so
the grid finds nothing (only `log 1.17` is above the stop at log sd
5.9e-8). The far point, log sd + 10, did raise the log-likelihood, but
it entered only `moved` ("the likelihood reads the sd"), so the rule
read "no gain, and read": at zero. Now a gain at the far point counts.
The 16 fits now warn "the fit stopped short of the maximum: this
standard deviation is near zero ... a larger value raises the
log-likelihood, so it is not at a boundary". At spread 1e-4 the grid's
top point already gave the warning before the change.

The remaining 11 silent fits are 6 poisson `(1 + x | g)` at 1e-4,
which reach the maximum (repro: short <= 8e-9), and 5 that end with a
boundary message about the correlation, filed below. All 11 are
`autoscale = FALSE` only; under the default the slope is rescaled.

### Step 4: lme4 and glmmTMB

`dev/slopescale-lme4.R`: lme4 2.1.0 (`lmer(REML = FALSE)`, `glmer`)
and glmmTMB 1.1.15.2 (`REML = FALSE`) on the same 48 data sets of step
3, each at spread s and at spread 1. Each column is the shortfall
against the best of the six fits of that data set:

```
before (dev/slopescale-lme4-before.log)
  frm_s : 42 of 48 short > 1e-6 (max 111, NA 0)
  frm_1 : 0 of 48 short > 1e-6 (max 1.08e-08, NA 0)
  lme4_s : 14 of 48 short > 1e-6 (max 16.7, NA 0)
  lme4_1 : 0 of 48 short > 1e-6 (max 1.76e-08, NA 0)
  tmb_s : 4 of 48 short > 1e-6 (max 5.24, NA 37)
  tmb_1 : 0 of 48 short > 1e-6 (max 9.87e-11, NA 0)
after (dev/slopescale-lme4-after.log)
  frm_s : 0 of 48 short > 1e-6 (max 2.86e-09, NA 0)
  frm_1 : 0 of 48 short > 1e-6 (max 7.79e-07, NA 0)
  lme4_s : 14 of 48 short > 1e-6 (max 16.7, NA 0)
  lme4_1 : 0 of 48 short > 1e-6 (max 7.81e-07, NA 0)
  tmb_s : 4 of 48 short > 1e-6 (max 5.24, NA 37)
  tmb_1 : 0 of 48 short > 1e-6 (max 7.79e-07, NA 0)
```

So frmtmb on the small-spread data now reaches the best optimum within
2.86e-9 on all 48; on spread 1 within 7.79e-7, the one data set
(poisson `(0 + x | g)` seed 2) where all three packages stop 7.8e-7
below frmtmb's rescaled fit. lme4 on the small spread is short on 14
(up to 16.7) and glmmTMB errors on 37 and is short on 4: the same data
in other units, which these two do not rescale for a slope either.

### Step 5: cost

`dev/slopescale-cost.R`. Every time is a block grown past 1.2 s, the
minimum of 5 rounds (3 for a whole fit in (a); the 1e5-row fit once),
because proc.time ticks at 10 ms. Two runs, the second with
block-timed fits in (a), while 14 lanes shared the machine:

```
dev/slopescale-cost.log (first run; (a) fit times were min of 3 single runs)
  sleepstudy fit 0.010  plan_z 8.42e-05  note 6.77e-05  (plan_z+note)/fit 1.52e-02
  pois5k     fit 0.300  plan_z 3.37e-04  note 3.44e-04  (plan_z+note)/fit 2.27e-03
  gau1e5     fit 17.970  plan_z 9.61e-03  note 8.36e-03  (plan_z+note)/fit 1.00e-03
  sleepstudy default 0.0522  autoscale=FALSE 0.0756  spread 1 0.0325  default/spread1 1.61
  pois5k     default 0.5400  autoscale=FALSE 0.6950  spread 1 0.3275  default/spread1 1.65
dev/slopescale-cost2.log
  sleepstudy fit 0.0800  plan_z 2.83e-04  note 1.71e-04  (plan_z+note)/fit 5.68e-03
  pois5k     fit 0.6650  plan_z 8.01e-04  note 7.86e-04  (plan_z+note)/fit 2.39e-03
  gau1e5     fit 34.1500  plan_z 1.46e-02  note 1.55e-02  (plan_z+note)/fit 8.83e-04
  sleepstudy default 0.1163  autoscale=FALSE 0.1025  spread 1 0.0550  default/spread1 2.11
  pois5k     default 0.8800  autoscale=FALSE 1.2300  spread 1 0.6400  default/spread1 1.38
```

The absolute times doubled between the runs (machine load); the ratios
are the result. The first run's sleepstudy fit time is one tick, so its
ratio is not usable. Every fit pays the planning: 0.09 to 0.6 % of the
fit (sleepstudy `Reaction ~ Days + (Days | Subject)`; a poisson GLMM,
5000 rows, `(1 + x | g)`; a gaussian, 1e5 rows, 1000 groups, three
correlated terms). A fit where the default engages (a slope column
spread 1e-4) pays the pre-fit: 1.4 to 2.1 times the fit of the same
data at spread 1, about the cost of the stalled `autoscale = FALSE` fit
(0.7 to 1.1 times it).

### Acceptance checks

- Tests seen failing on rellib-r21, the final test file
  (`dev/slopescale-tests/autoscale-before-final.log`):
  `RESULT test-autoscale.R pass=61 fail=68 err=1 skip=0 warn=23`, by
  test: 41 FAIL "rescaled in every covariance structure", 21 FAIL "no
  fixed effect on its column", 3 FAIL "not called at zero", 2 FAIL
  "names every structure", 1 FAIL "spread 1 leaves the fit" (the plan
  had no z entry), 1 ERROR "each structure ... scales"
  (`autoscale_exact_core` absent: the weak form; the others are the
  behavioral failures).
- After: `dev/slopescale-tests/autoscale-after.log`, `RESULT
  test-autoscale.R pass=176 fail=0 err=0 skip=0 warn=0`.
- Guard absent (sd(x) 1): test "a slope column of spread 1 leaves the
  fit as autoscale = FALSE" (`opt$par` identical), and 20 of 20 ulp
  draws identical in `dev/slopescale-ulp.R`.
- Ulp robustness (`dev/slopescale-ulp.R`, y or x times
  `1 + 2e-15 k`, k = 1 to 20; each number is the largest ratio of the
  asserted quantity to its bound):

  ```
  dev/slopescale-ulp-ref.log:
  slope cases: 80 ; worst ratio to bound: ll 0.137 sd 0.0232 eff 0.0158 fitted 0.0157 ; not autoscaled 0 ; with warnings 0
  guard identical: 20 of 20 ; verdict short not zero: 20 of 20
  ```

- OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
  lapack`; Rblas 74 routines to OpenBLAS, Rlapack 616):
  `dev/slopescale-ulp-ob0326.log` the same line as above;
  `dev/slopescale-tests/autoscale-ob0326.log` `RESULT test-autoscale.R
  pass=176 fail=0 err=0 skip=0 warn=0`;
  `dev/slopescale-repro-gau-ob0326.log` `short > 1e-6: 0 of 36 ... max
  short 2.861441e-09`. 0.3.32 was not run: no assertion here sits at a
  rounding-decided branch (the default's thresholds are three orders
  of magnitude from the test spreads).
- All nine suites, one file per process, `NOT_CRAN=true`,
  `FRMTMB_BRMS_FIT_TESTS=true`, lane core ahead of rellib-r21 (no
  extension code changed, so the extensions are rellib-r21's, each
  loading this core; the frmtmb path of every `lib:` line, 529 of 529,
  is the lane library): `dev/release/slopescale-suite-all.log`,
  `SLOPESCALE-SUITE-ALL ran 529 of 529`, pass 37557, fail 0, err 0,
  skip 29, warn 0. The 29 skips are opt-in tiers: drmTMB agreement
  (13), the scale tier (15), the fuzz tier (1). Run apart:
  `FRMTMB_DRMTMB_FIT_TESTS=true` test-drmtmb-agreement.R pass=131
  fail=0 err=0 skip=0 warn=0; `FRMTMB_FUZZ=true` test-fuzz.R pass=2
  fail=0 err=0 skip=0 warn=0. The scale tier was not run (cost).
- One session per suite, three orders (`dev/release/run-sessions.sh`,
  N = 2, same environment): `dev/release/slopescale-sessions.log`,
  `SLOPESCALE-SESSIONS ran 27 of 27`. Every session fail 0 and warn 0;
  every one err 0 except frmtmb.learn preload, err 1 warn 1, in
  test-stan-identity.R: "TRE pattern compilation error 'Out of memory'"
  in `sub()` on a Stan DSO name, while 14 lanes shared the machine.
  Rerun alone (`dev/slopescale-tests/learn-preload-rerun.log`):
  `SESSION frmtmb.learn preload files=16/16 pass=545 fail=0 err=0
  skip=2 warn=0`. 194 LEAK lines, every one tagged `[with load of
  ...]`; none untagged. Core: forward pass=27568, reverse 27733,
  preload 27836 (skip 15, 15, 17): the counts depend on the order on
  this build; whether they do on 0.82.0 too was not measured.
- R CMD check --as-cran of core (`dev/slopescale-check.sh`, built with
  vignettes inside `dev/slopescale-check/`, since deleted; log kept as
  `dev/slopescale-check.log`; `R_LIBS` = lane library, rellib-r21, user
  library; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; pandoc and TinyTeX on
  PATH): `Status: 2 NOTEs`, tests `[100m] OK`. The two NOTEs are the
  expected ones: examples over 5 s (`frm_hazard_reads` 2.83 s user,
  29.10 s elapsed; `dharma_residuals` 1.19 s user, 6.39 s elapsed;
  elapsed far above CPU is machine load, dev/lane-rules.md) and the
  HTML manual's skipped math rendering (V8 unavailable). No extension
  was checked: no extension file changed.

### Acceptance checks, fix round

The checks above are the first pass. After the F1 and F2 changes:

- test-autoscale.R, reference BLAS: `pass=182 fail=0 err=0 skip=0
  warn=0` (`dev/slopescale-tests/autoscale-fix-after.log`).
- test-se-limit.R, reference BLAS: `pass=240 fail=0 err=0 skip=0
  warn=0` (`dev/slopescale-tests/selimit-ref.log`).
- OpenBLAS 0.3.32 with LAPACK (`dev/ciharden-openblas.sh 0.3.32
  lapack`: Rblas 76 routines to OpenBLAS, Rlapack 617;
  `OPENBLAS_CORETYPE=Haswell`, 4 threads; the copy is deleted):
  test-autoscale.R `pass=182 fail=0 err=0 skip=0 warn=0`, test-se-limit.R
  `pass=240 fail=0 err=0 skip=0 warn=0`
  (`dev/slopescale-tests/autoscale-ob0332.log`, `selimit-ob0332.log`).
  Whether this Windows OpenBLAS build honors `OPENBLAS_CORETYPE` is
  not established (dev/ci-findings.md).
- Ordinary fits bit for bit against 0.82.0: 21 of 21 (above).
- Core suite, one file per process (`NOT_CRAN=true`,
  `FRMTMB_BRMS_FIT_TESTS=true`, 2 at a time):
  `dev/release/slopescale-suite-core2.log`, `SLOPESCALE-SUITE-CORE2 ran
  330 of 330`, pass 27574, fail 0, err 0, skip 15 (the opt-in tiers, as
  before), warn 0; frmtmb from the lane library in 330 of 330 `lib:`
  lines. The extension suites and session runs were not repeated: the
  fix changes only when the default engages on a slope column and the
  default note, and the reviewer's 21 ordinary fits are bit-identical
  to 0.82.0.
- R CMD check --as-cran of core, same script and settings
  (`dev/slopescale-check-fix.log`; the check tree is deleted):
  `Status: 1 NOTE`, the HTML manual's skipped math rendering (V8
  unavailable). Examples `[65s] OK`, tests `[46m] OK`. The examples
  timing NOTE of the first pass did not recur (machine load).
- Process cap: 3 R processes at most until the machine-wide notice,
  then 2. Every process count in the fix round stayed within the cap in
  force (the check ran alone after the notice). The first pass exceeded
  3 once, briefly (R CMD check beside two scripts).

## Not done, and why

- The 0.05 threshold (`autoscale_slope_sd`) is unchanged. The backlog
  item's last sentence (at spread 0.06 the plain fit fell short by at
  most 2.2e-4 over 64 fits) was calibrated with a fixed twin; it was
  not remeasured for a slope without one.
- The large side was not measured for a random slope (a column spread
  1e3 or more); the default does not engage there, as before.
- `se_sd_scale()` was not made column-aware. Scaling the upward probe's
  grid by the column's spread would also let `sd_up_escape()` restart a
  stall under `autoscale = FALSE`, but it moves the grid of every
  random slope with spread far from 1, in selimit2's file. The far-point
  hunk is the smallest change that removes the wrong verdict.
- The message of the "short" verdict says "this standard deviation is
  near zero" for a stop at sd 1, which is near zero only against the
  slope's scale. Left to selimit2's wording.

## Defects found and not fixed (filed in `dev/test-backlog.md`)

- HIGH (fix round, F1b, the reviewer's finding): the default never
  engages on a poorly centered slope column (a calendar year, x = 5 +
  0.1 z), which stalls silently; 0.82.0 the same. Numbers in the fix
  round section above. Needs a user decision on the trigger.
- MEDIUM: under `autoscale = FALSE`, a random slope on a column spread
  1e-4 or 1e-5 whose sd does move can leave the correlation parameter
  run to the limit at a short stop (`(1 + x | g)`, gaussian, seed 4 at
  1e-4: theta = (-0.588, 6.649, 265.6), 6.08 units short; seed 5:
  3.90; seeds 1 and 2 at 1e-5: 2.54 and 4.51), and the fit says
  "Boundary (singular) fit: cor_g__Intercept__x is at the limit of its
  range" with no warning (`dev/slopescale-scratch.R`, last version).
  `se_edge_cor()` asks only whether moving further changes nothing; a
  move back toward 0 would show the gain.
- LOW, pre-existing, not investigated: `y ~ x + z + (1 + x + z | g)`,
  1000 groups of 100, gaussian, seed 2 (`dev/slopescale-gau1e5.R`)
  warns "Large maximum absolute gradient at the optimum (0.00372 at
  theta_6) ... may not have converged" on 0.82.0 and on this build
  alike (identical `opt$par`); `theta_3` (the sd of `z`, true value 0)
  is at log -6.51. Whether the warning is wrong was not checked
  against lme4.

## Files touched

- `R/autoscale.R`, `R/fit.R` (call site, roxygen), `R/se-check.R` (one
  hunk), `man/frmtmb_control.Rd` (generated), `NEWS.md`,
  `tests/testthat/test-autoscale.R`, `dev/test-backlog.md` (CLOSED line
  under 1018, three filed items: HIGH, MEDIUM, LOW).
- New in `dev/`: `slopescale-findings.md`, `slopescale-install.sh`,
  `slopescale-repro.R`, `slopescale-structs.R`, `slopescale-verdict.R`,
  `slopescale-lme4.R`, `slopescale-cost.R`, `slopescale-ulp.R`,
  `slopescale-gau1e5.R`, `slopescale-scratch.R`, `slopescale-check.sh`;
  `slopescale-suite.lst` (a job list, gitignored as `*.lst`).

## Release

NEWS has a BREAKING bullet (autoscale now changes fits it did not
change before) and a behavior bullet (the verdict). A minor bump; no
new floor (no new dependency, no extension change).
