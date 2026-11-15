# Short check of lane slopescale's fix round (round of 2026-11-13, base ec8178d3)

Checks only the fixes of `dev/reviews/2026-11-13-slopescale.md` (F1,
F1b, F2, docs), against the "Fix round" section of
`dev/slopescale-findings.md`. Lane library `wt-slopescale-lib`, base
`rellib-r21`. Scripts and logs: `dev/reviews/slopescale-short-*`. The
reviewer's own scripts were rerun unchanged, with output to
`slopescale-short-*.log`. At most 2 R processes at a time. Free space
on C: 146,582,106,112 bytes at start, 146,343,698,432 at end.

## Verdict: SHIP

Every fix holds. One LOW observation on the new trigger (below, 1c)
belongs with the filed F1b item. It does not block.

## 1. F1: the default engages on the spread, the scale stays RMS

Code (`R/autoscale.R`):

- `autoscale_block_cols()` returns `rms = sqrt(mean(nz^2))` and
  `spread = stats::sd(nz)` over the nonzero rows. It returns both only
  when the column has more than 2 distinct values and `rms > 0`.
- `autoscale_plan_z()` keeps `scale = cpos[pos]` (from `rms`) and adds
  `spread = cl$spread` to each entry.
- `autoscale_decide()`: `isTRUE(zz$spread < autoscale_slope_sd)`. This
  is a plain comparison with `5e-2`, with no rounding or tolerance
  branch. 0.82.0 compared `zz$scale` (the fixed twin's sd).
- The installed lane library reads `spread` in `autoscale_decide()`
  (`slopescale-short-trigger.log` line 1: `lane decide reads spread:
  TRUE`).

### 1a. The reviewer's 80 uncentered fits

`slopescale-review-uncentered.R`, log `slopescale-short-uncentered.log`:

```
fits: 80 ; short > 0.001 : base 58 (silent 34 ) ; lane 0 (silent 0 ) ; lane TRUE 0
lane short where base was not: 0 ; base short where lane is not: 58 ; max lane short 5.94e-09 ; max TRUE short 5.94e-09
```

The three former regressions, `y ~ x + (1 + x | g)`, lane short: c5
s0.03 seed 2 2.61e-10; c5 s0.03 seed 4 0; c7 s0.03 seed 4 0. The
counts match the findings exactly.

### 1b. The reviewer's 21 ordinary fits, bit for bit

`slopescale-review-ordinary.R`, log `slopescale-short-ordinary.log`:

```
base deterministic: 21 of 21 ; lane identical to base (all fields): 21 of 21
```

### 1c. Attack on the trigger

`slopescale-short-trigger.R`, log `slopescale-short-trigger.log`.
Gaussian, 20 groups of 15, seeds 1 to 5, `y ~ x + (1 + x | g)` and
`y ~ 1 + (1 + x | g)`. Arms: base default, lane default, lane FALSE,
lane TRUE, in one process. "short" is against the best of the 4 arms.

| design | distinct | spread (nonzero rows) | lane engaged | base engaged | max short lane | max short base | max short FALSE |
|---|---|---|---|---|---|---|---|
| dose01: x 0, or 1 + 0.01 z | 151 | 0.0093 to 0.0105 | 10 of 10 | 0 | 4.37e-08 | 1.41e-06 | 1.41e-06 |
| twonz: x 0, 1 or 1.01 | 3 | 0.00495 to 0.00501 | 10 of 10 | 0 | 7.03e-10 | 2.44e-09 | 2.44e-09 |
| ulpind: x 0, 0.3 or 0.1 * 3 | 3 | 3.6e-17 to 3.9e-17 | 10 of 10 | 0 | 4.95e-07 | 1.3e-06 | 1.3e-06 |
| indconst: x 0 or 3 | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| mostzero: 90 % 0, else N(0, 1) | 24 to 37 | 1 to 1.15 | 0 | 0 | 1.12e-10 | 1.12e-10 | 1.12e-10 |

```
fits per arm: 50 ; lane errors: 0 ; lane default engaged: 30 ; base engaged: 0 ; engaged on lane not base: 30
short > 0.001 : base 0 ; lane default 0 ; lane FALSE 0 ; lane TRUE 0 ; lane short where base not: 0 ; max lane short 4.95e-07
not engaged on lane: identical to base (par and conditions) 20 of 20
```

- No division by zero and no engagement at sd 0. A column whose
  nonzero rows are all equal has at most 2 distinct values, so it gets
  `rms = spread = NA` (indconst: not engaged, bit-identical to base).
  With more than 2 distinct values, the nonzero rows hold at least 2
  distinct values, so `spread > 0`.
- A column with exactly two distinct values (zero or not) never
  qualifies. Two distinct NONZERO values plus zero (twonz) qualify and
  engage.
- A mostly-zero column with spread near 1 is not engaged and is
  identical to base.
- LOW observation. The default engages on 30 of 30 well-scaled designs
  where 0.82.0 did not and where `autoscale = FALSE` reaches the
  maximum (dose01, twonz, ulpind: in effect an indicator times a
  near-constant). The result is never worse: every lane fit is within
  4.95e-7 of the best, with 0 warnings. Base gives 2 convergence
  warnings (ulpind seed 5, both forms, "singular convergence (7)",
  judged at the maximum). The cost is the pre-fit
  (`slopescale-short-cost.log`, dose01 seeds 1 and 2, both forms):
  the default takes 1.59 to 1.68 times as long as FALSE (0.052 to
  0.064 s against 0.032 to 0.038 s). In ulpind, engagement turns on
  ulp noise (0.1 * 3 against 0.3), a rounding-decided branch. Both
  sides give the same answer within 1.3e-6, and 0.82.0 already had
  this branch for a fixed-effect column (sd 1e-17 below 1e-3). This is
  cost only, from the reviewer's own recommended rule, not a defect.
  It is relevant to the filed F1b item: the ratio of RMS to sd over
  nonzero rows is about 100 on dose01, so a ratio trigger also fires
  there. Its calibration should include such a design.

## 2. F1b: filed HIGH, trigger unchanged for it

`dev/test-backlog.md`, "Filed in the round of 2026-11-13", "From lane
slopescale": "HIGH. The default never engages on a poorly centered
random-slope column, which stalls silently". It carries the reviewer's
numbers (year: 10 of 10 short, 7 silent, 3.3 to 41.1; 5 + 0.1 z: 8 of
10, 7 silent, 3.1 to 15.3; TRUE within 1.1e-8 on 30) and calls the
trigger a user decision. The threshold is still
`autoscale_slope_sd <- 5e-2` on the spread. No ratio rule was added.

## 3. F2: the default-mode note arm is removed

`autoscale_note_unscaled()` names a slope when
`asked || cl$cnm %in% cn`, where `cn` holds the planned X column names
of the decided plan. Under the default with nothing engaged, the plan
is NULL, so nothing is named. Under the default with an engaged plan,
a slope is named only when its fixed twin is rescaled, the condition of
0.82.0 (`comp$cnms[k] %in% cn` over `plan`). The test "autoscale names
every structure whose slope it leaves" now asserts
`expect_false(any(grepl("random slopes of", msgs)))` for
`y ~ 1 + homdiag(1 + x | g)` under the default (test-autoscale.R lines
295 to 300). The findings record it failing on the first-pass build.

## 4. Docs

- `autoscale_block_exact()` roxygen: "checked at ONE point" and "One
  point is not a proof", with the structure of review item 3.
- Findings, "R/se-check.R" section: names `sd_up_escape()`
  (R/fit.R:4432) and `headroom_drop_rows()` (R/fit.R:6177) as the two
  further callers.
- NEWS: the top `# frmtmb (development version)` section has a
  BREAKING bullet, and no lane names, `dev/` paths or version numbers
  (searched above line 38, the `# frmtmb 0.82.0` heading).
- Nit, not blocking: NEWS and the `frmtmb_control(autoscale =)` roxygen
  say "the column's standard deviation". The code uses the sd over the
  nonzero rows, which differs for a column with zeros (dose01: 0.01
  against 0.50). "the standard deviation of its nonzero values" is
  exact.

## 5. Test files on the lane library

`slopescale-short-testfile.R`, `NOT_CRAN=true`,
`FRMTMB_BRMS_FIT_TESTS=true`, one file per process. The `lib:` line is
the lane library on all three:

```
RESULT test-autoscale.R pass= 182 fail= 0 err= 0 skip= 0 warn= 0
RESULT test-se-limit.R pass= 240 fail= 0 err= 0 skip= 0 warn= 0
RESULT test-se-check.R pass= 134 fail= 0 err= 0 skip= 0 warn= 0
```

The first two match the findings' fix-round counts.
