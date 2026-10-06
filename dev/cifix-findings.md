# Lane cifix: the Ubuntu failures of 0.68.0

Base: d24f7b86 (frmtmb 0.68.0, frmtmb.spline 0.10.0), built in
`C:/Users/adf44/source/r/rellib-r6`. Worktree `frmtmb-wt-cifix`, private
library `wt-cifix-lib`. Logs are in `dev/cifix-log/` (gitignored); every
number below is copied from them, with the script that made it.

## What failed on CI

Run 37502518163 and its siblings on d24f7b86:

| Workflow | Job | Result |
|---|---|---|
| check-frmtmb.spline | ubuntu-latest, release | FAIL 2, WARN 0, SKIP 3, PASS 565 |
| R-CMD-check | ubuntu-latest, release and devel | FAIL 2, WARN 2, SKIP 198, PASS 14664 |
| R-CMD-check | windows-latest, macos-latest | success |
| test-coverage | ubuntu | the same 2 failures |

The task named only the spline job. The core R-CMD-check and
test-coverage jobs fail too, on `test-se-check.R:383` and `:384`
("separated data are named as separation"), and the Ubuntu core run lets
2 warnings escape (`test-aliased-grouping.R:138`,
`test-open-issues.R:52`). This lane fixes all of them.

The runners install OpenBLAS 0.3.26 (pthread) as BLAS and LAPACK
(`dev/cifix-log/ci-full.log`, setup-r step). Windows R uses the
reference BLAS.

## The instrument: OpenBLAS 0.3.26 inside a copy of R on Windows

`dev/cifix-openblas.sh` copies R 4.6.1 to `dev/cifix-out/Rob` and puts
OpenBLAS 0.3.26 (the official x64 release, the runners' version) behind
its `Rblas.dll`. R's `Rblas.dll` exports `dgemmtr_` and `zgemmtr_`,
which 0.3.26 lacks, so a forwarding `Rblas.dll` sends those two to the
reference BLAS and every other symbol to `libopenblas.dll`. LAPACK stays
R's reference build. A 2000 x 2000 product takes 0.05 s against 4.4 s
with the reference BLAS (`dev/cifix-blascheck.R`). Runs use
`OPENBLAS_NUM_THREADS=4`, the runners' core count.

It reproduces every CI symptom on rellib-r6, with the test files of
d24f7b86 (`dev/cifix-log/diff-base-ob4.log`,
`t-base-ob-sechk.log`, `t-base-ob-alias.log`, `t-base-ob-openiss.log`,
`probe354.log`):

- `test-difference.R:354`: the standard error spread is 6.245e-17 = 9
  ulps of 0.0618, row 1 apart and the other 88 rows equal to each
  other. CI shows exactly that shape: row 1 0.0618114539778425939,
  rows 2 to n 0.0618114539778426564, 9 ulps.
- `test-se-check.R:383` and `:384` fail.
- `test-aliased-grouping.R` and `test-open-issues.R` each let one
  "Standard errors are not available" warning escape (warn=1; 0 with
  the reference BLAS).

It also reproduces failure 2's mechanism, below, though not on the CI
fit itself. It is not the runner: four other files fail or warn under
it on rellib-r6 and do not on CI (see "Seen with the emulator only").

## Failure 1: `test-difference.R:354`, a 9-ulp spread

Cause. 0.68.0 reads a difference curve from ONE `frm_lp_basis()` call
on the stacked grid `rbind(newdata, contrast)`, where 0.67.0 made two
calls. A position at row i of `newdata` sits at row n + i of the stack.
`pred_design()` kriged every unseen row by `backsolve()` over the
columns of `t(Ks)`, and OpenBLAS's triangular solve rounds a column by
where it sits in the right-hand side; the reference BLAS solves every
column with the same loop. So the two copies of one position got
kriging weights and variances a few ulps apart, `A1 - A2` was not
exactly zero in the gp() columns, and `gp_krig_cov()`'s exact
cancellation (it sets the covariance of two rows at one position to
the first row's variance, and then overwrote the diagonal with each
row's own) left the difference of the two variances behind.

Decision: a guarantee, not an accident. `gp_krig_cov()` already claims
that two rows at one position "cancel exactly rather than to a
rounding", and the claim costs nothing to make true.

Fix, `R/predict.R`:

- `pred_design()` kriges each distinct unseen position once (by
  `pos_rowkey()`, the key the observed-position match already uses) and
  copies the row of `Ks`, `P`, `Xw` and `rs` to every row at it.
- `gp_krig_cov()` forms the covariance over the distinct positions and
  expands it, so rows at one position have identical rows and columns
  and covary by exactly their variance. The duplicate-group loop and
  the diagonal overwrite go.

Measured (`dev/cifix-probe354.R`, `sp_gp_fit()`'s data, seed 21, n = 90):

| Build | BLAS | SE spread | rows off row 1 | estimate spread / (eps max abs eta) |
|---|---|---|---|---|
| rellib-r6 | reference | 0 | 0 | 0.6677 |
| rellib-r6 | OpenBLAS | 6.245e-17 (9 ulps) | 88 | 5.175 |
| fix | reference | 0 | 0 | 0.6677 |
| fix | OpenBLAS | 0 | 0 | 1.002 |

The test keeps `expect_identical()` on the standard error, with its
comment now naming the mechanism. The estimate's bound was 16 eps of
the difference's mean; it is now 16 eps of the largest `|eta|`
subtracted (1.498), because cancellation rounding scales with the terms,
not with their difference. Measured 0.67 and 1.0 of that unit.

New core test, `test-gp-multidim.R` "rows at one unseen position krige
alike wherever they sit": one fit, a grid of 119 unseen positions twice
with one extra row between them; design rows, `extra_var`, and the
extra covariance blocks must be identical. rellib-r6 with OpenBLAS: 4
failures (`t-base-ob-gpmd.log`, pass=58 fail=4). It passes on
rellib-r6 with the reference BLAS, which is why it exists: it pins the
property where it can fail.

The same class, in the spline: `test-difference.R:130`, "a grid
differenced with itself is exactly zero", failed with OpenBLAS at 4
threads (row 25 estimate -2.2e-16; `diff-base-ob4.log`) and passed at
the default thread count and on CI. It is the eta of identical rows at
rows i and n + i. Fix, `extensions/frmtmb.spline/R/curve-cov.R`: the
stacked difference evaluates each row that both grids share exactly
once (`sp_row_key()`, which compares doubles with `match()` because
`duplicated()` on a data frame compares at 15 significant digits) and
reads it twice. A self-difference is then exactly zero on every BLAS,
and so is the extra covariance of a row against itself.

## Failure 2: `test-difference.R:490`, "disagrees ... by 1 relative"

The task's hypothesis was a disagreement about which rows are lost.
That is not the mechanism. On this fit `se_nonest` is all FALSE on both
routes: the fit has random effects, so `get_joint_cov()` inverts the
joint precision and `jc$null` is NULL (`dev/cifix-log/diag2-base.log`).

Cause, in two parts.

1. The spline formed a standard error as `sqrt(max(q + e, 0))` from the
   coefficient variance `q = diag(A V A')` and the kriging variance
   `e`; `frm_linpred(se.fit = TRUE)` forms `sqrt(max(q, 0) + e)`. They
   differ only when `q < 0`, and when `q < -e` the spline gives 0 and
   core gives `sqrt(e)`: `|0 / sqrt(e) - 1|` is exactly 1, the number
   on CI.
2. `q < 0` happens because the joint covariance is not one. Both
   smoothing sds run to zero and lose their standard errors ("flat"),
   so the outer Hessian is indefinite, and so is the joint precision
   TMB builds on it: smallest eigenvalue -9.9e-08 on the test's seed-8
   fit, -5.0e-09 on seed 21 (`dev/cifix-diag3.R`,
   `diag3-ref.log`). `solve(Q)` then has eigenvalues of -6.7e+07 to
   1.4e+08 over its theta rows (`diag2-base.log`), and the coefficient
   block it hands to predictions is wrong too. 0.68.0's standard-error
   repair (lane nanse) reaches predictions only on a fit without random
   effects; its warning says "a prediction that moves along these
   directions gets none", which was not true with random effects.

Reproduced on Windows by construction. `dev/cifix-scan2.R` fits the
test's model at n = 100, 120, 140 and seeds 1 to 40 and calls
`frm_curve(contrast = )` (`scan2-*.tsv`):

| Build | BLAS | fits | lost an SE | some grid row q < 0 | refused |
|---|---|---|---|---|---|
| rellib-r6 | reference | 120 | 64 | 3 | 3 |
| rellib-r6 | OpenBLAS | 120 | 59 | 4 | 4 |
| fix | reference | 120 | 66 | 0 | 0 |
| fix | OpenBLAS | 120 | 59 | 0 | 0 |

Every refusal was "by 1 relative", and every refusal had a row with
`q < 0`. The refusing fits on rellib-r6 with the reference BLAS are
n = 120 seed 4 (on `contrast`), n = 100 seed 21 and n = 140 seed 23 (on
`newdata`). The CI fit is seed 8, n = 120, which does not refuse on
Windows with either BLAS; the fit itself differs by platform (theta_1
-11.34 here, -11.03 with OpenBLAS), and the scan shows the refusal on
both. The smallest grid-row q in each arm: -0.1041 with OpenBLAS (n =
120, seed 14; the "-0.10" in NEWS and in `joint_cov_repair()`'s
comment) and -0.01558 with the reference BLAS (n = 140, seed 23). The
"lost an SE" counts 64 and 66 on one BLAS differ because some fits
differ between processes on this machine (see "Fits that differ
between processes", below), not because of the patch.

How wrong 0.68.0's covariance is, on the A grid (`dev/cifix-qrange.R`,
`qrange-ref.log`, reference BLAS):

| n, seed | lost | q, 0.68.0 | q, fix | se.fit, 0.68.0 | se.fit, fix |
|---|---|---|---|---|---|
| 120, 4 | theta_2 | 0.004563 to 0.02293 | 0.005043 to 0.02514 | 0.06755 to 0.1514 | 0.07102 to 0.1586 |
| 100, 21 | theta_1, theta_2 | -0.002577 to 0.009168 | 0.004114 to 0.01857 | 0.0007579 to 0.09575 | 0.06415 to 0.1363 |
| 140, 23 | theta_1, theta_2 | 0.0004262 to 0.02126 | 0.00358 to 0.02763 | 0.02066 to 0.1458 | 0.05984 to 0.1662 |
| 120, 8 | theta_1, theta_2 | 0.004754 to 0.02127 | 0.004758 to 0.02556 | 0.06896 to 0.1459 | 0.06899 to 0.1599 |

On seed 21, `frm_linpred(se.fit = TRUE)` reported 0.00076 on a row
whose standard error is 0.064 or more. The "fix" columns equal the
covariance conditional on the lost parameters computed independently
by `dev/cifix-diag3.R` (removing their rows of `Q` before the solve),
range for range.

Fix, core (`R/predict.R`, `joint_cov_repair()`, called by
`get_joint_cov()`): with random effects and a lost standard error,
the joint covariance is built from the identity TMB's joint precision
encodes,

    V = [Qrr^-1 + J Vf J'   J Vf]     J = -Qrr^-1 Qrf
        [Vf J'              Vf  ]

with `Vf` the outer covariance the rescue already keeps
(`cov_fixed_prop`: the pseudo-inverse over the determined directions,
zero along the lost ones) in place of `H^-1`. `Qrr` and `Qrf` are the
inner Hessian's blocks and do not depend on `H`. So the random block is
the inner covariance plus the outer uncertainty propagated along the
determined directions, which is what the no-random-effect branch
already did. It is conditional on the lost DIRECTIONS (corrected after
review, m3). When each lost parameter is a direction of its own, as on
the fits here, that is the covariance with the lost parameters held,
mgcv's convention for a smoothing parameter. When one lost direction
mixes several parameters, a combination orthogonal to it keeps its
variance: on the review's `gr(g, by = f)` fixture (seed 11, one
direction over theta_2 and theta_3) a slope effect's variance is
7.06e-04 against 1.72e-06 with both held, and the prediction standard
errors differ by at most 0.4 percent (`dev/cifixrev-r1.R`). `jc$null`
gets the
lost directions at the outer rows, so a prediction that loads a lost
parameter directly gets NaN through the existing `jc_nonest()`; a
linear predictor loads none. `frm_joint_cov()` marks the lost rows NaN
from a position list `get_joint_cov()` stores (`lost_pos`); it used to
match outer names against row numbers, which is right only without
random effects, the only case where it ran. A fit that lost nothing
takes the old `solve(Q)` path, and its joint covariance is unchanged
bit for bit.

The kriging change is bit-for-bit neutral for a grid without repeated
positions, except under `gp(x, by = <numeric>)`: `gp_krig_cov()` now
scales by the multipliers after symmetrizing rather than before, which
can move entries of its diagonal blocks by an ulp.

Fix, spline (`R/curve-cov.R`): one arithmetic, `sp_se()`, clamps
`q` and `e` separately, core's order. It serves the curve, the
difference, `frm_curve_deriv()` and both standard errors of
`frm_curve_feature()`. The spline's clamp order would have kept
refusing on any fit where core's covariance went negative; with the
core fix, `q` is a variance again and the order matters only at
rounding.

Lost rows, one source. `sp_lost_rows()` reads both of core's marks
from the one `frm_lp_basis()` call: `se_nonest` (NaN) and `nonest`, a
rank-deficient design's unestimable row (NA). The spline ignored
`nonest`, so a rank-deficient fit with one such grid row was refused
whole with "disagrees ... by NA relative" (`dev/cifix-nonest.R`:
rellib-r6 refuses; the fix returns the curve with NA at row 4, equal
to `frm_linpred()` elsewhere). `sp_cov_check()` now compares WHICH
rows have no standard error before it compares values, with its own
message. Every route (single grid, stacked difference, derivative,
feature, simultaneous band) goes through `sp_curve_parts()`, so all
read the rows from there; the derivative and the feature carry both
marks through their stencils, NaN for a lost direction (with the
existing warning) and NA for an unestimable row, as core reports them.

Tests that fail on rellib-r6 on Windows with the reference BLAS:

- spline `test-difference.R` "smoothing sds at zero leave the check
  and core agreeing" (n = 100 seed 21, n = 120 seed 4): error, "by 1
  relative" (`t-base-ref-diff2.log`, pass=61 error=1).
- spline `test-curve.R` "a row a rank-deficient fit cannot estimate is
  NA, not a refusal": error, "by NA relative" (`t-base-ref-curve.log`,
  pass=78 error=1).
- core `test-se-check.R` "a fit with random effects and a lost sd gets
  a covariance": `min(q / q_in)` is -0.3 against a bound of 1 - 1e-8,
  where `q_in` is the variance with every outer parameter held (the law
  of total variance); `se^2` is 6.28e-07, the kriging variance alone,
  where `q + e` is -0.000503 (`t-base-ref-sechk.log`, fail=3 error=1;
  the last two are the weak form, `lost_pos` not existing).
- core `test-se-lost-re-predict.R` (added after review, m4): on the
  OLRE fit of `test-se-check.R` (`y ~ x + (1 | id)`, seed 5, sigma and
  the id sd lost), 0.68.0 gives sigma predictions a standard error of
  `0 0 0` with no warning under ML and REML; the patch gives NaN and one
  warning, and mu keeps its standard errors.

On rellib-r6 with OpenBLAS, the CI platform, two of these behavioral
pins do not fail, because the fits differ by platform (review, m7): the
core law-of-total-variance test fails only in its weak `lost_pos` part,
and the spline "smoothing sds at zero" test passes. The CI fit (seed 8)
and the property assertions still catch a regression there. The CI fit
keeps its test unchanged.

## The core failures CI showed and the task did not name

`test-se-check.R:383`. The separation fit runs its estimates toward
infinity, so where nlminb stops is rounding. At the default budget
(`eval.max = iter.max = 1000`) it stopped at code 0 after 1939
objective and gradient calls on Windows and hit the limit, code 9, with
OpenBLAS after 1997; the convergence warning then explains the fit and
the SE check stays silent, so the expected warning never comes
(`dev/cifix-diagsep.R`, `diagsep-ref.log`, `diagsep-ob.log`). At 4000
both stop at code 0, after 1939 and 2090. The test now sets 4000 and
asserts code 0, so a platform that still stops short fails on that
line, by name.

Defect found, not fixed (review, m5): separation is not detected as
such. At the default budget with OpenBLAS, nlminb stops at code 9 after
1997 evaluations with the generic "function evaluation limit reached"
warning, and the word separation never reaches the user, while `glm()`
reports fitted probabilities of 0 or 1 after 25 iterations. Filed in
`dev/test-backlog.md`.

The two escaped warnings. `(1 | Subject/a)` on sleepstudy runs the
Subject:a sd to exp(-29). Its finite-difference Hessian row is
+7.1e-12 with the reference BLAS and -7.1e-12 with OpenBLAS
(`dev/cifix-diagnest.R`, `diagnest-*.log`), noise either way. Positive,
`solve(H)` has a positive diagonal and the check keeps it (tier 1, no
warning, a standard error of about 3.7e5 on the log sd); negative, tier
3 calls the parameter flat and warns. Both tests now allow that warning
(`allow_warnings()`, not required) with the reason in a comment.

Defect found, not fixed: tier 1 of `cov_from_hessian()` (and
`sdr_rescue()`'s first exit) accepts an inverse built on a Hessian row
of pure noise whenever the noise is positive. Tier 3's empty-row rule
(row maximum within its noise, or under `se_empty_row`) would call it
flat on every platform. Applying that rule before tier 1 is a
behavior change with reach across every suite, and it belongs to the
round's owner, not to a CI fix.

## Seen with the emulator only

On rellib-r6 AND on the fix, identically, so not introduced here and
not seen on CI (`fixob.sum`, `t-base-ob-*.log`):

- `test-cumulative-cs.R:132`: `mixture(cumulative, sratio)` stops with
  "NA/NaN gradient evaluation" (error).
- `test-ordinal-mixture.R:751`: "singular convergence (7)" escapes.
- frmtmb.eam `test-sampling.R:29`: two ESS warnings escape.
- frmtmb.sample `test-brms-shapes-draws.R:196`: an R-hat warning
  escapes.

The first two are 0.68.0 tests. They pass on the runner, whose CPU
kernels and LAPACK differ from this build, so they are fragile rather
than failing. Not changed.

Two more warnings escape on the runner and not here, in the
frmtmb.sample job (FAIL 0, WARN 2, SKIP 18, PASS 2274), and neither is
a rounding matter: `test-brms-shapes-draws.R:106` ("Method
'posterior_samples' is deprecated", from the runner's brms) and
`test-compat-preflight.R:95` ("no DISPLAY variable so Tk is not
available", a headless runner). Not changed; listed for the round's
owner.

## Other platform-fragile exactness in the tests 0.68.0 added

`git diff 4b903fa9 d24f7b86` over every `tests/` directory, grepped
for `expect_identical()` and for the SE check's warning; the grep is a
list to read, not evidence. The evidence is a run: all eight suites
ran with the emulator and the fix (table below), and beyond the cases
above the only differences from the reference BLAS are the four files
in the last section.

## Fits that differ between processes

WITHDRAWN: an earlier version of this section said that a fit depends
on what ran before it in the same process. The review measured that
this is wrong (dev/reviews/2026-10-06-cifix.md, item 6, m1).

What I saw: with 0.68.0 and the reference BLAS, the scan's fit at
n = 140 seed 11 gave theta_1 = -11.4241 in one process and -11.3562
in others (`dev/cifix-determ.R`, `det-*.tsv`). I read that as history
dependence, from too few processes.

What the review measured: in fresh processes with identical inputs,
`dev/cifix-determ.R` gave -11.424124791143393 eight times and
-11.356207797215273 four times in 12 runs. Six fits in one process
always agree (48 of 48). Pinned to one logical CPU each, the result
follows the core type: P-cores (Windows CPUs 0-1, 10-13, 22-23 on this
Intel Core Ultra 9 285K) give one value and E-cores the other. So some
compiled code chooses its kernel once per process from a CPUID-based
query, and the AD gradient's rounding follows that choice; the library
was not identified. It is not state in frmtmb, and the CI runners are
not hybrid. The rule it supports stands: on a flat-direction fit,
assert a property that holds on any fit, not a value of one.
`dev/rtmb-pitfalls.md` item 21 records the trap.

## Verification

Installed core and frmtmb.spline from the worktree into wt-cifix-lib;
every other extension loads from rellib-r6 behind it. Every log's
`lib:` line shows wt-cifix-lib's frmtmb. One test file per process
(`dev/cifix-pardrv.sh`, `dev/cifix-run1.R`), NOT_CRAN=true, all eight
suites, ungated:

| Package | files | reference BLAS (pass/fail/error/skip/warn) | OpenBLAS |
|---|---|---|---|
| frmtmb | 214 | 14739/1/0/177/0 | 14734/1/1/177/1 |
| frmtmb.coupling | 11 | 542/0/0/5/0 | 542/0/0/5/0 |
| frmtmb.eam | 29 | 1743/0/0/3/0 | 1743/0/0/3/2 |
| frmtmb.latent | 10 | 360/0/0/2/0 | 360/0/0/2/0 |
| frmtmb.learn | 15 | 430/0/0/13/0 | 430/0/0/13/0 |
| frmtmb.ode | 11 | 549/0/0/1/0 | 549/0/0/1/0 |
| frmtmb.sample | 48 | 2400/0/0/4/0 | 2400/0/0/4/1 |
| frmtmb.spline | 15 | 590/0/0/1/0 | 590/0/0/1/0 |

The one failure in both columns is `test-perf.R`, a wall-clock bound
(1.21 against 1.00) while 22 R processes ran; alone it passes 3/0
(`t-fix-ref-perf.log`). The OpenBLAS extras are the four files above.
After these runs `test-curve.R` gained its test and the spline split
its NaN and NA rows; the spline suite then ran again on the final code,
590/0/0/1/0 with each BLAS (`splref.sum`, `splob.sum`). The only core
change after the runs is a comment.

The files the task named, reference BLAS: test-se-check.R 123/0,
test-gp-by.R 67/0 (2 skips), test-gp-multidim.R 62/0, test-fd-chain.R
11/0, test-predict-re-uncertainty.R 78/0; spline test-difference.R
71/0; the frmtmb.sample suite 2400/0.

Gated tier, reference BLAS, `FRMTMB_BRMS_FIT_TESTS=true`, the 36 files
that skipped above (`gatedref.sum`): frmtmb 26 files 2609/0/0/14/0,
frmtmb.sample 3 files 324/0/0/1/0, frmtmb.learn 2 files 71/0/0/2/0;
the remaining skips are the scale and fuzz tiers and
`test-drmtmb-agreement.R`, which then ran with
`FRMTMB_DRMTMB_FIT_TESTS=true`: 131/0/0/0/0 (`drmref.sum`).

`R CMD check --as-cran`, built and checked in `dev/cifix-check/`, each
Status: 1 NOTE, the V8 math-rendering NOTE. frmtmb: tests FAIL 0,
WARN 0, SKIP 198, PASS 14678. frmtmb.spline, on the final code: tests
FAIL 0, WARN 0, SKIP 3, PASS 582.

What this could not verify: the Ubuntu runner itself. The emulator
matches it on every symptom it showed, but its CPU kernels and LAPACK
differ, so a green run there is the remaining evidence. A push of a
branch would give it.

## Not done, and why

- `vcov()`, `confint()` and `brms-shapes.R` under REML read the joint
  precision through their own `solve()`, and keep doing so; the REML
  warning says it does not repair them, which stays true. Routing them
  through `get_joint_cov()` would change vcov() on REML fits.
- The tier-1 noise-sign defect above, and the review's other
  found-not-fixed items (m2, m5, m6, m9, items 4 and 5): all filed in
  `dev/test-backlog.md` under "Filed at 0.68.1".

## After review: versions and the final pass

The review (`dev/reviews/2026-10-06-cifix.md`) found the patch
MERGEABLE. Applied after it: versions frmtmb 0.68.1 (DESCRIPTION,
codemeta.json) and frmtmb.spline 0.10.1 with `frmtmb (>= 0.68.1)`;
the m3 wording in NEWS, `?frm_joint_cov` and `joint_cov_repair()`;
`lost_pos` documented in `?frm_joint_cov` (m8); the OLRE test (m4);
this record's corrections (m1, m7, m10, m5).

No other extension needs the new floor. None calls `frm_joint_cov()`,
reads `lost_pos` or the kriging internals; frmtmb.sample's
`test-gp-by-draws.R` reads `frm_lp_basis()$extra_var` and `extra_cov`
at two distinct positions, where the dedup changes nothing. All seven
extension suites passed against this core (Verification), and their
floors stay at 0.68.0.

Final pass on 0.68.1 / 0.10.1, installed into wt-cifix-lib, one file
per process (`finalref.sum`, `finalob.sum`, `baseolre*.sum`):

| File | reference BLAS | OpenBLAS |
|---|---|---|
| test-se-check.R | 123/0/0/0/0 | 123/0/0/0/0 |
| test-gp-multidim.R | 62/0/0/0/0 | |
| test-aliased-grouping.R | 33/0/0/0/0 | |
| test-open-issues.R | 43/0/0/0/0 | |
| test-se-lost-re-predict.R | 14/0/0/0/0 | 14/0/0/0/0 |
| test-se-lost-re-predict.R on rellib-r6 | 6/8/0/0/0 | 6/8/0/0/0 |
| frmtmb.spline, 15 files | 590/0/0/1/0 | 590/0/0/1/0 |

`R CMD check --as-cran` frmtmb.spline 0.10.1 against frmtmb 0.68.1:
Status: 1 NOTE (V8); tests FAIL 0, WARN 0, SKIP 3, PASS 582.
