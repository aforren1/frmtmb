# Lane ci701: the macOS failure of 0.70.0

Base: d2ef6312 (frmtmb 0.70.0), built in
`C:/Users/adf44/source/r/rellib-r8`. Worktree `frmtmb-wt-ci701`,
private library `wt-ci701-lib`. Logs are in `dev/ci701-log/`
(gitignored); every number below is copied from them, with the script
that made it. "ref" is R's reference BLAS and LAPACK; "ob26" and "ob32"
are the copies of R that `dev/ciharden-openblas.sh 0.3.26 lapack` and
`0.3.32 lapack` build in this tree (Rblas 74 and 76 symbols to
OpenBLAS, Rlapack 616 and 617; `ob26-build.log`, `ob32-build.log`), run
with `OPENBLAS_NUM_THREADS=4`.

## Fix round 1 (review `dev/reviews/2026-10-13-ci701.md`)

This section supersedes the text below where they disagree.

**R1, D2 on wider blocks.** `se_edge_cor()` now asks only a `us` block
of dimension 2. In a wider block `us_chol_L()` normalizes rows, so
`|t| > 70.7` on one parameter need not put any correlation near
+/-1. `test-v16.R` row 3 ends at t = (1318, -708) with correlations
0.881 and -0.235. Pulling one parameter back also moves the interior
correlations of its row. In the review's case A, seed 3, corr(3, 2)
went from 0.19674 to 0.20613 (`dev/ci701-rev-d3c.R`). Every
measurement behind D2 is of dimension 2, so the restriction loses
nothing measured. A wider block keeps 0.70.0's handling.

New test, `test-review-v29.R` "a block of dimension 3 keeps the 0.70.0
check": case A, seed 3, started at t21 = 2000. It asserts that a
correlation parameter is past the clamp, that `se_edge_cor()` is all
FALSE, and that `cov_from_hessian()` gives an identical lost set and
covariance with `se_edge_cor()` mocked to FALSE. On the lane build
before this round it fails 3 assertions: `se_edge_cor()` fires on
theta_4, and the lost set and covariance differ
(`t3-prefix-ref-v29.txt`). On rellib-r8 it errors, because
`se_cor_clamp_theta()` does not exist there.

New case in "a boundary correlation reports NA bounds, not a clamp":
the review's natural code-0 stop (data moved with seed 8004, start
t = 2000). On rellib-r8 its sd rows are NA and the forced analysis
loses theta_1 and theta_2 as flat (`t3-base-ref-v29.txt`, lines 434
to 436 and 442, 445). On the lane build the sd intervals are finite
and the SEs are within 1 percent of the rank-one fit.

**R2, accuracy.** The "19 of 20 flat" figure came from forcing the
check on fits where it does not run. The natural case is the code-0
stop above. The 0.55 percent SE move is a correction toward an
independent reference (see "What the pull-back changes elsewhere").
NEWS says both. Its long sentences are split.

**The code-1 path.** D2 also runs when a fit with code 1 reaches the
analysis. That happens through `se_boundary_stop()` in
`check_convergence()` (gradient at most `grad_tol`) and through
`nl_flat_message()`. In the full scan the pull-back ran 66 times
against 10 code-0 fits past the clamp. So it runs on code-1 fits, and
more than once per fit: `cov_from_hessian()` takes a second full
Hessian (2 np gradients) on every call where `se_edge_cor()` holds.
That is 90 gradients on the review's RB4 fit. A test covers the
`se_boundary_stop()` path now. In "a lost correlation leaves the sd
intervals of its block", seeds 7 and 8 end with code 1 at |t| of 200
and 169. The pull-back runs there, and the test asserts that
`fit$cache$se_boundary_stop` is TRUE. In the review's runs no verdict
on a code-1 fit changed (R2, R6). The `nl_flat_message()` path has no
test.

**Rerun after the fix** (core reinstalled into `wt-ci701-lib`; all 96
files load it). The job list `dev/ci701-log/fix2.jobs` holds 32 files:
test-review-v29, test-se-check, test-se-tier, test-se-flags, test-v16,
and every test file whose formulas grep as a `us` block of 3 or more
terms or a shared `|p|`, `|q|` or `|ID|` block (27 core files, two of
frmtmb.learn, one of frmtmb.sample). `dev/ci701-scan.sh`, every gate
on, on ref, ob26 and ob32 (`fix2ref.sum`, `fix2ob26.sum`,
`fix2ob32.sum`): "ran 32 of 32" each, 32 of 32 files with fail, error,
skip and warn all 0. test-review-v29.R 186/0/0/0/0, test-se-check.R
127/0, test-se-tier.R 91/0, test-se-flags.R 68/0 and test-v16.R 24/0
on every arm. `dev/ci701-scan-diff.R` finds no difference between the
three arms in any RESULT, DETAIL or TEST line.

**Filed as pre-existing.** Under `check_se = "ignore"`, `vcov()` and
`summary()` on rs20 seeds 7 and 12 (`y ~ x + (1 + x | g)` with no
slope variance) let two "NaNs produced" warnings escape. This happens
on 0.70.0 and on the lane build (review R6). It is in "Found, not
fixed" below.

## What failed on CI

Run 38457152625, job 115502711709, macos-latest release (arm64), R CMD
check `tests`: FAIL 4, WARN 1, SKIP 200, PASS 16377
(`dev/ci701-macos-r-cmd-check.log`). Every other job passed.

1. `test-review-v29.R:376`: `expect_warning(confint_varcorr(fit),
   "boundary")` saw no warning, and `:387` to `:389`: both sd rows of
   `confint_varcorr()` were NA, on the `(1 | q | id)` fit whose
   correlation runs to +1 (`set.seed(61)`).
2. `test-ordinal-mixture.R:359`: "Optimizer did not report
   convergence: singular convergence (7)" escaped `allow_warnings()`.

## Item 1: the boundary correlation

### What the fit is

The correlation of the merged block has its maximum at +1, so the
optimizer runs its parameter `t` (`rho = t / sqrt(1 + t^2)`) out until
it stops. The model at the boundary is rank one,
`y_r = a_r + b_r x + s_r u_id + e_r`, `u_id ~ N(0, 1)`, which RTMB
fits with no near-singular covariance (`dev/ci701-bcor-rank1.R`,
`bcor-rank1-ref.log`):

- rank-one maximum: logLik -67.7319914882, code 0; frmtmb's default fit
  -67.7319918044 (3.2e-7 below), at the same log sds to 1e-5
  (-0.300257 against -0.300250, -0.797233 against -0.797243);
- rank-one standard errors of the log sds: **0.217708 and 0.222476**.
  So the two sds are identified at the boundary, and withholding their
  intervals is wrong.

frmtmb's objective with `t` held and the rest re-maximized, the
finite-difference Hessian of the two log sds there, and the largest
gradient of the two log sds over 6 ulp-moved points (same log):

| t | logLik | H11 | H22 | se(log sd1), se(log sd2) | max abs gradient |
|---|---|---|---|---|---|
| 10 | -67.8299420148 | 77.58 | 74.19 | 0.21455, 0.21949 | 0.000167 |
| 30 | -67.7425344021 | 80.8 | 77.39 | 0.21738, 0.22218 | 0.000239 |
| 70.7 | -67.7338823492 | 81.16 | 77.74 | 0.21765, 0.22243 | 0.000212 |
| 300 | -67.7320964166 | 81.23 | 77.82 | 0.21770, 0.22248 | 8.35e-05 |
| 1000 | -67.7320009299 | 80.66 | 77.97 | 0.21816, 0.22193 | 0.00124 |
| 3000 | -67.7319925379 | 139.9 | 69.98 | 0.27484, 0.39875 | 0.0768 |
| 10000 | -67.7319915182 | 580.1 | 3045 | 0.20183, 0.20221 | 7.9 |

Below about t = 1000 the objective has its digits; above it the
gradient is noise that grows fast with t. The cause is the dense
`RTMB::dmvnorm()` of the `us` block: it factorizes a covariance that is
singular to within `1 / t^2`, so the Cholesky cancels away the digits
the outer gradient is made of. (At t = 70.7, the clamp of
`confint_varcorr()`, a correlation of 0.9999, the standard errors are
within 1.3e-4 of the rank-one ones.)

### Where the fit stops, and what each side reports

`dev/ci701-bcor-probe.R`: the test's data, then every `x`, `y1`, `y2`
moved by 0 to 4 ulps (seed 1000 + k), 40 fits per arm on rellib-r8:

| arm | code 0 | code 1, false convergence (8) | code 1, singular convergence (7) |
|---|---|---|---|
| ref (`bcor-base-ref-*.log`) | 2 | 35 | 3 |
| ob26 (`bcor-base-ob0.3.26-*.log`) | 2 | 35 | 3 |
| ob32 (`bcor-base-ob0.3.32-*.log`) | 2 | 35 | 3 |

(ob26 and ob32 gave the same fits; ref differs from them.) All 120 fits
ended with a convergence or a gradient warning, so the standard-error
check did not run (`se_explained` "convergence"), and all 120 had
finite sd intervals: that is the side Windows and Ubuntu took. `t`
ended between 1665.9 and 13370.6, and the largest gradient between
0.00043 and 8.3 (`dev/ci701-bcor-alt.R`, `balt-*.log`).

The other side, where nlminb reports convergence and no gradient
warning fires, is the macOS runner's. It is reached here from a moved
start (`dev/ci701-bcor-start.R`, `bstart-base.log`): from
`theta = (0, 0, 500)` and `(0, 0, 1000)` the fit stays at that `t` with
code 0 and no warning, gives the boundary message, and
`confint_varcorr()` gives NA for both sd rows and no warning, which is
the CI failure exactly. From `t0` = 100 and 200 it is code 0 with
theta_3 "flat" and one warning, sd rows NA too; every other `t0` from
0.5 to 3000 ends on the warning side.

### Two defects on the checked side

D1, `confint_varcorr()`. `varcorr_trans_rows()` read `sdr$cov.fixed`,
whose lost rows are NaN, and the delta-method product spread the
correlation's NaN to every row of the block. `VarCorr()` had this fixed
at lane nanse (`hyp_par_cov()`, `hyp_prop_var()`): it reads
`cov_fixed_prop` and gives NaN only to rows that move along a lost
direction (`jc_nonest()`). This does not need the boundary correlation
fit: on `y ~ x + (1 + x | g)` with no slope variance, seeds 1 to 20
(`dev/ci701-vc-spread.R`, `vcspread-base-ref.log`), the check removed
only the correlation (`theta_3` "boundary") in 8 fits (seeds 2, 7, 8,
11, 12, 15, 17, 19), and in all 8 the intercept sd had no interval,
while `VarCorr()` gave it a standard error. That outcome is stable:
10 of 10 ulp-moved fits per seed for seeds 2, 7 and 8
(`dev/ci701-vc-check.R`, `vccheck-fix-ref.log`).

D2, the standard-error check. At the fit's own point the check reads
the Hessian at t in the thousands, which is noise (table above).

The natural case (found by the review, `dev/ci701-rev-tstart.R`,
`dev/ci701-rev-log/tstart-*.log`): the review-v29 data moved by a few
ulps (seed 8004) and started at t = 2000. The fit ends with code 0 at
t = 2000, gradient 0.0069, and the check runs. rellib-r8 loses
`theta_1:flat, theta_2:flat, theta_3:boundary`. It gives NA sd rows,
the boundary message and the flat-SE warning. That is the macOS
signature. The lane build loses `theta_3` alone, with finite sd
intervals and se 0.21768 and 0.22245.

How often the check would say this. The 20 ref fits on the lane's
data (seed 1000 + k) all end with a convergence warning, so the check
does not run on them. Forced to run (`se_fit_analysis(fit, force =
TRUE)`, `dev/ci701-bcor-force.R`, `bforce-base-ref-*.log`), it called
both sds flat in 19 of 20. That figure is a measure of the Hessian at
those points, not of what a user sees. There the unit-diagonal matrix
had a largest eigenvalue of 243 to 1.3e6 and a negative one of the
same size. theta_3's diagonal is about 1e-9, so the scaling blows its
noisy off-diagonal entries up. Its most negative direction loaded 0.26
to 0.66 on each log sd.

So on the macOS side the sd rows are NA through D2 when the check runs
at a t in the thousands. They are NA through D1 even where the check
does not call them flat (the `t0 = 500` fit).

Which of the two the runner hit cannot be read from its log; both are
fixed.

### The fixes

`R/confint.R`, `varcorr_trans_rows()`: the generic (jacobian) branch
reads `cov_fixed_prop` and gives NA to a row `jc_nonest()` says moves
along a lost direction, as `VarCorr()` does. The single-entry branches
(smooth, gp, car, spde) keep reading `cov.fixed`; their NaN is their
own. A block with nothing lost reads the same numbers as before.

`R/se-check.R`:

- `se_edge_cor()`: the correlation parameter of a `us` block of
  dimension 2 with `|t|` above the clamp's value
  (`se_cor_clamp_theta()`, 70.7), where `se_at_edge()` says the
  likelihood does not change further out. Not asked under quadrature,
  as `se_edge_sd()`. Dimension 2 only since fix round 1: see that
  section.
- `cov_from_hessian()`: when there is such a parameter, the Hessian is
  taken again with it pulled back to the clamp (`fit_outer_hessian(fit,
  pull =)`, `se_pull_back()`), and the parameter joins `edge`, so tier 3
  removes it before the decomposition; `se_boundary_names()` then calls
  it "limit", a boundary, as before.

Measured on the same fits (`balt-*.log`, candidate b, before the code
existed; then `bforce-fix-*.log` on the lane build):

- Hessian with `t` pulled to the clamp: se(log sd) 0.21767 or 0.21768
  and 0.22244 or 0.22245 in **120 of 120** fits (ref, ob26, ob32), that
  is within 1.3e-4 of the rank-one 0.217708 and 0.222476.
- For contrast, the fit-point Hessian with theta_3 only dropped
  (candidate a) gave NaN in 21 of 120 fits and values from 0.0579 to
  0.699 otherwise. Dropping the row is not enough.
- The full inverse at the fit point (what sdreport() shows on the
  warning side) gives 0.21203 and 0.21203 in 35 of 40 ref fits and
  0.21180 to 0.21206 in the rest (`balt2-ref-*.log`): finite, 2.6 and
  4.7 percent low, because the noisy row constrains the sds along one
  direction. That side is left alone: the fit carries a convergence
  warning there, and `sdr_rescue()` keeps sdreport()'s answer on such a
  fit by design.
- Lane build, forced analysis: `lost` is `theta_3 = "boundary"` alone
  in 120 of 120 fits, with shown se 0.2177 and 0.2224 or 0.2225
  (`bforce-fix-ref-*.log`, `bforce-fix-ob0.3.2*-*.log`).
- Lane build, moved starts (`bstart-fix.log`): every start gives finite
  sd intervals; on the checked side [0.4834, 1.1347] and [0.2913,
  0.6968], on the warning side [0.4888, 1.1222] and [0.2974, 0.6827].
- D1 (`vcspread-fix-ref.log`): the 8 seeds have finite intercept-sd
  intervals, for example seed 2 [0.6046, 1.2154], equal to `VarCorr()`'s
  delta method (Est.Error 0.1527 on 0.8572).

What the pull-back changes elsewhere. `dev/ci701-rs20-se.R` on the 20
`(1 + x | g)` fits, base against lane (`rs20se-base.log`,
`rs20se-fix.log`): the lost set is identical in 20 of 20. The 6 fits
whose correlation is past the clamp and at its edge (seeds 7, 8, 11,
12, 15, 19, |t| 169 to 354) see their standard errors move. The
log-sd of the slope moves by up to 0.55 percent (seed 12: 1.21476 to
1.22143), the others by 1e-5 relative or less. Seeds 2 and 17 (|t| 68
and 58) are unchanged bit for bit.

The move is a correction. The review fitted the rank-one model in
closed form, with no RTMB (`dev/ci701-rev-se.R`, `se-base.log`,
`se-lane.log`). Against it the largest relative SE error on these 6
fits is 7.6e-4 to 5.3e-3 on 0.70.0 and 1.2e-4 to 2.8e-4 on the lane.
Seed 12 goes from 5.3e-3 to 1.3e-4.

### The test

`test-review-v29.R`, "a boundary correlation reports NA bounds, not a
clamp", now asks what holds on both sides, for the default start and
for `theta = (0, 0, 500)` (the checked side on every arm here):

- the boundary is said exactly once, by the fit's boundary message or
  by `confint_varcorr()`'s warning, and nothing else warns there;
- the correlation row is NA, the sd rows finite and around the
  estimate (as before);
- `se_fit_analysis(fit, force = TRUE)` loses `theta_3 = "boundary"`
  alone, and its sd standard errors are within 1 percent of the
  rank-one model's, fitted in the test with RTMB (measured gap 1.3e-4).

New test, "a lost correlation leaves the sd intervals of its block":
seeds 2, 7, 8 of the `(1 + x | g)` shape; sd intervals finite and equal
to `VarCorr()`'s delta method (relative 1e-6); where the check removed
something, it is `theta_3 = "boundary"` alone.

Seen to fail on rellib-r8, all three arms, pass=155 fail=11
(`t2-base-{ref,ob0.3.26,ob0.3.32}-v29.txt`): the `t0 = 500` fit's sd
rows NA (`:426` to `:428`, the CI failure itself), the default fit's
forced analysis losing the sds and fixed effects (`:435`, `:438`), and
the three seeds of the new test (`:459`, `:463`, 3 each). Lane build:
pass=166 fail=0 in all three arms (`t2-fix-*-v29.txt`).

## Item 2: "singular convergence (7)" on the degenerate ordinal mixture

`dev/ci701-omx-probe.R`: the test's fit, `x` moved by 0 to 4 ulps
(seed 2000 + k), on rellib-r8:

| arm | fits | relative convergence (4) | singular convergence (7) |
|---|---|---|---|
| ref (`omx-ref-*.log`) | 20 | 13 | 7 |
| ob26 (`omx-ob0.3.26-0.log`) | 10 | 9 | 1 |
| ob32 (`omx-ob0.3.32-0.log`) | 10 | 9 | 1 |

Every fit warns that component 1 is a step function at a degenerate
boundary. The code-7 fits are not worse: logLik -380.18978 to
-380.18973 against -380.18996 to -380.18973 for code 0, largest
gradient 7.96e-05 to 0.000707 against 5.27e-05 to 0.000706. PORT calls
the Hessian singular because the degenerate component's directions are
flat, which the required warning names. Not a defect of the fit; the
same decision as lane ci691's for `test-nonfinite-gradient.R`.

Fix, test only: "did not report convergence" is allowed, not
required. `dev/ci701-omx-body.R` runs the revised test body on 20 ref,
10 ob26 and 10 ob32 ulp-moved fits (`omxbody-fix-*.log`): 31 code 0, 9
code 1, the required warning in 40 of 40, the frequency assertion in
40 of 40, nothing escaped.

## Siblings

The scans run every test file of all eight suites, one process per
file, every gate on (`dev/ci701-scan.sh`, a copy of
`dev/ciharden-scan.sh` that runs `dev/ci701-run1.R`, which is
`dev/ciharden-run1.R` with two traces that change no number):
`CI701COR` prints each fit whose `us` correlation parameter ended past
the clamp, with its optimizer code; `CI701EDGE` prints each time
`se_edge_cor()` removed one.

On the reference arm of the lane build (`fixref-suite/`), 73 fits in
18 files ended past the clamp, 10 with code 0 and 63 with code 1, and
the pull-back ran 66 times in 15 files. Only 5 of the 73 ended above
t = 1000, where the gradient is noise: the review-v29 fit itself, two
fits of the gated fuzz tier (2346.6, 2974.9; CI does not run it),
`test-v16.R`'s `us(spp + 0 | site)` fit (1318.4, under
`suppressWarnings()`, asserting a log-likelihood agreement that does
not depend on the code), and frmtmb.learn `test-stan-identity.R`'s
4 x 4 rlddm block (2413.0, gated).

The last allows only "false convergence (8)" and the gradient warning,
so a "singular convergence (7)" there would escape as it did on the
ordinal mixture. `dev/ci701-sib-rlddm.R` (rt moved by 0 to 4 ulps,
seed 5000 + k): 10 of 10 fits false convergence (8), t 1324.6 to
3099.8, logLik -731.359032 to -731.358989, nothing escaped
(`sib-rlddm-fix-ref.log`); ob26 10 of 10 false convergence (8), t 1599.9
to 2890.0 (`sib-rlddm-fix-ob26.log`); ob32 10 of 10, t 1599.9 to 2890.0
(`sib-rlddm-fix-ob32.log`). 30 of 30 on the allowed side, so the test
is left as it is; it is the one to widen if a runner ever reports
"singular convergence (7)" there.

The other ordinal-mixture fits whose `allow_warnings()` list lacks
"did not report convergence" (`dev/ci701-sib-probe.R`, x moved by a
few ulps, seed 4000 + k): `test-ordinal-mixture.R:401` (the `cs()`
mixture) 20 of 20 relative convergence (4), and frmtmb.sample
`test-ordinal-mixture-draws.R`'s "none" fit 20 of 20 relative
convergence (4) (`sib-cs-ref.log`, `sib-none-ref.log`). Neither is
near the branch.

`expect_warning(confint_varcorr(...))` appears nowhere else in the
suites.

The ob26 and ob32 arms give the same trace counts as the reference arm
in every package but core, where the pull-back ran 63 times against
62 (see Verification). Across the three arms one file differs, and it
is the coverage-only row lane ciharden recorded: `test-v11.R` "toep
matches a hand-rolled reference (glmmTMB when it converges)" passes 1
assertion on ref and 2 on ob26 and ob32 (glmmTMB converges or not).

## Verification

Core installed from the worktree into `wt-ci701-lib` (roxygenised; no
generated file changed). The two changed test files, one per process
(`t2-*-v29.txt`, `t-*-ref-test-ordinal-mixture.R.txt`;
pass/fail/error/skip/warn):

| File | base ref | base ob26 | base ob32 | lane ref | lane ob26 | lane ob32 |
|---|---|---|---|---|---|---|
| test-review-v29.R | 155/11/0/0/0 | 155/11/0/0/0 | 155/11/0/0/0 | 166/0/0/0/0 | 166/0/0/0/0 | 166/0/0/0/0 |
| test-ordinal-mixture.R | 99/0/0/0/0 | | | 99/0/0/0/0 | | |

(The base columns run the revised tests against 0.70.0. On this
machine the original ordinal-mixture test passes on both builds: the
default data end at code 0 here; the escape is shown on the moved data
above.)

All eight suites, every file in its own process, every gate on
(`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS`, `FRMTMB_FUZZ`,
`FRMTMB_SCALE_TESTS`), lane build, by `dev/ci701-scan.sh`; totals by
`dev/ci701-scan-table.R` (`fixref.sum`, `fixob26.sum`, `fixob32.sum`:
"ran 371 of 371" each; the core library on all 371 files of each arm
is `0.70.0 C:/Users/adf44/source/r/wt-ci701-lib`; the extensions load
from rellib-r8):

| Package | files | ref pass | ob26 pass | ob32 pass | fail, error, skip, warn (all arms) |
|---|---|---|---|---|---|
| frmtmb | 229 | 18974 | 18975 | 18975 | 0, 0, 0, 0 |
| frmtmb.coupling | 11 | 548 | 548 | 548 | 0, 0, 0, 0 |
| frmtmb.eam | 29 | 1765 | 1765 | 1765 | 0, 0, 0, 0 |
| frmtmb.latent | 10 | 366 | 366 | 366 | 0, 0, 0, 0 |
| frmtmb.learn | 15 | 515 | 515 | 515 | 0, 0, 0, 0 |
| frmtmb.ode | 11 | 551 | 551 | 551 | 0, 0, 0, 0 |
| frmtmb.sample | 51 | 2806 | 2806 | 2806 | 0, 0, 0, 0 |
| frmtmb.spline | 15 | 592 | 592 | 592 | 0, 0, 0, 0 |

The gated ported brms-suite files pass in all three arms, so the change
flips no row of `dev/brmsport-verdicts.tsv`. `dev/stan-cache` here was
seeded from the main checkout (121 programs; the seeded union was not
in this worktree), and the gated tiers compiled 70 programs fresh (191
in the cache after the scans);
`dev/brms-suite/` is brms 2.23.0's CRAN tarball, sha256
b5f5bb5604ec3f87b3ad99f0da4e6a37d8c5488221e8abe459171340992c37a5, the
one `dev/brmsport-blocks.R` asserts.

`R CMD check --as-cran` of frmtmb on the final tree, built and checked
in `dev/ci701-check/frmtmb/` (`dev/ci701-check.sh`; the Rcheck tree is
deleted, `check.log`, `keep/00check.log` and `keep/testthat.Rout` are
kept): Status: 1 NOTE, the V8 math-rendering NOTE on the HTML manual.
Tests: FAIL 0, WARN 0, SKIP 200, PASS 16409 (the macOS runner on
0.70.0: FAIL 4, WARN 1, SKIP 200, PASS 16377).

## Version and files

A patch: frmtmb 0.70.1. No extension changes and none needs a new
floor. Files touched: `R/confint.R` (`varcorr_trans_rows()`),
`R/se-check.R` (`fit_outer_hessian()` gains `pull`,
`cov_from_hessian()`, new `se_edge_cor()`, `se_cor_clamp_theta()`,
`se_pull_back()`), `NEWS.md`, `tests/testthat/test-review-v29.R`,
`tests/testthat/test-ordinal-mixture.R`, `dev/test-backlog.md` (fix round 1), and
the lane's `dev/ci701-*`
scripts and this file. No Rd file changed.

## What was not done, and why

- The warning side of the boundary-correlation fit keeps sdreport()'s
  standard errors from the noisy Hessian (2.6 and 4.7 percent low) and
  its convergence and gradient warnings, which are false alarms at a
  point 3.2e-7 below the maximum. Repairing either changes the
  convergence verdict, which is outside a patch.
- The root cause, the dense `dmvnorm()` of a `us` block losing digits
  as the correlation runs to +/-1, is not touched: evaluating the block
  from its Cholesky factor would change the rounding of every `us`
  fit. It is for the next minor; filed in `dev/test-backlog.md`
  ("Filed by lane ci701").
- `se_edge_cor()` covers `us` blocks of dimension 2 only, and the
  trace looked at `us` blocks only. A correlation of another structure
  that runs to +/-1 (`ar1`, `hetar1`, `homcs`) keeps 0.70.0's handling;
  whether its density loses digits the same way was not measured.
- `confint_varcorr()`'s warning says "on the boundary of the parameter
  space" for every NA row, including a row the check called flat
  rather than at a boundary. The wording predates this lane.
- Two "NaNs produced" warnings escape `vcov()` and `summary()` under
  `check_se = "ignore"` on rs20 seeds 7 and 12, on 0.70.0 and on this
  build (review R6). Filed in `dev/test-backlog.md`.

## Fix round 1: R CMD check and cleanup

`R CMD check --as-cran` of the fix-round tree (`dev/ci701-check.sh
frmtmb-r1`; logs kept in `dev/ci701-check/frmtmb-r1/keep/`): Status:
1 NOTE, the V8 math-rendering NOTE. Tests: FAIL 0, WARN 0, SKIP 200,
PASS 16429. `dev/ciharden-out` (the OpenBLAS copies of R, 925 MB) is
deleted from this worktree. To rerun an OpenBLAS arm, rebuild it with
`dev/ciharden-openblas.sh`.
