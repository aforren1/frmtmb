# Lane ci082: the red ubuntu CI of frmtmb 0.82.0

Round of 2026-11-13, base ec8178d3 (frmtmb 0.82.0), branch
`wt-r83-ci082`. Brief: `round-briefs/r083-ci082.md`. This file is the
record; logs are local (`dev/ci082-out/`, gitignored).

Status: FIX ROUND 2 DONE (short check `dev/reviews/2026-11-13-ci082-short.md`; review `dev/reviews/2026-11-13-ci082.md`,
verdict FIX; see "Fix round" at the end, which supersedes item 1's rule
and the TMB-7 entry below). Everything is uncommitted in the worktree;
the lane library holds the change (core only; no extension changed).

## Summary

| CI failure | reproduced on | cause | fix |
|---|---|---|---|
| covstruct-seam:126-128 | k = 0: every OpenBLAS setting tried (0.3.26 and 0.3.32, core type unset, Haswell, Zen, Sandybridge; `dev/ci082-out/base-xcs-k0-*.txt`), not the reference BLAS; 13 of 20 ulp perturbations under 0.3.26, 11 of 20 on the reference BLAS | zero-se rule read a rounding zero on a fixed, not estimated, correlation | code: rows the structure fixes (same value at its start and at theta +/- 1, to 1e-10 relative) get se 0 and a zero-width interval (fix round) |
| mixture-starts:333 | reference BLAS, x scaled by 1 + 48e-15; at k = 0 on no setting probed (reference; 0.3.26 unset, Haswell, Zen, Sandybridge; 0.3.32 unset) | the test pinned one of two valid endings | test: the signal matches its own run's record |
| se-limit:98 | reference BLAS k = 3, 4; OpenBLAS 0.3.26 k = 8, 18; 0.3.26 Sandybridge k = 6, 18; at k = 0 on no setting probed (reference; 0.3.26 unset, Haswell, Zen, Sandybridge) | nlminb code 8 at a supremum; "ignore" keeps the convergence warning by design | test: accepts the code-8 ending |
| mixture-starts:182 (devel) | not run (no R-devel here); R-devel source read | `RNGkind()` has a fourth kind | code and test: every kind, by position |

Found and filed, not fixed: a false gradient warning at the poisson
limit from TMB's `dnbinom_robust()` derivative (3b, MEDIUM, RTMB-4);
"ignore" and a code-8 stop (LOW, a user decision).

## The CI failures (`dev/ci082-ci-failed.log`, run 43407296535)

- ubuntu-24.04 release (R 4.6.2): FAIL 3, WARN 1.
  `test-covstruct-seam.R:126` (warning) and `:128`,
  `test-mixture-starts.R:333`, `test-se-limit.R:98`.
- ubuntu-24.04 devel: FAIL 4, WARN 1. `test-covstruct-seam.R:126/128`
  and `test-mixture-starts.R:182` three times (one per RNG case).
  `:333` and `test-se-limit.R:98` passed on devel.

## Files touched

- `R/confint.R`: `varcorr_trans_rows()` (three small hunks: the call
  to `varcorr_fixed_rows()` beside the jacobian, `se_g[fixed] <- 0`
  before the lost-direction rules, the final NA line), the new
  `varcorr_fixed_rows()` after `varcorr_us_row_edge()`, and the
  `confint_varcorr()` roxygen paragraph; `man/confint_varcorr.Rd`
  regenerated. Lane selimit2 owns this file's estimability; these hunks
  are in the covariance-rows code only.
- `R/mixture-starts.R`: `mixture_starts_draws()` restore and its doc.
- `tests/testthat/test-covstruct-seam.R`: one new test.
- `tests/testthat/test-mixture-starts.R`: the RNG test (lines 169-205),
  the `on.exit()` of the next test, the line-333 assertion. Lane
  priors2 edits the `oms_fit()` helper and the `step =` line of
  `wfit()`; this lane left both as they were.
- `tests/testthat/test-se-limit.R`: the `check_se = "ignore"` assertion.
- `NEWS.md`: a new `# frmtmb (development version)` section.
- `dev/test-backlog.md`: a new "Filed in the round of 2026-11-13"
  section, "From lane ci082", two items.
- `dev/upstream-bugs.md`: the existing RTMB-4 entry, with this lane's
  numbers and consequence added (the TMB-7 entry of the first round is
  gone; it duplicated RTMB-4).
- New: `dev/ci082-*.{R,sh,jobs,md}`.

Bump: patch (0.82.1). Floor: none new; R-devel's fourth RNG kind is
handled on any R version.

## Instruments

- `dev/ci082-probe.R <lib|base> <case> [k] [n]`: one case per process.
  Cases `xcs`, `xcsk` (k = 0..19), `nb`, `nbk` (k = 0..19), `mix`,
  `mixk` (k from `k`, `n` of them). k scales a column by
  `1 + 2e-15 * k`.
- `dev/ci082-coretypes.sh <arm> <case> <k> <ver...>`: the probe over
  `OPENBLAS_CORETYPE` default, Haswell, SkylakeX, Zen, Sandybridge.
- OpenBLAS builds: `bash dev/ciharden-openblas.sh 0.3.26 lapack` and
  `0.3.32 lapack` (logs `dev/ci082-ob326.log`, `dev/ci082-ob332.log`).
  `OPENBLAS_NUM_THREADS=4`, as on the 4-core runners.
- SkylakeX exits 132 (illegal instruction) on every run: this CPU
  (Core Ultra 9 285K) has no AVX-512, so that kernel cannot be tried
  here. A runner on an AVX-512 Intel host would use it; that arm is
  untested.
- `dev/ciharden-out/` (925 MB) and `dev/ci082-check/` were deleted at
  the end.

## 1. `test-covstruct-seam.R:126-128`: the rule decided on rounding

Reproduced on rellib-r21 with OpenBLAS 0.3.26 + LAPACK, core type
unset, k = 0 (`dev/ci082-out/ob326-xcs.txt`): the three cor rows have
`se_t` NA, NA bounds and the "No interval for 3 components" warning, as
on the runner. On the reference BLAS at k = 0 the same rows have
`se_t` 5.8392760179632705e-13 and a width of 1.7167933741291108e-12
(`dev/ci082-out/ref-xcs.txt`), which passes.

Cause. `xcs()` fixes its correlation at 0.5; it is not estimated. The
jacobian row of atanh(cor) is a central difference of two evaluations
of `cov2cor(exp(2 theta) C)`, which are equal or a few ulps apart by
how the products round at the two points. So the delta-method se is
either exactly 0 or about 1e-12, and 0.82.0's rule "a zero se is NA"
(lane boundfam, meant for an ESTIMATED component whose se rounds to 0)
read the exact 0 as a lost estimate.

Over 20 ulp perturbations of y (`dev/ci082-probe.R base xcsk`):

| build | se exactly 0, NA bounds, warning | finite tiny se |
|---|---|---|
| reference BLAS, rellib-r21 | 11 | 9 |
| OpenBLAS 0.3.26, rellib-r21 | 13 | 7 |
| OpenBLAS 0.3.26, lane | 0 (20 of 20 se set to 0, width 0, no warning) | 0 |

(`dev/ci082-out/{ref-base,ob326-base,ob326-lib}-xcsk.txt`.)

Decision: the test's expectation is right and the rule was wrong. A
component no parameter moves is known exactly, so its interval is its
value with no width; a zero se on a component that IS moved stays NA
(boundfam's mock test, `test-se-limit.R` "a zero standard error gives
no interval", is unchanged and passes). The rule now decides "not
moved" against rounding, not on `se > 0`: a row is fixed when, for every
theta, the difference `gfun(t + h) - gfun(t - h)` is within
`16 * eps * max(|g0|, 1)`, the tolerance the lost-direction test beside
it already uses. Such a row gets se 0 and keeps it; the boundary rules
(`bd`, a non-finite estimate) still take precedence. For the fixed cor,
g0 = atanh(0.5) = 0.549 and the two evaluations differ by a few ulps
(1.1e-16 each), against a tolerance of 3.6e-15; an estimated row moves
by about 1e-5 over the same step.

New test `test-covstruct-seam.R` "a fixed correlation has no width on
every rounding": 20 perturbations, exact zero width
(`identical(lwr, estimate)`), and the guard-absent case, an estimated
`(0 + t | g)` correlation that keeps a positive width. On rellib-r21 it
fails 40 expectations on the reference BLAS (11 warnings) and 40 under
0.3.26 plus line 128 (14 warnings) (`dev/ci082-base/covstruct-*.txt`).

## 2. `test-mixture-starts.R:333`: a test that pinned one ending

Reproduced on the reference BLAS, rellib-r21, x scaled by
`1 + 2e-15 * 24` (`dev/ci082-out/ref-mixk2.txt`):

    RESULT mixk k 24 scaled FALSE ... reached 9 gain 0 above 5 tier 0 step TRUE
    RESULT mixk k 24 scaled TRUE  ... reached 7 gain 0 above 0 tier 0 step FALSE

The step-start warning says that some of the 10 jittered starts reached
a higher point where a component is a step function. How many do is
each run's optimizer path, and the plain and autoscaled fits run on
different scales. Counts of `degenerate_above` over 240 fits (k = 0 to
39, plain and scaled; reference BLAS, OpenBLAS 0.3.26 default and
Haswell, 0.3.32 default and Zen; `dev/ci082-out/*mixk*.txt`):

    above:  0   1   2   3   4   5   6
    fits:   1  23  66  68  38  34  10

Every one of the 240 fits ends at the same maximum (logLik
-388.095489, tier 0, gain 0), and in each the warning is present
exactly when its own record has `degenerate_above > 0`. So the two
endings are both valid: the fit is the same, and the signal reports
what that run's starts found. The test now asserts that the signal
matches the record of the one run the fit carries (with `picks == 1`
already asserting one run). That still fails the defect the test was
written for: a pre-fit that muffled the warning leaves a record with
`degenerate_above > 0` and no warning. No code change.

## 3. `test-se-limit.R:98`: nlminb's code at a supremum

The warning is "Optimizer did not report convergence: false
convergence (8). A Newton step from the stop gains ... within
frmtmb_control(grad_tol = 0.001), so the estimates are at the maximum
...", class `frmtmb_conv_code`.

The limit rule's verdict does not flip. Over 20 ulp perturbations of x
(`dev/ci082-probe.R base nbk`), check_se = "warning" gave 0 warnings
and the boundary message on 40 of 40 fits (reference BLAS and OpenBLAS
0.3.26). What flips is nlminb's code at the supremum, where the log
shape runs to 15 to 21: "false convergence (8)" on k = 3, 4 (reference)
and k = 8, 18 (0.3.26), "relative convergence (4)" otherwise
(`dev/ci082-out/{ref,ob326}-nbk.txt`). Under check_se = "warning" the
SE check reads a code-8 stop whose lost SEs are all limits as a
boundary stop (`se_boundary_stop()`); under "ignore" that check is
off, so the convergence warning stays. That is the documented contract
(`se_boundary_stop()`'s doc) and `test-se-conv.R` pins it ("and so it
does where the standard errors are not asked about"). The test pinned
the code-4 ending. The same file's first fit (line 86) passed on the
runner, which matches: code 8 under "warning" gives 0 warnings and 1
message.

Fix: the "ignore" assertion accepts the two valid endings: code 0 and
silent, or a nonzero code with exactly one warning, the convergence
one; never a message (no boundary report). No code change. The old
`expect_silent()` fails on the code-8 ending: k = 3 on the reference
BLAS gives 1 warning (`ref-nbk.txt`); the new branch's expectations are
what that probe prints there (0 messages, 1 warning starting
"Optimizer did not report convergence: "). Whether "ignore" should
still ask the boundary question for a stop with a nonzero code is a
user decision (backlog, LOW, filed).

I first ran the "ignore" fit over 20 ulp perturbations in the test.
Under `OPENBLAS_CORETYPE=Sandybridge` (0.3.26 and 0.3.32) that failed
once, at k = 15, with a THIRD ending that is a defect, so the test
keeps the k = 0 data only, as before, and the defect is filed:

### 3b. A false gradient warning at the poisson limit (filed, MEDIUM)

`dev/ci082-probe.R lib nbk` under 0.3.26 Sandybridge, 1 thread
(`ob326sb-nbk.txt`): k = 15 stops at code 0 at log shape 24.479 with
"Large maximum absolute gradient at the optimum (0.00504 at
shape_(Intercept)): one Newton step ... would still gain 0.158", under
"warning" and "ignore" alike, and no limit message. k = 6 and 18 are
code 8 (the two endings above). `nbone 15` (`ob326sb-nbone15.txt`,
`ref-nbone15.txt`): the log-likelihood at 24.479 is -386.278985178424,
the poisson glm's -386.278985178388, so no step gains more than 3.6e-11;
the tape gradient in the log shape is noise past about 18 (0.0036 at
24, -0.005 at 24.48, 0.10 at 28, 218 at 35), the same on both builds.
On the reference BLAS the same data stop at 15.21 with the limit
message.

`dev/ci082-nbgrad.R` (RTMB alone, `dev/ci082-out/nbgrad.txt`) shows it
is TMB's `dnbinom_robust()` derivative, against the exact derivative:

    log shape  10.0000  tape grad  7.1714e-05  exact  7.1719e-05
    log shape  18.0000  tape grad -1.3212e-05  exact  2.3410e-08
    log shape  24.0000  tape grad -3.6231e-03  exact  5.8029e-11
    log shape  24.4791  tape grad  5.0411e-03  exact  3.5950e-11
    log shape  30.0000  tape grad -1.8418e-01  exact  1.6132e-13
    log shape  35.0000  tape grad -2.1800e+02  exact -3.2058e-15

Filed: `dev/upstream-bugs.md` TMB-7 and a MEDIUM backlog item. Not
fixed here: a stable negbinomial density touches every count family,
and a gradient verdict that measures the Newton step's gain touches the
fit-end checks of every fit; both are larger than this lane.

### A stall that is the machine's, not the lane's

With `OPENBLAS_NUM_THREADS=4` and the machine at 100 percent load (65 R
processes of all lanes), `test-se-limit.R` under OpenBLAS 0.3.26 stopped
making progress in "a family limit on a nonlinear fit still triggers
the restart" (the `gp()` nonlinear student fit, line 241): about 1900 s
of CPU in 17 minutes per process, on rellib-r21 and on the lane library
alike (`dev/ci082-out/selimit-loc-*.txt`, location reporter; last
expectation 238 in each). With `OPENBLAS_NUM_THREADS=1` the same file,
same build, passes 281 of 281 in 87.7 s. I stopped the stalled processes
(mine, matched on `ci082-run1.R`). So the se-limit matrix below runs
1 thread; 0.82.0's release ran this file under 4 threads to the end
(dev/round-20261111.md), so this is contention, not a hang in the code.

## The changed files on every BLAS setting tried (lane library)

`dev/ci082-matrix.sh lib <dir> <file>`; each log's `lib:` line reads
`wt-ci082-lib`. Settings: reference BLAS; OpenBLAS 0.3.26 and 0.3.32
with LAPACK, `OPENBLAS_CORETYPE` unset, Haswell, Zen, Sandybridge
(SkylakeX: illegal instruction on this CPU).

| file | threads | settings | each setting |
|---|---|---|---|
| test-covstruct-seam.R | 4 | 9 of 9 | pass=192 fail=0 error=0 skip=0 warn=0 |
| test-mixture-starts.R | 4 | 9 of 9 | pass=96 fail=0 error=0 skip=0 warn=0 |
| test-se-limit.R | 1 | 9 of 9 | pass=241 fail=0 error=0 skip=0 warn=0 |

(`dev/ci082-matrix1/`, `dev/ci082-matrix2/`, `dev/ci082-matrix4/`; the
summaries `dev/ci082-out/matrix{1,2,4}.txt`.) The se-limit row is the
final test file (k = 0 only); `dev/ci082-matrix3/` holds the run of the
20-perturbation version that found item 3b (Sandybridge 0.3.26 and
0.3.32: pass=280 fail=1, the other 7 settings pass=281).

## Full suites (reference BLAS, lane library, one file per process)

`NOT_CRAN=true FRMTMB_RELLIB=<lane lib> bash dev/release/run-par.sh
ci082-core dev/ci082-core.jobs 3` (log `dev/release/ci082-core.log`,
files `dev/release/ci082-core-files/`): CI082-CORE ran 330 of 330;
summed pass=24451 fail=0 err=0 skip=240 warn=0. All 330 `lib:` lines:
frmtmb 0.82.0 from `wt-ci082-lib`.

Extensions, the same way (`dev/release/ci082-ext.log`): CI082-EXT ran
199 of 199; all 199 logs load frmtmb from `wt-ci082-lib` (extensions
from rellib-r21, unchanged):

| package | files | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1851 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 13 | 697 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 91 | 4238 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 529 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 16 | 575 | 0 | 0 | 0 | 0 |

## R CMD check --as-cran, core (`dev/ci082-check.sh`)

`dev/ci082-out/check.log` (copied out; `dev/ci082-check/` deleted):
Status: 2 NOTEs. Both are the expected ones: examples over 5 s
(`frm_hazard_reads` 8.45 s, `dharma_residuals` 9.51 s elapsed, a
loaded machine) and "Skipping checking math
rendering: package 'V8' unavailable". Tests in one session, forward
order (`dev/ci082-out/check-testthat.Rout`): [ FAIL 0 | WARN 0 | SKIP
273 | PASS 24195 ].

## Core suite in one session

`dev/release/run-session.R frmtmb reverse` (`dev/ci082-out/
sess-frmtmb-reverse.txt`): SESSION frmtmb reverse files=317/330
pass=24616 fail=0 err=0 skip=240 warn=0. The 13 files without a row are
gated brms files whose file-level skip leaves no row (the same files
read pass=0 skip=0 one per process). Every LEAK line carries "[with
load of ...]" (an extension's own registrations), none is untagged.

`run-session.R frmtmb preload` (`dev/ci082-out/sess-frmtmb-preload.txt`):
SESSION frmtmb preload files=317/330 pass=24719 fail=0 err=0 skip=242
warn=0; 0 LEAK lines. The forward order is R CMD check's (above).

Extensions, `bash dev/release/run-sessions.sh ci082-extsess 2 <8
packages>` (`dev/release/ci082-extsess.log`): CI082-EXTSESS ran 24 of
24, every session fail=0 err=0 warn=0, no untagged LEAK line, all 24
`lib:` lines load frmtmb from `wt-ci082-lib`. Each package's three
orders give the one-file counts above; frmtmb.sample reads
files=86/91, its 5 files without a row being the 5 that read pass=0
skip=0 one per process.

The gated tier (`FRMTMB_BRMS_FIT_TESTS=true`) was not run. The
R-CMD-check job that failed does not set it; a gated test that reads
`confint_varcorr()` or `summary()` of a covariance block could reach the
changed rule, so the consolidation's gated run is the check of that.

## 4. `test-mixture-starts.R:182` on R-devel: a fourth RNG kind

R-devel's `RNGkind()` and `set.seed()` gained `binom.kind` (NEWS: "To
ensure reproducibility, RNGkind() gets new argument binom.kind"; kinds
"Buggy BTPE", "BTPE", read from R-devel's `src/library/base/R/RNG.R`
through `gh api repos/wch/r-source/contents/...`). `RNGkind()` returns
four values there. `mixture_starts_draws()` restored three by name and
the test compared three.

Fix: the code restores every kind by position,
`do.call(RNGkind, as.list(kinds))`, on any R version; the test sets and
compares the whole vector, and on a version with a fourth kind adds a
case with `binom.kind = "Buggy BTPE"`. Not run on R-devel: no R-devel
is installed here. On R 4.6.1 the test passes (the matrix above).

## Not done, and why

- No run on R-devel (none installed). Item 4 rests on the CI log and on
  R-devel's `RNG.R` read through `gh api`; the consolidation's CI is the
  check.
- No SkylakeX kernel (no AVX-512 on this CPU).
- `test-se-limit.R` under 4 OpenBLAS threads only to expectation 238;
  the full file ran under 1 thread on every setting (the stall above).
- No gated tier.
- Item 3b and the "ignore" question are filed, not fixed.

## For the user

1. LOW, a decision: should `check_se = "ignore"` still ask the boundary
   question for a stop with a nonzero optimizer code, so that a code-8
   stop at a limit is as quiet as a code-4 one? Recommendation: keep
   the present contract (pinned by `test-se-conv.R`): "ignore" exists
   to skip the SE analysis, and the warning it leaves is true.
2. MEDIUM, a defect: TMB's `dnbinom_robust()` derivative is noise at a
   large size (TMB-7), which can turn the limit message into a false
   "Large maximum absolute gradient" warning. Fixing it in frmtmb
   means a stable negbinomial density across the count families, a
   lane of its own.

## Fix round (review `dev/reviews/2026-11-13-ci082.md`, verdict FIX)

The review's scripts `dev/reviews/ci082-review-{vc,proto,nb8}.R` were
rerun on the fixed lane build; outputs `dev/ci082-out/fix-*.txt`.
Before = the first-round lane build (review logs `lib-*.txt`) or 0.82.0
(`base-*.txt`); after = this build.

### FIX 1 (MAJOR): the fixed-row rule, option (a)

The 16-eps derivative test is gone. `varcorr_fixed_rows(bk, t0)`
(`R/confint.R`) reads each row, unclamped (log sd, atanh correlation),
at the structure's `start(bk)` and at `t0` moved by -1 and +1 in each
theta coordinate; a row is fixed when every point gives the value at
`t0` to `1e-10` of the larger magnitude, no absolute floor; an error or
a non-finite value at a point counts as a move. No covstruct slot was
added. `se_g[fixed] <- 0` now runs right after the delta method and
BEFORE the nocurv and `jc_nonest()` NA rules, so a lost-direction NA
always stands; for a row not fixed every line is 0.82.0's.

| case | 0.82.0 | first round | fix round |
|---|---|---|---|
| xcs 0.5, 20 perturbations, reference | 11 NA + warning, 9 width ~1e-12 | 20 zero width | 20 zero width (`fix-ref-lib-xcsk.txt`: se 0 on 20 of 20, no NA) |
| fixed 0.999, 20 perturbations (`vc fix999`) | 6 NA + warning, 14 width 3.0e-12 to 9.1e-12 | 6 zero width, 14 width 3e-12 to 9e-12 | 20 zero width, se 0, no warning |
| ou independence, seeds 1-3 (`vc ouind`) | na 11/11/10, positive tiny se 8/8/9, 1 warning each | 19 se 0 each, no warning | identical to 0.82.0 |
| homcs started at -20, -26, -28, -30, -34 (`vc csstart`) | NA + warning at -20 to -30; 3 tiny positive at -34 | se 0 at -28, -30, -34 | identical to 0.82.0 |
| cs / homcs natural fits, seeds 1-3 (`vc cslow`) | all interior | the same | identical to 0.82.0 |

"Identical" is `diff` of the whole log against the review's
`base-*.txt` with the `lib:` line and arm name dropped: no difference
for ouind, csstart and cslow (every printed se and width). fix999
differs from 0.82.0 in every k, as intended.

Interior `us()` rows unchanged: the cost script below asserts
`identical()` output of `varcorr_trans_rows()` with the probe and with
the probe replaced by "nothing fixed" on a d = 6 us() fit ("same rows
with and without the probe: TRUE"), and the review's prototype found
no us row fixed (`proto2.txt`).

New test `test-covstruct-seam.R` "a fixed correlation has no width on
every rounding" now covers xcs at 0.5 and a registered `xch` at 0.999
(20 perturbations each, `identical(lwr, estimate)`), the estimated
`(0 + t | g)` block, and the ou independence fit (seed 1): every cor
row has se NA or > 0, at least one > 0, every finite interval positive
width, every row whose estimate is exactly 0 NA, and the "No interval"
warning required. On rellib-r21, reference BLAS
(`dev/ci082-base2/covstruct-ref.txt`): pass=177 fail=80 warn=17 (40
failures each at lines 161 and 162, 17 escaped "No interval"
warnings: 11 from xcs, 6 from xch). On the lane build
(`dev/ci082-out/fix-covstruct-ref.txt`): pass=257 fail=0 warn=0.

Cost (`dev/ci082-cost.R`, log `dev/ci082-out/cost.txt`; us() block
d = 6, 21 theta, 60 groups; interleaved arms, blocks of 128 calls, min
of 15 rounds): `varcorr_trans_rows()` 8.52 ms without the probe, 11.64
ms with it (ratio 1.367; control 0.987). The probe adds 44 of the 88
structure `vcov()` calls (1 at t0, 1 at the start, 42 at +/-1).
`confint_varcorr()` of that fit: 0.013 s per call.

### FIX 2 (MINOR): the "ignore" assertion of test-se-limit.R

The nonzero-code branch now also requires "false convergence (8)",
"so the estimates are at the maximum", and coefficients equal to the
"warning" fit's (ratio within 1e-6 of 1). `dev/ci082-nb8.R`
(`fix-nb8.txt`, reference BLAS, the assertion as the file has it):
k = 0 (code 0), k = 3 and 4 (code 8, "gains 9.9e-07 [3.51e-06] ... so
the estimates are at the maximum"): pass. The short-stop arm (nlminb
held to 3 iterations, 6 evaluations) fails on all three k: code 9,
"A Newton step from the stop would still gain 0.0657 ... a point short
of the maximum" plus a gradient warning; it fails at the warning count,
and its first warning carries neither of the two new phrases.

### FIX 3 (MINOR): RTMB-4, not TMB-7

The TMB-7 entry is deleted; its numbers and the consequence are added
to RTMB-4 (`dev/upstream-bugs.md`); the MEDIUM backlog item points at
RTMB-4. The part of `upstream-bugs.md` after RTMB-4 is byte-identical
to ec8178d3's.

### Optional: plain and scaled fits at one maximum

`test-mixture-starts.R` (autoscale test) also asserts that the plain
and scaled fits' log-likelihoods differ by less than 1e-7 of the
log-likelihood (3.9e-5): all 240 probe fits read -388.095489 to 6
decimals, so they agree within 1e-6. `oms_fit()` and `step =` are
untouched.

### Fix-round runs (lane library, final tree)

Matrix `bash dev/ci082-matrix2.sh lib dev/ci082-fixmx dev/ci082-fix.jobs`
(summary `dev/ci082-out/fixmx.txt`): 18 files x 7 settings = 126 runs,
every `lib:` line `wt-ci082-lib`, every `R:` line the build named.
Settings: reference BLAS; OpenBLAS 0.3.26 and 0.3.32 with LAPACK
(rebuilt with `dev/ciharden-openblas.sh`, logs `dev/ci082-ob326b.log`,
`dev/ci082-ob332b.log`), core type unset, Haswell, Sandybridge; 4
threads, `test-se-limit.R` 1 thread. The same counts on all 7 settings:

| package | file | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | test-covstruct-seam.R | 257 | 0 | 0 | 0 | 0 |
| frmtmb | test-mixture-starts.R | 97 | 0 | 0 | 0 | 0 |
| frmtmb | test-se-limit.R | 241 | 0 | 0 | 0 | 0 |
| frmtmb | test-se-conv.R | 55 | 0 | 0 | 0 | 0 |
| frmtmb | test-confint-anova.R | 35 | 0 | 0 | 0 | 0 |
| frmtmb | test-review-v29.R | 187 | 0 | 0 | 0 | 0 |
| frmtmb | test-v07.R | 27 | 0 | 0 | 0 | 0 |
| frmtmb | test-v18.R | 27 | 0 | 0 | 0 | 0 |
| frmtmb | test-multiple-pooling.R | 33 | 0 | 0 | 0 | 0 |
| frmtmb | test-double-bar.R | 154 | 0 | 0 | 0 | 0 |
| frmtmb | test-covstruct-floor.R | 935 | 0 | 0 | 0 | 0 |
| frmtmb | test-covstruct.R | 7 | 0 | 0 | 0 | 0 |
| frmtmb.spatial | test-car.R | 132 | 0 | 0 | 0 | 0 |
| frmtmb.spatial | test-guardrail.R | 148 | 0 | 0 | 0 | 0 |
| frmtmb.spatial | test-sar.R | 157 | 0 | 0 | 4 | 0 |
| frmtmb.spatial | test-spde.R | 48 | 0 | 0 | 0 | 0 |
| frmtmb.latent | test-hmm.R | 120 | 0 | 0 | 0 | 0 |
| frmtmb.latent | test-scale.R | 0 | 0 | 0 | 2 | 0 |

Full core suite, one file per process, reference BLAS
(`NOT_CRAN=true FRMTMB_RELLIB=<lane lib> bash dev/release/run-par.sh
ci082-core2 dev/ci082-core.jobs 2`, log `dev/release/ci082-core2.log`):
CI082-CORE2 ran 330 of 330; pass=24517 fail=0 err=0 skip=240 warn=0
(first round 24451: +65 in test-covstruct-seam.R, +1 in
test-mixture-starts.R); all 330 `lib:` lines `wt-ci082-lib`.

R CMD check --as-cran, core, final tree (`dev/ci082-check.sh`; logs
copied to `dev/ci082-out/check2.log`, `check2-testthat.Rout`;
`dev/ci082-check/` deleted): Status: 1 NOTE ("Skipping checking math
rendering: package 'V8' unavailable"). Tests in one session, forward
order: [ FAIL 0 | WARN 0 | SKIP 273 | PASS 24261 ]. A first attempt
failed at `R CMD build` ("copying to build directory failed") because
I deleted `dev/ciharden-out/` while the build was copying the tree;
the copies are deleted and the check was rerun from the start.

Not rerun in the fix round: the one-session runs and the extension
suites in full (only the extension covstruct files above). The fix
changes only which rows `varcorr_fixed_rows()` calls fixed, and for a
row it does not call fixed the code is 0.82.0's, line for line.

## Fix round 2 (short check `dev/reviews/2026-11-13-ci082-short.md`, FIX on one item)

### 1. An underflowed correlation called fixed

The short check's `oufar` (ou() on independent data, times 0, 0.5, 1,
2, 8, 20, 40, 1000, seeds 1 to 3): every correlation at a lag of 960 or
more is exactly 0 at the estimate, at the start and at both +/- 1
points, so the probe called those 7 rows per seed fixed. Before (fix
round 1, `dev/reviews/ci082-short-out/atk-lib-oufar.txt`): 7 rows se 0
and zero width per seed; seed 2 lost its warning (0.82.0: 7 NA + 1
warning). 0.82.0 (`atk-base-oufar.txt`): 17 / 7 / 18 NA, 1 warning
each.

Change (`varcorr_fixed_rows()`): a CORRELATION row at exactly 0 at
`t0` is never fixed. The coordinator's line was `fixed <- is.finite(g0)
& g0 != 0`; I applied the `!= 0` to correlation rows only, because a
log sd of exactly 0 is an sd of exactly 1, which no underflow
produces, and a structure fixing a unit sd would otherwise be sent to
the zero-se NA. For correlations the two forms are the same. The row's
jacobian row is exactly 0 there, so it keeps 0.82.0's NA and warning
deterministically.

After (`dev/ci082-out/fix3-lib-oufar.txt`, reference BLAS): 17 / 7 / 18
NA, 0 fixed, 1 warning each; the whole log, every printed se and
width, equals 0.82.0's apart from the probe column (absent on 0.82.0)
and trailing spaces.

Test: "a fixed correlation has no width on every rounding" now runs
the ou independence fit on two grids, the 7 times (seed 1) and the 8
times with 1000 (seeds 1 to 3), and requires on each: the "No interval"
warning, every cor se NA or > 0, at least one > 0, every finite
interval of positive width, the estimate-0 rows NA, and every
`tim(1000)` row NA. On the fix-round-1 build
(`dev/ci082-out/fix2-pre-covstruct.txt`): pass=263 fail=13 (the missing
warning for seed 2 at line 197, and lines 201, 203, 205, 207 on each
of the 3 far seeds). On this build: pass=276 fail=0 warn=0. On
rellib-r21 (`dev/ci082-base2/covstruct-ref-fix3.txt`): pass=196 fail=80
warn=17, the 80 being xcs and xch (lines 161, 162); the ou cases pass
there, as they must (0.82.0's behavior is the target).

### 2. Docs: zero width only inside the clamp

Decided: keep 0.82.0 at the clamp. NEWS and `?confint_varcorr` now say
that for a correlation the zero-width interval holds when the fixed
value is not 0 and lies inside -0.9999 to 0.9999, and that a fixed
correlation at or past that limit keeps NA bounds and the warning
(the boundary rule `bd`). `man/confint_varcorr.Rd` regenerated.

### 3. Recorded: a vcov read back through two inverses (not a blocker)

Short check `attack invgauss` (a registered structure with vcov
`solve(D^-1 solve(R) D^-1)`, R a gaussian kernel `exp(-lag^2 / 2)`,
d = 10, 10 ulp perturbations): 35 of 45 correlations fixed (zero
width), the 10 far rows (down to 2.6e-18) carry a read-back error near
eps, a relative move up to 1.96, so they fall back to 0.82.0's rule: a
positive width up to 2.0e-12 on the lane build; on 0.82.0 all 45
positive (up to 1.7e-11) and at k = 4 one row NA with the warning. So
on those rows the outcome is still rounding's, as on 0.82.0. Only a
registration slot naming the fixed rows (the review's option (b))
would make it exact.

### Runs

Reference BLAS, lane build, the review and short-check scripts rerun:
`ci082-short-attack.R lib oufar` as above; `ci082-review-vc.R lib`
`ouind`, `csstart`, `cslow` byte-identical to 0.82.0's `base-*.txt`
(`lib:` line dropped, arm normalized); `fix999` byte-identical to fix
round 1 (`fix-lib-fix999.txt`: 20 of 20 se 0, zero width, no warning).

Test files (`dev/ci082-fix3mx/`, every `lib:` line `wt-ci082-lib`;
OpenBLAS 0.3.32 + LAPACK rebuilt (`dev/ci082-ob332c.log`), 4 threads,
`OPENBLAS_CORETYPE=Sandybridge`, deleted afterwards):

| file | reference BLAS | 0.3.32 Sandybridge |
|---|---|---|
| test-covstruct-seam.R | 276/0/0/0 | 276/0/0/0 |
| test-covstruct-floor.R | 935/0/0/0 | 935/0/0/0 |
| test-se-conv.R | 55/0/0/0 | 55/0/0/0 |

(pass/fail/error/skip; warn 0, exit 0.)

R CMD check --as-cran, core, final tree (`dev/ci082-check.sh`; logs
`dev/ci082-out/check3.log`, `check3-testthat.Rout`; `dev/ci082-check/`
deleted): Status: 2 NOTEs, the expected ones (examples over 5 s on a
loaded machine; V8 unavailable). Tests in one session: [ FAIL 0 |
WARN 0 | SKIP 273 | PASS 24280 ].

Machine notice during this round (memory commit limit, disk): from then
on at most 2 R processes; only the check ran after it. No run failed on
space or allocation.
