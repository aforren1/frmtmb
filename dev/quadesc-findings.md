# Lane quadesc (round of 2026-11-16): the rule's stationary escape

Base b255350e (frmtmb 0.83.0, rellib-r22). Worktree
`frmtmb-wt-quadesc`, library `wt-quadesc-lib`. This file is the record
if the report does not arrive.

## Summary

- The CI failure (ubuntu-next, 26.04) and item 5015 are one defect. A
  restart of the rule's stationary escape runs alpha toward its end,
  the rule's gradient turns NaN, `optimize_obj()` raises, and
  `escape_stationary()` dropped the restart. Whether the optimizer
  stopped before the NaN was a rounding branch: 0.83.0 stopped on
  alpha = 0 on 3 to 10 of 21 roundings per case and BLAS, the
  reference BLAS included (section 4).
- The CI's exact configuration (Cooperlake kernel, gcc 15) cannot run
  here (AVX-512); its ending reproduces locally under every BLAS at
  other roundings, e.g. k = 9 (sections 1, 2).
- Fix (`R/fit.R`): an abandoned restart keeps its lowest
  finite-gradient point (`escape_track()`), used only when no restart
  that ended on its own gained. On 189 fits (3 cases, 21 roundings, 3
  BLAS builds) none stops on the point; per rounding the lane is never
  below 0.83.0 in exact marginal log-likelihood and is above it exactly
  where 0.83.0 stopped (section 5).
- Test rewritten to assert what a user sees, on a rounding where 0.83.0
  fails under all three BLAS builds (section 7).
- Filed: the rule's log-likelihood at a run-off alpha reads 0.27 to
  0.54 above the exact value (MEDIUM); the LOW of item 5015 refiled;
  the NaN gradient of skew_normal()'s density below alpha * z of about
  -3e9, from RTMB's `pnorm(log.p = TRUE)` (MEDIUM, corrected in fix
  round 1; RTMB-3 in `dev/upstream-bugs.md` extended).
- Fix round 1 (section 10): a restart counts as a gain only when it
  left the point's window; an abandoned point carries `evals` and
  `nonfinite_trials`; the test adds a mock of the veto and the CI's
  k = 0.
- Bump: patch. Floors: none change (no extension needs this core).

## Instruments

- `dev/quadesc-run.sh <ver> <coretype|default> <log> <script> [args]`:
  <ver> is ref, 0.3.26 or 0.3.32;
  one R process under the reference BLAS or an OpenBLAS build of
  `dev/ciharden-openblas.sh <ver> lapack`; log `dev/quadesc-log/<log>.txt`.
- `dev/quadesc-probe.R <lib|base> <case> [k...]`: the fit of each case,
  every `escape_stationary()` call logged. Cases: `alphax` (the fit of
  `test-quad-escape-record.R`, seed 1011), `symsmall` (SYMSMALL seed 1 of
  `dev/reviews/quadseam2-review-esc2.R`), `weak` (WEAK seed 4). k scales
  the response by `1 + 2e-15 * k`.
- `dev/quadesc-cands.R`: the same, with every calibration candidate of
  `quad_fit()` and every `quad_settle()` call logged.
- `dev/quadesc-restarts.R`: inside the rule's escape, each start
  optimized on the frozen rule and on a rule re-anchored at the start,
  each endpoint scored by the exact marginal log-likelihood (one-dimensional
  integrals by `stats::integrate()`, no rule).

## 1. The CI failure (run 43878491850, ubuntu-26.04)

The runner: Intel Xeon 6975P-C, OpenBLAS 0.3.32+ds-5 (Ubuntu amd64v3
build) with LAPACK, kernel `Core: Cooperlake`, gcc 15.2.0. The failing
fit stopped with alpha_x at 3.862e-15 and the warning "alpha_x ... where
the log-likelihood rises by a further 4"; the record kept was the Laplace
one (`identical(esc, lap)`).

Local, at k = 0 (`dev/quadesc-log/cands-base-*-alphax.txt`), base build:

    reference BLAS          escapes, ll -84.5783086103
    OpenBLAS 0.3.32 default escapes, ll -84.5602809547
    OpenBLAS 0.3.32 Haswell escapes, ll -84.5602809547
    OpenBLAS 0.3.32 Zen     escapes, ll -84.5602809547
    OpenBLAS 0.3.32 SkylakeX, Cooperlake: exit 132 (illegal instruction;
      this CPU, a Core Ultra 9 285K, has no AVX-512)

So the CI failure does not reproduce at k = 0 on any core type this
machine can run. It is a rounding branch: see section 2.

## 2. The cause (item 5015 and the CI failure are one defect)

Base build, k = 0..20 on ALPHAX (`dev/quadesc-log/ulp-base-*-alphax.txt`,
`dev/quadesc-probe.R base alphax 0:20`): the fit stops on alpha = 0
(logLik -91.2452361629, the record kept is the Laplace one, the CI's
failure) at

    reference BLAS   k = 1, 7, 9, 16          (4 of 21)
    OpenBLAS 0.3.32  k = 4, 9, 13, 14, 15, 17, 19   (7 of 21)

So the CI's ending is a rounding branch that the local reference BLAS
takes too. The test at k = 0 passed here by rounding.

What happens (`dev/quadesc-cands.R`, `dev/quadesc-restarts.R`; logs
`cands-base-*`, `rs-base-*`, `rs2-lib-*`):

1. `quad_fit()` escapes on the Laplace objective (ALPHAX: alpha from
   1e-15 to 4.07, -8.53; objective 91.245 to 86.224).
2. The rule anchored there is optimized; alpha runs off (to 1.4e6), the
   tape's gradient fails, and the re-anchor at that point does not tape
   (`CAND retry`, then `CAND NULL`).
3. The next calibration point is the cold template (alpha = 0). Its
   optimum is on the stall and is "stationary" (gradient 1.7e-4 below
   `grad_tol`), so `quad_keep_best()` keeps it outright (objective
   91.245, 5.02 worse than the Laplace escape's point).
4. The rule's escape in `fit_assembled()` restarts from alpha = +-2 on
   that rule. Each restart runs alpha toward its end, where the rule's
   gradient becomes NaN ("NA/NaN gradient evaluation"); `optimize_obj()`
   restarts once from the tape's best point, whose gradient is NaN too,
   and raises; `escape_stationary()` drops the restart. Whether the
   optimizer stops (false convergence, code 1) before it reaches a NaN
   gradient is the rounding branch.

In the failing run (reference BLAS, k = 1, `rs2-lib-ref-alphax-k1.txt`)
the tape's best point after the first restart is alpha (2404321,
-6929248), objective 84.560031, exact marginal log-likelihood -85.105051,
gradient not finite. The stall's exact value is -91.245236. On SYMSMALL
under OpenBLAS 0.3.32 (`rs2-lib-ob332-symsmall-k0.txt`): best point
alpha 13388137, objective 84.925389, exact -85.191921, gradient not
finite.

The item's inference ("the restarts run on the rule frozen at alpha =
0") is half right: they do, but the frozen rule is not what sends them
back. They raise.

WITHDRAWN (fix round 1): this section said RTMB's `pnorm(log.p = TRUE)`
was not the source, because its value and derivative are exact to
x = -1000 on the tape (`dev/quadesc-pnorm.R`). That probe stopped far
too early. The review showed it IS the source (section 10): the
derivative is non-finite below about -3.1e9, which the rule's far nodes
reach. The fact was already in `dev/upstream-bugs.md` as RTMB-3 (NaN
near -4.46e9, lane optima), and I did not search for it there.

Correction to step 4 for SYMSMALL under OpenBLAS 0.3.32 (k = 0): the
restart from alpha = +2 raised; the one from -2 came back to alpha =
-0.0017 (objective unchanged). So one start raising is enough when the
other start's side holds no gain.

## 3. The change

`R/fit.R`: `escape_stationary()` optimizes each restart on a wrapper of
the objective, `escape_track()`, which remembers the lowest point where
the objective and its gradient were both finite (read at the points
where nlminb asks for the gradient after reading the objective). When
`optimize_obj()` raises, the restart returns that point with code 1 and
the message "the optimizer stopped on an error (...) on its way from a
restart off a stationary point; the lowest point with a finite gradient
is kept", so the fit warns `frmtmb_conv_code`. Before, the restart was
dropped (`error = function(e) NULL`).

Why not the tape's own best point (`last.par.best`, the first build of
this change): it is where the gradient first failed. In
`rs2-lib-ref-alphax-k1.txt` and `rs2-lib-ob332-symsmall-k0.txt` the
best point's gradient is not finite in all six restarts, so a salvage
of it gave nothing (that build's sweep still stopped on alpha = 0 at
k = 1, 7 and 9 under the reference BLAS, `ulp-lib-ref-alphax.txt`).

Not done, and why:

- Re-anchoring the rule at each restart's start (`quad_fit()`'s
  `restart`): measured in `rs-base-*`; the re-anchored restarts ran off
  too and returned NULL where the frozen rule raised (ALPHAX k = 1: all
  three starts NULL). It does not address the error.
- Changing `quad_keep_best()` so that the cold template's stationary
  candidate does not beat the Laplace escape's point: the Laplace-
  anchored rule itself runs off and its re-anchor does not tape
  (`CAND NULL`), so there is no candidate from that point to keep. The
  escape on the frozen rule is what has to work.

## Permission denied (reported, not worked around)

`git diff --no-index` of a `git show HEAD:R/fit.R` copy against
`R/fit.R` was denied by the permission system. No other diff tool was
tried. The lane's hunks are listed in section 3 and in the report.

## 4. Sweep: 0.83.0 against the lane's first build

k = 0..20 (response times 1 + 2e-15 k), three cases, three BLAS builds:
`dev/quadesc-sweep.sh`, `dev/quadesc-sweep-b.sh` (logs
`dev/quadesc-log/sw-a-<arm>-<ver>-<case>.txt`), summary
`dev/quadesc-sweep-sum.R a` (`sweep-sum-a.txt`); 0.83.0's ALPHAX under
the reference BLAS and 0.3.32 from `ulp-base-*-alphax.txt` (same fits,
no exact score). "On point": every alpha coefficient below 1 in
absolute value at the end. Exact: the exact marginal log-likelihood
(trapezoid rule, 20001 points over +-12 group sds; matches the
`stats::integrate()` score of `dev/quadesc-restarts.R` to 1e-6, e.g.
-85.961442 at SYMSMALL's Laplace point).

    0.83.0, fits on alpha's point (of 21):
      ALPHAX   ref 4  (k 1,7,9,16)   0.3.32 7 (k 4,9,13,14,15,17,19)
               0.3.26 7 (same k as 0.3.32)
      SYMSMALL ref 10 (k 3,4,6,8,11,12,16,17,18,20)
               0.3.32 10, 0.3.26 10 (k 0,1,2,4,7,9,13,15,19,20)
      WEAK     ref 3 (k 2,13,19)   0.3.32 3, 0.3.26 3 (k 0,14,18)
    first lane build: 0 of 21 on every case and BLAS (189 fits)

The stall's exact value and the lane's endings (first build, all
BLAS): ALPHAX -91.2452 against -85.1032 to -85.0361; SYMSMALL -87.2168
against -85.1956 to -85.1926; WEAK -75.5907 against -74.5794 to
-74.5791. In rule log-likelihood the gains are those of the item: 6.67
to 6.69, 2.29, 1.45.

The Laplace fit from the same start, scored exactly: ALPHAX -86.009814,
SYMSMALL -85.961442, WEAK -75.420328, the same on every k and BLAS. The
lane's fits are above it by at least 0.9066, 0.7658 and 0.8409.

A regression the first build had: on ALPHAX, where 0.83.0 escaped by a
restart that ended on its own at alpha about (5920, -16445) and one that
raised, the first build kept the raised restart's point (rule objective
0.48 lower, alpha about (2e6, -6e6)), whose exact log-likelihood is
0.0633 to 0.1234 LOWER than 0.83.0's: 12 of 21 roundings under 0.3.26
(per-k comparison in the lane's notes from `sw-a-base-0.3.26-alphax.txt`
and `sw-a-lib-0.3.26-alphax.txt`; e.g. k = 1: 0.83.0 rule ll -85.0393,
lane -84.5601, exact 0.0633 the other way). The rule anchored at alpha
= 0 overstates the likelihood far from its anchor (section 6). On
SYMSMALL and WEAK the lane was below 0.83.0 by more than 1e-3 (exact) on
0 of 21 under every BLAS.

So the second build keeps an abandoned restart's point only when no
restart that ended on its own gained (`escape_stationary()`).

## 5. Sweep: the lane's second build (the change as delivered)

`dev/quadesc-sweep-c.sh` (logs `sw-c-lib-<ver>-<case>.txt`), summary
`dev/quadesc-sweep-sum.R c` (`sweep-sum-c.txt`), per-k exact comparison
against 0.83.0's `sw-a-base-*` logs:

    fits on alpha's point: 0 of 21 on every case and BLAS (189 fits)
    smallest |alpha coefficient| over the 189 fits: 5917.7
    exact log-likelihood, lane minus 0.83.0, per k:
      below by more than 1e-3: 0 on every case and BLAS compared
      above by more than 1e-3: SYMSMALL 10, 10, 10 (ref, 0.3.32, 0.3.26);
        WEAK 3, 3, 3; ALPHAX 7 (0.3.26); exactly the k where 0.83.0
        stopped on the point; the other k are identical
      (ALPHAX ref and 0.3.32 have no 0.83.0 exact score; there the
      0.83.0 fits on the point are the k listed in section 4 and the lane
      escapes on every k)
    exact range: ALPHAX -85.1020 to -84.9742, SYMSMALL -85.1956 to
      -85.1926, WEAK -74.5797 to -74.5790
    over the Laplace fit's exact value (ALPHAX -86.009814, SYMSMALL
      -85.961442, WEAK -75.420328): at least 0.9078, 0.7658 and 0.8406

The gains of item 5015 under each BLAS, rule log-likelihood (0.83.0 ->
lane at the k where 0.83.0 stopped): SYMSMALL 2.29 (-87.2168 ->
-84.925), WEAK 1.45 (-75.5907 -> -74.138), ALPHAX 6.69 (-91.2452 ->
-84.560); in exact log-likelihood 2.02, 1.01 and 6.14. At k = 0 the
item's cases: SYMSMALL and WEAK stopped on the point under 0.3.26 and
0.3.32 (as the item says) and escape on the lane.

Warnings on the lane's 189 fits (classes): conv_code 98, conv_gradient
42, quad_unsettled 35, frmtmb_warning (standard errors) 5. The
unsettled and gradient warnings are on the ALPHAX endings near (5920,
-16445) that 0.83.0 reaches too.

## 6. Accuracy of a run-off fit (filed, not fixed)

At alpha near 1e6 the estimate is read on the rule anchored at the
stall (`quad_settle()` cannot re-anchor: the rule taped there is not
finite, `SETTLE ... rounds NULL` in the cands logs). It reads high
against the exact score: ALPHAX -84.5601 against -85.1009 (0.54,
`sw-a-lib-ref-alphax.txt` k = 0), SYMSMALL -84.9254 against -85.1952
(0.27, k = 0), WEAK -74.1378 against -74.5794 (0.44). 0.83.0's own
escaped fits read the same way (SYMSMALL ref k = 0: -84.9253951576,
exact -85.195206). Where an ending near alpha 5920 is re-anchored, the
rule and the exact score agree to 2e-3 (ALPHAX k = 10: -85.0393 against
-85.0375). Filed MEDIUM ("Filed in the round of 2026-11-16").

## 7. The test, and whether ubuntu-next would pass

`tests/testthat/test-quad-escape-record.R`: the record test is kept; the
second test ("the rule's escape runs after the Laplace escape ran") is
replaced by two:

- "a restart that raises keeps its last finite-gradient point": a unit
  test of `escape_track()` on a mock objective (gradient NaN beyond a =
  5): the kept point is a = 4, code 1; a gradient read at a point other
  than the last objective's is not remembered; nothing remembered gives
  NULL.
- "the rule's escape takes alpha off its point on every rounding": the
  ALPHAX fit with the response times 1 + 2e-15 * 9, a rounding on which
  0.83.0 stops on alpha = 0 under the reference BLAS, 0.3.26 and 0.3.32.
  It asserts what a user sees, and accepts each valid record: the
  Laplace record is kept with a positive gain; the kept record is the
  Laplace one or a rule record with a positive gain; every alpha
  coefficient above 1 in absolute value; and the fit's exact marginal
  log-likelihood above that of the Laplace fit from the same start (the
  higher optimum found by other means; trapezoid rule in the test).
  Allowed warnings: the optimizer's code, a large gradient, an unsettled
  rule, standard errors not available: the four that the 63 lane fits
  of this design give (section 5). No absolute tolerance.

Seen failing on 0.83.0 (rellib-r22), `dev/quadesc-tests.sh <ver> base`,
logs `t-base2-<ver>-test-quad-escape-record.txt`, under each of ref,
0.3.32, 0.3.26: pass=7 fail=2 err=1 (escape_track missing; `all(abs(a) >
1)` FALSE; exact -91.2 <= -86.0). On the lane, same three builds
(`t-lib2-<ver>-...`): pass=14 fail=0 err=0 skip=0 warn=0.

Every rounding: the lane's 63 ALPHAX fits of section 5 (k = 0..20, three
builds) satisfy each assertion: 0 on alpha's point, smallest |alpha|
5917.7, exact at least 0.9078 above the Laplace fit's, every record a
rule record with gain > 0 or the Laplace one (checked from the FIT
lines), warnings only of the four allowed classes.

ubuntu-next: its two failures (FAIL 2, both in this test, `:61` and
`:64`) are the ending this test now pins. The runner's kernel
(Cooperlake, AVX-512) and gcc 15 cannot be run here (SkylakeX and
Cooperlake exit 132), so its exact rounding is not reproduced. On the
evidence, 189 lane fits over three cases, 21 roundings and three BLAS
builds, none stops on the point, and the new test passes on every
rounding of its design measured, so I expect ubuntu-next to pass this
file. That is an inference, not a run: only a workflow run on 26.04
shows it, and this lane changes no workflow and starts none.

## 8. Regression checks and suites

The reviewer's 45-fit set (`dev/quadesc-esc2.R`, a copy of
`dev/reviews/quadseam2-review-esc2.R` with this round's libraries; logs
`esc2-*.txt`), lane against 0.83.0, logLik per fit:

    reference BLAS:  45 fits, lane below by > 1e-8: 0, above: 0
    OpenBLAS 0.3.32: 45 fits, lane below by > 1e-8: 0, above: 2
                     (SYMSMALL 1: +2.291, WEAK 4: +1.453; item 5015)

OpenBLAS with LAPACK (`dev/quadesc-tests.sh <ver> lib`, logs
`t-ob-<ver>-*.txt`), each fail 0 err 0 skip 0 warn 0 under both 0.3.26
and 0.3.32: test-quad-escape-record.R 14, test-stationary-passes.R 4,
test-quadrature-defects.R 89, test-skew-normal-start.R 70,
test-optima.R 107, test-quadrature-nested.R 69, frmtmb.latent
test-hmm-stationary.R 12.

Every suite, one file per process, two at a time (`dev/release/run-par.sh
quadesc-all dev/quadesc-all.jobs 2`, `NOT_CRAN=true
FRMTMB_BRMS_FIT_TESTS=true`, `FRMTMB_RELLIB` = the lane library, so each
extension of rellib-r22 runs on the lane's core; every per-file log
names `wt-quadesc-lib`): `QUADESC-ALL ran 549 of 549`.

    frmtmb           349 files pass 28762 fail 0 err 0 skip 15 warn 0
    frmtmb.coupling   11 files pass   551 fail 0 err 0 skip  5 warn 0
    frmtmb.eam        33 files pass  1857 fail 0 err 0 skip  3 warn 0
    frmtmb.latent     13 files pass   700 fail 0 err 0 skip  2 warn 0
    frmtmb.learn      16 files pass   545 fail 0 err 0 skip  2 warn 0
    frmtmb.ode        11 files pass   549 fail 0 err 0 skip  1 warn 0
    frmtmb.sample     92 files pass  4651 fail 0 err 0 skip  1 warn 0
    frmtmb.spatial     8 files pass   634 fail 0 err 0 skip  0 warn 0
    frmtmb.survival   16 files pass   575 fail 0 err 0 skip  0 warn 0

The skips are the opt-in scale tier (`FRMTMB_SCALE_TESTS`) and the
drmTMB agreement tier. The drmTMB tier, which fits a skew normal, run on
its own with `FRMTMB_DRMTMB_FIT_TESTS=true`
(`t-drm-test-drmtmb-agreement.txt`): pass=131 fail=0 err=0 skip=0
warn=0.

One session per suite, three orders (`dev/release/run-sessions.sh
quadesc-sess 2 frmtmb frmtmb.latent frmtmb.sample`, the suites whose
fits reach the escape): `QUADESC-SESS ran 9 of 9`, fail 0 err 0 warn 0
in each; frmtmb forward pass 28762 skip 15, reverse 28927 skip 15,
preload 29031 skip 17; frmtmb.latent 700 each, skip 2; frmtmb.sample
4651 each, skip 1. LEAK lines without a load tag or the extension-state
tag: 0.

Why only those three suites in one session: `escape_stationary()` acts
only when a family declares `post$stationary` and the fit lands on it;
only `skew_normal()` declares one (and mixtures and hmm states built
from it). No extension's R code names the escape, and only
frmtmb.latent and frmtmb.sample's tests name `skew_normal` (grep of
`extensions/*/R` and `extensions/*/tests`). All nine suites ran one file
per process above.

## 9. R CMD check

`dev/quadesc-check.sh` (build and `--as-cran` check inside
`dev/quadesc-check/`, R_LIBS = lane library, rellib-r22, user library;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; pandoc and TinyTeX on PATH):
frmtmb `Status: 1 NOTE`, the V8 math-rendering NOTE of the HTML
manual; tests 48 min, OK. Only core changed, so no extension was
checked. `dev/quadesc-check/` and `dev/ciharden-out/` were deleted
afterwards.

## Files touched

(As of fix round 1.)

- `R/fit.R`: `escape_stationary()` (the restart loop: `escape_track()`
  wrapper; per-restart tally; abandoned restarts kept only when no
  restart that ended on its own gained and left the point; their
  counts); new `escape_track()` after it. No other hunk.
- `tests/testthat/test-quad-escape-record.R`: second test replaced by
  three (sections 7 and 10); the record test unchanged.
- `NEWS.md`: new top heading `# frmtmb (development version)`, section
  "Convergence and standard errors", one BREAKING bullet.
- `dev/test-backlog.md`: CLOSED line under item 5015; section "Filed in
  the round of 2026-11-16", subsection "From lane quadesc" (2 MEDIUM, 1
  LOW).
- `dev/upstream-bugs.md`: RTMB-3, one "Also" paragraph and its Record
  line.
- New under `dev/`: `quadesc-findings.md`, `quadesc-run.sh`,
  `quadesc-probe.R`, `quadesc-cands.R`, `quadesc-restarts.R`,
  `quadesc-pnorm.R`, `quadesc-sweep.sh`, `quadesc-sweep-b.sh`,
  `quadesc-sweep-c.sh`, `quadesc-sweep-d.sh`, `quadesc-sweep-sum.R`,
  `quadesc-tests.sh`, `quadesc-install.sh`, `quadesc-check.sh`,
  `quadesc-esc2.R`, `quadesc-veto-ulp.R`, `quadesc-fix3-seen.R`,
  `quadesc-core.jobs`, `quadesc-ext.jobs`, `quadesc-all.jobs`. Logs in
  `dev/quadesc-log/` and `dev/release/quadesc-*.log` (gitignored).

No compat row, brms-suite ledger row, workflow file or extension
changed. Roxygen output unchanged (the new function is `@noRd`).

## 10. Fix round 1 (review `dev/reviews/2026-11-17-quadesc.md`, FIX)

### The cause, corrected (review finding 1)

RTMB 2.0's `pnorm(log.p = TRUE)` derivative is the NaN source
(`dev/reviews/quadesc-review-iso2.R`, log `-iso2-ref.txt`): non-finite
at 174 of 4001 log-spaced points of x in [-1e14, -1e2], all in
[-9.93e13, -3.10e9]; relative error against -x 3.7e-5 at -1e6, 3.2e-3
at -1e7, 0.34 at -1e8. A plain RTMB replica of ALPHAX at k = 9 confines
the NaN to group 11, whose rule value (6.0384) is within 0.022 of its
exact value (6.0606): the value is sound, only the derivative fails.
There alpha * z is -3.0e9 to -3.9e9 at the rule's far nodes, whose
weight exp(-5e18) is 0, and 0 times Inf is NaN. Section 2's claim is
withdrawn there. Changes to the record:

- `dev/test-backlog.md`: the LOW "where the NaN arises" is replaced by
  a MEDIUM owned by the skew normal's density, with the review's soft
  floor `-1e6 + logspace_add(0, w + 1e6)` as the candidate, NOT built
  (it moves every skew-normal density at the ulp level).
- `dev/upstream-bugs.md`: the fact is RTMB-3 ("Derivatives of `pnorm(x,
  log.p = TRUE)` lose accuracy in the lower tail", status not reported,
  draft ready), which lane optima had already recorded with NaN near
  -4.46e9. I added the review's numbers and the skew-normal route to it
  as an "Also" paragraph rather than a second entry for the same
  defect.

### The veto counts only a restart that left the point (finding 2)

`escape_stationary()`: a restart that ended on its own still replaces
the kept optimum when its objective is lower (as in 0.83.0), but it
vetoes an abandoned point only when it left the point's window: the
number of predictors on their point (`attr(stationary_escapes(obj, o,
frame), "n_on")`, the escape's own window, `tol` of the declaration)
fell below the number at the start of the pass.

Test (`test-quad-escape-record.R`, "a restart back on the point does
not veto an abandoned one"), the review's construction made
deterministic: a skew normal `y ~ 1`, n = 200, seed 4041, true alpha 4,
started at alpha = 0. `optimize_obj()` is mocked: the restart from
alpha = +2 runs to the maximum and then raises; the one from -2 returns
the stall with an objective 1e-12 lower. Seen failing
(`t-fr-prefix-...` on the second build, `t-fr-base-...` on 0.83.0, ref
BLAS): the fit stays at alpha 6e-15, code 0, "relative convergence
(4)", only a standard-error warning; second build pass=23 fail=4,
0.83.0 pass=16 fail=6 err=1. Fix-round build: pass=27 fail=0 err=0
skip=0 warn=0. Over 21 roundings (y times 1 + 2e-15 k, k = 0..20,
`dev/quadesc-veto-ulp.R`, logs `veto-lib-{ref,ob332,ob326}.txt`): PASS
21 of 21 under each of ref, 0.3.32 and 0.3.26 (alpha 3.44123, gain
11.309, evals 43 on every one).

### Counts on an abandoned point (finding 5)

The abandoned result now carries `evals` (the tally's count) and, for
nlminb, `nonfinite_trials`. When the caller passes no tally (the
Laplace escape inside `quad_fit()`), each restart gets a fresh one,
which `optimize_obj()` would create itself, so a finished restart is
unchanged. Seen failing: `dev/quadesc-fix3-seen.R` assigns the lane's
`escape_stationary()` without the two count lines into the installed
namespace and runs the test file (log `fix3-seen.txt`): the mock test
fails 2 of 8 (`evals` NULL, `nonfinite_trials` NULL); with them, 0.

### Test, optional items

- The fit test now also runs k = 0, the data the 26.04 runner failed
  on, beside k = 9. 0.83.0 passes k = 0 locally (it escapes there under
  every local BLAS) and fails k = 9.
- NEWS: the bullet moved to "Convergence and standard errors" and says
  that the change acts wherever the escape runs (Laplace fits, and both
  escapes of a quadrature fit).
- Review minor (section 2 of the review), recorded and not fixed: the
  kept abandoned point is the one the frozen rule ranks higher, which
  is not always the higher in exact terms. ALPHAX k = 9: restart 1 rule
  objective 84.560080, exact -85.105859; restart 2 84.578059, exact
  -85.024758 (0.081 higher). Valid endings span -85.106 to -84.974
  exact (0.13). No simple sound fix: the rule is the only objective the
  fit has, and re-anchoring at a run-off alpha does not tape. Added to
  the MEDIUM on the rule's accuracy in the backlog.

### Fix-round runs

Sweep of the fix-round build (`dev/quadesc-sweep-d.sh`, logs
`sw-d-lib-*`, summary `sweep-sum-d.txt`): 189 fits, 0 on alpha's point,
smallest |alpha| 5917.7, and logLik identical to 10 decimals to the
second build's sweep (section 5) on 189 of 189. So every ALPHAX ending
of k = 0..20 under the three builds meets the fit test's assertions as
before (exact at least -85.1020 against the Laplace fit's -86.009814).

Escape-related files, lane library, each fail 0 err 0 skip 0 warn 0
under ref, 0.3.32 and 0.3.26 with LAPACK (`t-fr1-<ver>-*.txt`):
test-quad-escape-record.R 27, test-stationary-passes.R 4,
test-quadrature-defects.R 89, test-skew-normal-start.R 70,
test-optima.R 107, test-quadrature-nested.R 69, frmtmb.latent
test-hmm-stationary.R 12. `dev/ciharden-out/` deleted afterwards.

One session per suite, three orders, fix-round build
(`dev/release/run-sessions.sh quadesc-sess2 2 frmtmb frmtmb.latent
frmtmb.sample`, every session log names the lane library):
`QUADESC-SESS2 ran 9 of 9`, fail 0 err 0 warn 0 in each; frmtmb
forward pass 28775 skip 15, reverse 28940 skip 15, preload 29044 skip
17 (13 more than section 8: the test file's new expectations);
frmtmb.latent 700 each, skip 2; frmtmb.sample 4651 each, skip 1. LEAK
lines without a load tag or the extension-state tag: 0.

R CMD check `--as-cran` of the fix-round build (`dev/quadesc-check.sh`):
frmtmb `Status: 1 NOTE`, the V8 math-rendering NOTE of the HTML manual;
tests 41 min, OK. `dev/quadesc-check/` deleted afterwards.
