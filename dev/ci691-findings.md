# Lane ci691: the Ubuntu and macOS failures of 0.69.0

Base: 48c1e6c9 (frmtmb 0.69.0), built in
`C:/Users/adf44/source/r/rellib-r7`. Worktree `frmtmb-wt-ci691`, private
library `wt-ci691-lib`. Logs are in `dev/ci691-log/` and
`dev/ci691-suite-*/` (gitignored); every number below is copied from
them, with the script that made it. "ref" is R's reference BLAS and
LAPACK; "ob26" and "ob32" are the copies of R that
`dev/ciharden-openblas.sh 0.3.26 lapack` and `0.3.32 lapack` build,
run with `OPENBLAS_NUM_THREADS=4`.

## What failed on CI

Runs on 48c1e6c9, read with `gh run view --job <id> --log`
(`dev/ci691-log/job-<id>.log`):

| Workflow, job | Result | Failures |
|---|---|---|
| R-CMD-check ubuntu-24.04 release, devel | FAIL 3, WARN 2 | xbeta :296, :311; optima :163 |
| R-CMD-check macos-latest | FAIL 2, WARN 1 | xbeta :296, :311 |
| R-CMD-check windows-latest | success | |
| test-coverage (ubuntu) | FAIL 3, WARN 2 | as ubuntu-24.04 |
| ubuntu-next ubuntu-24.04 (2026-10-12) | FAIL 1, WARN 1 | optima :163 |
| ubuntu-next ubuntu-26.04 (2026-10-12) | FAIL 2, WARN 1 | optima :163; se-check :546 |

The two escaped warnings (the task's item 3) are
`test-xbeta-zibb-hurdle-cum.R:290` ("Standard errors are not available
for 3 of 3 parameters") and `test-nonfinite-gradient.R:58` ("Optimizer
did not report convergence: singular convergence (7)"). macOS has only
the first. ubuntu-26.04's `test-se-check.R:546` error was not in the
task; it is fixed here too.

The ubuntu-next ubuntu-24.04 job ran the same commit with the same R
(4.6.1), RTMB (2.0), TMB (1.9.25), Matrix (1.7-5) and OpenBLAS
(0.3.26) as the R-CMD-check ubuntu-24.04 job, and the xbeta test
PASSED there. So the platform does not decide that failure; rounding
does, and it flips between runners of one image.

## Why the round's OpenBLAS emulation did not catch them

None of the three is a BLAS effect. Each test asserted the outcome of
a branch that a few ulps decide: where nlminb stops on a ridge, and on
which side of the SE check's threshold the Hessian then falls; whether
PORT ends on a trial it rejected; whether PORT calls a singular
Hessian "singular convergence". These fits barely touch BLAS, and the
emulator keeps this machine's compiler, libm and CPU, so it gives one
more draw of the same coin, not the runner's draw. macOS, which has no
OpenBLAS, failing the xbeta test is the direct evidence. The emulator
passes all four files on rellib-r7 (`before-ob-*.log`, `t-base-ob26-*`
and `t-base-ob32-*`: xbeta 606/0, nonfinite-gradient 5/0, se-check
126/0, the old optima tests all pass).

The instrument that does reproduce them is a perturbation: move the
start (or the data) by a few ulps and refit. Item 21 of
`dev/rtmb-pitfalls.md` (P and E cores) is the same coin on this
machine: two processes with identical inputs gave different phi 2e4
xbeta fits (`xbeta-perturb-ref.log` against `xbeta-cert-ref-2e4.log`,
k = 0: 308 iterations to logLik 4049.107337 in one, 4049.107365 in the
other).

## Failure 1: xbeta, `test-xbeta-zibb-hurdle-cum.R:296` and `:311`

CI: `all(is.finite(sqrt(diag(vcov(fit)))))` FALSE, then
`max(abs(g) * se)` NaN at phi 2e4, after the SE check warned that all
three parameters are flat.

The likelihood (`dev/ci691-xbeta-profile.R`, `xbeta-profile.log`,
the reference written from `stats::dbeta()` and `pbeta()` alone; seed
4200, n = 1000). Profiled over log kappa, at phi 2e4 the maximum is
4049.10738442 at log kappa -1.5 (logit mu -1.692, log phi 11.53) and
the profile falls to 4049.10142593 at log kappa 3, toward the normal
limit kappa = Inf: 0.0060 in all. At phi 1e4 the maximum is 3371.32315712
at log kappa -3.5 and the fall to log kappa 3 is 0.52. So at phi 2e4
mu, phi and kappa lie on a curved ridge that the data barely bend.

What the fit does on Windows (`dev/ci691-xbeta-perturb.R`, 25 fits per
arm, start k = 0 the default and k = 1..24 moved by 1 to 4 ulps):

| Arm, phi | fits | code 0 | SEs lost | lost with code 0 |
|---|---|---|---|---|
| ref, 2e4 (`xbeta-perturb-ref.log`) | 25 | 8 | 8 | 8 |
| ob26, 2e4 (`xbeta-perturb-ob-2e4.log`) | 25 | 11 | 9 | 9 |
| ref, 1e4 (`xbeta-perturb-ref-1e4.log`) | 25 | 4 | 0 | 0 |

That is the CI failure, reproduced: a start a few ulps away loses all
three standard errors in about a third of fits at phi 2e4, and never at
1e4. Every lost case had code 0. A fit with code 1 (false convergence)
is explained by the convergence warning, so the SE check does not run
and its standard errors are shown unchecked (`se_explained` is
"convergence").

Why the check removes them (`dev/ci691-xbeta-tier3.R`, k = 2,
`xbeta-tier3-k2.log`): at the code-0 estimate (-1.705, 11.51, -1.516)
the unit-diagonal AD Hessian has eigenvalues 2.013, 0.9868 and
1.076e-9, a ratio of 5.34e-10, under `se_flat_tol` (1e-9). Then
`se_curvature_real()` steps straight along that eigenvector: the second
difference is 4851 at the full step and 220.1 at the half step, where
the quadratic model predicts 0.004, a ratio of 22.0 outside (3, 5.5).
That is the rule for a straight step off a curved ridge, so the
direction is called flat, and it loads on all three parameters. The
plain inverse of the same Hessian, which is positive definite, gives
standard errors 6.07, 8.75 and 7.48 on the link scale. Across 25 fits
the ratio ranges 1.758e-10 to 3.198e-09 at phi 2e4, astride 1e-9, and
is 1.362e-06 at phi 1e4 in all 25 fits of both arms
(`dev/ci691-xbeta-ratio.R`, `xbeta-ratio-*.log`).

Verdict: the fit is not wrong. The estimates match the reference's
log-likelihood (the test's `:295` passed on every runner), and the
code-0 fits sit within 3.4e-5 of the profile maximum (lowest
4049.107351, ob26). Whether the standard errors are shown is decided
by nlminb's code and by the side
of 1e-9 that rounding puts the eigenvalue ratio on. With a fall of
0.006 from the maximum to kappa = Inf, "the data do not determine
them" is a fair description, and Wald standard errors of 6 to 9 on the
link scale say the same. So the test asserted something not guaranteed.

The test's second assertion was fragile as well. `max(abs(g) * se)`,
with `g` the reference's central-difference gradient
(`dev/ci691-xbeta-cert.R`, `xbeta-cert-*.log`), reached 0.246 (ob26)
and 0.196 (ref) at phi 2e4 on fits with standard errors, against a
bound of 0.01, at points whose logLik is at most 4.2e-4 below the
profile maximum (lowest fit 4049.106969, ob26). At phi 1e4 it reached
0.00902, 0.9 of the bound. The product of a gradient with each
coordinate's standard error is not invariant on a ridge; the Newton
gain `g' H^-1 g / 2`, with `H` the fit's AD Hessian, is the
log-likelihood the quadratic model says a step would still gain. Its
largest value was 2.91e-4 (ob26) and 2.65e-4 (ref) at phi 2e4, and
1.2e-9 at phi 1e4.

Fix, test only (`tests/testthat/test-xbeta-zibb-hurdle-cum.R`):

- the SE warning is allowed, not required;
- at phi 1e4 every standard error must still be finite (ratio 1.4e-6,
  far from the threshold);
- at phi 2e4 a standard error is either shown or named by the check
  (`sdr_of(fit)$se_lost`), never NaN without a reason;
- the optimality check is the Newton gain, below `fit$control$grad_tol`.

Seen to hold on both outcomes (`dev/ci691-xbeta-testbody.R` on the
fix, `xbeta-testbody-*.log`): ref 2e4, 25 of 25 fits pass every
assertion, 8 with standard errors lost, largest gain 0.000265; ob26
2e4, 25 of 25, 9 lost, 0.000291; ref 1e4, 25 of 25, 0 lost, 1.2e-09.

## Failure 2: `test-optima.R:163`, a code defect

CI: `fit$obj$fn(fit$opt$par)` was 374.16747089090501 and
`fit$opt$objective` 374.16747089090489 (seed 12), on every Ubuntu job.
frmtmb did not error there: the fit returned with the expected
warnings, and only the identity failed. The test name ("warns, not
errors") describes 0.68.1's defect, which stays fixed.

Cause. PORT can stop on a trial point it rejected and return that
trial as `par` with the best value as `objective`. 0.69.0's
`nlminb_best_par()` (lane optima) replaced `par` with the best point
only when that trial's objective was not finite. PORT also does it on
a FINITE rejected trial. `dev/ci691-port-probe.R` runs nlminb on the
cs() mixture's dry-run objective, seeds 1 to 20, each from the default
start and 4 starts moved by a few ulps (`port-probe-ref.log`,
`port-probe-ob.log`):

| Arm | runs | fn(par) is objective | last trial not finite | finite, rejected |
|---|---|---|---|---|
| ref | 100 | 57 | 33 | 10 |
| ob26 | 100 | 55 | 35 | 10 |

Every finite rejected trial had code 1 and `fn(par) - objective` of
5.7e-14 to 5.1e-13 (1 to 9 ulps), the CI shape. In 100 of 100 runs per
arm the best value tracked by `nlminb_trial_fn()` was the reported
objective exactly, and `fn()` at the best point reproduced it. Seed
12 was a finite rejected trial from one moved start in each arm.

Decision: a code defect. `fit$opt$par` and `fit$opt$objective` are
documented (and asserted) to describe one point, and every caller
reads the estimates from `par`.

Fix, `R/fit.R`, `nlminb_best_par()`: the best finite point replaces
`par` whenever `par` is the last trial and that trial's value is not
the reported objective. A run that ended on the point it reports is
returned as it was, so a fit that PORT did not stop on a rejected
trial is unchanged bit for bit.

New test, `test-optima.R`, "a finite trial PORT rejected is not
reported either": run_optimizer() on seeds 1, 9 and 13 from the
default start, and the rule on a made-up run (a finite last trial
above the objective is replaced; a last trial at the objective is left
alone, the guard's absent case). On rellib-r7 it fails in every arm
(`t-base-*-test-optima.log`): ref pass=102 fail=5 (seeds 1, 9, 13, for
example seed 1 `fn(par)` 394.679712369822425 against 394.679712369822369,
and the two rule assertions), ob26 and ob32 fail=3 (seed 13 and the
rule). `dev/ci691-port-seeds.R` in 8 fresh processes per build:
rellib-r7 gave DIFF on all three seeds in 8 of 8, the fix "same" on
all three in 8 of 8 (`portseeds-*.log`), so the pin does not depend
on the core type.

## Failure 3: `test-nonfinite-gradient.R:58`, an escaped warning

The degenerate three-component ordinal mixture can end with "singular
convergence (7)", which the test did not allow. That is PORT's honest
report on a singular Hessian, and the degenerate-component warning
(required) is given beside it. `dev/ci691-mix3-probe.R` refits the
test's data with each x moved by a few ulps (`mix3-probe-*.log`):
ref 3 of 17 fits singular convergence, ob26 2 of 17; 17 of 17 gave
the degenerate-component warning and a finite gradient in each arm.
Fix, test only: "did not report convergence" is allowed, not
required.

## Failure 4: `test-se-check.R:546` on ubuntu-26.04

`startsWith(lost, "theta_")` stopped with "non-character object(s)":
`sdr_of(fit)$se_lost` was NULL, so the check removed nothing on that
runner. The property asked is that whatever is removed is a theta,
which an empty set satisfies. Fix, test only: `as.character()` around
the names. Not reproduced here: ob32 passes the file on rellib-r7
(`t-base-ob32-test-se-check.log`, 126/0), so whether the lme4
comparison that follows holds on 26.04 is not verified.

## Verification

Installed core 0.69.1 from the worktree into wt-ci691-lib (roxygenised;
no Rd changed). The four files, one per process, base and fix
(`t-<build>-<arm>-<file>.log`; pass/fail/error/skip/warn):

| File | base ref | base ob26 | base ob32 | fix ref | fix ob26 | fix ob32 |
|---|---|---|---|---|---|---|
| test-optima.R | 102/5/0/0/0 | 104/3/0/0/0 | 104/3/0/0/0 | 107/0/0/0/0 | 107/0/0/0/0 | 107/0/0/0/0 |
| test-xbeta-zibb-hurdle-cum.R | 606/0/0/0/0 | 606/0/0/0/0 | 606/0/0/0/0 | 606/0/0/0/0 | 606/0/0/0/0 | 606/0/0/0/0 |
| test-nonfinite-gradient.R | 5/0/0/0/0 | 5/0/0/0/0 | 5/0/0/0/0 | 5/0/0/0/0 | 5/0/0/0/0 | 5/0/0/0/0 |
| test-se-check.R | 126/0/0/0/0 | 126/0/0/0/0 | 126/0/0/0/0 | 126/0/0/0/0 | 126/0/0/0/0 | 126/0/0/0/0 |

(The base columns run the revised tests against 0.69.0; the base
failures are the new pin.)

All eight suites on the fix, reference BLAS, one file per process,
NOT_CRAN=true, ungated (`dev/ci691-par.sh fixref`,
`dev/ci691-suite-fixref/ci691-fixref.sum`: "SUITE fixref ran 364 of
364"; all 364 `lib:` lines read `wt-ci691-lib/frmtmb 0.69.1`).
`nlminb_best_par()` is on every nlminb fit, so every extension suite
can reach it:

| Package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 223 | 15179 | 0 | 0 | 177 | 0 |
| frmtmb.coupling | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 29 | 1746 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 10 | 360 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 430 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 50 | 2474 | 0 | 0 | 4 | 0 |
| frmtmb.spline | 15 | 590 | 0 | 0 | 1 | 0 |

The extensions load from rellib-r7 behind the private library; none
changed, and none needs a new floor. After this run only two comments
in `test-optima.R` changed.

`R CMD check --as-cran` of frmtmb 0.69.1 on the final tree, built and
checked in `dev/ci691-check/frmtmb/` (`dev/ci691-check.sh`): Status:
1 NOTE, the V8 math-rendering NOTE on the HTML manual. Tests: FAIL 0,
WARN 0, SKIP 198, PASS 15117 (CI on 0.69.0: FAIL 3, WARN 2, SKIP 198,
PASS 15108 on ubuntu-24.04).

## Version and files

`R/fit.R` changed, so frmtmb goes to 0.69.1 (DESCRIPTION,
codemeta.json) with a `# frmtmb 0.69.1` section in NEWS.md, after the
0.68.1 precedent. No extension needs a new floor. Files touched:
`R/fit.R` (`nlminb_best_par()` and its comment), `DESCRIPTION`,
`codemeta.json`, `NEWS.md`, `tests/testthat/test-optima.R`,
`test-xbeta-zibb-hurdle-cum.R`, `test-nonfinite-gradient.R`,
`test-se-check.R`, and the lane's scripts `dev/ci691-*.R`,
`dev/ci691-*.sh` and this file. No Rd file changed.

## Found, not fixed

- The SE check's curved-ridge rule cannot tell a ridge the likelihood
  does not move along from one it bends along slightly: here the
  profile falls 0.006 to kappa = Inf and the curvature is real, and
  the check calls the three parameters "flat" where the Hessian is
  positive definite. Defensible on this fit (see the verdict above),
  but the wording "flat along them at the estimates" is stronger than
  what was measured. It is the SE check's design, overlapping lane
  seflags, so it is for the round's owner.
- The SE check does not run on a fit with a convergence warning, so on
  this ridge a fit nlminb calls false convergence shows unchecked
  standard errors (0.33, 0.46 and 0.40 at the Windows default-start
  point, `xbeta-probe-ref.log`) while a converged fit a few ulps away
  shows none. The reported uncertainty flips with rounding.
- `test-se-check.R` "an identified remainder keeps its SEs" on
  ubuntu-26.04: why nothing was removed there is not known.

## What this could not verify

The runners themselves. Every mechanism above is reproduced on Windows
by construction (ulp-moved starts or data), and the CI values match
its shape (the 2-ulp objective gap, the three "flat" parameters at
code 0, singular convergence), but the runners' own draws can only be
seen on a push.
