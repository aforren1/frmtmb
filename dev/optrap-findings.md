# Lane optrap: fits that stop short of the maximum, or end on luck

Round of 2026-10-07, base `48c1e6c9` (frmtmb 0.69.0), worktree
`frmtmb-wt-optrap`, private library `wt-optrap-lib`, before arm
`rellib-r7`. Nothing is committed.

STATUS: the consolidation fixes done (first section below); punch round
2's and punch round 1's sections follow as written then. Round 0's
sections follow as written then, with corrections marked; their generated blocks
are from round 0's runs, of which the sweep, the `mo()` studies and
the scale study ran before round 0's last source edit (each
reproduced on round 0's final build in the review). The ledger row
`brmsfit-methods:326` flips back (see "Left open").

## Consolidation fixes

Two defects the 0.70.0 consolidation found on the merged build
(`frmtmb-wt-release`, `rellib-r8`; its `dev/round-20261013.md` and
`dev/test-backlog.md`, "Found at the consolidation"), and the medium
item filed with them. Final lane build: install 32
(`dev/optrap-log-install32.txt`). Scratch merge: the release tree
copied to `dev/optrap-relmerge/` (without its `dev/` logs) with this
round's diffs against build 27 applied to its `R/fit.R` and
`R/autoscale.R` (`dev/optrap-c3-apply-merge.sh`; the diffs are
`dev/optrap-log/c3/fit.patch` and `autoscale.patch`, the base
`dev/optrap-c3-base27-fit.R` and `dev/optrap-log/c3/base27-autoscale.R`)
and this tree's `tests/testthat/test-optrap.R`, installed into
`optrap-relmerge-lib` ahead of `rellib-r8`
(`dev/optrap-c3-install-merge.sh`). Logs: `dev/optrap-log/c3/`.

### Changed functions (for the re-merge)

- `R/fit.R`, `newton_polish()`: new argument `ill`; with random effects
  it runs when `ill` is true, on a Hessian from `stats::optimHess()`
  in natural units (`parscale = par_units`); a step is kept only if it
  gains at least 0.1 of its prediction; the steps are undone when they
  end where the Hessian has no Cholesky factor.
- `R/fit.R`, `fit_assembled()`: computes `ill` (random effects, at most
  100 outer parameters, code 0, `predictor_design_kappa(frame) > 1e5`)
  and passes it; the verbose line reads "Newton steps from the stop".
- `R/autoscale.R`, `predictor_design_kappa()` (new): the largest
  condition number of a frame's fixed designs, columns scaled to unit
  length, aliased columns left out, on evenly spaced rows up to
  `frmtmb.unit_ls_cells`.
- `tests/testthat/test-optrap.R`: five new tests (below).
- `NEWS.md`: the Newton-step bullet (random effects, the tenth, the
  undo). The release tree's NEWS is the integration's to merge.

### 1. A step that gains far less than it predicts

The integration's candidate gate was a quarter. Measured on the merge
(`dev/optrap-c3-gate.R`, which patches the gate in-process): it holds
the negative binomial data, but it refuses the first step of a degree-7
raw polynomial (seed 17, no group term), whose actual gain is 0.203 of
the prediction because log sigma is not a quadratic direction; the
next steps gain 1.036, 1.004 and 1.058 of theirs and reach lm()'s
maximum (`dev/optrap-c3-trace17.R`, `trace17-merge.txt`). With a
quarter that fit stays 1.16 short in silence on the merge. A tenth
keeps it and still refuses the negative binomial step (3e-8 of a
prediction above 1e-6). Two fits it keeps on the negative binomial
data (seeds 25 and 28, gain 2e-6) ended where the Hessian is slightly
indefinite (-2.8e-6), which turned 0.69.1's silent fit into an SE
warning on the shape; undoing steps that end without a Cholesky factor
removes that too (and is item 3's fix).

Negative binomial at the Poisson limit, 200 seeds
(`dev/optrap-c3-negbin.R`, the consolidation's data in its order):

| build | warned | coefficient SE lost | polished | largest gain |
|---|---|---|---|---|
| 0.69.1 (the consolidation's log) | 66 | 0 | 0 | - |
| rellib-r8 | 74 | 3 | 11 | 2.1e-6 |
| merge, gate 0.25 (build 31) | 67 | 0 | 1 | 1.7e-6 |
| merge, gate 0.1, no undo | 68 | 0 | 4 | 2.1e-6 |
| merge, final (0.1 and the undo) | 66 | 0 | 2 | 9.0e-7 |
| lane build 32 | 66 | 0 | 2 | 9.0e-7 |

New test "a Newton step that gains far less than it predicts is
refused" (seed 2): fails with punch round 2's `newton_polish()`
(`OPTRAP_ROUND2=1`, `t-optrap-round2.txt`), passes on build 32. New
test "a Newton step the quadratic model overrates is still taken"
(degree 7, seed 17): fails with the gate at 0.25 (`OPTRAP_GATE=0.25`,
`t-optrap-gate025.txt`), passes on build 32.

### 2. Raw polynomials with (1 | g): the cause and the fix

The cause is the coefficient units. Ablations on build 27
(`dev/optrap-c3-repoly.R` with `OPTRAP_ABLATE`, degrees 6 and 7, 20
seeds, against the exact profile maximum):

| build 27 with | deg 6 short (silent) | deg 7 short |
|---|---|---|
| everything | 13 (13), to 1.1 | 20, to 3.1 |
| no theta start | 14 (14), to 1.1 | 20, to 3.1 |
| no inner tolerance | 13 (13), to 1.1 | 20, to 3.1 |
| no coefficient units | 5 (5), to 0.023 | 15, to 1.78 |
| no units at all | 5 (4), to 0.023 | 16, to 1.78 |

So without the coefficient units these fits are 0.69.1's, which also
stops short (5 of 20 at degree 6, 4 silent). The units make every
response behave as if its residual spread were 1, and on this design
nlminb's relative function test stops short more often there, as on
the fixed-effect polynomials of punch round 1. Restricting the units
would bring back 0.69.1's shortfalls and its scale dependence, so the
fix is the Newton steps, extended to random effects where a design is
near collinear. The Hessian of the Laplace objective comes from central
differences of its exact gradient in natural units; the measured
prototype (`dev/optrap-c3-repolish.R`, `repolish-b27.txt`) took 13 and
20 short fits to the maximum at degrees 6 and 7, at most 0.03 s each.
The steps run only above a condition number of 1e5
(`predictor_design_kappa()`: 3e5 at degree 5, 4e6 at degree 6, 3.4 for
sleepstudy), so other fits with random effects pay one QR of their
design (at most 5e6 cells, 0.3 s at the cap) and are otherwise
unchanged.

The study (`dev/optrap-c3-repoly.R`, 20 seeds per degree, short = 1e-4
below the exact maximum; `repoly-*.txt`):

| | deg 5 | deg 6 | deg 7 | deg 8 |
|---|---|---|---|---|
| (1 \| g), 0.69.1 | 0 | 5 (4 silent), to 0.023 | 16, warned | 20, warned |
| (1 \| g), rellib-r8 | 0 | 13 (13 silent), to 1.1 | 20 (20 silent), to 3.1 | 20, warned |
| (1 \| g), merge | 0 | 0 | 0 | 20, warned |
| (1 \| g), lane build 32 | 0 | 0 | 0 | 20, warned |
| no group term, 0.69.1 | 2 (2 silent), to 4.7e-4 | 11 (10 silent), to 0.33 | 20, warned | 20, warned |
| no group term, rellib-r8 | 0 | 0 | 1 (silent), 4.0e-4 | 20, warned |
| no group term, merge | 0 | 0 | 1 (silent), 4.0e-4 | 20, warned |
| no group term, lane build 32 | 0 | 0 | 1 (warned), 4.0e-4 | 20, warned |

(lmer is short on 16 of 20 at degree 7 and 7 at degree 8, and fails on
some degree-8 seeds;
"warned" at degree 8 is the SE check's ill-conditioning or flat
verdict, so every degree-8 fit says it is unreliable.) The one
remaining silent shortfall at degree 7 without a group term (seed 4)
is 4e-4 below the maximum, as on rellib-r8.

New test "a near-collinear fit with a group term reaches the maximum"
(degree 6, seeds 1, 13, 14; within what moving log sigma by 1e-2 of
its SE costs, against lmer): fails 3 of 3 with punch round 2's
`newton_polish()`, passes on build 32.

### 3. The B5 design without random effects (medium; fixed)

`dev/optrap-c3-b5.R` (degree 8 plus z and w, seeds 4 and 5):

| build | steps | SE ratio to poly(x, 8): z, wb, sigma | warning |
|---|---|---|---|
| rellib-r8 | 2, 1 | NaN, NaN, NaN | large gradient, curvature unusable |
| merge | 0, 0 | 1.0018, 1.0034, 1.0057 and 1.0034, 1.0028, 1.0076 | the SE check's, for the 9 polynomial coefficients |
| lane build 32 | 0, 0 | 1.0011, 1.0029, 1.0000 and 1.0024, 1.0024, 1.0000 | the same |

The steps gained what they predicted and ended where the Hessian does
not factor; they are now undone, so the fit is lane seflags' again. The
end point's gradient is not the test: at the maximum of a degree-6
polynomial it is 1.6e-3 to 2.9e-3 in natural units, above `grad_tol`,
with no headroom left. New test "Newton steps ending where the Hessian
does not factor are undone": fails 4 assertions with punch round 2's
`newton_polish()`, passes on build 32.

### What did not move (build 32 and the merge)

- The high-R^2 design (`dev/optrap-rev-r2b.R`, `r2b.txt`): 0 / 0 / 0
  at k = 300, 1e3, 3e3, 1e4.
- The wide design (`dev/optrap-rev2-b1.R wide`,
  `dev/optrap-log/rev2/b1-wide-c3.tsv`): unit 1; seeds 2 and 3 at
  lmer's maximum; seed 1 9107.66 short at code 0, the Open-high item of
  punch round 2 (the design's condition number is small, so the steps do
  not run).
- The scale study (`lanescale.txt`): 0 short in every row, the same
  invariance as punch round 2.
- The `mo()` study (`dev/optrap-c3-mo-sum.R`, `mo-sum.txt`): 200 / 200
  and 200 / 200 at the exact maximum, 115647 and 102810 evaluations,
  as punch round 2.
- The polynomial study without a group term (`poly-b32.txt`,
  `poly-merge.txt`): 0 of 20 short at degrees 3 to 6, as before.

### Reruns

`dev/optrap-scan.sh` on the scratch merge (`TREE=dev/optrap-relmerge`,
`ROUND_BASE=rellib-r8`, so the merge's core from `optrap-relmerge-lib`
and the release's extensions from `rellib-r8`), one test file per R
process, every gate on: `c3m-ref` under the reference BLAS, `c3m-ob26`
under OpenBLAS 0.3.26 as BLAS and LAPACK, 4 threads. Generated by
`OPTRAP_COUNT_LIB=optrap-relmerge-lib dev/optrap-counts.R c3m-ref c3m-ob26`
(`dev/optrap-log/c3/counts.md`):

<!-- BEGIN GENERATED: dev/optrap-counts.R c3m-ref c3m-ob26 -->

### c3m-ref

371 jobs; 371 files with a RESULT line; 371 whose `lib:` line names the lane's library.

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 229 | 18946 | 0 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1765 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 51 | 2806 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **371** | **26089** | **0** | **0** | **0** | **0** |

Files with a failure, an error, an escaped warning or no result: 0


### c3m-ob26

371 jobs; 371 files with a RESULT line; 371 whose `lib:` line names the lane's library.

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 229 | 18947 | 0 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1765 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 51 | 2806 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **371** | **26090** | **0** | **0** | **0** | **0** |

Files with a failure, an error, an escaped warning or no result: 0

<!-- END GENERATED -->

The release's own tests pass with these changes under both BLAS. The
ledger row `brmsfit-methods:326`, which fails on this lane's tree
because its ledger still records it as pass, passes on the merge with
the release tree's ledger.
The five new tests are in the merge's `test-optrap.R`.

R CMD check of this lane's tree (`dev/optrap-check.sh`, after the last edit, `dev/optrap-check/check.log`): Status: 1 NOTE (`checking HTML version of manual`); tests `[ FAIL 0 | WARN 0 | SKIP 383 | PASS 11935 ]`.

## Punch round 2

The review's "Re-check after punch round 1" (`dev/reviews/2026-10-08-optrap.md`,
scripts `dev/optrap-rev2-*`). Final build: install 27
(`dev/optrap-log-install27.txt`). Logs are in `dev/optrap-log/p2/`, and
the reviewer's scripts ran through `dev/optrap-p2-redirect.R` so that
their tables land in `dev/optrap-log/rev/` and `dev/optrap-log/rev2/`
(runner `dev/optrap-p2-run.sh`, job list `dev/optrap-p2-jobs.txt`).
Round 1's outputs that had gone into the reviewer's
`dev/optrap-rev-log/` (240 files, `mo/` 64, and this lane's 161 lines
of `run-ends.txt`) are moved to `dev/optrap-log/rev-p1/`.

### What changed in this round

- `R/autoscale.R`, `predictor_resid_sd()`: the spread is the MAD of
  the least-squares residuals (the sd where the MAD is 0); a design of
  more than `cells` entries (`n * p`, default 5e6, internal option
  `frmtmb.unit_ls_cells`) is fitted on evenly spaced rows up to that
  size instead of being skipped; a design whose rows do not match the
  response gives no unit. The raw spread of y is never used.
- `R/fit.R`, `newton_polish()`: a step needs a Cholesky factor of the
  natural-unit Hessian; the gate also reads optim's `reltol`, `abstol`
  and L-BFGS-B's `factr`, `pgtol`, and skips a caller's optimizer
  function.
- Tests: `test-optrap.R` four new tests (outliers and the MAD, the size
  cap, the indefinite Hessian, optim's tolerances) and the unit test
  reads the MAD; `test-se-tier.R`'s healthy sleepstudy fit is in units
  of 32 ms (the MAD gives 64 ms a unit of 0.5).
- `NEWS.md`: the breaking bullet says MAD and evenly spaced rows; the
  Newton-step bullet says positive definite and names the gate.
- `dev/optrap-testfile.R`: `OPTRAP_ROUND1=1` patches in round 1's
  `predictor_resid_sd()` (sd, no least squares above the cap) and
  `newton_polish()` (`dev/optrap-p2-round1.R`).

### B1 above the cap (fixed)

`dev/optrap-rev2-lscost.R` (`p2/lscost.txt`): no unit (residual spread
1) on both sides of 5e6 cells, n 99000 and 101000 by 50, 990000 and
1.01e6 by 5, in 0.28 to 0.67 s. Far above it (`dev/optrap-p2-lscost.R`,
`p2/lscost-big.txt`): 1e6 by 50 dense, no unit, 0.54 s; 1e5 rows with
a 500-level factor (a 1e5 by 501 design), no unit, 2.28 s (1.03 s on
an idle machine), because the rows kept are 5e6 / p.

The review's wide fit (`dev/optrap-rev2-b1.R wide`, n 101000, 50
columns, seeds 1 to 3; `rev2/b1-wide-p2.tsv`): unit 1, code 0, 0
warnings, logLik - lmer -9107.66, -1.5e-9, -3.5e-8, which is 0.69.0's
fit and the reviewer's no-cap arm row for row (round 1: unit 8192, 3 of
3 at the evaluation limit 8488 to 8964 short). The narrow version
(101000 by 5): unit 1, 0 short. The reviewer's `nocap` and `mad`
wrappers patch lines that no longer exist (`<= 5e+06`; the sd), so the
cases they ran were rerun on the build itself, which now does both.

New test "a design above the least-squares size cap keeps the residual
unit" (cap lowered to 150 cells so 100 by 2 crosses it; the capped unit
within one power of two of the full design's): 3 of 3 fail with round
1's rule (`OPTRAP_ROUND1=1`, `p2/t-optrap-round1.txt`; unit 8192), pass
on build 27 (`p2/t-optrap-b27.txt`).

### r1: the MAD

Heavy tails (`dev/optrap-p2-tcmp.R`, `p2/tcmp.txt`; short = 1e-4
below the best any of build 27, build 26 or rellib-r7 reached at either
scale):

| data | scale | build 27 | build 26 | rellib-r7 |
|---|---|---|---|---|
| t1 | 1 | 0 | 1 (-29.8) | 0 |
| t1 | 1e3 | 0 | 1 (-39.6) | 8 (to -0.59) |
| t2 | 1 | 0 | 0 | 0 |
| t2 | 1e3 | 0 | 0 | 9 (to -2.83, 1 silent) |

New test "a few outliers do not set the unit" (seed 18 of t1: the
unit within a factor sqrt(2) of the residuals' MAD, and their sd more
than 8 times it): fails with round 1's sd, passes on build 27.

Unchanged by the MAD: the high-R^2 design (`dev/optrap-rev-r2b.R`,
`p2/r2b.txt`) 0 / 0 / 0 at k = 300, 1e3, 3e3, 1e4 and
`dev/optrap-rev-r2.R` 0 short; the scale study (`dev/optrap-scale.R`,
`p2/lanescale.txt`) 0 short in every row, invariance 1.68e-7 at most
(ri20) and 1.94e-8 (slope); the `mo()` study (`dev/optrap-p2-mo-sum.R`,
`p2/mo-sum.txt`) 200 / 200 and 200 / 200 at the exact maximum with
115647 and 102810 evaluations (the reviewer's MAD arm to the
evaluation); the unit-1 boundary grids (`bnd_lo`, `bnd_hi`) 0 warned,
0 nonzero codes, 0 short on all 22 grid points; `re10`, `re100` 0
everywhere; the one-level design 40 of 40 at code 0; lmprec as round 1
(median 2.0e-7 SE at sd 1, at most 8.45e-5); `dev/optrap-rev-bitwise.R`
41 of 41 fits without a unit identical in hex to rellib-r7's (the MAD
moves one student fit across the unit-1 boundary each way: 41 no-unit,
1 with a unit, against round 1's 40 and 2).

The scale attack set (`dev/optrap-rev-scale.R`, 16 cases, 20 seeds;
`dev/optrap-p2-scale-sum.R`, `p2/scale-sum.txt` against
`p2/scale-sum-p1i.txt`, the same reference): unchanged except huber
(short 9 / 10 / 8 / 7 / 7 at 1e-6 to 1e6, round 1 9 / 10 / 9 / 11 /
10; largest gap 0.363 at 1e3, elsewhere at most 0.147), exgaussian
(2 nonzero codes at 1e3 and 1e6, round 1 3) and skew_normal at 1e6
(seed 7 now lands on the local optimum it reaches at every other scale
on every build, 2.81 below, at code 0, where round 1 stopped there at a
nonzero code with a warning). r5 is gone: test-v14's `asym_laplace`
(`dev/optrap-rev2-ald.R`, `p2/ald.txt`) now gets no unit and stops
where 0.69.0 does, 0.000184 below the exact maximum.

### r3 and r4: the Newton steps

`dev/optrap-rev2-polish.R` (`p2/polish-*.txt`): the saddle is left at
(1, 0.1), 0 steps; the bound case takes 0 steps on 10 of 10 with the
objectives of the arm without steps; `bign` 0 of 10 at n = 1e2, 1e4,
1e6. Predfix (`dev/optrap-rev2-predfix.R`, `p2/predfix.txt`): 3 steps,
62.578 gained, 4.7e-9 from glm(). The polynomial study (`p2/poly.txt`,
`p2/poly-yk16.txt`): degrees 5 and 6, 0 of 20 short at y and y x 16;
degree 7 1 of 20 (1.3e-4 and 1.1e-4 below lm(), warned; round 1 20 of
20, all warned). New tests
"a Newton step needs a positive definite Hessian" (the review's
saddle) and "a caller's optimizer tolerance leaves the fit where it
stopped" (`optimizer = "optim"` with `reltol` and with `factr`: no
steps; without them the steps run): both fail with round 1's
`newton_polish()`, pass on build 27.

### Filed

- **r2**: group sd 1000 times the residual sd (`re1000`, 20 seeds):
  2 of 20 (seeds 6 and 14) at "false convergence (8)" with a warning,
  at lmer's maximum to 1e-12; 0.69.0 0. The MAD of the residuals holds
  the group variance (unit 1024), so the group sd starts right and the
  coefficients step right, but sigma starts at the family's own
  `sd(y)`, 1000 times its value. Not fixed.
- **Open, high: a silent wrong answer on 0.69.0 and on this build.**
  Seed 1 of the wide design ends 9107.66 below lmer at code 0 with no
  warning (`dev/optrap-p2-wide1.R`, `p2/wide1-r7.txt`, rellib-r7): 3
  iterations, 9 evaluations, the group sd at 6.5e-16 against lmer's
  0.486, sigma 1.11 against 0.999, a largest gradient of 0.0062, and a
  restart of nlminb from there stops again (relative convergence,
  -2.9e-9). The start is the cause: sigma at `log(sd(y))` = 9.21 (sd
  1e4) and every slope at 0, so the first steps put the group sd at the
  boundary, where the log-scale gradient vanishes. The units do not
  catch it (the residual spread is 1, no unit), and `newton_polish()`
  does not run (random effects). Measured candidates
  (`dev/optrap-p2-wide1-starts.R`, `p2/wide1-starts-r7.txt`, through
  `start =`): sigma started at the residual MAD reaches lmer (-4.7e-8)
  but at code 1 with a warning; sigma at the MAD and the mu
  coefficients at least squares reach it at code 0, silent (-1.3e-8, 9
  iterations). That is a change of starting values for every
  location-scale fit whose response is far from its residual spread,
  across ten families' sigma conventions; not cheap, so filed for its
  own round.
- m4 (y + 1e9), mixture components, huber()'s jumps and the ledger row
  :326 as in round 1.

### Overlap

- Lane ci691 composes with no textual conflict (the review's trial
  merge); nothing done here.
- The `par_units` semantics written for lane seflags in round 1 hold;
  only how `u` is measured changed (MAD, rows above the cap).

### Reruns on the final build

`dev/optrap-scan.sh`, one test file per R process, every gate on.
`p2-ref`: the eight suites, reference BLAS. Its `test-se-tier.R` ran
before that file's sleepstudy fit moved to 32 ms (2 failures: the
MAD gives 64 ms a unit) and was rerun alone after the edit
(`p2-refsetier`, 91 pass); `p2-final` is the two merged. `p2-ob26`:
core under OpenBLAS 0.3.26 as BLAS and LAPACK, 4 threads, started
after the edit. Generated by `dev/optrap-counts.R p2-final p2-ob26`
(`dev/optrap-log/p2/counts.md`):

<!-- BEGIN GENERATED: dev/optrap-counts.R p2-final p2-ob26 -->

### p2-final

365 jobs; 365 files with a RESULT line; 365 whose `lib:` line names the lane's library.

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 224 | 17799 | 1 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1764 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 50 | 2693 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **365** | **24828** | **1** | **0** | **0** | **0** |

Files with a failure, an error, an escaped warning or no result: 1

- frmtmb--test-brms-suite-methods.txt: pass 161 fail 1 err 0 warn 0
    - FAIL test-brms-suite-methods.R:429 fitted has reasonable outputs brms brmsfit-methods:326 is recorded as pass and no longer holds: Expected `all(fi > 0)` to be TRUE. Differences: `actual`:   <NA> `expected`: TRUE  

### p2-ob26

224 jobs; 224 files with a RESULT line; 224 whose `lib:` line names the lane's library.

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 224 | 17800 | 1 | 0 | 0 | 0 |
| **total** | **224** | **17800** | **1** | **0** | **0** | **0** |

Files with a failure, an error, an escaped warning or no result: 1

- frmtmb--test-brms-suite-methods.txt: pass 161 fail 1 err 0 warn 0
    - FAIL test-brms-suite-methods.R:429 fitted has reasonable outputs brms brmsfit-methods:326 is recorded as pass and no longer holds: Expected `all(fi > 0)` to be TRUE. Differences: `actual`:   <NA> `expected`: TRUE  
<!-- END GENERATED -->

Fits that changed their optimum against 0.69.0
(`dev/optrap-fitdiff.R base-final p2-ref`, `dev/optrap-log/p2/fitdiff.txt`):

<!-- BEGIN GENERATED: dev/optrap-fitdiff.R base-final p2-ref -->
```
fits logged: base-final 9810 | p2-ref 9822 
test files: 308 | 308 

             base only             error both    lane higher (worse) 
                   986                    123                     10 
   lane lower (better)              lane only           same optimum 
                    21                    998                   8487 
unfitted or non-finite 
                   183 

changed optima by convergence codes:
                     
                      base code != 0 only both code 0 neither code 0
  lane higher (worse)                   0           1              9
  lane lower (better)                   2          16              3

changed fits (objective lane - base; negative is better):
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +94.56 (base -199.5743685 code 1 | lane -105.0174293 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +94.56 (base -199.5743685 code 1 | lane -105.0174293 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +7.245 (base -142.9797326 code 1 | lane -135.7343254 code 1) count ~ Trt + Age + mo(Exp) + s(Age) + volume + (1 + Trt | visit) + ar
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +0.0001231 (base 128.8283478 code 0 | lane 128.8284708 code 0) count ~ Trt + Age + mo(Exp) + s(Age) + volume + (1 + Trt | visit) + ar
lane higher (worse)  frmtmb.sample/test-brms-suite-methods.R    +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-data-helpers.R      +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-emmeans.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane lower (better)  frmtmb/test-optrap.R                       -118.1 (base 2913.226547 code 0 | lane 2795.149471 code 0) y ~ t + (1 + t | s) gaussian
lane lower (better)  frmtmb/test-predfix.R                      -62.58 (base 512.543671 code 0 | lane 449.9653873 code 0) y ~ 0 + xs poisson
lane lower (better)  frmtmb/test-predfix.R                      -62.58 (base 512.543671 code 0 | lane 449.9653873 code 0) y ~ 0 + xs poisson
lane lower (better)  frmtmb/test-se-check.R                     -24.78 (base 443.4558358 code 1 | lane 418.678587 code 1) ls ~ mo(inc) + age gaussian
lane lower (better)  frmtmb/test-optrap.R                       -18.68 (base 1543.327912 code 0 | lane 1524.65226 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -14.99 (base 1532.567176 code 0 | lane 1517.573737 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-predfix.R                      -12.08 (base 325.3373564 code 0 | lane 313.2592905 code 0) y ~ xr skew_normal
lane lower (better)  frmtmb/test-optrap.R                       -10.68 (base 1532.355266 code 0 | lane 1521.676117 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.8245 (base -553.9283063 code 1 | lane -554.7528471 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base -547.0088252 code 0 | lane -547.674324 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base -1237.784353 code 1 | lane -1238.449852 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base 834.5422026 code 0 | lane 833.8767318 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-se-tier.R                      -0.6164 (base 326.4954019 code 0 | lane 325.878983 code 0) y ~ x1 + x2 poisson
lane lower (better)  frmtmb/test-optrap.R                       -0.4487 (base 324.463624 code 0 | lane 324.0149446 code 0) ls ~ mo(income) * age gaussian
lane lower (better)  frmtmb/test-emmeans.R                      -0.08293 (base 1.321783367 code 0 | lane 1.23885053 code 0) y ~ f + gp(x) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.04888 (base 331.5744044 code 0 | lane 331.5255205 code 0) ls ~ mo(income) * age gaussian
lane lower (better)  frmtmb.sample/test-brms-pins.R             -0.0009362 (base 463.5842302 code 1 | lane 463.583294 code 1) y ~ x + (1 + x + z + w | g) gaussian
lane lower (better)  frmtmb/test-ordinal-mixture.R              -1.141e-05 (base 678.3999931 code 0 | lane 678.3999817 code 0) ysh ~ 1 mixture(cumulative, cumulative)
lane lower (better)  frmtmb/test-gp-by.R                        -9.74e-06 (base 53.90114814 code 0 | lane 53.9011384 code 0) ysig ~ x gaussian
lane lower (better)  frmtmb/test-v12.R                          -3.865e-06 (base 87.80208891 code 0 | lane 87.80208504 code 0) y ~ t2(Tmat, Xmat, by = Wmat) gaussian
lane lower (better)  frmtmb/test-review-v29.R                   -2.035e-06 (base 67.73199382 code 1 | lane 67.73199179 code 1) mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(), bf(y2 ~ x + (1 | q | id))

worse by more than 1e-6 relative: 3 
```
<!-- END GENERATED -->

Against round 1's list: test-v14's `asym_laplace` is no longer worse
(r5); fixture 1's refits on changed data, which run toward an
unbounded objective at code 1 on both builds, stop elsewhere (+94.56
twice and +7.25 where round 1 had +40.6 twice and +32.05); the
`iter.max = 3` fit of `test-se-check.R` stops at another point (code
1 on both); `test-predfix.R`'s skew_normal is 12.08 better (round 1
11.32), and two fits gain 1e-5 and 4e-6 (`test-ordinal-mixture.R`,
`test-v12.R`). Every file but the ledger row passes.

R CMD check (`dev/optrap-check.sh`, after the last edit, `dev/optrap-check/check.log`): Status: 1 NOTE (`checking HTML version of manual`); tests `[ FAIL 0 | WARN 0 | SKIP 383 | PASS 11921 ]`.

## Punch round 1

Review `dev/reviews/2026-10-08-optrap.md` (NOT MERGEABLE on B1). Final
build: install 26 (`dev/optrap-log-install26.txt`); every number in
this section is from that build unless it names another, with logs in
`dev/optrap-log/p1*` and `dev/optrap-log/rev-p1/p1i-*` (moved there from
the reviewer's directory in punch round 2; tag `p1i`; `p1`,
`p1f`, `p1g`, `p1h` are the intermediate builds 21, 23, 24 and 25).

### What changed in this round

- `R/autoscale.R`: `predictor_unit_of()` takes the unit from the
  residual spread, `predictor_resid_sd()` (new): the sd of `y` (log y
  for the lognormal families, the positive rows of a hurdle) after the
  fixed design and offset are taken out by least squares, on rows of
  positive weight, rounded to a power of two; a residual sd below
  1e-10 of the response's magnitude (an exact fit) gives no unit.
  `predictor_col_spread()`
  keeps spread 1 for a `mo()` placeholder whose codes times multiplier
  do not line up with the rows or hold an NA. Doc comments: the
  powers-of-two claim (m5), why huber() stays (m6), why mixtures stay
  out (m3).
- `R/fit.R`: `grad_units()` removed, `grad_verdict()` and
  `check_convergence()` back to `fit$par_units` (m1), `fit$resp_units`
  removed; new `newton_polish()`, called in `fit_assembled()` after
  `mo_search()` (m2, m10; below).
- Tests: `test-optrap.R` two new tests (high R^2, ill-conditioned
  polynomial) and the unit test names the residual sd;
  `test-grad-verdict.R` is 0.69.0's again except its ridge test (below);
  `test-se-check.R`'s m3 data are 0.69.0's again; `test-autoscale.R`
  one comparison made SE-relative.
- `NEWS.md`: the breaking bullet says residual spread, "to about 1e-7
  in the log-likelihood" (m5) and the 1e-5 SE precision (m2); a new bug
  fix bullet for the Newton steps.

### B1: the unit is the residual spread

The review's design (`y = 1 + k x + b_g + e`, seeds 1 to 40,
`dev/optrap-rev-r2b.R`), fits warned / code != 0 / short of lme4:

| k | round 0 (`OPTRAP_ROUND0=1`, `dev/optrap-log/p1-r2-round0.txt`) | final (`p1i-r2b.txt`) |
|---|---|---|
| 300 | 4 / 4 / 0 | 0 / 0 / 0 |
| 1e3 | 13 / 13 / 0 | 0 / 0 / 0 |
| 3e3 | 26 / 26 / 0 | 0 / 0 / 0 |
| 1e4 | 34 / 34 / 0 | 0 / 0 / 0 |

The review's `dev/optrap-rev-r2.R` (k 10, 100, 1e4): 0 short, 0
warned, 0 nonzero codes (`p1i-r2.txt`). New test "a covariate that
explains most of y leaves the fit silent" (`test-optrap.R`, seeds 1 to
3, k = 1e4): fails 6 assertions with the round-0 rule patched in
(`dev/optrap-testfile.R` with `OPTRAP_ROUND0=1`,
`dev/optrap-log/p1-t-optrap-round0.txt`), passes on the final build.

The scale study (`dev/optrap-scale.R`, `p1i-lanescale.txt`, and the
review's copy `p1i-scalecopy.txt`) is round 0's row for row: 0 short
in every row, invariance 1.68e-07, 6.24e-08, 1.05e-07, 1.9e-07 (ri20)
and 1.1e-08 to 1.14e-08 (slope). The `mo()` study (seeds 1 to 400,
`dev/optrap-p1-mo-sum.R`, `p1i-mo-sum.txt`): 200 / 200 and 200 / 200
at the exact maximum, 116783 and 102579 evaluations, as round 0.

The one-level group with a bound (m7): 40 of 40 at code 0
(`p1i-onelevel.txt`; round 0 5 of 40 at the evaluation limit): the
design's residual spread rounds to 1, so it has no unit. That is why
`test-se-check.R` m3 no longer halves its response (its data now have
no unit, halved they would get 0.5; `dev/optrap-p1-m3unit.R`).

An exact fit has no residual scale. Build 25's fit diff against
0.69.0 (`dev/optrap-fitdiff.R base-final p1h-ref`,
`dev/optrap-log/fitdiff-p1h.txt`) showed test-edgecases.R's constant
response (`y = 2`, `y ~ x + (1 | g)`) with residuals near 1e-16, a unit
of 4.4e-16, stopping at -6342.06 at code 0 where 0.69.0 runs on to
-19272.11 (code 1), both toward a sigma of 0 (`dev/optrap-p1-flat.R`,
`p1-flat-*.txt`). `predictor_resid_sd()` now gives no unit below
1e-10 of `max(abs(y))`; the build-26 fit is 0.69.0's (-19272.11, code
1). Two assertions of the units test (a constant response and one the
design reproduces exactly) fail on build 25
(`p1h-t-optrap-exactfit.txt`) and pass on build 26.

### m1: `grad_units()` is gone

The trip-wire reads the gradient in the full `par_units`.
`test-grad-verdict.R` is 0.69.0's file apart from the ridge test.
Student `sigma ~ x` at y x 1e-6 (the review's `dev/optrap-rev-warn1.R`):
the raw gradient 7.03 at (Intercept) is 1.3e-5 in natural units, no
warning; the scale attack set at 1e-6 warns on 0 of 20 student fits
(round 0 20) and 3 of 20 exgaussian (round 0 20; the 3 are at a
nonzero code).

### m2 and m10: Newton steps on the exact Hessian

Measured first: raw polynomials `y ~ x + ... + x^deg`, x on [1, 2],
residual sd 0.2, seeds 1 to 20 per degree, against lm()
(`dev/optrap-p1-polystudy.R`; logs `p1-poly-*.txt`). Short means more
than 1e-4 below lm(); silent means code 0 and no warning.

| arm | deg 5 short (silent) | deg 6 short (silent) | deg 6 min |
|---|---|---|---|
| rellib-r7, y | 0 | 3 (2) | -0.0218 |
| rellib-r7, y x 4 | 6 (6) | 19 (17) | -1.87 |
| rellib-r7, y x 16 | 20 (20) | 19 (19) | -1.87 |
| round-0 unit rule (unit 0.5) | 2 (2) | 14 (14) | -0.944 |
| residual unit, no Newton steps (build 23, unit 0.25) | 3 (3) | 19 (19) | -2 |
| final, y and y x 16 | 0 | 0 | -7.1e-07 |

So the shortfall is nlminb's relative function test on a Hessian of
condition number near 3e14, and how often it bites depends on the
response's scale; 0.69.0 was lucky at the data's own scale, and the
units, which make every scale behave like residual spread 1, inherit
the rate 0.69.0 has there. Ablation (`OPTRAP_POLY_ARM`): with the
column spreads off 15 of 20 at degree 6 stay short; with the
coefficient units off it is 0.69.0. At the stop the coefficients are
0.5 to 1 SE from lm() and one Newton step predicts the shortfall to 1
percent (`dev/optrap-p1-polish2.R`).

`newton_polish()`: for a fit with no random effects, at code 0, with
at most 100 outer parameters and no tolerance set in `optCtrl`, Newton
steps on `obj$he()` in natural units while a step predicts a gain
above 1e-6, kept only if the objective falls and the point stays
strictly inside the bounds; at most five. A fit that takes no step has
its objective state restored, so it is the fit 0.69.0 made bit for bit:
`dev/optrap-rev-bitwise.R` (`p1i-bitwise.txt`), 40 of 40 fits without
a unit identical in objective and every parameter to rellib-r7's (P
or E core logs of the review), the 2 student fits with a unit differ.
Degree 7 stays 20 of 20 short, all warned (the Hessian is singular in
double precision). New test "an ill-conditioned fixed-effect fit
reaches lm()'s maximum" (seeds 6001 and 6002 at y and y x 16): with
`newton_polish()` taken out (`OPTRAP_NOPOLISH=1`) 4 of 4 fail, on
rellib-r7 2 of 4 (y x 16), on the final build 0
(`dev/optrap-log/p1g-t-optrap*.txt`).

The tolerance gate: three `test-grad-verdict.R` tests loosen
`rel.tol` to make a fit stop short and assert the warning; the steps
took the shortfall away (7 failures on build 24, `p1-ref` scan). A
caller who sets a tolerance chose where nlminb may stop, so the steps
skip that fit. The ridge test (`x2 = x1 + 1e-5 noise`) stopped short
at the default control; it now states `rel.tol = 1e-10` to keep
nlminb's stop for the verdict's assertions, and asserts that with the
default control the fit is silent and gains at least half the
headroom the verdict predicted.

m2's precision (`dev/optrap-rev-lmprec.R`, `p1i-lmprec.txt`): median
max |beta - lm| / SE 2.0e-7 at sd 1, 3.0e-6 to 1.45e-5 at sd 1.5 to
1e6, largest 8.45e-5; log-likelihood within 5.1e-9. The steps do not
move these (their predicted gain is below 1e-6), so NEWS states the
precision. `dev/optrap-rev-prec.R` (`p1i-prec.tsv`): IQ mu 2.2e-7 SE and
log sigma 1.3e-7 SE (round 0 6.5e-6 and 4.8e-5), sleepstudy fixef vs lmer 1.7e-5 SE, ar1 ML 1.7e-5 to
2.4e-5 SE, ma1 ML up to 1.5e-5, cosy ML up to 6.1e-6; REML fits equal.
(Not pinned to a core type this time; the review found P and E agree
to all printed digits.)

m10, the review's probe (`dev/optrap-rev-att/probe-optrap.R poly`,
`p1i-probe-poly.txt`): "Large maximum absolute gradient" on 0 of 15
fits (round 0 7, rellib-r7 8). Degree 6 without a group term: SEs
within 1.5e-3 (x on [1, 2]) and 5.1e-5 (x on [1e5, 2e5]) of lm()'s;
build 23 had 0.054 at a fit 3.33 short. Degrees 7 and 8 without a group
term now lose every coefficient's SE as "flat" (the SE check, lane
seflags' territory) where round 0 warned about the gradient; with
`(1 | g)` degree 7 loses 3 of 10, as round 0. The last row (degree 8
with `(1 | g)`) is lmer's "Downdated VtV is not positive definite" on
every build. Filed below for the remainder.

### m3: mixtures stay out, measured

Attempt (build 22, `predictor_unit_of()` for `mu<k>` of a mixture
component in the ten families, with `sigma<k>`): the review's mixture
scale case went to 0 short at every scale (rellib-r7 and round 0:
0/0/0/4/20 short at 1e-6 / 1e-3 / 1 / 1e3 / 1e6; `dev/optrap-p1-mixsum.R`,
`p1-mixsum.txt`; the final build 0/0/0/1/20, `p1i-scale-sum.txt`). The 55 test files that call `mixture()`
(`dev/optrap-scan.sh p1mix`, `dev/ciharden-log/p1mix.sum`) failed in 4:

- `test-v18.R:199`: the group-level mixture with `(1 | g)` stops 111.2
  lower at "false convergence (8)" (`dev/optrap-p1-mix111.R`,
  `p1-mix111.txt`). The unit is 2: the residual spread holds the two
  components' separation (means -1 and 2, sd 0.8). The group sd starts
  at log 2, never moves (theta 0.693 at the stop), the modes sit at
  +-1.7 and the gradient is 4.45e8: a random intercept that large
  explains a group's class, the inner problem has a mode per class and
  the Laplace objective jumps between them. rellib-r7 starts the sd at
  1, runs it to 0 (singular) and converges. frmtmb.sample's
  `test-laplace-draws.R` has the same fit 1.04 lower.
- `test-brms-likelihood.R:715` (gradient 0.0013 against `tol_grad`
  1e-3), `frmtmb.latent test-hmm.R:187` (mixture against HMM to 8e-7),
  `frmtmb.eam test-defects.R:92` (residuals not finite on
  `mixture(wiener, lognormal)`, at an objective 3.26 lower): the
  precision or the path of the stop.

Reverted; mixtures keep their starts. Filed. Note that the Newton steps
now take mixture fits without random effects to their maximum where
nlminb stops short: at y x 1e3 the case is 1 of 20 short (round 0 4).

### m4: y + 1e9, filed

Unchanged (6 of 20 short at 1e9, 4 warned, 2 silent, every build).
The intercept already starts at mean(y) (`init_dpars` of gaussian), so
centering the start does nothing; at 1e9 each residual `y - mu` loses
about 7 of 16 digits, and the cure is centering `y` inside the
objective (an internal offset), which touches every gaussian path.
Filed.

### m5, m8, m9: the record

m5: `predictor_units()`'s comment now says the fit of `y * 2^k` is not
bitwise the scaled fit (13 of 30 on the same path to 1 ulp, 17 on
another; `p1i-pow2.txt` reproduces: 0 of 30 bitwise, 13 same
evaluations, 6.15e-6 relative in beta, 1.45e-7 in the log-likelihood);
NEWS says "to about 1e-7 in the log-likelihood". m8: STATUS and the
"Tests seen to fail" paragraph are corrected (round 0 counted 97 on a
build before its last edit; on build 26 `test-optrap.R` passes 122,
`dev/optrap-log/p1i-t-optrap.txt`, and the suites' counts are below).
m9: the
fixture-1 bullet says the objective is unbounded within reach (the
review's 41 jittered starts reach -1308.3; rerun, `p1i-fix1-sum.txt`:
default 109.3401042156 at code 1, :326 fails, 33 of 41 hold at sd 0.3,
24 of 41 at sd 1, best -1308.2774, as round 0).

### m6: huber() stays in

Its `k` is in units of sigma (`huber_rho((y - mu) / s, k)`), so the
unit does not change its meaning. Per seed (`dev/optrap-p1-huber2.R`,
`p1i-huber.txt`; short = 1e-4 below the best any arm reached):

| scale | rellib-r7 short / silent | final short / silent |
|---|---|---|
| 1e-6 | 12 / 0 | 9 / 5 |
| 1e-3 | 13 / 0 | 10 / 4 |
| 1 | 7 / 4 | 9 / 2 |
| 1e3 | 17 / 2 | 11 / 3 |
| 1e6 | 20 / 0 | 10 / 4 |

rellib-r7's largest gaps are 87.8, 47.1, 0.04, 0.98, 12.4; the final
build's 0.363 at most. The silent shorts are the same points on every
build (seeds 4, 12, 13, 18 at 0.0285, 0.0019, 0.0041, 0.0077, equal on
rellib-r7 at scale 1) plus seed 19 (0.363). Seed 19
(`dev/optrap-p1-huber19.R`, `p1i-huber19.txt`, at y x 1e3): the
final build stops at a gradient of 8.4e-5 (a stationary point);
rellib-r7's stop at scale 1, 0.32 higher, is at false convergence with
a gradient of 3.7 (the best any arm reached, 0.363 higher, is also a
code 1 stop). On the segment between the two the log-likelihood
moves in jumps of 0.14, 0.09, 0.09 and 0.08 among steps of 4e-4: the
Laplace objective of huber() is discontinuous, because the inner
Hessian is a step function of whether a residual is within `k`. The
"best" is the high side of a jump, not a maximum. The unit removes the
large failures and does not cause these; kept, with that comment at
`predictor_unit_families`.

### `par_units`, for lane seflags

`fit$par_units` is one positive number per outer parameter, aligned
with `names(fit$opt$par)`, or NULL when every entry is 1. Entry `i` is
the size of one natural step of parameter `i`, in that parameter's own
units: the optimizer works on `par / par_units` (nlminb's `scale` is
`1 / par_units`), a gradient in natural units is `g * par_units`, and
a Hessian in natural units is `diag(par_units) H diag(par_units)`. It
is the product of three factors, each 1 where it does not apply:

1. autoscale (`frmtmb_control(autoscale =)`): a column's 1 / sd for a
   badly scaled column (`autoscale_units()`);
2. the smooth null-space units (`smooth_fx_units()`);
3. this lane (`predictor_par_units()`): for a coefficient of the
   identity-link `mu` of the ten families with an estimated `sigma`,
   `u / s_j`, `u` the response's unit (`predictor_units()`, a power of
   two near the residual sd, in the response's units) and `s_j` the
   column's sd rounded to a power of two (1 for the intercept and a
   binary column, `sd(codes * mult)` for a `mo()` placeholder); `u`
   alone where factor 1 or 2 already applies.

Every non-`mu` parameter (sigma's, theta, zeta, cutpoints) stays 1.
Readers: nlminb's scale, the restart trip-wire (`optimize_obj()`), the
convergence trip-wire (`check_convergence()`, `grad_verdict()`), the
SE check (`fit_outer_hessian()`, tier 3, the edge screen),
`autoscale_sdreport()`, predict's estimability, `newton_polish()`. The
semantics are those of round 0; this round changed only how `u` is
chosen. The review's O1 (seflags' degree-7 thresholds on the merged
metric) still applies; see m10 above for the lane's own SE behavior on
those designs.

### Overlap with lane moparity

- `mo_search()`'s tail: take moparity's cox step 3 first, then this
  lane's `obj_settle(obj, opt, best, saved)` and `if (!runs) return(opt)`
  after it; drop moparity's `probed` restore, which `obj_settle()`
  covers (it restores the saved state when `best` is `opt` and checks
  it with `obj_state_at()`).
- `predictor_col_spread()` against moparity's frame
  (`frmtmb-wt-moparity/R/frame.R`, read only): each entry keeps `col`,
  `codes` (one per row) and `mult` (NULL for a main term, else
  `mm[, k]` of `mo_sp_design()`, one column per coded factor level, so
  `mox:fb` is its 0/1 contrast). On entries of that shape
  (`dev/optrap-p1-mospread.R`, `p1-mospread.txt`): `mox` spread 2,
  `mox:fb` and `mox:fc` 1, matching the spreads by hand. Because
  `mo_sp_design()` builds its model frame with `na.pass`, a missing
  factor value gives an NA multiplier; the guard added this round keeps
  that column's spread at 1 (also for a multiplier that does not line
  up with the rows) instead of an NA or recycled spread. No other guard
  is needed.

### Filed (left open after this round)

- Mixture components get no unit (m3, measured above).
- `y + 1e9` (m4).
- Raw polynomials of degree 7 and up: the Hessian is singular in
  double precision, the steps cannot help, and the fits are warned
  (by the SE check without a group term); seflags' O1 on the merged
  build stands.
- huber()'s discontinuous Laplace objective (m6): stops on either side
  of a jump, on every build.
- The ledger row `brmsfit-methods:326` (unchanged).

### Reruns on the final build

The studies and the review's scripts (`dev/optrap-p1-jobs.txt`, 40
jobs through `dev/optrap-rev-run.sh`, tag `p1i`; logs
`dev/optrap-log/rev-p1/p1i-*` (moved there from the reviewer's directory in punch round 2), summaries `dev/optrap-log/p1i-*`): every
output equals build 25's (`p1h`) apart from the tag.
Two jobs end at rc 1 on every
build, the review's own: `dev/optrap-rev-healthy.R` at its last
`print()` (after its four comparisons, which match round 0), and the
polynomial probe at lmer's error on its last row.

`dev/optrap-scan.sh`, one test file per R process, every gate on.
`p1i-ref`: the eight suites under the reference BLAS; its four
`test-scale.R` files of the extensions were cut off by a dropped
session and rerun alone as `p1i-refscale`; `p1i-final` is the two
merged (`dev/ciharden-log/p1i-final-suite/`). `p1i-ob26`: core under
OpenBLAS 0.3.26 as BLAS and LAPACK, 4 threads. Generated by
`dev/optrap-counts.R p1i-final p1i-ob26` (`dev/optrap-log/counts-p1i.md`):

<!-- BEGIN GENERATED: dev/optrap-counts.R p1i-final p1i-ob26 -->

### p1i-final

365 jobs; 365 files with a RESULT line; 365 whose `lib:` line names the lane's library.

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 224 | 17790 | 1 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1764 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 50 | 2693 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **365** | **24819** | **1** | **0** | **0** | **0** |

Files with a failure, an error, an escaped warning or no result: 1

- frmtmb--test-brms-suite-methods.txt: pass 161 fail 1 err 0 warn 0
    - FAIL test-brms-suite-methods.R:429 fitted has reasonable outputs brms brmsfit-methods:326 is recorded as pass and no longer holds: Expected `all(fi > 0)` to be TRUE. Differences: `actual`:   <NA> `expected`: TRUE  

### p1i-ob26

224 jobs; 224 files with a RESULT line; 224 whose `lib:` line names the lane's library.

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 224 | 17791 | 1 | 0 | 0 | 0 |
| **total** | **224** | **17791** | **1** | **0** | **0** | **0** |

Files with a failure, an error, an escaped warning or no result: 1

- frmtmb--test-brms-suite-methods.txt: pass 161 fail 1 err 0 warn 0
    - FAIL test-brms-suite-methods.R:429 fitted has reasonable outputs brms brmsfit-methods:326 is recorded as pass and no longer holds: Expected `all(fi > 0)` to be TRUE. Differences: `actual`:   <NA> `expected`: TRUE  
<!-- END GENERATED -->

Both: the ledger row :326 and nothing else. Build 24's run (`p1-ref`,
before the tolerance gate) had 7 failures in `test-grad-verdict.R`
and 2 in `test-autoscale.R`; both are addressed above and pass.

Fits that changed their optimum against 0.69.0
(`dev/optrap-fitdiff.R base-final p1i-ref`, `dev/optrap-log/fitdiff-p1i.txt`):

<!-- BEGIN GENERATED: dev/optrap-fitdiff.R base-final p1i-ref -->
```
fits logged: base-final 9810 | p1i-ref 9813 
test files: 308 | 308 

             base only             error both    lane higher (worse) 
                   982                    123                     11 
   lane lower (better)              lane only           same optimum 
                    19                    985                   8494 
unfitted or non-finite 
                   181 

changed optima by convergence codes:
                     
                      base code != 0 only both code 0 neither code 0
  lane higher (worse)                   0           1             10
  lane lower (better)                   2          14              3

changed fits (objective lane - base; negative is better):
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +40.6 (base -199.5743685 code 1 | lane -158.9725173 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +40.6 (base -199.5743685 code 1 | lane -158.9725173 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +32.05 (base -142.9797326 code 1 | lane -110.9339239 code 1) count ~ Trt + Age + mo(Exp) + s(Age) + volume + (1 + Trt | visit) + ar
lane higher (worse)  frmtmb/test-v14.R                          +0.003463 (base 910.304292 code 1 | lane 910.3077552 code 1) y ~ x asym_laplace
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +0.0001231 (base 128.8283478 code 0 | lane 128.8284708 code 0) count ~ Trt + Age + mo(Exp) + s(Age) + volume + (1 + Trt | visit) + ar
lane higher (worse)  frmtmb.sample/test-brms-suite-methods.R    +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-data-helpers.R      +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-emmeans.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt
lane lower (better)  frmtmb/test-optrap.R                       -118.1 (base 2913.226547 code 0 | lane 2795.149471 code 0) y ~ t + (1 + t | s) gaussian
lane lower (better)  frmtmb/test-predfix.R                      -62.58 (base 512.543671 code 0 | lane 449.9653873 code 0) y ~ 0 + xs poisson
lane lower (better)  frmtmb/test-predfix.R                      -62.58 (base 512.543671 code 0 | lane 449.9653873 code 0) y ~ 0 + xs poisson
lane lower (better)  frmtmb/test-se-check.R                     -55.11 (base 443.4558358 code 1 | lane 388.3498424 code 1) ls ~ mo(inc) + age gaussian
lane lower (better)  frmtmb/test-optrap.R                       -18.68 (base 1543.327912 code 0 | lane 1524.65226 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -14.99 (base 1532.567176 code 0 | lane 1517.573737 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-predfix.R                      -11.32 (base 325.3373564 code 0 | lane 314.0195391 code 0) y ~ xr skew_normal
lane lower (better)  frmtmb/test-optrap.R                       -10.68 (base 1532.355266 code 0 | lane 1521.676117 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.8245 (base -553.9283063 code 1 | lane -554.7528471 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base -547.0088252 code 0 | lane -547.674324 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base -1237.784353 code 1 | lane -1238.449852 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base 834.5422026 code 0 | lane 833.8767318 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-se-tier.R                      -0.6164 (base 326.4954019 code 0 | lane 325.878983 code 0) y ~ x1 + x2 poisson
lane lower (better)  frmtmb/test-optrap.R                       -0.4487 (base 324.463624 code 0 | lane 324.0149446 code 0) ls ~ mo(income) * age gaussian
lane lower (better)  frmtmb/test-emmeans.R                      -0.08293 (base 1.321783367 code 0 | lane 1.23885053 code 0) y ~ f + gp(x) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.04888 (base 331.5744044 code 0 | lane 331.5255205 code 0) ls ~ mo(income) * age gaussian
lane lower (better)  frmtmb.sample/test-brms-pins.R             -0.0009362 (base 463.5842302 code 1 | lane 463.583294 code 1) y ~ x + (1 + x + z + w | g) gaussian
lane lower (better)  frmtmb/test-gp-by.R                        -9.74e-06 (base 53.90114814 code 0 | lane 53.9011384 code 0) ysig ~ x gaussian
lane lower (better)  frmtmb/test-review-v29.R                   -2.035e-06 (base 67.73199382 code 1 | lane 67.73199179 code 1) mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(), bf(y2 ~ x + (1 | q | id))

worse by more than 1e-6 relative: 4 
```
<!-- END GENERATED -->

Reading it against round 0's: of the 11 worse, 10 are fixture 1 and
its refits (9 at code 1 on both builds; one at code 0, 1.2e-4
higher), as before, and the 11th is `test-v14.R`'s `asym_laplace` `y ~ x` (+0.0035, code 1 on
both: a kinked likelihood). New among the 19 better, from the Newton
steps: `test-predfix.R`'s poisson `y ~ 0 + xs` with x at 1e-6 and
autoscale off, 62.58 better at code 0 on both builds (0.69.0 stops
there 62.58 short in silence), its skew_normal `y ~ xr` (-11.32, code
0 on both), and `test-se-tier.R`'s RB1 poisson pair at a correlation
near 1 - 5e-11 (-0.62, code 0 on both). Every one of these files
passes.

R CMD check (`dev/optrap-check.sh`, after the last edit, log `dev/optrap-check/check.log`): Status: 1 NOTE (`checking HTML version of manual`, the environmental note); tests `[ FAIL 0 | WARN 0 | SKIP 383 | PASS 11912 ]`.

## Round 0 (the record before the review; the sections below)

## What changed

Edited in `R/fit.R`: `fit_assembled()`, `quad_fit()`,
`fit_recovery_starts()`, `optimizer_from_best()`, `optimize_obj()`,
`obj_settle()`, `mo_search()`, `make_start()`; new in `R/fit.R`:
`obj_state_at()`. New in `R/autoscale.R`: `predictor_units()`,
`predictor_unit_families`, `predictor_unit_of()`,
`predictor_resid_sd()`, `unit_theta_start()`, `predictor_par_units()`,
`predictor_par_units_resp()`, `predictor_col_spread()`,
`predictor_inner_control()`. Nothing in `R/objective.R`. Tests: new
`tests/testthat/test-optrap.R`; edited `test-id-kron.R`,
`test-se-tier.R`, `test-se-check.R`, `test-bcm-gaussian.R`,
`test-glmm-gaussian.R`, `test-autocor.R`. `NEWS.md`: a development
section. (Punch round 1 took out `grad_units()` and the round-0 edits
of `grad_verdict()`, `check_convergence()` and `test-grad-verdict.R`.)

1. **The response's unit** (item 1). The identity-link `mu` of the ten
   core location-scale families (`predictor_unit_families`: gaussian,
   student, lognormal, shifted_lognormal, hurdle_lognormal,
   skew_normal, exgaussian, asym_laplace, zero_inflated_asym_laplace,
   huber) gets a unit `u`: since punch round 1 the sd of the response
   (log y for the lognormal families) after its fixed design and
   offset are taken out by least squares, rounded to a power of two,
   computed on the rows of positive weight (round 0 used the family's
   own start of `sigma`, `sd(y)` for gaussian). Three things use it:
   - `unit_theta_start()` (called by `make_start()`): every
     random-effect sd of a block whose columns enter that `mu` starts
     at `log(u)` instead of 0 (shared sds only when all columns agree;
     mapped entries never).
   - `predictor_par_units()` (in `fit_assembled()`): the coefficients
     of that `mu` get the natural unit `u / s_j`, `s_j` the column's
     sd rounded to a power of two (1 for the intercept and binary
     columns; `sd(codes * mult)` for a `mo()` coefficient; `u` alone
     for a column autoscale or the smooth null space already scales).
     Multiplied into `par_units`, so nlminb's `scale`, the restart
     trip-wire in `optimize_obj()`, the convergence check's trip-wire
     (`grad_verdict()`, `check_convergence()`) and the standard-error
     check's Hessian all read them.
   - `predictor_inner_control()`: TMB's inner Newton `grad.tol` is
     divided by the largest unit (passed to every `MakeADFun()` of the
     fit: `fit_assembled()`, `quad_fit()`, `fit_recovery_starts()`).
   A response with spread in [0.71, 1.41) has `u = 1`: no change.
   Mixtures, extension families and a `sigma` held at a constant
   (`se()`, a fixed `sigma`) are left out (measured below).
2. **The restart that raises** (item 2). `optimizer_from_best()`
   checks that the objective and gradient are finite at
   `last.par.best` before restarting there, and re-raises the original
   error otherwise; `optimize_obj()` catches a restart that raises,
   keeps the run before it, and settles the objective there.
3. **Objective state after a discarded run** (item 4).
   `optimize_obj()` saves the state before a restart and calls
   `obj_settle()` when the restart raises or ends worse; `obj_settle()`
   now verifies that a restored state's last and best points are the
   kept optimum (`obj_state_at()`) and evaluates there if not;
   `mo_search()` settles whether or not it refitted anything.
4. **The `mo()` cone search** (item 5). The other cone's search starts
   at least one of the coefficient's units from 0 (`par_units`), tries
   the next vertex by score while the coefficient ends on its bound,
   and searches both cones when the coefficient is within 1e-3 of its
   unit of 0.
5. **Item 3** is the P/E-core split; the test is made robust to it.

Tried and taken out (each measured; details in the sections below):
the objective's additive constant `n log u` subtracted for nlminb; a
check of every `mo()` face for a descent direction; every vertex as a
start of the cone search; a Newton-step start of the cone search; and
a rule that kept a converged run over a restart that ended at a
nonzero code without gaining more than 1e-10 of the objective. That
rule was written for `mo()` fits at their exact maximum that a restart
turned into "false convergence (8)" on an intermediate build (seed 3
of dev/optima-mo-study.R); on the final build it changes none of the
200 fits of that study (`OPTRAP_NOKEEP=1`, logs
`dev/optrap-log/mo-nokeep/`: 200 of 200 at the maximum, 0 nonzero
codes, the same 116783 evaluations), and under OpenBLAS's LAPACK it
kept test-bcm-simple.R's SIMPLE_1 at code 0 with a gradient of 0.00251
where rellib-r7 keeps the restart's 1.9e-4, which the file's Stan
identity check (`tol_grad = 1e-3`) needs (`dev/optrap-simple.R`).

## Item 1: the response's units

### The cause, measured on rellib-r7

`dev/optrap-diag1.R` and `dev/optrap-diag2.R` (ri20 design, seeds 1 to
40; logs `dev/optrap-log-diag1-r7.txt`, `-diag2-r7.txt`). Short means
more than 1e-4 below lme4's logLik.

| arm | x1e-3 short | x1 short | x1e3 short |
|---|---|---|---|
| default (rellib-r7) | 4 (min -0.8245) | 0 | 20 (min -0.6655) |
| group log sd started at log(sd(y)) only | 0 | 0 | 25 (min -1.677) |
| mu coefficients stepped in sd(y) only | 5 (min -0.7156) | 0 | 22 (min -1.677) |
| both | 0 (min -7.5e-08) | 0 | 0 (min -1.1e-07) |

Two things are scale-dependent and both are needed: the group sd
started at 1 in the response's units (1e-3 or 1e3 of the residual
sd), and nlminb stepped the mu coefficients in absolute units. A third
shows at 1e6 (below): TMB's inner Newton stops when the largest
absolute inner gradient is under 1e-8, which in units 2^20 leaves the
modes measurably unconverged. The objective's additive constant
(`n log u`) was also tried: subtracting it made the whole study
invariant to 1e-8 instead of 1e-7, but it makes nlminb's relative
function test stricter whenever `|f|` shrinks, and it turned converged
`mo()` fits into "false convergence (8)" (seed 54 of
`dev/optima-mo-study.R`); it was taken out.

### The scale study (dev/optrap-scale.R)

Seeds 1 to 40, two designs: `ri20` (`y ~ x + (1 | g)`, the setier
design) and `slope` (`y ~ t + (1 + t | s)`, 18 subjects of 10 days,
sleepstudy-like). Columns: fits more than 1e-4 below lme4 (min), fits
at a nonzero code, and the largest `|logLik(y*k) + n log k -
logLik(y)|`. Logs `dev/optrap-log-scale-r7.txt` (rellib-r7) and
`dev/optrap-log/scale-final.txt` (the lane's final build).

| design | scale | rellib-r7 short | r7 code != 0 | r7 invariance | lane short | lane code != 0 | lane invariance |
|---|---|---|---|---|---|---|---|
| ri20 | 1e-6 | 6 (-0.665) | 40 | 0.665 | 0 (-1.81e-07) | 0 | 1.68e-07 |
| ri20 | 1e-3 | 4 (-0.825) | 9 | 0.825 | 0 (-7.5e-08) | 0 | 6.24e-08 |
| ri20 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| ri20 | 1e3 | 20 (-0.665) | 0 | 0.665 | 0 (-1.22e-07) | 0 | 1.05e-07 |
| ri20 | 1e6 | 40 (-23.4) | 0 | 23.4 | 0 (-2.1e-07) | 0 | 1.9e-07 |
| slope | 1e-6 | 0 | 40 | 4.9e-09 | 0 (-1.01e-08) | 0 | 1.1e-08 |
| slope | 1e-3 | 0 | 13 | 4.3e-09 | 0 (-2.83e-09) | 0 | 6.14e-09 |
| slope | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| slope | 1e3 | 17 (-8.62) | 0 | 8.62 | 0 (-3.18e-09) | 0 | 5.49e-09 |
| slope | 1e6 | 40 (-127) | 0 | 127 | 0 (-7.97e-09) | 0 | 1.14e-08 |

Ablations on the final build (`OPTRAP_ABLATE=` replaces one part by
the identity; logs `dev/optrap-log/scale-final-abl-*.txt`):

| without | ri20 1e-6 | ri20 1e-3 | ri20 1e3 | ri20 1e6 | slope 1e-6 | slope 1e3 | slope 1e6 | code != 0 |
|---|---|---|---|---|---|---|---|---|
| the theta start | 12 (-0.825) | 8 (-0.716) | 22 (-1.68) | 22 (-1.68) | 2 (-92.8) | 0 | 40 (-104) | 0 |
| the coefficient units | 0 | 0 | 24 (-1.68) | 40 (-23.4) | 0 | 1 (-0.041) | 40 (-26.1) | 0 |
| the inner tolerance | 0 | 0 | 0 | 0 (-2.6e-06) | 0 | 0 | 0 (-1.28e-06) | 6 and 3 at 1e6 |

A plain linear model is in the trap too (`dev/optrap-se1e6.R`,
`y ~ x`, seeds 1 to 10; log `dev/optrap-log/se1e6-final.txt`):
rellib-r7 ends short of `lm()` on 1 of 10 at
x1e3 (by up to 8.6e-4) and on 10 of 10 at x1e6 (by up to 23.3), all
without a warning; the lane on 0 of 10 (largest gap 1.8e-8 at 1e6).

### Is a response near unit spread unchanged? (the absent case)

`predictor_units()` is empty for `u = 1`, for a log-link family, for a
lognormal response multiplied by 1e3 (its sigma is `sd(log y)`), and a
`sigma ~ x` predictor never gets one (test "only a predictor in the
response's units gets a unit"). With no unit, `par_units`,
`make_start()` and `inner.control` are what rellib-r7 used, bit for
bit. On the eight suites (below) every fit whose objective moved by
more than 1e-8 relative is of one of the ten families with a unit; no
poisson, binomial, Gamma, ordinal, mixture or extension-family fit
moved.

### Mixtures and extension families are left out

The first version gave a unit to every `mu*` with a `sigma*`. On the
eight suites (`dev/optrap-fitdiff.R base-ref lane-ref`, log
`dev/optrap-log/fitdiff-ref1.tsv`) it moved 68 fits to a lower maximum
and 38 to a higher one; among the lower: test-v18's
`mixture(gaussian, gaussian)` ran a component's sigma to 0 (111 units
lower, code 1), frmtmb.latent's hidden Markov models reached the
optimum `hmm_starts()` is tested to escape from (+22.7 and +14.5 on 31
fits of `test-hmm-starts.R`, +17.1 in `test-hmm.R`), frmtmb.learn's
families (+4.4, +1.16). Restricted to the ten core families, all of
these return to rellib-r7's objective (the rerun
`dev/optrap-scan.sh lane-sub1`, 17 files).

A `sigma` held at a constant is left out too. With a unit, the fuzz
specs `student()` with `se(sdy)` (`test-fuzz.R`; identified with
`dev/optrap-scan.sh base-fuzz` and `lane-fuzz`, which log each fit's
formula and family) stopped elsewhere, all without converging on one
build or both: `+ gp(xs, k = 10) + (1 + x || g)` at 301.68 against
237.01 (code 1 on both), `+ mo(xo) * z + (1 | ga/gb)` at 168.65
against 162.56 (code 1 on both), and `+ (1 + x | g)` at 190.05 (code
1) against 194.51 (code 0). With `se()` the residual scale is given,
not estimated, so there is no estimated sigma whose start is the
unit.

## Item 2: the restart that raises

Fuzz seed 20379118 (`dev/optrap-fuzz1.R`, logs
`dev/optrap-log-fuzz1-{ref,ob26}-{r7,lane}.txt`):

| build | reference BLAS | OpenBLAS 0.3.26 BLAS and LAPACK, 4 threads |
|---|---|---|
| rellib-r7 | code 1, -4074.5359144421 | ERROR "NA/NaN gradient evaluation" |
| lane | code 1, -4074.5359144421 | code 1, -49.1818324061 (the first run), 2 warnings |

Why `last.par.best` is not a usable start (`dev/optrap-fuzz2.R`, log
`dev/optrap-log-fuzz2-ob26.txt`, lane build under OpenBLAS): the tape
recorded `value.best` -1675.96303402 during the restart that raised;
evaluated again at that point the objective is NaN and all 8 gradient
entries are non-finite, because the inner problem is solved from
another start. The REML likelihood here has no bound (the reference
build runs it to -4074.5). TMB keeps the lowest objective it has seen,
not a point the optimizer can restart from. The objective at the kept
optimum is not reproducible either: -49.0733 evaluated after the
failure against -49.1818 reported, so the lane restores the state
saved after the first run (`obj_settle()`) rather than re-evaluating.
The release review's candidate (`obj$fn(opt$par)` after a failed
restart) left `last.par.best` at the NaN point, where `sdreport()`
reads it; it gave -54.0241 on this seed through a different mo()
search path.

The 300-seed sweeps: see "Fuzz sweep" below.

## Item 3: test-id-kron.R's two optima

It is the core type (`dev/rtmb-pitfalls.md` item 21). The merged-dpar
fit pinned with `start /affinity` (`dev/optrap-idk.R`,
`dev/optrap-idk-run.cmd`; logs `dev/optrap-log/idk-*.txt`):

| build | P cores (mask C03C03) | E cores (mask 3FC3FC) |
|---|---|---|
| rellib-r7 | 6 of 6: 128.07435688, boundary message | 6 of 6: 128.074356883, "singular convergence (7)" warning |
| lane | 6 of 6: 128.07435688, boundary message | 6 of 6: 128.074356883, "singular convergence (7)" warning |

The objective and gradient at the start are identical to the bit on
both core types (`f0 0x1.4195731e98f9cp+8`); at four points off the
start they differ in the last ulps, and so does the joint
negative log-likelihood with the random effects held, which no inner
solve touches. The difference is in the taped evaluation itself, which
for this `gr(cov = A)` block runs dense products over the 128-level
relationship matrix; a kernel chosen per process from the core's
reported cache sizes (Eigen's blocking does this) fits every
observation. The end points differ along the flat direction the
correlation runs on (its unconstrained parameter 499.17 on P, 499.80
on E); on P the SE check calls it a boundary, on E it finds no lost SE
and the code 7 warning stands. frmtmb cannot pin a kernel below RTMB,
so the result is made robust instead: the test asserts that the fit
says so once, by the boundary message or the singular-convergence
warning (`allow_boundary()` and `allow_warnings()` around a counting
handler), and passes 3 of 3 pinned to P cores and 3 of 3 pinned to E
cores (`dev/optrap-pin-test.cmd`, logs `dev/optrap-log/pin-idk*.txt`).

## Item 4: the objective's state after a discarded run

Where a run is discarded and what now settles it:

| path | before | now |
|---|---|---|
| `escape_stationary()` keeps another run | `obj_settle()` | `obj_settle()` |
| `mo_search()` keeps another run | `obj_settle()` | `obj_settle()` |
| `mo_search()` evaluates and refits nothing | not settled | `obj_settle()` |
| `optimize_obj()` restart ends worse | not settled | `obj_settle()` |
| `optimize_obj()` restart raises | the error propagated | `obj_settle()` on the kept run |
| restored state not at the kept optimum | trusted | evaluated there (`obj_state_at()`) |

Measured on rellib-r7 and the lane (`dev/optrap-state.R`): after an
ordinary nlminb run the last and best points are the reported optimum
on all five families tried (sleepstudy, poisson and bernoulli GLMMs,
ri20, a skew normal with an escape), and `ranef()` equals the modes
at the optimum to 0; so no settle is added at the end of every run.

## Item 5: mo() seed 194, fit_loss2, F11 seed 8

`dev/optrap-mo-study.R` (dev/rel069-mo-study.R with arms r7 and
optrap), `ls ~ mo(income) * age`, exact maximum by quadratic
programming over the sign patterns. Logs `dev/optrap-log/mo-r7-*.tsv`,
`dev/optrap-log/mo-holdout/`, `dev/optrap-log/mo-final*/`.

| seeds | build | at the exact max (gap < 1e-6) | short | code != 0 | warned | evaluations |
|---|---|---|---|---|---|---|
| 1 to 200 | rellib-r7 | 199 | 194 (0.0489) | 0 | 0 | 85701 |
| 1 to 200 | lane | 200 | none (max 1.2e-07) | 0 | 0 | 116783 |
| 201 to 400 | rellib-r7 | 194 | 212, 214, 273, 279, 354 (0.449), 394 | 0 | 1 | 87398 |
| 201 to 400 | lane | 200 | none (max 9.1e-09) | 0 | 0 | 102579 |

Seeds 201 to 400 were not used to build the search. The path to it,
all measured on seeds 1 to 200 (logs in `dev/optrap-log/mo-*/`, one
directory per build): the response's unit alone on the coefficients
left 6 short, up to 0.95, and 20 at a nonzero code (`mo-units`); a
check of every face for a descent direction found none on those seeds
(the cone's maximum is unique and the faces were cone maxima) and was
removed (`mo-face`); a start at each vertex's own Newton step in b
left 7 short (`mo-newton`), since with the other parameters held every
vertex favors the current sign; the column spread in the units left 10
short and none at a nonzero code (`mo-spread-newton`); starting the
cone search one unit from 0 left seeds 184 and 194 (`mo-unitstart`);
the next vertex on a collapse reached 200 of 200; rounding the column
spreads to powers of two then left seed 80 short (`mo-pow2`: an
interior simplex, coefficient 6e-5); searching both cones near 0
reached 200 of 200 again (`mo-final2`). Every vertex as a start (D
refits per face fit, `mo-vloop`) also fixed seed 194 but left 37 and
189 and cost 6.1 times the evaluations; it was taken out.

fit_loss2 (`dev/optrap-loss2.R`, nonlinear, no unit) and F11 seed 8
(`dev/optrap-f11.R`) are unchanged: the default fits are -358.545112
(code 1) and -589.0237261 on both builds, 0.032 below the best
jittered start that reaches a maximum (-358.5132; one start ends at
-355.47 with code 1 and a gradient of 1e8, a degenerate point, which
the script's own "below the best start by 3.074" line counts) and
0.1418 below the best of F11's multi-start, as lane optima recorded
(logs `dev/optrap-log/loss2-*.txt`, `f11-*.txt`, measured with the
keep rule in; neither model has a unit or a `mo()` term). Their
cause is not shared: fit_loss2's better optimum is a jittered start of
a nonlinear model, F11's a smoothing sd started at the boundary. Filed
again below.

## Found while doing this

- **The one-level group with a bound** (`test-se-check.R` m3,
  `dev/optrap-onelevel.R`, logs `dev/optrap-log/onelevel-*.txt`):
  `y ~ x + (1 | g1)`, g1 one level, `ub = 0.1` on x, unit 2. nlminb
  zigzags between log sigma and the creeping log sd to its evaluation
  limit on 5 of 40 seeds on the lane, 0 of 40 on rellib-r7 (0 of 40 on
  both without the bound). PORT's own combination (`dev/optrap-zigzag.R`
  straight through nlminb): only the log-sd start log 2 together with
  the coefficient scale 1/2 zigzags; a restart in raw units recovers 1
  of the 5. The model is degenerate (the variance is not identified
  and frmtmb warns so); the test's data are halved so that its unit is
  1, and this is a regression of the lane on this design, not fixed.
  (Punch round 1: with the residual unit the test's own data have unit
  1, and the halving is taken out; see that section.)
- **Fixture 1 of the ported suite** (student, `mo()`, `s()`, `arma`,
  `sigma ~ Trt`, unit 8): both builds stop at "false convergence (8)"
  with an sd at -723 to -1253 on the log scale. The objective is not
  flat there but unbounded within reach: from one of 41 jittered starts
  the lane runs it to -1308.3 (`dev/optrap-rev-fix1.R`, the review's
  item 6); the lane's stop is
  4.5e-5 higher in the objective and its sigma coefficients' variances
  in `sdreport()` are negative, so `fitted(fit1, dpar = "sigma")` has
  NaN errors and the ledger row `brmsfit-methods:326` no longer holds
  (`dev/optrap-fix1.R`). Lane optima's change made that row hold; it
  follows the path of a non-converged fit.
- **Precision tests.** Five assertions compared estimates to an
  external optimum at an absolute 1e-6 or 1e-4 that neither optimizer
  promises; rellib-r7's path met them. On the lane the objectives
  agree with rellib-r7's to 1e-10 relative (79 fits of those files:
  20 lower, 28 higher, 31 equal) and the parameters move by up to
  1.1e-5 (log sigma of the IQ fit) and 1.1e-4 (sleepstudy's intercept
  against lmer). The assertions now ask for agreement within 1e-3 of
  the fit's own standard error (`test-bcm-gaussian.R`,
  `test-glmm-gaussian.R`, `test-autocor.R`).
- The restart trip-wire and nlminb's `scale` now read units on every
  location-scale fit with a unit; the convergence check's trip-wire
  does not (`grad_units()`): on a build whose units were the
  response's alone, the converged `mo()` fit of `helper-mo-flat.R`
  (unit 16, age near 40) read 0.00116 at age and warned about
  convergence in place of the SE warning (`test-se-check.R`, 9
  assertions). (Punch round 1: superseded. With the column spread in
  `par_units` the full units pass that file, and `grad_units()` gave
  false warnings at small units; it is gone.)

## Tests seen to fail on rellib-r7

`dev/optrap-testfile.R "" frmtmb tests/testthat/test-optrap.R`, log
`dev/optrap-log/t9-optrap-r7.txt` (lane: `t9-optrap-lane.txt`, pass 97,
fail 0, error 0, warn 0, on a build before round 0's last edit; the
review counts 99 on round 0's final build; punch round 1's numbers are
in its section); every test of the new file fails or
errors on rellib-r7 and passes on the lane:

| test | rellib-r7 |
|---|---|
| a random intercept fit is the same at every scale of y | 26 of 64 fail, 10 warnings |
| a random slope fit is the same at every scale of y | 4 of 20 fail, 2 warnings |
| only a predictor in the response's units gets a unit | error (no `predictor_units()`, the weak form; the two tests above are its behavior) |
| a restart that raises keeps the run before it | error "NA/NaN gradient evaluation" leaves `optimize_obj()` |
| a restart that ends worse leaves the objective at the kept run | 2 fail: last and best point are the restart's |
| the other cone's search reaches its maximum | 2 of 8 fail: seeds 194 and 354 short |

## Fuzz sweep (item 2)

`dev/optrap-sweep.sh 1 300 spec mo` (dev/relrev069-fuzzsweep.R's
generator from this worktree, `dev/optrap-fuzzsweep.R`; seeds
`20260901 + k * 977`, k = 1 to 300; summary `dev/optrap-sweep-sum.R`,
logs `dev/optrap-log/sweep/`):

<!-- BEGIN GENERATED: dev/optrap-sweep-sum.R -->
```
mo    ob26  seeds 300 | errors r7 35 lane 23 | r7 errors where lane fits 12 | lane errors where r7 fits 0
mo    ref   seeds 300 | errors r7 35 lane 23 | r7 errors where lane fits 12 | lane errors where r7 fits 0
spec  ob26  seeds 300 | errors r7 42 lane 31 | r7 errors where lane fits 11 | lane errors where r7 fits 0
spec  ref   seeds 300 | errors r7 37 lane 28 | r7 errors where lane fits 9 | lane errors where r7 fits 0
```
<!-- END GENERATED -->

`ob26` is OpenBLAS 0.3.26 as BLAS and LAPACK, 4 threads
(`dev/ciharden-openblas.sh 0.3.26 lapack`, build log
`dev/optrap-log-openblas-build.txt`: Rblas 74 symbols to OpenBLAS,
Rlapack 616). The lane removes 9 to 12 errors per 300 and adds none,
the release review's figures for its candidate (28 and 31 for the
spec, 23 and 23 for `mo(xo)`). The remaining errors are fits whose
first run itself fails: on the four traced (`dev/optrap-fuzz1.R <lib>
<seed>`, seeds 20266763, 20267740, 20303889 and 20335153) the first
run raises "NA/NaN gradient evaluation" before it returns, its best
point is no usable start, and the recovery starts raise too.

## The eight suites

`dev/optrap-scan.sh` (dev/ciharden-scan.sh with the fit logger
`dev/optrap-profile.R`): one test file per R process, every gate
(`NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS`,
`FRMTMB_FUZZ`, `FRMTMB_SCALE_TESTS`). Runs: `base-final` (rellib-r7,
reference BLAS, this tree's tests), `lane-final3` (the final build,
reference BLAS), `lane-ob26c` (the final build, OpenBLAS 0.3.26 as
BLAS and LAPACK, 4 threads). The extensions load from rellib-r7 with
the lane's core (no extension changed). Generated by
`dev/optrap-counts.R` (`dev/optrap-log/counts-final.md`); the base
run's per-expectation lines are left out here (35 failures and 5
errors, all in this lane's new or rewritten tests, which read lane
functions or assert what rellib-r7 does not do).

<!-- BEGIN GENERATED: dev/optrap-counts.R lane-final3 lane-ob26c -->
### lane-final3

365 jobs; 365 files with a RESULT line; 365 whose `lib:` line names the lane's library.

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 224 | 17766 | 1 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1764 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 50 | 2693 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **365** | **24795** | **1** | **0** | **0** | **0** |

Files with a failure, an error, an escaped warning or no result: 1

- frmtmb--test-brms-suite-methods.txt: pass 161 fail 1 err 0 warn 0
    - FAIL test-brms-suite-methods.R:429 fitted has reasonable outputs brms brmsfit-methods:326 is recorded as pass and no longer holds: Expected `all(fi > 0)` to be TRUE. Differences: `actual`:   <NA> `expected`: TRUE  

### lane-ob26c

365 jobs; 365 files with a RESULT line; 365 whose `lib:` line names the lane's library.

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 224 | 17767 | 1 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1764 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 50 | 2693 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **365** | **24796** | **1** | **0** | **0** | **0** |

Files with a failure, an error, an escaped warning or no result: 1

- frmtmb--test-brms-suite-methods.txt: pass 161 fail 1 err 0 warn 0
    - FAIL test-brms-suite-methods.R:429 fitted has reasonable outputs brms brmsfit-methods:326 is recorded as pass and no longer holds: Expected `all(fi > 0)` to be TRUE. Differences: `actual`:   <NA> `expected`: TRUE  
<!-- END GENERATED -->

base-final: 365 of 365 files, 24741 pass, 35 fail, 5 error, 12
escaped warnings, all in `test-optrap.R` (34 fail, 2 error, 12 warn,
its tests seen to fail) and `test-grad-verdict.R` (3 errors: the
rewritten tests call `grad_units()` and `predictor_units()`). Every
other rewritten test passes on rellib-r7 too.

Reading the lane's blocks:

- Reference BLAS: 24795 pass, 1 failure, 0 errors, skips or escaped
  warnings. The failure is the ledger row `brmsfit-methods:326`
  ("Left open").
- OpenBLAS 0.3.26 BLAS and LAPACK: the same failure and nothing else.
  The two results 0.69.0 fails there (`dev/round-20261007.md`) pass:
  `test-fuzz.R`'s seed 20379118 (item 2), and `data-helpers:7`, but
  the second only because fixture 1 no longer converges under
  OpenBLAS's LAPACK on this build, so the prediction test that gave
  the false alarm is not reached; the false alarm itself (the lost
  basis restricted to named parameters) is lane seflags' and not
  touched here.

## Which fits changed their optimum

`dev/optrap-fitdiff.R base-final lane-final3` (log
`dev/optrap-log/fitdiff-final3.txt`, table `fitdiff-final3.tsv`):
every `fit_assembled()` call of the eight suites, matched within its
test file by a digest of the response, the call, REML and quadrature
and by occurrence. A change is an objective difference above 1e-8
relative.

<!-- BEGIN GENERATED: dev/optrap-fitdiff.R -->
```
fits logged: base-final 9810 | lane-final3 9811 
test files: 308 | 308 

             base only             error both    lane higher (worse) 
                   781                    123                     11 
   lane lower (better)              lane only           same optimum 
                    24                    782                   8688 
unfitted or non-finite 
                   183 

changed optima by convergence codes:
                     
                      base code != 0 only both code 0 neither code 0
  lane higher (worse)                   0           1             10
  lane lower (better)                   2          16              6

changed fits (objective lane - base; negative is better):
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +89.44 (base -247.8022574 code 1 | lane -158.363546 code 1) count ~ Trt + mo(Exp) + s(Age) + volume + (1 +
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +40.6 (base -199.5743685 code 1 | lane -158.9725173 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume 
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +40.6 (base -199.5743685 code 1 | lane -158.9725173 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volume 
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +32.05 (base -142.9797326 code 1 | lane -110.9339239 code 1) count ~ Trt + Age + mo(Exp) + s(Age) + volume
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +0.0001231 (base 128.8283478 code 0 | lane 128.8284708 code 0) count ~ Trt + Age + mo(Exp) + s(Age) + volu
lane higher (worse)  frmtmb.sample/test-brms-suite-methods.R    +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volu
lane higher (worse)  frmtmb/test-brms-suite-data-helpers.R      +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volu
lane higher (worse)  frmtmb/test-brms-suite-emmeans.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volu
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volu
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volu
lane higher (worse)  frmtmb/test-brms-suite-methods.R           +4.469e-05 (base 109.3400595 code 1 | lane 109.3401042 code 1) count ~ Trt * Age + mo(Exp) + s(Age) + volu
lane lower (better)  frmtmb/test-optrap.R                       -118.1 (base 2913.226547 code 0 | lane 2795.149471 code 0) y ~ t + (1 + t | s) gaussian
lane lower (better)  frmtmb/test-se-check.R                     -55.11 (base 443.4558358 code 1 | lane 388.3498424 code 1) ls ~ mo(inc) + age gaussian
lane lower (better)  frmtmb/test-optrap.R                       -18.68 (base 1543.327912 code 0 | lane 1524.65226 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-grad-verdict.R                 -18.37 (base 739.7082038 code 1 | lane 721.3384413 code 1) y ~ x + z gaussian
lane lower (better)  frmtmb/test-grad-verdict.R                 -15.52 (base 26783.24826 code 1 | lane 26767.72855 code 1) y ~ x + z gaussian
lane lower (better)  frmtmb/test-optrap.R                       -14.99 (base 1532.567176 code 0 | lane 1517.573737 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -10.68 (base 1532.355266 code 0 | lane 1521.676117 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.8245 (base -553.9283063 code 1 | lane -554.7528471 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base -547.0088252 code 0 | lane -547.674324 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base -1237.784353 code 1 | lane -1238.449852 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.6655 (base 834.5422026 code 0 | lane 833.8767318 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.4487 (base 324.463624 code 0 | lane 324.0149446 code 0) ls ~ mo(income) * age gaussian
lane lower (better)  frmtmb/test-emmeans.R                      -0.08293 (base 1.321783367 code 0 | lane 1.23885053 code 0) y ~ f + gp(x) gaussian
lane lower (better)  frmtmb/test-optrap.R                       -0.04888 (base 331.5744044 code 0 | lane 331.5255205 code 0) ls ~ mo(income) * age gaussian
lane lower (better)  frmtmb/test-fuzz.R                         -0.002016 (base 176.8571893 code 0 | lane 176.8551738 code 0) y | weights(w) ~ 1 + x + mo(xo) * z + (1 | g
lane lower (better)  frmtmb/test-fuzz.R                         -0.002016 (base 176.8571893 code 0 | lane 176.8551738 code 0) y | weights(w) ~ 1 + x + mo(xo) * z + (1 | g
lane lower (better)  frmtmb/test-fuzz.R                         -0.002016 (base 176.8571893 code 0 | lane 176.8551738 code 0) y | weights(w) ~ 1 + x + mo(xo) * z + (1 | g
lane lower (better)  frmtmb/test-fuzz.R                         -0.002016 (base 176.8571893 code 0 | lane 176.8551738 code 0) y | weights(w) ~ 1 + x + mo(xo) * z + (1 | g
lane lower (better)  frmtmb.sample/test-brms-pins.R             -0.0009362 (base 463.5842302 code 1 | lane 463.583294 code 1) y ~ x + (1 + x + z + w | g) gaussian
lane lower (better)  frmtmb/test-v14.R                          -0.0006332 (base 1001.633972 code 1 | lane 1001.633339 code 1) y ~ x zero_inflated_asym_laplace
lane lower (better)  frmtmb/test-importance.R                   -8.707e-06 (base 193.0767336 code 0 | lane 193.0767249 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-importance.R                   -8.707e-06 (base 193.0767336 code 0 | lane 193.0767249 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-importance.R                   -8.707e-06 (base 193.0767336 code 0 | lane 193.0767249 code 0) y ~ x + (1 | g) gaussian
lane lower (better)  frmtmb/test-review-v29.R                   -2.035e-06 (base 67.73199382 code 1 | lane 67.73199179 code 1) mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian()

worse by more than 1e-6 relative: 4 
```
<!-- END GENERATED -->

Reading it:

- 8688 fits reach the same optimum; 123 error on both builds (the
  fuzz tier's expected refusals); 183 have no optimum to compare
  (unfitted objects, quadrature refusals); 781 and 782 are unmatched
  (`dev/optrap-unmatched.R`): 1554 of the 1563 in five files whose fits
  draw data without a seed (`test-ce-bands.R` 680, `test-ce-levels.R`
  450, `test-ce-options.R` 184, `test-simulate-newdata.R` 160,
  `test-portability.R` 80), and 9 in `test-verbose.R`,
  `test-sugar.R` and the two test files this lane rewrote after the
  base run (`test-se-tier.R` 2, `test-grad-verdict.R` 1).
- 35 fits changed, all of `gaussian`, `student` or
  `zero_inflated_asym_laplace` (the units). 24 are better: 16 where
  both builds report convergence, 2 that converge only on the lane
  (the trap fits of `test-optrap.R`), 6 where neither does (fits that
  stop at an iteration cap by construction: the far-start fits of
  `test-grad-verdict.R` and the `iter.max = 3` fit of
  `test-se-check.R`; and `test-v14.R`, `test-review-v29.R`,
  frmtmb.sample's `test-brms-pins.R`).
- 11 are worse, all fixture 1 of the ported suite and its refits in
  `test-brms-suite-methods.R` (student, `mo()`, `s()`, `arma`,
  `sigma ~ Trt`, unit 8). Ten stop at "false convergence (8)" on both
  builds; four of those are refits on changed data whose objective
  goes negative on 40 counts (-247.80 against -158.36, -199.57
  against -158.97 twice, -142.98 against -110.93), a descent neither
  build completes. One converges on both builds and is 1.2e-4 (9.6e-7
  relative) higher on the lane. These are not improvements, and they
  are the ones left: the fixture's sd runs to -723 to -1253 on the log
  scale on both builds, and where nlminb stops on that flat edge
  follows its path.

## R CMD check

`dev/optrap-check.sh` (core, built and checked in `dev/optrap-check/`,
`R_LIBS` the lane's library, rellib-r7, the user library; pandoc and
TinyTeX on PATH; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`): Status: 1
NOTE (`checking HTML version of manual`, the environmental V8 note);
tests OK, `[ FAIL 0 | WARN 0 | SKIP 382 | PASS 11897 ]`. No extension
changed, so no extension was checked.

## Left open

- **The ledger row `brmsfit-methods:326` no longer holds** on either
  BLAS (fixture 1's sigma coefficients get negative variances at the
  lane's stop). It is recorded as pass; the regeneration at
  consolidation should move it back (`dev/brmsport-verdicts.tsv` not
  edited, as the round asks).
- **The one-level group with a bound**: 5 of 40 seeds zigzag to the
  evaluation limit on the lane, 0 on rellib-r7 (above). Not fixed.
- **Fixture 1's degenerate refits** stop elsewhere (above). Not fixed.
- **fit_loss2 and F11 seed 8** stay short (0.032, 0.1418), as lane
  optima filed; their causes are not shared with items 1 and 2.
- **test-id-kron.R** still ends at one of two points by core type; the
  test no longer depends on which (item 3). The SE check's verdict at
  the two points differs (boundary on P, nothing lost on E), which is
  lane seflags' territory.
- **The objective's constant** (`n log u`) stays in the objective
  nlminb sees, so the invariance is 1e-7 in the log-likelihood rather
  than exact; removing it made relative convergence stricter and broke
  converged `mo()` fits (above).
- **A Newton polish** of fixed-effect fits (TMB's exact Hessian is
  available without random effects) would make the precision tests'
  external agreement exact; not built (out of the lane's scope, and it
  would touch every fixed-effect fit).

## Overlaps with the other lanes

- `R/se-check.R` is not edited, but `fit$par_units` now carries the
  response's unit for the ten families, so the SE check's Hessian
  scaling (`fit_outer_hessian()`, `se_tier3_could_act()`, the edge
  screen) and `autoscale_sdreport()` read those units on such fits;
  lane seflags should rebase on this.
- Tests edited that lane seflags may also edit: `test-se-tier.R` (the
  RB2 and RQ1 constructions now start the sd at exp(-20); the healthy
  sleepstudy fit uses Reaction / 64), `test-se-check.R` (the healthy
  `mo()` fit uses ls / 16; m3's data are 0.69.0's again since punch
  round 1). `test-grad-verdict.R` is 0.69.0's again.
- `R/fit.R`: `mo_search()` (lane moparity works on `mo()`; the search
  changed here), `make_start()`, `check_convergence()`,
  `grad_verdict()`, `fit_assembled()`, `optimize_obj()`,
  `optimizer_from_best()`, `obj_settle()`, `quad_fit()`,
  `fit_recovery_starts()`.
- `NEWS.md`: a development section with two headings.

## Version

A minor bump (breaking: starting values and optimizer units of the
ten location-scale families change, and with them fits that stopped
short). No floor moves: no extension changed, and frmtmb.sample reads
nothing this lane changed.
