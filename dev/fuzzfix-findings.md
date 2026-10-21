# Lane fuzzfix findings (round opened 2026-10-20, base 4854e425)

Status: FIX ROUND 2 DONE (see "Fix round 2"; fix round 1 below).
The round-0 record below is kept as it was;
where fix round 1 changed a decision, the sections "Fix round 1",
"To file" and "Fix round 1: suites" say so and win.

Fix round 1 in one paragraph: nested, crossed and mu-sigma quadrature
are refused (lme4's rule), single-block quadrature is bitwise as on
rellib-r10, and the exact iterated route is saved as
`dev/fuzzfix-nested-quad.patch` with a design note. The weights
refusal needs an intercept in the dispersion formula, ignores priors
that do not provably bound it, runs on frm_sample()'s formula route
as a posterior check, and now covers asym_laplace (Laplace, any
weights) and exgaussian. `R CMD check --as-cran`: core `Status: 1
NOTE` (V8 math rendering, `fr1-check-frmtmb.txt`), frmtmb.sample
`Status: OK` (`fr1-check-frmtmb.sample.txt`). The first core check of
this round stopped inside the tests with no status line and no error
in its log (`fr1-check-frmtmb-run1.txt`; cause not found, it ran
beside 4 session processes of this lane and other lanes' runs); the
rerun alone gave the status above. Bump: frmtmb minor; frmtmb.sample
minor, with a floor on this core (it calls the new export
`check_weights_unbounded()`).

The round-0 status line read: both fixes in core, all eight suites
green, `R CMD check --as-cran` of core `Status: 1 NOTE`. The review
(`dev/reviews/2026-10-20-fuzzfix.md`) returned MERGEABLE AFTER FIXES.

This file is the record if the report does not arrive.

Two specs that the gated fuzz plan draws once lane cdfs merges. Both
fail the same way on rellib-r10 (0.72.0). Neither involves a CDF.
Logs are in `dev/fuzzfix-out/` (local, gitignored); every number below
is pasted from them.

## Reproductions on rellib-r10

`dev/fuzzfix-repro-spec.R base cum|gam` rebuilds each spec from the
generator in this worktree's `helper-fuzz.R` (`fuzz_data()` and
`fuzz_design()` are the same as cdfs's for an aterm of `weights`),
saves the data to `dev/fuzzfix-out/data-{cum,gam}.rds`, and fits it.

1. cumulative (seed 20322452), log `repro-cum-r10.txt`:

       frm(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) + cumulative(),
           data = d, REML = FALSE, quadrature = TRUE, ...)
       SUMMARY ERROR: no 'dimnames' attribute for array
       In sdreport_xtra(obj, sdr, ..., what = "raneffvector") :
         sdreport: Cannot mix inner and outer integration methods (yet)

2. Gamma (seed 20502220), log `repro-gam-r10.txt`:

       frm(bf(y | weights(w) ~ 1 + x + mo(xo) * z + ar1(tt + 0 | g),
              shape ~ 1 + z) + Gamma(link = "log"), ...,
           control = frmtmb_control(profile = TRUE))
       The optimizer failed on this model (Gamma, ML, profile, nlminb):
       NA/NaN gradient evaluation.

## Defect 1: quadrature over nested blocks was half Laplace

### Root cause

`fit_assembled()` passed `integrate = list(b = GK(...))`, one TMBad
`marginal_gk` pass over every `b`. TMB (`MakeADFun`, TMB 1.9.25)
integrates the effects the transform can split into one-dimensional
factors, removes them, and keeps the rest as `random`, which the
Laplace approximation then integrates. Over `(1 | ga/gb)` it removed
the 12 `ga` effects and kept the 24 `ga:gb` effects
(`dev/fuzzfix-nested-exact.R`, `nested-exact-{cum,bern}.txt`):

    block 1 | ga group ga levels 12 b_idx 1 12
    block 1 | ga:gb group ga:gb levels 24 b_idx 13 36
    b total: 36  left random in quadrature tape: 24

So the fit was a mixed one, and

- its objective was the Laplace one (cum: quadrature 126.188224721,
  Laplace 126.18822478; bern: 57.9445979746 against 57.9445979762),
  while the iterated Gauss-Hermite integral written in R over frmtmb's
  own joint objective (80 nodes per level) is 126.086168724 and
  57.5587442351: errors of 0.1021 and 0.3859, the Laplace errors;
- `RTMB::sdreport()` calls `sdreport_xtra()`, which warns "Cannot mix
  inner and outer integration methods (yet)" and returns nothing, so
  `sdr$par.random` is empty while `sdr$env$random` lists all 36
  effects; `TMB:::as.list.sdreport()` then indexes an empty summary
  (`summary(x, select = "random")[, what]`) and fails with "no
  'dimnames' attribute for array". `par_est_se()` in `summary()` is
  the first caller; `vcov()` does not read the random part and works.

Neither the family nor the weights matter. `dev/fuzzfix-quad-survey.R`
on rellib-r10 (`quad-survey-r10.txt`), every multi-block arm keeps
random effects and fails summary():

    ARM cum_fuzz       nb  36 left_random  24 ... diff 5.92e-08 summary ERROR
    ARM cum_fuzz_now   nb  36 left_random  24 ... diff 1.08e-08 summary ERROR
    ARM cum_fuzz_ri    nb  12 left_random   0 ... diff 0.0168 summary ok
    ARM bern_single    nb  24 left_random   0 ... diff 0.114 summary ok
    ARM bern_nested    nb  72 left_random  48 ... diff 4.25e-08 summary ERROR
    ARM pois_nested    nb  72 left_random  48 ... diff 0.0288 summary ERROR
    ARM pois_nested3   nb 168 left_random 144 ... diff -0.0206 summary ERROR
    ARM pois_crossed   nb  32 left_random   8 ... diff 0.0851 summary ERROR
    ARM cum_nested     nb  72 left_random  48 ... diff 1.3e-08 summary ERROR
    ARM sratio_nested  nb  72 left_random  48 ... diff 7.49e-09 summary ERROR
    ARM acat_nested_w  nb  72 left_random  48 ... diff 3.12e-08 summary ERROR
    ARM gaus_mu_sigma  nb  48 left_random  24 ... diff 0.268 summary ERROR

(`diff` is logLik quadrature minus logLik Laplace; `cum_fuzz_now` is
the fuzz spec without weights, `cum_fuzz_ri` with `(1 | ga)` only.)
The other modes are unaffected: `dev/fuzzfix-cum-modes.R`
(`cum-modes-r10.txt`) runs summary() on the fuzz data under ML,
`profile = TRUE` and `autoscale = TRUE` for cumulative, sratio, cratio
and acat: "summary ok" in all 12 arms on rellib-r10. The compat row claimed "(1 |
ga/gb) becomes an iterated one-dimensional integral"; it was not.

### Fix (R/fit.R)

- `quad_passes()`: one `marginal_gk` pass per block, the block with the
  most effects first, each pass's `which` indexing the `b` entries
  earlier passes left. The fine block's effects are then
  one-dimensional integrals given the coarse effect, and the coarse
  block's integrand is a product over its own levels again. Measured
  before wiring it in (`nested-exact-*.txt`, "two-pass"): 0 effects
  left, objective minus the exact integral -3.73e-09 (bern) and
  5.728e-09 (cum), gradient finite.
- `quad_crossed_stop()`: crossed blocks cannot be split this way; the
  passes leave effects, and `quad_fit()` refuses when the tape keeps
  any `random` (a property of the tape's graph, so it holds for every
  calibration point).
- `quad_check_depth()`: refuses a row that reads more than two blocks,
  before any tape. `dev/fuzzfix-quad-depth.R` (Poisson, 32 rows, seed
  7, one process per depth):

      DEPTH 1 rows 32 nb 4 time 0.5s peakMB 355 logLik -48.352747
      DEPTH 2 rows 32 nb 12 time 3.0s peakMB 843 logLik -48.333403
      DEPTH 3 rows 32 nb 28 time 377.6s peakMB 49119 logLik -48.333403

  The first survey run on the lane build, without this guard, died at
  `std::bad_alloc` in `TransformADFunObjectTemplate` on the 168-effect
  three-level arm. The count is per row, so a by-split `gr(g, by = f)`
  with three blocks is still accepted (`quad-ratio-lane.txt`: "by-split
  blocks: 3, left random 0, logLik -41.708837", the same as r10).

After (`quad-survey-lane2.txt`):

    ARM cum_fuzz       nb  36 left_random   0 logLik_quad -126.08035126 ... diff 0.108 summary ok vcov ok (27.3s)
    ARM cum_fuzz_now   nb  36 left_random   0 ... diff 0.0306 summary ok (26.1s)
    ARM cum_fuzz_ri    nb  12 left_random   0 ... diff 0.0168 summary ok (0.3s)
    ARM bern_single    nb  24 left_random   0 ... diff 0.114 summary ok (0.3s)
    ARM bern_nested    nb  72 left_random   0 ... diff 0.717 summary ok (31.9s)
    ARM pois_nested    nb  72 left_random   0 ... diff 0.283 summary ok (43.6s)
    ARM pois_nested3   fit-error: ... at most two deep, and here a row reads 3 blocks
    ARM pois_crossed   fit-error: ... Here 7 random e[ffects are left]
    ARM cum_nested     nb  72 left_random   0 ... diff 0.617 summary ok (53.8s)
    ARM sratio_nested  nb  72 left_random   0 ... diff 0.4 summary ok (68.4s)
    ARM acat_nested_w  nb  72 left_random   0 ... diff 0.697 summary ok (97.4s)
    ARM gaus_mu_sigma  nb  48 left_random   0 ... diff 0.218 summary ok (26.8s)

On the fuzz spec the lane's quadrature objective minus the exact
integral is 4.794e-09 (`nested-exact-cum-lane.txt`).

The cost is time. The fuzz cumulative spec fits in 27.3 s against 1.4
s for the r10 mixed fit; `dev/fuzzfix-quad-time.R`
(`quad-time-lane.txt`): tape 3.95 s, optimize 5.97 s, fit 19.95 s,
summary() 15.52 s (RTMB's `sdreport_xtra()` tapes the integral twice
more).

## Defect 2: a weighted likelihood with no maximum

### Root cause

Not a step outside the support, and not the profile route.
`ar1(tt + 0 | g)` in the fuzz design has one row per (g, tt) cell
(`gam-weights-r10.txt`: "rows: 72 (g, tt) cells: 72"), so every row
has a random effect of its own in `mu`, and the weights
(`runif(0.5, 2)`) sum to n + 15.9748. Take each effect so that the
row's mean equals its observation and let the dispersion go to its
limit (Gamma `shape` to infinity). The density at its peak grows like
1 / width, and the curvature of the row's own effect like 1 / width^2,
so the Laplace (and the exact) marginal grows like
width^(n - sum(w)): without bound when sum(w) > n. The same holds
under REML and `profile = TRUE` (the beta block of the joint Hessian
has a width-free Schur complement when every row owns an effect).

Minimized one component at a time (`dev/fuzzfix-gam-arms.R`,
`gam-arms-r10.txt`):

    ARM full         ERROR ... NA/NaN gradient evaluation
    ARM no_weights   ok logLik -159.557357 | warnings: 0
    ARM no_profile   ok logLik 517083.644502 | warnings: 2 (false convergence, large gradient)
    ARM no_shape     ok logLik 59.448552 | warnings: 2
    ARM re_ri        ok logLik -196.929849 | warnings: 0
    ARM re_none      ok logLik -196.027523 | warnings: 0
    ARM mo_main      ok logLik 52.866061 | warnings: 2
    ARM mo_none      ok logLik 12.217626 | warnings: 2
    ARM min_pro      ok logLik -161.667767 | warnings: 0

So weights plus the observation-level ar1 are the cause; mo(), the
shape formula and profile only change how the runaway ends. The
weights alone (`dev/fuzzfix-gam-weights.R`, `gam-weights-r10.txt`,
`ar1` and `(1 | obs)` arms, ML and profile):

    ARM ar1 w=w profile=FALSE        logLik 139844810784.0173 betad 53.91
    ARM obs w=w profile=FALSE        logLik 7.35e50 betad 144.6
    ARM ar1 w=w profile=TRUE         logLik 5088.6100 betad 36.74
    ARM ar1 w=w/mean profile=FALSE   logLik -171.8420 betad 15.72  (1 warning)
    ARM ar1 w=0.9 profile=FALSE      logLik -158.9845 betad 0.748
    ARM ar1 w=1 profile=FALSE        logLik -174.2091 betad 1.066
    ARM ar1 w=1.1 profile=FALSE      logLik 8.2e62 betad 171.9
    ARM obs w=1.1 profile=TRUE       ERROR ... NA/NaN gradient evaluation
    ARM gaussian obs w=1.1           logLik 148.1991 betad -35.91
    ARM gaussian obs w=0.9           logLik -89.2837 betad -0.0411

`w/mean` (sum(w) = n) still drifts to a log shape of 15.7: the leading
term vanishes there and the lower-order terms decide; it is not
refused (not proved unbounded).

### Fix (R/fit.R, `check_weights_unbounded()`, called from frm())

Refuses before fitting, after the `dry_run = "frame"` return (the
fuzz brms oracle reads the frame) and not for `dry_run = "objective"`
(frmtmb.sample's priors can bound the growth). It fires only where the
argument above is a proof, and fails open otherwise: one response; a
family in `unbounded_spike_dpar` with that dpar free (not constant);
no `se()`, `cens()`, residual correlation, `mi()` or prior; every row
with a positive weight owns a column (one non-zero row) of the `mu` Z
restricted to full-rank covariance structures; sum(w - 1) over those
rows above 1e-8 n.

The family table is measured, not assumed (`dev/fuzzfix-spike-families.R`,
seed 20261020, n = 80, `(1 | obs)`, constant weights). On rellib-r10
(`spike-families-r10.txt`) every table family runs away at w = 1.2 and
fits at 0.9:

    FAM gaussian w=1.2 logLik 486.334 | w=0.9 -69.7294
    FAM student  w=1.2 256.895        | -68.7394
    FAM lognormal 487.574 | -67.7058     FAM Gamma 205.535 | -68.3301
    FAM weibull 170.012 | -70.241        FAM Beta 1696.28 | 49.9467
    FAM inverse.gaussian 480.766 | -67.2527
    FAM skew_normal 213.179 | -61.8233   FAM shifted_lognormal 484.174 | -67.0479
    FAM huber 3704.62 | -64.7878         FAM von_mises -3.32338 | -68.047
    FAM frechet 4.649e218 | -70.7233     FAM gen_extreme_value 5.233e192 | -47.2256
    FAM poisson -136.439 | -102.382      FAM exponential -123.444 | -92.5829

On the lane (`spike-families-lane.txt`) all 13 table families are
refused at 1.2 and give the same logLik at 0.9; poisson and
exponential (no spike) are unchanged. asym_laplace is left out: it
runs away at 0.9 too (logLik 1584.11, betad -23.75), because its log
density is linear on each side of the mode and the Laplace Hessian of
an observation-level effect is then the prior's alone. That is a
different defect, filed below, not fixed.

## brms 2.23.0 on the same data (`dev/fuzzfix-brms.R`, seed 20261020)

4 chains x 2000 iterations, brms defaults.

cumulative + weights + (1 | ga/gb), fuzz data (`brms-cum.txt`):

| parameter | frmtmb quadrature (SE) | frmtmb Laplace (SE) | brms mean (SD) |
|---|---|---|---|
| x | 0.71366 (0.21912) | 0.71010 (0.21830) | 0.72 (0.23) |
| sd(ga:gb) | 0.8184 | 0.7807 | 0.84 (0.42) |
| sd(ga) | 0.00032 | 0.00032 | 0.56 (0.39) |

The ML sd(ga) is at the boundary; the brms posterior mean of a sd
with mass near 0 is not comparable to it. x and sd(ga:gb) agree within
0.1 posterior SD. logLik: quadrature -126.0804, Laplace -126.1882.

Gamma + weights + latent AR(1) + mo(xo) * z + shape ~ 1 + z
(`brms-gam.txt`). brms spells frmtmb's `ar1(tt + 0 | g)` as
`ar(time = tt, gr = g)`, whose latent residuals are the same
observation-level AR(1) effects.

- With the fuzz weights (sum(w) - n = 15.97) frmtmb refuses. brms
  does not converge: Rhat 3.0 to 4.0 on every parameter, bulk ESS 4,
  `shape_Intercept` 23.63 (a shape of about 1.8e10). Its default prior
  on the dispersion intercept does not bound the growth either.
- With the weights times 0.7 (sum(w) - n = -10.42), frmtmb's ML fit
  against the brms posterior:

| parameter | frmtmb (SE) | brms mean (SD) |
|---|---|---|
| Intercept | 0.77 (0.17) | 0.80 (0.19) |
| shape_Intercept | 0.74 (0.17) | 0.77 (0.28) |
| x | 0.53 (0.09) | 0.52 (0.11) |
| z | 0.50 (0.15) | 0.41 (0.17) |
| shape_z | -0.47 (0.21) | -0.44 (0.26) |
| moxo | 0.36 (0.07) | 0.35 (0.09) |
| moxo:z | -0.09 (0.07) | -0.04 (0.10) |

  Every coefficient is within 0.6 posterior SD. The AR(1) sd is at 0
  in the ML fit (boundary, NaN SE) and 0.29 (0.18) in brms.

## Regression tests (plain NOT_CRAN), seen to fail on rellib-r10

`dev/fuzzfix-testfile.R <lib> frmtmb <file>` (logs `tf-*-{base,lane}.txt`):

    RESULT base test-quadrature-nested.R pass=2 fail=4 err=2 skip=0 warn=6
    RESULT lane test-quadrature-nested.R pass=10 fail=0 err=0 skip=0 warn=0
    RESULT base test-weights-unbounded.R pass=13 fail=2 err=1 skip=0 warn=4
    RESULT lane test-weights-unbounded.R pass=17 fail=0 err=0 skip=0 warn=0

The base failures are behavioral: `fq$obj$env$random` "Actual length:
24"; `err_q/err_l` "0.353 >= 0.001"; summary() "no 'dimnames'
attribute for array" (both nested tests); the crossed fit "to throw a
error" did not; the Gamma ar1 ML and profile fits did not throw; the
REML fit died at "NA/NaN gradient evaluation". The absent cases (sum
(w) < n, six rows per level, unit weights, poisson, dry_run = "frame")
pass on both builds.

Under OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`, then `dev/fuzzfix-testfile.R`; `tf-*-ob0326.txt`):

    RESULT ob0326 test-quadrature-nested.R pass=10 fail=0 err=0 skip=0 warn=0
    RESULT ob0326 test-weights-unbounded.R pass=17 fail=0 err=0 skip=0 warn=0

0.3.32 was not run: neither test sits near a rounding-decided branch
(the ratio is 1.9e-6 against a bound of 1e-3; the weight sums are n
- 0.3 n and n + 0.25 n against n).

`dev/fuzzfix-quad-ratio.R` prints the numbers test 1 asserts on:
lane err_q 1.61e-06, err_l 0.853, ratio 1.89e-06; r10 err_q 0.278,
err_l 0.789, ratio 0.353.

## Suites, one file per process (NOT_CRAN)

`dev/release/run-par.sh fuzzfix-notcran dev/fuzzfix-out/jobs-all.txt 4`
with `FRMTMB_RELLIB` = the lane library (core from this worktree, the
seven extensions from rellib-r10; every `lib:` line names the lane's
frmtmb). `FUZZFIX-NOTCRAN ran 412 of 412`, summed from
`dev/release/fuzzfix-notcran.log` (`notcran-summary.txt`):

    frmtmb files=255 pass=18279 fail=0 err=0 skip=191 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=32 pass=1838 fail=0 err=0 skip=5 warn=0
    frmtmb.latent files=11 pass=368 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=15 pass=430 fail=0 err=0 skip=13 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=61 pass=2849 fail=0 err=0 skip=8 warn=0
    frmtmb.spline files=16 pass=594 fail=0 err=0 skip=1 warn=0

One session per suite, three orders (`dev/release/run-sessions.sh
fuzzfix-sess 4`, NOT_CRAN, `dev/release/fuzzfix-sess.log`, "ran 24 of
24"); every LEAK line carries "[with load of ...]" (untagged: 0):

    SESSION frmtmb forward files=245/255 pass=18279 fail=0 err=0 skip=191 warn=0
    SESSION frmtmb reverse files=245/255 pass=18377 fail=0 err=0 skip=191 warn=0
    SESSION frmtmb preload files=245/255 pass=18473 fail=0 err=0 skip=192 warn=0
    SESSION frmtmb.coupling forward|reverse|preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam forward|reverse|preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.latent forward|reverse|preload files=11/11 pass=368 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn forward|reverse|preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.ode forward|reverse|preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample forward|reverse|preload files=56/61 pass=2849 fail=0 err=0 skip=8 warn=0
    SESSION frmtmb.spline forward|reverse|preload files=16/16 pass=594 fail=0 err=0 skip=1 warn=0

(the three orders of each extension print identical counts; the log
has one line each.)

Gated, core only: the 29 core files with a skip above were rerun with
`FRMTMB_BRMS_FIT_TESTS=true` (`dev/release/fuzzfix-gated.log`, "ran
29 of 29", summed: pass=2689 fail=0 err=0 skip=14 warn=0); the two that still skipped on their own
gates were rerun with `FRMTMB_FUZZ=true` and
`FRMTMB_DRMTMB_FIT_TESTS=true` (`fuzzfix-gated2.log`):

    RESULT test-fuzz.R pass=2 fail=0 err=0 skip=0 warn=0
    RESULT test-drmtmb-agreement.R pass=131 fail=0 err=0 skip=0 warn=0

That fuzz run is this worktree's plan (base families table plus the
`fuzz_design()` change), not cdfs's; see the next section for the
merged plan's subset.

## What the merged gated fuzz plan needs

`dev/fuzzfix-fuzz-subset.R` reads cdfs's `helper-fuzz.R` in place
(read-only) and draws its plan, `fuzz_plan(seed = 20260901L, size =
300L)`: 300 rows, 6 with ar1 + weights + a continuous family, 6 with
quadrature + nested. Those 12 run through `fuzz_run()`.

- The six ar1 + weights specs are models with no maximum likelihood
  under the generator's weights (`runif(0.5, 2)`, mean 1.25). With
  cdfs's generator, on rellib-r10 they gave 1 REAL-NEW (the Gamma
  spec) and 5 UNCONVERGED findings, among them a logLik identity off
  by 2.19 and NaN objectives (`fuzz-subset-r10-cdfsgen.txt`); on the
  lane build all 6 are REAL-NEW fit_error, the refusal
  (`fuzz-subset-lane-cdfsgen.txt`).
- So `helper-fuzz.R`'s `fuzz_design()` (not its known lists) scales
  the weights by 0.6 when `re == "ar1"`, after the draw, so every
  other column keeps its random stream. With that change (merged into
  cdfs's helper by the script) the 12 specs give
  (`fuzz-subset-lane.txt`, 311 s): 0 REAL-NEW from either fix, 1
  UNCONVERGED (student + ar1 + smooth, vcov_psd, "singular
  convergence (7)", which the tier does not fail on), and 2 REAL-NEW
  from Beta + cens + nested quadrature (seed 20535438): "cens()/trunc()
  need a family with a CDF", which is cdfs's to add. That spec meets
  nested quadrature with a censored Beta for the first time at the
  merge; run it on the merged build.
- No known-list entry is needed for either spec. The merge has to keep
  this lane's `fuzz_design()` hunk beside cdfs's `fuzz_data()` and
  family-table hunks; they do not overlap.

## Not done, and why

- No brmsport ledger rerun: neither path is a brms grammar question,
  and the gated core files (the generated brms-suite files among them)
  pass. Consolidation regenerates the ledger.
- OpenBLAS 0.3.32 not run (see the tests section).
- The merged-build run of Beta + cens + nested quadrature (seed
  20535438) needs cdfs's Beta CDF; it is left for consolidation.
- Extensions not reinstalled from this worktree: none changed; the
  rellib-r10 copies served them, each `lib:` line naming this lane's
  core.

## Defects found, not fixed

- A quadrature fit with a variance at zero gets a finite Wald SE for
  the sd and an interval below zero, where the Laplace fit gives NaN
  and the boundary message (`dev/fuzzfix-quad-boundary.R`, seed 3):
  one block on rellib-r10 and on the lane alike, "sd 4.773e-05
  Est.Error 0.1355 Q2.5 -0.2654"; the nested fuzz spec on the lane,
  "sd 0.0003206 Est.Error 0.9595 Q2.5 -1.88". The boundary and
  no-curvature rules (`sdr_nan_nonpos()`) do not run on a quadrature
  objective. Older than this lane; nested quadrature only now reaches
  summary() to show it.

- asym_laplace with an observation-level effect has no maximum of the
  Laplace objective at any weight (w = 0.9: logLik 1584.11 with sigma
  at exp(-23.75)); the exact marginal is bounded there. Refused in fix
  round 1.
- mvbf() with two gaussian responses, each with `(1 | obs)` and
  weights summing above n, still runs away without a refusal (logLik
  790.247 on both builds, review); the check fails open for several
  responses.
- `sum(w) = n` (weights normalized to mean 1) with an observation-level
  effect drifts to a log shape of 15.7 on the fuzz data; not proved
  unbounded, not refused.

## Fix round 1 (review `dev/reviews/2026-10-20-fuzzfix.md`)

Decisions are the coordinator's (the user may override at
consolidation). Logs are `dev/fuzzfix-out/fr1-*`.

### Quadrature: refuse, do not iterate (findings 2, 3, 4)

- `quad_passes()`, `quad_check_depth()` and `quad_crossed_stop()` are
  gone. The single `integrate = list(b = GK(...))` of rellib-r10 is
  back, and `quad_check_blocks()` (R/fit.R) refuses before any tape
  when a row reads more than one block, read through the `Z` of every
  linear predictor of the response. Nested, crossed and `(1 | g)` in
  `mu` and `sigma` are refused; a by-split `gr(g, by = f)` (three
  blocks, one per row) and `mvbf()` with one block per response pass.
  The message names the blocks and points to `quadrature = FALSE`. A
  backstop in `quad_fit()` still stops any tape that keeps `random`.
- `quad_breakdown_message()` lost its nested branch (now unreachable).
- Single-block quadrature is unchanged. `dev/fuzzfix-rev-quad-cost.R`
  under `dev/fuzzfix-rev-capped.ps1` (cap 3000 MB), two processes per
  arm, interleaved (`fr1-cost-single-*`):

      single 288 rows: quad base 0.33s 0.34s | lane 0.33s 0.38s | logLik -189.30039101 in all four
      single 576 rows: quad base 0.57s 0.57s | lane 0.56s 0.55s | logLik -336.02942142 in all four
      peak MB: 288 rows base 396 401 | lane 399 392; 576 rows base 402 438 | lane 447 431

  The nested 72-row arm on the lane stops in the refusal after its
  0.59 s Laplace step (`fr1-cost-nested-12-3-lane.txt`).
- Fuzz tier: a nested quadrature spec now meets a deliberate refusal.
  It is registered in `FUZZ_KNOWN_REFUSAL` as `quad-one-block-per-row`,
  matched on `mode=quadrature`, `re=nested`, `fit_error` and the fixed
  string "quadrature = TRUE integrates one scalar random effect per".
  Chosen over not drawing it because `fuzz_ok()` feeds the pairwise
  cover: excluding nested quadrature would redraw the whole 300-spec
  plan, including every spec cdfs measured, while the entry leaves the
  plan as it is and still exercises the refusal on the 6 nested specs.
  The stale nested sentence of the `quad-breakdown` entry's `why` was
  removed.
- NEWS: one BREAKING bullet for the refusal; the round-0 bug-fix
  bullet and the "more than two blocks" bullet are gone. The
  `quadrature` docs and the compat row say one block per row.

### Weights refusal (findings 1, 5, 6, 10)

Moved to the new `R/weights-unbounded.R`.

1. Finding 1: `unbounded_spans_constant()` requires every spike dpar's
   design to have a column and to span the constant (`qr.resid()` of
   `rep(1, n)` below 1e-8 sqrt(n)), and not to be constant or
   nonlinear. The reviewer's `dev/fuzzfix-rev-wu.R` on the lane
   (`fr1-rev-wu-lane.txt`) against its rellib-r10 log
   (`dev/fuzzfix-rev-out/wu-base.txt`):

       FP1_sigma_0+xs     base ok logLik -137.731 betad 0.5282 | lane ok logLik -137.731 betad 0.5282
       FP2_sigma_offset   base ok logLik -143.083              | lane ok logLik -143.083
       FP2b_sigma_const   base ok logLik -143.083              | lane ok logLik -143.083

   Round 0 refused FP1 and FP2 (`dev/fuzzfix-rev-out/wu-lane.txt`).
   The test block "a dispersion without an intercept cannot reach its
   limit" fits both and asserts FP2's logLik equals the constant-sigma
   fit's. Every other arm of that script matches the base log line for
   line, except the refusals listed here and in item 2.
2. Finding 5: a prior exempts the model only when
   `prior_spike_verdict()` proves it bounds the growth; a prior on
   anything but a spike dpar's intercept is ignored. On the natural
   scale toward 0 the density goes like x^order (0 for normal,
   student_t, logistic, exponential; shape - 1 for gamma; infinite for
   lognormal, inv_gamma or a lower bound above 0), and MAP is bounded
   when order > sum(w) - n. Toward infinity a tail faster than any
   power bounds it; student_t and inv_gamma there are "unknown" and
   fail open. On the log link, normal bounds and student_t does not.
   Lane (`fr1-rev-wu-lane.txt`):

       prior_sigma_t      student_t(3, 0, 2.5) on sigma  base ok logLik 533.348 | lane REFUSED
       prior_sigma_gamma2 gamma(2, 1) on sigma           base ok logLik 500.951 | lane REFUSED
       prior_b_only       normal(0, 5) on b              base ok logLik 567.274 | lane REFUSED
       prior_shape_gamma  gamma(0.01, 0.01) on shape     base ok logLik -92.5793 | lane ok -92.5793

   The message (`dev/fuzzfix-fr1-msg.R`, `fr1-msg.txt`) prints the
   term as written (`ar1(tt + 0 | g)`, `(1 | obs)`), drops "Rescale the
   weights", explains that a weight above 1 counts a row more than
   once, and states when the model is bounded: weights summing to less
   than the number of rows, no observation-level term, or a prior on
   sigma that vanishes at 0 faster than sigma^(sum(w) - n)
   (lognormal() and inv_gamma() do; student_t(), normal() and
   exponential() do not).
3. Finding 6: `check_weights_unbounded()` is exported to the
   extensions (sampling-API list, `man/frmtmb-sampling-api.Rd`), and
   `frm_sample()`'s formula route calls it after it resolves the priors
   the chains sample under, with `route = "sample"`: the posterior must
   integrate, which a prior finite at 0 allows only while
   sum(w) - n < 1. frm() still skips it for `dry_run = "objective"`,
   with the comment corrected. With brms's default student_t on sigma,
   where the reviewer measured Rhat 2.497 (`wu-sample-lane.txt`), the
   formula route now refuses.
4. Finding 10: asym_laplace is in the table. Under Laplace
   (`exact = FALSE`) its excess is sum(w), not sum(w) - n.
   `dev/fuzzfix-fr1-families.R` (seed 20261020, n = 80, cap 3000 MB,
   peak 460 MB) on rellib-r10 (`fr1-families-r10.txt`):

       asym_laplace obs unweighted       logLik 448.985 betad -7.10771 conv 1 (false convergence)
       asym_laplace obs w=0.9            logLik 403.443 betad -7.11291 conv 1
       asym_laplace obs quad unweighted  logLik -56.7035 (exact integral)
       asym_laplace obs quad w=1.2       logLik 24.3051 betad -8.0625
       exgaussian obs unweighted         logLik -70.6264 conv 0
       exgaussian obs w=0.9              logLik -63.5638 conv 0
       exgaussian obs w=1.2              ERROR NA/NaN gradient evaluation

   The proof: the asym_laplace log density is linear on each side of
   the mode, so a row's own effect gets no curvature from the data and
   the Laplace Hessian is the prior's alone, while the density at the
   mode grows like 1/sigma; the approximation grows like sigma^-sum(w).
   exgaussian joins with both scales (`sigma`, `beta`): with beta held,
   its density stays below 1/beta. On the lane (`fr1-families-lane.txt`)
   asym_laplace obs is refused under Laplace at every weight, its
   quadrature fits at unit weights and 0.9 are unchanged, and it is
   refused at 1.2; exgaussian is refused at 1.2 and unchanged
   unweighted and at 0.9; every group-of-4 arm is unchanged.

### Other findings

- Finding 7: NEWS and the `bf()` docs give the intercept condition,
  the prior rule and the families. The rewrite says "centered".
- Finding 8: the absolute `abs(betad) < 10` is replaced by convergence
  code 0, no warning, and the weighted log shape within
  |unit-weight log shape| + 1 of the unit-weight fit's. Both core
  files now run in 2 to 3 s, process included (most fits are
  refusals), so they take no `skip_on_cran()`; the frmtmb.sample block
  that samples has one.
- Finding 9: `fuzz_design()` scales by 0.5, so sum(w) < n holds on
  every draw (runif(0.5, 2) / 2 < 1).

### Regression tests, rellib-r10 and lane, one file per process

    RESULT base test-weights-unbounded.R      pass=22 fail=3 err=3 skip=0 warn=8
    RESULT lane test-weights-unbounded.R      pass=40 fail=0 err=0 skip=0 warn=0
    RESULT base test-quadrature-nested.R      pass=4 fail=4 err=1 skip=0 warn=1
    RESULT lane test-quadrature-nested.R      pass=12 fail=0 err=0 skip=0 warn=0
    RESULT base test-quadrature-defects.R     pass=42 fail=1 err=0 skip=0 warn=0
    RESULT lane test-quadrature-defects.R     pass=43 fail=0 err=0 skip=0 warn=0
    RESULT base frmtmb.sample test-weights-unbounded.R pass=3 fail=1 err=0 skip=0 warn=5
    RESULT lane frmtmb.sample test-weights-unbounded.R pass=4 fail=0 err=0 skip=0 warn=0

Base failures (`fr1-tf-*-base.txt`): the Gamma ar1 fit returns "a
<frmtmb_fit> object", not an error; the three prior arms, nested,
crossed and mu-sigma quadrature, and the frm_sample() formula route
fail "to throw a error". The asym_laplace and sampler blocks call
`check_weights_unbounded()` directly, so their base form is the weak
"could not find function"; the asym_laplace block's `frm()` arm is
behavioral. `test-quadrature-defects.R` lost its nested arms (59
assertions in round 0).

## To file: nested quadrature, a design note for a later lane

The exact route is saved as `dev/fuzzfix-nested-quad.patch` against
4854e425 (`git apply --check --cached` passes there):
`quad_passes()`, `quad_check_depth()`, `quad_crossed_stop()`, their
call sites in `fit_assembled()` and `quad_fit()`, the round-0
`quadrature` docs, and the round-0 test as
`tests/testthat/test-quadrature-nested-exact.R`. On top of this round,
remove the `quad_check_blocks()` call first, and fix the pass order.

- Accuracy. Lane: quadrature minus an iterated Gauss-Hermite integral
  in R over frmtmb's own joint objective, 4.794e-09 (cumulative fuzz
  spec) and 4.505e-08 (Bernoulli); the test's weighted Bernoulli,
  1.61e-06 against a Laplace error of 0.853. Reviewer, five models
  (Bernoulli, weighted Poisson, weighted cumulative, an inner random
  slope, `(1 | ga)` in mu and sigma): |quadrature - exact| from
  4.1e-09 to 1.0e-05 in the log-likelihood (the review's table; its
  summary says 1.7e-08 to 1.0e-05), SEs equal to a numerical Hessian
  of the exact
  integral (ratios 0.99998 to 1.00000).
- Memory at depth 2 is linear in rows (reviewer, `cost-nested-*`):
  Bernoulli 12.3 to 19.4 MB per row, 1356 MB at 72 rows and 11016 MB
  at 576; censored Beta on the cdfs merge 551 MB per row, 6612 MB at
  12 rows, killed at the 8 GB cap at 24. One block: 405 MB at 576
  rows. On the merged gated plan, Poisson + cens_chr and Beta + cens
  nested specs passed 8 GB (one past 21 GB uncapped). Depth 3: 378 s
  and 49 GB on 32 Poisson rows (lane).
- Ordering defect. The patch orders passes by level count. A by-split
  inner factor has fewer levels per block than the outer block, so the
  outer block went first: 5 GB on 12 rows, a false "crossed" refusal
  at 24, killed at 6 GB at 48. Inner blocks first (reviewer,
  `dev/fuzzfix-rev-order.R`): 1.9 s and 439 MB, 4.7 s and 534 MB, 12.1
  s and 694 MB. Order by nesting: a block goes before another when the
  rows of each of its levels lie inside one level of the other, a
  pairwise test on the blocks' `Z` columns.
- A memory guard needs a prediction before the iterated tape: the
  Laplace tape's operation count per row (or a probe tape on a few
  rows), a per-family factor per depth calibrated on the points above
  (Bernoulli about 19 MB per row, censored Beta 551 MB), times the
  rows, refused above a budget the user can set. Not built or measured
  here.
- Ideas, unmeasured. The integral factors over outer levels, so one
  tape per outer cluster would bound memory by the largest cluster
  rather than by all rows. The per-row growth suggests that the
  transform copies each row's subgraph for every pair of nodes, which
  a smaller inner rule would cut; TMB's `adaptive = TRUE` failed on
  negbinomial cases in an earlier lane (see the `quad-breakdown` entry
  of `helper-fuzz.R`).

## Fix round 1: suites

All with core and frmtmb.sample from this worktree in the lane
library, the other six extensions from rellib-r10 (each `lib:` line
names the lane's frmtmb).

One file per process, NOT_CRAN (`dev/release/run-par.sh
fuzzfix-fr1-notcran`, "ran 413 of 413", `fr1-notcran-summary.txt`):

    frmtmb files=255 pass=18288 fail=0 err=0 skip=191 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=32 pass=1838 fail=0 err=0 skip=5 warn=0
    frmtmb.latent files=11 pass=368 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=15 pass=430 fail=0 err=0 skip=13 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=62 pass=2853 fail=0 err=0 skip=8 warn=0
    frmtmb.spline files=16 pass=594 fail=0 err=0 skip=1 warn=0

Gated: the 29 core and 4 frmtmb.sample files with a skip above, with
`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_FUZZ` and `FRMTMB_DRMTMB_FIT_TESTS`
true (`fuzzfix-fr1-gated.log`, "ran 33 of 33"): pass=3192 fail=0
err=0 skip=1 warn=0. The one skip is frmtmb.sample's
`test-scale.R`, gated on `FRMTMB_SCALE_TESTS`. Among them:

    RESULT test-fuzz.R pass=2 fail=0 err=0 skip=0 warn=0
    RESULT test-drmtmb-agreement.R pass=131 fail=0 err=0 skip=0 warn=0

That gated fuzz run is this worktree's plan. The merged plan's 12
specs on these paths (`dev/fuzzfix-fr1-fuzz-subset.R`: cdfs's helper
read in place, with this lane's `fuzz_design()` and
`FUZZ_KNOWN_REFUSAL`; capped at 4000 MB; `fr1-fuzz-subset-lane.txt`)
ran in 5 s at a peak of 521 MB, against 8 to 21 GB per spec in the
review: 5 KNOWN-REFUSAL (the nested quadrature specs, including
Poisson + cens_chr, seed 20335153), the 6 ar1 + weights specs with no
finding, and 2 REAL-NEW for Beta + cens (seed 20535438), "cens()/trunc()
need a family with a CDF", because this build lacks cdfs's Beta CDF.
On the merged build that spec has a CDF and reaches the quadrature
refusal in `fit_assembled()`; run it there at consolidation.

One session per suite, three orders (`run-sessions.sh fuzzfix-fr1-sess
4`, "ran 24 of 24", untagged LEAK lines: 0):

    SESSION frmtmb forward files=245/255 pass=18288 fail=0 err=0 skip=191 warn=0
    SESSION frmtmb reverse files=245/255 pass=18386 fail=0 err=0 skip=191 warn=0
    SESSION frmtmb preload files=245/255 pass=18482 fail=0 err=0 skip=192 warn=0
    SESSION frmtmb.coupling (all three) files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam (all three) files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.latent (all three) files=11/11 pass=368 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn (all three) files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.ode (all three) files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample (all three) files=57/62 pass=2853 fail=0 err=0 skip=8 warn=0
    SESSION frmtmb.spline (all three) files=16/16 pass=594 fail=0 err=0 skip=1 warn=0

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`; each log prints `R: .../Rob0.3.26-lapack | LAPACK: 3.12.0`;
the copy of R was deleted afterwards):

    RESULT ob0326 test-weights-unbounded.R pass=40 fail=0 err=0 skip=0 warn=0
    RESULT ob0326 test-quadrature-nested.R pass=12 fail=0 err=0 skip=0 warn=0
    RESULT ob0326 test-quadrature-defects.R pass=43 fail=0 err=0 skip=0 warn=0
    RESULT ob0326-sample test-weights-unbounded.R pass=4 fail=0 err=0 skip=0 warn=0

## Fix round 2 (review, "Re-check after fix round 1")

Logs are `dev/fuzzfix-out/fr2-*`.

1. N1: `prior_spike_verdict()` now checks whether the spike dpar's
   design has an intercept column (a constant non-zero column). When
   it has none and still spans the constant (cell means,
   `sigma ~ 0 + fac`), `prior_cellmeans_verdict()` reads that dpar's
   class `b` entries. It applies the link-scale test of an Intercept
   prior to each column's entry (coef-specific first, then class-wide).
   - All columns bound: "bounds".
   - No `b` entry: "free" (refused).
   - Anything else: "unknown" (fails open).

   The reviewer's `dev/fuzzfix-rev-r1-extra.R` (`fr2-r1-extra-{base,lane}.txt`):

       no prior                       base ok logLik 415.712 betad -8.457 -37.06 conv 1 | lane REFUSED
       normal(0,1) on b, dpar sigma   base ok logLik 9.56287 betad -11.59 -12.28 conv 0 | lane ok logLik 9.56287 betad -11.59 -12.28 conv 0

   The round-1 lane refused the second arm
   (`dev/fuzzfix-rev-out/r1-extra-lane.txt`, line 4). Every other arm
   of that script is unchanged from the review's lane log. A new test
   arm in "a dispersion without an intercept cannot reach its limit"
   covers both: refused without a prior, and fits with conv 0 under
   `normal(0, 1)` on class b, dpar sigma. The NEWS bullet now says that
   a `normal()` prior on every coefficient of an intercept-free
   cell-means formula bounds the model, and a prior on the slopes of a
   formula with an intercept does not.
2. N2: the message and both NEWS files no longer name `lognormal()` or
   `inv_gamma()`.
   - The message names the remedies that give a usable fit: weights
     summing to the number of rows or fewer, a random effect that is
     not observation-level, or a fixed dispersion.
   - It adds that a prior can make the objective bounded but not the
     estimate informative (`fr2-msg.txt`).
   - frmtmb.sample's bullet cites the review's lognormal(0, 1) run:
     posterior median sigma 7.2e-05, Rhat 2.83.
   - New test assertions: no "lognormal" in the message, and "not the
     estimate informative" in it.
3. Finding 8: the slack is now `abs(ru)` alone. The weighted log shape
   must lie within the unit-weight log shape's own magnitude of it. On
   rellib-r10 with the test's data the values are 0.666 against 0.896
   (`dev/fuzzfix-fr2-shape.R`).
4. Cosmetic: a dpar's block is printed as `sigma: (1 | ga)`.
   - frmtmb labels such blocks "sigma: 1 | ga", and the prefix is now
     kept outside the parentheses (`block_term_as_written()`).
   - New test assertion in `test-quadrature-nested.R`.
5. N3 is filed below, not fixed.

Tests, one file per process (`fr2-tf-*.txt`):

    RESULT lane test-weights-unbounded.R pass=45 fail=0 err=0 skip=0 warn=0
    RESULT lane test-quadrature-nested.R pass=13 fail=0 err=0 skip=0 warn=0
    RESULT lane test-quadrature-defects.R pass=43 fail=0 err=0 skip=0 warn=0
    RESULT lane-sample test-weights-unbounded.R pass=4 fail=0 err=0 skip=0 warn=0
    RESULT base test-weights-unbounded.R pass=24 fail=4 err=3 skip=0 warn=10

## To file: student_t priors on an upward dispersion (N3)

A student_t prior on Gamma's shape or Beta's phi gets the verdict
"unknown", so the check fails open and the runaway is not refused.
Review, `dev/fuzzfix-rev-r1-priors.R`:

- Gamma: log shape 33.76, logLik 352.966.
- Beta: log phi 39.11, logLik 29663.6.

The growth is a power, about k^(excess / 2) for the shape families
(k^excess for weibull and frechet). A student_t tail is k^-(nu + 1).
So the verdict can be decided from nu, with the family's exponent:
"free" when nu + 1 is below the growth exponent, "bounds" when it is
above.

The fix-round-1 core check that stopped at "Running 'testthat.R'"
with no status line (`fr1-check-frmtmb-run1.txt`) was not this lane's
defect. The coordinator reports that lane cdfs stopped it by mistake at
about 19:51, when it matched processes on `frmtmb_0.72.0.tar.gz`; that
lane's matching is now limited to its own worktree's paths. From fix
round 2, `dev/fuzzfix-check.sh` passes the tarball by its absolute path
under this worktree (`.../frmtmb-wt-fuzzfix/dev/fuzzfix-check/`), so its
command line carries this lane's path. A substring match on the bare
tarball name would still find it; the protection is cdfs's narrowed
matching. The tarball keeps the name R CMD build gives it.

Fix round 2, one-session runs (`run-sessions.sh fuzzfix-fr2-sess 4
frmtmb frmtmb.sample`, "ran 6 of 6", untagged LEAK lines: 0):

    SESSION frmtmb forward files=245/255 pass=18294 fail=0 err=0 skip=191 warn=0
    SESSION frmtmb reverse files=245/255 pass=18392 fail=0 err=0 skip=191 warn=0
    SESSION frmtmb preload files=245/255 pass=18488 fail=0 err=0 skip=192 warn=0
    SESSION frmtmb.sample forward|reverse|preload files=57/62 pass=2853 fail=0 err=0 skip=8 warn=0

`R CMD check --as-cran`, each alone: core `Status: 1 NOTE` (the V8
math-rendering NOTE on the HTML manual, `fr2-check-frmtmb.txt`);
frmtmb.sample `Status: OK` (`fr2-check-frmtmb.sample.txt`). The check
directory is deleted.
