# Lane gpfield: gp fields and gp parity

Round opened 2026-11-13, base ec8178d3 (frmtmb 0.82.0). Worktree
`frmtmb-wt-gpfield`, branch `wt-r83-gpfield`, library `wt-gpfield-lib`;
"base" is `rellib-r21`. Scripts are `dev/gpfield-*.R` and
`dev/gpfield-*.sh`; logs are in `dev/gpfield-out/` (gitignored). Every
number below is pasted from the log named beside it. Nothing committed.

Status: fix round done (section below). Installed in wt-gpfield-lib (install14). Nothing committed.

## Fix round (review dev/reviews/2026-11-13-gpfield.md: FIX)

This section holds where it disagrees with the sections after it. Machine
notice in force: at most 2 R processes, no probe above 4 GB. No run failed
on disk or allocation in this round.

### F1. The gp() warnings on the fit frm() returns

`gp_end_warn(fit)` (R/gp-start.R) says the nugget warning and
`gp_end_messages()` once per fit (a flag in `fit$cache`, `gp_end_said`),
and returns the fit. R/fit.R calls it in two places: in the fit-end block
of `fit_assembled()` under `announce_start` (where the nugget warning and
the loop were), and around `nl_fit_probed()`'s result in `frm()`:

    gp_end_warn(nl_fit_probed(run, frame, start, prior))

So whatever fit `nl_fit_probed()` returns, the first fit, a saddle
restart's from `nl_fit_off_saddle()`, or the sign probe's refit, is
checked, and the first fit is not checked twice. R/par-template.R is not
touched. After lane nlstop2 merges, nothing moves: its new return points
(`pe$fit` from `nl_pole_escape()`, twice, beside `limit()`; `pr$fit` after
`limit()`) are all inside `nl_fit_probed()`, whose result this one call in
`frm()` reads. The gp() warnings then follow every message and warning
`nl_fit_probed()` says, in that order.

The relin is not run on a refit: the cold-start rule in `fit_assembled()`
(`relin = gp_body && announce_start && cold_start`) is the settled one
(lane-env: do not fork it), and a refit is handed a start. `refit()`
and the other refit paths (`influence()`, `anova()`'s refits) call
`fit_assembled()` unannounced and are not wrapped, so their output gets no
gp() warning, by design: only `frm()`'s returned fit is checked.

The reproducer (`dev/gpfield-smoke.R lane sweep:mm_V:bernoulli:40:1`,
`smoke-f1.txt`):

    SMOKE sweep:mm_V:bernoulli:40:1 ll -17.607794 conv 0 nwarn 1 wall 15.56 | frmtmb_gp_laplace
      W frmtmb_gp_laplace: ... its first-order correction is 60.2 in log-likelihood, more than 1 ...

(the review: code 0 and no warning on base and on install7). New test
"a refit the sign probe returns is checked too" (test-gp-field.R,
`skip_on_cran()`), with the said-once check. NEWS and the 5112 CLOSED
line say so.

The gp end checks evaluate the objective; `gp_end_messages()` now puts
the objective's state back (`obj_state_save()`/`obj_state_restore()`).

### F2. predict()'s cost on an exact field

`gp_field_draw()` takes the replicate's inner mode by Newton steps on
numbers (`gp_field_newton()`, R/predict-brms.R) where one exact gp() block
is the model's only random effect, in a linear predictor that is no
body's, each row reading one position: from the fit's field, `f <- C (C +
W^-1)^-1 (f + g / w)`, `g` and `w` each position's gradient and curvature
of the log-likelihood by central differences of the family's density, `C`
the replicate's kernel. One step on a gaussian or lognormal identity mean
(exact there); elsewhere until a step is under 1e-2 of the field (Newton
converges quadratically, so the field returned is about 1e-4 off), at most
8 steps, else the tape's solve. The factor of the kernel it takes for z is
left on the replicate's copy of the block (`gp_factor_cache`, read by
`gp_factor_ordered()` when the theta is identical), so the replicate's
field does not factor it again. Elsewhere (another random effect, a body,
a curvature not positive) the tape's solve stays, with the objective's
state saved and restored by the helpers (m1).

Accuracy against the tape's inner mode (`dev/gpfield-newton.R poisson
300`, 20 parameter draws from `vcov(full = TRUE)`, `set.seed(5)`):

    one Newton step only (newton-p300.txt):  z off by up to 1.748 posterior sds
    shipped rule, 1e-2 (newton5-p300.txt):   z off by about 2e-3 posterior sds at most (median 2e-5 or less; this run 5.5e-4), the field by 8.3e-5 of its sd

(a 1e-8 rule read 7.3e-6, the tape's own tolerance). The predictions are
the tape's: `sd` of the six columns 2.576 1.304 2.567 1.091 2.112 14.045
on both routes (the review's `cost-lane-p300.txt`, `cost-lane8-p300.txt`),
and on the gaussian 0.317 0.313 0.289 0.381 1.191 3.376 on both.

Cost (`dev/gpfield-cost.R`, the review's script; base and lane in two
processes at once, same seed; per draw from the second of two calls, s):

| fit | call | base | lane, review (install7) | lane, fixed | fixed / base |
|---|---|---|---|---|---|
| gaussian n 500 | newdata 6 rows | 0.0879 | 0.451 | 0.0764 | 0.87 |
| gaussian n 500 | no newdata | 0.0346 | 0.392 | 0.0452 | 1.31 |
| poisson n 300 | newdata 6 rows | 0.0234 | 0.280 | 0.0328 | 1.40 |
| poisson n 300 | no newdata | 0.0110 | 0.239 | 0.0200 | 1.82 |

(`cost-base7-g500.txt` with `cost-lane7-g500.txt`, `cost-base8-p300.txt`
with `cost-lane8-p300.txt`.) Before the factor cache, the poisson
replicates took 3.73 factors of n x n each, z's included
(`newton-bench.txt`), `gp_corr()` 5.6 ms and a factor 2.85 ms at 300
positions. Without the cache the no-newdata poisson row read 4.4 times
base (`cost-lane4-p300.txt` against `cost-base4-p300.txt`), with it 2.0
(`cost-lane7-p300.txt`) and 1.82 in the paired run above. The short check
(dev/reviews/2026-11-13-gpfield-short.md) found that on binomial fields
with many trials 20 to 38 percent of the replicates fall back (their
steps do not settle or a curvature is not positive) to the tape solve,
which costs what it did before
(the review's 3 to 21 times); `?predict.frmtmb_fit` says so. `?predict.frmtmb_fit` and NEWS
give the range and say where the tape's solve stays (the review's 3 to 21
times).

The draw's law (the review's `gpfield-review-predict.R`, rerun on the
fixed build): exp_quad, 4000 draws, seed 7 (`law-eq-4000-s7.txt`), z sd
0.992 to 0.998 and KS p 0.649 to 0.989 at x 10.5 to 60, pair correlation
0.997 against 0.997; matern32, 2000 draws (`law-m32-2000.txt`), z sd 0.966
to 0.998, KS p 0.111 to 0.781, pair 0.973 against 0.975. The same as the
review's, line for line.

### F3. NEWS

The `frmtmb_gp_laplace` and `frmtmb_gp_basis` bullets are BREAKING. The
basis text says "every basis function gets the same prior weight" (code,
`?frmtmb-gp`, NEWS), and NEWS gives the review's rate (13 of 30). Filed
MEDIUM: the SE check should treat that limit as a limit.

### Runs of the fix round (install14)

The gp test files, gated tier, one per process, 2 at a time
(`dev/release/gpfield-fx.log`, "GPFIELD-FX ran 9 of 9", every `lib:` line
wt-gpfield-lib):

    RESULT test-gp-field.R pass=22 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-kernels.R pass=39 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-multidim.R pass=67 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-predict-draw.R pass=7 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-by.R pass=86 fail=0 err=0 skip=0 warn=0
    RESULT test-marginaleffects-re.R pass=59 fail=0 err=0 skip=0 warn=0
    RESULT test-brms-lik-gp-exact.R pass=74 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-start.R pass=251 fail=0 err=0 skip=0 warn=0
    RESULT test-predict-re-uncertainty.R pass=71 fail=0 err=0 skip=0 warn=0

The new test seen failing on rellib-r21 (`dev/release/gpfield-fieldbase.log`,
test-gp-field.R `pass=3 fail=8 err=5`): "a refit the sign probe returns
is checked too" fails on the missing warning ("expected a warning
containing 'first-order correction'"), then on the missing symbol.

OpenBLAS 0.3.32 with LAPACK, `OPENBLAS_CORETYPE=Haswell`
(`bash dev/ciharden-openblas.sh 0.3.32 lapack`,
`openblas-build-0332.txt`: Rblas 76 symbols and Rlapack 617 to OpenBLAS;
logs `gpfield-out/ob0.3.32/`, every `blas:` line the copy and Haswell):

    RESULT test-gp-field.R pass=21 fail=0 err=0 skip=0 warn=1   (first run)
    RESULT test-gp-predict-draw.R pass=7 fail=0 err=0 skip=0 warn=0

The first run of test-gp-field.R there read `pass=21 fail=0 err=0 skip=0
warn=1`: the F1 row's refit stops at code 8 under that BLAS (a Newton step
from it gains 1e-9) and says so beside the Laplace warning. Both are
right; the test now allows the code-8 note. Rerun
(`ob0.3.32/test-gp-field-2.txt`): `pass=21 fail=0 err=0 skip=0 warn=0`
(21: st80 ends at its spike under OpenBLAS, so its correction assertion
does not run, as under 0.3.26).

### Minor items

- m1: `gp_field_draw()`'s tape route uses `obj_state_save()` /
  `obj_state_restore()`.
- m2, the gate: kept at 1, with the rate recorded (NEWS, `?frmtmb-gp`).
  The review's 20 bernoulli fits of `2 sin(3x)` (`lap-white.txt`): 2
  warn, tk 1.51 and 1.08 against importance sampling 0.90 and 0.70; five
  more read 0.78 to 0.90. At a gate of 1.5 the two would go quiet, and so
  would wig160 (tk 1.0015, importance sampling 0.99) and exp_b student
  n120 s2 (tk 1.37, importance sampling 1.11), where the Laplace value is a
  nat off. The white-field artifacts read 60 to 2660 and stay far above
  either gate. A warning on an error of 0.7 to 0.9 says nothing false
  ("can differ"); missing errors of 1 to 1.1 would.
- m3: the inner-Hessian branch now runs only where every random effect
  sits in a predictor that holds an exact field, and `?frmtmb-gp` states
  each branch's scope.
- m4: test-gp-by.R's `fm` allows only the basis text.
- m5: not changed (the review's suggestion): the relin's regular check
  does not ask the Laplace correction; no sweep row moved up and newly
  warned.
- m6: the filed 5319 item points at parity2's `draws_latres_cols()` and
  silent2's `draws_brms_order()`, to do after both merge.
- m7, the relin's cost as the review measured it (`relincost-*.txt`, one
  run per arm, separate processes; identical log-likelihoods):

  | fit | base s | lane s | ratio |
  |---|---|---|---|
  | gaussian `gp(x, k = 30)`, n 5000 | 11.19 | 11.06 | 0.99 |
  | gaussian `gp(x, x2, k = 12)`, n 5000 | 219.39 | 241.75 | 1.10 |
  | poisson `gp(x, k = 30)`, n 5000 | 10.18 | 19.36 | 1.90 |
  | poisson exact `gp(x)`, n 400 | 46.19 | 60.95 | 1.32 |

  The approximate poisson row nearly doubles for no gain: the grid reads up
  to 55 objectives, each an inner solve over 5000 rows. Not changed in
  this round.

### Core suite and check (fix round)

Core, one file per process, gated tier, `FRMTMB_RELLIB=wt-gpfield-lib`
(install14): 33 files on one process (`gpfield-out/fxcore1/`), then the
other 300, 2 at a time (`dev/release/gpfield-fxcore2.log`,
"GPFIELD-FXCORE2 ran 300 of 300"); the first driver was stopped mid-file to
run 2 at a time, and that file reran in the second. Tally
(`gpfield-out/fxcore-sum.txt`):

    frmtmb files 333 results 333 noresult 0 pass 27546 fail 0 err 0 skip 15 warn 0 badlib 0

The extensions were not rerun in the fix round: none calls a function
this round changed (`gp_end_warn()`, `gp_field_draw()`, the factor cache
read in `gp_factor_ordered()` only on a predict() replicate's copy of a
block), and their run of the first submission stands.

`R CMD check --as-cran` of core (`bash dev/gpfield-check.sh .`, install14's
tree; kept as `gpfield-out/check2-frmtmb.txt` and
`check2-testthat.Rout.txt`; the check tree and `dev/ciharden-out/` are
deleted):

    Running 'testthat.R' [30m]
    [ FAIL 0 | WARN 0 | SKIP 274 | PASS 24200 ]
    * checking HTML version of manual ... [28s] NOTE   (V8 math rendering)
    Status: 1 NOTE

### Files touched in the fix round

- R/gp-start.R: `gp_end_warn()`; `gp_end_messages()` puts the
  objective's state back; `gp_laplace_message()`'s inner-Hessian scope;
  the basis warning's text.
- R/fit.R: the fit-end block calls `gp_end_warn(fit)` (replacing the
  nugget warning and the loop there); `frm()` wraps `nl_fit_probed()` in
  `gp_end_warn()`.
- R/predict-brms.R: `gp_field_draw_plan()` (the Newton plan),
  `gp_field_newton_plan()`, `gp_field_newton()`,
  `gp_field_newton_steps`, `gp_field_draw()` (Newton first, the factor
  cache, the state helpers); the cost paragraph of `?predict.frmtmb_fit`.
- R/covstruct.R: `gp_factor_ordered()` reads `gp_factor_cache`.
- R/gp-by.R: `?frmtmb-gp` (the prior-weight wording, each branch's scope,
  the correction's rate).
- NEWS.md: BREAKING on the two warning bullets, the predict() cost, the
  basis rate, the refit check.
- tests/testthat/test-gp-field.R (the F1 test, `mm_V` in `gfd_sweep()`),
  test-gp-by.R (`fm`).
- dev/test-backlog.md: the 5112 CLOSED line, the 5319 item, a MEDIUM item
  (the SE check at the basis's white limit).
- man/frmtmb-gp.Rd, man/predict.frmtmb_fit.Rd (roxygenised).
- dev/gpfield-cost.R (the review's script), dev/gpfield-newton.R,
  dev/gpfield-newton-bench.R, dev/gpfield-newton-prof.R,
  dev/gpfield-predict-prof.R, dev/gpfield-fx*.jobs.

## Items

| backlog line | item | result |
|---|---|---|
| 5112 | white-noise field at code 0, silent | closed: the 2 artifacts (bernoulli) warn `frmtmb_gp_laplace`; the 3 poisson ones are correct fits and stay quiet |
| 5139 | approximate gp() in a body, degenerate, silent | closed: grid at the optimum (2 short rows reach their best start); the white maximum warns `frmtmb_gp_basis` |
| 5149 | `pois_mid two 80` 1.59 and `gsigma step 40` 1.45 short | pois_mid fixed (mean's working marginal read again); gsigma filed again (not a field start) |
| 5125, 5158 | student() step field: Laplace limit; spikes | closed as silent: a first-order Laplace correction warns (st80 3.46 against IS 2.72); spikes named |
| 1906 | Matern and exponential kernels refused | closed: brms's four kernels, exact and approximate, to 1.2e-10 against brms's code; gated test against brms's program |
| 1919 | core predict() holds a gp() curve at its mode | closed: drawn; within 13 percent of brms's sd at 5 positions, where 0.82.0 was 8x too wide inside and half outside |
| 5319 | zgp order against brms | not changed: measured (same field, z differs by up to 1.54); filed with a design |

## Before: the rows on rellib-r21

`dev/gpfield-probe.R base <tags> <ndraw>`: each row's default fit, then
15 theta starts (sd 0.3, 1.5, 6 by length scale 0.1, 0.3, 1, 3, 8 in x's
units), the best and the best "regular" one (code 0, no warning, every
exact length scale at least half the median nearest-neighbour gap), and
the marginal by importance sampling (`gf_is()`, lane gpstop's
`is_ll()`: multivariate t, df 4, proposal precision the inner Hessian's
eigenvalues floored at 1, 20000 draws, `set.seed(9001)`).

5112's rows (`probe-white-base.txt`):

    PROBE sweep:exp_a:bernoulli:120:2 default ll -63.973108 conv 0 nwarn 0 wall 15.26 inner_ev 1 | gp sd 0.9536 ls 1.338e-05 ls/mednn 0.000535 |
    PROBEREF sweep:exp_a:bernoulli:120:2 default -63.973108 best -63.973108 (-0.0000) | best regular -67.067961 (3.0949) | gp sd 1.176e-05 ls 0.5158 ls/mednn 20.6
    PROBEIS sweep:exp_a:bernoulli:120:2 default  laplace -63.9731 is -104.1235 (halves -103.4519 -107.2812) ess 2 of 20000 ev_min 1
    PROBEIS sweep:exp_a:bernoulli:120:2 bestreg  laplace -67.0680 is -67.0655 (halves -67.0603 -67.0708) ess 4856 of 20000 ev_min 1
    PROBE sweep:logis_a:poisson:40:1 default ll -76.324095 conv 0 nwarn 0 wall 0.53 inner_ev 1 | gp sd 0.3221 ls 0.003458 ls/mednn 0.031 |
    PROBEREF sweep:logis_a:poisson:40:1 default -76.324095 best -76.471708 (0.1476) | best regular -76.471708 (0.1476) | gp sd 0.2707 ls 0.3625 ls/mednn 3.25
    PROBEIS sweep:logis_a:poisson:40:1 default  laplace -76.3241 is -76.3447 (halves -76.3576 -76.3320) ess 7729 of 20000 ev_min 1
    PROBE sweep:logis_a:bernoulli:120:2 default ll -64.476101 conv 0 nwarn 0 wall 28.46 inner_ev 1 | gp sd 31.05 ls 0.001677 ls/mednn 0.0493 |
    PROBEREF sweep:logis_a:bernoulli:120:2 default -64.476101 best -64.476101 (-0.0000) | best regular -74.087044 (9.6109) | gp sd 3.245e-05 ls 47.2 ls/mednn 1.39e+03
    PROBEIS sweep:logis_a:bernoulli:120:2 default  laplace -64.4761 is -111.7197 (halves -111.0266 -124.1782) ess 1 of 20000 ev_min 1
    PROBEIS l96:B:wig:60:2 default  laplace -118.6294 is -118.6593 (halves -118.6601 -118.6584) ess 6602 of 20000 ev_min 1.01
    PROBEIS l96:B:wig:60:2 best     laplace -118.6111 is -118.6349 (halves -118.6358 -118.6341) ess 6594 of 20000 ev_min 1
    PROBE l96:B:wig:60:1 default ll -120.294447 conv 0 nwarn 0 wall 0.61 inner_ev 1 | gp sd 0.1163 ls 0.004018 ls/mednn 0.0738 |
    PROBEIS l96:B:wig:60:1 default  laplace -120.2944 is -120.3033 (halves -120.3042 -120.3024) ess 6675 of 20000 ev_min 1

So of the five white-noise stops, the two on bernoulli rows are Laplace
artifacts (importance sampling 40 and 47 below the Laplace value, with an
effective sample size of 1 or 2, so those two are not numbers, only "far
below"; a regular maximum reads the same by both), and the three on
poisson rows are not: a white field there is an observation-level effect,
and importance sampling agrees with the Laplace value to 0.009 to 0.030.
A rule on the length scale alone would warn on three correct fits.

5139's and 5149's rows (`probe-approx-base.txt`):

    PROBE gate:approx_a default ll 40.926376 conv 0 nwarn 0 wall 0.82 inner_ev 9.24 | hsgp sd 103.8 ls 3.995e-05 ls/mednn 0.00125 |
    PROBEREF gate:approx_a default 40.926376 best 40.926376 (-0.0000) | best regular 36.042782 (4.8836) | hsgp sd 0.6452 ls 9.136 ls/mednn 285
    PROBE gate:approx_b default ll 100.401246 conv 0 nwarn 0 wall 0.72 inner_ev 1.32e+03 | hsgp sd 0.04822 ls 2.211 ls/mednn 85 |
    PROBEREF gate:approx_b default 100.401246 best 101.584647 (-1.1834) | best regular 101.584647 (-1.1834) | hsgp sd 0.05127 ls 0.5634 ls/mednn 21.7
    PROBE sweep:pow_a:gaussian:120:2:12 default ll 2.257621 conv 0 nwarn 0 wall 0.21 inner_ev 36.3 | hsgp sd 59.21 ls 3.083e-05 ls/mednn 0.00121 |
    PROBEREF sweep:pow_a:gaussian:120:2:12 default 2.257621 best 3.193095 (-0.9355) | best regular 3.193095 (-0.9355) | hsgp sd 0.7877 ls 9.032 ls/mednn 354
    PROBE ng:pois_mid:two:80:1 default ll -202.602532 conv 0 nwarn 0 wall 0.48 inner_ev 1 | gp sd 0.3453 ls 0.8446 ls/mednn 17 |
    PROBEREF ng:pois_mid:two:80:1 default -202.602532 best -201.011189 (-1.5913) | best regular -201.011189 (-1.5913) | gp sd 0.3979 ls 0.1016 ls/mednn 2.05
    PROBE ng:gsigma:step:40:1 default ll 34.412387 conv 0 nwarn 0 wall 0.15 inner_ev 1 | gp sd 0.7925 ls 0.5653 ls/mednn 8.96 |
    PROBEREF ng:gsigma:step:40:1 default 34.412387 best 35.861222 (-1.4488) | best regular 35.861222 (-1.4488) | gp sd 0.7552 ls 0.7996 ls/mednn 12.7

approx_a's white field (gaussian, where the Laplace value is exact) IS the
likelihood's maximum: the basis of k = 20 cannot draw the wiggle at any
regular length scale. approx_b and pow_a k12 are short stops.

## A cheap detector for the Laplace limit (5125, 5158, and 5112)

`dev/gpfield-tk.R` (prototype, no frmtmb change): the first-order
correction to the Laplace approximation, latres_tk()'s formula, with
each row's derivatives of minus its log density in the field's linear
predictor `u` by a 5-point stencil of the family's own density
(row_lpdf() with the predictor's offset shifted) and `v = A H^-1 A'`.
Against importance sampling at the fit's own point
(`dev/gpfield-tkrun.R base ...`, `tk1-base.txt`, `tk2-base.txt`):

    TK ng:student:step:80:1 ll 68.5503 conv 0 nwarn 0 | tk 2.1736 | is 71.2654 ess 75 is-lap 2.7151 ev_min 0.999
    TK ng:student:step:80:1@1/6 ll 79.5450 conv 1 nwarn 2 | tk 2.6e23 | is 73.8216 ess 38 is-lap -5.7234 ev_min 4.17e-08
    TK ng:student:step:160:1 ll 201.2938 conv 0 nwarn 0 | tk 1.7227 | is 203.9653 ess 55 is-lap 2.6715 ev_min 1
    TK sweep:logis_a:poisson:40:1 ll -76.3241 conv 0 nwarn 0 | tk -0.0228 | is -76.3447 ess 7729 is-lap -0.0206 ev_min 1
    TK ng:pois_mid:two:80:1 ll -202.6025 conv 0 nwarn 0 | tk 0.0043 | is -202.5861 ess 5694 is-lap 0.0164 ev_min 1
    TK ng:bern:wig:80:1 ll -54.7853 conv 0 nwarn 0 | tk 0.2400 | is -54.5513 ess 5888 is-lap 0.2339 ev_min 1
    TK ng:negbin:mix:80:1 ll -559.7495 conv 0 nwarn 0 | tk 0.0228 | is -559.7035 ess 5298 is-lap 0.0460 ev_min 1
    TK sweep:exp_a:bernoulli:120:2 ll -63.9731 conv 0 nwarn 0 | tk 2660.4277
    TK sweep:logis_a:bernoulli:120:2 ll -64.4761 conv 0 nwarn 0 | tk -23.2396
    TK l96:B:wig:60:2 ll -118.6294 conv 0 nwarn 0 | tk -0.0230
    TK ng:bern:step:40:1 ll -26.9051 conv 0 nwarn 0 | tk 0.0207
    TK ng:binom5:wig:40:1 ll -64.7063 conv 0 nwarn 0 | tk 0.1857

(stencil step 1e-2 in `u` in the prototype; the package's step is 5
percent of each row's posterior sd of `u`, which reads st80 at 3.46.)
The correction tracks importance sampling where the Laplace value holds
(0.004 to 0.24 against 0.016 to 0.234), reads 1.7 and 2.2 where it is 2.7
off (st80, st160), and is far from 0 at the spike and at both bernoulli
white fields. It separates every case of the three items from every
correct fit measured, at the residual-term note's gate of 1.

So 5125's "no check at the stop separates it" is not so: the stop is a
regular maximum of the Laplace objective, but the objective there is
2.7 below the marginal, and the correction says so. It does not say
whether the marginal's maximum is elsewhere (st160's is not; st80's is),
only that the Laplace value is not the marginal there.

## 1906: brms's Matern and exponential kernels

`gp(cov = )` takes `"matern52"`, `"matern32"`, `"exponential"` (and
`"matern12"`) beside `"exp_quad"`, exact and approximate. brms 2.23.0's
definitions, read in `dev/brms-suite/brms/R/formula-gp.R` and
`inst/chunks/fun_gp_*.stan`, `fun_spd_gp_*.stan`:

- exact, isotropic: the kernel of the Euclidean distance over the
  dimensions; anisotropic: the PRODUCT of one kernel per dimension
  (`cov_gp()`, and the Stan chunk multiplies `gp_*_cov(x[, d], 1,
  lscale[d])`);
- approximate: the spectral densities `spd_gp_*()`, which for an
  anisotropic Matern field weigh each frequency by its own length scale
  (an ARD Matern), so brms's exact and approximate anisotropic Matern
  fields are different kernels. frmtmb follows brms in both.

Code: `R/covstruct.R` (`gp_cov_of()`, `gp_kern_s()`, `gp_kern_dlog()`,
`gp_kern_d2()`, `gp_kern_xy()`, `gp_corr()`, `gp_field_dtheta()`,
`gp_cross_cov()`, `hsgp_sds(cov =)`, `hsgp_log_spd_matern()`),
`R/predict.R` (the three kriging kernels: `gp_krig_cov()`,
`gp_krig_factor()`, `gp_krig_factor_dense()`), `R/frame.R` (`gp_cov` on
the block), `R/parse.R` (the refusal removed), `R/extra-deriv.R`
(`frm_extra_cov_deriv()` refuses a Matern term by name: its derivatives
are not the Hermite polynomials it takes), `R/compat.R` (the `gp_pred()`
family row's text). The squared exponential keeps 0.82.0's arithmetic on
every path (a branch, not a rewrite), so its fits do not move.

On the tape a Matern kernel reads `r / rho` with the distance `r` taken
on numbers: a root on the tape at distance 0 has no derivative.

**Against brms's own code** (`dev/gpfield-kernels.R lane`,
`kernels1.txt`): frm()'s gaussian fit, then the model's marginal at
frm()'s estimates computed by hand from brms's `cov_gp()`, or from
`spd_gp()` at `standata()`'s `slambda_1` and basis `Xgp_1`, on brms's
scaled inputs; and the kernel itself against brms's at one theta. Seed
`set.seed(9100 + case)`, 40 rows:

    KERN exact1    exp_quad    frm ll -37.26460657 brms-kernel ll -37.26460657 diff 9.45e-13 | kernel rel diff 3.98e-16 | conv 0 nwarn 0
    KERN exact1    matern52    frm ll -28.33024786 brms-kernel ll -28.33024786 diff -2.13e-14 | kernel rel diff 5.06e-16 | conv 0 nwarn 0
    KERN exact1    matern32    frm ll -38.12988069 brms-kernel ll -38.12988069 diff 4.83e-13 | kernel rel diff 7.24e-16 | conv 0 nwarn 0
    KERN exact1    exponential frm ll -40.13234144 brms-kernel ll -40.13234144 diff -2.84e-14 | kernel rel diff 7.11e-16 | conv 0 nwarn 0
    KERN exact2iso exp_quad    frm ll -20.86499579 brms-kernel ll -20.86499579 diff -2.13e-14 | kernel rel diff 3.74e-16 | conv 0 nwarn 0
    KERN exact2iso matern52    frm ll -15.28334960 brms-kernel ll -15.28334960 diff -2.49e-14 | kernel rel diff 4.81e-16 | conv 0 nwarn 0
    KERN exact2iso matern32    frm ll -21.78009696 brms-kernel ll -21.78009696 diff -1.42e-14 | kernel rel diff 4.49e-16 | conv 0 nwarn 0
    KERN exact2iso exponential frm ll -26.44536499 brms-kernel ll -26.44536499 diff -1.16e-10 | kernel rel diff 1.57e-16 | conv 0 nwarn 0
    KERN exact2ani exp_quad    frm ll -9.21646742 brms-kernel ll -9.21646742 diff -3.55e-15 | kernel rel diff 3.57e-16 | conv 0 nwarn 0
    KERN exact2ani matern52    frm ll -8.98036783 brms-kernel ll -8.98036783 diff -1.24e-13 | kernel rel diff 5.97e-16 | conv 0 nwarn 0
    KERN exact2ani matern32    frm ll -15.95026172 brms-kernel ll -15.95026172 diff 7.07e-13 | kernel rel diff 4.25e-16 | conv 0 nwarn 0
    KERN exact2ani exponential frm ll -18.44508252 brms-kernel ll -18.44508252 diff 8.17e-14 | kernel rel diff 2.29e-16 | conv 0 nwarn 0
    KERN hsgp1     exp_quad    frm ll -38.62286498 brms-kernel ll -38.62286498 diff -7.11e-15 | kernel rel diff 4.19e-16 | conv 0 nwarn 0
    KERN hsgp1     matern52    frm ll -32.97500984 brms-kernel ll -32.97500984 diff 0.00e+00 | kernel rel diff 7.23e-16 | conv 0 nwarn 0
    KERN hsgp1     matern32    frm ll -36.26477846 brms-kernel ll -36.26477846 diff -7.11e-15 | kernel rel diff 2.60e-16 | conv 0 nwarn 0
    KERN hsgp1     exponential frm ll -39.52201092 brms-kernel ll -39.52201092 diff -7.11e-15 | kernel rel diff 4.70e-16 | conv 0 nwarn 0
    KERN hsgp2ani  exp_quad    frm ll -23.80058938 brms-kernel ll -23.80058938 diff 9.80e-12 | kernel rel diff 2.39e-16 | conv 0 nwarn 0
    KERN hsgp2ani  matern52    frm ll -25.12091663 brms-kernel ll -25.12091663 diff 1.72e-12 | kernel rel diff 1.95e-15 | conv 0 nwarn 0
    KERN hsgp2ani  matern32    frm ll -29.16127470 brms-kernel ll -29.16127470 diff 2.58e-12 | kernel rel diff 8.46e-15 | conv 0 nwarn 0
    KERN hsgp2ani  exponential frm ll -27.88991859 brms-kernel ll -27.88991859 diff 1.83e-12 | kernel rel diff 3.75e-15 | conv 0 nwarn 0

20 of 20 agree to 1.2e-10 in the log-likelihood and 8.5e-15 in the
kernel (relative). (The kernel comparison for the 2-D approximate rows
sorts the two sets of spectral values: brms orders its multi-indices
otherwise.)

## 5319: the order of an exact gp()'s z against brms

`dev/gpfield-zgp.R` (`zgp.txt`): one field, frmtmb's fit of `y ~ gp(x)`
on test-brms-lik-gp-exact.R's data (13 integer positions, unsorted), its
z against brms's z for the same field (`L_b^-1 f / sdgp`, `L_b` the
factor over brms's positions in `standata()`'s order):

    ZGP positions 13 | frmtmb order of the positions 7,1,13,4,10,2,3,5,6,8,9,11,12
    ZGP brms order of the positions (sorted index) 7,6,4,9,10,8,11,5,1,3,2,13,12
    ZGP max |frmtmb zgp - brms zgp| 1.5390 | the two fields' max difference 4.44e-16 | sum of squares 7.792759 7.792759

(per-element lines in the log). Same field, same prior term (an
orthogonal map between the two z), elementwise different by up to 1.54.

Not changed. Two routes, each larger than this item's size:

- Fit in brms's order. The farthest-first order exists because sorted
  or appearance orders round the field by 1.7e-5 at n 15 where it
  rounds by 6.1e-11, which a sampler cannot mix through
  (`gp_farthest_first()`, dev/ncgp-findings.md). Not an option.
- Report brms's z as a derived column of frmtmb.sample's draws: per
  draw, `L_b(theta_draw)^-1 sdgp L_ff(theta_draw) z / sdgp`, one
  Cholesky per draw, with the sampled z kept under an internal name that
  every reader of the draws (`draws_internal_matrix()` and the chart
  inverses) takes back. frmtmb.sample's derived columns read only the
  outer parameters (`draws_cov_cols()` via `hyp_env_vals()`), so this
  is a new kind of derived column, in the draws layer lanes parity2 and
  silent2 change in this round. Filed with this design.

The mapping itself exists and is tested: `stan_pars_from_fit()` in
`tests/testthat/helper-brms.R` builds brms's zgp this way for every
gated likelihood check.

## 1919: predict() draws a gp() field

`predict_b_drawer()` (`R/predict-brms.R`, predraw2's file; my hunks are
the `gov` line, the `xn` line and the `gp_field_draw()` call in the
drawer, the `krige_draw` flag, and three new functions after it:
`gp_drawn_b()`, `gp_field_draw_plan()`, `gp_field_draw()`; plus two doc
paragraphs) now draws:

- every `gp()` block's coefficients, exact and approximate, at every
  `re_formula` (a gp() term is a population term, which `re_formula`
  does not drop), from the conditional law the group effects take;
- an exact field's z at the inner mode of the replicate's own
  parameters (one inner solve per replicate, `obj$fn()` with the
  objective's state put back), plus the drawn deviation. The first-order
  mode shift the group effects take is not enough for a field, below;
- at a new position, the field given the replicate's field
  (`krige_draw`, as frmtmb.sample's draws already do).

REML, profile, quadrature and importance fits keep the first-order z
draw (gp_field_draw_plan() returns NULL; filed LOW).

**Against brms 2.23.0** (`dev/gpfield-predict-brms.R`, 30 rows,
`set.seed(9200)`, `y = 2 sin(x) + N(0, 0.2)` on 0 to 10; frmtmb
`predict(ndraws = 4000)` after `set.seed(1)`, arms in one process, base
= `gp_drawn_b()` mocked to return nothing; brms default priors, 2 chains
of 2000, seed 1, 6 divergent transitions; logs `predict-brms.txt` for
the first draft and `predict-brms2.txt` for the shipped code):

    frm sdgp 2.800249 lscale(data units) 2.180201 sigma 0.1987712
    brms sdgp 2.923566 lscale 0.2007489 sigma 0.2096869     (lscale on brms's scaled x)
    PRED x   5.5 | sd lane 0.246 base 1.887 brms 0.242 | 95% width lane 0.968 base 7.419 brms 0.933
    PRED x  10.5 | sd lane 0.671 base 2.034 brms 0.633 | 95% width lane 2.636 base 8.091 brms 2.494
    PRED x  12.0 | sd lane 2.277 base 3.018 brms 2.188 | 95% width lane 8.904 base 12.052 brms 8.676
    PRED x  15.0 | sd lane 3.794 base 3.186 brms 3.373 | 95% width lane 15.302 base 12.845 brms 13.711
    PRED x  60.0 | sd lane 3.913 base 1.786 brms 3.612 | 95% width lane 15.364 base 7.081 brms 14.666

0.82.0 was not only narrow far out (1.786 against 3.612 at x 60, the
item); inside the positions it was 8 times too wide (1.887 against
0.242 at x 5.5). The field held at its mode z was read through each
replicate's own sd and length scale while the intercept, which the field
can absorb (its standard error 1.75), was drawn on its own.

The first draft took the field to first order, `c0 + Jt dtheta + Jz dz`
read back through the replicate's factor (`predict-brms.txt`): right at
the positions (0.224 at x 5.5) and far away (3.913 at x 60), but 52.2 at
x 15, where kriging from a field that the replicate's kernel does not
draw amplified it. The shipped draw takes the replicate's own inner
mode, which that kernel does draw.

The approximate field (`dev/gpfield-dbg4.R`, same data, `k = 15`;
`dbg4-base.txt`, `dbg4-lane.txt`): variance at x 5.5, 0.82.0 3.208,
lane 0.0516, the noise alone (`propagate_error = FALSE`) 0.0414. The
same cause: its basis weights held while the intercept was drawn.

Cost: one inner solve per replicate for an exact field; the replicate
already factors the kernel once (coef_b() at the drawn theta).

## The fit items: what changed (R/gp-start.R, two lines of R/fit.R)

1. **5149: the mean's field read again at the fit's optimum**
   (`gp_mean_relin()`, gpstop's candidate built). Where an exact gp()
   sits in a linear predictor that is no body's and the model's marginal
   did not apply (any family but the gaussian with a constant sigma),
   each row's log density is replaced by its quadratic in its predictor
   `u` at the inner mode: `y_w = u + g / w` with known noise `1 / w`, `g`
   and `w` by central differences of the family's own density
   (`gp_row_ll_shift()`: `row_lpdf()` with the predictor's offset
   shifted, on a fit-like list from `solved_par_list()`), and the
   marginal of that (`gp_marginal_opt(lrs = 0)`: the noise is known, so
   the noise-ratio grid has one value). A start where it reads another
   length scale (a factor 1.1) goes through gpstop's loop: a descent,
   kept only at a higher regular maximum (`gp_relin_regular()`), up to
   3 rounds. Same gate as the body's: the user's own cold start, no
   `start`, natural beta. It does not reach `gsigma step 40` (filed
   again): there the working marginal is the model's at the fit's sigma
   coefficients, and it reads the fit's own theta (`dbg1.txt`: relin
   theta -0.2325667 -0.5704147, the fit's); the higher maximum also moves
   sigma's z coefficient, -2.13 against -0.77, and the other theta with
   the fit's other values reads 6.19, 28 below (`dbg2.txt`, `dbg3.txt`).
2. **5139: the approximate field's grid at the optimum**
   (`gp_hsgp_relin()`), for a model with approximate fields and no exact
   one: length scales from the positions' span (`gp_span`, new on the
   block) down to half the basis's spacing (`gp_spacing = 2 L / k`, new
   on the block), sds the fit's own times 0.3, 1, 3 and the predictor's
   unit times 0.1, 1, read from the objective with every other value
   where the fit put it; a descent from the best, kept only at a higher
   regular maximum, where regular for this block also means not white
   over the basis. A first draft read sds from the unit only and missed
   `b ~ gp(x, k = 20)` (the field's sd there is 0.05; `dbg1.txt` GRID:
   ls 0.05 sd 0.05 reads 101.08 against the fit's 100.40).
3. **5139: `frmtmb_gp_basis`** (`gp_basis_message()`): an approximate
   field whose length scales all lie under a tenth of the basis's
   spacing (`gp_basis_white()`, `gp_basis_white_frac = 0.1`; at a tenth
   the squared exponential spectral density falls 5 percent across the
   basis). approx_a's maximum is there (gaussian, Laplace exact; the best
   regular start is 4.88 below), so k = 20 cannot draw the field. The
   gated test-brms-lik-gp-exact.R row `y ~ gp(x, k = 10)` on 13 integer
   positions ends there too, and the profile (`dev/gpfield-dbg6.R base`,
   `dbg6-base.txt`) is flat: -39.8719 at ls 1e-4, -39.8973 at 0.4 of the
   spacing, -40.0330 at 0.8, -42.894 at 1.6. So do test-brms-likelihood.R's
   rows 10 and 10c, on 80 and 60 distinct positions (`dev/gpfield-dbg7.R
   lane`, `dbg7-lane.txt`): row 10 -29.5025 at the limit, -29.8737 at 0.4
   of the spacing, -40.2231 at 1.6; row 10c -27.7552, -28.3808,
   -42.2202. A modest k ends at the limit often: in the first suite pass
   (`gpfield-out/gated-part1-files/`) the warning escaped 9 times from 4
   files (test-brms-likelihood.R 2, test-gp-by.R 4, test-gp-multidim.R 2,
   test-marginaleffects-re.R 1; k from 4 to 15), besides the approximate
   row of test-brms-lik-gp-exact.R. Every one profiled is a
   likelihood maximum at that limit (Laplace exact on the gaussian), so
   the warning says only what holds there: sd and the length scales are
   not determined apart (only `sdgp^2` times the product of the length
   scales is), the fitted field and its predictions hold, and a larger k
   or the exact field reads the smoothness. It does not claim the field
   is rough. It is not said where the field's sd ran to 0 (the field's
   fitted values under 1e-6 of its predictor's): the field is absent
   there, the standard-error check names that, and the length scale is
   moot. Those test sites allow it, each with the reason.
4. **5112, 5125, 5158: `frmtmb_gp_laplace`** (`gp_laplace_message()`):
   (a) where the inner Hessian at the optimum has an eigenvalue below
   0.01 (`gp_inner_pd()`, the relin rule's cut), naming it (5158); (b)
   else where `|gp_laplace_tk()| > 1`. The check runs on one response
   whose random effects all sit in the field's predictor, with at most
   2000 distinct rows of its design (the sums run over distinct rows),
   and not where every field is in the mean of a gaussian or lognormal
   response with the identity link and no censoring or truncation, where
   the correction is 0 exactly. Stencil step: 5 percent of each row's
   posterior sd of `u`.

Package tk against importance sampling at the fit's own point, rows the
150 sweep newly warned and two that stay quiet (`dev/gpfield-tkrun.R
lane`, `tk3-lane.txt`, 20000 draws, `set.seed(9001)`):

    TK ng:student:two:160:1 ll 185.4971 conv 0 nwarn 1 | pkg tk 1.8973 | is 187.1699 ess 912 is-lap 1.6729 ev_min 0.953
    TK ng:student:wig:160:1 ll 195.9667 conv 0 nwarn 1 | pkg tk 1.0015 | is 196.9571 ess 1938 is-lap 0.9905 ev_min 1
    TK ng:student:step:80:1 ll 68.5503 conv 0 nwarn 1 | pkg tk 3.4581 | is 71.2654 ess 75 is-lap 2.7151 ev_min 0.999
    TK ng:student:mix:160:1 ll 260.7212 conv 0 nwarn 0 | pkg tk 0.2606 | is 261.0303 ess 1072 is-lap 0.3092 ev_min 1
    TK ng:student:smooth:160:1 ll 321.0538 conv 0 nwarn 0 | pkg tk 0.0834 | is 321.1468 ess 4052 is-lap 0.0930 ev_min 1

Every warned row is one where the Laplace value is about a nat or more
off the marginal (wig 160 sits at the gate: tk 1.0015, importance
sampling 0.9905). No row warned where importance sampling and Laplace
agree within 0.31.

## Tests and seen failing

New files `tests/testthat/test-gp-kernels.R` (7 tests),
`test-gp-field.R` (6), `test-gp-predict-draw.R` (2); test-gp-by.R's
refusal of `cov = "matern32"` became a fit; test-brms-lik-gp-exact.R
gains a gated test (the three kernels, exact iso, exact anisotropic and
approximate, against brms's compiled program, `brms_lp_check()`), and
its approximate row allows the basis warning (above).

Lane (`dev/release/gpfield-new.log`, `FRMTMB_RELLIB=wt-gpfield-lib`):

    RESULT test-gp-kernels.R pass=39 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-field.R pass=17 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-predict-draw.R pass=7 fail=0 err=0 skip=0 warn=0

rellib-r21 (`dev/release/gpfield-newbase.log`):

    RESULT test-gp-kernels.R pass=4 fail=0 err=7 skip=0 warn=0
    RESULT test-gp-field.R pass=3 fail=6 err=4 skip=0 warn=0
    RESULT test-gp-predict-draw.R pass=2 fail=5 err=0 skip=0 warn=0

The behavioral failures on rellib-r21: -202.6 < -201.0 (pois_two), 100.4
< 101.6 and 2.26 < 3.19 (the approximate fields), no basis warning, no
Laplace warning on st80 nor on the bernoulli white field, predict far sd
1.81 <= 2.25, at a fitted position 1.9 >= 0.6 (3 sigma), the approximate
field's variance 3.2 >= 3 times the noise's; the kernel tests stop at
the refusal ("not implemented here yet").

## Sweeps: the lane against 0.82.0's route, interleaved

`dev/gpfield-sweep.R <set> <part>`: one process on wt-gpfield-lib, arms
`lane` and `base` row by row, the order alternating; `base` mocks
`gp_mean_relin_applies()` and `gp_hsgp_relin()` to do nothing and
`gp_end_messages()` to say nothing, which is 0.82.0's route (the
squared exponential kernel keeps 0.82.0's arithmetic). Summaries by
`dev/gpfield-sweepsum.R` (`sum-*.txt`); MOVE lists every row whose
log-likelihood, code or warning set differs. The sweeps ran on the
installs install4/install5 (the 150 rows on install3); the code since
changes only predict()'s draw, the basis warning's text and the docs
(no value on these rows), and between install3 and install4 only the
Laplace check's sums over distinct rows (an identity) and its skip of a
gaussian mean (where it is 0).

Lane arstep's 150 rows (`sum-ng150.txt`; references are lane gpstop's
best of 15 theta starts, `bestc` the best that ended at code 0 with no
warning):

    SUM ng150 rows 150 (lane 150 base 150) | up >1e-6 1 >1e-3 1 | down < -1e-6 0 < -1e-3 0 | max up 1.5913 min 0.0000
    SUM ng150 conv changes 0 | warning-set changes 5 | wall lane 566.1 base 507.5 ratio 1.116
    MOVE ng150 pois_mid two 80 1 | lane -201.011189 c0 w0 [] | base -202.602532 c0 w0 [] | d +1.591343
    MOVE ng150 student step 160 1 | lane 201.293775 c0 w1 [frmtmb_gp_laplace] | base 201.293775 c0 w0 [] | d +0.000000
    MOVE ng150 student step 40 1 | lane 47.593930 c1 w3 [frmtmb_gp_laplace;frmtmb_conv_code;frmtmb_conv_gradient] | base 47.593930 c1 w2 [frmtmb_conv_code;frmtmb_conv_gradient] | d +0.000000
    MOVE ng150 student step 80 1 | lane 68.550327 c0 w1 [frmtmb_gp_laplace] | base 68.550327 c0 w0 [] | d +0.000000
    MOVE ng150 student two 160 1 | lane 185.497063 c0 w1 [frmtmb_gp_laplace] | base 185.497063 c0 w0 [] | d +0.000000
    MOVE ng150 student wig 160 1 | lane 195.966651 c0 w1 [frmtmb_gp_laplace] | base 195.966651 c0 w0 [] | d +0.000000
    REF ng150 arm l rows 150 | short of best 4 (silent 1) worst -10.995 | short of bestc 1 (silent 1) worst -1.449
    REF ng150 arm b rows 150 | short of best 5 (silent 4) worst -10.995 | short of bestc 2 (silent 2) worst -1.591
    REFSHORT ng150 gsigma step 40 1 | lane 34.4124 w0 [] | base 34.4124 | bestc 35.8612 best 35.8612
    REFSHORT ng150 pois_mid two 80 1 | lane -201.0112 w0 [] | base -202.6025 | bestc -201.0112 best -201.0112

The review's 122-row sweep, pb part (`sum-pb.txt`), gauss part
(`sum-gauss.txt`):

    SUM pb rows 56 (lane 56 base 56) | up >1e-6 0 >1e-3 0 | down < -1e-6 0 < -1e-3 0 | max up 0.0000 min 0.0000
    SUM pb conv changes 0 | warning-set changes 8 | wall lane 2129.5 base 2112.2 ratio 1.008
    MOVE pb exp_a_bernoulli_n120_s2 | lane -63.973108 c0 w1 [frmtmb_gp_laplace] | base -63.973108 c0 w0 [] | d +0.000000
    MOVE pb exp_b_bernoulli_n120_s2 | lane -64.797158 c1 w3 [frmtmb_gp_laplace;frmtmb_conv_code;frmtmb_conv_gradient] | base -64.797158 c1 w2 [frmtmb_conv_code;frmtmb_conv_gradient] | d +0.000000
    MOVE pb logis_a_bernoulli_n120_s2 | lane -64.476101 c0 w1 [frmtmb_gp_laplace] | base -64.476101 c0 w0 [] | d +0.000000
    MOVE pb mm_K_bernoulli_n120_s1 | lane -79.035250 c1 w3 [frmtmb_gp_laplace;frmtmb_conv_code;frmtmb_conv_gradient] | base -79.035250 c1 w2 [frmtmb_conv_code;frmtmb_conv_gradient] | d +0.000000
    MOVE pb mm_K_bernoulli_n120_s2 | lane -78.473077 c1 w3 [frmtmb_gp_laplace;frmtmb_conv_code;frmtmb_conv_gradient] | base -78.473077 c1 w2 [frmtmb_conv_code;frmtmb_conv_gradient] | d +0.000000
    MOVE pb mm_K_bernoulli_n40_s2 | lane -27.301011 c1 w3 [frmtmb_gp_laplace;frmtmb_conv_code;frmtmb_conv_gradient] | base -27.301011 c1 w2 [frmtmb_conv_code;frmtmb_conv_gradient] | d +0.000000
    MOVE pb mm_K_poisson_n40_s2 | lane -78.828493 c1 w3 [frmtmb_gp_laplace;frmtmb_conv_code;frmtmb_conv_gradient] | base -78.828493 c1 w2 [frmtmb_conv_code;frmtmb_conv_gradient] | d +0.000000
    MOVE pb pow_a_bernoulli_n40_s2 | lane -21.263602 c0 w2 [frmtmb_gp_laplace;frmtmb_conv_gradient] | base -21.263602 c0 w1 [frmtmb_conv_gradient] | d +0.000000
    SUM gauss rows 35 (lane 35 base 35) | up >1e-6 1 >1e-3 1 | down < -1e-6 0 < -1e-3 0 | max up 0.9355 min 0.0000
    SUM gauss conv changes 0 | warning-set changes 0 | wall lane 306.4 base 306.5 ratio 1.000
    MOVE gauss pow_a_gaussian_n120_s2_k12 | lane 3.193095 c0 w0 [] | base 2.257621 c0 w0 [] | d +0.935474

fitedge's 192 rows, fresh set (`sum-fresh.txt`):

    SUM fresh rows 96 (lane 96 base 96) | up >1e-6 0 >1e-3 0 | down < -1e-6 0 < -1e-3 0 | max up 0.0000 min 0.0000
    SUM fresh conv changes 0 | warning-set changes 1 | wall lane 564.0 base 549.3 ratio 1.027
    MOVE fresh J sin 60 1 | lane -37.501990 c1 w2 [frmtmb_gp_laplace;frmtmb_conv_code] | base -37.501990 c1 w1 [frmtmb_conv_code] | d +0.000000

The `base` arm is 0.82.0: on the rows probed on rellib-r21 itself
(`probe-white-base.txt`, `probe-approx-base.txt`) the two read the same
log-likelihood to the printed 1e-6: exp_a bernoulli n120 s2 -63.973108,
logis_a poisson n40 s1 -76.324095, logis_a bernoulli n120 s2 -64.476101,
pois_mid two 80 -202.602532, gsigma step 40 34.412387, pow_a k12
2.257621.

## Roundings

`dev/gpfield-ulps.R` (`ulps.txt`): each test-gp-field.R assertion on 21
copies of its data, `y * (1 + 2e-15 k)` for k = 0 to 20 (`x * (1 + 2e-15
k)` for the count rows, whose response has no ulps):

    ULPSUM pois_two draws 21 fail 0
    ULPSUM approx_b draws 21 fail 0
    ULPSUM pow_a draws 21 fail 0
    ULPSUM st80 draws 21 fail 0
    ULPSUM st40 draws 21 fail 0
    ULPSUM logis_pois draws 21 fail 0

fitedge's lane96 set (`sum-lane96.txt`):

    SUM lane96 rows 96 (lane 96 base 96) | up >1e-6 0 >1e-3 0 | down < -1e-6 0 < -1e-3 0 | max up 0.0000 min 0.0000
    SUM lane96 conv changes 0 | warning-set changes 0 | wall lane 379.4 base 375.7 ratio 1.010

The 122 sweep's student part, with the 4 n = 400 rows (`sum-student.txt`):

    SUM student rows 32 (lane 32 base 32) | up >1e-6 0 >1e-3 0 | down < -1e-6 0 < -1e-3 0 | max up 0.0000 min 0.0000
    SUM student conv changes 0 | warning-set changes 1 | wall lane 3346.1 base 3184.2 ratio 1.051
    MOVE student exp_b_student_n120_s2 | lane 55.945521 c0 w1 [frmtmb_gp_laplace] | base 55.945521 c0 w0 [] | d +0.000000

The one row of the 122 that newly warns at code 0 besides 5112's two,
exp_b student n120 s2, against importance sampling (`tk4-lane.txt`):

    TK sweep:exp_b:student:120:2 ll 55.9455 conv 0 nwarn 1 | pkg tk 1.3678 | is 57.0599 ess 943 is-lap 1.1144 ev_min 0.999

The Laplace value is 1.11 below the marginal there: the warning holds.
Every other newly warned row (6 pb rows, 1 fresh row) already ended at
code 1 or with a gradient warning.

## Cost

Wall time of the interleaved sweeps, lane over base (one process each,
under the load of 14 lanes): the 150 rows 1.116 (566.1 s against 507.5
s; the mean's second reading runs on every non-gaussian `y ~ gp(x)`
fit), the 122 sweep's pb part 1.008, gauss 1.000, student 1.051,
fitedge's fresh 1.027 and lane96 1.010. The Laplace check took 0.02 to
0.22 s per fit (`tk1`, `tk3`, `tk4` logs, column after vmax).

## Not done, and why

- 5319: see its section (filed with a design).
- `gsigma step 40` (5149's second row): not a field start (filed).
- The Laplace check's scope (one response, all random effects in the
  field's predictor): filed LOW. A field beside a group term elsewhere
  would need the cross terms of `v`.
- predict() under REML or profile keeps the first-order z draw: filed
  LOW, not measured.
- Derivatives of Matern fields (`frm_extra_cov_deriv()`): refused by
  name, filed LOW.
- 0.3.32: the new assertions are not on a rounding-decided branch (21
  roundings each, above; the student rows' tk sits at 3.46 against a
  gate of 1, the bernoulli row's at 147 or more), so only 0.3.26 ran.

## Files touched

- `R/gp-start.R` (lane-owned): `gp_marginal_opt(lrs =)`;
  `gp_start_scan()` (the approximate field's closure, the mean's
  reading beside the body's); new `gp_mean_relin_applies()`,
  `gp_mean_relin_keys()`, `gp_relin_pfit()`, `gp_row_ll_shift()`,
  `gp_mean_relin()`, `gp_relin_start()`, `gp_hsgp_relin()`,
  `gp_basis_white()`, `gp_basis_white_frac`, `gp_basis_message()`,
  `gp_laplace_tk()`, `gp_laplace_max_rows`, `gp_laplace_tk_max`,
  `gp_laplace_exact()`, `gp_laplace_message()`, `gp_end_messages()`;
  `gp_relin_regular()` reads an approximate block.
- `R/gp-by.R` (lane-owned): `?frmtmb-gp` sections "Kernels" and "The
  Laplace check", the `k` paragraph, two aliases.
- `R/covstruct.R`: the kernels (1906), `hsgp_sds(cov =)` and its three
  callers.
- `R/predict.R` (shared): the three kriging kernels, one hunk each in
  `gp_krig_cov()`, `gp_krig_factor()`, `gp_krig_factor_dense()` (new
  argument `krig`).
- `R/predict-brms.R` (predraw2's): `predict_b_drawer()`'s `gov` line,
  `xn` line, the field draw call and `krige_draw`; new `gp_drawn_b()`,
  `gp_field_draw_plan()`, `gp_field_draw()` after it; two roxygen
  paragraphs (`man/predict.frmtmb_fit.Rd` regenerates).
- `R/fit.R` (shared): three lines after `gp_nugget_message()` in the
  fit-end block (`gp_end_messages()`).
- `R/frame.R` (shared): `gp_cov`, `gp_spacing`, `gp_span` on the gp and
  hsgp components and the merged block.
- `R/parse.R`: the kernel refusal removed.
- `R/extra-deriv.R`: the Matern refusal.
- `R/compat.R`: the `gp_pred()` `kind:family` row's text (one line). No
  other compat row changes.
- `NEWS.md`: a development section, "Gaussian processes", 6 bullets.
- `man/frmtmb-gp.Rd`, `man/predict.frmtmb_fit.Rd` (roxygenised).
- `tests/testthat/`: new `test-gp-kernels.R`, `test-gp-field.R`,
  `test-gp-predict-draw.R`; `test-gp-by.R` (the refusal became a fit;
  its two gated tests allow the basis warning);
  `test-brms-lik-gp-exact.R` (a gated test; the approximate row allows
  the basis warning); `test-brms-likelihood.R` (rows 10 and 10c),
  `test-gp-multidim.R` and `test-marginaleffects-re.R` allow the basis
  warning, each with its reason.
- `dev/test-backlog.md`: CLOSED lines under 1906, 1919, 5112, 5125,
  5139, 5149, 5158; section "Filed in the round of 2026-11-13", lane
  gpfield, 5 items.
- `dev/gpfield-*` (scripts, this file).

## Bump and floor

Minor: three BREAKING bullets (fit results on non-gaussian exact fields
and approximate fields, predict() draws) and a new feature (kernels). No
extension calls a changed function's signature in a way that moves (the
new `hsgp_sds()` argument has a default; the extensions call none of the
changed functions), so no floor moves.
    ULPSUM exp_bern draws 21 fail 0

## OpenBLAS 0.3.26 with LAPACK

`bash dev/ciharden-openblas.sh 0.3.26 lapack` (`openblas-build-0326.txt`:
Rblas 74 symbols and Rlapack 616 to OpenBLAS), then
`OB_N=1 bash dev/gpfield-openblas-run.sh 0.3.26` (`ob0326-sum.txt`,
logs `gpfield-out/ob0.3.26/`; every `blas:` line names the copy, every
`lib:` line wt-gpfield-lib):

    RESULT test-gp-ar-start.R pass=49 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-by.R pass=71 fail=0 err=0 skip=2 warn=0
    RESULT test-gp-field.R pass=16 fail=1 err=0 skip=0 warn=3
    RESULT test-gp-inner-pd.R pass=10 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-kernels.R pass=39 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-multidim.R pass=67 fail=0 err=0 skip=0 warn=2
    RESULT test-gp-noncentered.R pass=76 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-predict-draw.R pass=7 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-start.R pass=251 fail=0 err=0 skip=0 warn=0
    RESULT test-predict-re-uncertainty.R pass=71 fail=0 err=0 skip=0 warn=0

Two findings, both fixed in the tests:

- Under OpenBLAS st80's default fit ends at the spike (code 8, gradient
  1.49e10), where on reference BLAS it ends at the regular point. The
  warning is there (the inner-Hessian form), but the test required the
  correction form. The test now accepts either end and requires the form
  that end calls for (`spike <- !gp_inner_pd(...)`). Rerun:
  `test-gp-field.R pass=18 fail=0 err=0 skip=0 warn=0` under OpenBLAS
  (`ob0.3.26/test-gp-field-2.txt`), `pass=19` on reference BLAS
  (`dev/release/gpfield-field2.log`; one more assertion where the row
  is regular).
- test-gp-multidim.R's `k = 4` over two dimensions ends at the basis's
  white limit under OpenBLAS and warns `frmtmb_gp_basis`; that test reads
  the argument evaluation, so it now allows the warning. Rerun: `pass=67
  fail=0 err=0 skip=0 warn=0` (`ob0.3.26/test-gp-multidim-2.txt`).

## Suites

Every test file of core and the eight extensions, one file per process,
`NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`, `FRMTMB_RELLIB=wt-gpfield-lib`
ahead of rellib-r21, in two passes: `dev/release/run-par.sh
gpfield-gated dev/gpfield-all.jobs 2` on install6 (stopped after 187
files, kept as `gpfield-out/gated-part1-files/`, to reinstall with the
basis warning's fixes), then `run-par.sh gpfield-gated2
dev/gpfield-rest.jobs 3` on install7 (`dev/release/gpfield-gated2.log`,
"GPFIELD-GATED2 ran 354 of 354"): the 345 files not yet run, and 9 rerun
(every file of the first pass that warned, and the gp files). Between
install6 and install7 only the basis warning changed (its text, and not
said where the field is absent). The tally, each file's last result
(`gpfield-out/suite-sum.txt`; `badlib` counts files whose `lib:` line
does not name wt-gpfield-lib):

    frmtmb files 333 results 333 pass 27543 fail 0 err 0 skip 15 warn 0 badlib 0
    frmtmb.coupling files 11 results 11 pass 551 fail 0 err 0 skip 5 warn 0 badlib 0
    frmtmb.eam files 33 results 33 pass 1857 fail 0 err 0 skip 3 warn 0 badlib 0
    frmtmb.latent files 13 results 13 pass 697 fail 0 err 0 skip 2 warn 0 badlib 0
    frmtmb.learn files 16 results 16 pass 545 fail 0 err 0 skip 2 warn 0 badlib 0
    frmtmb.ode files 11 results 11 pass 549 fail 0 err 0 skip 1 warn 0 badlib 0
    frmtmb.sample files 91 results 91 pass 4603 fail 0 err 0 skip 1 warn 0 badlib 0
    frmtmb.spatial files 8 results 8 pass 612 fail 0 err 0 skip 0 warn 0 badlib 0
    frmtmb.survival files 16 results 16 pass 575 fail 0 err 0 skip 0 warn 0 badlib 0

The 15 core skips are the gates as before: test-drmtmb-agreement.R 13,
test-curve-scale.R 1, test-fuzz.R 1; each extension's skips are its
test-scale.R. No extension was reinstalled: none changed.

The files this lane added or changed, on install7 (from the tally above):

    RESULT test-gp-field.R pass=19 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-kernels.R pass=39 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-predict-draw.R pass=7 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-by.R pass=86 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-multidim.R pass=67 fail=0 err=0 skip=0 warn=0
    RESULT test-brms-lik-gp-exact.R pass=74 fail=0 err=0 skip=0 warn=0
    RESULT test-brms-likelihood.R pass=544 fail=0 err=0 skip=0 warn=0
    RESULT test-marginaleffects-re.R pass=59 fail=0 err=0 skip=0 warn=0

test-brms-lik-gp-exact.R's gated test is the kernels against brms
2.23.0's compiled program (`brms_lp_check()`, joint density and its
gradient, 9 programs compiled fresh: three kernels, exact isotropic,
exact anisotropic, approximate); its first run on install6 stopped at
`object 'cv' not found` (brms's helpers rebuild a formula away from its
environment), fixed by writing the kernel into each formula.

## R CMD check --as-cran (core)

`bash dev/gpfield-check.sh .` on install7's tree (built and checked in
`dev/gpfield-check/`, `R_LIBS` = wt-gpfield-lib, rellib-r21, the user
library; pandoc and TinyTeX on PATH; `_R_CHECK_CRAN_INCOMING_REMOTE_=
FALSE`). Kept as `gpfield-out/check-frmtmb.txt` and
`check-testthat.Rout.txt`; the check tree and `dev/ciharden-out/` are
deleted.

    Running 'testthat.R' [77m]
    [ FAIL 0 | WARN 1 | SKIP 274 | PASS 24197 ]
    * checking examples ... [101s] NOTE   (frm_hazard_reads 6.08 s elapsed)
    * checking HTML version of manual ... [34s] NOTE   (V8 math rendering)
    Status: 2 NOTEs

The test WARN is the machine's, not the package's: test-importance.R:922,
`cannot create dir '...\Rtmp...' reason 'No space left on device'` in
testthat's own temp directory while 14 lanes ran; the disk read 136 GB
free afterwards, and test-importance.R passes one file per process
(`pass=344 fail=0 err=0 skip=0 warn=0`, the suite above). The examples
NOTE is load (lane-rules: an examples-timing NOTE on this box measures
load); `frm_hazard_reads()` is not this lane's.
No extension changed, so no extension check ran.

## Sessions

`bash dev/release/run-sessions.sh gpfield-sess 3 frmtmb frmtmb.sample`
(gated tier, `FRMTMB_RELLIB=wt-gpfield-lib`, install7; logs kept as
`gpfield-out/sess1/`):

    SESSION frmtmb forward files=333/333 pass=27543 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb reverse files=333/333 pass=27670 fail=0 err=3 skip=15 warn=0
    SESSION frmtmb.sample forward files=91/91 pass=4603 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample reverse files=91/91 pass=4603 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample preload files=91/91 pass=4603 fail=0 err=0 skip=1 warn=0
    GPFIELD-SESS ran 5 of 6

Core's reverse session's 3 errors were `std::bad_alloc` in
test-brms-lik-cdfs.R, and core's preload session died with a
segmentation fault at test-v14.R while the driver's own forks failed
with 0xC000012D (the commit limit: `fork: retry: Resource temporarily
unavailable`), the machine at 8.5 GB available. Neither file is this
lane's and both pass one per process. Lane-rules: a `bad_alloc` under
memory pressure is not evidence until it reproduces alone. The two
sessions ran again, two at a time beside the check
(`dev/release/run-session.R frmtmb <order>`, `gpfield-out/sess2/`):

    SESSION frmtmb preload files=333/333 pass=27811 fail=0 err=0 skip=17 warn=0
    SESSION frmtmb reverse files=333/333 pass=27708 fail=0 err=0 skip=15 warn=0

LEAK lines: core forward 48, reverse 49, preload 38; frmtmb.sample 15,
15, 8; untagged 0 in each of these six. The first, crashed preload
session had one untagged line, `LEAK test-se-tier.R search path:
+AD-overloads -`, which the full preload rerun does not show.
