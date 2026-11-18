# Lane nbrobust: a stable negbinomial density at large size (RTMB-4)

Round of 2026-11-16, base b255350e (frmtmb 0.83.0), worktree
`frmtmb-wt-nbrobust`, library `wt-nbrobust-lib`. Backlog items 4903
(MEDIUM) and 5224 (LOW). Status: DONE; nothing runs.

Logs are under `dev/nbrobust-out/` and `dev/nbrobust-tests/`
(gitignored); every number below is pasted from them.

## Summary

- New `nb_lpmf(y, lmu, lsh)` and `nb_rising(y, lsh)` in `R/families.R`:
  the negative-binomial log mass from the log mean and the log shape,
  `R + y lmu - (shape + y) log1p(mu / shape) - lgamma(y + 1)` with
  `R = lgamma(y + shape) - lgamma(shape) - y log(shape)`, which is
  `sum(log1p(j / shape))` over `0 < j < y` for a count up to 24 (one
  partial-sum vector when the shape is one number) and
  `lgamma_shift_diff(shape, y) - y log(shape)` above.
  `log1p(mu / shape)` is `logspace_add(0, lmu - lsh)`.
- Every family that read the negative-binomial density reads it:
  `negbinomial()` (so `negbinomial2`, through `fam_reparam()`),
  `nbinom1()`, `geometric()`, `zero_inflated_negbinomial()`,
  `hurdle_negbinomial()`. No extension has its own negative-binomial
  density (grep of `extensions/*/R` for `nbinom`/`negbinom`: only the
  frmtmb.sample prior table); frmtmb.survival, latent, learn and sample
  read core's families, so their suites ran (below).
- The non-log mean-link branches of negbinomial, nbinom1, geometric and
  zero_inflated_negbinomial used `RTMB::dnbinom2(y, mu, var)`, which
  forms `log(var - mu)` by subtraction and keeps the same lgamma
  cancellation; they now use `log(mu)` and the same `nb_lpmf()`, as
  hurdle_negbinomial already did. zero_inflated_negbinomial's non-log
  `lp0 = shape * (log(shape) - log(shape + mu))` cancelled at a large
  shape too; it now uses the log-link expression for every link.
- Derivative in the log shape, 3888 points (log shape 0 to 40, counts 0
  to 1e4, means 0.1 to 1e4) against a 256-bit reference: worst error
  6.63 times `eps * (y + mu + 1)` (old: 4.6e15 times, 672 absolute).
- Fits at the poisson limit: 0 code-1 stops and one stop point per data
  set over 126 fits on both BLAS (old: 26 code-1 stops, log shapes 14.99
  to 24.54); the limit message on all 126 under check_se "warning"; no
  warning under "ignore" (old: 2 gradient warnings; from a start at log
  shape 30, 21 of 21).
- Cost: faster. A negbinomial GLMM's objective takes 0.41 times, its
  gradient 0.33 times the 0.83.0 time.
- Items 4903 and 5224 closed. RTMB-4 in `dev/upstream-bugs.md` updated.

## 1. Reproduction on rellib-r22

Script `dev/nbrobust-repro.R <base|lib> <nb200|nb300> 0 20 <check_se>`,
driver `dev/nbrobust-repro-all.sh <arm> <ref|ob326>`, table
`dev/nbrobust-repro-sum.R` (log `dev/nbrobust-out/repro-sum.txt`). Case
nb200 is item 4903's data (seed 2, n = 200, `dev/ci082-probe.R`); case
nb300 is item 5224's data (seed 4, n = 300,
`dev/reviews/famend-review-nbpert.R`), fitted as `negbinomial()` and
`brmsfamily("negbinomial2")`. `x * (1 + 2e-15 * k)`, k = 0 to 20.
OpenBLAS: `dev/ciharden-openblas.sh 0.3.26 lapack`, one thread; nb200
with `OPENBLAS_CORETYPE=Sandybridge` (item 4903's setting), nb300 with
it unset (each log's `R:` and `CORETYPE:` lines).

| log | lib | fits | code 0 | code 1 | limit msg | grad warn | any warn | log shapes (min, max) |
|---|---|---|---|---|---|---|---|---|
| base-ob326-nb200-ignore.txt | rellib-r22 | 21 | 19 | 2 | 0 | 1 | 3 | 14.9899, 24.4791 |
| base-ob326-nb200-warning.txt | rellib-r22 | 21 | 19 | 2 | 21 | 0 | 0 | 14.9899, 24.4791 |
| base-ob326-nb300-ignore.txt | rellib-r22 | 42 | 30 | 12 | 0 | 1 | 12 | 15.9484, 24.2991 |
| base-ob326-nb300-warning.txt | rellib-r22 | 42 | 30 | 12 | 42 | 0 | 0 | 15.9484, 24.2991 |
| base-ref-nb200-ignore.txt | rellib-r22 | 21 | 19 | 2 | 0 | 0 | 2 | 15.0339, 21.1495 |
| base-ref-nb200-warning.txt | rellib-r22 | 21 | 19 | 2 | 21 | 0 | 0 | 15.0339, 21.1495 |
| base-ref-nb300-ignore.txt | rellib-r22 | 42 | 32 | 10 | 0 | 0 | 10 | 15.8933, 24.5351 |
| base-ref-nb300-warning.txt | rellib-r22 | 42 | 32 | 10 | 42 | 0 | 0 | 15.8933, 24.5351 |
| lib-ob326-nb200-ignore.txt | wt-nbrobust-lib | 21 | 21 | 0 | 0 | 0 | 0 | 17.2212, 17.2212 |
| lib-ob326-nb200-warning.txt | wt-nbrobust-lib | 21 | 21 | 0 | 21 | 0 | 0 | 17.2212, 17.2212 |
| lib-ob326-nb300-ignore.txt | wt-nbrobust-lib | 42 | 42 | 0 | 0 | 0 | 0 | 18.9741, 18.9741 |
| lib-ob326-nb300-warning.txt | wt-nbrobust-lib | 42 | 42 | 0 | 42 | 0 | 0 | 18.9741, 18.9741 |
| lib-ref-nb200-ignore.txt | wt-nbrobust-lib | 21 | 21 | 0 | 0 | 0 | 0 | 17.2212, 17.2212 |
| lib-ref-nb200-warning.txt | wt-nbrobust-lib | 21 | 21 | 0 | 21 | 0 | 0 | 17.2212, 17.2212 |
| lib-ref-nb300-ignore.txt | wt-nbrobust-lib | 42 | 42 | 0 | 0 | 0 | 0 | 18.9741, 18.9741 |
| lib-ref-nb300-warning.txt | wt-nbrobust-lib | 42 | 42 | 0 | 42 | 0 | 0 | 18.9741, 18.9741 |

What reproduced, and what did not:

- Item 4903 as filed (0.82.0, check_se "warning" and "ignore" alike):
  on 0.83.0, nb200 k = 15 under OpenBLAS Sandybridge stops at the same
  point, log shape 24.4791, logLik -386.278985178424, code 0, but under
  "warning" it gets the limit message and no warning. 0.83.0's
  BREAKING change ("A gradient above grad_tol along the direction of a
  coefficient ... at the end of its range ... no longer hides the
  limit", NEWS 0.83.0) already reads that stop. Under "ignore" (no SE
  analysis) the false warning remains: "Large maximum absolute gradient"
  on nb200 k = 15 and nb300 k = 9 (OpenBLAS), 2 of 126.
- Item 4903 under "warning" does reproduce on 0.83.0 from a start far
  out (`dev/nbrobust-start.R`, `start = list(betad = 30)`, nb200 data,
  reference BLAS): 15 of 21 roundings warn "Large maximum absolute
  gradient at the optimum (0.191 at shape_(Intercept)) ... would still
  gain 0.105" together with "false convergence (8)"; under "ignore",
  21 of 21.
- Item 5224 as filed (code 1, convergence warnings, no limit message
  under "warning"): not on 0.83.0. nb300 k = 9 under OpenBLAS stops at
  log shape 24.1914, code 1, logLik -555.262460822484, exactly as in the
  famend review, and gets the limit message with 0 warnings: every
  code-1 stop in the table above gets the message under "warning". The
  same 0.83.0 change reads them. Under "ignore" a code-1 stop keeps the
  convergence warning, which is the documented contract of
  `se_boundary_stop()` and the separate user decision filed by lane
  ci082 (backlog, "should check_se = 'ignore' still ask the boundary
  question").
- What 0.83.0 still has is the cause: where the fit stops is decided by
  the noise of the derivative (log shapes from 14.99 to 24.54 across
  roundings) and 26 of 126 fits stop at code 1.

| start log | lib | fits | code 0 | code 1 | limit msg | grad warn | any warn |
|---|---|---|---|---|---|---|---|
| start-base-ref-ignore.txt | rellib-r22 | 21 | 2 | 19 | 0 | 21 | 21 |
| start-base-ref-warning.txt | rellib-r22 | 21 | 2 | 19 | 6 | 15 | 15 |
| start-lib-ob326-ignore.txt | wt-nbrobust-lib | 21 | 21 | 0 | 0 | 0 | 0 |
| start-lib-ob326-warning.txt | wt-nbrobust-lib | 21 | 21 | 0 | 21 | 0 | 0 |
| start-lib-ref-ignore.txt | wt-nbrobust-lib | 21 | 21 | 0 | 0 | 0 | 0 |
| start-lib-ref-warning.txt | wt-nbrobust-lib | 21 | 21 | 0 | 21 | 0 | 0 |

(`dev/nbrobust-out/start-sum.txt`; the lane's start fits stop at log
shape 30.0000, logLik -386.278985178452 or ...453, on every rounding.)

## 2. The density, and its jacobian against the exact derivative

`dev/nbrobust-jac.R`, log `dev/nbrobust-out/jac.txt`, table
`dev/nbrobust-out/jac.tsv`. On the tape (`RTMB::MakeTape`, `jacfun()`),
log shape 0 to 40 by 0.5 (81), counts 0, 1, 2, 3, 5, 10, 24, 25, 30,
100, 1000, 1e4, means 0.1, 3, 100, 1e4: 3888 points per arm. Reference
at 256 bits (Rmpfr): the value from `lgamma()`, the derivative in the
log shape `n (digamma(y + n) - digamma(n) - log1p(mu / n) + (mu - y) /
(n + mu))`, in the log mean `y - (n + y) mu / (n + mu)`. "units" is the
error over `eps * (y + mu + 1)`. Arms: robust = `RTMB::dnbinom_robust(y,
lmu, 2 lmu - lsh)` (0.83.0); scalar and vector = `nb_lpmf()` with one
shape for all rows and with a shape per row (the two branches of
`nb_rising()`).

    robust lsh 0-10     n 1008  max|d dsh| 4.159e-11  max units   4.9e+04  max|d dmu| 2.068e-12  max rel value 5.599e-13
    robust lsh 10.5-18  n  768  max|d dsh| 1.291e-07  max units  1.61e+08  max|d dmu| 3.638e-12  max rel value 6.953e-12
    robust lsh 18.5-24  n  576  max|d dsh| 4.488e-05  max units  2.46e+10  max|d dmu| 1.819e-12  max rel value 4.775e-12
    robust lsh 24.5-30  n  576  max|d dsh| 1.958e-02  max units  2.62e+13  max|d dmu| 3.638e-12  max rel value 9.830e-12
    robust lsh 30.5-35  n  480  max|d dsh| 5.269e+00  max units  4.39e+15  max|d dmu| 2.231e-12  max rel value 1.354e-11
    robust lsh 35.5-40  n  480  max|d dsh| 6.715e+02  max units  4.58e+15  max|d dmu| 2.231e-12  max rel value 1.051e-11
    scalar lsh 0-10     n 1008  max|d dsh| 2.208e-12  max units       2.6  max|d dmu| 2.068e-12  max rel value 5.632e-12
    scalar lsh 10.5-18  n  768  max|d dsh| 4.206e-12  max units      2.55  max|d dmu| 3.638e-12  max rel value 6.241e-12
    scalar lsh 18.5-24  n  576  max|d dsh| 1.346e-11  max units      6.46  max|d dmu| 1.819e-12  max rel value 6.383e-12
    scalar lsh 24.5-30  n  576  max|d dsh| 1.131e-11  max units      5.86  max|d dmu| 3.638e-12  max rel value 8.661e-12
    scalar lsh 30.5-35  n  480  max|d dsh| 1.089e-11  max units      6.63  max|d dmu| 2.231e-12  max rel value 8.334e-12
    scalar lsh 35.5-40  n  480  max|d dsh| 7.535e-12  max units      5.91  max|d dmu| 2.231e-12  max rel value 8.133e-12
    vector (identical to scalar, every band)

    nb_lpmf (scalar) against dnbinom_robust where it is right: 3805 of 3888 points, max rel 4.610e-12
    scalar against vector: max |value| diff 0.000e+00, max |dsh| diff 3.553e-15

    log shape, y = 3, mu = 3: robust, nb_lpmf, reference (d / d log shape)
      0.0  4.47039e-01  4.47039e-01  4.47039e-01
     10.0  6.80916e-05  6.80917e-05  6.80917e-05
     14.0  1.24669e-06  1.24729e-06  1.24729e-06
     18.0 -6.67438e-08  2.28450e-08  2.28450e-08
     20.0  5.81693e-07  3.09173e-09  3.09173e-09
     24.0 -2.22391e-05  5.66270e-11  5.66270e-11
     28.0 -2.52922e-03  1.03786e-12  1.03716e-12
     30.0 -1.35716e-02  1.39824e-13  1.40364e-13
     35.0 -4.00000e+00  7.73000e-16  9.45768e-16
     40.0 -4.00000e+00 -1.27451e-17  6.37253e-18

Reading it:

- The derivative's error no longer grows with the shape: at most 6.63
  units at every band, against 4.6e15 for `dnbinom_robust()`. An
  absolute error of order `eps * (y + mu)` per row is the floor of any
  double formula: `(shape + y) log1p(mu / shape)` and its derivative
  `shape log1p(x) - shape x / (1 + x)` each are of order `mu`, and their
  difference, the true derivative, is of order `mu^2 / shape`. At a log
  shape past about 30 the tape reads that floor (1e-15 against 9e-16
  above), far below any gradient tolerance.
- The value: the largest relative error, 8.7e-12 (old 1.4e-11), is at a
  count and a mean of 1e4, where the log mass is -5.5 and is formed from
  terms of about 9e4 (`y lmu` = 92103, `lgamma(10001)` = 82109); the
  absolute error there is 4.8e-11 (old 7.5e-11). For counts up to 30 the
  relative error is below 1e-14 (old 2.7e-14); at 100 and 1000, 2.3e-13
  and 2.9e-13 (old 2.1e-13 and 4.9e-14). Against `dnbinom_robust()`
  where that one is right (within 1e-13 of the reference, 3805 of 3888
  points) the new value agrees to 4.6e-12 relative, the largest again at
  a count of 1e4; the brief's 1e-12 holds for every count up to 1000.
- The derivative in the log mean is unchanged (same maxima, both arms).

Why not `lgamma_shift_diff()` alone, as the brief suggested: its value
`lgamma(a + s) - lgamma(a)` is right, and so is its derivative, but the
negative binomial then subtracts `y log(shape)` and `y log(shape + mu)`;
writing the mass as `R + y lmu - (shape + y) log1p(mu / shape)` takes
those out before they are formed. `lgamma_shift_diff()` serves the
counts above 24 (`R = lgamma_shift_diff(shape, y) - y lsh`, absolute
error of order `eps * y`, the same floor), as in `lgamma_count_shift()`.

The `osa` route (count on the tape) uses `lgamma_shift_diff()` for every
row; the `simref` route (simulation) goes to `dnbinom_robust()`, which
draws. `dev/nbrobust-osa.R` (logs `osa-base.txt`, `osa-lib.txt`): for
negbinomial, nbinom1, geometric, zinb, hurdle nb and negbinomial with
the identity link on 60 rows (seed 3, two counts of 31 and 45), the
logLik agrees to all 12 printed digits, `simulate()` is identical, and
the one-step residuals ("cdf" and "oneStepGeneric" with
`discreteSupport = 0:120`) agree to 1e-10 in their sum
(-0.9225771187 base, -0.9225771188 lane, negbinomial).

## 3. Cost

`dev/nbrobust-cost.R`. Data: 2000 rows, 50 groups, shape 2, counts up
to 59 (18 rows above 24), seed 1.

Density alone (`dens`, one process, arms interleaved in blocks of at
least 1.2 s, minimum of 5 rounds; control = the robust arm timed again;
log `cost-dens.txt`):

    nb        fn 6.371e-06 s gr 4.517e-04 s  ratio to robust: fn 1.084 gr 0.356
    nb_v      fn 6.561e-06 s gr 4.541e-04 s  ratio to robust_v: fn 1.089 gr 0.387
    robust2   fn 6.142e-06 s gr 1.235e-03 s  ratio to robust: fn 1.045 gr 0.973

The objective's value costs about the same (1.08 against a control of
1.05); its gradient 0.36 times.

Whole GLMM, `y ~ x + (1 | g)` (`glmm`, base and lane alternating, three
processes each; per process the minimum of five rounds over blocks of at
least 1.2 s, against a fixed arithmetic control in the same process;
logs `cost-glmm-{base,lib}-{1,2,3}.txt`; both arms logLik
-4576.3261915086, log shape -0.576956, code 0):

    RESULT glmm base control 2.6562e-03 fn 2.4805e-03 gr 1.2031e-02 fn/control 0.9338 gr/control 4.5294 fit 1.90
    RESULT glmm base control 2.5781e-03 fn 2.4414e-03 gr 1.2813e-02 fn/control 0.9470 gr/control 4.9697 fit 2.38
    RESULT glmm base control 2.5586e-03 fn 2.4414e-03 gr 1.3906e-02 fn/control 0.9542 gr/control 5.4351 fit 2.25
    RESULT glmm lib control 2.5391e-03 fn 1.0156e-03 gr 3.8281e-03 fn/control 0.4000 gr/control 1.5077 fit 1.09
    RESULT glmm lib control 2.5000e-03 fn 9.6680e-04 gr 4.0625e-03 fn/control 0.3867 gr/control 1.6250 fit 1.23
    RESULT glmm lib control 2.6562e-03 fn 1.1914e-03 gr 4.8047e-03 fn/control 0.4485 gr/control 1.8088 fit 1.36

Best over best: objective 0.3867 / 0.9338 = 0.41, gradient 1.5077 /
4.5294 = 0.33. Why the Laplace objective gains more than the density's
own value (1.08) is not measured; it differentiates the density in the
random effects, which the plain value does not. The whole `frm()` call:
1.09 to 1.36 s against 1.90 to 2.38 s.

## 4. The fit at the poisson limit

The repro table in section 1: on the lane, every rounding of both data
sets on both BLAS stops at code 0 at one point (nb200 log shape 17.2212,
logLik -386.278985229396; nb300 log shape 18.9741 for negbinomial and
-18.9741 for negbinomial2, logLik -555.262460892335 for both), with the
limit message under "warning" (126 of 126) and nothing under "ignore"
(0 warnings, 0 messages on 126). The noise had also pushed the two
spellings apart; they now stop at mirrored points.

5224 (code 1): the stable density removes the code-1 stops (0 of 126,
against 26), so the question of the brief's step 4, whether the SE
check's limit rule should read a code-1 stop on that ridge, does not
arise for this case; on 0.83.0 the rule already reads them (every
code-1 stop got the message under "warning"). No change to
`R/se-check.R`.

## 5. Tests

- `tests/testthat/test-numerical-robustness.R`:
  - new "negative-binomial log masses survive a large or small shape":
    `expect_robust_at()` sweeps the shape at linear predictor -30, 30
    and 40 for negbinomial, zero_inflated_negbinomial,
    hurdle_negbinomial, nbinom1 (`phi`, -e0) and negbinomial2 (`sigma`,
    -e0), counts 0, 1, 3, 40 (40 takes the `lgamma_shift_diff()`
    route). The comment that said the shape is not swept is replaced.
  - new "the negative-binomial gradient in the shape reads the poisson
    limit": item 4903's data, negbinomial and zero_inflated_negbinomial,
    the tape's derivative at log shape 20, 24.479, 30, 35, 40 against
    the exact derivative in a cancellation-free form (series in
    `mu / shape`), bound 64 `eps * sum(y + mu + 1)`. Margin
    (`dev/nbrobust-gradmargin.R`, 21 ulp perturbations of eta, 210
    points per build): worst 2.882 units of 64 on the reference BLAS
    and 2.882 on OpenBLAS 0.3.26 Sandybridge.
  - Seen to fail on rellib-r22 (`dev/nbrobust-tests/base-*.txt`):
    `RESULT test-numerical-robustness.R pass=758 fail=20`: the sweep at
    30 (`0.00189 >= 0.00001`, 5 families) and 40 (`47.00000 >=
    0.00001`, 5 families); the gradient test at every log shape, e.g.
    `0.18417777626 >= 0.00000000002` at 30, `804` at 40. Lane:
    `pass=778 fail=0 err=0 skip=0 warn=0`.
- `tests/testthat/test-family-fitend.R`:
  - new "a negbinomial fit started far along its poisson limit stops
    there": item 4903's data from `start = list(betad = 30)`, the limit
    message and no warning; under "ignore", no warning. Seen to fail on
    rellib-r22: `pass=51 fail=3` (no message; "false convergence (8) ...
    would still gain 0.105" under both settings). Lane `pass=54 fail=0`.
  - "negbinomial and negbinomial2 both name the poisson limit": the
    branch that accepted a code-1 stop with the convergence warnings and
    no message is removed; both spellings must give the message and no
    warning. It passes on rellib-r22 as well (the 0.83.0 limit rule
    reads code-1 stops), so it is a tightening, not a pin; the lane's 84
    nb300 fits on both BLAS all end so.

(Suite, OpenBLAS and check results: section 6.)

## 6. Suites, OpenBLAS, check

Every suite, reference BLAS, one file per process, two at a time
(`NOT_CRAN=true FRMTMB_RELLIB=<lane lib> bash dev/release/run-par.sh
nbrobust-core dev/nbrobust-core.jobs 2`, the same with
`dev/nbrobust-ext.jobs`; logs `dev/release/nbrobust-{core,ext}.log`,
files `dev/release/nbrobust-{core,ext}-files/`). Every per-file log's
`lib:` line names `wt-nbrobust-lib` for frmtmb. Summed from the RESULT
lines:

    NBROBUST-CORE ran 349 of 349: pass=25706 fail=0 err=0 skip=243 warn=0
    NBROBUST-EXT ran 200 of 200:
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1851 fail=0 err=0 skip=5 warn=0
    frmtmb.latent files=13 pass=700 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=92 pass=4286 fail=0 err=0 skip=18 warn=0
    frmtmb.spatial files=8 pass=551 fail=0 err=0 skip=6 warn=0
    frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0

Gated (`FRMTMB_BRMS_FIT_TESTS=true`), every test file of every suite
that names negbinomial, nbinom1 or geometric (`dev/nbrobust-gated.jobs`,
43 files; log `dev/release/nbrobust-gated.log`):

    NBROBUST-GATED ran 43 of 43: pass=8039 fail=0 err=0 skip=13 warn=0

The 13 skips are all of `test-drmtmb-agreement.R`, behind its own gate;
with `FRMTMB_DRMTMB_FIT_TESTS=true` as well (it has an NB2 agreement
test; `dev/nbrobust-tests/lib-drmtmb-gated.txt`):
`RESULT test-drmtmb-agreement.R pass=131 fail=0 err=0 skip=0 warn=0`.

OpenBLAS 0.3.26 with LAPACK, one thread (`dev/nbrobust-tests/ob326/`;
the first four with `OPENBLAS_CORETYPE=Sandybridge`):

    RESULT test-numerical-robustness.R pass=778 fail=0 err=0 skip=0 warn=0
    RESULT test-family-fitend.R pass=54 fail=0 err=0 skip=0 warn=0
    RESULT test-se-limit.R pass=379 fail=0 err=0 skip=0 warn=0
    RESULT test-families.R pass=252 fail=0 err=0 skip=0 warn=0
    RESULT test-negbinomial2.R pass=69 fail=0 err=0 skip=0 warn=0
    RESULT test-zi.R pass=11 fail=0 err=0 skip=0 warn=0
    RESULT test-hurdle-nb-zoib.R pass=89 fail=0 err=0 skip=0 warn=0

The fit-level perturbation runs (section 1, 21 roundings each, both
BLAS) and the gradient margin (section 5) are the 20-plus-draw evidence
for the new assertions. OpenBLAS 0.3.32 was not run: with the exact
derivative the poisson-limit fits end at one point on every rounding,
so no new assertion sits near a rounding-decided branch.

Every suite in one session, three orders (`NOT_CRAN=true
FRMTMB_RELLIB=<lane lib> bash dev/release/run-sessions.sh nbrobust-sess
2`, log `dev/release/nbrobust-sess.log`): NBROBUST-SESS ran 27 of 27,
every session fail=0 err=0 warn=0; all 27 load frmtmb from
`wt-nbrobust-lib` (extensions from rellib-r22, unchanged). No untagged
LEAK line: every one carries "[with load of ...]".

    SESSION frmtmb forward files=336/349 pass=25706 fail=0 err=0 skip=243 warn=0
    SESSION frmtmb reverse files=336/349 pass=25871 fail=0 err=0 skip=243 warn=0
    SESSION frmtmb preload files=336/349 pass=25975 fail=0 err=0 skip=245 warn=0
    SESSION frmtmb.coupling {forward,reverse,preload} files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam {forward,reverse,preload} files=33/33 pass=1851 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.latent {forward,reverse,preload} files=13/13 pass=700 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn {forward,reverse,preload} files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.ode {forward,reverse,preload} files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample {forward,reverse,preload} files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0
    SESSION frmtmb.survival {forward,reverse,preload} files=16/16 pass=575 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial {forward,reverse,preload} files=8/8 pass=551 fail=0 err=0 skip=6 warn=0

(The 13 core and 5 sample files without a row are the gated files whose
file-level skip leaves none, as in `dev/ci082-findings.md`.)

`R CMD check --as-cran` of core (`dev/nbrobust-check.sh`, built and
checked in `dev/nbrobust-check/`, `R_LIBS` = lane lib, rellib-r22, user
library; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; pandoc and TinyTeX on
PATH; log kept as `dev/nbrobust-out/check.log`, the check tree deleted):
`Status: 1 NOTE`, the expected "Skipping checking math rendering:
package 'V8' unavailable". No extension changed, so no extension check.

## Bump and floors

Minor (a BREAKING change in every negative-binomial fit's stop at the
poisson limit, and in the value to about 1e-12). No dependency floor
changes: `nb_lpmf()` uses `RTMB::logspace_add()`, `log1p()`, `cumsum()`
and `lgamma_shift_diff()`, all used by 0.83.0 already.

## Files touched

- `R/families.R`: `nb_lpmf()`, `nb_rising()` (new, after
  `beta_binomial_lpmf()`); the `lpdf` slots of `fam_negbinomial()`,
  `fam_nbinom1()`, `fam_zi_negbinomial()`, `fam_hurdle_negbinomial()`,
  `fam_geometric()`. No `lcdf`/`cdf_log` line touched (lane cdfname).
- `R/links.R`: two comments that named `dnbinom_robust()` (the roxygen
  of the robust-link section and the softplus `log_eta` comment).
- `man/frmtmb-links.Rd`: regenerated from that roxygen.
- `NEWS.md`: one bullet under `# frmtmb (development version)`,
  `## Families`.
- `tests/testthat/test-numerical-robustness.R`,
  `tests/testthat/test-family-fitend.R`.
- `dev/upstream-bugs.md` (RTMB-4), `dev/test-backlog.md` (two CLOSED
  lines), this file, and the scripts `dev/nbrobust-*.R`,
  `dev/nbrobust-repro-all.sh`, `dev/nbrobust-check.sh`,
  `dev/nbrobust-{core,ext,gated}.jobs`.

## Not done, and why

- The negative-binomial CDF (`nbinom_cdf_log()`, `negbinomial_cdf_pair()`
  and the censored and truncated rows) is not changed: the brief is the
  density, and the CDF slots belong to lane cdfname this round. Whether
  its derivative in the shape has the same large-shape noise is not
  measured.
- brms was not run: the density's value is unchanged to 5e-12 relative,
  so no brms comparison moves.

## Fix round 1 (review `dev/reviews/2026-11-17-nbrobust.md`, verdict FIX)

The review found the code right; its fixes are records. No package code
changed in this round, so nothing was rerun.

- **What the lane takes away and does not replace.** On 0.83.0 the
  noise of `dnbinom_robust()`'s shape derivative gave, by accident, a
  "Large maximum absolute gradient" warning on 5 of 6 fits started at
  log shape 25 or 30 on overdispersed data (seed 5, n = 200,
  `rnbinom(size = 2)`, maximum at log shape 0.7613, logLik
  -424.5690959863). Those fits stop at their start, 46 log-likelihood
  units below the maximum, on both builds. With the exact derivative
  (6e-9 at 25, 4e-11 at 30) no gradient warning comes, and the only
  signal left is the SE check's "flat ... the data do not determine
  it", which is false there. This lane removes an accidental alarm and
  nothing replaces it. Filed MEDIUM in `dev/test-backlog.md` ("Filed in
  the round of 2026-11-16", "From lane nbrobust") with a candidate fix
  (the flat verdict testing the objective at the data's moment estimate
  of the shape), not built here; said in NEWS as well.
- Filed LOW: `nbinom_cdf_log()` is NaN from log shape about 37.5 for
  counts of 1 or more (review finding 4), so a censored fit started at
  `betad = 40` errors with "NA/NaN gradient evaluation". LOW because
  only a far start reaches it and it fails loudly; owner
  `R/families-cdf.R` (cdfname's area this round, not edited).
- Filed LOW: `nb_rising()` stops on an NA count (review finding 9.4);
  no package path reaches it, so no guard was added in this round.
- NEWS: a second bullet with the value fixes the review measured (zinb
  under a non-log link off by up to 0.061 on the tape; off the tape,
  the values `frmtmb.sample::log_lik()`/`loo()` read, NaN at log shapes
  35.7 to 39.9 for zinb log/sqrt/softplus/squareplus, negbinomial
  softplus with `rate()` and nbinom1 sqrt, and zinb off by up to 0.899
  on every link), and a paragraph on the far start above.
- `dev/upstream-bugs.md`: new section "R (stats)", R-1,
  `stats::dnbinom(mu = )` losing digits at a large size (9.8e-8
  relative at log size 22; 0.184 absolute at log mean 11.5, log size
  24), which brms's `log_lik()` inherits; not reported.

Files touched in this round: `NEWS.md`, `dev/test-backlog.md`,
`dev/upstream-bugs.md`, this file.
