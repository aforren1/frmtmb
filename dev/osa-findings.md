# Lane osa: one-step-ahead residuals on discrete responses

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1), worktree
`frmtmb-wt-osa`, library `wt-osa-lib`. "Before" is `rellib-r18`;
"after" is this worktree installed into `wt-osa-lib`. RTMB 2.0, TMB
1.9.25 and RTMBdist 1.2.0 from the user library in both arms.

Every block marked "generated" is filled from its log by
`dev/osa-findings-gen.R`; the logs are under `dev/osa-out/` (local,
gitignored). Seeds and scripts are named beside each block.

## 0. Fix round after the review (`dev/reviews/2026-11-06-osa.md`)

Sections 1 to 8 are the first pass. Where this section and they
disagree, this section holds. The reviewer's scripts (`dev/review-osa-*`)
were rerun on the fix build; their logs here are `dev/osa-out/fix-*.log`.

### 0.1 B1: weights under a group-level term

The density term keeps the weight `w`, and each tail gets `(w - 1)`
copies of the density (`osa_cdf_loglik()`), the reviewer's form. TMB's
row mass is then `p` and its CDF `F`: the residual of the last of `w`
copies, with the later rows conditioned on all `w`. A row of weight 0 is
`NA`, and so is a row of weight between 0 and 1, with a warning. The
tail factor is clamped at `pmax(w, 1) - 1`: with `w - 1` below 1 the
tails were `S(y | u) / p(y | u)^(1 - w)`, which a group effect's prior
does not integrate (a weight-0 row read `nlcdf.upper` NaN), and with the
clamp they are finite. The reason for `NA` below 1 is not the integral
(the first version of this paragraph said so; that held only without the
clamp). With the clamp, `Fx` is the one-step CDF of a fresh copy of the
row given the earlier rows at their weights, but TMB's row mass `px`
comes from the density term `w ll`, so it is the mean of `p^w`, not of
`p`, and a count's randomization over it is wrong. For a continuous
response `U = 0`, so the fresh-copy residual `qnorm(Fx)` would be sound.
Open option, for the user: return that fresh-copy residual for a
continuous response with a weight in (0, 1), as the reviewer suggested
(`dev/reviews/2026-11-06-osa.md`, short check, item 1); for now every
family gives `NA` there. On those rows TMB's own `qnorm(Fx - U px)` warned
"NaNs produced"; that one warning is muffled around the call when such
rows exist, and only the classed `frmtmb_osa_weight` warning reaches the
user (`test-osa-cdf.R`). `dev/review-osa-weights.R` on the
fix build (seed 301, weight 2 on rows t <= 3 of 20 groups of 10, against
the rows written twice; seed 302 for gaussian). The 8.27e-07 is poisson's:
the short check measured a weighted gaussian with `(1 | g)` 0.0018 from
the rows written out, which keep the exact `fullGaussian`, the Laplace
error of the `"cdf"` route (`dev/review-osa-out/sc-wgauss.log`); NEWS
now says so.

<!-- gen:fix-weights -->
```
BEGIN review-osa-weights
estimates weighted  0.15940010  0.41420402 -0.00994222 | duplicated  0.1593995  0.4142042 -0.0099423
rows t >= 4 (weight 1, after the weighted rows): max|qnorm(hi) weighted - duplicated| 8.27e-07, median 1.26e-07
weighted fit, rows t >= 4: max|lane - exact weighted| 0.00236, max|lane - exact unweighted| 0.33, max|exact weighted - exact unweighted| 0.328 (qnorm of the upper end)
weighted fit residual vs exact weighted-conditioning residual, rows t >= 4: max 0.00244 median 0.000267; rows over 0.05: 0 of 140
gaussian weight 3 rows: mean r_fullGaussian / r_cdf 1.0000 (sqrt(3) = 1.7321); weight 1 rows: 1.0000. student() same data, sd of weight-3 rows 1.023 vs gaussian default 1.027
END review-osa-weights
```
<!-- /gen -->

The reviewer's m5 read 8.27e-07 and 0.00244 on the same lines; this
build reads the same. Weight 0 (`dev/review-osa-edge.R`, seed 407): the
weight-0 rows are `NA` and the script's mean of the other rows reads
`NA` for that reason; `test-osa-cdf.R` checks the other rows against the
fit without them instead (section 0.5, m6).

m1: a weighted `gaussian()` now takes `"cdf"`, so its weight is a
frequency, as on every other family and in brms. The gaussian line of
the block above reads `mean r_fullGaussian / r_cdf 1.0000` on the
weight-3 rows (it was 1.7321, `sqrt(3)`), and student() and gaussian()
on the same data agree (sd 1.023 against 1.027). `fullGaussian`, asked
for by name, still reads a weight as a precision; the help says so.

### 0.2 M4 and m5: a residual from the smaller tail

`osa_tail_residual()` forms the residual from the `nll`, `nlcdf.lower`
and `nlcdf.upper` columns `oneStepPredict()` returns: `qnorm(log(Fx - U
px), log.p = TRUE)` where `Fx <= 1/2`, and `-qnorm(log(Sx + U px))`
above, with TMB's own uniforms drawn again (`osa_uniforms()`: the same
`set.seed(seed)`, or the stream as it stood before the call for
`seed = NULL`; the session stream is left as TMB left it). The tails are
floored at -1e4 instead of `log(1e-300)`. `dev/review-osa-edge.R` and
`-edge2.R` on the fix build:

<!-- gen:fix-edge -->
```
BEGIN review-osa-edge
poisson outliers y = 14, 22, 60              row  10 r     4.5577  exact step [    4.4368,     4.7696]
poisson outliers y = 14, 22, 60              row  20 r     6.8031  exact step [    6.7974,     7.0772]
poisson outliers y = 14, 22, 60              row  30 r    15.2428  exact step [       Inf,        Inf]  OUTSIDE
                                             upper tails S(y - 1): 4.56e-06 5.33e-12 4.71e-52
                                             true residual -qnorm(S(y-1) - (1-U) p): row 10 4.5577, row 20: 6.8031, row 30: 15.243
Gamma outlier y = 30                         row 5 r 6.256541; shape 2.828, pgamma upper tail 1.97e-10 -> true 6.257
Gamma outlier                                y 40: r 6.58663, true 6.587 (S 2.25e-11)
Gamma outlier                                y 50: r 6.77672, true 6.777 (S 6.15e-12)
Gamma outlier                                y 60: r 6.89314, true 6.893 (S 2.73e-12)
Gamma outlier                                y 80: r 7.01461, true 7.015 (S 1.15e-12)
                                             WARN 
gaussian outlier at about 9.5 sd             row 7: fullGaussian 7.537, trunc(lb = -50) via cdf 7.537399
poisson y = 0 at large mu                    mu 702.0, log F(0) = -702.0; r -37.36237; true in [-Inf, -37.35]
binomial support edges                       rows at 0 or 6: 48, outside their step: 0, nonfinite 0
poisson trunc(lb = 3, ub = 7)                rows 92 (at 3: 18, at 7: 16), outside step 0, nonfinite 0
Warning message:
Optimizer did not report convergence: false convergence (8). A Newton step from the stop gains 1.95e-13 in log-likelihood, within frmtmb_control(grad_tol = 0.001), so the estimates are at the maximum and their standard errors are judged as for a converged fit 
poisson trunc(lb = 3), mu about 0.05         window mass 1.2e-42; outside step 0 of 80, nonfinite 0, y range 3-3
weight 0 on rows t <= 4 (y + 6)              estimates equal: TRUE; rows t > 4: mean r weighted fit -0.005, mean r fit without the rows 0.014
END review-osa-edge
```
<!-- /gen -->

<!-- gen:fix-edge2 -->
```
BEGIN review-osa-edge2
outlier  +7.0: fullGaussian   6.9385 | trunc() via cdf  6.93855 | (y - mu) / sigma   6.9385
outlier  +8.0: fullGaussian   7.8973 | trunc() via cdf   7.8973 | (y - mu) / sigma   7.8973
outlier  +8.5: fullGaussian   8.3726 | trunc() via cdf  8.37258 | (y - mu) / sigma   8.3726
outlier  +9.0: fullGaussian   8.8449 | trunc() via cdf  8.84494 | (y - mu) / sigma   8.8449
outlier +12.0: fullGaussian  11.6103 | trunc() via cdf  11.6103 | (y - mu) / sigma  11.6103
outlier  -9.0: fullGaussian  -8.7980 | trunc() via cdf -8.79797 | (y - mu) / sigma  -8.7980
outlier -12.0: fullGaussian -11.5654 | trunc() via cdf -11.5654 | (y - mu) / sigma -11.5654
END review-osa-edge2
```
<!-- /gen -->

Poisson y = 60 at a mean of about 3 reads 15.2428 (true 15.243; the
script's "OUTSIDE" is its own reference, which is `Inf` in double
precision there). Poisson y = 0 at mu 702 reads -37.362, inside its step
(at most -37.35), where the floor clamped it at -37.05. The gaussian
outliers at +-12 read the exact `(y - mu) / sigma`.

### 0.3 M3: the cost

First, the closed forms the coordinator named, on the one-step tape
(`dev/osa-tapecost.R`, 500 rows, 300 gradient-and-Hessian evaluations of
a tape that carries the density and both tails at every row, seed 1):

<!-- gen:tapecost -->
```
BEGIN osa-tapecost
pois_density   300 gradients and Hessians   0.030s
pois_summed    300 gradients and Hessians   4.580s
pois_pgamma    300 gradients and Hessians   4.920s
binom_density  300 gradients and Hessians   0.050s
binom_ibeta    300 gradients and Hessians   4.390s
binom_pbeta    300 gradients and Hessians   3.360s
bern_ibeta     300 gradients and Hessians   3.360s
bern_closed    300 gradients and Hessians   0.420s
nb_density     300 gradients and Hessians   0.470s
nb_ibeta       300 gradients and Hessians   5.940s
nb_pbeta       300 gradients and Hessians   4.880s
END osa-tapecost
```
<!-- /gen -->

`RTMB::pgamma()` (poisson) costs more than the summed pair on the tape,
and `RTMB::pbeta()` less. Inside whole `residuals()` calls, poisson with
`pgamma()` took 0.63 of the summed pair's time, negbinomial with
`pbeta()` 0.67, binomial with `pbeta()` 1.15, and bernoulli in closed form
0.45 (500 rows, `(1 | g)`, in one process each; `dev/osa-time2.R` at
that stage, logs `time2b-*-25.log`). None of them came near 2 times the
old route for poisson (4.31 with `pgamma()`, 6.80 summed).

The cost was the tails of every row on every evaluation: the indicator
of a tail is a parameter of TMB's tape, so nothing is skipped, and TMB
makes three evaluations per row. `residual_values()` now hands the rows
to `oneStepPredict()` in blocks of about `sqrt(n)` (at least 20): block
`a..b` with rows `1..a-1` as `conditional` and the rows after it switched
off, which is the same one-step predictive, and the tape of that call
carries the tails of the block alone (`osa_cdf_blocks()`,
`osa_cdf_loglik()`, option `frmtmb.osa_tail_rows`, set and restored
around each block). With blocks the summed pairs and the closed forms
were within 15 percent of each other (poisson 7.5 s against 6.6 s,
negbinomial 13.6 s against 11.8 s, `time2c-*-25.log`), so the summed pairs
stay: they keep every digit in both tails, which `pgamma()`'s complement
did not on poisson's lower tail. Only bernoulli takes a closed form
(`osa_bernoulli_pair()`), which is exact.

`dev/osa-time2.R` on the final route: "cdf" is the default (blocks),
"old" the 0.79.1 call on the same objective, "one block" the default
route over every row at once (option `frmtmb.osa_block = Inf`), "control"
the default again. Seed 601, groups of 20. Generated:

<!-- gen:time3-poisson-25 -->
```
BEGIN osa-time2
poisson      rows  500 (1 | g): cdf     8.6s old    28.1s one block   144.0s control     8.6s | cdf/old 0.31 one block/old 5.13 control/cdf 1.00
END osa-time2
```
<!-- /gen -->
<!-- gen:time3-negbinomial-25 -->
```
BEGIN osa-time2
negbinomial  rows  500 (1 | g): cdf    13.9s old   126.7s one block   136.8s control    12.2s | cdf/old 0.11 one block/old 1.08 control/cdf 0.88
END osa-time2
```
<!-- /gen -->
<!-- gen:time3-bernoulli-25 -->
```
BEGIN osa-time2
bernoulli    rows  500 (1 | g): cdf     3.2s old    57.8s one block     7.3s control     3.3s | cdf/old 0.06 one block/old 0.13 control/cdf 1.03
END osa-time2
```
<!-- /gen -->
<!-- gen:time3-poisson-100 -->
```
BEGIN osa-time2
poisson      rows 2000 (1 | g): cdf    59.1s old   303.2s one block      NAs control    53.9s | cdf/old 0.20 one block/old NA control/cdf 0.91
END osa-time2
```
<!-- /gen -->
<!-- gen:time3-negbinomial-100 -->
```
BEGIN osa-time2
negbinomial  rows 2000 (1 | g): cdf   117.7s old  1800.9s one block      NAs control   101.3s | cdf/old 0.07 one block/old NA control/cdf 0.86
END osa-time2
```
<!-- /gen -->
<!-- gen:time3-bernoulli-100 -->
```
BEGIN osa-time2
bernoulli    rows 2000 (1 | g): cdf    38.0s old   771.1s one block      NAs control    34.0s | cdf/old 0.05 one block/old NA control/cdf 0.89
END osa-time2
```
<!-- /gen -->

The uniforms of a block come from one draw over all rows, so the
residuals do not depend on the blocks: `test-osa-cdf.R` compares every
count family with its exact step on 150 rows, which run in 8 blocks.

### 0.4 M2: the Laplace error, stated

NEWS and `?residuals.frmtmb_fit` now give the reviewer's measurement
(`dev/review-osa-out/exact-summary.txt`, seeds 201 to 203, 20 groups of
10): a median error below 0.001 under `(1 | g)` on seven count families,
0.006 to 0.016 under a latent `ar()` on poisson and the zero-inflated
and hurdle counts, a largest error of 0.15 under `(1 | g)` and 0.13 under
`ar()`, carried by the first row of each group, and the KS p-value within
0.13 of the exact residuals' on 42 fits.

### 0.5 The rest

- M1: `residuals_osa x von_mises` is `works`, with a test. Every compat
  row this lane changed was reread against its note (section 5 has the
  list); the note of `group:autocor x residuals_osa` said bernoulli was
  broken and the route exact, and now gives the Laplace error instead;
  `cens()`/`trunc()` x bernoulli no longer say "The family has no CDF";
  `cens()` and `trunc()` x nbinom1 have rows of their own, which said
  "needs a family with an AD log-CDF" through `kind:family`.
- m2: a family with a CDF that leaves out `cens()` or `trunc()`
  (bernoulli, nbinom1) is refused with "is not supported for this
  family" instead of the generic "would change nothing about the fit"
  (`check_accepted_aterms()`, `R/families.R`).
- m3: NEWS and the help say `u` is uniform on `[F(y - 1), F(y)]`, give
  1.2e-12 for the count families, say bernoulli was wrong with no group
  term too, mark the new refusals BREAKING, and document `seed = 123`
  and `seed =`. The "on bernoulli() they are faster" line is gone.
- m4: a `quadrature = TRUE` fit is refused up front, class
  `frmtmb_osa_quadrature`, with the reason.
- m6: new tests: weight 2 against the duplicated rows and weight 0
  against the dropped rows under `(1 | g)`; the truncated gaussian
  `(1 | g)` test now compares with the exact predictive (ratio 0.01);
  step checks for zero_inflated_negbinomial, hurdle_negbinomial,
  zero_inflated_binomial and zero_inflated_beta_binomial; the poisson
  y = 60 and gaussian +-12 outliers; the RNG state; the quadrature
  refusal; the von_mises row.
- m7: TMB-6 in `dev/upstream-bugs.md`, and RTMBdist-5 points to it.

### 0.6 The tests, on four builds

`tests/testthat/test-osa-cdf.R` (final), one process each, on the fix
build, on rellib-r18, and on the reviewer's `rev-osa-lib/m4` (the first
pass with the `trunc()` window ignored) and `m5` (the first pass with
B1's two lines). Every new test fails on a build without its fix: the
weight-2 test on m4 and passes on m5, which is B1's fix; the truncated
gaussian `(1 | g)` test now fails on m4; the outlier, weight-0 and
quadrature tests fail on m5. Two tests pass on rellib-r18 as well: under
the 0.79.1 route a weight scaled the density in the conditioning, as a
duplicated row does, and weight-0 rows read NaN there. Generated:

Fix build:

<!-- gen:t0-fix -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-osa-cdf.R pass=58 fail=0 err=0 skip=0 warn=0
```
<!-- /gen -->

rellib-r18:

<!-- gen:t0-before2 -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18; frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18
RESULT test-osa-cdf.R pass=17 fail=19 err=7 skip=0 warn=0
      1 ERROR [a quadrature fit and the compat rows say what is true]
      1 ERROR [a truncated count is read inside its window, row by row]
      1 ERROR [a truncated gaussian with a group-level term has every row]
      1 ERROR [a zero-inflated or hurdle count reads its point mass as data]
      1 ERROR [continuous residuals are qnorm of the fitted CDF]
      1 ERROR [count residuals are randomized quantile residuals of the fit]
      1 ERROR [the zero-inflated and hurdle counts sit in their own steps]
      1 FAIL [a case weight does not reshape a row's own step]
      1 FAIL [a non-finite residual is reported, not passed on silently]
      2 FAIL [an outlier's residual is read in its own tail]
      6 FAIL [bernoulli residuals are calibrated under group-level terms]
      2 FAIL [continuous residuals are qnorm of the fitted CDF]
      4 FAIL [count residuals are randomized quantile residuals of the fit]
      3 FAIL [the routes that cannot read a count refuse it by name]
```
<!-- /gen -->

m4:

<!-- gen:t0-m4 -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/rev-osa-lib/m4; frmtmb 0.79.1 C:/Users/adf44/source/r/rev-osa-lib/m4
RESULT test-osa-cdf.R pass=46 fail=9 err=1 skip=0 warn=1
      1 ERROR [a quadrature fit and the compat rows say what is true]
      2 FAIL [a row of weight 0 is NA and leaves the other rows alone]
      2 FAIL [a truncated count is read inside its window, row by row]
      1 FAIL [a truncated gaussian with a group-level term has every row]
      1 FAIL [a weight is a frequency under a group-level term]
      3 FAIL [an outlier's residual is read in its own tail]
      1 WARN [an outlier's residual is read in its own tail]
```
<!-- /gen -->

m5:

<!-- gen:t0-m5 -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/rev-osa-lib/m5; frmtmb 0.79.1 C:/Users/adf44/source/r/rev-osa-lib/m5
RESULT test-osa-cdf.R pass=51 fail=4 err=1 skip=0 warn=4
      1 ERROR [a quadrature fit and the compat rows say what is true]
      1 FAIL [a row of weight 0 is NA and leaves the other rows alone]
      3 FAIL [an outlier's residual is read in its own tail]
      3 WARN [a row of weight 0 is NA and leaves the other rows alone]
      1 WARN [an outlier's residual is read in its own tail]
```
<!-- /gen -->

The other changed files, and `test-compat.R` and `test-arg-refusal.R`,
one per process on the fix build (`dev/osa-out/t-fix`):

<!-- gen:t-fix -->
```
RESULT test-osa-cdf.R pass=58 fail=0 err=0 skip=0 warn=0 | lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-beta-binomial-density.R pass=96 fail=0 err=0 skip=0 warn=0 | lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-hurdle-nb-zoib.R pass=89 fail=0 err=0 skip=0 warn=0 | lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-xbeta-zibb-hurdle-cum.R pass=606 fail=0 err=0 skip=0 warn=0 | lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-cdfs.R pass=234 fail=0 err=0 skip=0 warn=0 | lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-cens-trunc.R pass=72 fail=0 err=0 skip=0 warn=0 | lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-compat.R pass=759 fail=0 err=0 skip=0 warn=0 | lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-arg-refusal.R pass=122 fail=0 err=0 skip=0 warn=0 | lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
```
<!-- /gen -->

The whole suites on the fix build:

<!-- gen:suites-fix -->
```
BEGIN suites (fix build)
SKIPS frmtmb RESULT test-bcm-binomial.R pass=17 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb RESULT test-bcm-esp.R pass=5 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-bcm-gcm.R pass=9 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-bcm-model-selection.R pass=13 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-bcm-psychophysics.R pass=7 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-bcm-signal-detection.R pass=6 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-brms-agreement.R pass=172 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-brms-lik-cdfs.R pass=0 fail=0 err=0 skip=5 warn=0
SKIPS frmtmb RESULT test-brms-likelihood.R pass=5 fail=0 err=0 skip=55 warn=0
SKIPS frmtmb RESULT test-ce-aterms.R pass=19 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-curve-scale.R pass=0 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-drmtmb-agreement.R pass=0 fail=0 err=0 skip=13 warn=0
SKIPS frmtmb RESULT test-emmeans.R pass=92 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-fuzz.R pass=0 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-gp-by.R pass=71 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-latres-brms.R pass=30 fail=0 err=0 skip=5 warn=0
SKIPS frmtmb RESULT test-mi-mm-brms-l1.R pass=0 fail=0 err=0 skip=10 warn=0
SKIPS frmtmb RESULT test-offset-grid.R pass=27 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-refcat-brms.R pass=0 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb RESULT test-bcm-bart.R pass=8 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-bcm-data-analysis.R pass=20 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb RESULT test-bcm-gaussian.R pass=11 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-bcm-latent-mixtures.R pass=30 fail=0 err=0 skip=5 warn=0
SKIPS frmtmb RESULT test-bcm-mpt.R pass=15 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-bcm-retention.R pass=15 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-bcm-simple.R pass=8 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-brms-lik-gp-exact.R pass=5 fail=0 err=0 skip=4 warn=0
SKIPS frmtmb RESULT test-brms-methods.R pass=18 fail=0 err=0 skip=47 warn=0
SKIPS frmtmb RESULT test-brms-priors.R pass=0 fail=0 err=0 skip=12 warn=0
SKIPS frmtmb RESULT test-ce-predict-median.R pass=11 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-fcor.R pass=50 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-mo-factor-re.R pass=34 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-nb2dm-brms.R pass=0 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb RESULT test-nl-brms.R pass=0 fail=0 err=0 skip=10 warn=0
SKIPS frmtmb RESULT test-nl-stan-funs.R pass=122 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-ordinal-mixture-brms.R pass=0 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-rl-example.R pass=96 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-simplex-brms.R pass=0 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb.coupling RESULT test-scale.R pass=0 fail=0 err=0 skip=5 warn=0
SKIPS frmtmb.eam RESULT test-sample-ndt-default.R pass=7 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.eam RESULT test-scale.R pass=0 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb.latent RESULT test-scale.R pass=0 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.learn RESULT test-scale.R pass=0 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.learn RESULT test-stan-identity.R pass=0 fail=0 err=0 skip=11 warn=0
SKIPS frmtmb.ode RESULT test-scale.R pass=0 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb.sample RESULT test-default-priors-fit.R pass=0 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb.sample RESULT test-loglik-autocor.R pass=110 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.sample RESULT test-loo-predict.R pass=114 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb.sample RESULT test-loo.R pass=138 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.sample RESULT test-mi-draws-brms.R pass=0 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb.sample RESULT test-mo-coding.R pass=9 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.sample RESULT test-sampling-ported.R pass=207 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb.sample RESULT test-scale.R pass=0 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb.spatial RESULT test-brms-likelihood-car.R pass=4 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.spatial RESULT test-sar.R pass=106 fail=0 err=0 skip=4 warn=0
frmtmb.spatial     files    8 pass    478 fail 0 err 0 skip 6 warn 0
frmtmb.eam         files   33 pass   1848 fail 0 err 0 skip 5 warn 0
frmtmb.sample      files   85 pass   3962 fail 0 err 0 skip 18 warn 0
frmtmb.coupling    files   11 pass    551 fail 0 err 0 skip 5 warn 0
frmtmb             files  309 pass  22574 fail 0 err 0 skip 239 warn 0
frmtmb.survival    files   15 pass    552 fail 0 err 0 skip 0 warn 0
frmtmb.ode         files   11 pass    549 fail 0 err 0 skip 1 warn 0
frmtmb.latent      files   12 pass    683 fail 0 err 0 skip 2 warn 0
frmtmb.learn       files   16 pass    474 fail 0 err 0 skip 13 warn 0
files with no RESULT line: 0
frmtmb loaded from (files per lib line):
   500 frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
END suites
```
<!-- /gen -->

<!-- gen:check-fix -->
```
BEGIN check (fix build)
* checking whether package 'frmtmb' can be installed ... [28s] OK
* checking R code for possible problems ... [57s] OK
* checking examples ... [46s] OK
* checking examples with --run-donttest ... [48s] OK
* checking tests ...
  Running 'testthat.R' [29m]
 [29m] OK
* checking re-building of vignette outputs ... [163s] OK
* checking HTML version of manual ... [28s] NOTE
Skipping checking math rendering: package 'V8' unavailable

Status: 1 NOTE

END check
```
<!-- /gen -->

### 0.7 After the short check

The short check (end of `dev/reviews/2026-11-06-osa.md`) read SHIP with
small items: TMB's "NaNs produced" on fractional-weight count rows is
muffled; NEWS gives 8e-7 for poisson and 0.0018 for gaussian `(1 | g)`;
the fractional-weight rationale is rewritten (0.1); `test-osa-cdf.R`
pins the `set.seed(123)` uniforms directly against
`qnorm(ppois(y) - U dpois(y))`, and `seed = 7` likewise; the tapecost
block in 0.3 is regenerated from the corrected script (it ended in an R
error on the negbinomial arms before). `test-osa-cdf.R` on the final
build and on the reviewer's n2 mutant (`seed + 1`), generated:

<!-- gen:t0-final -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-osa-cdf.R pass=63 fail=0 err=0 skip=0 warn=0
```
<!-- /gen -->

<!-- gen:t0-n2 -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/rev-osa-lib/n2; frmtmb 0.79.1 C:/Users/adf44/source/r/rev-osa-lib/n2
RESULT test-osa-cdf.R pass=58 fail=5 err=0 skip=0 warn=6
      1 FAIL [a row of weight 0 is NA and leaves the other rows alone]
      1 FAIL [a weight between 0 and 1 gives NA and one classed warning]
      1 FAIL [a weight is a frequency under a group-level term]
      2 FAIL [the residuals leave the session's random stream alone]
      3 WARN [a row of weight 0 is NA and leaves the other rows alone]
      3 WARN [a weight is a frequency under a group-level term]
```
<!-- /gen -->

## 1. What was wrong

### 1.1 The cause, corrected

The backlog entry (`dev/test-backlog.md`, "Filed at the 0.79.1 release")
says the default `oneStepGeneric` integrates a discrete response "with no
`discrete = TRUE`". That is not what the code did: `residual_values()`
passed `discrete = TRUE` for every family of type `"discrete"`, with
`range = c(0, Inf)`. What it did not pass was `discreteSupport`. Without
one, TMB's `oneStepGeneric` takes its continuous branch even with
`discrete = TRUE`: it runs `tmbprofile()` on the observation with an
adaptive step, so it evaluates the density at NON-integer counts (the
mass function's continuous `lgamma()` extension), fits a spline through
that profile, and sums the spline on the integer grid of the profile's
range. Past a binomial's trials the extension is meaningless and the
spline reads NaN; elsewhere the sum is an approximation. Only the
ordinal families, which set `discreteSupport = 1:K`, took TMB's exact
discrete branch.

The beta-binomial NaN that 0.79.1 refused was the same route, not
RTMBdist 1.2.0. On RTMB and RTMBdist alone, with no frmtmb, the generic
call reads NaN on `dbetabinom()` and on `dbinom()` both; with
`discreteSupport = 0:10`, or with `method = "cdf"` (RTMBdist 1.2.0's new
`pbetabinom()`), it reads none. Generated, `dev/osa-rtmbdist.R`, data
seed 9, n = 150, 10 trials:

<!-- gen:rtmbdist -->
```
BEGIN osa-rtmbdist
dbetabinom  generic  nan 17 of 150
dbetabinom  support  nan 0 of 150
dbetabinom  cdf      nan 0 of 150
dbinom      generic  nan 2 of 150
dbinom      support  nan 0 of 150
dbinom      cdf      nan 0 of 150
END osa-rtmbdist
```
<!-- /gen -->

### 1.2 Every family, before and after

`dev/osa-audit.R`: every core family on one fixed design with no random
effect (data seed 9, n = 150; designs per support in the script). The
route is traced from the `oneStepPredict()` call itself. Truth is formed
from the family's DENSITY on doubles, never from its CDF slot or TMB: a
sum of the mass for a count, `stats::integrate()` for a continuous
response, and the ordinal families from `ord_probs()`. `outside` counts
rows outside `[qnorm(F(y - 1)), qnorm(F(y))]` (counts only); `maxdiff`
is the largest `|r - truth|`, where a count's truth uses TMB's own
uniforms (`set.seed(123); runif(n)`, `oneStepPredict()`'s default seed)
and a continuous one is `qnorm(F(y))`. `nan` is the number of non-finite
residuals, all silent before. `REFUSED` is an error; before, eight of
them were RTMB's S4 type errors rather than refusals.

Before (generated, `dev/osa-out/audit-before.log`):

<!-- gen:audit-before -->
```
BEGIN osa-audit before 
   family                      method         discrete support range   nan outside maxdiff  secs note                                                                                               
1  gaussian                    fullGaussian      NA    -       -        0  NA      6.66e-10 0.05                                                                                                    
2  student                     oneStepGeneric FALSE    -       -       NA  NA            NA 0.06 REFUSED: Not compatible with requested type: [type=S4; target=double].                             
3  skew_normal                 oneStepGeneric FALSE    -       -        0  NA      7.90e-06 1.33                                                                                                    
4  exgaussian                  oneStepGeneric FALSE    -       -       NA  NA            NA 0.00 REFUSED: Not compatible with requested type: [type=S4; target=double].                             
5  asym_laplace                oneStepGeneric FALSE    -       -       NA  NA            NA 0.02 REFUSED: Not compatible with requested type: [type=S4; target=double].                             
6  huber                       oneStepGeneric FALSE    -       -       NA  NA            NA 0.00 REFUSED: Not compatible with requested type: [type=S4; target=double].                             
7  gen_extreme_value           oneStepGeneric FALSE    -       -        0  NA      1.59e-04 1.37                                                                                                    
8  Gamma                       oneStepGeneric FALSE    -       -        0  NA      8.61e-02 0.84                                                                                                    
9  lognormal                   oneStepGeneric FALSE    -       -       NA  NA            NA 0.00 REFUSED: error in evaluating the argument 'x' in selecting a method for function 'dnorm': non-numer
10 inverse.gaussian            oneStepGeneric FALSE    -       -        3  NA      6.77e+00 0.99                                                                                                    
11 exponential                 oneStepGeneric FALSE    -       -       NA  NA            NA 0.01 REFUSED: Not compatible with requested type: [type=S4; target=double].                             
12 weibull                     oneStepGeneric FALSE    -       -       NA  NA            NA 1.03 REFUSED: $ operator is invalid for atomic vectors                                                  
13 shifted_lognormal           oneStepGeneric FALSE    -       -       NA  NA            NA 0.02 REFUSED: error in evaluating the argument 'x' in selecting a method for function 'dnorm': Not compa
14 frechet                     oneStepGeneric FALSE    -       -        3  NA      7.29e+00 1.20                                                                                                    
15 tweedie                     oneStepGeneric FALSE    -       -       71  NA      4.20e+00 4.17                                                                                                    
16 beta                        oneStepGeneric FALSE    -       -        0  NA      6.91e-01 1.19                                                                                                    
17 xbeta                       <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'xbeta': its density puts a point mass
18 zero_inflated_beta          <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_inflated_beta': its density puts
19 zero_one_inflated_beta      <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_one_inflated_beta': its density 
20 poisson                     oneStepGeneric  TRUE    -       [0,Inf]  0   5      3.35e-01 0.70                                                                                                    
21 negbinomial                 oneStepGeneric  TRUE    -       [0,Inf]  0  12      2.92e-01 0.71                                                                                                    
22 negbinomial2                oneStepGeneric  TRUE    -       [0,Inf]  0  12      2.92e-01 0.58                                                                                                    
23 nbinom1                     oneStepGeneric  TRUE    -       [0,Inf]  0  10      2.82e-01 0.60                                                                                                    
24 compois                     oneStepGeneric  TRUE    -       [0,Inf]  0  13      3.20e-01 0.60                                                                                                    
25 geometric                   oneStepGeneric  TRUE    -       [0,Inf]  0   0      3.36e-05 0.57                                                                                                    
26 discrete_weibull            oneStepGeneric  TRUE    -       [0,Inf]  0  12      3.60e-01 0.83                                                                                                    
27 binomial                    oneStepGeneric  TRUE    -       [0,Inf]  2   0      8.32e-03 0.91                                                                                                    
28 beta_binomial               oneStepGeneric  TRUE    -       [0,Inf] NA  NA            NA 0.00 REFUSED: beta_binomial() does not support OSA residuals. Use another `type` of residuals()         
29 bernoulli                   oneStepGeneric  TRUE    -       [0,Inf]  0  57      6.18e+00 0.75                                                                                                    
30 zero_inflated_poisson       <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_inflated_poisson': its density p
31 zero_inflated_negbinomial   <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_inflated_negbinomial': its densi
32 hurdle_poisson              <NA>              NA    <NA>    <NA>    NA  NA            NA 0.01 REFUSED: residuals(type = "osa") is not available for family 'hurdle_poisson': its density puts a p
33 hurdle_negbinomial          <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'hurdle_negbinomial': its density puts
34 zero_inflated_binomial      <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_inflated_binomial': its density 
35 zero_inflated_beta_binomial <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_inflated_beta_binomial': its den
36 hurdle_gamma                <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'hurdle_gamma': its density puts a poi
37 hurdle_lognormal            <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'hurdle_lognormal': its density puts a
38 zero_inflated_asym_laplace  <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_inflated_asym_laplace': its dens
39 cumulative                  oneStepGeneric  TRUE    1:4     -        0  NA            NA 0.04 truth: argument "extra" is missing, with no default                                                
40 sratio                      oneStepGeneric  TRUE    1:4     -        0  NA            NA 0.03 truth: argument "extra" is missing, with no default                                                
41 cratio                      oneStepGeneric  TRUE    1:4     -        0  NA            NA 0.01 truth: argument "extra" is missing, with no default                                                
42 acat                        oneStepGeneric  TRUE    1:4     -        0  NA            NA 0.04 truth: argument "extra" is missing, with no default                                                
43 hurdle_cumulative           <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'hurdle_cumulative': its density puts 
44 von_mises                   oneStepGeneric FALSE    -       -        0  NA      4.45e+00 1.08                                                                                                    
45 categorical                 <NA>              NA    <NA>    <NA>    NA  NA            NA 0.00 REFUSED: residuals() is not defined for a categorical family: the categories carry no order, so a r
46 multinomial                 <NA>              NA    <NA>    <NA>    NA  NA            NA   NA matrix/simplex response: not probed                                                                
47 dirichlet_multinomial       <NA>              NA    <NA>    <NA>    NA  NA            NA   NA matrix/simplex response: not probed                                                                
48 dirichlet                   <NA>              NA    <NA>    <NA>    NA  NA            NA   NA matrix/simplex response: not probed                                                                
49 dirichlet2                  <NA>              NA    <NA>    <NA>    NA  NA            NA   NA matrix/simplex response: not probed                                                                
50 logistic_normal             <NA>              NA    <NA>    <NA>    NA  NA            NA   NA matrix/simplex response: not probed                                                                
END osa-audit
```
<!-- /gen -->

After (generated, `dev/osa-out/audit-final.log`):

<!-- gen:audit-final -->
```
BEGIN osa-audit after 
   family                      method         discrete support range nan outside maxdiff  secs note                                                                                               
1  gaussian                    fullGaussian   FALSE    -       -      0  NA      6.66e-10 0.04                                                                                                    
2  student                     cdf            FALSE    -       -      0  NA      6.66e-10 0.47                                                                                                    
3  skew_normal                 cdf            FALSE    -       -      0  NA      7.50e-13 0.19                                                                                                    
4  exgaussian                  cdf            FALSE    -       -      0  NA      6.66e-10 0.09                                                                                                    
5  asym_laplace                cdf            FALSE    -       -      0  NA      2.53e-06 0.08                                                                                                    
6  huber                       oneStepGeneric FALSE    -       -     NA  NA            NA 0.02 REFUSED: residuals(type = "osa", osa_method = "oneStepGeneric") failed in oneStepPredict(): Not com
7  gen_extreme_value           cdf            FALSE    -       -      0  NA      2.99e-09 0.06                                                                                                    
8  Gamma                       cdf            FALSE    -       -      0  NA      2.05e-12 0.20                                                                                                    
9  lognormal                   cdf            FALSE    -       -      0  NA      8.01e-13 0.06                                                                                                    
10 inverse.gaussian            cdf            FALSE    -       -      0  NA      7.66e-13 0.08                                                                                                    
11 exponential                 cdf            FALSE    -       -      0  NA      8.90e-14 0.06                                                                                                    
12 weibull                     cdf            FALSE    -       -      0  NA      2.22e-14 0.07                                                                                                    
13 shifted_lognormal           cdf            FALSE    -       -      0  NA      8.02e-13 0.07                                                                                                    
14 frechet                     cdf            FALSE    -       -      0  NA      8.19e-13 0.06                                                                                                    
15 tweedie                     <NA>              NA    <NA>    <NA>  NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'tweedie': it puts a point mass on an 
16 beta                        cdf            FALSE    -       -      0  NA      6.11e-12 0.42                                                                                                    
17 xbeta                       <NA>              NA    <NA>    <NA>  NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'xbeta': it puts a point mass on an ex
18 zero_inflated_beta          <NA>              NA    <NA>    <NA>  NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_inflated_beta': it puts a point 
19 zero_one_inflated_beta      <NA>              NA    <NA>    <NA>  NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'zero_one_inflated_beta': it puts a po
20 poisson                     cdf             TRUE    -       -      0   0      2.65e-14 0.07                                                                                                    
21 negbinomial                 cdf             TRUE    -       -      0   0      3.33e-14 0.36                                                                                                    
22 negbinomial2                cdf             TRUE    -       -      0   0      3.36e-14 0.33                                                                                                    
23 nbinom1                     cdf             TRUE    -       -      0   0      3.25e-14 0.32                                                                                                    
24 compois                     cdf             TRUE    -       -      0   0      2.75e-14 0.18                                                                                                    
25 geometric                   cdf             TRUE    -       -      0   0      2.73e-14 0.27                                                                                                    
26 discrete_weibull            cdf             TRUE    -       -      0   0      2.02e-14 0.06                                                                                                    
27 binomial                    cdf             TRUE    -       -      0   0      1.87e-14 0.24                                                                                                    
28 beta_binomial               cdf             TRUE    -       -      0   0      1.18e-12 0.11                                                                                                    
29 bernoulli                   cdf             TRUE    -       -      0   0      6.02e-15 0.05                                                                                                    
30 zero_inflated_poisson       cdf             TRUE    -       -      0   0      2.15e-14 0.11                                                                                                    
31 zero_inflated_negbinomial   cdf             TRUE    -       -      0   0      2.36e-14 0.37                                                                                                    
32 hurdle_poisson              cdf             TRUE    -       -      0   0      2.01e-14 0.09                                                                                                    
33 hurdle_negbinomial          cdf             TRUE    -       -      0   0      7.22e-14 0.39                                                                                                    
34 zero_inflated_binomial      cdf             TRUE    -       -      0   0      4.20e-14 0.25                                                                                                    
35 zero_inflated_beta_binomial cdf             TRUE    -       -      0   0      6.63e-13 0.11                                                                                                    
36 hurdle_gamma                <NA>              NA    <NA>    <NA>  NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'hurdle_gamma': it puts a point mass o
37 hurdle_lognormal            <NA>              NA    <NA>    <NA>  NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'hurdle_lognormal': it puts a point ma
38 zero_inflated_asym_laplace  <NA>              NA    <NA>    <NA>  NA  NA            NA 0.01 REFUSED: residuals(type = "osa") is not available for family 'zero_inflated_asym_laplace': it puts 
39 cumulative                  oneStepGeneric  TRUE    1:4     -      0   0      1.67e-14 0.09                                                                                                    
40 sratio                      oneStepGeneric  TRUE    1:4     -      0   0      1.74e-14 0.08                                                                                                    
41 cratio                      oneStepGeneric  TRUE    1:4     -      0   0      1.74e-14 0.05                                                                                                    
42 acat                        oneStepGeneric  TRUE    1:4     -      0   0      2.24e-14 0.06                                                                                                    
43 hurdle_cumulative           <NA>              NA    <NA>    <NA>  NA  NA            NA 0.00 REFUSED: residuals(type = "osa") is not available for family 'hurdle_cumulative': it puts a point m
44 von_mises                   cdf            FALSE    -       -      0  NA      4.22e-14 0.34                                                                                                    
45 categorical                 <NA>              NA    <NA>    <NA>  NA  NA            NA 0.00 REFUSED: residuals() is not defined for a categorical family: the categories carry no order, so a r
46 multinomial                 <NA>              NA    <NA>    <NA>  NA  NA            NA   NA matrix/simplex response: not probed                                                                
47 dirichlet_multinomial       <NA>              NA    <NA>    <NA>  NA  NA            NA   NA matrix/simplex response: not probed                                                                
48 dirichlet                   <NA>              NA    <NA>    <NA>  NA  NA            NA   NA matrix/simplex response: not probed                                                                
49 dirichlet2                  <NA>              NA    <NA>    <NA>  NA  NA            NA   NA matrix/simplex response: not probed                                                                
50 logistic_normal             <NA>              NA    <NA>    <NA>  NA  NA            NA   NA matrix/simplex response: not probed                                                                
END osa-audit
```
<!-- /gen -->

Read together:

- Counts: every count family now takes `"cdf"` and matches the
  randomized quantile residual to 1e-12 or better. Before, binomial was
  NaN on 2 rows, bernoulli put 57 rows outside its own step, and
  poisson, negbinomial, negbinomial2, nbinom1, compois and
  discrete_weibull put 5 to 13 rows outside it (geometric happened to
  read within 3.4e-5). beta_binomial was refused.
- The six zero-inflated and hurdle count families were refused by name;
  they now work, because the `"cdf"` route hands the density the count
  as data.
- Continuous families were wrong too, which the brief did not name:
  Gamma 0.086 off, Beta 0.69, inverse.gaussian NaN on 3 rows and 6.77
  off, frechet NaN on 3 and 7.29 off, von_mises 4.45 off (its compat
  row said "refused upstream", which was false: it ran), tweedie NaN on
  71 of 150, and student, exgaussian, asym_laplace, huber, exponential,
  weibull, lognormal and shifted_lognormal stopped in RTMB with an S4
  type error. After: every one with a CDF matches `qnorm(F(y))` to
  3e-9 or better, except asym_laplace at 2.5e-6, which is the reference
  integral across the density's kink, not the residual: against the
  closed-form CDF on the same fit it is 1.5e-14 (`dev/osa-alap.R`,
  generated):

<!-- gen:alap -->
```
BEGIN osa-alap
max|r - qnorm(F(y))| closed form: 1.5e-14
END osa-alap
```
<!-- /gen -->
 tweedie is refused by name (a mass at zero), and huber,
  which has no CDF, stops with a refusal that says why.
- Unchanged: gaussian (`fullGaussian`, 6.7e-10), the four ordinal
  families (`oneStepGeneric` over `1:K`, 2.2e-14 or better), and the
  refusals for the continuous point-mass families, hurdle_cumulative,
  categorical, and the matrix and simplex responses.

### 1.3 Defect 3: bernoulli near -3.7

`dev/osa-latent.R` reruns `dev/review-latres-osa2.R`'s construction
(seeds 1 to 3, 15 groups of 10, latent AR(1) at 0.5 with sd 0.6). It is
the same route as 1.1 on a two-point support. Generated:

Before (`dev/osa-out/latent-before.log`):

<!-- gen:latent-before -->
```
BEGIN osa-latent before 
seed 1 bern_latent   mean -3.707 sd 0.335 min  -5.08 max  -3.04 nan 0 ks.p 2.11e-130  9.24s
seed 1 bern_re       mean -4.019 sd 0.224 min  -5.09 max  -3.66 nan 0 ks.p 1.11e-130  10.53s
seed 1 bern_ar1re    mean -3.707 sd 0.335 min  -5.08 max  -3.04 nan 0 ks.p 2.11e-130  10.37s
seed 1 bern_plain    mean -3.072 sd 1.621 min  -5.06 max   0.10 nan 0 ks.p 2.17e-63  3.14s
seed 1 binom_latent  mean -0.086 sd 0.991 min  -2.99 max   2.60 nan 0 ks.p 0.588  9.69s
seed 1 binom_re      mean -0.106 sd 1.064 min  -3.12 max   2.60 nan 0 ks.p 0.263  8.87s
seed 1 pois_latent   mean -0.028 sd 1.018 min  -2.69 max   2.69 nan 0 ks.p 0.917  10.40s
seed 1 nb_latent     mean  0.007 sd 0.978 min  -2.21 max   2.86 nan 0 ks.p 0.694  24.61s
Boundary (singular) fit: sderr = 0.0011, ar[1] = 1 in ar(t, g, cov = TRUE) are at the limits of their ranges. The data read as cosy() or (1 | group): the marginal sd of the latent, 0.272, is the number they identify, and neither of the two is known. Est.Error and the interval are NA in summary() and confint_varcorr(); diagnose() lists the limit
seed 1 ord_latent    mean -0.026 sd 0.961 min  -2.32 max   2.17 nan 0 ks.p 0.815  1.40s
seed 1 ord_re        mean -0.026 sd 0.961 min  -2.32 max   2.17 nan 0 ks.p 0.815  0.95s
seed 2 bern_latent   mean -3.786 sd 0.265 min  -4.82 max  -3.30 nan 0 ks.p 1.37e-130  11.04s
seed 2 bern_re       mean -4.143 sd 0.201 min  -5.02 max  -3.87 nan 0 ks.p 1.06e-130  12.78s
seed 2 bern_ar1re    mean -3.786 sd 0.265 min  -4.82 max  -3.30 nan 0 ks.p 1.37e-130  12.08s
seed 2 bern_plain    mean -4.323 sd 0.199 min  -5.08 max  -4.00 nan 0 ks.p 1.05e-130  2.43s
seed 2 binom_latent  mean -0.058 sd 0.923 min  -2.87 max   1.94 nan 0 ks.p 0.896  13.43s
seed 2 binom_re      mean  0.020 sd 1.127 min  -3.07 max   2.26 nan 0 ks.p 0.523  10.57s
seed 2 pois_latent   mean -0.002 sd 0.973 min  -2.48 max   2.68 nan 0 ks.p 0.611  13.13s
seed 2 nb_latent     mean  0.031 sd 0.966 min  -2.30 max   3.09 nan 0 ks.p 0.922  27.12s
Boundary (singular) fit: sderr = 0.00096, ar[1] = 1 in ar(t, g, cov = TRUE) are at the limits of their ranges. The data read as cosy() or (1 | group): the marginal sd of the latent, 0.652, is the number they identify, and neither of the two is known. Est.Error and the interval are NA in summary() and confint_varcorr(); diagnose() lists the limit
seed 2 ord_latent    mean -0.016 sd 0.940 min  -2.20 max   1.95 nan 0 ks.p 0.569  1.19s
seed 2 ord_re        mean -0.016 sd 0.940 min  -2.20 max   1.95 nan 0 ks.p 0.569  1.25s
Boundary (singular) fit: sderr = 0.002, ar[1] = 1 in ar(t, g, cov = TRUE) are at the limits of their ranges. The data read as cosy() or (1 | group): the marginal sd of the latent, 0.74, is the number they identify, and neither of the two is known. Est.Error and the interval are NA in summary() and confint_varcorr(); diagnose() lists the limit
seed 3 bern_latent   mean -3.757 sd 0.284 min  -4.67 max  -3.21 nan 0 ks.p 1.54e-130  12.89s
seed 3 bern_re       mean -3.757 sd 0.284 min  -4.67 max  -3.21 nan 0 ks.p 1.54e-130  11.94s
Boundary (singular) fit: cor_g__ft1__ft2 is at the limit of its range. The likelihood does not change as it moves further that way, so it has no standard error (NaN in summary(), vcov() and confint()); the other standard errors are kept. lme4 reports such a fit as singular; diagnose() lists the components at the boundary
seed 3 bern_ar1re    mean -3.757 sd 0.284 min  -4.67 max  -3.21 nan 0 ks.p 1.54e-130  11.93s
seed 3 bern_plain    mean -3.349 sd 1.457 min  -5.05 max   0.76 nan 0 ks.p 1.89e-72  2.86s
seed 3 binom_latent  mean -0.061 sd 0.998 min  -2.60 max   1.97 nan 0 ks.p 0.501  11.63s
seed 3 binom_re      mean -0.059 sd 1.060 min  -2.73 max   2.28 nan 0 ks.p 0.473  11.57s
seed 3 pois_latent   mean -0.029 sd 1.024 min  -2.87 max   3.00 nan 0 ks.p 0.937  12.64s
seed 3 nb_latent     mean -0.007 sd 1.010 min  -2.86 max   2.47 nan 0 ks.p 0.941  30.61s
Boundary (singular) fit: sderr = 0.0013, ar[1] = 1 in ar(t, g, cov = TRUE) are at the limits of their ranges. The data read as cosy() or (1 | group): the marginal sd of the latent, 0.309, is the number they identify, and neither of the two is known. Est.Error and the interval are NA in summary() and confint_varcorr(); diagnose() lists the limit
seed 3 ord_latent    mean -0.033 sd 1.058 min  -2.91 max   3.50 nan 0 ks.p 0.972  1.31s
seed 3 ord_re        mean -0.033 sd 1.058 min  -2.91 max   3.50 nan 0 ks.p 0.972  1.18s
END osa-latent
```
<!-- /gen -->

After (`dev/osa-out/latent-after.log`):

<!-- gen:latent-after -->
```
BEGIN osa-latent after 
seed 1 bern_latent   mean -0.027 sd 0.943 min  -2.63 max   2.02 nan 0 ks.p 0.997  2.60s
seed 1 bern_re       mean -0.025 sd 0.945 min  -2.65 max   2.05 nan 0 ks.p 0.994  2.42s
seed 1 bern_ar1re    mean -0.027 sd 0.943 min  -2.63 max   2.02 nan 0 ks.p 0.997  2.43s
seed 1 bern_plain    mean -0.023 sd 0.945 min  -2.66 max   2.06 nan 0 ks.p 0.994  0.06s
seed 1 binom_latent  mean -0.060 sd 1.014 min  -2.99 max   2.95 nan 0 ks.p 0.724  5.96s
seed 1 binom_re      mean -0.089 sd 1.081 min  -3.12 max   2.90 nan 0 ks.p 0.264  5.25s
seed 1 pois_latent   mean -0.004 sd 1.028 min  -2.68 max   2.71 nan 0 ks.p 0.837  15.53s
seed 1 nb_latent     mean  0.010 sd 0.977 min  -2.21 max   2.87 nan 0 ks.p 0.729  13.12s
Boundary (singular) fit: sderr = 0.0011, ar[1] = 1 in ar(t, g, cov = TRUE) are at the limits of their ranges. The data read as cosy() or (1 | group): the marginal sd of the latent, 0.272, is the number they identify, and neither of the two is known. Est.Error and the interval are NA in summary() and confint_varcorr(); diagnose() lists the limit
seed 1 ord_latent    mean -0.026 sd 0.961 min  -2.32 max   2.17 nan 0 ks.p 0.815  1.62s
seed 1 ord_re        mean -0.026 sd 0.961 min  -2.32 max   2.17 nan 0 ks.p 0.815  1.41s
seed 2 bern_latent   mean -0.012 sd 1.018 min  -2.47 max   3.38 nan 0 ks.p 0.981  2.47s
seed 2 bern_re       mean -0.008 sd 1.016 min  -2.48 max   3.38 nan 0 ks.p 0.978  2.21s
seed 2 bern_ar1re    mean -0.012 sd 1.018 min  -2.47 max   3.38 nan 0 ks.p 0.981  2.30s
seed 2 bern_plain    mean -0.009 sd 1.016 min  -2.49 max   3.38 nan 0 ks.p 0.994  0.06s
seed 2 binom_latent  mean -0.001 sd 0.971 min  -2.87 max   2.26 nan 0 ks.p 0.999  6.17s
seed 2 binom_re      mean  0.039 sd 1.150 min  -3.07 max   2.31 nan 0 ks.p 0.316  5.64s
seed 2 pois_latent   mean  0.019 sd 0.986 min  -2.48 max   2.68 nan 0 ks.p 0.573  16.28s
seed 2 nb_latent     mean  0.032 sd 0.966 min  -2.30 max   3.09 nan 0 ks.p 0.924  14.24s
Boundary (singular) fit: sderr = 0.00096, ar[1] = 1 in ar(t, g, cov = TRUE) are at the limits of their ranges. The data read as cosy() or (1 | group): the marginal sd of the latent, 0.652, is the number they identify, and neither of the two is known. Est.Error and the interval are NA in summary() and confint_varcorr(); diagnose() lists the limit
seed 2 ord_latent    mean -0.016 sd 0.940 min  -2.20 max   1.95 nan 0 ks.p 0.569  1.21s
seed 2 ord_re        mean -0.016 sd 0.940 min  -2.20 max   1.95 nan 0 ks.p 0.569  1.53s
Boundary (singular) fit: sderr = 0.002, ar[1] = 1 in ar(t, g, cov = TRUE) are at the limits of their ranges. The data read as cosy() or (1 | group): the marginal sd of the latent, 0.74, is the number they identify, and neither of the two is known. Est.Error and the interval are NA in summary() and confint_varcorr(); diagnose() lists the limit
seed 3 bern_latent   mean -0.013 sd 1.061 min  -2.52 max   3.45 nan 0 ks.p 0.924  2.50s
seed 3 bern_re       mean -0.013 sd 1.061 min  -2.52 max   3.45 nan 0 ks.p 0.924  2.82s
Boundary (singular) fit: cor_g__ft1__ft2 is at the limit of its range. The likelihood does not change as it moves further that way, so it has no standard error (NaN in summary(), vcov() and confint()); the other standard errors are kept. lme4 reports such a fit as singular; diagnose() lists the components at the boundary
seed 3 bern_ar1re    mean -0.013 sd 1.061 min  -2.52 max   3.45 nan 0 ks.p 0.924  2.87s
seed 3 bern_plain    mean -0.008 sd 1.063 min  -2.52 max   3.36 nan 0 ks.p 0.76  0.06s
seed 3 binom_latent  mean -0.013 sd 1.040 min  -2.60 max   2.21 nan 0 ks.p 0.636  6.81s
seed 3 binom_re      mean -0.021 sd 1.096 min  -2.72 max   2.30 nan 0 ks.p 0.473  6.24s
seed 3 pois_latent   mean  0.014 sd 1.040 min  -2.86 max   3.01 nan 0 ks.p 0.648  17.27s
seed 3 nb_latent     mean  0.004 sd 1.014 min  -2.86 max   2.47 nan 0 ks.p 0.958  19.19s
Boundary (singular) fit: sderr = 0.0013, ar[1] = 1 in ar(t, g, cov = TRUE) are at the limits of their ranges. The data read as cosy() or (1 | group): the marginal sd of the latent, 0.309, is the number they identify, and neither of the two is known. Est.Error and the interval are NA in summary() and confint_varcorr(); diagnose() lists the limit
seed 3 ord_latent    mean -0.033 sd 1.058 min  -2.91 max   3.50 nan 0 ks.p 0.972  0.97s
seed 3 ord_re        mean -0.033 sd 1.058 min  -2.91 max   3.50 nan 0 ks.p 0.972  1.79s
END osa-latent
```
<!-- /gen -->

Before reproduces the review's figure (seed 1, `bern_latent`: mean -3.707,
sd 0.335, against the review's -3.71 and 0.34). The review also asked why
the ordinal residuals read the same three decimals under a latent term
and under `(1 | g)`: in both arms each `ord_latent` fit stops at the
boundary `ar[1] = 1` with `sderr` near 0.001, and frmtmb's own message
says that model reads as `(1 | group)`, so the two fits are the same model
(the `Boundary (singular) fit` lines in both logs). Not a defect.

## 2. What changed

- `R/predict.R`, `residual_values()`: a core family with a CDF
  (`cdf_log`, or `lcdf`), uncensored, not ordinal, and not an
  untruncated gaussian, defaults to `osa_method = "cdf"`. A truncated
  gaussian takes it too (section 7). An extension's family keeps its own
  route: frmtmb.eam's Wiener `lcdf` is a first-passage CDF per boundary,
  not a CDF of the response alone.
- New in `R/predict.R`: `osa_has_cdf()`, `osa_is_cdf()`,
  `osa_log_pair()`, `osa_cdf_tails()`, `osa_cdf_loglik()`. Under
  `method = "cdf"` TMB keeps the observation as data and hands the
  objective an `osa` object with a three-column indicator (density,
  lower tail, upper tail). `R/objective.R` gets a six-line hook before
  `row_lpdf()`: on that object it adds `keep1 * ll + keep2 * log F(y) +
  keep3 * log S(y)`, from `fam_log_pair()` (or `lcdf`/`lccdf`), inside
  the `trunc()` window by `log_window()`, each tail floored at
  `log(1e-300)` as RTMB's `dGenericOSA()` floors them. The floor is
  needed: the indicator is a parameter of TMB's tape, and `0 * -Inf`
  (the survivor at a binomial's last count) is NaN there.
- Case weights are left out of the `"cdf"` terms. TMB reads a row's own
  mass from the same density term that conditions the later rows, so a
  weight `w` made that mass `p^w` (section 3.3).
- `oneStepGeneric` on a count now needs a `discreteSupport`: the ordinal
  `1:K` and a finite `trunc()` window still set one; otherwise it is
  refused with the remedy. Every method is told `discrete`, so TMB
  refuses `fullGaussian` and the other gaussian methods on a count;
  before, they returned a gaussian approximation of the count with no
  word.
- A non-finite residual warns (class `frmtmb_osa_nonfinite`), and an
  error inside `oneStepPredict()` is re-raised naming the method and,
  for a family with no CDF, why.
- `tweedie` joins `osa_point_mass_families` (refused by name).
- `R/families.R`: `bernoulli()` and `nbinom1()` declare `cdf_log`
  (`binomial_cdf_log(q, 1, .)` and `nbinom_cdf_log()` at shape
  `mu / phi`). Neither family accepts `cens()` or `trunc()`, so nothing
  else reads it. `beta_binomial_lpmf()` no longer refuses the `osa`
  object: under `oneStepGeneric` it unwraps it and forms the shifts with
  `lgamma_shift_diff()`, which takes a count on the tape.
- `R/compat.R`: the rows in section 5.

## 3. Evidence after the change

### 3.1 Random effects, against the exact one-step predictive

`dev/osa-re.R`, 10 refits per family. Base data seed 11 (20 groups of
10, random-intercept sd 0.7), then 10 responses simulated from that fit
with the group effects redrawn (`simulate(re_formula = NA, seed = 21)`),
each refitted. "exact" is the one-step predictive by quadrature over
the random intercept (801 points over +-8 sd) at the refit's estimates;
`maxdiff` and `meddiff` compare the OSA residual with it, using TMB's
uniforms. `ks.p` pools the 2000 residuals of each family. Generated:

Before (`dev/osa-out/re-before.log`):

<!-- gen:re-before -->
```
BEGIN osa-re before reps 10 
poisson        rows 2000 nan 0 err 0 | osa mean -0.010 sd 0.978 ks.p 0.661 | exact ks.p 0.829 | maxdiff 0.305 meddiff 0.000329 | 53.2s
bernoulli      rows 2000 nan 0 err 0 | osa mean -3.633 sd 0.366 ks.p 0 | exact ks.p 0.0822 | maxdiff 6.84 meddiff 3.59 | 99.3s
binomial       rows 2000 nan 0 err 0 | osa mean -0.004 sd 0.977 ks.p 0.587 | exact ks.p 0.468 | maxdiff 0.966 meddiff 0.000141 | 128.3s
negbinomial    rows 2000 nan 0 err 0 | osa mean -0.021 sd 0.976 ks.p 0.346 | exact ks.p 0.616 | maxdiff 0.292 meddiff 0.00021 | 204.1s
   beta_binomial ERR beta_binomial() does not support OSA residuals. Use another `type` of residuals() 
beta_binomial  refused on 10 of 10 refits
END osa-re
```
<!-- /gen -->

After (`dev/osa-out/re-after.log`):

<!-- gen:re-after -->
```
BEGIN osa-re after reps 10 
poisson        rows 2000 nan 0 err 0 | osa mean  0.002 sd 0.987 ks.p 0.802 | exact ks.p 0.829 | maxdiff 0.0575 meddiff 0.000365 | 185.7s
bernoulli      rows 2000 nan 0 err 0 | osa mean -0.026 sd 0.963 ks.p 0.0828 | exact ks.p 0.0822 | maxdiff 0.0197 meddiff 0.000412 | 35.5s
binomial       rows 2000 nan 0 err 0 | osa mean  0.003 sd 0.985 ks.p 0.465 | exact ks.p 0.468 | maxdiff 0.0207 meddiff 0.000197 | 107.7s
negbinomial    rows 2000 nan 0 err 0 | osa mean -0.007 sd 0.985 ks.p 0.62 | exact ks.p 0.616 | maxdiff 0.0343 meddiff 0.000231 | 195.8s
beta_binomial  rows 2000 nan 0 err 0 | osa mean -0.008 sd 1.000 ks.p 0.807 | exact ks.p 0.811 | maxdiff 0.0185 meddiff 0.000263 | 106.0s
END osa-re
```
<!-- /gen -->

The remaining `maxdiff` after (0.02 to 0.06) is the Laplace
approximation of the one-step integral, which both routes make; the
median is 2e-4 to 4e-4. bernoulli before was 3.59 off at the median.
The seconds in these two blocks were taken with both arms running at
once on a machine shared with eleven other lanes; section 3.4 is the
controlled measure.

### 3.2 trunc(), weights, the explicit methods, RTMBdist-5

`dev/osa-edge.R`, data seed 9, n = 150. Generated:

Before (`dev/osa-out/edge-before.log`):

<!-- gen:edge-before -->
```
BEGIN osa-edge before 
poisson trunc(lb = 1)              nan 0 outside 2 maxdiff 0.254
binomial trunc(ub = 7)             ERROR 'from' must be a finite number
poisson trunc(lb = row-varying)    ERROR residuals(type = "osa") needs trunc() bounds that are the same for every observation; got row-varying bounds
poisson weights(1, 3)              nan 0 outside 42 maxdiff 1.85
gaussian weights(1, 3)             max|r - z| 1.97  max|r - z sqrt(w)| 4.44e-16
gaussian weights, cdf              max|r - z| 2.65
binomial oneStepGeneric 0:10       nan 0 outside 0 maxdiff 3.47e-14
binomial oneStepGeneric no support nan 3
binomial fullGaussian              nan 0
binomial oneStepGaussianOffMode    nan 0
beta_binomial default              ERROR beta_binomial() does not support OSA residuals. Use another `type` of residuals()
beta_binomial oneStepGeneric 0:10  ERROR beta_binomial() does not support OSA residuals. Use another `type` of residuals()
beta_binomial oneStepGaussianOffMode ERROR beta_binomial() does not support OSA residuals. Use another `type` of residuals()
huber default                      ERROR Not compatible with requested type: [type=S4; target=double].
inverse.gaussian oneStepGeneric    nan 1
inverse.gaussian default           nan 1
RTMBdist-5 dbetabinom tape  log phi 10: -5.689958   binomial -5.690003
RTMBdist-5 dbetabinom tape  log phi 39: -1536.000000   binomial -5.690003
RTMBdist-5 dbetabinom tape  log phi 45: 0.000000   binomial -5.690003
END osa-edge
```
<!-- /gen -->

After (`dev/osa-out/edge-final.log`):

<!-- gen:edge-final -->
```
BEGIN osa-edge after 
poisson trunc(lb = 1)              nan 0 outside 0 maxdiff 2.23e-14
binomial trunc(ub = 7)             nan 0 outside 0 maxdiff 2.09e-14
poisson trunc(lb = row-varying)    nan 0 outside 0 maxdiff 2.99e-14
poisson weights(1, 3)              nan 0 outside 0 maxdiff 3e-14
gaussian weights(1, 3)             max|r - z| 1.97  max|r - z sqrt(w)| 4.44e-16
gaussian weights, cdf              max|r - z| 1.35e-14
binomial oneStepGeneric 0:10       nan 0 outside 0 maxdiff 3.47e-14
binomial oneStepGeneric no support ERROR residuals(type = "osa", osa_method = "oneStepGeneric") on a discrete response sums the mass over `discreteSupport`, and none is known for this fit. Pa
binomial fullGaussian              ERROR residuals(type = "osa", osa_method = "fullGaussian") failed in oneStepPredict(): method='fullGaussian' is not for discrete observations
binomial oneStepGaussianOffMode    ERROR residuals(type = "osa", osa_method = "oneStepGaussianOffMode") failed in oneStepPredict(): method='oneStepGaussianOffMode' is not for discrete observa
beta_binomial default              nan 0 outside 0 maxdiff 4.26e-14
beta_binomial oneStepGeneric 0:10  nan 0 outside 0 maxdiff 4.26e-14
beta_binomial oneStepGaussianOffMode ERROR residuals(type = "osa", osa_method = "oneStepGaussianOffMode") failed in oneStepPredict(): method='oneStepGaussianOffMode' is not for discrete observa
huber default                      ERROR residuals(type = "osa", osa_method = "oneStepGeneric") failed in oneStepPredict(): Not compatible with requested type: [type=S4; target=double]. The '
inverse.gaussian oneStepGeneric    nan 1 WARN: 1 of 150 one-step-ahead residuals are not finite: oneStepPredict(method = "oneStepGeneric") could not form the one-step 
inverse.gaussian default           nan 0
RTMBdist-5 dbetabinom tape  log phi 10: -5.689958   binomial -5.690003
RTMBdist-5 dbetabinom tape  log phi 39: -1536.000000   binomial -5.690003
RTMBdist-5 dbetabinom tape  log phi 45: 0.000000   binomial -5.690003
END osa-edge
```
<!-- /gen -->

### 3.3 Case weights

Under `"cdf"` with the weight on the density term (an intermediate
version of this lane), the poisson rows at weight 3 stayed inside their
step but sat up to 1.45 from the randomized quantile residual, because
TMB's mass for the row was `p^3`. With the weights left out it is 3e-14
(block above). Before, the generic route put 42 of 150 rows outside
their step at the same weights. `gaussian()`'s `fullGaussian` reads a
weight as a precision, `(y - mu) sqrt(w) / sigma`, to 4.4e-16 in both
arms, so the old compat note ("the residuals ignore the case weights")
was false there; under `osa_method = "cdf"` a gaussian residual ignores
them (1.35e-14).

### 3.4 Cost

`dev/osa-time.R`, one process of the lane build: the old route (the
0.79.1 call, made on the fit's own objective) against `"cdf"`,
interleaved, minimum of 5 rounds, with a control that runs `"cdf"`
twice. Seed 11, 20 groups of 10. Generated:

<!-- gen:time -->
```
BEGIN osa-time
poisson       old   1.42s  cdf   0.03s  control   0.03s  cdf/old  0.02  control/cdf  1.00
poisson_re    old   8.41s  cdf  22.06s  control  19.86s  cdf/old  2.62  control/cdf  0.90
bernoulli     old   2.88s  cdf   0.06s  control   0.06s  cdf/old  0.02  control/cdf  1.00
Boundary (singular) fit: sd_g__Intercept is at zero. The likelihood does not change as it moves further that way, so it has no standard error (NaN in summary(), vcov() and confint()); the other standard errors are kept. lme4 reports such a fit as singular; diagnose() lists the components at the boundary
bernoulli_re  old  13.14s  cdf   2.49s  control   2.58s  cdf/old  0.19  control/cdf  1.04
END osa-time
```
<!-- /gen -->

The poisson `(1 | g)` fit is 2.6 times slower: the cdf tape carries
every row's CDF pair (`poisson_cdf_log()` sums up to 70 terms a
row), and each of the 3 evaluations a row needs runs the inner
Laplace optimization on it. The bernoulli `(1 | g)` design here
fitted at the boundary (sd 0, the message in the block), so its
0.19 is close to the no-group case; the bernoulli refits of 3.1
(35.5 s against 99.3 s, both arms at once) are the other
reading.

## 4. Tests

### 4.1 The pins, seen failing on rellib-r18

`tests/testthat/test-osa-cdf.R` (new). On rellib-r18 (`OSA_ARM=before
bash dev/osa-par.sh`), behavioral failures, not missing symbols:

<!-- gen:t0-before -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18; frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18
DETAIL FAIL [count residuals are randomized quantile residuals of the fit] Expected `sum(!ok)` to be identical to 0L. | Differences: |   `actual`: 3 | `expected`: 0 |  | poisson
DETAIL FAIL [count residuals are randomized quantile residuals of the fit] Expected `sum(!ok)` to be identical to 0L. | Differences: |   `actual`: 3 | `expected`: 0 |  | binomial
DETAIL FAIL [count residuals are randomized quantile residuals of the fit] Expected `sum(!ok)` to be identical to 0L. | Differences: |   `actual`: 53 | `expected`:  0 |  | bernoulli
DETAIL FAIL [count residuals are randomized quantile residuals of the fit] Expected `sum(!ok)` to be identical to 0L. | Differences: |   `actual`: 8 | `expected`: 0 |  | negbinomial
DETAIL ERROR [count residuals are randomized quantile residuals of the fit] <frmtmb_error/error/condition> | Error: beta_binomial() does not support OSA residuals. Use another `type` of residuals()
DETAIL FAIL [continuous residuals are qnorm of the fitted CDF] Expected `max(abs(r - ref))/max(abs(ref))` < 1e-08. | Actual comparison: 0.04297798 >= 0.00000001 | Difference: 0.04297797 >= 0
DETAIL FAIL [continuous residuals are qnorm of the fitted CDF] Expected `max(abs(r - ref))/max(abs(ref))` < 1e-08. | Actual comparison: 0.10462070 >= 0.00000001 | Difference: 0.10462069 >= 0
DETAIL ERROR [continuous residuals are qnorm of the fitted CDF] Error in `h(simpleError(msg, call))`: error in evaluating the argument 'x' in selecting a method for function 'dnorm': non-numeric argument to mathematical function
DETAIL ERROR [a zero-inflated or hurdle count reads its point mass as data] <frmtmb_error/error/condition> | Error: residuals(type = "osa") is not available for family 'zero_inflated_poisson': its density puts a point mass on an e
DETAIL ERROR [a truncated count is read inside its window, row by row] <frmtmb_error/error/condition> | Error: residuals(type = "osa") needs trunc() bounds that are the same for every observation; got row-varying bounds
DETAIL ERROR [a truncated gaussian with a group-level term has every row] Error in `pred$nlcdf.lower`: $ operator is invalid for atomic vectors
DETAIL FAIL [a case weight does not reshape a row's own step] Expected `max(abs(u1 - u2))/max(u1)` < 1e-06. | Actual comparison: 1.068804 >= 0.000001 | Difference: 1.068803 >= 0
DETAIL FAIL [bernoulli residuals are calibrated under group-level terms] Expected `abs(mean(r))/(stats::sd(r)/sqrt(length(r)))` < 4. | Actual comparison: 28.7 >= 4.0 | Difference: 24.7 >= 0
DETAIL FAIL [bernoulli residuals are calibrated under group-level terms] Expected `stats::sd(r)` < 1.25. | Actual comparison: 1.58 >= 1.25 | Difference: 0.33 >= 0
DETAIL FAIL [bernoulli residuals are calibrated under group-level terms] Expected `abs(mean(r))/(stats::sd(r)/sqrt(length(r)))` < 4. | Actual comparison: 168.8 >= 4.0 | Difference: 164.8 >= 0
DETAIL FAIL [bernoulli residuals are calibrated under group-level terms] Expected `stats::sd(r)` > 0.75. | Actual comparison: 0.31 <= 0.75 | Difference: -0.44 <= 0
DETAIL FAIL [bernoulli residuals are calibrated under group-level terms] Expected `abs(mean(r))/(stats::sd(r)/sqrt(length(r)))` < 4. | Actual comparison: 168.8 >= 4.0 | Difference: 164.8 >= 0
DETAIL FAIL [bernoulli residuals are calibrated under group-level terms] Expected `stats::sd(r)` > 0.75. | Actual comparison: 0.31 <= 0.75 | Difference: -0.44 <= 0
DETAIL FAIL [the routes that cannot read a count refuse it by name] Expected `osa_r(fit, osa_method = "oneStepGeneric")` to throw a error.
DETAIL FAIL [the routes that cannot read a count refuse it by name] Expected `osa_r(fit, osa_method = "fullGaussian")` to throw a error.
DETAIL FAIL [the routes that cannot read a count refuse it by name] Expected `osa_r(ft)` to throw a error.
DETAIL FAIL [a non-finite residual is reported, not passed on silently] Expected `any(grepl(r, seen, fixed = TRUE))` to be TRUE. | Differences: | `actual`:   FALSE | `expected`: TRUE  |  | expected a warning containing '1 of 150 o
RESULT test-osa-cdf.R pass=10 fail=17 err=5 skip=0 warn=0
```
<!-- /gen -->

On the lane build:

<!-- gen:t0-after -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-osa-cdf.R pass=41 fail=0 err=0 skip=0 warn=0
```
<!-- /gen -->

The two existing tests whose expectations changed (beta_binomial and the
hurdle counts now work), on rellib-r18:

<!-- gen:t3-before -->
```
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18; frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18
DETAIL ERROR [simulate() and the OSA residuals keep their routes] <frmtmb_error/error/condition> | Error: beta_binomial() does not support OSA residuals. Use another `type` of residuals()
RESULT test-beta-binomial-density.R pass=92 fail=0 err=1 skip=0 warn=0
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18; frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18
DETAIL ERROR [one-step-ahead residuals read a count's point mass as data] <frmtmb_error/error/condition> | Error: residuals(type = "osa") is not available for family 'hurdle_poisson': its density puts a point mass on an exact resp
RESULT test-hurdle-nb-zoib.R pass=82 fail=0 err=1 skip=0 warn=0
```
<!-- /gen -->

`test-cdfs.R` and `test-cens-trunc.R` changed only the expected message
of a refusal that still happens (bernoulli and nbinom1 with `cens()`
are now refused by their list of addition terms, since they have a
CDF); `test-xbeta-zibb-hurdle-cum.R` changed one compat expectation.

Under OpenBLAS with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
`OPENBLAS_NUM_THREADS=4`, the final build). 0.3.32 was not run: no test
here sits near a rounding-decided branch (every bound is a ratio to the
row's own step or to the largest residual). Generated:

<!-- gen:openblas -->
```
OpenBLAS 0.3.26 with LAPACK (dev/ciharden-openblas.sh 0.3.26 lapack), OPENBLAS_NUM_THREADS=4;
dgemm 2000x2000: 0.2 s in this R, 10.02 s in the reference R (same machine, one run each)
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-beta-binomial-density.R pass=96 fail=0 err=0 skip=0 warn=0
lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib; frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
RESULT test-osa-cdf.R pass=41 fail=0 err=0 skip=0 warn=0
```
<!-- /gen -->

### 4.2 Suites

Every core file and every file of the eight extension suites, one file
per process (`dev/osa-par.sh`, `dev/osa-test1.R`), `NOT_CRAN=true`, the
lane's core in front of rellib-r18 (the extensions from rellib-r18,
built from the base commit). The core suite ran on the build before the
last install, which changed only the wording of one refusal; the new
test file was rerun on the final build (above). The skips are the gated
tier's; `test-brms-lik-cdfs.R` and `test-bcm-binomial.R` were also run
with `FRMTMB_BRMS_FIT_TESTS=true` (the `gated-after` rows, counted in
the core line). Generated (`dev/osa-suites-sum.sh`):

<!-- gen:suites -->
```
BEGIN suites
SKIPS frmtmb RESULT test-bcm-binomial.R pass=17 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb RESULT test-bcm-esp.R pass=5 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-bcm-gcm.R pass=9 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-bcm-model-selection.R pass=13 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-bcm-psychophysics.R pass=7 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-bcm-signal-detection.R pass=6 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-brms-agreement.R pass=172 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-brms-lik-cdfs.R pass=0 fail=0 err=0 skip=5 warn=0
SKIPS frmtmb RESULT test-brms-likelihood.R pass=5 fail=0 err=0 skip=55 warn=0
SKIPS frmtmb RESULT test-ce-aterms.R pass=19 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-curve-scale.R pass=0 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-drmtmb-agreement.R pass=0 fail=0 err=0 skip=13 warn=0
SKIPS frmtmb RESULT test-emmeans.R pass=92 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-fuzz.R pass=0 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-gp-by.R pass=71 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-latres-brms.R pass=30 fail=0 err=0 skip=5 warn=0
SKIPS frmtmb RESULT test-mi-mm-brms-l1.R pass=0 fail=0 err=0 skip=10 warn=0
SKIPS frmtmb RESULT test-offset-grid.R pass=27 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-refcat-brms.R pass=0 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb RESULT test-bcm-bart.R pass=8 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-bcm-data-analysis.R pass=20 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb RESULT test-bcm-gaussian.R pass=11 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-bcm-latent-mixtures.R pass=30 fail=0 err=0 skip=5 warn=0
SKIPS frmtmb RESULT test-bcm-mpt.R pass=15 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-bcm-retention.R pass=15 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-bcm-simple.R pass=8 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-brms-lik-gp-exact.R pass=5 fail=0 err=0 skip=4 warn=0
SKIPS frmtmb RESULT test-brms-methods.R pass=18 fail=0 err=0 skip=47 warn=0
SKIPS frmtmb RESULT test-brms-priors.R pass=0 fail=0 err=0 skip=12 warn=0
SKIPS frmtmb RESULT test-ce-predict-median.R pass=11 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-fcor.R pass=50 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-mo-factor-re.R pass=34 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-nb2dm-brms.R pass=0 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb RESULT test-nl-brms.R pass=0 fail=0 err=0 skip=10 warn=0
SKIPS frmtmb RESULT test-nl-stan-funs.R pass=122 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb RESULT test-ordinal-mixture-brms.R pass=0 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb RESULT test-rl-example.R pass=96 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb RESULT test-simplex-brms.R pass=0 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb.coupling RESULT test-scale.R pass=0 fail=0 err=0 skip=5 warn=0
SKIPS frmtmb.eam RESULT test-sample-ndt-default.R pass=7 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.eam RESULT test-scale.R pass=0 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb.latent RESULT test-scale.R pass=0 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.learn RESULT test-scale.R pass=0 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.learn RESULT test-stan-identity.R pass=0 fail=0 err=0 skip=11 warn=0
SKIPS frmtmb.ode RESULT test-scale.R pass=0 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb.sample RESULT test-default-priors-fit.R pass=0 fail=0 err=0 skip=6 warn=0
SKIPS frmtmb.sample RESULT test-loglik-autocor.R pass=110 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.sample RESULT test-loo-predict.R pass=114 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb.sample RESULT test-loo.R pass=138 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.sample RESULT test-mi-draws-brms.R pass=0 fail=0 err=0 skip=3 warn=0
SKIPS frmtmb.sample RESULT test-mo-coding.R pass=9 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.sample RESULT test-sampling-ported.R pass=207 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb.sample RESULT test-scale.R pass=0 fail=0 err=0 skip=1 warn=0
SKIPS frmtmb.spatial RESULT test-brms-likelihood-car.R pass=4 fail=0 err=0 skip=2 warn=0
SKIPS frmtmb.spatial RESULT test-sar.R pass=106 fail=0 err=0 skip=4 warn=0
frmtmb.spatial     files    8 pass    478 fail 0 err 0 skip 6 warn 0
frmtmb.eam         files   33 pass   1848 fail 0 err 0 skip 5 warn 0
frmtmb.sample      files   85 pass   3962 fail 0 err 0 skip 18 warn 0
frmtmb.coupling    files   11 pass    551 fail 0 err 0 skip 5 warn 0
frmtmb             files  311 pass  22664 fail 0 err 0 skip 239 warn 0
frmtmb.survival    files   15 pass    552 fail 0 err 0 skip 0 warn 0
frmtmb.ode         files   11 pass    549 fail 0 err 0 skip 1 warn 0
frmtmb.latent      files   12 pass    683 fail 0 err 0 skip 2 warn 0
frmtmb.learn       files   16 pass    474 fail 0 err 0 skip 13 warn 0
files with no RESULT line: 0
frmtmb loaded from (files per lib line):
   502 frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
END suites
```
<!-- /gen -->

### 4.3 One session, three orders

The core suite in one R session, forward, reverse and with every
extension preloaded (`dev/osa-session.R`, a copy of
`dev/release/run-session.R` reading rellib-r18), on the final build.
`files=296/309` in each: 13 files hold only gated brms tests and source
to no test row. Every LEAK line is tagged with the load of an
extension. The extension suites were not run in one session (section
6). Generated:

<!-- gen:sessions -->
```
BEGIN sessions
SESSION frmtmb forward files=296/309 pass=22557 fail=0 err=0 skip=239 warn=0
SESSION frmtmb preload files=296/309 pass=22800 fail=0 err=0 skip=241 warn=0
SESSION frmtmb reverse files=296/309 pass=22707 fail=0 err=0 skip=239 warn=0
forward: LEAK lines 10, untagged 0; DETAIL FAIL/ERROR/WARN 0
reverse: LEAK lines 10, untagged 0; DETAIL FAIL/ERROR/WARN 0
preload: LEAK lines 0, untagged 0; DETAIL FAIL/ERROR/WARN 0
      3 lib: frmtmb 0.79.1 C:/Users/adf44/source/r/wt-osa-lib
END sessions
```
<!-- /gen -->

### 4.4 R CMD check

`bash dev/osa-check.sh`: built and checked in `dev/osa-check/`, `--as-cran`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, `R_LIBS` = wt-osa-lib, rellib-r18,
the user library. Generated (the Status line and the NOTEs):

<!-- gen:check -->
```
BEGIN check
* checking whether package 'frmtmb' can be installed ... [90s] OK
* checking R code for possible problems ... [140s] OK
* checking examples ... [93s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
                 user system elapsed
frm_hazard_reads 3.01   0.17    5.16
* checking examples with --run-donttest ... [99s] OK
* checking tests ...
  Running 'testthat.R' [43m]
 [43m] OK
* checking re-building of vignette outputs ... [359s] OK
* checking HTML version of manual ... [36s] NOTE
Skipping checking math rendering: package 'V8' unavailable
Status: 2 NOTEs
END check
```
<!-- /gen -->

## 5. Compat rows changed (`R/compat.R`)

- `residuals_osa x kind:family` (conditional): note rewritten for the
  `"cdf"` route.
- New, `works`: `residuals_osa` x `zero_inflated_poisson`,
  `zero_inflated_negbinomial`, `zero_inflated_binomial`,
  `zero_inflated_beta_binomial`, `hurdle_poisson`, `hurdle_negbinomial`
  (the `group:point_mass` row stays `refused` for the rest).
- New, `refused`: `residuals_osa x tweedie`. (huber is not a feature of
  the compat vocabulary, so it gets no row; `test-compat.R` refused the
  one first written.)
- `residuals_osa x von_mises`: `refused` to `works`.
- `residuals_osa x trunc()` (conditional): row-varying bounds work under
  `"cdf"`.
- `residuals_osa x weights()` (conditional): note rewritten (3.3).

## 6. Not done, and why

- The continuous point-mass families (hurdle_gamma, hurdle_lognormal,
  zero_inflated_beta, zero_inflated_asym_laplace, zero_one_inflated_beta,
  xbeta) stay refused. TMB's `"cdf"` method can randomize at a point
  mass of a continuous response (`discrete = FALSE` with
  `discreteSupport = 0`), and four of them have a CDF, but that is new
  scope and no reference here measures it.
- `cens()` on a count stays refused: the window there is `[lo + 1, hi]`
  and nothing measures it (the existing refusal says so).
- First pass: with group-level terms the `"cdf"` route was slower on
  poisson (section 3.4). Fixed in the fix round by blocks (section 0.3).
- Fix round: the Laplace error of the one-step integrals (up to 0.15 on
  a group's first row) is stated, not reduced; `quadrature = TRUE`
  cannot help, since TMB's one-step machinery needs the random effects
  on the Laplace tape.
- Fix round: the core suite was not rerun in one session, three ways, on
  the fix build; it ran one file per process, with every extension
  suite, and in R CMD check's single session (section 0.6).
- Fix round: a weight between 0 and 1 gives `NA` with a warning, not a
  residual, on every family. For a continuous response the fresh-copy
  residual is available (section 0.1); the user decides.
- huber has no CDF; it is refused with the reason rather than given one.
- Extension families keep their routes; their compat rows are
  unchanged, and their OSA behavior was not audited (no extension code
  changed, and every extension suite passes on this core).
- The extension suites ran one file per process only, not in one
  session; the core suite ran in one session three ways.
- The ported brms-suite ledger and the equivalence record were not
  rerun. A grep for "osa" in `dev/brmsport-*` finds nothing, and in
  `dev/equiv-data/results.tsv` one row, p093 (`hurdle_cumulative` x
  `residuals`), whose measured call is `type = "response"`, which this
  change does not reach. brms has no OSA residual, so no pair compares
  one.

## 7. trunc() with a group-level term, and what is left

Under `oneStepGeneric`, `row_lpdf()` subtracts each row's truncation
normalizer `log Z_j(u)` whatever its data-term indicator, so a later row
whose density is switched off still reweights the random-effect integral
by `1 / Z_j(u)`; and the per-row quadrature failed on rows outright
("non-finite function value" from `integrate()`). Under `"cdf"` the
normalizer sits inside the density term, which the indicator switches
off. So a truncated `gaussian()`, which kept `oneStepGeneric` by
default, now takes `"cdf"` too. `dev/osa-trunc-re.R` (seed 31, 12 groups
of 8, the test's design), against the exact one-step predictive by
quadrature. Generated:

<!-- gen:trunc-re -->
```
BEGIN osa-trunc-re
cdf          rows 1:96 nan 0 max|r - exact| 0.00158 median 0.000222 mean 0.017 sd 1.041
0.79.1 call  rows 1:24 integrate() failed on 24 rows; oneStepPredict() stopped: $ operator is invalid for atomic vectors  300s
END osa-trunc-re
```
<!-- /gen -->

Not fixed: the same normalizer term under `cens()` with `trunc()` and a
group-level term (that pair keeps `oneStepGeneric`, the only route that
conditions on censoring events), and under an extension's family with
`trunc()` and `oneStepGeneric`. Not measured.

## 8. Bump and floor

Minor (BREAKING: the default OSA route and every residual it gives on a
non-gaussian core family change). No extension's code calls the changed
functions; no floor moves.
Every extension suite passes on this core with the extensions of the
base commit (section 4.2), so no extension needs a new frmtmb either.

`roxygen2::roxygenise()` with roxygen2 8.1.1 rewrites
`Config/roxygen2/version` in DESCRIPTION from 8.1.0; the line was set
back by hand, so DESCRIPTION carries no change from this lane.
