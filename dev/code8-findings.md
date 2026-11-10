# Lane code8: code 8 at the maximum, fixture 1, the SE verdict on a ridge

Round opened 2026-11-09, base e0cce906 (frmtmb 0.80.0). Brief:
`C:/Users/adf44/source/r/round-briefs/r081-code8.md` (Parts A and B,
and Part C added on 2026-11-09 by the user's decision).

Status: change complete; final suite runs and R CMD check at the end
of this file. Numbers below are pasted from the logs named
beside them (`dev/code8-out/`, gitignored). "before" is rellib-r19
(0.80.0); "after" is the lane library.

## Scripts

- `dev/code8-common.R`: shared setup. `CODE8_ARM=before` puts rellib-r19
  first; `after` puts `wt-code8-lib` first. `stop_info()` measures each
  stop with the package's own `conv_headroom()` and `sdr_of()`.
- `dev/code8-battery.R`: the bigx review's battery (products and powers
  of shifted columns at 1e4 to 1e6, seeds 41 to 45; raw polynomials on
  [1000, 1010], seeds 501 to 503, degrees 5, 6 and 8; `y ~ w` at 1e7 to
  1e9, gaussian and poisson, seeds 901 to 910), and an extension `ywx`
  (seeds 911 to 940 at 1e7 to 1e9, seeds 901 to 920 at 1e10). Each fit
  is compared with the same model on centered or orthogonal columns.
- `dev/code8-short.R`: code-8 stops that are short. Forced: nlminb's
  `xf.tol` made coarse, `restarts = 0`, on a gaussian LMM with a
  correlated slope, a poisson GLMM, a student location-scale model, a
  nonlinear exponential decay and a bernoulli GLMM. Natural:
  `asym_laplace()` (the check-loss kink) against quantreg's exact
  maximum.
- `dev/code8-summarize.R`: generates the counts below from the TSVs.
- `dev/code8-fix1*.R`: fixture 1 (Part B). `dev/code8-ridge*.R`: the
  ridge (Part C). `dev/code8-mm22.R`, `dev/code8-mpt.R`: two fits whose
  verdict moved during development.

## Part A: a classed message for "false convergence (8)" at the maximum

### Battery on rellib-r19 (before) and on the lane (after)

`Rscript dev/code8-summarize.R` (generated, pasted):

    battery before-base: 99 fits, 32 at code 8; gain 1.95e-09 to 9.92e-06; ...
      se all finite 32; warned (conv_code) 32; message 0; with another warning 3
    battery before-ywx: 220 fits, 77 at code 8; gain 2.42e-11 to 1.98e-06; ...
      se all finite 77; warned (conv_code) 77; message 0; with another warning 0
    battery after-base: 99 fits, 32 at code 8; ... warned (conv_code) 0; message 32
    battery after-ywx: 220 fits, 77 at code 8; ... warned (conv_code) 0; message 77
    before, all: 109 code-8 fits, gain 2.42e-11 to 9.92e-06

The 3 "with another warning" are the raw polynomials of seed 503, which
also carry `frmtmb_design_capped` (the model is not the one written);
the capped warning stays on the lane, the code-8 warning becomes the
message. Every code-8 fit is within 6.1e-6 below to 2.3e-5 above the
centered fit (`dll -6.08e-06 to 2.23e-05` over the 77 ywx fits; the
capped polynomials are 0.73 to 3.26 below the uncapped model, as the
capped warning says).

0.80.0's design basis did not change the counts against the review's
(29 of 90): 32 of 99 here, with the polynomials included.

### Short code-8 stops

    short before: 64 fits, 21 at code 8, 19 of them short by more than 1e-6
      (dll -1.64e+00 to -1.02e-06), gain on those 0.000068 to 1.785775;
      warned 19 of 19 short, message on 0 short
    short after: 64 fits, 21 at code 8, 19 of them short by more than 1e-6
      ... warned 19 of 19 short, message on 0 short;
      code-8 not short: 2 (gain 0.00e+00 to 3.34e-08, message 1)

On the forced stops the Newton gain tracks the true shortfall closely
(gain 4.93e-4 against dll -4.94e-4, 2.11e-2 against -2.16e-2, 0.983
against -0.920, ...; `dev/code8-out/short-before.log`). The one forced
stop that gets the message (student, seed 1, `xf.tol = 1e-4`) is 3.3e-8
below the default fit.

### The bound: 2e-5

A second forced battery (`dev/code8-short.R after-x 1,...,8
1e-3,3e-3,0.03`, `dev/code8-short-x-sum.R`) fills the range between
1e-6 and 1e-4:

    fits 124 code 8 78 short (> 1e-6) 74
    short: gain 0.0 24.2 ...
    short with gain in (1e-6, 1e-4]: 14

and on every one of those 14 the gain equals the shortfall to three
digits (1.665e-5 against -1.665e-5, 3.632e-5 against -3.629e-5, 3.864e-6
against -3.863e-6, ...). There is no gap between "at the maximum" and
"short": the Newton gain IS the shortfall on these models. So the bound
is a tolerance on the shortfall that a code-8 stop may carry without a
warning, and it has to clear what rounding alone reads at a true
maximum. On the 109 code-8 stops at the maximum, rounding reads up to
9.92e-6 (gaussian `w1 * w2` at 1e5, seed 42, a fit 9.4e-6 ABOVE the
centered one). 1e-5 (the review's proposal) clears that by 0.8 percent;
2e-5 by a factor of 2, and a shortfall of 2e-5 in log-likelihood is
statistically nothing. The bound is 2e-5 (`conv_at_max_tol`).

On the final build (bound 2e-5), both forced batteries together
(`dev/code8-short-bound.R`, generated):

    code-8 stops 99 ; short by more than 2e-5: 80 warned 80 message 0
    short by 1e-6 to 2e-5: 13 message 7 warned 6 ; their gain/shortfall
      ratio 0.001569 1.000205
    gain/shortfall on the far ones: 0.8562 1.6363

The 7 messages are on stops short by 3.6e-6 to 1.95e-5, and the message
states that gain. The ratio 0.0016 is asym_laplace seed 15, whose kink
hides the shortfall from the quadratic model (it keeps the warning:
its SEs are lost).

The natural short stops (asym_laplace seeds 14 to 16, 1.0e-6 to 1.6e-6
short of quantreg's exact maximum) keep the warning, because their
standard errors are lost (the kink).

### How it composes with the 0.80.0 classes

- The warning keeps class `frmtmb_conv_code` wherever the stop is not
  shown to be the maximum; `frmtmb_conv_gradient` and
  `frmtmb_conv_gradient_nonfinite` are unchanged. The message is
  `frmtmb_converged_at_max` (plus `frmtmb_message`, `message`).
- The message is decided after the gradient verdict: if anything else
  in `check_convergence()` warned, the code-8 warning is given instead,
  first in order.
- `nl_stop_warning()` reads warning classes only, so a message does not
  start the nonlinear ridge pass.
- Tests: `allow_warnings()` matches warnings by text and does not see a
  message; a test that wraps a fit that may stop at code 8 in
  `expect_no_message()` uses the new helper `allow_at_max()`
  (`tests/testthat/helper-warnings.R`), which muffles that class only.
- Users: `suppressMessages()`, or
  `withCallingHandlers(frmtmb_converged_at_max = function(m)
  invokeRestart("muffleMessage"))`. `fit$opt$convergence` keeps nlminb's
  code (1).
- `check_se = "ignore"` keeps the warning: the SEs are not asked about,
  so "every SE is finite" is not known.
- An exact `gp()` term keeps the warning: lane ncgp measured 18 of 73
  stops whose Newton clause said "at the maximum" short by up to 1 nat
  on the gp marginal's ridge; its warm-start proof (`gp_ok`) is the only
  way such a fit is quiet, as before.

## Part B: brms fixture 1

### Diagnosis (rellib-r19)

`dev/code8-fix1.R`, `dev/code8-fix1-probe.R`, `dev/code8-fix1-probe2.R`,
`dev/code8-fix1-jitter.R` (logs `fix1-*-before.log`).

- The "sd at -723" is not an sd. It is `theta_3`, the correlation
  parameter of `(1 + Trt | visit)` (theta_idx 1, 2, 3 of a 2-by-2 `us`
  block): the correlation is -0.999999 at the stop. With 4 visits this
  is lme4's singular fit. `theta_4`, the sd of `s(Age)`, is at -354.9 on
  the log scale (the smooth reduced to its linear part); its Hessian row
  is NaN.
- Profile of `theta_3` (others re-optimized): objective 113.585 at
  `theta_3 = -0.5` (cor -0.447), 111.399 at -2, 110.161 at -5, 109.627 at
  -10, 109.378 at -30, 109.3436 at -100, 109.3401041 at -723.6. Monotone
  toward cor = -1: a maximum at the boundary, not a flat direction and
  not an unbounded one along this parameter.
- Reductions that lme4 and glmmTMB can express (`count ~ Trt * Age +
  volume + offset(Age) + (1 + Trt | visit)`, with and without `Exp` as a
  number): `lmer()` gaussian is singular with cor -1 on both
  (logLik -112.6331557 and -112.6113166); frmtmb gaussian matches to
  1.4e-6 and 4.7e-6 with the correlation parameter at -3335 and -1773;
  glmmTMB `t_family()` with `dispformula = ~ Trt` stops at "false
  convergence (8)" with cor -0.999999 and pdHess FALSE, its logLik NA.
- The fixed design is rank 5 of 7 (`X` as frmtmb builds it):
  `s(Age).fx1` is linear in Age (R^2 = 1 to the printed digits), an
  exact alias the SE check names (se_design_null() was written for this
  very fixture). The other deficiency is the `mo()` placeholder column.
- From 40 jittered starts (sd 1, seed 20261008) the base build stops in
  several basins: 109.3400 (the default's, cor at -1 and the smooth's sd
  at 0) on most, but also 105.3995 to 105.4067 (lower objective, higher
  likelihood) with `nu` run to 1 + exp(438 to 583) (the gaussian limit)
  or to 1 + exp(-14) (the Cauchy limit), 106.62, 107.95 and 109.08 with
  no parameter at an edge. So the student likelihood of this 40-row
  model has several local maxima, each at a boundary of a different
  parameter or in the interior; the default stop is not the global one.
  The -1308.3 that lane optrap's review recorded (its lane build, sd 1)
  was not reached on rellib-r19 from these 40 starts (lowest 105.3995).
  brms's fixture is a Bayesian fit whose priors (lkj on the correlation,
  gamma(2, 0.1) on nu, the sd priors) hold all three edges away.
- Why the SEs were NaN: the stop is at code 8, and `conv_headroom()`
  returned NA (the curvature "unusable"), because the correlation's row
  is noise at -723. So the fit was called short, `se_explained` was set
  to "convergence", and `sdreport()`'s own inverse was kept, with
  negative variances for all 19 parameters ("nocurv"). The SE check's
  own analysis at the same stop was right: flat for `(Intercept)`,
  `Age`, `s(Age).fx1` (the alias), boundary for `theta_3`, nonfinite for
  `theta_4`, and finite SEs for everything else, sigma's included.

### What frmtmb reports now

`conv_headroom()` holds a covariance parameter at the edge of its range
as the SE check holds it (`headroom_edge()`: se_edge_cor() and
se_edge_sd()). The gain at the stop is then 8.1e-8, the stop is at the
maximum, and the SE check judges it as a converged fit: the boundary
message (class `frmtmb_boundary_fit`) names `cor_visit__Intercept__Trt1`,
the flat warning names the alias, and the sigma coefficients keep SEs
of 0.2848 and 0.4831. The code-8 warning stays (it says the stop is at
the maximum), because SEs are lost; the gain test does not take this
boundary stop for an interior maximum.

`dev/code8-fix1-draws.R` (count times 1 + 2e-15 j, j = -10 to 10),
`fix1-draws-after.log`: `h326 TRUE 21 of 21 ; grad warning 0 ; short
warning 0 ; boundary message 21`. The stop does not move under these
draws (theta_3 -723.6169 on all 21; count is an integer, so the
perturbation is 4e-14).

Not changed: `theta_4` (the smooth's sd at exp(-354.9)) stays
"nonfinite". se_at_edge() returns FALSE for it (the likelihood moves
when it goes further down from there, likely underflow), so it is not
called a boundary. Its row is the only part of the report that does not
say "boundary"; it is not a convergence matter.

### Ledger row brmsfit-methods:326

`dev/code8-fix1-jitter.R <arm> 40 1` (default start plus 40 jittered
starts, sd 1, seed 20261008):

- before (rellib-r19): `h326 TRUE: 25 of 41`; the default start fails
  (40 NaN errors), and so do 4 more starts that end in the default's
  basin (109.34007 to 109.34010).
- after (lane): `h326 TRUE: 34 of 41`. Every start that ends in the
  default's basin (30 of them, 109.3400 to 109.3402) holds, and so do
  the 3 that end at 105.4067 (nu run to the gaussian limit). The 7 that
  fail end elsewhere: 106.62, 107.95 and 109.08 (no parameter at an
  edge), 110.43 (both sds at 0 and the ARMA at its unit root), 106.996,
  109.41 and 105.3995 (nu at the Cauchy limit).

The row is asserted from the default start, which now holds on every
draw. It moves from "cannot transfer" to a pass: its line is removed
from `dev/brmsport-verdicts-manual-fit.tsv`, as 0.69.0's
`dev/rel069-verdicts.R` did when it last held. The generated
`dev/brmsport-verdicts.tsv`, the ledgers and
`tests/testthat/test-brms-suite-methods.R` (generated) are not edited:
until the consolidation regenerates them, the gated
`test-brms-suite-methods.R` fails at :326 with "now HOLDS but is
recorded as 'cannot transfer'", which is the expected signal.

## Part C: a stable SE verdict on the ridge of c0 + exp(a)^k

### The tally on rellib-r19

`dev/code8-ridge.R before trace` (y times 1 + 2e-15 j, j = -20 to 20;
`ridge-before.tsv`). The objective is 26.95097020 on all 41 (the ridge
is exactly flat: the invariance (a0, k, theta_1) -> (s a0, k / s,
theta_1 + log s)). a ranges from -29.6 to 25.4, k from -0.197 to 0.109.

    ridge before : 41 draws; se(a) finite 9
       3 : c0_(Intercept):flat a_(Intercept):flat k_(Intercept):flat sigma_(Intercept):flat theta_1:flat
       29 : c0_(Intercept):flat a_(Intercept):flat k_(Intercept):flat theta_1:flat
       9 : c0_(Intercept):nocurv k_(Intercept):nocurv theta_1:nocurv

The brief's "6 of 41" was CI's count; on this machine it is 9 of 41
(the draws are a few ulps apart and the count follows the platform).

### What decided each draw (three mechanisms)

1. The 9 draws with a finite SE for a (9.17 to 9.31) carried a gradient
   warning ("Large maximum absolute gradient", 0.0012 to 0.0202 at k).
   The headroom's quadratic model sees the gradient across the curved
   ridge as a gain along it, or cannot factor the Hessian; the fit was
   called short and kept sdreport()'s inverse. Every one of these draws
   is at the maximum (objective 26.95097020).
2. The 3 draws that lost sigma: at a = -26 (k = -0.067) the
   finite-difference step of 1e-3 on k raises the objective by 1.7, far
   beyond the quadratic region, and the asymmetry of the differences put
   a noise of 17.2 in sigma's row, whose largest entry is 144.5; the
   empty-row rule (10 times the noise) called the row empty
   (`dev/code8-ridge-trace2.R`, `ridge-trace2-before.log`).
3. The rest: the two smallest unit-diagonal eigenvalues (about +4e-5 and
   -4e-5) sit far inside the Hessian's noise (1e-3), both removed, and
   the flat subspace spans c0's direction too.

### The truth

`dev/code8-ridge-truth.R` (`ridge-truth.log`). With k held at its
estimate the ridge is gone and nothing else changes: objective
26.9509704 against 26.9509702, SEs c0 69.29, c0_x 0.0331, a 120.0,
sigma 0.0833, theta_1 12.07. The profile of c0 in the full model
(others re-optimized) loses 1.08e-6, 4.39e-6 and 1.87e-5 at steps of
0.1, 0.2 and 0.4: quadratic, with the fixed-k SE. So the ridge reaches
exactly a, k and theta_1; c0 is weakly but really identified (SE 69 on
an estimate of -3.4), and sigma and c0_x are well identified.

### The fix

- `fd_hessian()` keeps a retaken column whenever its rise was above the
  band (it used to keep the retake only where the retaken matrix
  factored, which a ridge never does). With the retake the unit-diagonal
  spectrum is the same on every draw: 3.00, 1.06, 1.00, 0.940, 8.8e-7,
  1.5e-8, with noise about 2e-7 (`ridge-tau.log`).
- `se_ridge_ritz()`: the small directions that load on nonlinear
  parameters are measured on the ridge's floor (Newton steps on the
  larger directions, then central gradient differences at 1e-3 and 1e-2,
  extrapolated) and re-diagonalized. On draws 0, -4, -7
  (`dev/code8-ridge-rr.R`, `ridge-rr.log`): eigenvalues 8.86e-7 and
  8.8e-11, 8.86e-7 and 7.5e-11, 8.85e-7 and -1e-12; the flat vector
  matches the invariance's tangent (0, 0, -0.7067, 0.7073, 0, 0.0186) to
  1e-6, where the Hessian's own was 7e-3 off on c0. The noise of the
  measurement is 7.3e-11, so the subspace bound tau is small and
  theta_1 (loading 0.0186) is named; c0 (loading 1e-6) is not.
- The gradient verdict and `conv_headroom()` (`headroom_checked()`):
  WITHDRAWN as first built (`headroom_probe()`, the objective along the
  Newton step and down the gradient in place of any prediction above
  `grad_tol`); see "Fix round" below. Now the prediction stands unless
  the Ritz step measured a direction flat; then the quadratic model is
  taken over the other directions and the objective is asked along the
  flat ones only (`headroom_off_flat()`). On the ridge every such probe
  gains 0.

### The tally on the lane

`dev/code8-ridge-exp.R none after4` (`ridge-exp-after4.tsv`):

    41 : a_(Intercept):flat k_(Intercept):flat theta_1:flat
    41 : W[frmtmb_warning] SE-NA 3 of 6 parameters.

c0's SE 69.09 to 76.99 on the 41 draws (truth 69.29), sigma's 0.083332
on every draw (truth 0.0833), no gradient or convergence warning.

### Under OpenBLAS (LAPACK routed too)

`dev/ciharden-openblas.sh 0.3.26 lapack` and `0.3.32 lapack`, then
`dev/code8-ridge-exp.R none ob<ver>` with that R copy:
`41 : a_(Intercept):flat k_(Intercept):flat theta_1:flat` on both, c0's
SE 70.49 to 70.68.

## Tests seen to fail on rellib-r19

`dev/release/run-tests.R` with no `FRMTMB_RELLIB` (`base-test-se-*.txt`):

- test-se-conv.R `pass=23 fail=4 err=2`: "false convergence (8) at the
  maximum is a message of its class" (a warning, no message) and "a
  stop at a correlation's limit is measured as the maximum" (the
  "curvature at the stop" warning, no boundary message, NaN sigma SEs).
- test-se-tier.R `pass=96 fail=8`: the rewritten ridge test fails on
  each of its three draws (c0 lost on draw 0; on draw -19 the gradient
  warning, a finite SE for a, c0 lost; on draw -4 sigma and c0 lost).
- Under OpenBLAS 0.3.26 the base fails the ridge test 6 times
  (`ob-0326-before.log`).

## OpenBLAS runs of the changed test files (lane)

`bash dev/code8-ob-run.sh <ver> after ...` (`ob-after.log`), 0.3.26 and
0.3.32 alike: test-se-reparam.R 79/0/0/0 warn 0, test-se-tier.R 104/0,
test-se-conv.R 41/0, test-design-offset.R 80/0. The CI WARN at
test-se-reparam.R:554 (code 8 at the maximum on ubuntu) did not appear
on these OpenBLAS copies for the base either (`ob-0326-before.log`,
warn 0), so the WARN itself is not reproduced on this machine; the
code-8 stop at the maximum that it reports is now a message when every
SE is finite. `dev/code8-aexp.R` (the fit at :554): on both builds,
reference BLAS and OpenBLAS 0.3.26, it stops at code 0 ("relative
convergence (4)") with a gain of 3.6e-9 (reference) or 2.6e-8
(OpenBLAS) and no SE lost, so on ubuntu's code-8 stop it should give
the message; that needs the CI run to confirm.

## Other verdicts that moved (core suite, one file per process)

Each was run on both builds with a script under `dev/`:

- test-mo-coding.R, `(mo(x) | g)` with x in 1, 2, 5, seed 23
  (`dev/code8-mogap.R`): before, "singular convergence (7)", "the
  curvature at the stop is unusable", every simplex SE and the
  correlation's NaN ("nocurv"), and get_vcov() NaN on `zeta1_cum1`.
  After, the correlation (theta_3 at 549) is held at its edge, the gain
  is 5.7e-11, the boundary message names `cor_g__Intercept__mox`, the
  simplex coordinates are "flat" (the gap), and get_vcov() returns the
  determined directions, `zeta1_cum1` included. The same mechanism as
  Part B on a code-7 stop. The test now asserts that, and keeps its old
  claim on the same model stopped short by `iter.max = 3`
  (`dev/code8-mogap2.R`: the "lost another direction" branch on both
  builds).
- test-mm-id-cov.R, a merged `mm()` block at the edge of its
  correlation space: the SE check now runs there and names theta_5 and
  theta_6 flat; the test allows that warning (its point is the design).
- test-brms-names-specials.R, a tight sd prior: the gradient warning
  ("a correlation's gradient short of the tolerance", in the test's
  own words) is gone where the objective gains nothing from the stop,
  and VarCorr() then warns that 16 of 22 entries move along a lost
  direction; the test allows that warning (its point is the names).
- test-bcm-mpt.R: see Part C's history; its gradient warning stays.
  With the retaken differences the group means lose their SEs too (the
  objective is spiked there); the test now asserts that the covariance
  block is among the lost parameters, not that it is the only one.
- test-se-exact-hessian.R: its reference "0.71.0's differences" is now
  spelled `fd_hessian(tau_band = Inf)` (no retake); and asym_laplace's
  kink fits, where the probe is not asked (gradient 0.66 to 4.5 in
  unit-diagonal units), are unchanged, with the SE warning allowed for
  the draws where the probe is asked.
- test-design-offset.R: four `expect_no_message()` sites around fits
  that stop at code 8 at the maximum wrap them in `allow_at_max()`.
- test-latres.R's cancelling ARMA roots: unchanged in the end. A
  version of the Ritz step that also ran on a single small direction
  read a curvature of 2.7e-4 on its floor and gave ar and ma SEs of
  1.1e4 and 147, so the limit message ("Their roots cancel") was lost
  (`dev/code8-latres.R`); the step now runs on two or more small
  directions only.

## Final suite runs (lane library, one file per process, NOT_CRAN=true)

- Core, `bash dev/release/run-par.sh code8-core2 ... 2`
  (`dev/release/code8-core2.log`): ran 317 of 317, every `lib:` line
  the lane library; pass 23282, fail 1, error 0, skip 240, warn 0. The
  one failure was test-mo-coding.R's short-stop half, which asserted the
  coefficient block of a stop short of the maximum finite; that
  assertion moved to the converged fit, and the file reran 166/0/0/0
  warn 0 (`dev/code8-out/rerun-mo-coding.txt`). test-se-tier.R 104/0
  (the base's known red at :323/:327 is the rewritten test).
- Extensions, `code8-ext2` (`dev/release/code8-ext2.log`): ran 192 of
  192, pass 9106, fail 0, error 0, skip 50, warn 0; every extension from
  rellib-r19 with the lane's core.
- Gated: test-brms-suite-methods.R with `FRMTMB_BRMS_FIT_TESTS=true`
  (`gated-brms-suite-methods.txt`): pass 161, fail 1, the expected
  "brms brmsfit-methods:326 now HOLDS but is recorded as 'cannot
  transfer'", which the consolidation's ledger regeneration clears.
- Not run: the three-way one-session runs (`run-sessions.sh`). A core
  session is about six hours of one process on this machine; R CMD
  check below runs the core suite in one session in alphabetical order.
  The change registers nothing and the new helper `allow_at_max()`
  leaves no state.

### Development history that the record must keep

- With the retake alone (`dev/code8-ridge-exp.R retake`, on rellib-r19
  by assigning the function): 37 of 41 lost a and k only, theta_1 kept
  an SE (the flat vector's theta loading 0.019 was under the capped
  threshold, since c0's real direction sat 8.6e-7 above it), and 4 still
  took the gradient-warning path.
- A first version of the Ritz step finalized "flat" on its own error;
  it called a raw polynomial's weakest real direction flat
  (test-se-reparam.R, 3 of 10 SEs lost) and named the correlated pair
  of test-se-flags.R. Now it only finalizes what it keeps; the rest goes
  to the step test (se_curvature_real()) and the re-measurement as
  before, with the vector it measured.
- A first version of the probe had no gradient condition. On
  test-bcm-mpt.R's hierarchical MPT (gradient 6.5e7, 2800 in
  unit-diagonal units) it found no gain, but the objective there has a
  spike: it falls by 0.015 within 1e-9 of the stop along one coordinate
  and rises like the logarithm of the distance beyond
  (`dev/code8-mpt.R`, `mpt-after3.log`). The gradient condition keeps
  that fit's gradient warning. On the binomial Michaelis-Menten ridge of
  test-nl-start.R the Ritz step at first took a direction of negative
  curvature (an ascent to K = 1e8) and gave finite SEs; it now asks only
  directions within their noise of zero or above.

## Fix round (review dev/reviews/2026-11-10-code8.md: B1, M1, m1 to m3)

### B1: the probe no longer certifies a maximum

The first build's `headroom_checked()` replaced any quadratic prediction
above `grad_tol` by the best of a few straight probes, a lower bound on
what is left; on a curved valley (y = A exp(B x), x on [30, 31] or
[100, 101]) it called stops 0.003 to 10.5 short the maximum, three with
the classed message. Withdrawn. The prediction (`newton_room()` or the
Cholesky one) now stands, except where the SE check's Ritz step measured
a direction flat on the floor (`ritz_flat`, new in se_tier3()'s result):
then the Newton gain is taken over the other directions and the
objective is asked only along those flat ones (`headroom_off_flat()`).

- The directions the model makes exactly flat (se_design_null(): an
  alias, a mo() simplex split by a category no row is in) join the
  Ritz-flat ones. Found by the fix build's suite run: test-mo-coding.R's
  gap fit (seed 23) has two such directions whose finite-difference
  eigenvalues read -4.1e-6 and -5.9e-6 (gradient along them 1e-11), so
  `newton_room()` called the curvature unusable and the fit was taken
  as short again (`dev/code8-mogap3.R`). With them taken out the gain is
  1.1e-9 and the file passes 166/0. On that build the reviewer's
  valley scripts give the same zeros as below (`valley8-fix2.log`,
  `valley-opt-fix2.log`) and the ridge tally is 41 of 41
  (`ridge-exp-after7.log`).
- BLOCKER from the reviewer's check of that change
  (`dev/review-code8-fix-dnull.R`, `-dnull2.R`): a design null is flat
  in the likelihood only. y ~ x + s(x), gaussian, n 200, seed 5,
  normal(0, psd) on x and s(x).fx1, stops stepped along the null: at psd
  1e4 code 8 gave the message at 0.00117 and 0.0187 short and "at the
  maximum" up to 11.7 short; code 0 was silent where the base warns.
  Two changes:
  1. The design-null directions are taken only where the fit has no
     prior with a density (`fit_prior_density()`); with one, the
     quadratic model, which has the prior's curvature, decides.
  2. The message's bound now also counts what the directions too flat
     for the quadratic model could hold (`newton_room()` and
     `headroom_off_flat()` return it as attribute `flat_bound`, the sum
     of g^2 / (2 bound) over them). Without it, at psd 1e4 the prior's
     curvature along the alias (1e-8 of the largest) is under
     `se_flat_tol` and the stops 0.00117 and 0.0187 short read a gain of
     0 and still got the message.
  On the final build `-dnull2.R` gives the base's verdicts exactly (the
  columns says_max, grad_w, code_w, at_max_m and nwarn on all 20 stops
  identical, `dnull2-review-fix4.log` against `dnull2-review-base.log`);
  `-dnull.R` (psd 1 to 100): "short > 1e-3: 28 ; called at max: 0".
  The two code-8 stops still called "at the maximum" in the warning
  text (0.00117 and 0.0187 short at psd 1e4) and the code-0 stops up to
  0.117 short with no warning are the base's own verdicts: the quadratic
  model's allowance of up to `grad_tol` along directions below
  `se_flat_tol` (pre-existing, not changed here). New test in
  test-se-conv.R, "an alias under a prior is not taken as flat" (psd
  1e4, a step of 2e4 along the null, 1.87 short, codes 8 and 0): fails
  twice on the build before this change (`dnull-before-fix.log`),
  passes after; test-mo-coding.R still passes 166/0
  (`rerun3-mo-coding.txt`).
- Reviewer's scripts on the fix build: `dev/review-code8-valley8.R`
  (288 stops) "short > 2e-5 with the converged_at_max message: 0 of
  286; short > 1e-3 whose warning says 'at the maximum': 0; short >
  1e-3 with no warning: 0" (`dev/code8-out/valley8-fix.log`);
  `dev/review-code8-valley-opt.R` "short by > 1e-3 with no warning: 0
  of 23" (`valley-opt-fix.log`). The x0 = 100 rows where optim()'s BFGS
  or CG returns code 0 while 13.8 to 225 short warn only about lost SEs,
  as on base; an optimizer that reports success is not questioned
  (pre-existing, not this lane's).
- New test, test-se-conv.R "a stop short on a curved valley is not
  called the maximum" (forced code 8 from BFGS cut at 40 iterations,
  seed 1; BFGS's own code-1 stop, seed 2): fails 3 times with the first
  build's `headroom_checked()` assigned into the namespace
  (`dev/code8-mutant-tests.R oldprobe`, `mutant-oldprobe-seconv.log`),
  passes on the fix build.

### Batteries rerun on the fix build

`dev/code8-summarize.R` and `dev/code8-short-bound.R after3 after-x3`
(generated):

    battery after-base2: 99 fits, 32 at code 8; gain 1.95e-09 to 9.92e-06;
      ... warned (conv_code) 0; message 32; with another warning 3
    battery after-ywx2: 220 fits, 77 at code 8; gain 2.42e-11 to 1.98e-06;
      ... warned (conv_code) 0; message 77
    code-8 stops 99 ; short by more than 2e-5: 80 warned 80 message 0
    short by 1e-6 to 2e-5: 13 message 7 warned 6 ; their gain/shortfall
      ratio 0.9995 1.0002

No fit at a true maximum now keeps the warning (109 of 109 give the
message, as before). The 7 messages are stops 3.6e-6 to 1.95e-5 short,
whose stated gain is that shortfall. asym_laplace's kink stops no
longer read a gain at all (the probe that read 0 there is gone); they
keep the warning, as on base.

### M1: a wider block's correlation at its limit

`(1 + poly(x, 2) | g)` with `normal(2, 0.001)` on one sd: the third
Cholesky row stops at (808, 1814), the correlations of `polyx22` with
the other two at 0.41 and 0.994. Moving one entry alone changes the
log-likelihood (theta_5 by +1: -3.6e-5), so the single-parameter edge
test said no; moving the row as a whole does not (doubling: -8.4e-8;
x10: -5.3e-7; pulled back to the clamp: +8.0e-5;
`dev/code8-tightsd-pull.R`). At the stop the finite-difference
Hessian's asymmetry is 254 in the theta rows and 3.3 in the intercept's
row (`dev/code8-tightsd.R`), and every eigen-direction of the kept block
was within three times its noise (eigenvalues 2.36 to -0.51, noise 0.07
to 1.27), so tier 3 removed them all and named the fixed effects. At
the row pulled back the asymmetry is 4.5e-4.

`se_edge_cor()` now treats a ROW of a block of three or more as at its
edge when its largest entry is past the clamp and doubling the row
changes the log-likelihood by at most `se_edge_tol`
(`se_row_at_edge()`); it is pulled back along itself (`se_pull_back()`
reads a `scale` attribute). `conv_headroom()` and `grad_headroom()` hold
it and read the pulled Hessian too, so this fit has no gradient warning.
On the fix build only theta_5 and theta_6 are lost
(`dev/review-code8-tightsd.R`, `-tightsd2.R`; `tightsd2-fix.log`):

| parameter | lane SE | reviewer's Richardson SE at the stop |
|---|---|---|
| (Intercept) | 0.283 | 0.2846 |
| poly(x, 2)1 | 1.464 | 1.453 |
| poly(x, 2)2 | 1.333 | 1.297 |
| sigma | 0.0522 | 0.0535 |
| theta_1 | 0.241 | 0.2145 |
| theta_2 | 0.497 | 0.01716 |
| theta_4 | 0.660 | 0.2027 |

theta_2 and theta_4: which is right? The profile decides
(`dev/code8-tightsd-prof.R`, every other parameter re-optimized,
theta_5 and theta_6 free; `tightsd-prof.log`). A quadratic profile
loses 0.5 at one SE:

    theta_2 clamp      SE 0.497:  loss at -1, +1 SE: 0.332 0.670
    theta_2 richardson SE 0.0172: loss at -1, +1 SE: 5.35e-04 6.92e-04
    theta_4 clamp      SE 0.66:   loss at -1, +1 SE: 0.717 0.219
    theta_4 richardson SE 0.203:  loss at -1, +1 SE: 0.0556 0.0389
    theta_1 clamp      SE 0.241:  loss at -1, +1 SE: 0.550 0.448

The values at the clamp are the profile's (theta_2 averages 0.50 at one
SE; theta_4 0.47, skewed); the Richardson references are 30 and 3 times
too small, because that Hessian is itself noise there (asymmetry 3843,
eigenvalues -5.55 and -15.5 whose probes rise both ways). So the clamp
values are reported as finite SEs and need no further qualifier: a user
sees `sd(polyx21)` and `cor(Intercept, polyx21)` with these SEs, and
the two correlations of `polyx22` without one ("Standard errors are not
available for 2 of 10 parameters ... theta_5, theta_6").

The comment in test-brms-names-specials.R now says what is lost and why.

Known limit (fix-round short check, note 1; filed in
`dev/test-backlog.md`, "Filed in the round of 2026-11-09"): a row held
at its edge loses BOTH entries' SEs, theta_5 and theta_6 here, although
the row's direction is determined (it carries cor(Intercept, polyx22) =
0.41); only its length, the third component's own share, is at the
boundary.

The weak-direction set of the review (`dev/code8-weak.R`, a copy of
`dev/review-code8-weak.R`, 120 fits): lost sets identical to the
review's lane and base runs on all 120 (`weak-fix.tsv` against
`dev/review-code8-out/weak-after.tsv` and `weak-before.tsv`).

### m1: the retake's cost

Kept, because it is what removes the noise in sigma's row on the ridge
(3 of 41 draws lost sigma without it). Measured by the reviewer
(`dev/review-code8-cost.R`): gradients for the SE Hessian 16 / 20 (LMM
1e5 rows, 8 parameters), 8 / 12 (poisson GLMM 1e5, 4), 10 / 16 (nl 2e4,
5); seconds 1.25 / 1.52, 1.14 / 1.65, 0.09 / 0.20; SEs move by at most
3e-7 relative. `se_check_at_fit()` now budgets the worst case, 4 np
gradients, against `se_check_share` of the optimizer's evaluations.
Which test fits that moves is logged in the suite runs below.

### m2, m3: the ridge test

Draws 0, -20, -3, 6, 17: the mutant without the ratio test
(`dev/code8-ridge-mutant.R`, a copy of the review's script) keeps a
finite SE for a on 12 of 41 draws (and loses "nocurv" on j = -11), and
fails the test 10 times on these 5 draws (`mutant-noratio-setier.log`);
the fix build passes. Gradient and convergence warnings are allowed;
the assertions are the lost set, its reason "flat", NaN for the three
and finite SEs for the rest.

### Ridge, fixture 1 on the fix build

- Ridge (`ridge-exp-after6.tsv`): 41 of 41 lose exactly a, k, theta_1;
  c0's SE 70.48 to 70.68, sigma's 0.083332 on every draw; no gradient
  or convergence warning.
- Fixture 1 (`fix1-draws-after2.log`): `h326 TRUE 21 of 21 ; grad
  warning 0 ; short warning 0 ; boundary message 21`.

### Final build (design-null prior gate and `flat_bound`), rechecks

Generated from the logs named:

- At-max battery (`battery-after-base3.tsv`, `battery-after-ywx3.tsv`):
  32 of 32 and 77 of 77 code-8 stops give the message, none warned.
- Short stops (`dev/code8-short-bound.R after4 after-x4`): "code-8
  stops 99 ; short by more than 2e-5: 80 warned 80 message 0 ; short by
  1e-6 to 2e-5: 13 message 7 warned 6 ; their gain/shortfall ratio
  0.9995 1.0002".
- Ridge (`ridge-exp-after8.log`): 41 of 41 lose exactly a, k, theta_1;
  c0's SE 70.48 to 70.68; no condition but the lost-SE warning.
- Fixture 1 (`fix1-draws-after3.log`): "h326 TRUE 21 of 21 ; grad
  warning 0 ; short warning 0 ; boundary message 21".
- Reviewer's valley8 (`valley8-fix4.log`): "short > 2e-5 with the
  converged_at_max message: 0 of 286 ; short > 1e-3 whose warning says
  'at the maximum': 0 ; short > 1e-3 with no warning: 0".
- Design null (`dnull2-review-fix4.log`): verdicts identical to the
  base build on all 20 stops; `dnull-review-fix3.log`: 0 of 28 called
  the maximum.

Not run by this lane, at the coordinator's instruction: the full core
and extension suites, the OpenBLAS runs of the SE files and R CMD check
on this build; a separate verification agent runs them.

## Final verification

Run on 2026-11-10 by the verification agent; no package code or test
changed.

### Build

The installed frmtmb in `wt-code8-lib` carries the worktree's code:
all 83 functions of `R/fit.R` and all 69 of `R/se-check.R` deparse the
same as the installed namespace (`dev/code8-verify-build.R`), among
them `headroom_ritz_flat`, `newton_room` and `headroom_off_flat`. No
reinstall. The extensions come from rellib-r19; the lane changes none.

### Suites

`dev/release/run-par.sh code8-final dev/release/code8-final.jobs 3`,
one file per process, `NOT_CRAN=true`,
`FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-code8-lib`,
`R_PROFILE_USER=dev/code8-budget-profile.R`; log
`dev/release/code8-final.log`, per-file logs
`dev/release/code8-final-files/`. 509 of 509 ran; every `lib:` line
names `frmtmb 0.80.0 C:/Users/adf44/source/r/wt-code8-lib`.

| package | files | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 317 | 23305 | 1 | 0 | 239 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1848 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 687 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 86 | 3967 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |

The known-red base test (test-se-tier.R:323, :327) passed, and
test-se-reparam.R had no WARN.

**The one failure is a lane regression.** test-se-check.R:472, "whether
the check waits is decided by counted work":
`expect_true(at_fit(50, 400))` gets FALSE. The test pins the 0.80.0
budget, 2 np <= 0.25 evals (100 <= 100); the lane's
`se_check_at_fit()` asks 4 np <= 0.25 evals (200 > 100), so the edge is
now 800 evaluations at np 50. The lane did not update the test. It is
deterministic: alone, without the budget profile, the lane build fails
it the same way (`dev/code8-out/final/lane-se-check.txt`: pass=133
fail=1), and rellib-r19 passes (`dev/code8-out/final/base-se-check.txt`:
pass=134 fail=0). Not fixed here.

### Fits whose fit-time SE check the 4 np budget defers

`dev/code8-out/budget-defer.log`, this run only (the earlier partial
log is kept as `budget-defer-prefinal.log`). Each line is a call that
0.80.0's 2 np rule checks at fit time and the lane's 4 np rule defers.
24 fits in 15 files; the 25th line (test-se-check.R, np 50, evals 400)
is the synthetic call of the failing test above, not a fit.

| file | np | evaluations |
|---|---|---|
| test-autocor.R | 11 | 94 |
| test-brms-port.R | 11 | 123 |
| test-curve.R | 13 | 159 |
| test-double-bar.R | 13 | 113 |
| test-double-bar.R | 12 | 135 |
| test-group-labels.R | 11 | 107 |
| test-me.R | 11 | 130 |
| test-me.R | 11 | 102 |
| test-optstart.R | 13 | 179 |
| test-optstart.R | 13 | 189 |
| test-optstart.R | 13 | 190 |
| test-se-check.R (synthetic) | 50 | 400 |
| test-se-flags.R | 11 | 104 |
| test-se-flags.R | 13 | 199 |
| test-se-flags.R | 13 | 151 |
| test-se-flags.R | 13 | 180 |
| test-se-flags.R | 13 | 125 |
| test-se-reparam.R | 11 | 98 |
| test-se-reparam.R | 11 | 127 |
| test-smooths.R | 16 | 212 |
| test-sparsear1.R | 12 | 108 |
| test-sparsex.R | 12 | 151 |
| test-sparsex.R | 12 | 151 |
| frmtmb.coupling test-coherence.R | 12 | 105 |
| frmtmb.latent test-hmm.R | 13 | 144 |

### OpenBLAS

`bash dev/code8-ob-run.sh <ver> after test-se-reparam.R test-se-tier.R
test-se-conv.R test-design-offset.R`, OpenBLAS with LAPACK, lane build
(logs `dev/code8-out/ob/`). Identical counts on both versions:

| file | 0.3.26 | 0.3.32 |
|---|---|---|
| test-se-reparam.R | pass=79, 0 fail/err/skip/warn | pass=79, 0 fail/err/skip/warn |
| test-se-tier.R | pass=112, 0 fail/err/skip/warn | pass=112, 0 fail/err/skip/warn |
| test-se-conv.R | pass=55, 0 fail/err/skip/warn | pass=55, 0 fail/err/skip/warn |
| test-design-offset.R | pass=80, 0 fail/err/skip/warn | pass=80, 0 fail/err/skip/warn |

`dev/ciharden-out/` deleted afterwards.

### R CMD check --as-cran (core)

`bash dev/code8-check.sh`, `dev/code8-check/check.log`:
`Status: 1 ERROR, 1 NOTE`.

- ERROR, checking tests: the same failure,
  `test-se-check.R:472:3` "Expected `at_fit(50, 400)` to be TRUE";
  `[ FAIL 1 | WARN 0 | SKIP 505 | PASS 18962 ]`.
- NOTE, checking HTML version of manual: "Skipping checking math
  rendering: package 'V8' unavailable" (environment, not the lane).

### Reviewer's last MINOR (accepted)

From `dev/reviews/2026-11-10-code8.md`, third fix-round check: through
the Ritz-flat exception, `c0 + exp(a)^k` with `normal(0.5, psd)` on k,
at psd 1000, ridge stops 1.9e-4 and 1.23e-3 short lose the gradient
warning the base gives. From 0.00316 short upward every stop warns on
both builds, and no stop gets the converged-at-maximum message. The
difference sits at the `grad_tol` scale (1e-3)
(`dev/review-code8-fix-ridgeprior.R`, `-ridgeprior2.R`). Filed in
`dev/test-backlog.md`.

### The one regression, fixed by the coordinator (2026-11-10)

`tests/testthat/test-se-check.R:472` ("whether the check waits is
decided by counted work") pinned the 0.80.0 budget, 2 np <= 0.25 evals:
`at_fit(50, 399)` FALSE and `at_fit(50, 400)` TRUE. Under the lane's
4 np budget the edge at np = 50 is 800 evaluations. The test now asserts
`at_fit(50, 799)` FALSE and `at_fit(50, 800)` TRUE, with a comment
giving the arithmetic. No package code changed.

- `test-se-check.R` on the lane build: pass 134, fail 0, err 0, skip 0,
  warn 0 (lib: wt-code8-lib frmtmb 0.80.0).
- `bash dev/code8-check.sh` rerun: `Status: 1 NOTE` (V8 math rendering,
  environment).
