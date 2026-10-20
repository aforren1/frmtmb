# Lane adhess: findings

Round of 2026-10-14, base 47cf659d (frmtmb 0.71.0). Defect filed by
lane illse2 (`dev/illse2-findings.md` in its tree, "Defects found, not
fixed") and measured by its reviewer (`dev/reviews/2026-10-14-illse2.md`
in that tree, findings 7 and 8). Every number below comes from a log in
`dev/adhess-log/` (local, gitignored), written by the script named
beside it. "base" is rellib-r9 alone; "lane" is
`C:/Users/adf44/source/r/wt-adhess-lib` over rellib-r9.

STATUS: DONE. Every item of the lane's test plan ran on the install of
2026-10-20 (see "Resumed, second part" and sections 7 to 9). The
subsections "What is not done", "Changed on resume, not run" and "The
blocker" below record the state on the morning of 2026-10-19 and are
superseded.

## Resumed, 2026-10-19

The first agent was stopped by the round's pause on 2026-10-14. A
second agent resumed on 2026-10-19 and was blocked before it could run
R. The text of this section up to "Resumed, second part" is as written
then.

### What was done before the pause

Established from the code, the diff and the logs in `dev/adhess-log/`.

- The change (sections 3 and 4) in R/se-check.R and R/autoscale.R, the
  coupling refusal in extensions/frmtmb.coupling/R/coupling.R, the
  NEWS entries, the vignette paragraph, the two re-pinned tests and
  the new test-se-exact-hessian.R. Core was last installed into the
  lane library at 2026-10-14 16:13 (install-5.txt) and no file under
  `R/` changed after that, so every log from 16:14 on is of the
  current core. Coupling was installed at 16:37 (install-coupling.txt).
- Sections 1 to 6: measured and recorded, base and lane.
- The fragility scan, reference BLAS, all eight suites one file per
  process (`dev/adhess-scan.sh`, `scan/lane-ref.sum`, started 16:14).
  Generated from the summary: 380 files with a RESULT line, pass 27172
  fail 2 err 0 skip 0 warn 0, and the driver's own line "ran 382 of
  383". The two failures are stale: test-se-tier.R ran at 16:18, before
  its re-pin at 16:26 (the failure is the old `RTMB::sdreport(obj)`
  pin, line 104), and frmtmb.coupling's test-coherence.R ran at 16:19
  with rellib-r9's coupling, before the coupling change. Both were
  rerun on the lane library at 16:38 (`t2-*.txt`): test-se-tier.R pass
  91, test-coherence.R pass 66, test-se-check.R pass 129,
  test-se-exact-hessian.R pass 25, each fail 0 err 0 skip 0 warn 0. The
  missing file is frmtmb.ode's test-scale.R: no RESULT line after
  1896 s, killed by the pause. The other frmtmb.coupling test files
  were not rerun with the changed coupling.
- The SE false-alarm battery (lane illse's five groups, 41 models,
  seeds 1 to 20) and lane illse's bitwise set, lane against rellib-r9,
  on all three BLAS (`dev/adhess-sf-run.sh`; `sf-*/falarm-cmp.txt`,
  `sf-*/bitcmp-*.txt`), generated from the logs:

      sf-ref:  fits 820 | rows differing in code 0, lost 0, se_warn 0,
               other_warn 0, max_rel (string) 377
      sf-ob26: fits 820 | rows differing in code 0, lost 0, se_warn 0,
               other_warn 0, max_rel (string) 384
      sf-ob32: fits 820 | rows differing in code 0, lost 0, se_warn 0,
               other_warn 0, max_rel (string) 384

  So no new warning and no lost SE on any BLAS. Of the 820 reference
  rows, `max_rel` (the largest |SE / ref - 1| against the battery's
  reference) moved on 377: lower on the lane on 305 and higher on 72
  (`sf-ref/sum.txt`). The higher rows are the poisson `pois_cor*`
  models (16 to 19 of 20 seeds each, by at most about 1e-7 in a value
  near 1e-6), `logit_uncentered` (10 of 20, below 1e-9 absolute), and
  three random-effects models (poly4_lmer 7, poly3_lmer 1,
  lmer_cor0.9999999 1). Not yet explained: the battery's reference for
  those models is not exact, so "farther" here may be the reference's
  own error. To check against an exact reference when R is back.

### What is not done

- Section 7 (cost): `dev/adhess-cost.R` is written and was never run.
- Section 8: the OpenBLAS 0.3.26 and 0.3.32 scans of the eight suites,
  frmtmb.ode's test-scale.R on every configuration, the coupling suite
  with the changed coupling, the one-session runs
  (`dev/adhess-onesession.R`, written, never run), and R CMD check of
  frmtmb and frmtmb.coupling.
- A fresh-install re-verification of sections 1 to 6.

### Changed on resume, not run

- test-se-exact-hessian.R carried absolute tolerances (1e-8, 1e-6,
  1e-5, `tolerance = 1e-12`), which lane-env.md forbids. Now: the two
  GLM tests bound the error by `1e3 * eps * kappa` of the reference's
  correlation matrix (the rounding of two inverses); the nonlinear test
  by 1e-3 times, and the GLMM test by 1e-2 times, the error that
  optimHess()'s step (fd_hessian() at its default) has on the same fit
  at the same point; the step tests by 64 eps. NOT RUN: see the
  blocker. The margins come from the logged values: on the lane the
  GLM errors are 1e-16 to 5e-14 (fdcov-lane2.txt, glm-lane.txt) against
  a bound of about 2e-11 for the poisson design (correlation -0.98,
  condition number about 100); nl 4.9e-7 against 0.133 times 1e-3;
  GLMM 2.6e-7 against 3.7e-3 times 1e-2.

### The blocker

The user library `C:/Users/adf44/AppData/Local/R/win-library/4.6` is
hollow: 223 of its 405 package directories hold no files (no
DESCRIPTION), all with a modification time of 2026-10-16 23:30, the
signature `dev/machine-library.md` describes. The hollow set includes
RTMB, TMB, Rcpp, RcppEigen, rlang, testthat, lme4, brms, rstan,
StanHeaders and numDeriv, so neither frmtmb nor roxygen2 loads
(`install-6.txt`, `install-7.txt`: "there is no package called
'rlang'"). The lane may not install into the user library, and a
private copy of 223 CRAN packages would put the lane on dependency
versions the round's records were not made with. The restore is the
user's: `dev/machine-library.md`, "Restoring it", from
`dev/win-library-4.6-manifest-2026-09-28.csv`.

## Resumed, second part (2026-10-19, after the restore)

The round's coordinator restored the user library in place (223 of 223
hollow packages, every comparable version that of the 2026-09-28
manifest). Core and frmtmb.coupling were roxygenised and reinstalled
into the lane library (`dev/adhess-install.sh frmtmb.coupling`,
`install-8.txt`, built 2026-10-20 01:34 UTC); roxygenise changed no
tracked file.

### Sections 1 to 6 re-verified on the fresh install

`dev/adhess-reverify.sh` reran every measurement script of sections 1
to 6 on both arms (logs in `rv/`). Each log, header line aside, is
identical to the one the section cites (`diff`): fdcov-base,
fdcov-lane (against fdcov-lane2.txt), glm-base, glm-lane, nls-base,
nls-lane, shift-base, shift-lane, routes-base, restep-lane (against
restep-4.txt) and movecmp. `rv/restep-base.txt` confirms the base's
lost set of section 4: pois_xbig_re, gamma_re and pois_xbig_n1e4 each
"lost flat 1" with the warning; every other model "lost none".

### The new test file, lane and base

`dev/ciharden-run1.R`, one process each (`rv/t3-exact-*.txt`):

    lane (wt-adhess-lib): RESULT test-se-exact-hessian.R pass=25 fail=0 err=0 skip=0 warn=0
    base (rellib-r9):     RESULT test-se-exact-hessian.R pass=5 fail=12 err=2 skip=0 warn=1

On the base every exactness bound fails behaviorally: the poisson GLM's
SE error is 6.17e-2 against a bound of 2.76e-11, the bernoulli's
6.70e-2 against 2.96e-13, the nonlinear body's 8.66e-3 against 8.66e-6
(the base's covariance IS the old route there, so the ratio is exactly
1e3), the AD-Hessian test finds `exact` FALSE and a Hessian that is not
the AD one, and the GLMM loses theta_1 as "flat" with the warning. The
two step tests error on the missing `outer_fd_steps()` (the weak form;
the unit-scale one is the guard's absent case and must hold on both).

Margins on the lane (`dev/adhess-tolmargin.R`, which runs the file with
each `expect_lt()` printing its value and limit;
`rv/tolmargin-lane.txt`):

    poisson GLM    se  1.577e-14 / 2.761e-11 = 5.7e-4
                   cor 3.147e-14 / 2.761e-11 = 1.1e-3
    bernoulli GLM  se  6.661e-16 / 2.956e-13 = 2.3e-3
                   cor 1.090e-15 / 2.956e-13 = 3.7e-3
    nonlinear      se  1.296e-12 / 8.664e-06 = 1.5e-7
                   cor 4.866e-07 / 1.330e-04 = 3.7e-3
    GLMM           se  2.736e-07 / 5.951e-04 = 4.6e-4
                   cor 5.473e-07 / 3.427e-03 = 1.6e-4

### The 72 battery rows where the lane looked farther

The battery's reference on those models is glm()'s or lmer()'s
`vcov()`, not an exact one. `dev/adhess-sfexact.R` reran the battery's
"lin" and "poly" groups with the battery's own loop (same data, seeds 1
to 20, both arms) and measured each fit's largest |SE / ref - 1|
against an exact reference: without random effects the inverse of
X' W X formed and inverted in 256 bits at frmtmb's estimate (pois_cor,
logit_uncentered) or the inverse of the AD Hessian (lm, poly lm); with
random effects a Richardson Jacobian of the AD gradient
(`dev/adhess-sfexact-cmp.R`, `sfx/cmp.txt`). 360 fits per arm; 0
differ in warnings or lost set. Comparing the unrounded values, the
lane is farther from the battery's reference on 132 rows (the battery's
72 compared values printed to 4 digits). Against the exact reference,
of those 132 the lane is closer on 107 and farther on 25:

    model              n  lane max   base max   lane closer  farther
    pois_cor0.999      20 4.28e-13   4.56e-07   20           0
    pois_cor0.99999    20 2.96e-11   5.47e-07   20           0
    pois_cor0.9999999  20 4.55e-09   4.23e-07   20           0
    lm_cor (3 models)  60 0          1.36e-07   60           0
    poly3_lm           20 7.70e-11   3.49e-09   20           0
    logit_uncentered   20 7.39e-10   7.72e-10   8            12
    lmer_cor (3)       60 1.47e-05   1.47e-05   8            12
    poly3_lmer         20 2.82e-07   2.82e-07   10           10
    poly4_lmer         20 3.97e-05   3.92e-05   7            11

The lm_cor and poly3_lm rows are an identity on the lane (its
covariance is the inverse of the AD Hessian, the reference used there);
section 3's 256-bit references are the independent check for that
route. So the pois_cor rows are glm()'s own error (glm() stops IRLS at a
relative deviance change of 1e-8 and evaluates its vcov() at that
estimate; its SEs are 1e-7 to 1e-5 from the exact ones), and the lane
is exact there to 4.6e-9 where the base was 2e-7 to 5e-7 off. The 25
rows where the lane is farther from the exact reference are noise:
logit_uncentered (x ~ N(1000, 1), a design whose correlation matrix
has a condition number near 1e6, so both arms sit at the rounding of a
double Hessian, at most 7.7e-10); and gaussian models with random
effects, where both arms are central differences, their errors (up to
4e-5) are equal to two or three digits, and the lane's step for a
column wider than one unit moves the rounding, not the truncation (the
gradient of a gaussian model is linear in its coefficients). poly6_lmer
has errors of 0.707 on both arms identically and its Richardson
reference is not invertible on 5 seeds: an ill-conditioned design that
is lane illse2's, not this lane's.

### Summary for the brief's item 2: which steps changed

Every finite-difference Hessian or covariance site whose step is not
scaled to its parameter is listed in section 2. Changed:

- The covariance of a fit without random effects: the AD Hessian
  (`outer_hessian_ad()` through `fit_outer_hessian()`), handed to
  sdreport() by `autoscale_sdreport()`. Measured best: exact to the
  rounding of a double Hessian (sections 3 and 7), and it costs no new
  Hessian, because the fit already holds it.
- The covariance of a fit with random effects: `fd_hessian()` takes a
  step per coordinate, `outer_fd_steps()`, scaled to each coefficient's
  design column. Chosen over the curvature retake and over Richardson
  by the measurements of section 4 (the retake moves an SE-check
  verdict; Richardson doubles the gradients for an accuracy nothing
  needs). The QR-coordinate route (`reparam_hessian()`) already scales
  its steps to the design and is unchanged.

Bumps and floors: frmtmb minor (every fit without random effects
changes its SEs; BREAKING in NEWS). frmtmb.coupling patch, with a
floor on this frmtmb: its new refusal answers the finite, huge SE the
exact Hessian gives at a coherence of 1, where 0.71.0 gave NaN and the
existing refusal fired. The other six extensions need no change (their
suites pass on this core unchanged, section 8).

Left, each with its reason in section 2: the nonlinear and variance
parameters with random effects (section 9), `newton_polish()`'s
optimHess() (tried, no measured effect), the headroom fallbacks, the
autoscale verdicts, and tmbroot()'s bracket.

## 1. The cause

On a fit without random effects the covariance came from
`autoscale_sdreport()`, which hands `RTMB::sdreport()` the Hessian of
`fit_outer_hessian()`. That was `fd_hessian()`: optimHess()'s central
differences of the AD gradient, column `i` as
`(g(p + eps e_i) - g(p - eps e_i)) / (2 eps)` with `eps = 0.001`
ABSOLUTE in the coordinates it is given (natural units, or `par / u`
where autoscale `par_units` are in play). `RTMB::sdreport()` itself
does the same when it gets no Hessian
(`hessian.fixed <- optimHess(par.fixed, obj$fn, obj$gr)`, read from
`TMB:::sdreport` as installed). `par_units` widen only a column whose
spread is below 1e-3 (autoscale), so a column of spread 100 keeps the
step 1e-3 and the step moves the predictor by 0.1 to 0.5 per row.

Why the error is large: the truncation of a central difference is
`h^2 / 6` times the third derivative of the gradient along the step.
For a GLM that is about `(h * x)^2 / 6` relative: with `x` near 500,
`h = 1e-3` gives `(0.5)^2 / 6 = 0.04`, which is the 6.2 percent SE error
measured below. A gaussian identity-link model is immune in its
coefficients because its gradient is linear in them.

The SE check (`se_fit_analysis()`) already read the exact AD Hessian
(`obj_he_cached()`) on a fit without random effects, found nothing, and
then the covariance was taken from the finite differences anyway
(`sdr_rescue_hessian()` keeps sdreport()'s inverse when it is usable).

With random effects there is no AD Hessian of the Laplace objective,
and `fd_hessian()`'s absolute step is the same defect. The retake in
`fd_hessian()` (a step that gives a rise of 1e-4) runs only when the
matrix has no Cholesky factor.

Reproduction (`dev/adhess-fdcov.R`, a copy of the reviewer's
`dev/illse2-rev-fdcov.R`; log `fdcov-base.txt`), reference the inverse
of the AD Hessian without random effects:

    model           base SE err   base cor err   warning
    exp_x20         8.66e-03      1.33e-01       -
    pois_xbig       6.17e-02      1.20e-01       -
    bern_xbig       6.74e-04      1.35e-03       -
    pois_xbig_re    NaN           NaN            lost, flat
    well scaled     4.6e-8 to 9.1e-7

## 2. Every finite-difference Hessian or covariance site

Inventory of core and the seven extensions (read-only sweep, then
checked by hand at each site). ABS is an absolute step.

Sites that build a covariance or a Hessian the covariance is read from:

- `fd_hessian()` (R/se-check.R), through `fit_outer_hessian()`: ABS 1e-3.
  CHANGED (sections 3 and 4).
- `autoscale_sdreport()` (R/autoscale.R): hands sdreport() the Hessian
  above; its fallbacks are sdreport()'s own optimHess() (ABS) when
  `fit_outer_hessian()` errors, and optimHess() in `par / u` when the
  kept Hessian is not at the fit's point. Reached only when
  `fit_outer_hessian()` fails. Not changed.
- `reparam_hessian()` (R/se-reparam.R, tier 4 and the Newton steps on
  an ill-conditioned design): differences along the design's whitened
  columns (scaled to the design), Richardson; ABS 1e-3 for the
  parameters outside the design blocks (log sd, sigma, nonlinear).
  Not changed: illse2's file, and a log sd's own scale is its unit.
- `newton_polish()` with random effects (R/fit.R, only when
  `predictor_design_kappa() > 1e5`): `optimHess(parscale = u)`, ABS in
  the units u. Tried with the new steps; it changed nothing measured
  (finding 8, section 6), so it was reverted. Not changed.
- `conv_headroom()`, `grad_headroom()` (R/fit.R): read
  `fit_outer_hessian_at()` with random effects, so they get the new
  steps; their fallbacks (optimHess, ABS) run only when that fails.
  Not changed.
- `TMB::tmbroot()` in `confint(method = "uniroot")` (R/confint.R) and
  in the profile intervals of `conditional_effects()`
  (R/conditional-effects.R): calls `sdreport(obj)` without a Hessian,
  so optimHess(), ABS. It uses the SD only to bracket the root at
  +-7 SD; uniroot's root does not depend on the bracket while the root
  lies inside it, and a 6 percent error in the SD does not move a
  bound at about 2 SD out of a 7 SD bracket. Not changed.
- `autoscale_prefit_converged()`, `autoscale_fit_sound()`
  (R/autoscale.R): optimHess() in `par / u`, ABS where there are no
  units. Verdicts (positive definite, condition ratio above 1e-4 or
  1e-8), not covariances. Not changed; a verdict on a badly scaled
  column could move if they read the AD Hessian, and nothing measured
  says it is wrong. Filed in section 9.

Sites that are not a Hessian of the objective, or are scaled already:

- `se_obj_curv()`, `se_line_probe()`, `se_curvature_real()`,
  `se_ridge_curvature()`, `se_span_remeasure()` (the SE check's second
  differences): along unit-diagonal (curvature-scaled) directions, with
  adaptive lengths. Scaled.
- `se_at_edge()`, `se_block_ridge_edge()` (`se_edge_step = 2` on the
  log scale), `se_sd_gain_up0()` (a grid of sds): probes of the
  likelihood, not derivatives.
- `nl_flat_message()`, `diagnose_flat()`: relative steps
  (`1e-4 * max(|p|, 1)`).
- Delta-method Jacobians in R/confint.R, R/autocor.R, R/methods-fit.R,
  R/predict.R, R/families-simplex.R, R/brms-shapes.R: relative steps;
  R/me.R's `me_hyper_table()` an absolute 1e-6 on a transformation.
  Not Hessians.
- `TMB::tmbprofile()` (confint profile): its own adaptive step from
  h = 1e-4. Not a covariance.
- frmtmb.sample: `mode_inits()` jitters the mode by N(0, 0.25) on the
  unconstrained scale; no Hessian, no mass matrix from one;
  `check_laplace()` reads core's covariance. frmtmb.latent's
  `hmm_starts()` jitters by SEs from core's `vcov()`. frmtmb.spline
  reads core's covariance. coupling, eam, learn and ode form none. No
  extension calls sdreport(), optimHess() or numDeriv.

## 3. Fix without random effects: the AD Hessian

`fit_outer_hessian()` now returns the AD Hessian, symmetrized, when the
objective has no random effects (`outer_hessian_ad()`, which reads
`obj_he_cached()`, the Hessian the SE check already took), with
`exact = TRUE`; `Hq = H * outer(u, u)` under `par_units`, an exact
diagonal scaling. It falls back to `fd_hessian()` only when RTMB
cannot build a finite Hessian. `autoscale_sdreport()` hands sdreport()
that matrix (in natural units under `par_units`, the one the check
kept), so `vcov()` is `solve(H)` of the exact Hessian, or the
well-conditioned inverse of `Hq` mapped back.
`cov_from_hessian()` keeps a Hessian taken at a correlation's clamp
(`se_edge_cor()`) as exact rather than swapping it for the one at the
fit's own point, and `sdr_rescue_hessian()` passes `exact` through.

### Exact references

GLMs (`dev/adhess-glm.R`; logs `glm-base.txt`, `glm-lane.txt`): ten
designs, poisson and bernoulli, n = 500, seeds `100 k + family`. The
reference is the inverse of `X' W X` formed and inverted in 256 bits
(Rmpfr) at each fit's own estimate; glm()'s vcov() against the same at
glm()'s estimate (epsilon 1e-12). SE err is max |SE / SE_exact - 1|,
cor err the largest covariance error over sqrt(V_ii V_jj).

    design      family     base SE / cor        lane SE / cor        glm SE / cor
    d01_unit    poisson    2.19e-07 / 4.39e-07  8.88e-16 / 1.24e-15  3.94e-07 / 7.87e-07
    d01_unit    bernoulli  2.56e-08 / 5.12e-08  4.44e-16 / 9.43e-16  2.73e-07 / 5.45e-07
    d02_x500    poisson    6.35e-02 / 1.23e-01  2.45e-14 / 4.92e-14  2.73e-07 / 5.46e-07
    d02_x500    bernoulli  1.06e-02 / 2.13e-02  7.33e-15 / 1.44e-14  4.66e-07 / 9.31e-07
    d03_sd100   poisson    3.68e-03 / 7.35e-03  5.55e-16 / 7.73e-16  2.83e-13 / 5.66e-13
    d03_sd100   bernoulli  3.59e-04 / 7.19e-04  4.44e-16 / 7.91e-16  3.63e-07 / 7.25e-07
    d04_sd1e3   poisson    2.70e-07 / 5.39e-07  6.66e-16 / 8.06e-16  3.91e-13 / 7.81e-13
    d04_sd1e3   bernoulli  8.07e-08 / 1.61e-07  4.44e-16 / 7.69e-16  7.04e-08 / 1.41e-07
    d05_mixed2  poisson    3.34e-07 / 6.68e-07  8.75e-14 / 1.75e-13  9.93e-13 / 1.99e-12
    d05_mixed2  bernoulli  2.03e-07 / 4.81e-07  7.34e-14 / 1.47e-13  2.91e-07 / 5.82e-07
    d06_shift   poisson    2.52e-08 / 5.04e-08  2.56e-08 / 5.13e-08  2.54e-12 / 5.09e-12
    d06_shift   bernoulli  2.30e-08 / 4.60e-08  2.22e-08 / 4.43e-08  4.65e-07 / 9.30e-07
    d07_three   poisson    2.07e-01 / 3.71e-01  9.99e-16 / 1.96e-15  5.88e-13 / 1.18e-12
    d07_three   bernoulli  6.88e-02 / 1.42e-01  5.55e-16 / 1.77e-15  3.26e-07 / 6.53e-07
    d08_factor  poisson    1.07e-01 / 2.03e-01  4.44e-16 / 3.37e-15  6.46e-07 / 1.29e-06
    d08_factor  bernoulli  1.13e-01 / 2.38e-01  4.22e-15 / 8.52e-15  5.34e-08 / 1.07e-07
    d09_poly3   poisson    8.68e-08 / 1.74e-07  3.73e-12 / 7.47e-12  2.71e-09 / 5.42e-09
    d09_poly3   bernoulli  8.68e-02 / 1.66e-01  7.92e-13 / 1.58e-12  5.69e-07 / 1.14e-06
    d10_skew    poisson    2.33e-15 / 4.59e-15  2.11e-15 / 4.22e-15  7.94e-08 / 1.59e-07
    d10_skew    bernoulli  2.92e-07 / 5.84e-07  1.78e-15 / 3.50e-15  1.62e-07 / 3.24e-07

No row warned on either arm. The coefficients agree with glm()'s to at
most 1.1e-5 SE. d06 (x on [1000, 1001], the design's condition number
about 1e7) is at the rounding of a double Hessian on both arms; glm()'s
QR keeps more there. The base's d06 and d10 rows are good because
fd_hessian()'s retake ran there: optimHess() at 1e-3 had no Cholesky
factor and the covariate's column rose by 53 to 1.2e30 over the step,
so it was retaken at a rise of 1e-4 (`dev/adhess-routes.R`,
`routes-base.txt`). A design whose 1e-3 step still factors keeps the
wrong matrix, which is every other bad row above.

Nonlinear (`dev/adhess-nls.R`; `nls-base.txt`, `nls-lane.txt`):
`y ~ a * exp(b * x)`, n = 200, x on x0 + [0, 1], seeds 1 to 3.
Reference: the observed information of (a, b, log sigma) in closed
form, inverted in 256 bits at frmtmb's estimate. nls() against its own
formula `(RSS / (n - 2)) (J'J)^-1` in 256 bits at its estimate.

    x0   seed  base SE / cor          lane SE / cor          nls SE / cor
    20   1     8.66e-03 / 1.33e-01    4.15e-13 / 4.87e-07    2.22e-07 / 4.44e-07
    20   2     6.11e-03 / 1.11e-01    3.45e-13 / 4.00e-07    2.93e-07 / 5.85e-07
    20   3     1.14e-12 / 5.58e-08    1.14e-12 / 5.58e-08    2.59e-07 / 5.18e-07
    0    1     6.34e-07 / 5.36e-06    5.55e-16 / 1.08e-09    5.63e-09 / 1.13e-08
    0    2     6.32e-07 / 4.97e-06    5.55e-16 / 9.53e-11    6.73e-09 / 1.35e-08
    0    3     5.81e-07 / 4.55e-06    1.55e-15 / 1.54e-08    3.82e-09 / 7.63e-09
    5    1     5.08e-05 / 2.67e-03    2.71e-14 / 7.59e-08    7.38e-08 / 1.48e-07
    5    2     5.17e-05 / 2.25e-03    7.84e-14 / 4.65e-08    7.74e-08 / 1.55e-07
    5    3     5.07e-05 / 2.38e-03    1.09e-13 / 2.00e-07    7.20e-08 / 1.44e-07

The lane's correlation errors of up to 4.9e-7 at x0 = 20 are the
rounding of a double Hessian whose scaled condition number is 2.7e4
(the point is the same: `opt$par` and `last.par.best` are identical,
`dev/adhess-nlspt.R`). frmtmb's SEs over nls()'s are 0.9929 to 1.0058
on the base and 0.9929 to 0.9974 on the lane: the definitions differ
(nls divides by n - 2 and drops the residual curvature term), so this
ratio is context only. x20 seed 3 is the same on both arms because the
fit-time SE check reached tier 3 there and kept its analysis, which is
read off the AD Hessian (`routes-base.txt`: check tier 3, kept
analysis TRUE).

## 4. Fix with random effects: steps scaled to the design

`fd_hessian()` takes a step per coordinate (`eps`, default 1e-3, the
old arithmetic bit for bit). `fit_outer_hessian()` passes
`outer_fd_steps(frame, pn, u)`: for a coefficient of a dense design
column (`beta`, `betad`, not a nonlinear body or a `mo()` term), 1e-3
over the column's root mean square in the coordinates `par / u`, where
that is above one; 1e-3 otherwise. No step moves a linear predictor by
more than about 1e-3 per row. The rule reads the design only: no
measured number decides it, so no rounding can move it.

Candidates measured (`dev/adhess-restep.R`; `restep-2.txt`,
`restep-3.txt`, `restep-4.txt`) against a Richardson Jacobian of the AD
gradient in standard-error units (numDeriv, six halvings), at two base
steps 0.1 and 0.01 whose agreement is the "refself" column (at most
3.3e-10 in SE): `abs` 0.71.0; `des` the rule above; `curv` every column
whose rise over 1e-3 exceeds 1e-2 retaken at a rise of 1e-4
(fd_hessian()'s retake, unconditional); `desR` `des` with Richardson
(twice the gradients); `desC` `des` plus `curv` on the columns outside
the designs. SE error:

    model             scale     abs       des       curv      desR      desC
    pois_xbig_re      512       3.72e-03  2.64e-07  2.63e-07  1.81e-11  2.64e-07
    pois_re           1.07      2.59e-07  2.61e-07  2.59e-07  7.67e-12  -
    bern_xbig_re      101       3.35e-04  4.74e-07  4.81e-07  1.33e-11  -
    negbin_re         50.7      1.96e-05  1.80e-07  1.80e-07  1.09e-11  -
    gamma_re          206       7.22e-05  8.18e-08  8.18e-08  1.31e-10  -
    pois_slope_re     1         2.65e-07  2.65e-07  2.65e-07  1.32e-11  -
    pois_n1e4         1.01      6.90e-08  6.90e-08  6.90e-08  5.19e-11  -
    pois_xbig_n1e4    510       2.98e-04  8.16e-08  8.16e-08  1.73e-10  -
    sleep_lmm         5.34      3.41e-07  3.41e-07  3.41e-07  3.53e-12  3.41e-07
    gauss_sigma_re    101       9.06e-03  9.09e-07  4.75e-07  9.07e-11  9.09e-07
    nl_re             1         2.05e-04  2.05e-04  5.35e-07  4.54e-09  5.35e-07
    pois_xsmall_re    0.01      2.25e-07  2.25e-07  2.25e-07  4.29e-12  2.25e-07
    bern_x10_re       10        2.83e-07  6.75e-08  2.83e-07  4.61e-12  -

On 0.71.0, gamma_re, pois_xbig_re and pois_xbig_n1e4 lost the group sd
as "flat" with a warning; with `des` none does (`restep-4.txt`, "lost
none", the fit's own vcov equal to `des`). lme4's glmer() SEs differ
from the same reference by 9e-5 to 5e-3 on the five GLMMs it fits
(`lme4 fixed SE / ref - 1` in the logs), so glmer() is not a reference
at this precision.

Chosen: `des`. It removes the scale defect at no cost, is deterministic,
and leaves a well-scaled design's Hessian bit for bit (tested). `desR`
is 2e4 times more accurate everywhere but costs twice the gradients on
every fit with random effects and changes the noise estimate `E` that
the SE check's tiers compare with (not measured against the
false-alarm battery). `curv` and `desC` fix the nonlinear body
(`nl_re`), but a retake on a measured rise is what moved the SE
check's verdict on test-se-tier.R's curved ridge in lane optstart
(`dev/optstart-findings.md`, "The step rule"), so it was not adopted.
What `des` leaves: a nonlinear body with random effects (nl_re,
2.05e-4) keeps 0.71.0's error (section 9).

## 5. How far ordinary fits move

`dev/adhess-move.R` and `dev/adhess-movecmp.R` (`movecmp.txt`): 17
model specs times seeds 1 to 3, unit-scale covariates (lm, gamma,
poisson, bernoulli, student, cumulative, Beta, sigma ~ z, mvbf with
rescor; and with random effects ri, rs, s(z) + (1 | g), sigma ~ z,
bernoulli, poisson, negbinomial, poly(x, 3)). Estimates, warnings and
lost sets are identical on all 51 fits.

    without random effects (27 fits): SE moved by at most 5.23e-07
      (median 3.33e-07), correlation by at most 1.05e-06; all 27 toward
      the exact inverse (an identity for the lane, whose covariance IS
      that inverse; sections 3's references are the independent ones)
    with random effects (24 fits): SE moved by at most 9.71e-09
      (median 0), correlation by at most 1.94e-08; vcov() bitwise
      identical on 17 of 24 (every design column at unit scale or
      below). The 7 that moved are seed 2's, whose x has a root mean
      square just above 1; the move is 1e-8 against an error of 3e-7
      either way (lane closer on 1 of 24)

## 6. Finding 8: ML on a shifted covariate (not fixed, filed)

`dev/adhess-shift.R` (`shift-base.txt`, `shift-lane.txt`): sleepstudy,
120 rows (seed 1), `Reaction ~ t + (1 | Subject)` by ML with
`t = Days + s`, against lme4 by ML on Days. Base and lane identical:

    s       code  logLik - lme4   slope - lme4   lost      warning
    0       0     -1.6e-09        -2.6e-05       none      -
    1e3     0     -9.9e-11        -4.5e-06       none      -
    1e4     0     +8.0e-13        +2.0e-07       none      -
    1e5     0     -2.452e-02      +1.96e-01      none      gradient
    1e6     0     -4.398e+01      -1.04e+01      flat 2    flat

The covariance treatment does not cover it: this is where the
optimizer stops (`dev/adhess-shiftdbg.R`, `shiftdbg-1e5.txt`,
`shiftdbg-1e6.txt`, `shiftdbg-1e6b.txt`).

- s = 1e5: nlminb stops at "X-convergence (3)", code 0, with a gradient
  of 1.3e-3 in the slope. `predictor_design_kappa()` is 6.9e4, below
  the 1e5 at which `newton_polish()` runs on a fit with random effects,
  so no Newton step is tried. One Newton step in the design's
  coordinates (`reparam_basis()`, `reparam_hessian_refined()`) from the
  stop gains 0.02454, which is the shortfall, and lands on lme4's
  objective 591.7714.
- s = 1e6: kappa 6.9e5, so `newton_polish()` runs; the natural Hessian
  calls for the design's coordinates (`hess_needs_reparam()` TRUE), but
  there the Hessian is indefinite at the stop (eigenvalues 204.6, 23.1,
  0.0113, -0.00928), the stop is near the start (slope 1.4e-6,
  intercept 293.5, the response mean), and no damped step passes the
  gate. With the new step rule in `newton_polish()`'s natural Hessian
  the result was the same (tried, then reverted).

Fix direction: optimize ML's `beta` in the design's coordinates (an
outer `gamma` with `beta = T gamma`, `T` from `design_whitener()`), as
lane illse2's REML does for `par.random`, and map back exactly; or, a
smaller change, start the steps of newton_polish() at a design kappa of
about 1e4 rather than 1e5. Both are in `R/fit.R` (newton_polish and the
fit path), where lane illse2 is working; neither is a covariance or
step-rule change.

## 7. Cost

Profiled first (`dev/adhess-cost.R`: per model three timed repeats of
`frm()` and the first `vcov()`, then one repeat under Rprof at 10 ms;
seeds 1 to 6 by model; logs `cost/cost-<model>-<arm>.txt`, the two arms
run at the same time so that they share the load). On the base the
outer Hessian (`fd_hessian()` under `fit_outer_hessian()`) is 43.9
percent of the profiled time of the 50-coefficient bernoulli GLM at
n = 1e5, and 17.6 and 18.1 percent of the 10-coefficient poisson and
gaussian ones. On the lane it is gone from the profile: the AD Hessian
it reads (`obj_he_cached()`) is the one the fit already took.

Whole fit and first `vcov()`, minimum of 3 (seconds):

    model                 np  n     fit base  fit lane  vcov base  vcov lane
    glm_bern_p50_n1e5     50  1e5   22.18     22.89     20.41      2.20
    glm_pois_p10_n1e5     10  1e5    3.67      3.79      1.63      0.78
    lm_p10_n1e5           11  1e5    3.22      3.26      1.41      0.69
    glm_pois_p3_n1e3       3  1e3    0.01      0.02      0.00      0.00
    glmm_pois_n1e3         3  1e3    0.10      0.09      0.00      0.02
    glmm_pois_n1e5         4  1e5   16.15     15.50      1.28      1.18

The fit itself does not change (the differences are within the spread
of the repeats, printed in the logs). The first `vcov()` of a fit
without random effects is 2.0 to 9.3 times faster, because it no longer
pays 2p gradient evaluations for a Hessian the fit already holds. With
random effects the steps change and their number does not, so the cost
is the same. The small models are below the 10 ms tick and say nothing.

The changed component alone, arms interleaved in one process with a
control (`dev/adhess-cost2.R`: blocks of at least 1.2 s, minimum of 5
interleaved rounds; `cost/cost2-<model>.txt`):

    model              RE    cached at fit  old s     new s     new/old  ctrl/old
    glm_bern_p50_n1e5  FALSE TRUE           19.3200   13.1000   0.678    0.999
    glm_pois_p10_n1e5  FALSE TRUE            0.9100    0.4300   0.473    0.978
    lm_p10_n1e5        FALSE TRUE            0.7500    0.4900   0.653    1.007
    glm_pois_p3_n1e3   FALSE TRUE            0.0008    0.0002   0.229    1.358
    glmm_pois_n1e5     TRUE  FALSE           1.3000    1.1600   0.892    0.931
    glmm_pois_n1e3     TRUE  FALSE           0.0061    0.0084   1.379    1.576

"new" without random effects is an uncached `obj$he()`, the worst case,
and still cheaper than the differences; "cached at fit" says the fit
had already cached the AD Hessian at the point the covariance is taken,
on all four, so in a real fit the new route costs no Hessian at all.
With random effects both routes are 2p gradients; the control's ratio
(0.931 and 1.576) shows these two rows are within the noise.

The cost models carry a covariate x ~ N(500, 100), so they show the
defect too (the `SE:` lines of the logs): on both GLMMs the base's
group sd SE is NaN and the lane's is 0.1374 (n = 1e3) and 0.03295
(n = 1e5); glm_pois_p3_n1e3's intercept SE is 0.12577 on the base and
0.13009 on the lane, 3.3 percent apart.

## 8. Tests, BLAS, false-alarm battery

All on the install of 2026-10-20 01:34 UTC (core and frmtmb.coupling
from the lane library, the other six extensions from rellib-r9 over
it; the `lib:` line of every per-file log says so).

Fragility scans, every test file of all eight suites, one R process
per file, every gate on (`dev/adhess-scan.sh`, which is
`dev/ciharden-scan.sh` with its logs in `adhess-log/scan/`). Generated
from the `.sum` files:

    r2-ref:  383 files | pass 27306 fail 0 err 0 skip 0 warn 0 | ran 383 of 383
    r2-ob26: 383 files | pass 27307 fail 0 err 0 skip 0 warn 0 | ran 383 of 383
    r2-ob32: 383 files | pass 27307 fail 0 err 0 skip 0 warn 0 | ran 383 of 383

r2-ref is the reference BLAS and LAPACK; r2-ob26 and r2-ob32 are
OpenBLAS 0.3.26 and 0.3.32 for BLAS and LAPACK, 4 threads
(`dev/ciharden-openblas.sh <ver> lapack`). frmtmb.ode's test-scale.R,
the file the pause killed, ran and passed on each (3144, 3177 and
2658 s). The one
per-file count that differs between configurations is test-v11.R (22
against 23), the pre-existing glmmTMB branch of section 9; between the
2026-10-14 scan and r2-ref the only differences are the two files
rerun after their edits (test-se-tier.R 90 to 91, test-coherence.R 67
to 66 with the changed coupling, whose eta = 16 step now refuses).

None of this lane's tests depends on a rounding-decided branch: the
bounds of test-se-exact-hessian.R hold with margins of 270 to 6e6
(section "Resumed, second part"), and its `identical()` assertions
compare a matrix with the same computation of itself.

The SE false-alarm battery on all three BLAS: section "What was done
before the pause" (820 fits per configuration, 0 rows differing in
code, lost set, SE warnings or other warnings). That ran on the
2026-10-14 install, whose core code is the code installed now (no
file under `R/` changed after 16:13 that day); the exact-reference
rerun of its "lin" and "poly" groups on the fresh install (section
"Resumed, second part", `sfx/cmp.txt`) again finds 0 rows differing in
warnings or lost set.

One session per package, files in alphabetical order, every extension
on the library path (`dev/adhess-onesession.R`, `one/<pkg>.txt`):

    ONE SESSION frmtmb.coupling: files 11 tests 76 pass 545 fail 0 error 0 skip 5 warn 0
    ONE SESSION frmtmb.eam: files 30 tests 325 pass 1808 fail 0 error 0 skip 3 warn 0
    ONE SESSION frmtmb.latent: files 10 tests 74 pass 360 fail 0 error 0 skip 2 warn 0
    ONE SESSION frmtmb.learn: files 15 tests 112 pass 430 fail 0 error 0 skip 13 warn 0
    ONE SESSION frmtmb.ode: files 11 tests 139 pass 549 fail 0 error 0 skip 1 warn 0
    ONE SESSION frmtmb.sample: files 48 tests 391 pass 2697 fail 0 error 0 skip 4 warn 0
    ONE SESSION frmtmb.spline: files 16 tests 115 pass 594 fail 0 error 0 skip 1 warn 0
    ONE SESSION frmtmb: files 227 tests 2524 pass 17462 fail 9 error 0 skip 185 warn 0

The skips are the gated tiers (only NOT_CRAN is set, as in R CMD
check). Core's 9 failures are all in test-compat-register.R (5) and
test-compat.R (4): the known one-session failure of 47cf659d with
extensions loaded (section 9).

R CMD check --as-cran (`dev/adhess-check.sh`, in
`dev/adhess-check/<pkg>/`):

    frmtmb:          Status: 1 NOTE (Skipping checking math rendering: package 'V8' unavailable)
    frmtmb.coupling: Status: OK

## 9. Not done, and defects found but not fixed

- A nonlinear body with random effects keeps 0.71.0's step for its
  parameters: nl_re (`a ~ 1 + (1 | g)`, x on 5 + [0, 1]) has SEs
  2.05e-4 off. Measured remedies: the curvature retake on those
  columns (5.4e-7) or Richardson (4.5e-9); see section 4 for why
  neither was adopted. A design-like scale for a nonlinear parameter,
  the root mean square of the body's Jacobian column at the estimate,
  would fit the rule of section 4 without a measured branch; not built.
- Finding 8 (section 6).
- `autoscale_prefit_converged()` and `autoscale_fit_sound()` decide with
  optimHess() (ABS where there are no units); on a fit without random
  effects they could read the AD Hessian. Not changed: they are
  verdicts, and changing them changes which start autoscale keeps.
- `TMB::tmbroot()`'s internal sdreport() (section 2): bracket only.
- Pre-existing, not this lane's: test-v11.R's toep test asserts against
  glmmTMB only "when it converges", and whether glmmTMB converges
  depends on the BLAS (22 expectations with the reference BLAS, 23 with
  OpenBLAS; already recorded in `dev/ciharden-findings.md`). The branch
  is glmmTMB's alone; no frmtmb covariance is involved.
- Pre-existing: poly6_lmer of the false-alarm battery (raw x^1 to x^6
  on [1, 2] with `(1 | g)`) has SEs 0.707 from a Richardson reference
  on both arms identically, and that reference is not invertible on 5
  of 20 seeds (`sfx/cmp.txt`). An ill-conditioned design with random
  effects: lane illse2's area.
- Pre-existing: core's one-session run fails 9 expectations in
  test-compat-register.R and test-compat.R, the failure lane-env.md
  records for 47cf659d with extensions loaded (fixed by lane testiso,
  not in this base). Reported, not fixed.
- Pre-existing, found in fix round 1 (`dev/adhess-kinksweep.R`,
  `dev/adhess-kinkpeek.R`, `f1/kink-*.txt`, `f1/kinkpeek-*.txt`; the
  same on both arms): no asym_laplace() or zero_inflated_asym_laplace()
  fit has a valid Hessian standard error (the log density is piecewise
  linear in the location; every fit fails the gradient check). It
  should refuse by name or use a sandwich estimate. A poisson fit to an
  all-zero response gives SEs of 1.5e136 with only the convergence
  warning (no separation verdict for poisson). `y | cens(cen) ~ x` with
  every row right-censored gives SEs near 4e9 on seed 1 with no
  warning at all (seeds 2 and 3 lose two SEs as "flat").

## Fix round 1 (review `dev/reviews/2026-10-19-adhess.md`)

Built 2026-10-20 (`install-9.txt`). Logs in `dev/adhess-log/f1/` and
`dev/adhess-log/scan/f1-*`.

### Finding 1: no variance of 0, negative or NaN is shown as itself

`sdr_nan_nonpos()` (R/se-check.R), called by `sdr_rescue_hessian()` on
both paths that keep sdreport()'s own inverse (the `cov_usable()`
return, where it changes nothing because `cov_usable()` already
requires a positive diagonal, and the `short` return): a parameter
whose variance is 0, negative or not finite gets a NaN row and column
in `cov.fixed` (and in `cov_fixed_prop` when there is one), and joins
`se_lost` with the new reason "nocurv", which `se_lost_clauses()`
renders for summary(), diagnose() and the non-finite warning. On these
fits the fit-time warnings are the gradient check's ("the curvature
there is unusable") and, on some seeds, the optimizer's; vcov() and
confint() add none (the short path has already said it).

The test ("a kinked log density shows no standard error of 0",
test-se-exact-hessian.R; asym_laplace with quantile free and fixed,
seed 2, n = 400) failed on the build before the fix
(`f1/t-exact-prefix.txt`: pass 32 fail 3, the three kink assertions,
plus an escaped "NaNs produced" from confint()) and passes after it.

The sweep (`dev/adhess-kinksweep.R`, seeds 1 to 3, n = 400): per case
the variances that are 0, negative or NaN, and confint() rows of zero
width.

    case                lane after the fix                 base (0.71.0)
    asym_laplace free   NaN 4 of 4, nocurv, fit warns      finite (noise)
    asym_laplace q 0.5  NaN 3 of 3, nocurv, fit warns      finite (noise)
    asym_laplace q 0.9  NaN 3 of 3, nocurv, fit warns      s3: 2 negative
    zi_asym_laplace     NaN 4 of 5, nocurv, fit warns      s2: 2 negative
    huber, k 1.345      finite, no warning                 same
    huber, k 0.01       finite, no warning                 same
    poisson, all y 0    finite (1.5e136), conv warning     same
    bernoulli, all y 1  NaN 2 of 2, separation             same
    cens, all right     s1 finite (4e9), no warning;       same
                        s2, s3 NaN 2 of 3, flat

On the lane no variance is 0 or negative and no confint() row has zero
width on any of the 27 fits. huber has a positive curvature wherever
some rows lie inside its quadratic core, so it is unaffected. The
poisson and censored cases are pre-existing and the same on both arms
(section 9).

### Finding 2: quadrature and importance fits

They have no inner random effects in TMB's sense, so they take the AD
Hessian. New test "an importance fit's covariance does not depend on
x's units" (poisson `y ~ x + (1 | g)`, x ~ N(500, 100), 30 groups,
seed 31, importance = 200): the covariance mapped back from the x / 100
fit is the reference, the bound is 1e-2 times the error of 0.71.0's
differences on the same fit, and the quadrature fit of the same model
must report `exact`. On rellib-r9 (`f1/t-exact-base.txt`) it fails: SE
error 0.848 and correlation error 2.41 against bounds of 0.008 and
0.02, and `exact` FALSE on both fits. NEWS and the vignette say so.

### Finding 3: the step rule's wording

NEWS, the vignette and `outer_fd_steps()`'s roxygen now say "about 1e-3
in root mean square over the rows", with the heavy-tail case. The rule
is unchanged. The reviewer's comparison (`dev/adhess-rev-step.R`,
poisson GLMM, n = 400, 20 groups, seeds 1 to 3, largest SE error
against Richardson):

    column                    base       lane       1e-3/max|x|
    x ~ N(0, 1e-2)            1.3e-07    1.3e-07    1.3e-07
    1000 + N(0, 1)            NaN, lost  1.3e-07    1.3e-07
    5 + N(0, 0.1)             4.8e-07    1.3e-07    1.3e-07
    lognormal, RMS 0.9        5.5e-06    5.5e-06    2.5e-07
    lognormal, RMS 3          6.1e-05    6.8e-06    2.5e-07
    lognormal, RMS 50         1.6e-02    6.8e-06    2.5e-07
    N(500, 100)               2.2e-03    1.3e-07    1.3e-07
    N(500, 100) * dummy       1.7e-02    1.5e-07    1.4e-07
    dummy + N(0, 1)           1.3e-07    1.3e-07    1.3e-07
    bernoulli, N(500, 100)    2.8e-02    1.5e-06    1.5e-06

`1e-3 / max|x|` is 30 times more accurate on heavy tails and equal
elsewhere, but it changes the steps of a standardized covariate (max
|x| about 3), so a design at unit scale would no longer keep 0.71.0's
Hessian bit for bit. Not adopted; the error it would remove is 7e-6.
The vignette also says now that nonlinear, `mo()` and variance
parameters keep 1e-3, and that an ill-conditioned design takes its
covariance from the design's QR coordinates.

### Finding 4: the coherence refusal

One bound on the logit scale, `abs(eta) + z * se > qlogis(1 - eps)`
(36.04), at both ends (`cp_logit_edge`). The old test of the rounded
ends refused only near 1: a lower end at a logit of -36.24 has
plogis() 1.8e-16 > 0 and was returned (`f1/edge.R`). New test "the
coherence refusal is the same at both ends" mocks `cp_link_se()`: both
ends a tenth of an SE past the bound are refused, and the same
intervals a tenth of an SE inside are returned inside (0, 1). With
rellib-r9's coupling over the lane's core (`f1/t-coh-basecoupling.txt`)
it fails (the lower end comes back as a data frame), as does the eta =
16 loop step. The coupling NEWS bullet is under "Breaking changes".

### Finding 5: memory

NEWS and the vignette's cost paragraph give the reviewer's numbers
(`dev/adhess-rev-cost.R`): under `check_se = "ignore"` the first
vcov() builds the AD Hessian, with a peak working set 350 MB (p = 200)
and 870 MB (p = 500) above 0.71.0's on a bernoulli GLM with 1e4 rows.
`frmtmb_control()`'s `check_se` documentation is in R/fit.R, lane
illse2's file, so it is not changed here; the consolidation may add
the same sentence there.

### Finding 6: stale comments

`fd_hessian()`'s roxygen now says the bit-for-bit claim holds at its
default `eps` only; test-se-check.R's comment says the same.

### Tests of fix round 1

One file per process: the SE files of core, every file that fits
asym_laplace (core and frmtmb.sample) and frmtmb.coupling's
test-coherence.R (`dev/adhess-log/f1/jobs.txt`; `scan/f1-ref.sum`,
`scan/f1-ob26.sum`). Reference BLAS and OpenBLAS 0.3.26 with LAPACK, 4
threads, give the same counts:

    frmtmb RESULT test-se-check.R pass=129 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-se-conv.R pass=21 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-se-exact-hessian.R pass=35 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-se-flags.R pass=68 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-se-lost-re-predict.R pass=20 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-se-reparam.R pass=57 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-se-tier.R pass=91 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-brms-suite-brmsterms.R pass=5 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-families.R pass=252 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-numerical-robustness.R pass=723 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-optstart.R pass=69 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-simulate-density.R pass=591 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-v14.R pass=81 fail=0 err=0 skip=0 warn=0
    frmtmb.sample RESULT test-brms-suite-brmsterms.R pass=4 fail=0 err=0 skip=0 warn=0
    frmtmb.coupling RESULT test-coherence.R pass=72 fail=0 err=0 skip=0 warn=0
    ran 15 of 15

Against the pre-fix scan r2-ref the only counts that moved are the new
tests' (test-se-exact-hessian.R 25 to 35, test-coherence.R 66 to 72).

R CMD check --as-cran (`dev/adhess-check.sh`):

    frmtmb:          Status: 1 NOTE (Skipping checking math rendering: package 'V8' unavailable)
    frmtmb.coupling: Status: OK

The first frmtmb check of this round was stopped: its tarball held a
stray 90-byte file `0` at the worktree root, written at 00:11 by an
inline `Rscript -e` whose printed text "(> 0, ...)" the shell took as a
redirect into a file named `0` (a one-off command, not a script; the
script `f1/edge.R` that replaced it has no redirect). That check got
"Non-standard file/directory found at top level: '0'". The file is
deleted, and the rerun above (`check3-frmtmb.out`, tests 534 s) has the
V8 NOTE only.

## Merge with illse2

The reviewer's trial merge (`dev/reviews/2026-10-19-adhess.md`,
"Composition with lane illse2"), built outside both worktrees and
deleted; nothing is merged here.

- Three-way `git merge-file` of the shared files against 47cf659d,
  line endings normalized to LF first: R/se-check.R clean;
  tests/testthat/test-se-tier.R clean; R/autoscale.R one conflict on
  the last line of `autoscale_sdreport()`; NEWS.md one conflict (both
  lanes add "Breaking changes").
- The autoscale conflict keeps both:
  `sdr_rescue(fit, sdr_unwhiten(fit, sdr_name_joint(sdr, obj)), Hn)`.
- Precedence: no change needed. Where a design without random effects
  is ill-conditioned (unit-diagonal condition at or above
  `se_reparam_kappa`, 1e8), illse2's QR factor (tier 4) is stored in
  `fit$cache$se_analysis` and `sdr_rescue_hessian()` puts it in the
  sdreport, so illse2's route wins there; elsewhere the AD inverse. The
  plain AD inverse is not accurate above about 1e8 (SE error 1.9e-4 to
  1.1e-3 at a condition of 1.3e13 to 1.7e13; `solve()` fails at 3e15),
  so the factor is needed there, and the route supplies it. The
  reviewer measured the merged build identical to illse2 on every row
  of illse2's accuracy battery (gaussian, poisson and bernoulli without
  random effects; gaussian with `(1 | g)`), and the merged build passed
  the six SE test files of both lanes.

Recipe:

1. Take illse2's versions of the files only it changes (R/fit.R,
   R/se-reparam.R, R/brms-shapes.R, R/confint.R, R/influence.R,
   R/interop.R, R/methods-fit.R, R/multiple.R, R/predict.R,
   man/vcov.frmtmb_fit.Rd, tests/testthat/test-se-reparam.R) and
   adhess's for the files only it changes (vignettes/diagnostics.Rmd,
   tests/testthat/test-se-check.R, tests/testthat/test-se-exact-hessian.R,
   frmtmb.coupling).
2. Add illse2's two new test files, untracked in its worktree:
   tests/testthat/test-se-reparam-inv.R and test-se-reparam-reml.R.
3. R/se-check.R and tests/testthat/test-se-tier.R: three-way merge;
   both were clean. Fix round 1 since added `sdr_nan_nonpos()`, edited
   the two early returns of `sdr_rescue_hessian()` and
   `se_lost_clauses()`, and the roxygen of `outer_fd_steps()` and
   `fd_hessian()`; recheck the R/se-check.R merge (illse2 edits
   `se_bound_held()`, the tier 3 tail of `cov_from_hessian()`,
   `se_sd_gain_up0()` and `se_sd_scale()`).
4. R/autoscale.R: resolve the one conflict as above.
5. NEWS.md: one "# frmtmb (development version)" heading and one
   "## Breaking changes" list with both lanes' bullets.
6. Rerun the six SE test files of both lanes, frmtmb.coupling's
   test-coherence.R, and the OpenBLAS 0.3.26 and 0.3.32 scans of core
   (neither lane has run on the other's code under OpenBLAS).
7. If wanted, add the memory sentence of finding 5 to
   `frmtmb_control()`'s `check_se` documentation in R/fit.R.
