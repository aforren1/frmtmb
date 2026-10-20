# Lane illse2: findings

Lane `illse2` of the round opened 2026-10-14. Worktree
`frmtmb-wt-illse2`, branch `wt-illse2`, base 47cf659d (frmtmb 0.71.0).
Lane library `C:/Users/adf44/source/r/wt-illse2-lib` over `rellib-r9`;
"base" is `rellib-r9` alone. Logs are local, in `dev/illse2-log/`;
every number below is copied from them, with the script that wrote it.
Shared setup and 256-bit helpers: `dev/illse2-common.R`
(`ILLSE2_ARM=lane|base`). Install script: `dev/illse2-install.sh`.

Final build: install 21 (`dev/illse2-log/install-21.txt`, 2026-10-19,
on the restored user library), core only; no extension changed. Bump:
minor (new behavior on REML, nonlinear and bounded fits, marked
BREAKING in NEWS); the extensions need no new floor. Every accuracy
number below is from install 21's logs in `dev/illse2-log/f4/`,
summarized by `dev/illse2-f4-sum.R` into `f4/summary.txt` (section 12,
"Fix round 1, resumed"), unless it names another build.

Summary:

1. `cooks.distance()` and D1 and D2 (Wald) read the factor: within
   2.3e-6 of exact at degree 8, where 0.71.0 stopped with
   "computationally singular" at degrees 7 and 8 (section 1).
2. Degree 7 seed 4 was really short (1.2e-3, by an exact profile); an
   sd escape takes it, and two other short fits, to the maximum, and
   the route then takes the interior sd (section 2).
3. REML is taped in the design's coordinates, with the Jacobian
   `-log |det T|` under REML: the log-likelihood is exact to 6.6e-10 and
   the covariance to 6.1e-9; nonlinear bodies take their Jacobian's
   coordinates; a guard refuses an asymmetric Hessian (section 3).
4. Newton steps hold the coordinates a bound stops; fits up to 400
   outer parameters take the route (section 4).
5. Half the gradients per Hessian on identity-link gaussian designs:
   0.87 and 0.89 of 0.71.0's time at n = 1e4, 0.88 and 0.88 at n = 1e5,
   without and with `(1 | g)`; other fits keep 0.71.0's differences
   (section 5).
6. Assessed, not implemented (section 6).

## 1. Item 1: routes that invert the covariance

### What changed

- `R/se-reparam.R`: `root_inv_quad(K, D)` forms `d' V^-1 d` for each
  column of `D` with `V = K K'` held as its factor: rows of `K` scaled
  to unit length, a pivoted LAPACK QR of `K'`, `|R^-T P' d|^2`. NULL
  when the rows are dependent to working precision (smallest diagonal
  of `R` at or below `p` epsilon of the largest); the caller then gives
  NaN with a warning (`warn_cov_singular()`, class
  `frmtmb_cov_products`).
- `R/interop.R`: `coef_vcov_root(model, extra)`, the factor's rows for
  `vcov_estimated()` (`extra = FALSE`) or `interop_vcov()` (`TRUE`).
  `cov_routes_right()` now names `cooks.distance()` and the Wald rules
  of `anova()` on `frm_multiple()` fits; `warn_cov_products()` no longer
  lists them among the routes that invert the whole matrix.
- `R/influence.R`: `cooks.distance.frmtmb_influence()` reads the factor
  when the fit holds one and `V` is finite; otherwise `solve(V)` exactly
  as before.
- `R/multiple.R`: `pooled_wald_roots()` gives a factor of each
  imputation's `U_j` when at least one fit holds a factor (a fit without
  one gets the pivoted Cholesky factor of its `U_j`); D2's per-imputation
  chi-squares and D1's trace and form (`ubar = Kb Kb'`, `Kb` the factors
  side by side over `sqrt(m)`, `bvar = C C'`) read them. Without a
  factor the code is the old one.

### Accuracy

`dev/illse2-p1.R <re> <deg> 1:3` (gaussian raw polynomial on [1, 2],
n 200, 20 groups of 10; three imputations that differ by noise of sd
0.05 in y). Reference: each fit's exact covariance at its own point
(256-bit Hessian, closed form without random effects, second
differences of the closed-form marginal with `(1 | g)`; inverted in 256
bits); `cooks.distance()` against `d' (V_sub)^-1 d / p` with the
group-deletion `d` the influence table holds; D1 and D2 (Wald, every
power) against the same rules evaluated with the exact `U_j`. Largest
|ratio - 1| over seeds 1 to 3 (`dev/illse2-p1sum.R`,
`p1-summary.txt` for the base; lane rows from `f4/p1-lane-*.txt`,
install 21, `f4/summary.txt`, identical to install 20's; 0.0e+00
means below 5e-10, the printed digits):

    arm  design  deg | cooks         | D1 all        | D2 all
    base no RE  deg 5 | 1.4e-06      | 4.3e-08       | 6.1e-08
    base no RE  deg 6 | 5.8e-04      | 4.0e-06       | 7.6e-06
    base no RE  deg 7 | 3 of 3 error | 3 of 3 error  | 3 of 3 error
    base no RE  deg 8 | 3 of 3 error | 3 of 3 error  | 3 of 3 error
    base (1|g)  deg 5 | 4.5e-07      | 6.1e-08       | 7.9e-08
    base (1|g)  deg 6 | 5.6e-05      | 4.4e-06       | 3.0e-05
    base (1|g)  deg 7 | 3 of 3 error | 3 of 3 error  | 3 of 3 error
    base (1|g)  deg 8 | 3 of 3 error | 3 of 3 error  | 3 of 3 error
    lane no RE  deg 5 | 7.5e-10      | 0.0e+00      | 0.0e+00
    lane no RE  deg 6 | 7.6e-09      | 0.0e+00      | 3.0e-09
    lane no RE  deg 7 | 7.1e-08      | 2.0e-09      | 6.0e-09
    lane no RE  deg 8 | 9.6e-07      | 4.3e-08      | 2.1e-07
    lane (1|g)  deg 5 | 1.1e-09      | 0.0e+00      | 0.0e+00
    lane (1|g)  deg 6 | 1.7e-08      | 1.0e-09      | 3.0e-09
    lane (1|g)  deg 7 | 1.1e-07      | 1.0e-08      | 5.0e-09
    lane (1|g)  deg 8 | 2.3e-06      | 2.2e-07      | 5.1e-07

Under OpenBLAS 0.3.32 LAPACK (install 17, `p1-summary-ob32.txt`), the
lane at degrees 7 and 8: Cook's 1.3e-7 to 1.5e-6, D1 5.0e-9 to 3.4e-8,
D2 3.1e-8 to 1.5e-7.

The errors on the base are `solve()`'s "system is computationally
singular" (reciprocal condition numbers 1e-20 to 1e-19). At degree 6
the base was silently off by up to 5.8e-4 (Cook's distance). The
riv of D2 (a variance of square roots of nearly equal chi-squares) is
cancellation-limited on both arms (up to 2.6e-5 on the lane at degree
8, 3.5e-5 on the base for the two top powers); the statistic is not.
The test of the two top powers alone (a well-conditioned 2 x 2 block)
gives the same digits on both arms.

The double-precision reference of the new tests (`d' H d` through
`X d`; `|Xc q|^2 / s2` with centered columns) agrees with the 256-bit
one within 2.9e-9 (Cook's) and 1.4e-10 (the Wald form), seeds 1 to 3
(`dev/illse2-p1ref.R`, `p1ref.txt`).

### Tests

`tests/testthat/test-se-reparam-inv.R` (new, 13 expectations): on the
base 2 pass and 4 error (`t-inv-base.txt`): `solve.default(V)` and
`solve.default(U[, , j], Q[, j])` "computationally singular" (the
behavioral failures), and the two unit tests of the new helpers. On the
lane 13 pass.

## 2. Item 2: degree 7, seed 4 of the 80 sd-at-zero fits

### The fit was short

`dev/illse2-p2.R 7 4` (`p2-base-7-4.txt`): frmtmb 0.71.0 stops at code
0 with the group sd at 1.877e-6 and log-likelihood -61.7572196279. The
exact profile of the log-likelihood in tau = sd / sigma (256 bits, beta
by GLS and sigma^2 = Q / n) has its maximum at tau 0.0331201,
-61.7559980476: the fit is 1.222e-3 short, above `grad_tol`, as its
warning said. lme4 (ML) ends at -61.7598409365 (sd 0.01274), 3.8e-3
below the maximum and below frmtmb's stop.

### What changed

- `R/fit.R`: `sd_up_escape()`, called in `fit_assembled()` after the
  fit is assembled (not under importance sampling or quadrature). For
  each log sd below `se_sd_escape_frac` (0.01) of its term's scale
  whose upward probe (`se_sd_gain_up()`, the rest held) gains more than
  `grad_tol`, the optimizer restarts with that sd at the probe's best
  value; the result is kept when its objective is lower, then
  `newton_polish()` runs on it as on the first stop. Another optimizer,
  or an `optCtrl` that asks for an earlier stop (optctrl_stops_earlier(),
  section 12, finding 2; install 17 used any non-default setting), turns
  it off (test-se-conv.R's `iter.max = 0` stop stays short).
  `opt$sd_escape` records the gain.
- `R/se-check.R`: `se_sd_gain_up()` returns attribute `at`, the log sd
  of its best gain; `se_sd_scale()` is the scale it was computed from,
  factored out.
- `R/se-check.R`, `cov_from_hessian()`: when tier 3 removes nothing and
  nothing is held or at an edge, the design's route is asked again
  without `se_reparam_defer()` (`se_reparam_cov(defer = FALSE)`). After
  the escape the sd of seed 4 is interior (SE of log sd 10.2) and below
  `se_edge_screen`, so the route deferred, tier 3 removed nothing, and
  `solve()` of the natural Hessian was 1.37e-3 from the exact
  covariance on the scale of the correlation; through the route it is
  5.63e-6 (`f4/p2-lane-7-4.txt`, install 21, as on install 20; this
  line said 2.38e-6, from a log that was overwritten and that no
  surviving log holds).

### Result

Lane, `dev/illse2-p2.R 7 4` (install 21, `f4/p2-lane-7-4.txt`): code
0, log-likelihood -61.7559980552, 7.6e-9 below the exact maximum, sd
0.009853; nothing lost, no warning; `vcov(full = TRUE)` within 5.63e-6
of the exact covariance at the fit's point (cor), 6.19e-7 on the
diagonal; SE(log sd) 10.2 as exact.

The 80 fits (`dev/illse2-edge4.R`, reference BLAS,
`edge4-base-ref.txt` for the base, `f4/edge4-ref-<deg>.txt` for the
lane on install 21, the same counts as installs 17 and 20; "short" is
the exact profile maximum minus the fit's log-likelihood):

    arm   short > 1e-3   1e-6 to 1e-3   <= 1e-6   theta_1 verdicts
    base  3              3              74        20 SE, 50 boundary,
                                                  6 nonfinite, 3 short,
                                                  1 flat
    lane  0              3              77        22 SE, 50 boundary,
                                                  6 nonfinite, 2 short

Three fits differ, each a base fit short by more than `grad_tol` that
the lane takes to the maximum through the escape:

    fit         base                                   lane
    deg 6 s16   0.17 short, theta_1 "flat"             2.0e-8, SE 0.962
    deg 7 s4    1.2e-3 short, 8 coefficient SEs NaN    7.6e-9, SE 10.2
    deg 8 s3    1.6e-2 short, theta_1 "short"          6.4e-9, SE 2.92

Under OpenBLAS 0.3.32 (`edge4-lane-ob32.txt`, install 17) the counts
are the same. Two fits differ between the BLAS in nlminb's code and so
in the "did not report convergence" warning, with the same theta_1
verdict and SE: degree 7 seed 11 (as on 0.70.0 and 0.71.0) and degree
8 seed 3, whose restart from the escape stops at "false convergence
(8)" under the reference BLAS and "X-convergence (3)" under OpenBLAS;
one Newton step then reaches the maximum under both (gain 5.4e-4 and
5.2e-4), and the warning says the stop is the maximum (a Newton step
gains 5.7e-9; `dev/illse2-p2d8s3.R`). Degree 6 seed 11 prints SE 400
against 399.

The lane fits still short by 1e-6 to 1e-3: degree 5 seed 4 (1.1e-4)
and degree 6 seed 4 (1.0e-4), theta_1 "short", as on the base (below
`grad_tol`, so no escape); degree 6 seed 11 (1.6e-6, exact maximum at
tau 0, SE of log sd 400 on both arms).

### Tests

- `test-se-reparam.R`, "an sd stopped near zero below a small interior
  maximum moves": degree 7 seed 4 against the profile maximum in double
  (whitened GLS through QR); base 3 failures (a warning, 1.98e-5
  relative short, NaN SEs), lane passes.
- `test-se-tier.R`, "an sd the fit left near zero short of the maximum
  is no boundary": its premise (a start of log sd -20 leaves the fit
  0.825 and 0.665 short) is what the escape removes. It now asks that
  the fits reach lme4's maximum (5e-11 and 4.8e-10,
  `dev/illse2-p2tier.R`) with no boundary message and no warning; base
  6 failures, lane passes. The short-sd stop stays pinned by
  test-se-conv.R (a caller's `iter.max = 0`).

## 3. Item 3: REML and nonlinear formulas

### Nonlinear formulas

What changed (`R/se-reparam.R`): `reparam_basis(frame, pn, cells,
est)` takes, before the designs, a block for each nonlinear body from
its Jacobian at `est` (`reparam_nl_jacobian()`: central differences of
the body's predictor on its link scale with respect to the outer
coefficients of the nonlinear parameters, a step of 1e-3 of each
coefficient's size or of 1). A body linear in its parameters has its
design as Jacobian. Every caller passes the estimates:
`se_reparam_cov()` and `headroom_reparam_hessian()` `fit$estimates`,
`newton_polish()` `obj$env$parList(p)`. A body with a `ps()` term, or a
parameter whose predictor is not a plain design, gets no block.

`dev/illse2-p3nl.R` (reference: the exact Hessian at the fit's point in
256 bits; the objective agrees with the written-out likelihood up to a
constant, which moves by at most 1e-10 between two points):

    case  seed  base                                   lane
    poly  1     0.41 short, 9 "illcond"                3.9e-8 short, cor 4.7e-6
    poly  2     0.047 short, gradient warning, NaN     1.1e-8, cor 3.8e-6
    poly  3     0.15 short, 9 "illcond"                3.5e-7, cor 1.3e-6
    sin   1     cor 3.0e-3, silent                     cor 4.9e-9
    sin   2     cor 1.8e-3, silent                     cor 1.8e-8
    sin   3     cor 4.9e-4, silent                     cor 4.7e-9

`poly` is the degree-8 raw polynomial as a body with one parameter per
coefficient (kappa 5.4e17 to 5.9e17 at unit diagonal, from the exact
Hessian); "short" is against lm(). `sin` is a degree-6 raw polynomial
beside `A sin(w z)`, z on [0, 10], n 300 (kappa 1.5e13 to 1.7e13).
`nls()` stops with an error on every `sin` seed from lm()'s
coefficients and the true A, w, so it gives no reference there.

A rational body, a degree-6 polynomial over `(1 + dd x)` (`rat`), is
weakly identified along dd (kappa 6.6e11 to 9.5e13). On seeds 1 and 2
the inverse of the exact Hessian at frmtmb's stop has negative diagonal
entries (the stop is not a maximum), and both arms lose 7 coefficients
as "flat" and 3 as "concave". Seed 3 is the case that made a guard:
a code-1 stop at the maximum (exact Newton gain 5.0e-9; BFGS on the
256-bit likelihood gains 5.1e-10), with the smallest unit-diagonal
eigenvalue 7.8e-14 of the largest (`dev/illse2-p3ratdbg.R`,
`p3ratdbg-lane-3.txt`). The Jacobian's coordinates there are not
linear enough for the differences: a step of 1e-3 along a column moves
dd far, and the Hessian in those coordinates came out with its entries
48 percent off the exact `M' H M` while its eigenvalues looked fine
(condition number 1.09); the route took it and gave SEs 0.82 to 0.90
of the exact ones, silently (install 7 to 16, central differences on
those columns). With Richardson the diagonal turned negative. Its
asymmetry, `max |H - H'| / 2` over its largest entry, was 1.7e-2 then
(1.0 on seed 1), and is 0.165 under the final scheme (`fx/p3asym.txt`,
install 20), where every fit the route should take has at most 1.71e-5 (27 fits:
gaussian, bernoulli, poisson, cumulative raw polynomials of degree 7
and 8 with and without `(1 | g)`, and the mixed-scale binomial;
`dev/illse2-asym2.R`, `f4/asym2.txt`). So se_reparam_cov() refuses a
Hessian whose asymmetry is above `se_reparam_asym` (1e-3). Under the
final scheme the route also refuses seed 3 with the guard off (its
other checks; `f4/mock-noguard.txt`), so this body no longer shows the
guard's effect; the body that does is `A exp(-k z)` with A = 0.05
(section 12, finding 4). Seed 3 gets 0.71.0's covariance: 0.474 from
exact on the scale of the correlation (`f4/p3ratdbg-lane-3.txt`), the
pre-existing defect filed below (sdreport()'s finite-difference
covariance on a fit without random effects), not a new one.

Test: `test-se-reparam.R`'s "a nonlinear body keeps the ill-conditioned
verdict" becomes "a nonlinear body takes the coordinates of its
Jacobian": nothing lost, no warning, lm()'s log-likelihood to 1e-8 and
its ML covariance to 1e-4; base 5 failures.

### REML (and `profile = TRUE`)

What "accurate" means. Under REML `beta` is an inner parameter,
integrated with the group effects, and frmtmb's coefficient covariance
is the `beta` block of the joint covariance its joint precision
carries: the conditional covariance at the variance parameters,
`sigma^2 (X' V^-1 X)^-1` (what lme4's `vcov()` reports), plus what the
variance parameters' uncertainty adds, `J V_theta J'`, with `J` the
derivative of the GLS estimate in (log sigma, log sd) and `V_theta` the
inverse Hessian of the restricted log-likelihood. That is frmtmb's
definition on every REML fit; the second part was 2.1e-9 to 6.2e-3 of
the variances here (degrees 5 to 8, `f4/p3reml-lane-<deg>.txt`).
`dev/illse2-p3reml.R` computes both parts exactly in 256 bits at
frmtmb's own REML estimates (the restricted log-likelihood in closed
form with a 256-bit log determinant; `J` by 256-bit central differences
of the GLS estimate), and finds the exact REML optimum (BFGS on the
256-bit criterion).

Why the design's route did not reach it. The outer Hessian has no
`beta` row, so se_reparam_cov() never saw the design. The inner Hessian
carries `X' W X`: at degree 7 to 8 the Laplace solve and its log
determinant lose their digits, and the REML objective itself is wrong
there, so the fit stops short.

What changed (`R/se-reparam.R`, `R/fit.R`, `R/autoscale.R`, consumers):

- `beta_whitener(frame)`: for each linear predictor whose coefficients
  are `beta` and enter through an ill-conditioned dense design, the
  block `T` of design_whitener() (which now also returns `T^-1`, from
  `R`). `fit_assembled()` tapes the objective in `gamma`, `beta = T
  gamma` (`beta_whiten_nll()`; the wrapper forces its `nll`, since the
  caller rebinds the name), for REML and `profile = TRUE`, not under
  quadrature, importance sampling or an objective-only build; the
  start is mapped in and the estimates out (`est_at()`), and the fit
  carries `beta_whiten`.
- `sdr_unwhiten()` (in `autoscale_sdreport()`): `par.random`'s `beta`
  becomes `T gamma`; the joint precision `Q_gamma = R'R` gives the
  factor `joint_root` of the joint covariance in natural coordinates
  (`R^-1` with its `beta` rows times `T`), marked `lossy` as
  se_reparam_cov() marks its factor; the joint precision is carried
  back, `W^-T Q_gamma W^-1`, for readers of its entries. At most
  `joint_root_max` (4000) joint parameters.
- `joint_cov_whole()` is read by `vcov()`, `hyp_par_cov()`,
  `fit_draw_space()` and `get_joint_cov()` instead of
  `solve_joint_precision()`; `square_root_rows()` reduces the factor's
  rows to a square factor for the draws and differences. The
  `frmtmb_cov_products` warning, `get_vcov()`'s NaN and
  `coef_vcov_root()` (cooks.distance, D1) read `joint_root` too.

`dev/illse2-p3reml.R <deg> 1:3`, gaussian with `(1 | g)`, n 200; base
rows from `p3reml-summary.txt`, lane rows from
`f4/p3reml-lane-<deg>.txt` (install 21; `f4/summary.txt`); "SE dev"
the largest |SE / exact - 1| of the
coefficients, "cor" every entry of `vcov()` over sqrt(Ve_ii Ve_jj),
"short" the exact REML optimum minus the fit's):

    arm   deg  short            SE dev           cor              warnings
    base  5    5.4e-13-1.4e-10   2.0e-6-7.5e-6     4.0e-6-1.5e-5     none
    base  6    2.6e-9-4.6e-4     1.0e-3-2.8e-3     2.1e-3-5.6e-3     2 of 3 conv
    base  7    1.8e-3-1.0e-1    0.50-0.88        1.2-2.6          conv, grad,
                                                                  short
    base  8    0.65, 16; s3 error  NaN             NaN               conv, grad
    lane  5    0-2.0e-12        9.5e-12-2.2e-10  1.9e-11-4.4e-10  none
    lane  6    0-8.5e-14        2.0e-12-4.7e-10  4.1e-12-9.3e-10  none
    lane  7    1.4e-14-3.1e-13  3.1e-10-3.2e-10  6.1e-10-6.5e-10  none
    lane  8    0-9.9e-14        3.4e-10-3.1e-9   6.8e-10-6.1e-9   none

(base degree 8 seed 3: "The optimizer failed on this model ... NA/NaN
gradient evaluation".) `fitted(re_formula = NA)`'s SE on the lane: at
most 1.8e-7 from exact (degree 8). lme4's REML `vcov()` against the
exact conditional covariance at lme4's own estimates: within 1.5e-6 at
degree 5, 1.5e-4 to 8.6e-4 at degree 6, 2.2 to 3.3 percent at degree 7
(where lme4's stop is 3.6e-4 to 3.0e-2 below the REML optimum), and
"Downdated VtV is not positive definite" at degree 8. So on these
designs lme4 is not a reference; the 256-bit computation is.

The `vcov()` products warning fires on the lane fits at degree 6 to 8
(the joint factor is lossy there), as on an ML fit.

Non-gaussian REML (Laplace over `beta` and the group effects), degree 7
with `(1 | g)`, n 400, seeds 1 to 3, against the same model written with
`poly(x, 7)` (an orthogonal basis; REML's variance parameters do not
depend on the basis of the fixed effects, and neither do the fitted
values' SEs) (`dev/illse2-p3remlglm.R`, `p3remlglm-*.txt`):

    family     arm   code     variance pars - poly  fitted SE / poly
                                                     (warnings)
    bernoulli  base  1, 1, 0  0.41, 0.12, 0.053     0 to 1.55 (conv, grad, sh.)
    bernoulli  lane  0, 0, 0  6.4e-6, 4.3e-6, 1e-10 within 1.4e-6 (none)
    poisson    base  0, 1, 1  0.038, 0.040, 0.024   1.00-1.32 (grad; conv, sh.)
    poisson    lane  0, 0, 0  1.8e-6, 7.6e-7, 5e-9  within 5.3e-7 (none)

("sh." is the "short of the maximum" clause of the convergence warning.)

Tests: `tests/testthat/test-se-reparam-reml.R` (new, 19 expectations):
degree 7 seed 1 against a double-precision reference (whitened GLS by
QR, `J` by central differences, `V_theta` from the fit), SEs and fitted
SEs within 1e-4; the guard's absent case (degree 2, no whitener, no
factor); and the two log-likelihood tests of section 12, finding 1.
Base (`f4/t-base-test-se-reparam-reml.R.txt`): 8 pass, 10 fail (2
warnings, no factor, SE 0.497 off, fitted SE 0.164 off, no whitener on
the four log-likelihood fits, logLik 4.6e-8 and 2.3e-5 relative off the
restricted likelihood at degrees 5 and 6), 1 error (the missing
helper); lane 19 pass.

## 4. Item 4: bounds, and fits with more than 100 outer parameters

### A bound that holds a coefficient

`dev/illse2-p4bound.R 5:8 1:5`: gaussian raw polynomial, n 200, an
upper bound one SE below lm()'s estimate of the top coefficient, so the
constrained maximum is lm() with that coefficient as an offset at the
bound (`p4bound-base.txt`, `p4bound-lane.txt`).

Three things stopped these fits on the base, each changed:

- `newton_polish()` refused every step that crossed a bound, and every
  full step did. Now a step that leaves the box is taken again with the
  coordinates it crosses held at their bounds (the Newton step of the
  quadratic model under those equality constraints, `s = -H^-1 (g + C
  lambda)`, `C' s = r`), in natural units or the design's; the line
  search runs along it. Without a held coordinate the step and its gain
  are computed exactly as before.
- On an ill-conditioned design the gradient's sign does not say whether
  a bound holds: at the constrained maximum of degree 7 seed 3 (7.5e-8
  short) the gradient along the top coefficient was +6.6e-4, rounding,
  so the headroom counted the step past the bound (0.524) and warned
  that the fit was short. `bound_held_by_step()` also holds a parameter
  on its bound that the Newton step over the others would push out
  (conv_headroom(), grad_headroom(), se_bound_held() when given `H`).
- The route refused a removed parameter inside a design's block.
  `reparam_basis(drop = )` leaves it out of its design's block, so the
  others' covariance is taken in the coordinates of the design without
  it; se_reparam_cov() and headroom_reparam_hessian() pass the removed
  set, and the headroom no longer skips the design's coordinates when a
  bound is active.

Result, 20 fits per arm:

    arm   at the constrained max   silent and wrong         warned, short
    base  10 (SE dev up to 2.9e-4)  2 (deg 6 s1, s2: 0.12   8
                                    and 0.60 short, SEs
                                    5.8 to 39 times lm's)
    lane  17 (SE dev up to 2.1e-4)  0                       3 (deg 8 s1-3)

The lane's 17 include degree 7 seed 1, 1.7e-5 short at nlminb's
evaluation limit (code 1), SE dev 2.1e-4; the others are within 8.8e-7
of the constrained maximum with SE dev at most 5.1e-5. The three
degree-8 fits that stay short (2580 to 2700 below) start badly: the
design's condition number is above the least-squares start's limit, so
the default start puts the top coefficient at its bound with the others
at zero, nlminb stops at once at a point where the Hessian is
indefinite (unit-diagonal eigenvalues -0.73 to 9.97), and the Newton
steps do not run from a saddle (`dev/illse2-p4bdbg.R 8 1`). They warn
(large gradient, "the curvature there is unusable"), as on the base.
Not fixed: a start that solves the least squares with the bound
coefficient held would be the fix, and it is lane optstart's code.

Test: `test-se-reparam.R`, "an active bound on an ill-conditioned
design is reached" (degree 6 seed 1): base 4 failures (no bound
message, nothing lost, 2.7e-3 relative short, SEs 35 times off); lane
passes.

### More than 100 outer parameters

`se_reparam_np` is 400 (was 100); `newton_natural_np` (100) keeps
newton_polish()'s natural-unit steps where they were, so a fit above
100 parameters that does not need the design's coordinates takes no
step, as before. The `ill` flag of fit_assembled() and newton_polish()'s
own cap read `se_reparam_np`.

`dev/illse2-p4wide.R <K> 2000 1:2` (degree-7 raw polynomial beside K
standard normal covariates, gaussian, n 2000; install 21,
`f4/p4wide-<K>.txt` and `f4/p4wide-<K>-base.txt`, 2026-10-19; the arms
ran under different machine loads, so the times are not a timing
comparison):

    outer  arm   time (s)    logLik - lm()       SE / lm()
    159    base  13.6, 13.1  -2.2, -0.63         0.52-1.0, NaN (warned)
    159    lane  24.7, 27.5  -6.0e-8, -2.4e-7    within 5.3e-6
    301    base  31.4, 28.8  -0.49, -0.80        NaN, NaN (warned)
    301    lane  80.8, 75.7  -1.2e-7, -5.2e-7    within 2.0e-6

The lane fit at 159 parameters prints, on install 21, three lines
"Error in (function (ptr) : promise already under evaluation" (three
again in each of two reruns, `f4/p4wide-150-re{a,b}.txt`; none with the cap at
100, `f4/p4wide-150-cap100.txt`; not on install 20). They come from
TMB's ADFun finalizer, `function(ptr)` in `TMB:::registerFinalizer()`,
run by a garbage collection while the lazy-load promise of
`TMB:::isNullPointer` is being forced: with that promise forced at
TMB's load, the messages are gone (`f4/p4wide-150-forced.txt`, through
`R_PROFILE_USER=f4/force-tmb.Rprofile`). The fit and its numbers are
the same; the three ADFun objects are not freed. It is TMB's, and the
route's many tapes only make the collection likelier; not changed.

The cap at 400 is a cost judgment: the route's two Hessians cost 2 k
gradients each, and the eigen decomposition k^3. Above 400 the fits keep
0.71.0's behavior.

## 5. Item 5: the cost

Profiled first (`dev/illse2-cost.R`, Rprof at 0.02 s, and
`dev/illse2-cost2.R`, which times every Hessian with its caller; the
reviewer's design: n rows, 100 outer parameters, with and without
`(1 | g)`), at n = 1e4 against 0.70.0 (`rellib-r8`) and 0.71.0:

- 0.70.0 spends 85 percent of the fit in nlminb and 8 percent in the
  gradient check. 0.71.0 adds newton_polish(): 35 percent of the fit
  without random effects (Rprof 22.5 s of 64.4) and 29 percent with
  them (33.2 of 114.3); the SE check itself is 0.1 percent (it reuses
  the cached Hessians).
- Inside it, two Hessians in the design's coordinates, one at the stop
  and one at the end, which conv_headroom() and the SE route then
  reuse from the cache: 10.4 s and 9.9 s without random effects, 25.2
  s and 27.1 s with them (`cost2-lane-*-1e4.txt`). Each is 4 k
  gradients (Richardson over steps 1e-3 and 2e-3).

So the change is in the difference scheme, not the allocations: the
time is in `obj$gr()`, 26 ms per gradient (62 ms with `(1 | g)`), and
the R work around it is a 100 x 100 product per column.

Schemes, accuracy first (`dev/illse2-acc.R`; A to E on installs 7 and 8
with a temporary option, `acc-summary.txt`; F on the final build,
`f2/acc-summary.txt`;
every entry of `vcov(full = TRUE)` against the exact covariance at the
fit's point, the largest over seeds 1 to 3, n 300):

    scheme                       family         deg 6    deg 7    deg 8
    A  Richardson, every column  gaussian       1.7e-8   1.6e-7   3.9e-6
                                 gaussian(1|g)  2.3e-8   4.6e-7   5.0e-6
                                 bernoulli      8.0e-8   1.2e-6   1.6e-5
                                 poisson        4.6e-8   6.5e-7   5.3e-6
    B  central, step 1e-3        gaussian       6.7e-7   6.7e-7   1.5e-6
                                 gaussian(1|g)  6.7e-7   6.7e-7   4.6e-6
                                 bernoulli      3.8e-7   5.5e-7   1.0e-5
                                 poisson        6.6e-7   1.0e-6   5.0e-6
    C  central, step 1e-4        gaussian       1.7e-7   2.3e-6   1.3e-5
                                 gaussian(1|g)  3.5e-7   4.0e-6   6.4e-5
                                 bernoulli      7.8e-7   4.8e-6   3.1e-4
                                 poisson        3.3e-7   5.4e-6   4.9e-5
    E  Richardson only for a     gaussian       8.7e-9   9.0e-8   1.8e-6
       column of one parameter   gaussian(1|g)  2.3e-8   7.2e-7   3.8e-6
                                 bernoulli      2.9e-7   4.5e-7   1.7e-5
                                 poisson        6.6e-7   8.1e-7   4.5e-6
    F  Richardson but on an      gaussian       8.8e-9   3.0e-7   3.2e-6
       identity-link gaussian   gaussian(1|g)  4.8e-8   1.6e-6   1.2e-5
       design (final)           bernoulli      8.0e-8   1.2e-6   1.6e-5
                                 poisson        4.6e-8   6.5e-7   5.3e-6

B's floor at 6.7e-7 is the truncation along log sigma; C is limited by
the rounding of the gradient. E keeps A's accuracy on gaussian fits
(the gradient is linear along a design's column, so a central
difference is exact there) and B's on the others, at half the
gradients. E was the scheme of install 17; the review showed it loses
digits along a nonlinear body's Jacobian and on non-gaussian designs
(section 12, finding 3), so the final scheme F extrapolates every
column but those of an identity-link gaussian design: bernoulli and
poisson are A's numbers exactly. Gaussian without random effects is
within 3.2e-6; with `(1 | g)` it is 1.6e-6 and 1.2e-5 at degrees 7 and
8 against A's 4.6e-7 and 5.0e-6 (the cause was not examined further),
the price of the gaussian saving with random effects.

Timing of the final scheme F (install 21, 2026-10-19;
`dev/illse2-time.sh <re> 1e4 2 f4 "base lane"`, one fit at a time,
arms interleaved, two rounds, n = 1e4, `time-f4-<re>-1e4.txt`; wall
seconds, CPU within 1 percent; two other test files of the full scan
ran beside it):

    design     0.71.0 (Richardson)   lane (F)
    no RE      80.1, 67.5            61.4, 58.4
    (1 | g)    149.4, 132.2          118.3, 120.8

Taking the faster round of each arm, the lane takes 0.87 of 0.71.0's
time without random effects (58.4 against 67.5 s) and 0.89 with
`(1 | g)` (118.3 against 132.2 s).

n = 1e5, the reviewer's own size (`dev/illse2-time1e5b.sh f4`, the four
fits at once beside the core one-session run and R CMD check, one run
each, `t1e5-f4-*.txt`):

    design     0.71.0               lane                 lane / 0.71.0
    no RE      855.8 s (cpu 851)    754.3 s (cpu 750)    0.88
    (1 | g)    1673.4 s (cpu 1664)  1471.2 s (cpu 1463)  0.88

Every fit ends at code 1 with the convergence warning that says the
stop is the maximum, nothing lost, the factor held, and the same number
of evaluations on both arms.

0.70.0 is no longer on the machine (`rellib-r8` was removed during the
round), so the final scheme is compared with 0.71.0 only. The earlier
runs with 0.70.0 (install 17, scheme E, which takes the same
differences as F on these gaussian designs: every design column
central, log sigma and the sd extrapolated) gave, at n = 1e4 and on the
faster round, 0.70.0 38.6 and 85.3 s, the lane build with 0.71.0's
scheme A (every column extrapolated; arm "lane" of that script, not
rellib-r9) 61.7 and 135.3 s, the lane with scheme E 53.4 and 117.0 s
(`time-e1-<re>-1e4.txt`): 1.38 and 1.37 times 0.70.0 where scheme A
took 1.60 and 1.59. At n = 1e5 the install-17 run gave
0.77 and 0.89 of 0.71.0 (`t1e5-f1-*.txt`), against 0.88 and 0.88 here;
both were taken on a shared machine, one run each.

What was not changed, and why: the natural Hessian newton_polish()
takes at the stop to decide whether the design's coordinates are needed
(an AD Hessian without random effects, which the SE check reuses;
optimHess(), 2 k gradients, with them) also finds the rows a step must
leave alone (an sd at zero); deciding from the design alone would skip
that, and the sd-at-zero fixes of lane illse rest on it.

## 6. Item 6: a correlation's boundary verdict on a flat ridge

Assessed, not implemented. On `rs_sd0.02` seed 4 the two stops are at
the same maximum to 6e-8 (release review, section 2); the verdict
(boundary message, or theta_3's SE 1.86e4 on t = 19.66, a delta-method
SE of rho of about 2.4) depends on where along the ridge the stop
lands. Both answers say the correlation is not estimated; neither is a
finite wrong number of the kind this lane fixes. A pull in correlation
space would change verdicts across the false-alarm battery's `rs_sd`
models, and there is no measured case where the current verdict misleads
a user, so the evidence the brief asks for is not there.

## Defects found, not fixed

- newton_polish() turns off whenever `optCtrl` names any tolerance, so
  a tighter `rel.tol = 1e-12` leaves degree 7 seed 4 of the 80
  sd-at-zero fits 0.211 short, with a warning and every coefficient SE
  NaN, on 0.71.0 and on the lane (review finding 2,
  `f4/esc2-lane.txt`). A rule like optctrl_stops_earlier() would fit
  there too; filed, not changed.
- TMB's ADFun finalizer errors with "promise already under evaluation"
  when a garbage collection runs it while `TMB:::isNullPointer`'s
  lazy-load promise is being forced; the ADFun is then not freed. Met
  on the lane's 159-parameter fit (section 4); not frmtmb's code.
- (Being handled by lane adhess this round, findings 7 and 8 of the
  review.) A no-RE fit whose SE check finds nothing keeps sdreport()'s
  finite-difference covariance (optimHess, absolute step 1e-3). On
  `y ~ a * exp(b * x)` with x on 20 + [0, 1] (kappa 2.7e4, not
  ill-conditioned) it is 0.133 from the exact covariance on the scale of
  the correlation and the SEs 0.6 to 0.9 percent off, silently, where
  solve() of the exact AD Hessian is within 8.3e-13
  (`dev/illse2-p3exp-dbg.R`). On the rational body of section 3 (seed
  3) the shown covariance is 0.47 off where solve() of the AD Hessian
  is 1.2e-2 off (`p3ratdbg-lane-3.txt`). A fix would read the AD
  Hessian into sdreport() on every fit without random effects, which
  moves every such fit's covariance in its last digits; that is a
  lane of its own: lane adhess (review finding 7, which measures it on
  many more models). ML on a shifted covariate stopping short (review
  finding 8) is also lane adhess's.
- An active bound on an ill-conditioned design whose least-squares
  start is off (condition number above start-ls's 1e5) starts with the
  bound coefficient at its bound and the others at zero, at a saddle,
  and stops there 2580 to 2700 short with a warning (section 4, degree
  8 seeds 1 to 3). Lane optstart's code.
- The REML `beta` whitener needs the joint precision's dense factor;
  above `joint_root_max` (4000) joint parameters it is not formed, and
  such a fit keeps 0.71.0's inverse of the whole joint precision.
- `joint_cov_repair()` (a REML fit whose variance parameter is lost,
  `cov_fixed_prop` set) still inverts the inner block `Qrr`, which under
  REML holds `beta`; not met on the designs here.

## 7. Tests

New files and tests, each seen to fail on rellib-r9 (one file per
process, `dev/ciharden-run1.R base`, 2026-10-19 on the final test
files; logs `f4/t-base-<file>.txt`; lane counts from `scan/f4-ref.sum`,
install 21):

    file / test                                    base             lane
    test-se-reparam-inv.R (new, 4 tests)           2 pass, 4 error  13 pass
    test-se-reparam-reml.R (new, 4 tests)          8 pass, 10 fail, 19 pass
                                                   1 error
    test-se-reparam.R: nonlinear body (rewritten), 56 pass, 16 fail, 76 pass
      sd below a small interior maximum, active    2 error
      bound, raised caps, asymmetry pass and
      refusal (new)
    test-se-tier.R: an sd near zero (rewritten)    89 pass, 6 fail  95 pass

The base errors on test-se-reparam-inv.R: `solve.default(V)` and
`solve.default(U[, , j], Q[, j])` "system is computationally singular"
(the behavioral failures), and the two unit tests of the new helpers
(missing symbols, the weak form). The REML file's base error is the
missing `beta_whitener` of its guard-absent test; its ten failures are
behavioral (2 warnings, no factor, SEs 0.497 and fitted SEs 0.164 off;
no whitener on the four log-likelihood fits, whose logLik is 4.6e-8 and
2.3e-5 relative off the restricted likelihood at degrees 5 and 6).
test-se-reparam.R's 16 base failures: 5 for the nonlinear body, 3 for
the sd (a warning, 1.98e-5 relative short, NaN SEs), 4 for the bound
(no bound message, nothing lost, 2.7e-3 relative short, SEs 35 times
off), 2 for the raised caps (a warning, no escape), 2 for the
asymmetry pass (no factor, SEs 3.5e-4 off); its 2 errors are
`optctrl_stops_earlier` and `reparam_basis(est =)` missing (the weak
form; the refusal test's behavioral failure is the guard-off run of
section 12, finding 4). test-se-tier.R's 6: the two seeds, each a
warning, 6.1e-3 and 4.7e-3 relative short of lme4 (`expect_lt` on
1e-8), no escape.

Each fix of the review undone in the namespace of install 21
(`dev/illse2-fx-mock.R`, `dev/illse2-fx-nojac.R`; logs `f4/mock-*.txt`,
`f4/t-reml-nojac.txt`):

    undone                       test that fails                 how
    REML Jacobian (finding 1)    REML logLik (degrees 4 to 6),   6 failures:
                                 shifted covariate               0.180, 0.286,
                                                                 0.440 off the
                                                                 reference,
                                                                 1.18e-3 between
                                                                 the shifts
    escape rule of install 17    raising the optimizer's caps    4 failures: a
      (finding 2)                                                warning, no
                                                                 escape, 1.98e-5
                                                                 short, the rule
    scheme E (finding 3)         identified nonlinear body       SEs 1.69e-4 off
                                                                 (limit 1e-4)
    asymmetry guard, Inf         nearly unidentified body        route taken
      (finding 4)                (A = 0.05, seed 2)              (factor held)
    asymmetry guard, 0           every route test (allguard)     44 failures
                                                                 and errors
                                                                 (testthat
                                                                 prints 10)

Guards' absent cases: the well-conditioned design of
test-se-reparam-inv.R keeps `solve(V)` bit for bit (`expect_identical`);
the degree-2 REML fit has no whitener and no factor; the asymmetry and
the sd escape are measured on the batteries (section 8).

## 8. Suites, OpenBLAS, batteries

### Final build, install 21 (2026-10-19)

All eight packages, one file per process, every gate on, reference
BLAS (`bash dev/illse2-scan.sh f4-full ref <lane lib> 6`,
`scan/f4-full.sum`): 384 of 384 files ran, and every file's `lib:` line
names the lane's frmtmb. frmtmb.ode's `test-scale.R` ran (3131 s).

    package          files   pass  fail  error  skip  warn
    frmtmb             238  20013     0      0     0     0
    frmtmb.coupling     11    548     0      0     0     0
    frmtmb.eam          30   1826     0      0     0     0
    frmtmb.latent       10    366     0      0     0     0
    frmtmb.learn        15    515     0      0     0     0
    frmtmb.ode          11    551     0      0     0     0
    frmtmb.sample       53   2916     0      0     0     0
    frmtmb.spline       16    596     0      0     0     0

The 44 files that read the covariance, the SE check, the optimizer's
stop, REML, nonlinear bodies or bounds (`ob-subset.jobs`;
`dev/illse2-f4-scans.sh`, scans `f4-ref`, `f4-ob26`, `f4-ob32`, 4
threads under OpenBLAS):

    BLAS                   files   pass   fail  error  skip  warn
    reference                 44   3697      0      0     0     0
    OpenBLAS 0.3.26 LAPACK    44   3697      0      0     0     0
    OpenBLAS 0.3.32 LAPACK    44   3697      0      0     0     0

These ran before the refusal test of `test-se-reparam.R` moved to the
A = 0.05 body (section 12, finding 4). The final file then passed 76
of 76 under all three (`dev/illse2-f4-reparam.sh`, scans `f4b-ref`,
`f4b-ob26`, `f4b-ob32`), and the full scan above ran the final file.

False-alarm battery of lane illse (5 model groups, 41 models, seeds 1
to 20) and its bitwise set, lane against rellib-r9, reference BLAS
(`TAG=f4 CFGS=ref N=2 SEEDS=1:20 bash dev/illse2-sf-run2.sh`,
`sf-f4/falarm-cmp.txt`, `sf-f4/bitcmp-ref-lane-vs-base.txt`): 820 fits
per arm, 0 rows differ in code, lost set, SE warnings, other warnings
or the largest |SE / ref - 1| as a string; the bitwise set is 30 of 30
fits `identical()` to rellib-r9 in `fixef()`, `vcov()`, `fitted()` and
predictions.

Core's whole suite in one session (`dev/illse2-onesession.R`,
`f4/onesession.txt`; NOT_CRAN, gated tiers off, every extension
installed): 228 files, 17490 pass, 9 fail, 0 error, 185 skip, 0
warning. The 9 failures are test-compat-register.R (5) and
test-compat.R (4), the one-session failures known on 0.71.0 that lane
testiso fixes; none is in a file this lane touched.

### Install 17 (the record before the review)

All eight packages, one file per process (`dev/illse2-scan.sh`, every
gate on, `scan/s1-ref.sum` on install 15 and `scan/s2-ob26.sum` on
install 16; `scan-summary-s1-s2.txt`; every log's `lib:` line names the
lane's frmtmb):

    package          files  ref pass  fail  OB 0.3.26 pass  fail
    frmtmb             238     19990     2           19993     0
    frmtmb.coupling     11       548     0             548     0
    frmtmb.eam          30      1826     0            1826     0
    frmtmb.latent       10       366     0             366     0
    frmtmb.learn        15       515     0             515     0
    frmtmb.ode          11       551     0             549     0
    frmtmb.sample       53      2916     0            2916     0
    frmtmb.spline       16       596     0             596     0

No error, skip or warning in any file of either run. The two reference
failures were test-bracket-access.R reading `bk$idx` in the new helpers
(hazard containers); install 16 uses `[[ ]]` and the file passes 33 of
33. Under OpenBLAS 0.3.26 frmtmb.ode's `test-scale.R` (49 min under the
reference BLAS, 2 pass) was stopped, by me, to install the asymmetry
guard; the other 383 files ran.

Install 17 adds `se_reparam_asym`. The 44 files that read the
covariance, the SE check, the optimizer's stop, REML, nonlinear bodies
or bounds (`dev/illse2-log/ob-subset.jobs`), on install 17 (scans
`s3-ref`, `s3-ob26`, `s3-ob32`):

    BLAS                   files   pass   fail  error  skip  warn
    reference                 44   3676      0      0     0     0
    OpenBLAS 0.3.26 LAPACK    44   3676      0      0     0     0
    OpenBLAS 0.3.32 LAPACK    44   3676      0      0     0     0

The false-alarm battery of lane illse (5 model groups, 41 models,
`dev/illse2-sf-run.sh`, `dev/illse2-sf-run2.sh`) and its bitwise set,
lane against rellib-r9, install 17 (`sf-a1/`, seeds 1 to 20):

- Reference BLAS, 820 fits: 0 rows differ in code, lost set, SE
  warnings, other warnings or the largest |SE / ref - 1| as a string.
- OpenBLAS 0.3.32: 0 differ in code, lost set or warnings; the largest
  |SE / ref - 1| differs as a string on 2 rows (`poly5_lmer` seeds 5 and
  9: 4.934e-6 against 4.935e-6, 1.767e-6 against 1.768e-6).
- Bitwise set: 30 of 30 fits `identical()` to rellib-r9 in `fixef()`,
  `vcov()`, `fitted()` and predictions, under both BLAS.

Seeds 21 to 100 (`sf-a2/`, 3280 fits per BLAS; the release review's
100 seeds), reference BLAS and OpenBLAS 0.3.26 LAPACK: 0 rows differ in
code, lost set, SE warnings or other warnings under either; the largest
|SE / ref - 1| differs as a string on 5 (reference) and 6 (0.3.26) rows,
all `poly6_lmer`, in the fourth significant digit (for example
8.451e-4 against 8.450e-4). Bitwise set: 30 of 30 under both. So over
4100 fits per BLAS the lane adds no warning and loses no SE against
rellib-r9.

Core's whole suite in one session (`dev/illse2-onesession.R`,
`onesession.txt`; NOT_CRAN, gated tiers off, every extension
installed): 228 files, 17469 pass, 9 fail, 0 error, 185 skip, 0
warning. The 9 failures are test-compat-register.R and test-compat.R,
the one-session failures known on 0.71.0 that lane testiso fixes; none
is in a file this lane touched.

(The install-17 record above stays as the record of that build; the
final build's runs are at the top of this section.)

## 12. Fix round 1 (review `dev/reviews/2026-10-14-illse2.md`)

The reviewer's scripts (`dev/illse2-rev-*.R`) were rerun unchanged on
the lane library; logs in `dev/illse2-log/fx/`. Every accuracy number in
sections 1 to 5 and in NEWS was regenerated on the final build, install
21, by `dev/illse2-f4-jobs.sh` (50 jobs, logs in `dev/illse2-log/f4/`,
summary `f4/summary.txt` by `dev/illse2-f4-sum.R`); they equal install
20's (`dev/illse2-f2-jobs.sh`) to every printed digit.

### Finding 1 (BLOCKING): the REML log-likelihood missed `log |det T|`

Laplace over `gamma` integrates the density over `beta` divided by
`|det T|`. `beta_whiten_nll(nll, bw, jacobian = REML)` now carries
`-log |det T|` (summed over the blocks) in the tape; `profile = TRUE`
maximizes `beta` and takes no term (the reviewer: -591.77142224 on both
arms at every shift).

`dev/illse2-fx-remlll.R` (`f4/remlll-lane.txt`, `fx/remlll-base.txt`):
each REML log-likelihood against the exact restricted log-likelihood at
its own estimates (256 bits, with the `(n - p) log(2 pi)` constant), on
test-se-reparam-reml.R's data, seed 1:

    deg   lane - exact   0.71.0 - exact   lme4 - exact   exact max - lane
    3     -3.8e-10       -3.8e-10         +4.9e-11       2.7e-13
    4     -4.7e-13       +4.6e-8          +6.3e-9        1.4e-14
    5     -6.2e-12       +3.3e-6          -8.0e-7        1.0e-13
    6     +4.7e-11       +1.5e-3          +1.6e-4        8.5e-14
    7     +6.6e-10       +0.33            +3.1e-2        1.4e-14

So the lane's REML log-likelihood is the restricted log-likelihood to
6.6e-10; 0.71.0's own value was off by 1.5e-3 and 0.33 at degrees 6 and
7, and lme4's by 1.6e-4 and 3.1e-2. The reviewer's
`dev/illse2-rev-remlll2.R` (`fx/remlll2-lane.txt`): lane minus lme4 at
degrees 3 to 6 is 0, 0, +1e-6 and -1.4e-4 (lme4's own error at 6).

Consumers of the REML objective value: `anova()` of REML fits is allowed
only between fits with one fixed design (`R/confint.R`), where a
constant cancels; profile intervals require ML; the sd probe and the
headroom read differences. All are exact now anyway.

Tests (`test-se-reparam-reml.R`, 19 expectations on the lane): "REML
logLik of a whitened fit is the restricted likelihood" (degrees 4 to 6,
against the restricted log-likelihood in double through the whitened QR
within 1e-9 relative, and lme4 at degrees 4 and 5 within 1e-7), and "a
shifted covariate leaves the REML logLik where it was" (sleepstudy,
`Days` against `Days + 1e5`). With the Jacobian removed
(`dev/illse2-fx-nojac.R`, `f4/t-reml-nojac.txt` on install 21, as
`fx/t-reml-nojac.txt` on install 17's tape) they fail 6 times:
relative errors 0.180, 0.286, 0.440 against the reference, 0.180 and
0.286 against lme4, and 1.18e-3 between the two shifts.

### Finding 2 (MAJOR): asking for more work turned the sd escape off

`optctrl_stops_earlier()` (`R/fit.R`): the escape is off only for an
`iter.max` or `eval.max` below 1000, a `rel.tol`, `x.tol`, `abs.tol` or
`xf.tol` above nlminb's (1e-10, 1.5e-8, 0, 2.2e-14), or another
optimizer. The reviewer's `dev/illse2-rev-esc2.R` (`f4/esc2-lane.txt`,
degree 7 seed 4):

    control                        short      escape    SE finite   warning
    default                        7.6e-9     0.00119   8 of 8      none
    iter.max = eval.max = 1000     7.6e-9     0.00119   8 of 8      none
    iter.max = eval.max = 5000     7.6e-9     0.00119   8 of 8      none
    rel.tol = 1e-12                0.211      no        0 of 8      conv, short
    restarts = 3                   7.6e-9     0.00119   8 of 8      none
    grad_tol = 1e-4                7.6e-9     0.00119   8 of 8      none

Test: "raising the optimizer's caps keeps the sd escape"
(`test-se-reparam.R`). With install 17's rule put back in install 21's
namespace (`dev/illse2-fx-mock.R esc17`, `f4/mock-esc17.txt`) it fails
4 times: a warning, no escape, logLik 1.98e-5 relative short, and the
rule's own assertion; on rellib-r9 2 failures and 1 error (section 7).
NEWS now says what the rule is. The `rel.tol = 1e-12` row stays 0.211
short, as on 0.71.0: the escape may run, but newton_polish() is off
whenever a tolerance is named, and the escape's probe gains nothing at
that stop. Filed below, not fixed.

### Finding 3 (MAJOR): the cheaper scheme lost digits off gaussian designs

`reparam_hessian()` now extrapolates every column except one that moves
only the coefficients of an identity-link gaussian design without an
aterm (`reparam_lin_lp()`; reparam_basis() returns their positions as
`lin_pos`, and `reparam_cols()` keeps the attribute when callers drop
columns). The refining pass extrapolates every column, as before.

The reviewer's cases on the final build (install 21, `f4/`):

- `dev/illse2-rev-asym.R` (degree-6 polynomial beside `A exp(-k z)`, n
  2000; `f4/asym-1-lane.txt`): A = 1, cor 3.17e-7 and 3.99e-7 (install
  17: 3.38e-4 and 3.85e-4; 0.71.0: 7.09e-4 and 1.78e-3), asymmetry
  2.46e-8 and 3.39e-8.
- `dev/illse2-rev-hess.R` (`f4/hess-*-lane.txt`; 0.71.0 from
  `fx/hess-*-base.txt`, 2026-10-14): the lane equals 0.71.0 to every
  printed digit (cor, SE dev, gain): poisson degree 6 at n = 1e5 cor
  5.72e-8 and 5.23e-8 (install 17: 5.2e-7), bernoulli degree 7 at n =
  1e5 4.67e-7 and 6.19e-7, poisson degree 8 at n = 1e4 3.90e-6 and
  9.98e-7. The script still labels the lane arm "scheme E"; its label
  is a constant. On these fits no column is in an identity-link
  gaussian design, so the lane takes 0.71.0's differences and no time
  is saved; the times in the logs (lane 21.3 to 50.3 s, written beside
  seven other jobs, against 0.71.0's 18.9 to 48.8 s of 2026-10-14) are
  not a comparison. On gaussian fits time is saved (section 5).
- `dev/illse2-p3nl.R sin 1:3` (`f4/p3nl-lane-sin.txt`): 4.92e-9,
  1.81e-8, 4.67e-9 (the reviewer's install-17 rerun gave 5.2e-7,
  4.0e-7, 4.5e-7; the earlier findings table was written before the
  scheme changed and happened to match the final build).

Test: "an identified nonlinear body passes the asymmetry guard" (finding
4) pins the scheme: with install 17's scheme E put back in install 21's
namespace (`dev/illse2-fx-mock.R schemeE`, `f4/mock-schemeE.txt`) its
SEs are 1.69e-4 off the reference, above its 1e-4; with F it passes.
The margin is a factor of 1.7 on E's side; on F's side the route's
covariance is within 3.2e-7 (cor), so the bound of 1e-4 is far from it.

### Finding 4 (MINOR): the asymmetry guard, tested

After finding 3 (`f4/asym-*-lane.txt`, and with the guard off,
`I2_ASYM=off`, `f4/asym-*-lane-off.txt`):

    A     seed  t(k)  asymmetry  guard on              guard off
    1     1, 2  459,  2.5e-8,    taken, cor 3.2e-7,    (same)
                452   3.4e-8     4.0e-7
    0.2   1, 2  92,   1.5e-5,    taken, cor 1.9e-4,    seed 1 the same
                90    2.1e-5     2.5e-4
    0.05  1, 2  23,   4.3e-3,    refused, 0.71.0's     taken, cor 5.9e-2,
                22    6.6e-3     covariance, 2.5e-3,   9.1e-2
                                 5.9e-4 off

At A = 0.2 the route is 16 and 14 times better than 0.71.0's covariance
(3.15e-3 and 3.55e-3, the review); at A = 0.05 the guard keeps a
covariance 24 and 156 times better than the route's. 0.71.0's
covariance there is still off silently: the pre-existing defect of
finding 7. The 27 fits of `dev/illse2-asym2.R` read at most 1.71e-5
(`f4/asym2.txt`), 58 times below the threshold; A = 0.05 is 4.3 and
6.6 times above it.

New tests in `test-se-reparam.R`: "an identified nonlinear body passes
the asymmetry guard" (A = 1, seed 1: the route is taken and the
coefficient SEs are within 1e-4 of a double-precision reference,
`(J'J - S) / s2` inverted through the QR of the Jacobian) and "a nearly
unidentified nonlinear body is refused by the guard" (A = 0.05, seed 2:
no factor, asymmetry above `se_reparam_asym`). Seen to fail:

- With the guard off (`dev/illse2-fx-mock.R noguard`,
  `f4/mock-noguard-b.txt`) the refusal test fails on its behavioral
  assertion (the factor is held) as well as on the threshold.
- With the guard refusing everything (`allguard`, `f4/mock-allguard.txt`)
  the pass test and every other route test fail (44 failures and
  errors).

The refusal test was first written on the rational body of section 3,
seed 3. With the guard off that test still passed its behavioral
assertion (`f4/mock-noguard.txt`: the route refuses that body by its
other checks under the final scheme), so it did not pin the guard, and
it let a warning escape under OpenBLAS 0.3.26 (scan `f2-ob26`). It is
replaced by the A = 0.05 body, which needs no allowed warning; it
passes 76 of 76 under the reference BLAS and OpenBLAS 0.3.26 and 0.3.32
(`scan/f4b-*.sum`). No test sits near the threshold.

### Finding 5 (MINOR): what the REML covariance is

`?vcov.frmtmb_fit` now says: the conditional covariance plus `J
V_theta J'`, glmmTMB's definition (the reviewer: 0.1566099622 on both),
and how lme4's differs (by up to 1e-3 of an SE on an unbalanced design,
not at all on a balanced one). NEWS says it in the REML bullet.

### Finding 6 (MINOR): NEWS

The nonlinear route, the 400 cap and the scheme are under "Breaking
changes" with the sentence that they change SEs. The REML numbers are
install 21's (log-likelihood to 6.6e-10, optimum to 2e-12, covariance
6.1e-9, from `f4/remlll-lane.txt` and `f4/p3reml-lane-<deg>.txt`; NEWS
had 7e-10 and 6.2e-9 until the resumed round); the escape bullet
states the rule; the nonlinear and Cook's-distance numbers are
regenerated, and the timing is the final scheme's (section 5). Two
comments that said how now say why. 81 sentences, none over 25 words,
the longest 24 (`dev/illse2-news.R`, on the final NEWS).

### Findings 7 and 8

Filed by the review and handed this round to lane adhess: sdreport()'s
finite-difference covariance on fits without random effects (wider than
the defect entry below said: any non-gaussian or nonlinear model with a
covariate on a scale above about 10), and ML on a shifted covariate
stopping short. Not worked on here.

### Fix round 1, resumed (2026-10-19)

The round was paused on 2026-10-14 after install 20 (16:29). On
resuming, items 1 to 6 were all in the code, NEWS and `?vcov`, but:
the tests of items 2, 3 and 4 had no recorded failing run; the timing
was install 17's; `p3reml-summary.txt` predated the fix round; NEWS
carried 7e-10 and 6.2e-9 for the REML numbers (the logs: 6.6e-10 and
6.1e-9); and no full scan existed for install 20 or later (`f3-ref`
ran 383 of 384, frmtmb.ode's `test-scale.R` stopped by the pause, and
found 2 test-bracket-access.R failures on `family$family` in
reparam_lin_lp(), which the source already read as `[["family"]]`).

The first attempt at install 21 failed: the user library was hollow
(223 of 405 package directories empty, modified 2026-10-16 23:30:11 to
23:30:17). After the user library was restored in place (223 packages
from CRAN binaries, versions matching the 2026-09-28 manifest), install
21 built (`install-21.txt`), and:

- `dev/illse2-f4-jobs.sh` reran every accuracy job and the four
  mock modes (`dev/illse2-fx-mock.R`) and the no-Jacobian run on
  install 21 (`f4/`); `dev/illse2-f4-sum.R` writes `f4/summary.txt`.
  Every accuracy number equals install 20's to every printed digit;
  the findings and NEWS now cite `f4/`. Three findings numbers that no
  surviving log held were corrected (section 2's 2.38e-6 is 5.63e-6;
  the deg 7 s4 and deg 8 s3 "short" are 7.6e-9 and 6.4e-9; the
  variance share of `J V_theta J'` is 2.1e-9 to 6.2e-3), and three
  NEWS numbers were set to the logs (5.6e-6, 1.8e-8, 5.8e-4).
- The missing failing runs: finding 2's test fails 4 times under
  install 17's rule; finding 3's scheme E puts the A = 1 body's SEs
  1.69e-4 off (limit 1e-4); finding 4's refusal test, moved to the
  A = 0.05 body because the rational body did not pin the guard, fails
  with the guard off (section 7 table, section 12 findings 2 to 4).
- Timing on the final scheme against 0.71.0 (section 5): 0.87 and 0.89
  at n = 1e4, 0.88 and 0.88 at n = 1e5. NEWS now gives these and no
  longer the install-17 comparison with 0.70.0.
- Test plan on install 21 (section 8): the 44 files under three BLAS,
  3697 pass each; all eight suites, 384 of 384 files, 0 fail, error,
  skip or warning; the battery, 820 fits, 0 rows differ, 30 of 30
  bitwise; core in one session, only the 9 known compat failures.
- New observation: TMB's finalizer message on the lane's 159-parameter
  fit (section 4, defects list); not frmtmb's.
- After install 21 two text-only edits were made: NEWS (the numbers
  above and the timing bullet) and one `@noRd` comment in
  `R/se-reparam.R` (reparam_hessian()'s timing, now 58 s against
  0.71.0's 68 s). Install 23 (`install-23.txt`) carries them; R CMD
  check ran on the tarball built from the tree before these two edits
  (section 10).

## 9. Overlaps, files, and what is not done

Files changed: `R/se-reparam.R` (root_inv_quad(), warn_cov_singular(),
reparam_nl_jacobian(), `est` and `drop` of reparam_basis(),
design_whitener()'s inverse, the beta whitener and sdr_unwhiten(),
joint_cov_whole(), square_root_rows(), the difference scheme,
`se_reparam_asym`, `se_reparam_np` 400, `newton_natural_np`),
`R/se-check.R` (se_sd_gain_up()'s `at`, se_sd_scale(), the route
without defer after tier 3, se_bound_held(H)), `R/fit.R`
(sd_up_escape(), est_at(), the beta whitener in fit_assembled(),
bounded Newton steps, bound_held_by_step(), the headroom under bounds),
`R/interop.R` (coef_vcov_root(), the gates read the joint factor,
cov_routes_right()), `R/influence.R`, `R/multiple.R`, `R/confint.R`,
`R/brms-shapes.R`, `R/methods-fit.R`, `R/predict.R`, `R/autoscale.R`,
`NEWS.md`, `man/vcov.frmtmb_fit.Rd` (from `R/methods-fit.R`); tests
`test-se-reparam-inv.R`, `test-se-reparam-reml.R` (new),
`test-se-reparam.R`, `test-se-tier.R`. Scripts and records under
`dev/illse2-*` (the resumed round added `dev/illse2-fx-mock.R`,
`dev/illse2-f4-jobs.sh`, `dev/illse2-f4-sum.R`,
`dev/illse2-f4-scans.sh`, `dev/illse2-f4-reparam.sh`,
`dev/illse2-time1e5b.sh`, `dev/illse2-fin.R`, `dev/illse2-fin2.R`,
`dev/illse2-fin3.R`).

Likely overlaps with other lanes: `R/fit.R` (newton_polish(),
fit_assembled(), the headroom functions), `R/se-check.R`
(cov_from_hessian(), se_sd_gain_up()), `R/predict.R`
(get_joint_cov()), `R/confint.R` (hyp_par_cov()), `R/interop.R`,
`R/autoscale.R`, `NEWS.md`. Lane optstart's start code is not touched;
a start that solves the least squares with a bound coefficient held
would fix section 4's three degree-8 fits.

Not done, and why:

- marginaleffects and a caller's `a' vcov() a` (not this lane's).
- sdreport()'s finite-difference covariance on fits without random
  effects (defects, above).
- Item 6 (section 6).
- The REML factor above 4000 joint parameters; `joint_cov_repair()`
  under REML with a lost variance parameter.
- Newton steps from a saddle under a bound (section 4, degree 8).

## 10. R CMD check

`R CMD check --as-cran` of core, once, on the final source of the
resumed round (tarball built 2026-10-19 23:08 from the tree of install
21 plus the final test-se-reparam.R; `dev/illse2-check.ps1`,
`dev/illse2-check/frmtmb-wt-illse2/frmtmb.Rcheck/00check.log`):
Status: 1 NOTE, "checking HTML version of manual ... Skipping checking
math rendering: package 'V8' unavailable", the known NOTE; tests 16
min, OK; examples OK, with `--run-donttest` OK. No extension changed,
so none was checked. After the build, NEWS text (three numbers and the
timing bullet) and one `@noRd` comment changed; no Rd and no code
changed (install 23 carries them).

## 11. Scratch

The OpenBLAS copies of R (`dev/ciharden-out/`, 925 MB) are deleted
again after the resumed round; `bash dev/ciharden-openblas.sh 0.3.26
lapack` and `0.3.32 lapack` rebuild them. The check tarball is deleted;
`00check.log` stays. No R process of this lane is running. The gated
brms-suite files pass on the lane (scan `f4-full`, every gate on), so
no row of the ported-suite ledger flips.
