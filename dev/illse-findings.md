# Lane illse: findings

Lane `illse` of the round opened 2026-10-13. Worktree
`frmtmb-wt-illse`, branch `wt-illse`, base d2ef6312 (frmtmb 0.70.0).
Lane library `C:/Users/adf44/source/r/wt-illse-lib` over `rellib-r8`;
base build `rellib-r8`. Final build: install 9
(`dev/illse-log/install-9.txt`, core and frmtmb.spline; the accuracy
scripts give the same output on install 8). Logs are local,
in `dev/illse-log/`; every number below is copied from them, with the
script that wrote it. OpenBLAS: `dev/ciharden-openblas.sh 0.3.26
lapack` and `0.3.32 lapack`, in `dev/ciharden-out/`.

Three items: the main task (accurate standard errors on ill-conditioned
designs), item c (the SE check on fits with a convergence warning) and
the optional item b (the flat simplex direction off the coordinates).

## 1. The measurement that decided the design

`dev/illse-floor.R` (`floor.txt`; base R, Rmpfr at 256 bits, a raw
polynomial on [1, 2], n 200, seed 1; kappa is the condition number of
`X'X` at unit diagonal). Against the exact `(X'X)^-1`:

    design          kappa    lm vcov()     a'Va from the    a'Va from  a'Va from
                             entry error   exact V rounded  lm vcov()  a factor
    degree 5        9.46e10  3.0e-12       1.7e-06          4.4e-07    7.5e-12
    degree 6        1.61e13  6.5e-12       8.2e-05          1.1e-04    5.8e-11
    degree 7        3.18e15  4.7e-10       0.0019           0.0044     1.2e-09
    degree 8        4.18e16  3.6e-09       2.8              3.5        1.6e-08
    degree 9        lm drops a column (rank 9 of 10)
    x on 1e6+[0,1]  4.32e13  6.4e-11       0.0017           0.0011     9.1e-10
    x on 1e7+[0,1]  lm drops a column

`a'Va` is the variance of a fitted value at x = 1, 1.25, ..., 2 (the
largest relative error of the five). So:

- At degree 8 `lm()`'s covariance is accurate entry by entry (3.6e-9),
  but NO double-precision matrix can be multiplied by a row of that
  design: the exact covariance rounded to double gives 2.8 times the
  variance of a fitted value, and `lm()`'s own `vcov()` 3.5 times. A
  factor `T` with `V = T T'` keeps it to 1.6e-8.
- This is why punch round 2b of lane seflags failed through every
  consumer that read `V` whole (review B6, B7): a property of double
  precision, which `lm()` shares through emmeans or marginaleffects.
- The brief's "degree 8 and up is singular in double precision" is not
  what was measured. `lm()`'s QR keeps degree 8 at 3.6e-9; degree 9 is
  where `lm()` and frmtmb's frame drop a column, and a degree-9 fit is
  then a degree-8 one. So degree 8 now gets standard errors.

Hence the design: compute the covariance in coordinates where the
Hessian is well conditioned, hold it as a factor, route every consumer
frmtmb controls through the factor, and for a consumer that forms
`J V J'` itself (marginaleffects, a caller's `a' vcov() a`) say so.

Prototype before any package change (`dev/illse-proto.R`,
`proto1-rich.txt`, gaussian without random effects, seeds 1 to 3): the
outer Hessian by Richardson-extrapolated central differences of the
gradient along the columns of `M` (the inverse R factor of the design's
QR) has a condition number of 1.00 to 1.36 at unit diagonal at every
degree from 3 to 8, and the covariance mapped back is within 1.5e-12
(degree 3) to 3.4e-6 (degree 8) of the exact one on the correlation
scale, where `solve()` of the exact AD Hessian was 1.1e-11 to 3.0e-11
(degree 3), 7.8e-6 to 1.1e-5 (degree 5), 7.5e-5 to 6.4e-4 (degree 6)
and singular from degree 7. With `(1 | g)` (`dev/illse-proto3.R`,
`proto3.txt`, before Richardson): 6.7e-7 to 2.9e-6, where 0.70.0's
`vcov(full = TRUE)` was 3.3e-3 to 1.8e-2 at degree 6 and NaN at
degree 8.

## 2. What changed

### The covariance (`R/se-reparam.R`, new; `R/se-check.R`)

- `reparam_basis(frame, pn)`: a basis `M` of the outer parameters. For
  each linear predictor whose coefficients are outer parameters and
  enter through a dense design (not a nonlinear body, not `mo()`), the
  block at its positions is `T` with `X T = Q sqrt(n)`, from the QR of
  the column-scaled design (`design_whitener()`). An ordinal predictor
  has no intercept, and a common shift of its thresholds is one; the
  shift (read off the family's own threshold maps, `ord_shift()`) joins
  the design. Without it a cumulative degree-7 fit kept a curvature at
  1/1.6e11 of the largest in the new coordinates; with it the first
  pass is at 2.2 (`dev/illse-cumdbg.R`).
- `reparam_hessian()`: `M' H M` by central differences of the gradient
  along each column of `M`, Richardson-extrapolated (steps 1e-3 and
  2e-3). `reparam_hessian_refined()` measures once more along the
  eigenvectors of the first pass when that is still above 1e4 at unit
  diagonal, scaling a direction at its own noise by that noise
  (`dev/illse-refine.R`, `refine.txt`: one step from a cumulative stop
  the first pass gave -2.1e-12 against a noise of 3.3e-7).
- `se_reparam_cov()`, asked by `cov_from_hessian()` before any tier.
  It applies when the fit has at most 100 outer parameters, no bound
  holds a parameter and no sd is at its edge, the Hessian at unit
  diagonal has a direction at or below 1e-8 of its largest curvature,
  and every such direction loads at least half on the covered
  coefficients. It accepts when every curvature in the new coordinates
  is above 1e-6 of the largest and above 10 times its own noise, and a
  Newton step gains at most `se_ill_gain` (0.005). It returns tier 4,
  `V = root root'` with `root = M R^-1`, and marks `attr(root,
  "lossy")` when the natural condition number is above 1e12 (products
  of `V` then keep fewer than about four digits, per the floor table).
- `sdr_rescue_hessian()` stores `sdr$cov_fixed_root`; `sdr_sep_lost()`
  zeroes its separated rows.
- Two guards found by the OpenBLAS runs of the SE test files (install
  8): `test-se-tier.R`'s nonlinear alias `a_k + b` (`a ~ 0 + f`, `b ~ 1
  + (1 | g)`) lost nothing under OpenBLAS 0.3.26 and 0.3.32, 4 failures
  each. The design `0 + f` has condition number 1, yet the fit's
  ill-conditioned Hessian sent it to the new coordinates, and there the
  refinement scaled the exactly flat alias by the gradient's noise and
  read that noise as curvature: smallest over largest 0.997 against
  -1.8e-15 under the reference BLAS (`dev/illse-ridgeob.R`). So a
  design whose own condition number (columns scaled) is below 1e4 keeps
  its coordinates (`design_whitener()`), and a direction the refinement
  rescaled must curve the objective itself (`se_obj_curv()`) before the
  covariance is accepted. Both OpenBLAS builds pass the file with them
  (install 9).

### The consumers

Through the factor (`jc_rowvar()`, `jc_quad()`, `lb_rowvar()`): the
joint covariance (`get_joint_cov()`, `joint_cov_repair()`, which forms
the group effects' `J Vf J'` as `(J L)(J L)'` and carries `root`,
`Vrest`), the predict path, `predict_mean_se()`, `ord_prob_se()`,
`conditional_smooths()`, emmeans (which takes the grid route on such a
fit, `G V G'` through the factor), `hyp_prop_var()` (`VarCorr()`,
summary deltas), `hypothesis()` and `frm_multiple()`'s pooled
hypotheses (`hyp_delta_var()`: the expression differenced along each
column of the factor, since a coordinate gradient's rounding is
magnified by the product: 1.8e-4 from `fitted()` before, 3.6e-6 after,
on test-se-reparam.R's polynomial at x = 1.5), `brms_fixef_extra_vcov()`,
`fit_fd_se()` (the output differenced along each column,
`fd_joint_root()`), the parameter draws of `predict()` (`root z`),
`ranef(condVar = TRUE)` and `smooth_edf()` (which read `sdreport()`'s
`diag.cov.random`, NaN on these fits), and frmtmb.spline (`sp_quad()`
in `curve-cov.R`, `curve-deriv.R`, `curve-feature.R`).
`frm_lp_basis()` returns `V_root` and `V_rest`
(`V = V_rest + V_root V_root'`), `frm_joint_cov()` `root` and `Vrest`.
`vcov()` and `frm_joint_cov()` warn once per lossy fit (class
`frmtmb_cov_products`); marginaleffects' `get_vcov()` returns NaN with
a warning of that class.

### Newton steps and the gradient verdict (`R/fit.R`)

`newton_polish()` takes its steps in the design's coordinates when the
natural Hessian at the stop needs them (`hess_needs_reparam()`: a
direction below 1e-8 of the largest, and none curving downward by more
than 1e-3 of it, which is a saddle; test-grad-verdict.R's start far
off has -0.42 against 2.46, and a finite-difference Hessian's noise on
poisson degree 8 with `(1 | g)`, seed 2, stays below the bound). There
a step is damped (Levenberg and Marquardt) while the Hessian is
indefinite and shortened until it gains a tenth of its prediction. It
also runs from a stop nlminb did not call converged, but only in those
coordinates. `grad_headroom()` measures the headroom there too. The
random-effect gate no longer needs code 0.

## 3. Accuracy, final build

`dev/illse-check.R` (gaussian; exact covariance at the fit's own point
in 256 bits: `gauss_exact_cov()`, and `lmm_exact_cov()` for `(1 | g)`,
whose sufficient-statistic form agrees with the per-row form to 9.5e-36
and with base R's dense marginal likelihood to 1.7e-6,
`dev/illse-lmmref.R`, `lmmref.txt`). "cor" is the largest entry error
of `vcov(full = TRUE)` over `sqrt(Vii Vjj)`, every entry; "diag" over
the variances; the fitted columns `fitted(newdata, re_formula = NA)`'s
Est.Error over the exact SE at x = 1, 1.25, ..., 2. Seeds 1 to 5
(`check4-lane.txt`, `check2-base.txt`, summarized by `dev/illse-sum.R`):

    design                 lane cor  lane diag  lane fitted  base cor  base fitted
    degree 5                1.2e-09   6.8e-10   1.00000      1.1e-04   1.00000-1.00002
    degree 5 (1|g)          2.9e-09   1.5e-09   1.00000      7.9e-05   0.99999-1.00001
    degree 6                1.6e-08   1.1e-08   1.00000      2.5e-02   0.99682-1.00471
    degree 6 (1|g)          3.7e-08   2.2e-08   1.00000      1.8e-02   0.99733-1.00179
    degree 7                2.8e-07   1.4e-07   1.00000      1.1e-06   0.99342-1.00214, 1 fit NaN
    degree 7 (1|g)          2.0e-07   1.5e-07   1.00000      1.8e-02   0.99598-1.00985
    degree 8                3.2e-06   3.2e-06   1.00000      5 of 5 lost, NaN
    degree 8 (1|g)          3.8e-06   2.2e-06   1.00000      4 of 5 lost; 1.0, 0.89-1.004
    B5 (8 + z + w)          6.5e-06   6.5e-06   1.00000      5 of 5 lost, NaN
    B5 (8 + z + w) (1|g)    5.4e-06   3.0e-06   1.00000      5 of 5 lost, NaN

No lane fit lost an SE or warned. The lane fits reach the maximum: the
shortfall to `lm()` is at most 1.0e-6 without a group term; with
`(1 | g)` frmtmb ends at or above `lmer()`'s log-likelihood (9.4e-3
above it at degree 7, seed 3), and `lmer()` fails at degree 8 on four of
five seeds.

Covariates far from zero, near-collinear pairs, mixed scales (seeds 1
to 3, `check-other2-lane.txt`, `check-other-base.txt`):

    design                      lane cor  lane fitted  base
    x on 1e5 + [0, 1]           6.2e-09   1.00000      s1 at code 1 NaN; s2, s3 158 and 153 short
    x on 1e6 + [0, 1]           4.8e-08   1.00000      144 to 158 short, NaN
    x on 3e6 + [0, 1]           6.7e-07   1.00000      s2 158 short; s1, s3 x dropped
    x2 = x1 + 1e-5 noise        5.0e-10   1.00000      1.9e-04
    x2 = x1 + 1e-6 noise        1.5e-08   1.00000      1.2e-02
    mixed scales                4.9e-09   1.00000      6.3e-04
    x on 1e5 + [0, 1] (1|g)     1.6e-08   1.00000      2.1e-03 (s3); s1, s2 code 1, NaN
    x on 1e6 + [0, 1] (1|g)     1.1e-07   1.00000      lost flat 2, or 134 to 140 short
    x2 = x1 + 1e-6 noise (1|g)  3.8e-08   1.00000      2.3e-02
    mixed scales (1|g)          1.8e-08   1.00000      1.1e-03

(6.7e-7 is the natural route's own truncation on a fit whose frame
dropped x or x2; those fits do not take the new route.) On 0.70.0 a
covariate at 1e5 to 3e6 is an optimizer failure, not only an SE one:
nlminb stops at "false convergence" at its start slope.

Bernoulli, poisson, cumulative, degrees 6 to 8, with and without
`(1 | g)`, seeds 1 to 3 (`dev/illse-glm.R`, `glm5-lane.txt`,
`glm1-base.txt`). "exact" is the Hessian at the fit's point in 256 bits
(no group term); "poly" the same likelihood fitted with `poly(x, deg)`,
its covariance carried to the raw coefficients through the exact
change of basis; "fitted" is `fitted()` (cumulative: every category)
over the poly fit's:

    family      lane exact cor    lane poly cor     lane fitted   lost: lane / base
    bernoulli   3.8e-08 to 1.7e-05  1.0e-05 to 7.4e-04  0.9999-1.0001   0 / 9 of 18
    poisson     1.6e-08 to 4.9e-06  8.2e-06 to 6.0e-05  1.0000-1.0001   0 / 11 of 18
    cumulative  9.4e-08 to 6.4e-06  5.3e-06 to 2.2e-04  1.0000          0 / 9 of 18

On 0.70.0 the cumulative category route gave P(y = k) SEs 1.04 to 6164
times the poly fit's at degree 6 (1.05 to 7043 with `(1 | g)`) and
10009 to 272087 times at degree 7, on fits that lost nothing and warned
nothing: a silent defect of 0.70.0 below the "illcond" verdict (B7 of
the seflags review at lower degree). The poly reference is limited by
its own fit (its point differs from the raw fit's by the optimizers'
tolerances); the exact column is the measurement where it exists.

The design families' verdicts (`dev/illse-survey.R`, 50 designs, seeds 1
to 3, `survey2-lane.txt`, `survey2-base.txt`): on the base 48 of 150
fits lost some SE and 86 had every SE finite; on the lane 0 and 150.

## 4. Every reader of the covariance

`dev/illse-routes.R` (`routes2-lane.txt`, `routes-base.txt`): degree 8,
seed 3, against the poly(x, 8) fit, at x = 1.1, 1.5, 1.9.

    route                         lane, no RE / (1|g)          base
    fitted(re_formula = NA)       1.00000 / 1.00000            NaN, warned
    fitted() at a group           - / 1.00000                  NaN, warned
    frm_linpred(se.fit = TRUE)    1.00000 / 1.00000            NaN, warned
    frm_lp_basis() via V_root     1.00000 / 1.00000            0.88 to 1.0007 from V, silent
    predict(), 4000 draws         0.996-1.003 / 0.997-1.000    plug-in, warned
    conditional_effects() band    1.00000 / 1.00000            NaN, warned
    emmeans() at x                1.00000 / 1.00000            NaN, warned
    emtrends(x)                   1.00000 / 1.00000            NaN, warned
    sigma's SE (summary)          0.99999 / 1.00000            1.004 / 1.063
    VarCorr() sd's SE             - / 1.00000                  0.999
    ranef(condVar = TRUE) condSD  - / 1.00000                  1.009
    marginaleffects predictions() NaN, warned                  NaN, warned

The seflags review's own scripts, copied by `dev/illse-sf-make.R` and
run by `dev/illse-sf-run.sh` on the final build (`dev/illse-log/sf-f3/`,
reference BLAS and OpenBLAS 0.3.26 LAPACK; `sf-f3ob32/`, OpenBLAS 0.3.32
LAPACK; every lane output is byte-identical to the `sf-f2` run of the
build before):

- `whole` (B6; degree 8 through `I(x^k)`, seeds 1 to 12): `fitted()`,
  `emmeans()` and `emtrends()` 0.9999 to 1 on all 12 under both BLAS;
  marginaleffects NaN with a warning on all 12; `A vcov() A'` formed by
  the caller 0.55 to 1.26 (reference) and 0.41 to 1.35 (OpenBLAS), each
  with the products warning: the double-precision limit of section 1.
- `cum` (B7; cumulative seeds 2, 5, 9): `P(y = 1)` SE over the poly
  fit's 1 on every row, both BLAS; nothing lost, no warning.
- `illpred` (B4): coefficient SEs over the exact ones 1 for all nine on
  the five degree-8 fits. The degree-9 rows are the review's erratum
  (x9 is dropped by the frame).
- `illid` (B5): every hypothesis and sigma SE equal to the poly fit's
  to the printed digits on seeds 3, 8 and 10 with `(1 | g)` and 4 and 5
  without.
- bitwise (`bitcmp-*-lane-vs-base.txt`): 30 of 30 identified fits
  `identical()` to 0.70.0 in `fixef()`, `vcov()`, `fitted()` and
  predictions, under all three BLAS: the guard's absent case.
- false-alarm battery (41 models, 820 fits per BLAS,
  `dev/illse-sf-falarm-cmp.R`): 0 differences in code, lost set, SE
  warnings or other warnings under any of the three BLAS. The largest
  |SE / reference - 1| differs as a string on 80 rows under each BLAS.
  One is `negbin_large_shape` seed 8, 0 on the lane against NaN on the
  base: an item c effect (a stop at the maximum now judged, which the
  base skipped). The other 79 (reference), 78 (OpenBLAS 0.3.26) and 78
  (0.3.32) are raw polynomials of degree 4 to 6
  (`dev/illse-falarm-rel.R`). The lane is closer on 51, 58 and 58 of
  them, and its largest values are lower in every model (degree 6 `lm`:
  0.0177 against the base's 0.0249 and 0.0370; degree 6 `lmer`: 0.0044
  against 0.0120 and 0.0133). Against `lm()` the lane's value is the
  ML-against-REML factor, constant per model (0.0152 at degree 5,
  0.0177 at degree 6); the base's scatter around it is its own error.

## 5. Item c: the SE check on fits with a convergence warning

What changed (`R/fit.R`, `R/se-check.R`):

- `check_convergence()`: on a fit nlminb did not call converged,
  `conv_headroom()` measures what a Newton step from the stop gains over
  the parameters no bound holds (one Hessian, natural or in the design's
  coordinates). `newton_room()` steps over the directions above
  `se_flat_tol` and counts a flatter one only through what it could
  gain, so an exact alias does not make the answer unusable. The warning
  ends with one of three clauses (`conv_headroom_clause()`): within
  `grad_tol`, "the estimates are at the maximum and their standard
  errors are judged as for a converged fit"; above it, "those of a point
  short of the maximum"; no usable Hessian, "the standard errors are
  those of the stop".
- A stop at the maximum whose only warning is that one is judged by the
  SE check as a converged fit is, messages and all. Every other fit with
  a convergence warning keeps `sdreport()`'s standard errors with the
  reason stated, as before (`fit$cache$se_explained` "convergence").
  That is the "kept SEs with a stated reason" answer to the brief's
  question; the next paragraph says why not NaN.

Why a short fit is not judged: on test-review-v29.R's correlation run
to 1 (code 1, curvature unusable) the check called every parameter of
an identified model flat, from finite differences whose asymmetry was
320 (`fd_hessian()` noise); a curvature verdict at such a point is not
about the data. The first form of item c judged every fit and failed
four test files that way (`test-se-check.R`, `test-review-v29.R`,
`test-bcm-mpt.R`, `test-brms-suite-methods.R`'s ledger row) and let two
warnings escape (`test-brms-families.R`, `test-optima.R`), in scan
`f1-ref`, which was stopped and discarded; all pass unchanged under the
final form.

Measured first. The probe `dev/illse-run1scan2.R` (the suite runner with
a hook on `check_convergence()`) records, for every fit of the eight
suites that warns, what a restart of nlminb from the stop gains and
what `conv_headroom()` says. Final build, reference BLAS, all 374
files (`dev/illse-log/convscan-f3/`, 177 probe rows from 54 files, 0
probe errors; `dev/illse-convsum.R`, `convsum-f3.txt`; the build before
gave the same within two fits, `convsum-f2.txt`):

    warning kind      fits   at max by restart   short   no restart
    code != 0          129   104                 24      1
    large gradient      13   5                   7       1

    code != 0, conv_headroom()   restart: at max   short   none
      within grad_tol            34                0       1
      above it                   3                 6       0
      not measured (NA)          67                18      0

So 104 of the 129 stops nlminb did not call converged are at the
maximum by a restart (gain at most 1e-4), and the headroom calls 35 of
them at the maximum with no short fit among them. It is conservative:
3 fits at the maximum read above `grad_tol`, and 67 are not measurable
(mostly singular Hessians on ridges, which keep 0.70.0's treatment).

The rounding flip (`dev/illse-xbeta-perturb.R`, lane ci691's xbeta fit,
phi 2e4, 25 starts a few ulps apart; counted from the logs
`xbeta-*-ref.txt`, `xbeta-base-ob0.3.2*.txt`, `xbeta2-lane-*.txt`):

    build                   code 0: lost / kept   code 1: lost / kept
    base, reference         8 / 0                 0 / 17
    base, OpenBLAS 0.3.26   9 / 2                 0 / 14
    base, OpenBLAS 0.3.32   9 / 2                 0 / 14
    lane, reference         8 / 0                 15 / 2
    lane, OpenBLAS 0.3.26   9 / 2                 11 / 3
    lane, OpenBLAS 0.3.32   9 / 2                 11 / 3

On the base the code decided the verdict: no code-1 fit was judged. On
the lane every code-1 fit is at the maximum (the headroom is 3.3e-6 on
start 6, `dev/illse-xbeta-k.R`) and judged. The two that keep their SEs
stopped where the unit-diagonal eigenvalue ratio is 3.2e-9, above
`se_flat_tol` (1e-9): the remaining split is the check's threshold on
where along the ridge the fit stopped, which ci691 documented, not the
code: under OpenBLAS the lane keeps 2 of 11 code-0 fits and 3 of 14
code-1 fits, the base 2 of 11 code-0 fits and all 14 code-1 fits.

The guard's absent case: converged fits keep their verdicts bit for bit
(the bitwise set and the false-alarm battery, section 4).

Newton steps from a coded stop: `y ~ x` with x on [1e5, 1e5 + 1]
stopped at "false convergence" 157 below `lm()`'s log-likelihood at its
start slope on 0.70.0 and reaches it to 5.4e-14 now
(`dev/illse-offdbg.R`; the check-other table, seeds 1 to 3).

## 6. Item b: the flat simplex direction off the coordinates

`R/se-check.R`: `se_mo_gaps()` lists the interior categories of each
`mo()` term that no row it reads is in; `se_mo_null()` carries each
split `w_c - w_(c + 1)` to the chart's coordinates through the chart's
Jacobian at the estimates and adds it to `se_design_null()`'s exact
flat directions; `se_tier3()` computes them whenever a gap exists and
removes a kept eigen-direction that lies in them. On the residue fits
the Hessian curved by 3.1e-9 of its largest along the split (seed 15 of
`ls ~ mo(inc)`, `dev/illse-segap1.R`), above `se_flat_tol`, and the
ridge rule read 4.6e-9 against its threshold 2.9e-9
(`dev/illse-segap2.R`).

The consolidation's four designs, 40 seeds each (`dev/illse-segap.R`,
`segap-base.txt`, `segap3-lane.txt`):

    design             base: named / silent / other    lane
    ls ~ mo(inc)       38 / 2 (seeds 15, 19) / 0       40 / 0 / 0
    ls ~ mo(inc) + age 38 / 2 (19, 26) / 0             40 / 0 / 0
    ls ~ mo(inc) * age 37 / 0 / 3                      37 / 0 / 3
    y ~ mo(x2)         35 / 4 (3, 10, 18, 31) / 1      39 / 0 / 1

The four "other" lane fits carry the large-gradient warning (short of
the maximum: seeds 2, 30, 37 and 17) and keep their SEs with that
reason, as on the base. Seed 10 of `y ~ mo(x2)` had zeta SEs 547.9 and
949.0 in silence on the base (`dev/illse-segap1.R base 10 4`).

## 7. Tests

New: `tests/testthat/test-se-reparam.R` (7 tests), `test-se-conv.R` (3),
`extensions/frmtmb.spline/tests/testthat/test-illcond.R` (1). Edited:
`test-se-flags.R` (the "illcond" test is now "a direction below double
precision is resolved, not lost"; the B5 test asks for nothing lost;
the simplex-face test gains seeds 15, 19 and `y ~ mo(x2)` seed 10),
`test-optrap.R` ("Newton steps ending where the Hessian does not factor
are undone" becomes "Newton steps on an ill-conditioned design reach
the maximum": in natural units the steps were undone and the fits
stayed 0.674 and 1.73 short, the filed 0.70.0 item; they now reach the
poly fit's log-likelihood).

Seen to fail on 0.70.0 (`rellib-r8`, the same test files, one file per
process; logs `dev/illse-log/t-*-base.txt`, `spl-oldext.txt`):

    file                                   pass  fail  error
    test-se-reparam.R                         8    25      1
    test-se-conv.R                           14     3      0
    test-se-flags.R                          54    14      0
    test-optrap.R                           143     4      0
    frmtmb.spline test-illcond.R              0     0      1

The spline row is the lane core with the 0.70.0 frmtmb.spline
(`dev/illse-spl-oldext.R`): the old extension forms `C V C'` from the
whole matrix and refuses the fit ("disagrees with frm_linpred(se.fit =
TRUE) by 1 relative").

Item c, the tests whose warnings change. In the final suite scan 129
fits in 42 test files stop with a nonzero code; each of their
convergence warnings now ends with one of the three clauses. 35 of
those fits, in 18 files, are at the maximum by `conv_headroom()` and
now go through the SE check: `test-rp-floored.R` (frmtmb.spline, 6),
`test-lca.R` (frmtmb.latent, 4), `test-bcm-simple.R`,
`test-numerical-robustness.R`, `test-thres-refit.R`,
`test-xbeta-zibb-hurdle-cum.R` (3 each), `test-ordinal-mixture.R` (2),
and one each in `test-allfit-start.R`, `test-frailty.R`, `test-fuzz.R`,
`test-hmm.R`, `test-mixture-start.R`, `test-ode.R`,
`test-ordinal-mixture-brms.R`, `test-scale.R` (frmtmb.ode),
`test-se-conv.R`, `test-se-reparam.R` and `test-sparsex.R`. The probe
does not record whether the check then names a parameter; no test
expectation outside the new and edited files changed, and no warning
escapes in any file (below).

## 8. Suites, OpenBLAS and R CMD check

Final build (install 9, `dev/illse-log/install-9.txt`).

Reference BLAS, all eight packages, one file per process, five at a
time (`dev/illse-scan.sh`, `dev/illse-log/scan/f3-ref.sum`,
`dev/illse-scansum.R f3-ref`): 374 of 374 files, every one loading
frmtmb from the lane library.

    package          files   pass  fail  error  skip  warn
    frmtmb             231  18999     0      0     0     0
    frmtmb.coupling     11    548     0      0     0     0
    frmtmb.eam          29   1765     0      0     0     0
    frmtmb.latent       10    366     0      0     0     0
    frmtmb.learn        15    515     0      0     0     0
    frmtmb.ode          11    551     0      0     0     0
    frmtmb.sample       51   2806     0      0     0     0
    frmtmb.spline       16    596     0      0     0     0
    total              374  26146     0      0     0     0

OpenBLAS 0.3.26 and 0.3.32 LAPACK builds of R (`dev/ciharden-out/`),
the 40 files that read the covariance, the SE check or the optimizer's
stop (`dev/illse-log/ob-subset.jobs`; scans `f3-ob26`, `f3-ob32`): 40
of 40 files, 4003 pass, 0 fail, error, skip or warning on each build.
The SE-check batteries (`test-se-check.R`, `test-se-flags.R`,
`test-se-tier.R`, `test-se-conv.R`, `test-se-reparam.R`,
`test-grad-verdict.R`) are among them and fire as on the reference
build.

R CMD check --as-cran (`dev/illse-check.ps1`, `dev/illse-check/`):

    frmtmb          Status: 1 NOTE   (HTML manual: V8 unavailable;
                                      tests 21 min, OK)
    frmtmb.spline   Status: 1 NOTE   (the same NOTE; tests OK)

## 9. Overlaps and what is not done

Overlaps with other lanes: `R/fit.R` (`newton_polish()`,
`check_convergence()`, `grad_headroom()`, `fit_assembled()`),
`R/se-check.R`, `R/predict.R`, `R/interop.R`, `R/confint.R`,
`R/brms-shapes.R`, `R/sandwich.R`, `R/insight.R`, `NEWS.md`,
frmtmb.spline's `R/curve-*.R` and
`test-optrap.R`. Lane optstart changes starting values: a changed start
changes which fits stop with a nonzero code and where, so the item c
counts above (129 coded, 35 judged) move with it, and a fit it moves
from "short" to "at the maximum" becomes judged. The policy itself does
not depend on the start.

Not done:

- `marginaleffects` and any `a' V a` formed by the caller from `vcov()`:
  impossible in double precision on these designs (section 1); NaN or
  the products warning.
- The influence measures (`cooks.distance()`, `solve(V)`) and the D1
  and D2 tests of `anova()` on `frm_multiple()` fits (`solve(U)`) read
  the whole matrix; at degrees 7 and 8 they stop with an error (review,
  attack 3). `vcov_cluster()` is converted in fix round 1 (section 10).
- REML; nonlinear formulas keep the "ill-conditioned" verdict.
- frmtmb.sample's jitter is not scaled by the covariance (it uses only
  the diagonal SEs).

## 10. Fix round 1 (review `dev/reviews/2026-10-13-illse.md`)

The rule of the round: no route gives a finite wrong answer where
0.70.0 gave NaN. The reviewer's scripts (`dev/illse-rev-*.R`) were
rerun on the lane library unchanged (`REV_ARM=lane`); their logs are in
`dev/illse-log/fx/rev-*.txt`. Final build: install 13
(`dev/illse-log/install-13.txt`).

### Item 1: `vcov_cluster()`

`R/sandwich.R`: `cluster_root()` forms `F = sqrt(cf) root K'` with
`K = S root`, reduced to a square factor by a QR of `F'` when there are
more clusters than parameters, and `V = F F'`; the bread is never
multiplied whole. The matrix carries `F` as attribute `root`.
`resolve_vcov_arg()` returns it, and `hypothesis(vcov = )` reads it
through hyp_delta_var(). A matrix without its factor, on a fit held as
a factor, gives NaN with a warning to every hypothesis that combines
parameters; a single coefficient keeps its SE. `vcov_cluster()` reads
the bread through `vcov_estimated()`, so the `vcov()` warning no longer
fires for it; it warns for itself (`warn_cov_products()`, now once per
fit and call).

`dev/illse-rev-lmm.R 101:103 8` (degree 8, `(1 | g)`), against the
reviewer's 256-bit sandwich with the same scores:

    seed   entries, cor scale     hypothesis(vcov = ) SE / exact
           ref        OB 0.3.32   ref          OB 0.3.32
    101    1.46e-5    8.24e-6     1.00000008   0.99999967
    102    1.73e-5    1.36e-5     1.00000017   0.99999919
    103    7.76e-6    1.46e-5     0.99999877   1.00000019

Before the fix (`dev/illse-rev-log/lmm8-lane.txt`): entries 2.33 to
12.3, SE 0, 1228 and 1511 times. Degree 7 (`rev-lmm7-lane.txt`):
entries 5.4e-7 to 1.8e-6, SE within 7e-8 (before: 0.5 to 2.6 percent,
SE 19.9, 0 and 23.5 times). Test: `test-se-reparam.R`, "vcov_cluster()
forms the sandwich through the factor", against lm()'s QR.

### Items 2 and 3: interop methods that hand out the whole matrix

`R/interop.R`, `R/insight.R`. The methods that hand a covariance to
another package: `marginaleffects::get_vcov()`, `insight::get_varcov()`
and `vcov()` itself (read by car, lmtest and sandwich's own callers);
there are no broom methods, and emmeans takes the grid route
(section 4). `get_vcov()` is NaN with a warning on every fit held as a
factor, since marginaleffects' numerical Jacobian is multiplied by the
cancellation even where the products keep their digits.
`get_varcov()` is NaN with a warning where products lose their digits
(`cov_products_lossy()`): its readers form analytic products, and those
keep four digits below `se_product_kappa`. `vcov()` keeps its entries
and warns. Each warning names the routes that are right
(`cov_routes_right()`) and the ones that read the whole matrix.

- `dev/illse-rev-insight.R` (`rev-insight-lane.txt`): degree 8, lossy:
  the matrix is NaN, 1 warning (before: finite, 0 warnings, 0.779 of
  lm()'s SE).
- `dev/illse-rev-glm.R binom 201:204` (`rev-glm-binom-lane.txt`,
  `-ob32.txt`): binomial with mixed scales, kappa 3.7e11 to 4.2e11, not
  lossy: `predictions()` NaN with the warning on all 4 seeds under both
  BLAS (before: 4.8 to 15.1 times the exact SE, silent). Every other
  route there is unchanged (within 7e-8).

### Item 4: parameters outside the designs

`R/se-reparam.R`, `R/se-check.R`. Two parts.

- `se_reparam_defer()`: the route at the top of cov_from_hessian() is
  refused when a parameter outside the designs' blocks has a row that
  tier 3 empties, or is a log sd whose curvature is below
  `se_edge_screen`. Degree 5 seed 4 now loses theta_1 as "short", as
  0.70.0 does (before: SE(log sd) 4.75e5, silent).
- `se_reparam_cov(fit, H, p, free)`: tier 3 removes such a parameter
  ("boundary", "nonfinite", an empty row, a bound) and then hands the
  rest to the route, conditional on it. On the prior build tier 3 in
  natural coordinates called x3 to x6 "flat" on the degree-8 fits with
  the sd at zero and kept the other coefficients' SEs at 0.05 to 0.10
  of lm()'s (`dev/illse-fx-edge8.R lane 1 2 4 6` on install 12, printed
  only): finite and wrong. Now, `dev/illse-fx-edge8.R lane 1:20`
  (`fx/edge8-lane.txt`): on the 16 seeds with the sd at zero only
  theta_1 is lost, and every coefficient SE is within 2.5e-5 of lm()'s
  (ML scale, the covariance with the sd at its boundary).

`dev/illse-rev-edge4.R` (`fx/rev-edge4d-lane.txt`, `-ob32.txt`), 80
fits, degrees 5 to 8, seeds 1 to 20: theta_1's verdict equals 0.70.0's
on 78. The two that differ: degree 5 seed 13 (0.70.0: SE NaN with a
large-gradient warning; lane: boundary, at the maximum) and degree 8
seed 4 (0.70.0: short, with a large-gradient warning; lane: boundary
at the maximum). No fit has a finite SE for a variance parameter that
0.70.0 lost. The coefficients that 0.70.0 lost on degree 8 now come
from the route: lm()'s on the 16 fits with the sd at zero (above), and
the route's accuracy of section 3 on the 4 with the sd away from zero.

### Item 5: the false large-gradient warning, and BLAS-decided verdicts

`R/fit.R`. `headroom_flat_rows()`: a row that tier 3 empties, or one
that is not finite (taken out one at a time, the worst first), is left
out of conv_headroom(), grad_headroom() and the Newton steps when a step
along it can gain at most `grad_tol`; `headroom_reparam_hessian()` then
takes the rest in the design's coordinates. Degree 6 seeds 4 and 7
(`fx/rev-edge5-lane.txt`): no warning about convergence, coefficient
SEs finite; theta_1 short and boundary, as 0.70.0 names them. Test:
"a fit at the maximum with an sd at zero is not unconverged".

Verdicts that differ between the reference BLAS and OpenBLAS 0.3.32 on
the 80 fits (theta_1's verdict, the lost count and the warning tags):

    build                          differing fits
    0.70.0                         2 (deg 5 s13, deg 7 s11)
    lane, before the fix           5 (+ deg 6 s4, deg 7 s14, deg 8 s5)
    lane, install 13               1 (deg 7 s11)

The three the lane introduced: degree 6 seed 4 was the false gradient
warning (item 5); degree 7 seed 14 took the route under one BLAS and
the boundary verdict under the other (the defer rule sends both to
tier 3); degree 8 seed 5 had a theta row that was empty under one BLAS
only, so the Newton steps ran under one and not the other (the steps
now leave that row out). Degree 7 seed 11 is nlminb's own: under
OpenBLAS it stops with a nonzero code and a convergence warning, under
the reference BLAS at code 0, with theta_1 nonfinite under both; 0.70.0
differs on it too.

### Item 7: cost

`R/se-reparam.R`, `R/se-check.R`, `R/fit.R`. Kept on the objective's
environment, keyed on the point: the Hessian in the design's
coordinates (`reparam_hessian_refined()`), the AD Hessian
(`obj_he_cached()`) and the finite-difference outer Hessian
(`fit_outer_hessian()`, which conv_headroom() and grad_headroom() now
read at the fit's point instead of their own optimHess()). The Newton
steps in the design's coordinates reuse a Hessian while its quadratic
model holds (a full, undamped step whose gain is within half of its
prediction), with at most five Hessians and ten steps. The fits end a
little further from the maximum than with five fresh Hessians (a
Newton gain of 1e-7 to 8e-7 instead of 1e-9 to 1e-13 on the degree-8
fits of `rev-lmm8-lane.txt`); every accuracy number in this section is
at that point.

`dev/illse-rev-cost.R` (`fx/rev-cost-*.txt`), n = 1e5, 100 outer
parameters, the four runs at once:

    design      0.70.0    lane before fix     lane, install 13
    no RE       442 s     1344 s (2.5x)*      641 s (1.45x)
    (1 | g)     863 s     2381 s (2.1x)*      1309 s (1.52x)

\* the reviewer's run, on a differently loaded machine; the ratios are
the comparable numbers. The Hessians in the design's coordinates now
cost 188 s and 354 s in all (two computations each), from 672 s and
1006 s. The rest of the difference is 11 s and 92 s. 0.70.0 stops on
both fits with a large-gradient warning (short of the maximum); the
lane reaches it. At n = 1e4
(`dev/illse-fx-cost.R`, `dev/illse-fx-cost2.R`, printed): 52.9 s
against 0.70.0's 32.8 s, and 119 s against 74.8 s with `(1 | g)`.

### Item 8: merge preparation

- `test-se-conv.R`'s short fit is now a poisson `y ~ x` from the
  explicit start `beta = c(-3, 3)` with one iteration, which no start of
  frm()'s own reaches; it still pins a coded stop short of the maximum
  that keeps finite, unjudged SEs and says so.
- The ci701 conflict in cov_from_hessian(): the route at the top runs
  only when `!any(held) && !any(edge)`, so it depends on `edge`. With
  ci701's `edge <- se_edge_sd(fit, H) | ec`, a correlation at the clamp
  refuses the top route; tier 3 then removes the edge parameters and
  hands the rest to the route through `free` (item 4), which needs only
  that they lie outside the designs' blocks. The composition holds.

### Batteries after the fix (install 13)

- Seflags batteries, reference BLAS and OpenBLAS 0.3.32
  (`dev/illse-log/sf-g2/`): bitwise 30 of 30 `identical()` to 0.70.0
  under both. False-alarm battery: 820 fits per BLAS, 0 differences in
  code, lost set, SE warnings or other warnings; the largest
  |SE / ref - 1| differs on 80 rows each (79 and 78 raw polynomials of
  degree 4 to 6, plus `negbin_large_shape` seed 8), lane lower in every
  model (degree 6 `lm`: 0.0177 against 0.0249 and 0.0370). The lane's
  other battery outputs moved from `sf-f3` by at most 1e-4 relative
  (the Newton stops above).
- The 40 files of `dev/illse-log/ob-subset.jobs`, reference BLAS and
  OpenBLAS 0.3.32 (scans `g2-ref`, `g2-ob32`): 40 of 40, 4026 pass, 0
  fail, error, skip or warning under each.
- Seen to fail on 0.70.0 (`dev/illse-fx-tbytest.R`,
  `fx/t-se-reparam-base-bytest.txt`, `fx/t-se-reparam-base.txt`):
  `test-se-reparam.R` 16 pass, 40 fail, 1 error. The new tests fail
  there: vcov_cluster 6, insight 2, marginaleffects 3, the coefficients
  with the sd at its boundary 3, the sd at zero at the maximum 1 (its
  coefficient SEs). The test of the sd near zero passes on 0.70.0: it
  pins 0.70.0's verdict, which the build before this fix broke
  (`dev/illse-rev-log/edge4-lane.txt`, `edge5-lane.txt`).
  `test-se-conv.R` 14 pass, 3 fail (`fx/t-se-conv-base.txt`).
- All eight packages, reference BLAS, one file per process (scan
  `g2-full`, `dev/illse-log/scan/g2-full.sum`): 374 of 374 files, 26169
  pass, 0 fail, error, skip or warning (frmtmb 231 files and 19022,
  coupling 11 and 548, eam 29 and 1765, latent 10 and 366, learn 15 and
  515, ode 11 and 551, sample 51 and 2806, spline 16 and 596).
- R CMD check --as-cran (`dev/illse-check/`): frmtmb Status: 1 NOTE
  (HTML manual, V8 unavailable; tests 16 min, OK); frmtmb.spline
  Status: 1 NOTE (the same NOTE; tests OK).

## 11. Fix round 2 (review round 2, R2-1 and R2-2)

Install 14 (`dev/illse-log/install-14.txt`). No Rd changed since the R
CMD check of section 10 (every `man/*.Rd` is byte-identical to the
checked tarball's), so the check was not rerun.

R2-1, `R/fit.R`. headroom_flat_rows() no longer decides a verdict; it
only says which rows the Newton steps leave alone. The new
headroom_drop_rows() leaves a row out of conv_headroom() and
grad_headroom() only when it is a log sd, and it asks the likelihood
what larger values of that sd gain, the rest held (se_sd_gain_up(), the
SE check's upward probe). That gain is added to the headroom of the
other parameters before the comparison with `grad_tol`. Any other empty
or non-finite row stays in, and its unusable curvature says so. The
threshold is `grad_tol`, not se_sd_at_zero()'s `se_edge_tol`: degree 6
seed 4 is "short" by the SE check's 1e-6 rule, but a larger sd gains
less than `grad_tol` there, and at `se_edge_tol` the false gradient
warning of item 5 came back.

- The reviewer's real case (`dev/illse-fx-r2flat.R lane 1 4`):
  `y ~ x + (1 | g)`, true sd 0.5, start log sd -15, `iter.max = 3`.
  Seed 1 stops 91.6 short and seed 4 139.0 short. The warning now says
  "would still gain 54.2" and "71.8" (the probe is a lower bound),
  "short of the maximum", and the large-gradient warning "may not have
  converged" is back. Install 13 said "at the maximum". Seeds 2 and 3:
  "would still gain 51.6" and "69.4".
- New test in `test-se-conv.R` on seed 1: on install 13, 19 pass and 3
  fail (`fx/t-se-conv-i13.txt`); on install 14 all 22 pass.
- Degree 6 seeds 4 and 7 (`dev/illse-fx-edge.R lane 6 4 7`): no
  convergence warning; theta_1 short and boundary; coefficient SEs
  finite.
- The 80 sd-at-zero fits (`fx/r2-edge4-lane.txt`, `-ob32.txt`): one
  verdict differs between the BLAS (degree 7 seed 11, nlminb's code, as
  before). Against install 13 one fit changed, the same under both
  BLAS: degree 7 seed 4. A larger sd gains 1.1e-3 there, above
  `grad_tol`, so the fit now warns that it may not have converged, and
  as a stop short of the maximum it keeps sdreport()'s standard errors,
  which are NaN. 0.70.0 named theta_1 "short" with finite coefficient
  SEs. So theta_1's verdict equals 0.70.0's on 77 of the 80 fits. The
  16 degree-8 fits with the sd at zero still lose only theta_1.
- NEWS now says what the rule measures.

R2-2. Core NEWS: 106 sentences, none over 25 words (longest 25);
frmtmb.spline NEWS: 12 sentences, longest 21 (`dev/illse-rev2-news.R`).
The cost sentence now gives both timing runs (1.45 to 1.66, and 1.47 to
1.52 with `(1 | g)`).

Tests on install 14, one file per process (scans `r2-ref`, `r2-ob32`),
the same counts under the reference BLAS and OpenBLAS 0.3.32, 0 fail,
error, skip or warning: `test-se-reparam.R` 57, `test-se-conv.R` 22,
`test-grad-verdict.R` 74, `test-se-check.R` 127, `test-se-flags.R` 68,
`test-optrap.R` 147.
