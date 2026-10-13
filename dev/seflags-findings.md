# Lane seflags: findings

Lane `seflags` of the round opened 2026-10-07. Worktree
`frmtmb-wt-seflags`, branch `wt-seflags`, base 48c1e6c9 (frmtmb
0.69.0). Lane build: `C:/Users/adf44/source/r/wt-seflags-lib` over
`rellib-r7`. Base build: `rellib-r7`. OpenBLAS: `bash
dev/ciharden-openblas.sh 0.3.26 lapack` (Rblas 74 symbols and Rlapack
616 symbols to OpenBLAS, 4 threads), in
`dev/ciharden-out/Rob0.3.26-lapack`. Logs are local, in
`dev/seflags-log/` and `dev/ciharden-log/seflags-*`; every number
below was copied from them.

## Punch round 2c (scope-down, 2026-10-13)

The check of punch round 2b ("Check of punch round 2b" in
`dev/reviews/2026-10-08-seflags.md`) found B6 and B7: every consumer
that reads the covariance whole (emmeans' design route,
marginaleffects through `get_vcov()`, a caller's own `vcov()` or
`frm_joint_cov()`, `vcov_cluster()`, the influence measures) and the
ordinal category route got 0.3 to 255929 times the true SE in silence
where 0.69.0 gave NaN. The coordinator scoped the lane down: an
"illcond" direction gets no SE again, through the lost-direction
machinery every consumer already handles, with the improved message.

Build: install 29. Every run below started after it.

### What changed

- `R/se-check.R` is punch round 2's file (install 22) with these
  changes:
  - `se_tier3()`: an "illcond" parameter (a parameter whose variance a
    confirmed direction below machine epsilon of the largest carries
    more than `se_ill_share`, 1 percent, of) is lost, and its unit row
    joins `null`, as a flat, non-finite or concave row does. So a
    prediction or hypothesis that reads it at all gets NaN with the
    prediction warning, and one that does not keeps its SE. The
    direction alone is not in `null`: inside the data's range a row's
    cosine with it is 1e-10, below `se_pred_tol`, though it carries
    most of the row's variance (B4). The `flat$ill` record is gone.
  - `jc_nonest()`: the cosine and subspace tests only. B4's share rule
    is gone; with the parameters lost it is not needed, and it took
    the SEs of terms the data determine (B5).
  - The "illcond" reason (`se_lost_clauses()`): "the design is
    ill-conditioned: the likelihood curves along them by less than
    machine precision of its largest curvature, so one double-precision
    covariance cannot hold their standard errors beside the others.
    Centering the covariate, or poly() for its powers, removes this",
    never "the data do not determine". When a Newton step gains more
    than `se_ill_gain` (0.005), `ill_short_clause()` adds that the
    optimizer stopped short of the maximum, with the gain. The gain,
    `g' V g / 2` over the kept covariance and the direction, is
    computed in `se_tier3()` for the message only (in
    `fit$cache$se_ill_gain`); no consumer reads it.
- `R/predict.R`, `R/brms-shapes.R`, `R/confint.R`: punch round 2's
  files, without the `flat$ill` lines. `R/conditional-smooths.R`,
  `R/interop.R`, `R/methods-fit.R`, `R/extra-deriv.R`,
  `R/predict-brms.R`, `man/frm_lp_basis.Rd` and
  `man/frm_extra_cov_deriv.Rd`: back to 48c1e6c9. Gone: `ill = list(D,
  lam)`, `ill_var()` and its siblings, the per-route additions, the
  `frmtmb_se_illcond` warning, the draw change and the BREAKING note.
- Kept as reviewed: B1 (`se_design_null()`), B2 (no joining, the span
  re-measured by second differences), B3 (the sd probe at + 10), the
  curved-ridge rule and its cost screen, the unit-diagonal metric
  (m5), the re-measured covariance when nothing is lost (the second
  defect), the ci691 wording, the frmtmb.eam test edit.
- `NEWS.md`: the raw-polynomial bullet is degree 7 only and not
  BREAKING; a new bullet says how an ill-conditioned design is named.
- `dev/test-backlog.md`: "Accurate SEs on ill-conditioned designs"
  under "Open - high priority (lane seflags)", with the 2b measurements
  and B6 and B7 as the consumers it must cover; the QR-route and
  stopping-short entries updated, the "covariance read whole" entry
  folded into it.

### Tests

`test-se-flags.R`, 11 tests:

- "a direction below double precision is ill-conditioned, not flat"
  (seed 3 with (1 | g)): the nine coefficients lost as "illcond", one
  lost warning that says "the design is ill-conditioned", names
  `poly()` and the stopping short, never "do not determine"; the four
  `fitted(newdata)` rows NaN with the prediction warning.
- "terms beside an ill-conditioned polynomial keep their SEs" (B5):
  z and wb are not lost; seed 3 `hypothesis(b_wb = 0)` and seed 10
  `fitted(dpar = "sigma")` finite with no warning, within 2 percent of
  the `poly(x, 8)` fit, or 3 where the message says the fit stopped
  short; seed 4 asks both on every build. Under OpenBLAS the lane stops
  seeds 3 and 10 at false convergence and the test asks only that the
  fit said so.

Final build: 68 of 68 (reference BLAS), 62 of 62 (OpenBLAS), 0
warnings. Install 22 (overlay): 59 pass, 9 fail (the new wording and
every B5 assertion).

### The review's scripts, final build

`dev/seflags-p2c-run.sh` runs `dev/seflags-rev4-{illpred,illid,whole,
cum,routes,gain}.R` on the lane, the base and the lane's merge with
optrap's current tree (`R/fit.R` 9e0afeba, `R/autoscale.R` 517e41ed),
both BLAS (`dev/seflags-log/p2c/rev4/`).

- `whole` (B6; `I(x^k)`, no random effects, 12 seeds): on every
  checked fit `fitted()`, `A vcov() A'`, marginaleffects
  `predictions()`, `slopes()`, `comparisons()`, `emmeans()` and
  `emtrends()` are NaN, `fitted()` and emmeans with a warning; the base
  gave the same NaN except `fitted()`, which kept two rows finite at
  0.66 to 0.89 of the reference. On the merge, seeds 6 and 9 (one
  BLAS each) stop with the large-gradient warning, so the check does
  not run and the rows are 0.32 to 0.92 with no warning at the call,
  as on op4 in the review (the filed convergence-warning policy).
- `cum` (B7, cumulative seeds 2, 5, 9): every `eta` and `P(y = 1)` SE
  NaN with the prediction warning on lane and merge, both BLAS. The
  base kept two rows finite at 0.86.
- `illpred` (B4): every row of every checked degree-8 fit NaN with the
  warning. The one finite set is lane degree 8 with (1 | g) seed 4, code
  0 with the large-gradient warning (filed), unchanged since 0.69.0.
- `illid` (B5), SE over the `poly()` fit's: lane seed 3
  `hypothesis(b_wb = 0)` 0.988, `b_z` 0.996, sigma 1.000; seed 10
  `b_wb` 1.003, sigma 1.004; seeds 4 and 5 without random effects
  1.001 to 1.018. Merge: 0.995 to 1.023 (sigma at 1.023 where the fit
  stops 1.5 short). Without the ill-conditioned direction these SEs
  miss its share of their variance (2.5 percent of b_wb's on seed 3,
  hence 1.2 percent low). The grid rows over the polynomial are NaN
  with the warning on every checked fit; the base kept 24 of 27 finite
  on seeds 4 and 5.
- `routes`: on every checked fit each route is NaN with a warning
  (`fitted()`, `frm_linpred()`, emmeans and its contrast,
  marginaleffects, `conditional_effects()`, the frmtmb.spline curves)
  or finite on terms that do not read the lost coefficients:
  `VarCorr()`, `ranef()`, `hypothesis(b_wb = 0)` and `vcov()` of z and
  wb at 0.988 to 1.023. The `frm_lp_basis()` lines form `A V A'` from
  the seam's `V` without its `se_nonest` flags, as on the base. The
  fits not checked carry a convergence warning (lane OpenBLAS B5 seed
  3 at code 1; merge seed 4 without random effects with the
  large-gradient warning, rows 0.28 to 0.92), as on op4 in the review.
- `gain`: unchanged from 2b (the gain is now in the lost message).

### Reruns

- False-alarm battery (`dev/seflags-p2c-falarm.sh`, the review's
  `dev/seflags-rev3-falarm.R`, 41 models, 1640 fits, both BLAS, against
  the review's base run): 0 differences in the lost set, the SE
  warnings or any other warning (`p2c/falarm-sum.txt`).
- Suites, one R process per file, every gate on (`seflags-ref9`,
  `seflags-ob9`; 318 of 318 logs name `wt-seflags-lib`, the OpenBLAS
  ones `Rob0.3.26-lapack`):

        package           files  pass ref  pass OB  fail ref  fail OB
        frmtmb              224     17733    17727         0        1
        frmtmb.eam           29      1765     1765         0        0
        frmtmb.sample        50      2693     2693         0        0
        frmtmb.spline        15       592      592         0        0

  0 errors, skips or warnings. The OpenBLAS failure is `test-fuzz.R:59`
  (seed 20379118, Gamma REML), as on every build; `test-se-flags.R`
  asks 6 fewer there (above). frmtmb.eam's `test-defects.R` passes with
  its `require =`.
- Merge with optrap's current tree: the review's scripts above; the
  merge build is `dev/seflags-merge-out/lib` (md5s in
  `p2c/merge-optrap-md5.txt`).
- R CMD check --as-cran of frmtmb after the last edit
  (`dev/seflags-check.ps1`): `Status: 1 NOTE` (the HTML manual's math
  check skipped for want of V8, as before); tests 13 minutes, OK.

Every process this round started has ended.

## Punch round 2b (2026-10-13)

Superseded by punch round 2c (above) where it kept the ill-conditioned
direction beside the covariance; kept as the record of the measurements.

The final check after punch round 2 ("Final check after punch round 2"
in `dev/reviews/2026-10-08-seflags.md`) found B5: the 1 percent
variance-share rule of B4 took the SEs of quantities the data
determine (a covariate z, a factor w, sigma) beside a degree-8 raw
polynomial. The review measured that the re-measured curvature is
accurate and asked for the "illcond" directions to go back into the
covariance at that curvature, for parameters and predictions alike.

"Install 22" below is the lane before this round. For the tests it is
reconstructed by `dev/seflags-prepunch.R` with `SEFLAGS_OVERLAY` set to
`dev/seflags-out/install22-overlay.R`: `R/se-check.R`, `R/predict.R`,
`R/brms-shapes.R` and `R/confint.R` as install 22 had them, evaluated
over the final lane library. Final build: install 27. Every run below
started after it, except where a run is named as earlier. Install 28
changed one comment of `se_tier3()`; `test-se-flags.R` (66 and 62)
and the R CMD check ran again on it.

### What changed

1. **An "illcond" direction keeps its SE.** A curvature below machine
   epsilon of the largest that the objective confirmed is kept, at the
   curvature the objective measured. No parameter loses its SE for it
   and no prediction is refused; reason "illcond" is gone from
   `se_lost`. The parameters whose variance such a direction carries
   more than `se_ill_share` (1 percent) of are named in a warning
   (class `frmtmb_se_illcond`, also under `check_se = "stop"`, since
   the SEs are there): the design is ill-conditioned, the SEs use the
   curvature the log-likelihood itself gives, and centering or `poly()`
   removes it.
2. **Trust criterion.** A direction below machine epsilon counts as
   measured only when `se_span_remeasure()` confirmed it: the symmetric
   second difference of the objective is quadratic over a factor of two
   in step (`se_quad_lo` to `se_quad_hi`), and more than `se_span_keep`
   (half) of the curvature is left once the mixing with every other
   direction (from the Hessian times the direction, by gradient
   differences) is taken out. A direction that fails is not kept: it
   is "flat", with NaN, as before, and a prediction on it gets NaN by
   the cosine and subspace tests. A single-parameter prediction reads
   the same covariance as `fixef()`, so it gets the parameter's verdict.
   On the 73 confirmed directions of the scan below on install 22
   (lane and merge, reference BLAS; `dev/seflags-log/p2b/scan1-*.txt`),
   the curvature at 0.25, 0.5 and 2 times the step is within 0.03
   percent of the step's own.
3. **The direction is kept beside the covariance, not in it.** Folded
   into the one matrix (install 23; its `illpred` log was overwritten
   by the final run, and the numbers are copied from it): the
   coefficient SEs were right, but `fitted(newdata)` in the data's
   range was 1.16 and 1.47 times the exact SE on seed 8 (x = 1.5, 2)
   and 0 on seed 4 without (1 | g). A matrix in double precision cannot
   hold a variance 1e17 times another: `g' V g` loses every digit to
   cancellation. The direction is now a separate record (`ill =
   list(D, lam)`, natural units), and a variance adds `(g . d)^2 /
   lambda` (`ill_var()`):
   - `vcov()`, `summary()`, `fixef()` and `frm_joint_cov()` show the
     sum, `V + D diag(1 / lam) D'`.
   - `get_joint_cov()` carries `ill` (moved to the group effects
     through `J`); `hyp_par_cov()` carries it over its positions.
   - Consumers: `frm_lp_basis()` puts it in `extra_var` and
     `extra_cov`, where every consumer of the seam already adds them
     (`frm_linpred()`, `conditional_effects()`, frmtmb.spline);
     `frm_extra_cov_deriv()` gives its derivative from the differenced
     design rows; `predict_mean_se()`, the categorical SEs, the
     emmeans bridge, `conditional_smooths()`, `ranef(condVar = TRUE)`,
     `hyp_prop_var()` (`hypothesis()`, `VarCorr()`, multiple
     imputation), `brms_fixef_extra_vcov()`; `fit_fd_se()` takes
     `J d` as one central difference along d (`se_ill_step`, a
     thousandth of the direction's SE), since a column-wise J errs by
     more than `J d` is in the data's range; `predict()`'s draws add
     the direction as its own normal draw.
4. **The warning says when the fit stopped short.** The gradient test
   cannot see a shortfall where the curvature is 1e-17. On these fits
   a Newton step (`g' V g / 2` with the kept covariance and the
   direction) gains up to 4.5 of log-likelihood, and the warning says
   so above `se_ill_gain` (0.005, a tenth of an SE). On the merge with
   lane optrap's units the B5 fits stop 1.1 to 1.5 short, mostly along
   the weak kept directions beside the ill-conditioned one, and the
   Newton model gives the shortfall to 3 digits
   (`dev/seflags-p2b-gain.R`):

        fit     build  short   Newton (kept + ill)   after 3 steps
        seed 3  lane   0.0464  0.0000 + 0.0464       0.0000
        seed 10 lane   0.5906  0.0005 + 0.5965       0.0000
        seed 3  merge  1.4808  1.4820 + 0.0340       0.0000
        seed 4  merge  1.0808  1.0973 + 0.0030       0.0000

   Seed 7 without (1 | g): three Newton steps along the direction gain
   3.63, and the coefficient SEs go from 1.036 to 1.000 times exact
   (`dev/seflags-p2b-span2.R`). Polishing the estimates is filed.
5. **Draws.** Without the direction the kept covariance is singular,
   and `chol()` refused it under OpenBLAS, so `predict()` fell back to
   the plug-in draws with its warning (seed 8 with (1 | g), install
   25; that log was overwritten by the final run). An ill-conditioned
   fit now draws from the symmetric square root.

`jc_nonest()` is back to the cosine and subspace tests; the B4 share
rule and `flat$ill` are gone. One defect found on the way: `sdr$se_ill`
read as `$` partial-matched `sdr$se_illcond` when the first was NULL;
the names no longer share a prefix (`se_cond_named`).

### B4 and B5, final build

B4 rows (`dev/seflags-rev3-illpred.R`, fitted(newdata, re_formula =
NA) at x = 1, 1.5, 2, 2.5 over the exact SE from a QR of the whitened
design at the fit's sigma and sd), every checked degree-8 fit, lane
and merge, both BLAS (`run-{lane,m9}/illpred-{ref,ob}-*.txt`):

    fit              lane ref                 lane OB
    seed 3, (1|g)    1.001 1.001 1.002 1.006  1.001 1.001 1.010 1.006
    seed 8, (1|g)    1.002 1.001 1.002 1.016  1.002 1.001 1.002 1.016
    seed 10, (1|g)   1.002 1.000 1.001 1.004  1.004 1.000 1.000 1.003
    seed 4, no RE    1.000 1.000 1.000 1.000  1.000 1.000 1.000 1.000
    seed 4, (1|g)    not checked (code 0, large-gradient warning, filed)

    fit              merge ref                merge OB
    seed 3, (1|g)    1.003 1.004 1.011 1.020  1.000 1.004 1.014 1.019
    seed 8, (1|g)    1.002 1.001 1.001 1.015  1.002 1.001 1.001 1.015
    seed 10, (1|g)   1.005 1.001 0.999 1.003  1.003 1.000 1.001 1.004
    seed 4, no RE    1.000 1.000 1.000 1.000  1.000 1.000 1.000 1.000
    seed 4, (1|g)    1.001 1.000 0.999 1.000  0.999 1.000 1.001 1.000

Stated tolerance: 3 percent (all 0.999 to 1.020). No warning on any
row. Degree 9 drops x9 from the frame (the review's erratum), so its
rows are not compared.

B5 (`dev/seflags-p2b-illid.R`, the review's `illid` with its share
instrumentation only where the rule still exists), SE over the
`poly(x, 8)` fit's:

    fit             query                  lane ref  merge ref  merge OB
    seed 3, (1|g)   hypothesis(b_wb = 0)   1.0007    1.0076     1.0075
    seed 3, (1|g)   fitted(dpar = sigma)   1.000     1.023      1.023
    seed 10, (1|g)  hypothesis(b_wb = 0)   1.0037    1.0073     1.0073
    seed 10, (1|g)  fitted(dpar = sigma)   1.011     1.023      1.023
    seed 4, no RE   hypothesis(b_wb = 0)   1.0029    not checked  0.9976
    seed 5, no RE   fitted(dpar = sigma)   1.028     not checked  not checked

No warning on any query. The merge's 2.3 percent on sigma and lane
seed 5's 2.8 percent are the stopping point: those fits stop 1.1 to
1.5 short of the maximum, and lane seed 5 by 1.9 (sigma 0.38 against
0.37, as the review found); the warning says so. Under OpenBLAS the
lane alone stops seeds 3, 8 and 10 with (1 | g) at false convergence (code 1, as
on install 22), and the merge stops seeds 4 and 5 without random
effects with the large-gradient warning; the check does not run on
those (filed policy).

Every route on two fits (`dev/seflags-p2b-paths.R`, seed 8 with
(1 | g) and seed 4 without, x = 1, 1.5, 2, 2.5; lane and merge, both
BLAS): `fitted()`, `frm_linpred(se.fit = TRUE)`, `frm_lp_basis()`
(`A V A' + extra_var` and `+ extra_cov`), `hypothesis()` of the
polynomial at x = 1.5, and `fit_fd_se()` are 1.000 to 1.016 times
exact, with no warning; `frm_joint_cov()`'s diagonal is the coefficient
SEs (1.018 on seed 8, 1.000 on seed 4); `predict()`'s parameter draws
(2000, Monte Carlo 1.6 percent) 0.990 to 1.036.

Coefficient SEs over the exact ones (`dev/seflags-p2b-scan.R`: 12
seeds, degree 8 and the B5 design, with and without (1 | g)):

    build        checked  SE / exact      not checked
    lane ref     25       0.9999 to 1.040 23 (19 false conv., 4 large grad.)
    lane OB      26       1.0000 to 1.040 22 (19, 3)
    merge ref    40       0.9837 to 1.049 8 (large gradient)
    merge OB     38       0.9989 to 1.050 10 (large gradient)

The excess over 1 is where the optimizer stopped short (item 4).

### Tests

`test-se-flags.R` (11 tests, 66 assertions under the reference BLAS;
62 under OpenBLAS, where seeds 3 and 10 of the B5 test stop at false
convergence and the test asks only that the fit said so):

- "an ill-conditioned design keeps its SEs and says why" replaces the
  degree-8 test: seed 8 with (1 | g); nothing lost; one warning, "is
  ill-conditioned", naming `poly()` and the stopping short, never the
  lost-SE phrase; coefficient SEs and the four `fitted()` rows within 3
  percent of exact, with no warning; `predict()`'s draws within 10
  percent (six Monte Carlo SEs).
- "terms beside an ill-conditioned polynomial keep their SEs" (B5):
  seed 3 `hypothesis(b_wb = 0)` and seed 10 `fitted(dpar = "sigma")`
  within 2 percent of the `poly(x, 8)` fit, or 3 where the warning says
  the fit stopped short (the merge); seed 4, which converges on every
  build, asks both.

Seen to fail on install 22: reference BLAS 51 pass, 15 fail (every new
assertion of both tests that pins a defect: SEs and rows NaN, the lost
warning, the prediction warning, the draws at the estimates);
OpenBLAS 53 pass, 9 fail. Final build: 66 of 66 (reference), 62 of 62
(OpenBLAS), 0 warnings.

### Reruns on the final build

- Merge with optrap's current tree (`R/fit.R` 9e0afeba, `R/autoscale.R`
  517e41ed, unchanged since op3; `dev/seflags-log/p2b/merge-optrap-md5.txt`),
  both BLAS, the merge tree's own test files: test-se-flags 66,
  test-se-check 127, test-se-tier 91, test-optrap 131, test-optima
  100, test-diagnostics-ux 130, test-autoscale 46,
  test-se-lost-re-predict 20, test-brms-postfit 76, frmtmb.eam
  test-defects 60. All pass, 0 fail, error, skip or warning.
- The review's `illpred`, `illid`, `illfull`, `deg7pred` and
  `deg7exact`, lane and merge, both BLAS (`p2b/run-{lane,m9}/`).
  Degree 7: `deg7exact` identical to install 22; `deg7pred` rows 0.999
  to 1.004 times exact.
- False-alarm battery (the review's `dev/seflags-rev3-falarm.R`, 41
  models, 1640 fits, both BLAS, against the review's base run):
  0 differences in the lost set, the SE warnings or the other warnings
  (`p2b/falarm-sum.txt`). The lane's own battery and review scripts
  (`dev/seflags-revrun.sh`, `rev-p2b`, `rev-p2bob`, `falarm-p2b`,
  `falarm-p2bob`): identical to install 22's runs (`rev-ref7`,
  `rev-ob7`) but for glmer's seventh digit in `year-*.txt`.
- Attack set and studies (`dev/seflags-attacks.sh`, `att-p2b`,
  `att-p2bob`): lane against base as on install 22, except
  `setier-rev-fixes-table`, whose `a * b` and `a + exp(b)` rows change
  between runs of the base too (SE against NLFLAT on the same build).
- Scans of all eight suites, one R process per file, every gate on
  (`seflags-ref8`, `seflags-ob8`; 365 of 365 logs name
  `wt-seflags-lib`, and the OpenBLAS ones `Rob0.3.26-lapack`):

        package           files  pass   ref fail  OB fail  err skip warn
        frmtmb              224  17731         0        1    0    0    0
        frmtmb.coupling      11    548         0        0    0    0    0
        frmtmb.eam           29   1765         0        0    0    0    0
        frmtmb.latent        10    366         0        0    0    0    0
        frmtmb.learn         15    515         0        0    0    0    0
        frmtmb.ode           11    551         0        0    0    0    0
        frmtmb.sample        50   2693         0        0    0    0    0
        frmtmb.spline        15    592         0        0    0    0    0
        total               365  24761         0        1    0    0    0

  The OpenBLAS failure is `test-fuzz.R:59` (seed 20379118, Gamma REML,
  NA/NaN gradient), as on every earlier build. OpenBLAS passes 24757:
  the count differs where a test branches on the build, as
  `test-se-flags.R` does (above).
  frmtmb.eam's `test-defects.R` passes with the `require =` of punch
  round 2: its directions are flat, not ill-conditioned.
- R CMD check --as-cran of frmtmb (`dev/seflags-check.ps1`, install 28,
  after the last edit): `Status: 1 NOTE`, the HTML manual's math
  check skipped for want of V8, as before; tests 11 minutes, OK.

Every process this round started has ended.

### Files

Functions edited: in `R/se-check.R`, `se_tier3()` (the direction kept
beside the covariance, the record `ill`, the names and the Newton gain
for the warning), `sdr_rescue_hessian()`, `sdr_sep_lost()`,
`se_check()`, `se_report()`, `se_deferred_report()`,
`se_lost_clauses()` (no "illcond" reason) and `jc_nonest()` (no share
rule); new `se_illcond_message()`, `ill_var()`, `ill_cov()`,
`ill_diag()`, `ill_add()`, `ill_at()` and the constants `se_ill_gain`
(0.005) and `se_ill_step` (1e-3); `se_ill_share` now only names. In
`R/predict.R`, `get_joint_cov()`, `joint_cov_repair()`,
`predict_mean_se()`, `ord_prob_se()`, `frm_joint_cov()`,
`frm_lp_basis()` (and its documentation) and `lp_basis_out()`. In
`R/extra-deriv.R`, `frm_extra_cov_deriv()` (and its documentation) and
new `ill_deriv_cov()`. In `R/brms-shapes.R`, `fit_draw_space()`,
`fit_fd_se()` and `brms_fixef_extra_vcov()`. In `R/predict-brms.R`,
`predict_par_drawer()`. In `R/confint.R`, `hyp_par_cov()`,
`hyp_prop_var()` and `hypothesis()`'s `vcov =` branch. In
`R/conditional-smooths.R`, `conditional_smooths.frmtmb_fit()`; in
`R/interop.R`, `emm_basis_grid()`; in `R/methods-fit.R`,
`ranef.frmtmb_fit()`. Regenerated: `man/frm_lp_basis.Rd`,
`man/frm_extra_cov_deriv.Rd`. Also `tests/testthat/test-se-flags.R`,
`NEWS.md` (the raw-polynomial bullet, now BREAKING for
`frm_lp_basis()`'s `V` on such a fit), `dev/test-backlog.md`. New dev
scripts: `dev/seflags-p2b-{scan,span,span2,gain,b5sig,paths,illid,
illfull,msg}.R`, `dev/seflags-p2b-{run,batch,final,falarm,merge-tests}.sh`,
`dev/seflags-p2b-falarm-sum.R`; `dev/seflags-prepunch.R` leaves active
bindings alone; overlay `dev/seflags-out/install22-overlay.R`.

Filed (`dev/test-backlog.md`): the optimizer stopping short where the
curvature is this small (a Newton polish, lane optrap's ground); a
covariance read whole (`vcov()`, `frm_joint_cov()`) keeps the
cancellation that the lane's own routes avoid; the QR route updated; a
mixing weight at its limit named "boundary" (the review's suggestion).

## Punch round 2 (2026-10-12)

The re-check ("Re-check after punch round 1" in
`dev/reviews/2026-10-08-seflags.md`) found one blocker (B4) and three
minors. "Punch 1" below is the lane as it was before this round,
reconstructed for the tests the same way as before
(`dev/seflags-out/punch1-se-check.R`, today's file with this round's
changes to `se_tier3()`, `jc_nonest()` and `jc_drop_boundary()` taken
out, evaluated over the lane library by `dev/seflags-prepunch.R` with
`SEFLAGS_OVERLAY`).

Functions edited: in `R/se-check.R`, `se_tier3()` (records the
unit-diagonal scaling of every outer parameter and the "illcond"
directions with their curvatures and the kept covariance; a
re-measured span is kept even when nothing is lost), `jc_nonest()`
(the variance share for "illcond" directions; every cosine in the
unit-diagonal metric), `jc_drop_boundary()` (keeps the scaling when no
flat column is left), `cov_from_hessian()` (returns tier 3 when a span
was re-measured) and `se_lost_clauses()` (the "flat" wording); in
`R/predict.R`, `joint_cov_repair()` (the repaired covariance whenever
the check made one, and the new fields moved to the joint positions);
in `R/brms-shapes.R`, the guard in `fit_fd_se()`. No change to
`R/fit.R` or `R/confint.R`. Also touched: `tests/testthat/test-se-flags.R`
(two tests), frmtmb.eam's `tests/testthat/test-defects.R` (m5,
below), `NEWS.md`, `dev/test-backlog.md`; new dev scripts
`dev/seflags-eamdbg.R`, `dev/seflags-check-eam.ps1`; `dev/seflags-prepunch.R`
and `dev/seflags-merge-run.sh` extended.

### B4: predictions on "illcond" fits got narrow SEs in silence

An "illcond" direction leaves the covariance, and the prediction test
read it only by cosine against `se_pred_tol`; at a curvature below
machine epsilon of the largest, a cosine far below 1e-6 carries most
of a prediction's variance. The parameter rule's variance share now
applies to predictions: `jc_nonest()` gives a prediction no SE (and
the prediction warning) when an "illcond" direction would add more
than `se_ill_share` (1 percent) of the variance the kept directions
give, `(g . d)^2 / lambda` at the re-measured curvature against
`g' V g`, both in natural units, so the ratio does not depend on any
scaling. This is the conservative answer; the QR route for gaussian
linear predictors (m6) is not built and is filed.

`dev/seflags-rev2-illpred.R` (fitted(newdata) at x = 1, 1.5, 2, 2.5,
against the exact SE from a QR of the whitened design), reference
BLAS, log `merge-m9/illpred-{lane,merge}.txt`:
every row of the six "illcond" fits (degree 8 with (1 | g), seeds 3,
8, 10; without, seed 4; degree 9 seeds 3 and 4) is NaN with the
prediction warning on the lane and on the merge with lane optrap,
under both BLAS builds (`merge-m9`, `merge-m9ob`); punch 1 gave 0.81
to 0.98 of the exact SE in range and 0.31 to 0.35 at x = 2.5, silent.
The one finite fit, degree 8 with (1 | g), seed 4 on the lane, is not
"illcond": it carries the large-gradient warning, so the check does
not run (filed below). `dev/seflags-illscan.R` (40 fits): lane 16
"illcond" and 24 with nothing lost (code 1 or a convergence warning;
reference BLAS), 20 and 20 (OpenBLAS); on the merge all 40 "illcond"
under both builds.

While checking the absent case, a second defect showed: on a fit with
a re-measured span and nothing lost, `cov_from_hessian()` returned
sdreport()'s inverse when it existed, and `joint_cov_repair()` built
the prediction covariance only when something was lost, so
predictions read `solve(Q)`, the inverse of the unresolved Hessian.
Degree 7 with (1 | g), seed 7: fitted() at x = 2.5 was 127.1 against
the exact 131.9 on the lane, and 0 on the merge with lane optrap
(`merge-m8/flags-merge.txt`, the new test failing 2 expectations on
that build). Both now use the re-measured covariance: 132.3 against
131.9, and the in-range rows within 0.3 percent, on both arms.

Tests (`test-se-flags.R`, 51 assertions, 10 tests): the degree-8 test
now also requires fitted() at x = 1, 1.5, 2, 2.5 to be NaN with the
prediction warning; a new test, "a prediction on a resolved near-null
direction keeps its SE", requires degree 7's four rows within 1
percent of the exact SE. Lane 51 pass; punch 1 49 pass, 2 fail (the B4
expectations; the degree-7 rows are fixed in `R/predict.R`, which the
overlay does not replace, and were seen to fail on the merge build
before the fix); `rellib-r7` 27 pass, 24 fail (`p2c/`).

### m5: the prediction test in the units-free metric

`jc_nonest()` measured the `se_pred_tol` columns as cosines in the
optimizer's units, which lane optrap's `par_units` change. It now
reads them in the unit-diagonal metric of the outer Hessian, whose
scaling `se_tier3()` records (`dscale`, moved to the joint positions
by `joint_cov_repair()`); the optimizer's units remain only as a
fallback for an analysis that recorded none (the separation path).
On the merge with optrap's current tree, the degree-8 rows that were
NaN at x = 2.5 and finite in range on the re-check's op2 build are NaN
on all four rows, as on the lane.

### m7 and lane ci691's two findings

- The two m7 defects (pre-existing) are filed in `dev/test-backlog.md`.
- The "flat" wording now says "flat along them, or too nearly flat for
  its curvature to be measured, at the estimates", which is what the
  curved-ridge rule measures (ci691: a ridge the profile bends along
  slightly). Cheap; done.
- The SE check not running on a fit with a convergence warning is a
  policy (the convergence warning explains the fit); changing it is
  not cheap. Filed, with degree 8 with (1 | g), seed 4, as an example
  (code 0 with the large-gradient warning, fitted() 0.37 of the exact
  SE at x = 2.5, silent; the base does the same).

### Reruns on the final build

Code: install 22 (`R/se-check.R`, `R/predict.R`, `R/brms-shapes.R`);
every run below started after it. One later edit is to a test only:
frmtmb.eam's `test-defects.R` (below), rerun alone under both builds.

Scans (`seflags-ref7`, `seflags-ob7`; all eight suites, one R process
per file, every gate on; 365 of 365 logs name `wt-seflags-lib`, and
the OpenBLAS logs `Rob0.3.26-lapack`):

    package           files  pass   ref fail  OB fail  err skip warn
    frmtmb              224  17716         0        1    0    0    0
    frmtmb.coupling      11    548         0        0    0    0    0
    frmtmb.eam           29   1765         0        0    0    0    0
    frmtmb.latent        10    366         0        0    0    0    0
    frmtmb.learn         15    515         0        0    0    0    0
    frmtmb.ode           11    551         0        0    0    0    0
    frmtmb.sample        50   2693         0        0    0    0    0
    frmtmb.spline        15    592         0        0    0    0    0
    total               365  24746         0        1    0    0    0

The scans themselves ran frmtmb.eam's `test-defects.R` before its
edit and showed one escaped warning there under both builds (59 pass,
1 warning); the eam row above has the edited file's own runs
(`seflags-tgt3`, `seflags-tgt3ob`: 60 pass, 0 fail, 0 warnings). That
warning is m5 at work: `test-defects.R:92` fits a wiener plus
lognormal mixture whose mixing weight runs to 6e-11, so the Wiener
parameters are flat, and `residuals(type = "response")` now gets NaN
SEs with the prediction warning on all 300 rows. In the unit-diagonal
metric the rows' largest cosine with the lost directions is 8.3e-5;
in the optimizer's units it was 3.8e-12 (`dev/seflags-eamdbg.R`). Every
row does move along a direction the fit does not determine, so the
warning is the invariant verdict; the test (which asserts only the
estimates) now allows and requires it. The OpenBLAS failure is
`test-fuzz.R:59`, seed 20379118, lane optrap's restart, as before.

The attack set, lane against base (`bash dev/seflags-sum.sh att ref7`,
`ob7`): every differing line is items 3 and 4, the B2 and degree-8
changes, and the new "flat" wording (cases 16, pred 24, trapq 4 are
the wording; b1-rho 4 are SE digits; fixes-table 4 and 6: `a * b`, `a
+ exp(b)` named, and one identified row's optimizer code changing
between processes). Singular (600), small-sd (200) and separation
(500) studies: 0 verdict differences under both builds. The reviewer's
scripts (`rev-{ref7,ob7}`): join, sd2, sd, deg7 as in punch round 1;
false-alarm battery 41 models, 820 fits per build, 0 verdict
differences (`falarm-7-sum.txt`). Extras (`ext-{ref7,ob7}`): fixture 1
under OpenBLAS LAPACK `fitted(newdata)` 0 warnings (base 1); segap 9
of 9; the alias construction silent with SEs 0.2 to 0.212; the
flat-peak sd "flat". Cost (`cost-p1`): `se_check()` evaluations equal
to the base's on every design (0/8, 2/8, 0/0, 0/0 objective/gradient)
and 0 ridge calls.

The merge with lane optrap (optrap's `R/fit.R` md5 98274516,
`R/autoscale.R` 226a95c9, read only; `merge-m9`, `merge-m9ob`):
`test-se-flags.R` 51 of 51 on lane and merge under both builds;
degree 7 SE/lmer 1.09 all; degree 8 "illcond" with every prediction
NaN and warned; sd2 9 of 9; nl `a + b` 15 of 15; b1ref 10 of 10 off
the alias NaN with a warning, 2 of 2 on it finite.

R CMD check --as-cran after the last code edit (install 22):
frmtmb `Status: 1 NOTE` (the V8 math-rendering NOTE); frmtmb.eam
(after its test edit) `Status: 1 WARNING, 1 NOTE`: the NOTE is V8's,
and the WARNING is "'::' or ':::' import not declared from: 'rstan'"
in its tests, which is pre-existing (`test-sampling.R` calls
`rstan::` and rstan is not in Suggests since 0.69.0's 5f39a6f5); filed.

## Punch round 1 (2026-10-08)

The review (`dev/reviews/2026-10-08-seflags.md`) found three blockers
in the first round's new rules. All three are fixed, each with a test
that fails on the first round's code; two first-round claims are
withdrawn in place below (marked WITHDRAWN). Only `R/se-check.R`
changed in code. "Pre-punch" below is the first round's build,
reconstructed for the tests by evaluating its `se-check.R`
(`dev/seflags-out/prepunch-se-check.R`, the file install 9 had, which
the review compared with the installed namespace) over the lane
library (`dev/seflags-prepunch.R`).

Functions edited in `R/se-check.R`: `se_tier3()` (the joining loop
removed; the design alias, the span re-measurement, the "illcond"
rule, the screen of the ridge correction, the cap regime's prediction
rule), `se_tier3_could_act()` (the screen only with a nonlinear
parameter), `se_sd_gain_up0()` (the far probe), `se_curvature_real()`
(quiet probes), `se_ridge_curvature()` (quiet gradients),
`se_lost_clauses()` (the "illcond" reason), and the documentation of
`jc_nonest()`. New: `se_lp_outer_pos()`, `se_linear_pars()`,
`se_design_null()`, `se_quiet()`, `se_obj_curv()`,
`se_span_remeasure()`; constants `se_ridge_nonlin` (0.1),
`se_span_keep` (0.5), `se_span_max` (4), `se_design_flops` (1e10),
`se_ill_share` (0.01), `se_sd_far` (10). No change to `R/fit.R`,
`R/predict.R`, `R/confint.R` or `R/brms-shapes.R` in this round.

`test-se-flags.R` (47 assertions, 9 tests): lane 47 pass; pre-punch 29
pass, 18 fail (B1 4, B2 2, B3 4, degree 7 3, degree 8 4, the alias
test 1); `rellib-r7` 25 pass, 22 fail (logs `dev/seflags-log/p1d/`).

A fourth change came from lane optrap's review: on a trial merge with
optrap's parameter units, degree 8 with (1 | g) reaches code 0 and the
first punch build gave it SEs 5.7 to 6.4 times lmer()'s, silently
(below, "Degree 8 and parameter units").

### B1: a hypothesis off an exact alias passed

The cap `min(tau, 0.5 * largest projection)` is a rule for naming
parameters, not an error bound, and the first round used it as the
prediction threshold. Two changes:

1. **The design's alias is taken from the design.**
   `se_design_null()` (new) finds, by QR with the frame's own rank
   test, the directions along which a linear predictor's fixed design
   is exactly flat (a smooth's fixed column joins the design after the
   frame's alias check). The part of the flat subspace in their span
   is replaced by them: every coefficient they move is named (on
   fixture 1 the intercept too; its unit-diagonal share of the exact
   alias is 0.11), and they go into the null basis as exact columns,
   tested at `se_pred_tol`. The covariance is the inverse on the
   complement of everything removed, so it excludes the exact
   direction and not the measured one.
2. **Where tau does not bound the subspace's error, a prediction that
   reads a named parameter gets no SE**, as in 0.69.0; the flat
   subspace with its threshold goes to `jc_nonest()` only where tau
   is the threshold.

Fixture 1 under OpenBLAS 0.3.26 with LAPACK (`dev/seflags-rev-pred3.R`,
`rev-ob6/pred3-lane.txt`): `b_Age + t * bs_sAge_1` is NaN with a
warning at t = 0.15 and 0.10 (pre-punch 0.5603 and 0.4777, silent) and
at every other t of the sweep; at t = 0.190041455 (full precision,
`dev/seflags-out/pred3x.R`) it is 0.629, silent. At t = 0.190041 (six
digits, as the review's script prints it) it is NaN with a warning:
that hypothesis is 4.6e-7 off the alias, which the exact test reads as
off. `fitted(fit1, newdata)` stays silent
(`ext-ob6/dh7diag-lane.txt`).

Reference BLAS (`dev/seflags-b1ref.R`, logs
`p1a/b1ref-{pre,lane}.txt`): the alias of `y ~ x + s(x) + (1 | g)`
beside a pair z2 = z1 + N(0, 1e-4) had a noise bound of 0.021, and the
hypothesis b_x + t bs_sx_1, estimable only at t = b, got finite SEs and
no warning for t up to 3 percent off b (pre-punch 0.40608, 0.405898,
0.404083, 0.40005 at 1e-4, 1e-3, 1e-2, 3e-2 off). The lane gives NaN
with a warning at every t off b and 0.4061 at t = b. This is the
suite's B1 test; the review had not found a reference-BLAS
construction.

### B2: confirmed directions joined to removed ones

The joining loop is gone. A direction the step test confirmed is never
removed on the Hessian's word. What the joining was meant to do, it
now does by measurement: a removed direction that no objective
evaluation has asked about (an eigenvalue below 100 epsilons of the
largest, or within its noise of zero and not positive, and not in a
design alias) is measured again together with the confirmed small
directions, by second differences of the objective on their span
(`se_span_remeasure()` and `se_obj_curv()`, new; the review's
`dev/seflags-rev-deg7fix.R` measurement). A re-measured direction is
kept, at the measured curvature, when the objective confirms it: a
quadratic change of 4 `grad_tol` within a step of 4^20, and at least
half of that curvature left once its mixing with the other directions
is taken out (the Hessian times the direction, from central differences
of the gradient at the same step). The second condition was needed:
the nl `a + b` alias read a curvature of 1.95e-19 through its mixing
with `a + b`, and 1.7e-23 after the span was re-measured, through a
residual mixing of 6e-8 with the near-collinear pair; without it 3 of
15 such fits kept the alias's SEs in silence (`p1a/join*.txt`).

- Fixture 1 under OpenBLAS LAPACK (`dev/seflags-rev-fx1theta.R`,
  `rev-ob6/fx1theta-lane.txt`): theta_1 and theta_2 SEs 0.4086 and
  0.4145 (base 0.4086 and 0.4144); lost (Intercept), Age,
  `s(Age).fx1` (the exact alias), theta_3 boundary, theta_4 flat.
- `dev/seflags-rev-join.R`, 15 fits per variant
  (`rev-{ref6,ob6}/join-{lane,base}.txt`): the pair lost as "flat" on
  1 of 15 (reference; seed 2 at s = 1e-6, the base's too) and 0 of 15
  (OpenBLAS) on lane and base alike (pre-punch 10 and 11); the nl
  `a + b` alias named on 15 of 15, lane and base, both BLAS builds; the
  pair's SEs are the no-alias fit's within 0.05 percent up to s = 3e-6
  (seed 1 at 3e-6: 7604.19 against 7605.29).
- Degree 7 with (1 | g), seed 7 (`rev-ref6/deg7-lane.txt`,
  `dev/seflags-deg7diag.R`): nothing lost, SE/lmer 1.09 for all eight
  coefficients (3032, 14800, 30720, 35110, 23890, 9671, 2158, 204.8
  against 2774 ... 187.4); base SE/lmer 0.119, 0.089, 0.060, NaN, NaN,
  NaN, 0.054, 0.081. The 9 percent is the review's own figure for the
  same measurement and is not settled (the outer Hessian carries the
  theta and sigma coupling that lmer's conditional covariance does
  not).

### Degree 8 and parameter units

Lane optrap scales the identity-link mean coefficients and group sds
of gaussian-type families by `fit$par_units`. The thresholds this lane
added are all in the unit-diagonal metric of the natural-unit Hessian
(eigenvalues, projections, tau, the curvature probes), which a
diagonal change of units leaves alone; what the units do change is
where the optimizer stops and the finite-difference Hessian, which is
taken in the optimizer's units (as sdreport() takes it). The empty-row
floor of 0.68.0 (`se_empty_row`, `se_noise_mult`) is in the
optimizer's units on purpose, since that is where the finite-difference
noise is, and is unchanged.

On the merge, degree 8 with (1 | g), seed 7, stops at code 0 (on the
lane alone at code 1, where the check is not asked). Its scaled design
has a condition number of 7.3e8 (a squared singular-value ratio of
1.9e-18), so the smallest curvature is below machine epsilon of the
largest and no double-precision inverse resolves it; lmer()'s own SEs
move by 13 percent between its two optimizers there. The first punch
build kept those directions as confirmed and gave SEs 5.68 to 6.35
times lmer()'s (`merge-m1/deg8-merge.txt`). Now a confirmed curvature
below `.Machine$double.eps` times the largest is "illcond": a
parameter whose variance such a direction carries more than 1 percent
of (`se_ill_share`), with a loading above sqrt(`se_flat_tol`), loses
its SE with a reason that says the design is ill-conditioned (raw
powers; center or use poly()), never that the data do not determine
it. Degree 7 (smallest 5.4 epsilons of the largest) is not affected.
`dev/seflags-illscan.R` on the lane alone (seeds 1 to 10, degree 8):
with (1 | g) seeds 3, 8, 10 reach code 0 and lose the nine
coefficients as "illcond" (base: six "flat"); without (1 | g) seeds 3
to 7 (base: nine "flat"); the others stop at code 1 or carry the
convergence warning, as on the base (`p1b/illscan-{lane,base}.txt`).
Test: "a direction below double precision is ill-conditioned, not
flat" (seed 3 with (1 | g)).

Both arms, reference BLAS and OpenBLAS LAPACK
(`dev/seflags-merge-run.sh`, `dev/seflags-merge.sh`: this lane's
package with optrap's `R/fit.R`, `R/autoscale.R` and changed tests
copied in, optrap's tree read only; logs `merge-m7`, `merge-m7ob`):

    check                         lane ref  merge ref  lane OB   merge OB
    test-se-flags.R               47/47     47/47      47/47     47/47
    degree 7, SE/lmer             1.09 all  1.09 all   1.09 all  1.09 all
    degree 8, verdict             code 1    illcond 9  illcond 9 illcond 9
    sd2 boundary (of 9)           9         9          9         9
    join: pair lost (of 15)       1         1          0         0
    join: nl a+b named (of 15)    15        15         15        15
    b1ref: off-alias t NaN+warn   10 of 10  10 of 10   10 of 10  10 of 10

The pair lost on one reference fit (seed 2, s = 1e-6) is the base's
too. Fixture 1 (fx1theta, pred3) stops at code 1 on the merge under
both BLAS builds (optrap's path), so its B1 and B2 checks do not apply
there; on the lane under OpenBLAS LAPACK they hold as above. The merge
was built from optrap's tree as it stood at the run (`R/fit.R` md5
98274516, `R/autoscale.R` d7cf916b, `dev/seflags-log/merge-optrap-md5.txt`;
optrap's tree as of 2026-10-07 16:20, the latest when this was rerun on
2026-10-12; an earlier merge with its 15:3x tree, `merge-m6`, gave the
same table).

### B3: an sd read only far above its estimate

`se_sd_gain_up()` now also evaluates the objective at the log sd + 10
(`se_sd_far`) and counts that change in `moved`. On the review's
`dev/seflags-rev-sd2.R` (9 fits, `rev-ref6/sd2-{lane,base}.txt`) the
lane gives the boundary message and theta_1 "boundary" on 9 of 9, as
the base does (pre-punch: "flat" on 5, nothing on 4 with SEs 778 to
995). Fit 58 stays "flat" (its objective is unchanged at every
value). `dev/seflags-rev-sd.R` (63 fits per build): lane and base
verdicts are identical (`rev-{ref6,ob6}/sd-*.txt`; the gain and moved
columns differ only because the base has none).

### Cost

The review's screen is in: `se_ridge_curvature()` is asked only of a
direction whose unit-diagonal loadings on parameters outside the
family's primary linear predictors (`se_linear_pars()`, new: nl
parameters, simplex coordinates, variance and distributional
parameters) have norm above 0.1 (`se_ridge_nonlin`), and a model with
no such parameter screens at `se_flat_tol` as 0.69.0 did. The
design-null QR runs only when a direction was removed.
`dev/seflags-cost-p1.sh` (lane and base side by side, minimum of 3,
`dev/seflags-log/cost-p1/`):

    design               se_check alone     check evals fn/gr  ridge
                         base     lane      base    lane       calls
    glmm n=1e5           6.17 s   5.72 s    0/8     0/8        0
    glmm n=1e5, year     9.08 s   9.11 s    2/8     2/8        0
    lmm, 400 coef.       0.00 s   0.00 s    0/0     0/0        0
    glm, 400 coef.       0.50 s   0.49 s    0/0     0/0        0

The 2051-coefficient bernoulli fit: 19.0 s against 18.8 s with the
check (base, lane, the faster of two rounds), 0 ridge calls. On the
review's false-alarm battery (below) ridge calls fall from 20 to 60
per linear model to 0; they remain on the nonlinear models (emax 8,
logistic_early 22, lapse 8). These are the final build's (install 19,
`cost-p1` rerun at the end of the chain); install 13 gave the same
counts.

### Probe warnings

The first full scan of the punch build let 8 warnings escape from
frmtmb.ode's `test-ode-nlf.R` ("Excessive precision requested",
"Returning early"): its exact ridge (lke + t, lV + t) is a removed
direction no objective evaluation had asked about, and the
re-measurement steps along it to 4^20, where the ODE solver complains.
The probes this lane added (`se_obj_curv()`, the gradients of
`se_ridge_curvature()` and `se_span_remeasure()`, the far sd probe,
and `se_curvature_real()`, which now gets steps sized from a
re-measured curvature) evaluate the objective through `se_quiet()`,
which muffles a warning: what the model's code says far from the
estimate is not about the fit. The file then passes with 0 escaped
warnings (`seflags-tgt2`).

### Minors

- m1 (a downward direction the line probe does not confirm keeps
  0.69.0's zeroed loadings): recorded, not changed. Such a direction
  comes from a finite-difference row of a near-zero sd, so its
  eigenvector is noise; read as a flat subspace, that noise would set
  the threshold. Filed in `dev/test-backlog.md`.
- m2: NEWS rewritten (the design alias bullet, the far probe, the
  screen; the degree-7 bullet now gives lmer's SEs).
- m3: the fixture 1 statement is withdrawn in place; the tables are
  within 80 columns.
- m4: filed in `dev/test-backlog.md`, "Filed by lane seflags".


### Reruns on the final build (install 19)

All on install 19 (the last code edit, `se_quiet()`), started after it.
Scans: every test file of the eight suites, one R process per file,
every gate on (`dev/ciharden-scan.sh`); every log names
`wt-seflags-lib` for frmtmb (365 of 365 per scan), and the OpenBLAS
one names `Rob0.3.26-lapack` (365 of 365).

    package           files  pass   ref fail  OB fail  err skip warn
    frmtmb              224  17712         0        1    0    0    0
    frmtmb.coupling      11    548         0        0    0    0    0
    frmtmb.eam           29   1764         0        0    0    0    0
    frmtmb.latent        10    366         0        0    0    0    0
    frmtmb.learn         15    515         0        0    0    0    0
    frmtmb.ode           11    551         0        0    0    0    0
    frmtmb.sample        50   2693         0        0    0    0    0
    frmtmb.spline        15    592         0        0    0    0    0
    total               365  24741         0        1    0    0    0

(`seflags-ref6`, `seflags-ob6`; the pass counts are the same under
both builds.) The one OpenBLAS failure is `test-fuzz.R:59`, seed
20379118, `fit_error`, the 0.69.0 review's B1 (1), lane optrap's
restart. `test-brms-suite-data-helpers.R` passes 6 of 6 under both.

The attack set, lane against base (`bash dev/seflags-sum.sh att
ref6` and `ob6`; differing output lines):

    script                                   ref   OB
    nanse cases (relabeled groups, item 3)    14   14
    b1-rho (SE ratio digits 2e-6 to 4e-6)      4    4
    spread2, grby, pred, b1, falseloss,
      ridge-re                                 0    0
    window designs                             0    0
    aexpb (item 3)                            70   70
    c0k-sweep (item 3)                        38   36
    c0k (seed 77, 3 reps)                      0    0
    deg7 (item 4)                             20   20
    probe (degrees 7 and 8, items 4, B2)       6    6
    fixes-table                                4    0
    smooth2 (item 3, OpenBLAS)                 0    2
    bnd, rs20, item8, trap, short, scale2,
      trapq, quad                              0    0
    sepmiss (budget, few; bigp)             0, 0 0, 0
    singular study, 600 fits (verdicts)        0    0
    small-sd study, 200 fits (verdicts)        0    0
    separation study, 500 fits (verdicts)      0    0

The b1-rho lines are the reported max |SE/ref - 1| in the sixth
digit, the verdicts identical. The fixes-table's two lines under the
reference BLAS are `a * b` (item 3) and `a + exp(b), b ~ 1 + x`, an
identified row that stopped at code 1 in the base process and code 0
in the lane's (finite SEs 3 of 3 either way; the optimizer, not the
check). The probe script's degree-7 rows now keep every SE: without
(1 | g) max |SE/lm - 1| 0.00209, with (1 | g) 0.0938 against lmer();
degree 8 with (1 | g) is "illcond".

The review's scripts, lane against base, both BLAS builds
(`dev/seflags-revrun.sh`, `rev-{ref6,ob6}/`): B1 to B3 and degree 7 as
above; `sd` 0 verdict differences (63 fits per build); `year` the
same as the base (theta_1 "flat" on the year design, the pre-existing
false alarm, filed); `neg` differs only by item 3 (the lane names the
c0k ridges), and its unconfirmed downward directions still have a
largest zeroed loading of 0. The false-alarm battery
(`dev/seflags-falarm-sum.R ref6 ob6`): 41 models, 820 fits per build,
0 verdict differences between lane and base on either build.

Extras (`ext-{ref6,ob6}/`): fixture 1 under OpenBLAS LAPACK
`fitted(newdata)` 0 warnings (base 1); the flat-peak sd "flat" with no
boundary message; the alias construction silent with SEs 0.2 to 0.212
(base 5 of 5 NaN); segap 9 of 9 named under both builds.

R CMD check --as-cran of core after the last edit (built and checked
in `dev/seflags-check/`, `R_LIBS` the lane library, `rellib-r7`, the
user library): `Status: 1 NOTE`, the V8 math-rendering NOTE of the
HTML manual that lane-rules.md expects; tests 15 min, OK.

## First round (as reviewed)

The theme: the standard errors and the warnings must be right. Items
1 to 4 are fixed, each with a test seen to fail on `rellib-r7`. Item 5
is partly fixed (the simplex face) and partly recorded.

## What changed

The first round's changes, as reviewed. Punch round 1 (at the top)
replaced some of them; each is marked.

- `se_tier3()` (R/se-check.R): flat directions kept whole in the null
  basis, with their scaling and threshold (`flat`); one threshold for
  parameters and predictions (punch round 1: only where tau is that
  threshold); tau from the noise along the removed directions;
  unresolved near-null directions joined the removed (withdrawn in
  punch round 1); kept small eigenvalues corrected for a curved ridge;
  `lazy`.
- `jc_nonest()`: reads a flat subspace as one subspace in the
  unit-diagonal metric against its threshold.
- `jc_drop_boundary()`, new: the boundary-column drop of two callers,
  keeping the flat subspace's columns in step.
- `se_ridge_curvature()`, new: the corrected curvature (item 3).
- `se_tier3_could_act()`: screens at `se_ridge_screen` (1e-3), not
  `se_flat_tol` (punch round 1: only for a model with a parameter
  outside its primary linear predictors).
- `cov_from_hessian()`: passes `lazy`.
- `se_edge_sd()`, `se_boundary_names()`: an sd is at zero only if the
  likelihood reads it (`se_sd_at_zero()`, new).
- `se_sd_gain_up()`: also returns the largest change either way
  (`moved`); cached per sd and point (`se_sd_gain_up0()` is the old
  body; punch round 1 added the far probe).
- `sdr_rescue_hessian()`: stores `se_null_flat`.
- R/predict.R, `get_joint_cov()` and `joint_cov_repair()`: carry
  `flat`, its scaling moved to the joint positions.
- R/confint.R, `hyp_par_cov()`, `hyp_prop_var()` and
  `se_null_boundary()`: carry `flat`; boundary drop through
  `jc_drop_boundary()`; a flat subspace is boundary by the parameters
  it names.
- R/brms-shapes.R, `fit_draw_space()` and `fit_fd_se()`: carry `flat`;
  `jc_drop_boundary()`.

No change to `R/fit.R`.

New constants: `se_ridge_screen` (1e-3), `se_ridge_iter` (20),
`se_ridge_step` (0.001), `se_ridge_shrink` (0.5).

## Item 1: an estimable prediction flagged (high)

Cause, as the 0.69.0 review diagnosed: tier 3 put a flat direction in
the null basis with the loadings of the parameters that kept their SE
set to zero (lane nanse's m1 follow-up). What was left of an exact
alias was not flat, and `jc_nonest()` tested it at a cosine of 1e-6.

The rule now. A flat direction is a direction along which the
likelihood is flat, all of it. It goes into the null basis whole,
marked as a subspace with the unit-diagonal scaling of the Hessian it
came from and the threshold that decided which parameters lost their
SE there: `flat_tol = min(tau, se_dominant * max(projection))`. A
prediction is non-estimable when its projection onto the subspace, in
that metric, is above `flat_tol` of its norm. `tau` bounds the angle
between the measured and the true subspace (below), so an estimable
prediction stays under it; and a prediction of a single parameter gets
exactly that parameter's verdict, since a parameter is lost when its
projection is above `flat_tol`. The `se_dominant` cap replaces the
per-column dominance rule of 0.69.0 for flat directions (one rule, so
the two verdicts agree); downward directions keep the 0.69.0 rule,
zeroed loadings and all, so lane nanse's m1 property and
`test-optima.R`'s loading test hold.

Measured on ported fixture 1 under OpenBLAS 0.3.26 with LAPACK
(`dev/seflags-dh7tau.R` on `rellib-r7`, its output on the console):
the flat direction's eigenvalue is -4.66e-5, its own noise 0.00727,
the gap 0.0225, so tau is 28.4 by 0.69.0's formula and 3.2 by the new
one: meaningless either way, and `flat_tol` is the cap,
0.5 * 0.707 = 0.354. The measured direction is within sin 0.0124 of the
exact alias. Its unit-diagonal loadings: Age 0.7074, `s(Age).fx1`
-0.6981, (Intercept) -0.1101, the rest at most 0.0070. The five
predictions' cosines with it: 0.00126, 5.82e-05, 0.000275, 0.00185,
0.00224, all far under 0.354. `fitted(fit1, newdata)` warnings: 1 on
`rellib-r7`, 0 on the lane (`dev/relrev069-dh7diag.R`, logs
`ext-ob3/dh7diag-{base,lane}.txt`). WITHDRAWN (punch round 1, review
B2 and m3): this paragraph said "The SE verdicts are unchanged: Age
and `s(Age).fx1` flat, theta_3 boundary, theta_4 flat". That was false
for the first lane build: its own log `ext-ob3/dh7diag-lane.txt` shows
theta_1 and theta_2 lost as "flat" too, by the joining rule of item 4.
The first round's cap rule also let hypotheses 12 to 30 percent off
the alias pass with finite SEs (review B1). Both are corrected in punch
round 1, at the top.

Built under the reference BLAS (`dev/seflags-alias.R`, seed 1, `y ~ x +
s(x) + (1 | g)`, n = 200, the mean of x set to m; the smooth's fixed
column is exactly `-0.1387 m + 0.1387 x`):

    mean 1e-05: lost x, s(x).fx1 (the intercept keeps its SE)
      rellib-r7: fitted(newdata): 5 of 5 predictions ... NaN NaN NaN NaN NaN
      lane:      fitted(newdata): silent, se 0.2 0.202 0.2 0.212 0.2

At means 0 and 1e-6 the base does not flag (the intercept's share is
below the 1e-6 cosine); from 1e-4 the intercept is lost too and neither
build flags. Test: `test-se-flags.R`, "a prediction along an alias
with a kept parameter is estimable". It checks the construction (x and
`s(x).fx1` lost, the intercept kept), no warning, finite SEs, the SEs
equal to those of `y ~ s(x) + (1 | g)` (the same likelihood without the
alias) within 1 percent, and the absent case: `hypothesis(fit, "x =
0")` is NaN with its warning. On `rellib-r7` it fails 3 expectations
(the warning, NaN SEs, the comparison: `NaN NaN NaN NaN NaN` against
`0.2 ...`).

## Item 2: a boundary message on an sd the likelihood does not read (high)

`test-diagnostics-ux.R` fit 58 (`dev/relrev069-flatsd.R`): the bump is
zero on every row, so the sd of pk's group effect (0.2199) does not
enter the likelihood. The edge test asked only whether a move toward
zero changes nothing, and whether a larger value gains; an sd the
likelihood ignores passes both. A boundary verdict now also needs the
likelihood to change somewhere on `se_sd_gain_up()`'s grid above the
estimate (`se_sd_at_zero()`): at a boundary it changes, by a loss;
an ignored sd does not change at all. The grid is the existing one, so
no new evaluations: the answer is cached per sd and point, which also
removes the duplicate grid that `se_edge_sd()` and
`se_boundary_names()` both asked.

    rellib-r7: MESSAGE Boundary (singular) fit:
               sd_id__pk_Intercept is at zero ...
               WARNING SE not available for 3 of 7: pk, lamp, lsig
    lane:      WARNING SE not available for 4 of 7: pk, lamp, lsig, theta_2
               (no boundary message; diagnose()$singular empty)

(logs `ext-ref3/flatsd-{base,lane}.txt`.) The lane's answer is 0.68.1's.
Test: `test-se-flags.R`, "an sd the likelihood does not read is flat,
not at zero"; it checks the construction (the objective is identical
with that log sd moved by -3 and +1), no `frmtmb_boundary_fit`
message, theta_2 lost as "flat", and the one SE warning naming it. On
`rellib-r7` it fails 3 expectations. The absent case (true boundary
sds still get the message) is the singular study below: 600 fits,
identical verdicts.

## Item 3: curved ridges keep silent finite SEs (medium)

Why 0.69.0 misses them. At a point on a ridge along which the
likelihood is flat, the curvature along the ridge is zero; a fit stops
within its gradient tolerance of the ridge, and there the gradient
across a curved ridge times the ridge's curvature leaks into the
Hessian along it. `dev/seflags-spectrum.R` (seeds 1 to 8, log
`ext-ref3/spectrum-base.txt`): the smallest unit-diagonal eigenvalue
on the curved ridges is 1.2e-10 to 4.1e-7 of the largest (a + exp(b)),
1.6e-7 to 3.7e-7 (its (1 | g) variant), 3.5e-7 to 1.4e-5 (c0k), 2.4e-9
to 2.1e-7 (a mo() simplex with an unobserved category, with and without
age): above `se_flat_tol` (1e-9) on most fits, so tier 1 stood.

The detector, `se_ridge_curvature()`, takes the leak out twice, for
each kept direction at or below `se_ridge_screen` (1e-3) of the
largest:

1. Newton steps on the other directions (those kept with a larger
   eigenvalue, the Hessian held, accepted while the gradient along
   them falls) move the estimate onto the ridge's floor.
2. The curvature along the direction is read there from central
   differences of the gradient at 0.001 and 0.01 unit-diagonal steps.
   A straight line leaves a curved ridge, so on a ridge that
   curvature grows as the square of the step; on real curvature it
   does not. The two are extrapolated to a step of zero.

When the polish leaves more than half the eigenvalue
(`se_ridge_shrink`), the direction curves as the eigenvalue says and
the small step's value is returned unextrapolated. When the corrected
curvature is at or below the flat threshold, the existing step test
(`se_curvature_real()`) decides, at the step that curvature calls for;
a curved ridge fails it (its straight step loses as the fourth power).

`dev/seflags-ridgecurv.R` (seeds 1 to 6, `rellib-r7`, gradient
curvature at steps 1e-3, 1e-2, 1e-1 at the estimate and after the
polish): on a + exp(b) seed 4 the eigenvalue 8.68e-7 reads 2.53e-9,
2.52e-7, 2.52e-5 after the polish (the square law, intercept about
1e-11); a raw polynomial of degree 5 reads 5.55e-11, 5.57e-11, 5.58e-11
(its eigenvalue 5.58e-11), and a poisson pair at x2 = x1 + N(0, 1e-4)
reads 5.61e-9 at all three steps (eigenvalue 5.61e-9).

Measured on the attack set, lane against base (logs `att-ref3/`,
reference BLAS, and `att-ob3/`, OpenBLAS 0.3.26 with LAPACK):

- `dev/setier-rev-aexpb.R`, flat-ridge code-0 fits with every
  coefficient SE finite, the same under both BLAS builds: a + exp(b)
  9 of 20 -> 0 of 20; a * b 3 of 20 -> 0 of 20; a + exp(b) with
  (1 | g) 20 of 20 -> 0 of 20. All of the base's were silent.
- `dev/setier-rev-c0k-sweep.R`, code-0 fits with a flat ridge: with a
  finite SE on a ridge parameter, of those silent: reference 18 fits,
  16 and 16 -> 2 and 0; OpenBLAS 19 fits, 16 and 15 -> 3 and 0. The
  ones left: seed 1 keeps theta_1 (2.88) and seed 7 keeps
  a_(Intercept) (17.8 and 18.1), and under OpenBLAS seed 17 keeps
  a_(Intercept) (3.3), each warned for the other two: their projection
  is under the cap of half the largest (the FD Hessian's tau is
  meaningless there). Recorded below.
- `dev/optima-segap.R` (item 5): 6 of 9 fits warned on the base, 9 of
  9 on the lane (`ext-ref3/segap-{base,lane}.txt`); `ls ~ mo(inc) *
  age` seed 1 now loses both simplexes (both share the unobserved
  category), where the base lost the interaction's only.
- Two curved ridges of the attack set that the base kept with finite
  SEs are now named. `dev/setier-rev-a-rev-cases.R` (lane nanse's
  cases), `y ~ x + (1 | g) + (1 | g2)` with g2 a relabeled g, gaussian
  and poisson: only sd1^2 + sd2^2 enters the likelihood, a circle in
  the two log sds; the base kept theta_1 and theta_2 at 321.8 and
  321.8 (gaussian) and 153.9 and 153.9 (poisson) in silence; the lane
  loses both with the SE warning, under both BLAS builds.
  `dev/setier-rev-smooth2.R`, `s(x) + s(x2)` with x2 = x: only the sum
  of the two penalty variances enters, the same circle; under OpenBLAS
  the lane loses theta_1 and theta_2 (base: SE 166 and 166), under the
  reference BLAS it does not (base and lane alike), so this one still
  depends on rounding.

False alarms. Every identified design of the attack set keeps its
verdicts: the nanse attack scripts other than the two relabeled-group
cases above (spread2, grby, pred, b1, b1-rho, falseloss-sweep,
ridge-re), the window designs, the bound, rs20, trap, short, scale,
item 8 and quadrature scripts, the singular study (600 fits), the
small-sd study (200 fits) and the separation study (500 fits) give identical
verdicts on lane and base under both BLAS builds (table below).
`dev/setier-rev-probe.R`'s raw polynomials of degree 5 and 6 with and
without (1 | g) keep every SE, max |SE/ref - 1| 5.7e-7, 9.2e-6,
2.8e-4 and 4.0e-3 (reference) against lm() and lmer(); degrees 7 and 8
with (1 | g) are item 4.

The first lane build extrapolated from steps 0.01 and 0.1 and took
all six SEs of frmtmb.eam's LBA fits (`test-lba.R:142`,
`test-rdm-gng.R`): the likelihood has a wall 0.3 units along the
smallest direction, and the curvature read 0.0020, 0.0024 and 5.2 at
0.001, 0.01 and 0.1 (`dev/seflags-lbadbg.R`). The full ref scan showed
it as 6 escaped warnings; the steps moved down by ten and the shrink
guard went in, and both files are clean.

Tests: `test-se-flags.R`, "a curved ridge above the flat threshold
loses its SEs" (a + exp(b) seed 4 exact, and seed 1 with (1 | g)), and
"a simplex face off the coordinates is found on every seed" (segap
`mo(inc)` seed 2 and `mo(inc) + age` seed 3, with the absent case:
every category observed keeps every SE in silence). On `rellib-r7`
they fail 4 and 4 expectations. The absent cases of the eigenvalue
gate stay in `test-se-tier.R` ("an ill-conditioned but identified
design keeps every SE", "a weak but identified direction keeps its SE
off stationarity"), which pass.

## Item 4: degree-7 polynomial with (1 | g) (medium)

`dev/seflags-deg7diag.R` (seed 7, log `ext-ref3/deg7diag-*.txt`). The
spectrum: 6.78 ... 3.73e-6, 9.52e-9, 1.37e-11, 8.9e-15; the noise of
the three smallest directions 9.5e-13, 1.8e-12, 2.9e-12; the noise
norm of the whole block 6.5e-8. 1.37e-11 is kept by the step test,
8.9e-15 (below 100 epsilons) is removed. 0.69.0's tau used the whole
block's noise over the gap (1.37e-11): far above 1, so only the
dominant loadings x3 to x5 were named, and the others kept SEs built
without either near-null direction: (Intercept) 329.7 against
lmer()'s 2774, x1 1210 against 13540, x6 105.8 against 1974.

Two changes:

1. tau is now the noise along the removed directions over the gap,
   `||Es V_rem|| / gap` (Davis and Kahan's sin-theta bound in its
   residual form), not the whole block's noise, which is ruled by the
   largest directions.
2. When that tau cannot tell the removed subspace's parameters from
   the rest (above half the largest projection), the kept directions
   the step test rescued join the removed, smallest first, until it
   can. Here the two smallest together have tau 0.0036 against the
   projections (Intercept) 0.022, x1 0.15, x2 0.48, x3 0.71, x4 0.65,
   x5 0.73, x6 0.52, x7 0.14, sigma 3.4e-8, theta_1 6.8e-8.

All eight coefficients then had NaN SEs; sigma and the group sd kept
theirs. WITHDRAWN (punch round 1, review B2 and claim 4): change 2 is
gone, and so is the claim that "no frmtmb SE of these coefficients
would be right". The joining rule removed directions the step test had
confirmed (fixture 1's theta_1 and theta_2 under OpenBLAS LAPACK, an
identified pair beside an alias on 10 of 15 fits), and the design is
identified: second differences of frmtmb's own objective resolve both
directions, and the SEs come within 10 percent of lmer()'s. See punch
round 1 at the top. Note also that this fit is 0.391 below
lmer()'s log-likelihood (5.685 against 6.076), so the two do not sit
at the same point (lmer itself warns "failed to converge" and "nearly
unidentifiable"); not this lane's to fix (optimizer paths are lane
optrap's).

`dev/setier-rev-probe.R`: degree 7 with (1 | g) "lost x3, x4, x5,
finite 5, max |SE/ref - 1| 0.946" -> "lost all 8 coefficients, finite
0"; degree 8 "lost x2 to x7, finite 3, 0.994" -> "lost all 9". Degrees
5 and 6 (with and without (1 | g)) and the 1e5-scaled rows unchanged.
Fixture 1's tau went from 28.4 to 3.2 (still the cap decides). Test:
`test-se-flags.R`, "a near-null direction next to a kept one keeps no
shrunken SE": every coefficient's SE is NaN or above half of lmer()'s,
and sigma and theta_1 keep finite SEs; it fails on `rellib-r7`.

## Item 5 (low)

- **The simplex face off the coordinates: fixed by item 3.**
  `dev/optima-segap.R` on 0.69.0 still missed 3 of 9 (seed 2 of
  `mo(inc)`, seeds 1 and 3 of `mo(inc) + age`); the lane misses 0 of
  9. Still missed, base and lane alike: on `ls ~ mo(inc) * age`, seeds
  2, 3 and 5 of 6 (lane), the interaction's simplex sits at a face
  (moinc:age1[1] = 7.7e-14 on seed 2) and its flat split along the
  face is not found: there the gradient across the face is real (a
  constrained optimum), the polish cannot remove it, and the corrected
  curvature stays 4.8e-7 (`dev/seflags-segapface.R`). Its SE is
  1846 on a weight in [0, 1], so it is not hidden, but it is not
  named. Recorded.
- **A mixture precision with an infinite MLE: not fixed.**
  `dev/seflags-bcmmal.R` on Malingering_2: code 0, no warning on
  either build; log phi2 = 20.39, its unit-diagonal eigenvalues are
  all 0.98 to 1.01 (the row is decoupled), its Hessian diagonal
  1.32e-5 (an SE of 275 on the log scale), and the objective's own
  rounding at that point is about 1e-3 (ciharden's bound 0.00109). No
  test on the Hessian can call it flat: the curvature is real to the
  objective's precision. The remedy would be a limit test per family
  precision (phi -> infinity is the binomial limit), the analogue of
  the sd edge test, against an objective whose rounding is larger than
  the changes it would read. Recorded, not attempted.
- **The cs() ordinal-mixture wall has no message naming its cause:
  not fixed.** Naming it needs a test for a category that closes in
  one component, in the ordinal-mixture code, which is not this lane's
  file set. Recorded.

## Composition with the existing messages

No new message. The ridge and item 4 losses go through `se_report()`
and its per-parameter precedence: on `dev/setier-rev-aexpb.R` (60 fits
per build) no fit has both the nl flat warning and the SE warning
(lane: 11 NLFLAT, 49 SE; base: 11 NLFLAT, 17 SE, 32 silent).
`dev/setier-rev-fixes-table.R`'s "a * b intercepts" row now gets the
nl flat warning (NLFLAT) and NaN SEs, where the base was silent with
finite SEs. A flat subspace whose named parameters are all at their
boundary counts as a boundary column, so a prediction along it does
not repeat the boundary message (`se_null_boundary()`).

## Cost

Counted, not timed (`dev/seflags-cost.R`, gradient and objective
evaluations of `se_check()` on a fit made with the check off; 
`cost-{lane,base}.txt`):

    design            base fn/gr   lane fn/gr
    glmm_uncentered   0 / 6        5 / 14
    lm_poly5          5 / 0        5 / 23
    lmm_poly5         5 / 16       5 / 33
    sleepstudy        0 / 12       0 / 12

The largest cost is on a big fit with a moderately small eigenvalue,
which tier 3 now decomposes: `dev/seflags-bigp.R`, the 2051-coefficient
sparse bernoulli fit of `dev/setier-rev2-sepmiss.R`, alternating builds
three rounds (`bigp-alt.txt`): fit with the check, quiet rounds, lane
24.6 and 25.4 s against base 18.4 and 18.5 s (the check alone is 10
times the fit on both; loaded rounds 32.8 to 57.2 s against 19.7 to
43.6 s). Two changes keep it there: a direction's own noise is
computed only where it can decide anything (`||Es v|| <= ||Es||_F`),
and a covariance that loses nothing is not built (`lazy`). The ode
scale row ran 3175, 3624 and 3439 s on lane builds (in three full
scans) and 4178 s on the base (alone, under load from the scans), so
no regression there.

The small fits of `dev/seflags-cost.R` take 0.01 to 0.06 s with or
without the check on both builds, at the clock's 10 ms resolution, so
their times say nothing; the counts above are the measurement.

## The runs (the first round's final build, install 9)

Every run below is on the final lane build (install 9 of the lane);
the scans run every test file of the eight suites, one R process per
file, with `NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`,
`FRMTMB_DRMTMB_FIT_TESTS`, `FRMTMB_FUZZ` and `FRMTMB_SCALE_TESTS` on
(`dev/ciharden-scan.sh`; `dev/seflags-sum.sh scan <run>` made the
tables). Every log names `wt-seflags-lib` for frmtmb (365 of 365 in
each scan); the extensions load from `rellib-r7`, built from the same
sources, since this lane changes none of them.

Reference BLAS (`seflags-ref3`):

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 224 | 17692 | 0 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1764 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 50 | 2693 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **365** | **24721** | **0** | **0** | **0** | **0** |

OpenBLAS 0.3.26 as BLAS and LAPACK, 4 threads (`seflags-ob3`; every log
names `Rob0.3.26-lapack`): the same counts, except 1 failure in
frmtmb (17692 pass, 1 fail). The failure is `test-fuzz.R:59`, seed
20379118, `fit_error` ("NA/NaN gradient evaluation" after the restart),
the 0.69.0 review's B1 (1), which the release scan had too; it is lane
optrap's restart. `test-brms-suite-data-helpers.R`, the review's other
OpenBLAS LAPACK failure, now passes (6 of 6). The release review's
scan of `rellib-r7` under the same configuration had 364 files, 24693
pass, 2 fail; the lane adds `test-se-flags.R` (27 assertions).

The first full reference scan of the lane (`seflags-ref1`, before the
step change of item 3) had 0 failures and 6 escaped warnings, all in
frmtmb.eam (`test-lba.R` 2, `test-rdm-gng.R` 4): the LBA false alarm
described under item 3, which the final build does not have.

The attack set (`dev/seflags-attacks.sh`, lane against base, `bash
dev/seflags-sum.sh att <tag>`: differing output lines, library lines
and timings left out):

| script | ref | OpenBLAS |
|---|---|---|
| nanse attack cases (`setier-rev-a-rev-cases`) | 14 | 14 |
| spread2, grby, pred, b1, b1-rho, falseloss-sweep, ridge-re | 0 each | 0 each |
| window designs (`setier-rev2-window`) | 0 | 0 |
| `setier-rev-aexpb` | 70 | 70 |
| `setier-rev-c0k-sweep` | 38 | 36 |
| `setier-rev-c0k` (seed 77, 3 reps) | 0 | 0 |
| `setier-rev-deg7` | 12 | 12 |
| `setier-rev-probe` | 4 | 4 |
| `setier-rev-fixes-table` | 2 | 2 |
| `setier-rev-smooth2` | 0 | 2 |
| bnd, rs20, item8, trap, short, scale2, trapq, quad | 0 each | 0 each |
| sepmiss (budget, few; bigp) | 0, 0 | 0, 0 |
| singular study, 600 fits (verdict columns) | 0 | 0 |
| small-sd study, 200 fits (verdict columns) | 0 | 0 |
| separation study, 500 fits (verdict columns) | 0 | 0 |

Every differing line is one of items 3 and 4 above: the curved ridges
now named (aexpb, c0k-sweep, fixes-table's `a * b`, the relabeled
groups of the cases script, smooth2 under OpenBLAS) and degrees 7 and
8 with (1 | g) (deg7, probe).

R CMD check --as-cran of core (built and checked in
`dev/seflags-check/`, `R_LIBS` the lane library, `rellib-r7`, the user
library): `Status: 2 NOTEs`. The examples-timing NOTE
(`dharma_residuals` 6.16 s elapsed against 1.51 s of CPU, measured
during the scans: load) and the V8 math-rendering NOTE of the HTML
manual, which lane-rules.md expects. Tests 25 min, OK.

## Defects found and not fixed

- c0k seeds 1 and 7 (and 17 under OpenBLAS) keep one ridge
  parameter's SE (2.88, 17.8; 3.3) while the other two are named: its
  projection is under the cap of half the largest. A sharper bound on
  a finite-difference Hessian's eigenvector would close it.
- `s(x) + s(x2)` with x2 = x (`dev/setier-rev-smooth2.R`) is named
  under OpenBLAS and not under the reference BLAS, where both builds
  keep SEs of 166 on the two log sds.
- `mo(inc) * age` with the interaction's simplex at a face: the flat
  split along the face is not named (seeds 2, 3, 5 of 6).
- Malingering_2's precision at infinity; the cs() mixture wall
  (item 5).
- The degree-7 fit stops 0.391 below lmer()'s log-likelihood (lane
  optrap's area).

## Brms-suite ledger rows this flips

`data-helpers:7` (fixture 1, `expect_silent(fitted(fit1, newdata))`)
passes under OpenBLAS LAPACK now; under the reference BLAS fixture 1
stops at code 1 and the row passed before and after. Every ported
brms-suite file passes in both final scans, so no gated file shows a
stale "now HOLDS" row; `dev/brmsport-verdicts.tsv` is not edited.

## Version

Behavior changes a user sees, two of them marked BREAKING in NEWS
(standard errors that are now NaN): a minor bump, as for the SE
verdict changes of 0.68.0 and 0.69.0. Floors: nothing exported
changed and no extension needs a new core floor.
