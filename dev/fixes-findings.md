# Lane fixes: emmeans transforms, ledger row 955, trunc()/se() in
# conditional_effects(), ordinal leftovers, the s() null space, the flat warning

Base: 9e902909 (frmtmb 0.67.0), branch `wt-fixes`, uncommitted.
Library: `C:/Users/adf44/source/r/wt-fixes-lib`; "before" is
`rellib-r5`. brms 2.23.0 from the user library. Logs are in
`dev/fixes-log/`, the suite logs in `dev/fixes-suite-p1/` (both local,
gitignored); every number below is pasted from them.

## Punch round 2 (re-check in `dev/reviews/2026-10-05-fixes.md`)

**B2, B3: the pre-fit check is gone; the fitted point warns.** Two
checks before the fit refused identified models: the symbolic one on
`pnorm()` bodies (B1) and the fixed-point Jacobian on natural-unit
bodies, where saturation underflowed columns to 0 (B2: logistic growth
over calendar years, a lapse psychometric at x 200 to 400). The second
also cost O(n^2) (B3: 54 s at n = 32000). Both failure modes come from
judging the model at points that are not the fit's. So
`check_nl_identified()` is deleted, with its call in `frm()`, and
`nl_flat_message()` (`R/fit.R`) reads the outer Hessian block of the
nonlinear parameters' coefficients at the optimum, by central
differences of the exact gradient (2 k gradients for k coefficients),
scaled to unit diagonal. An eigenvalue below `nl_flat_tol = 1e-9` of
the largest is a flat direction, and the fit warns, naming the
coefficients that load on it (|v| > 0.1). It runs on the user's own
call only (not the autoscale pre-fit, `refit()`, bootstrap replicates
or `allFit` members), not under `importance` or `quadrature`, and it
restores TMB's `last.par`, `last.par.best` and `value.best` so that
`sdreport()` reads the optimum. No model is refused before its fit any
more; a fit that cannot reach an optimum stops as before.

The tolerance, measured (`dev/fixes-p2-flat.R`,
`dev/fixes-log/p2-flat.txt`): exact linear ridges `a + b` with shared
`1 + x` 8.65e-17, `a - b` and `a + 2 * b` 2.07e-16, the partial
overlap 2.00e-16; `a * b` with two intercepts 3.90e-11. Curved ridges
whose differenced Hessian is noisier: `a + exp(b)` with intercept-only
`b` 2.31e-07, `c0 + exp(a)^k` 4.58e-07; the BCM product 2.50e-02 (its
ML point sits at a boundary). Every identified model: 3.37e-06 (`a +
exp(b)` with `b ~ 1 + x`, weakly identified) and up; near-collinear
designs (cor 0.9995) 2.64e-04; B2's logistic 9.15e-02 and lapse model
2.79e-01; drmTMB's `exp(lsd) * zz` 5.52e-01. So 1e-9 sits 3400 times
under the weakest identified case and 7 decades over the exact ridges.

False-alarm table on the models of both reviews (`dev/fixes-p2-table.R`,
`dev/fixes-log/p2-table-*.txt`; rellib-r5 beside the lane): on all 16
identified models (B1's three, B2's two, `a + b` disjoint, `a *
exp(b * x)`, `a + exp(b)` with `b ~ 1 + x`, psychometric at x 50 to
150, decay at t 1000 to 5000, Emax, a covariate at 1e3 to 2e3, n = 4,
near-collinear, drmTMB's model, an absent factor level) the lane fits
with no flat warning and with estimates and standard errors equal to
rellib-r5's (B2 logistic SE 0.391 0.361 0.0239, lapse 0.0185 2.25
0.0947). Of 8 unidentified models it warns on 5 (`a + b`, `a - b`, `a +
2 * b`, the partial overlap, `a * b`) and is silent on the 3 curved
ridges (`a + exp(b)` with intercept-only `b`, `exp(a)^k`, the BCM
product), which fit as on rellib-r5, their standard errors NaN or in
the thousands.

Cost (`dev/fixes-p2-cost.R`, `dev/fixes-log/p2-cost.txt`, the
reviewer's `a * exp(b * x) + c0` with `a, b ~ 1 + z`; the whole `frm()`
with the check and with it replaced by a no-op, interleaved in one
process, minimum of 3 rounds; the check alone on the fitted
objective):

    n =   2000: fit with check  0.06 s, without (control)  0.06 s, check alone 0.000 s
    n =   8000: fit with check  0.20 s, without (control)  0.19 s, check alone 0.010 s
    n =  32000: fit with check  1.37 s, without (control)  1.28 s, check alone 0.080 s
    n = 200000: fit with check 11.28 s, without (control) 10.56 s, check alone 0.610 s

Linear in n and about 5 to 7 percent of the fit (the clock ticks at
10 ms, so n = 2000 is below resolution).

Survey (`dev/fixes-nlwatch.sh`, `dev/fixes-nlwatch-run.R`, now wrapping
`nl_flat_message()`; `dev/fixes-log/nlwatch-results.txt`): 64 of 64
files RESULT, 204 calls with a body. Warnings: test-nl.R's own (8), and
2 in frmtmb.ode's test-ode-nlf.R, where `pk_dyn` reads ke and V only as
ke / V and the output is an amount, so `lke + t, lV + t` is exactly
flat (also unidentified on rellib-r5). That test checks two spellings
agree, which they do on the ridge; it now allows that warning, with the
reason in a comment. The recovery start of m1 was offered 0 times.

Seen to fail on the punch-round-1 build (still installed when the
tests were written): test-nl.R pass=60 fail=14 err=2, the two B2
must-fit models refused (`dev/fixes-log/p2-test-nl-p1build.txt`);
lane pass=78 (`dev/fixes-log/p2-test-nl.txt`). The round-1 seen-to-fail
runner `dev/fixes-p0-run-test.R` injects into `check_nl_identified()`,
which no longer exists; it is kept as the record of that round.

**Row :955 now** (`dev/fixes-p2-gam.R`, `dev/fixes-log/p2-gam.txt`):
the ported update stops at its default start ("The nonlinear fit
failed from its default starting values (NA/NaN gradient evaluation)"),
because the fixture has no priors and `mu = a + b` is 0 there on
`Gamma("identity")`. It is not refused for identifiability any more;
started anywhere usable it reaches the ridge and warns (test-nl.R's
Gamma update from the parent fit's estimates warns naming
`a_Intercept, a_x, b_Intercept, b_x`), and with fit2's priors it
converges to 7.149 / 11.585 / 1.647 / -0.793 (m1). The row stays a
defect: the fixture lacks fit2's priors, and `is(up, "brmsfit")` cannot
hold under rule 2 of 2026-09-17.

**m11.** NEWS now reports the reviewer's 1020-fit study
(`dev/fixes-rev2-smooth.R`, `-smooth-cross.R`,
`dev/fixes-rev2-log/smooth-cmp{,2}.txt`) in place of "1.1e-7 on 140
fits": convergence codes other than 0 in 11 (lane) against 16 (base),
non-finite fixef SEs 37 against 31, 9 fits on each build more than 1e-3
below the best of lane, base and mgcv, 6 of them shared; each fit only
one build got wrong reaches the other's optimum from the other's
smoothing SDs. Net neutral on wrong answers; a different set of local
optima. F11 seed 8 (`bf(y ~ s(x1), sigma ~ s(x2))`, lane 0.142 below,
converged, finite SEs) was looked at (`dev/fixes-p2-f11.R`,
`dev/fixes-log/p2-f11.txt`): the lane stops with sigma's smoothing SD
at a local optimum, theta -7.99 and -0.80, where rellib-r5 runs it to
the boundary (-28.88, -9.87; logLik -588.881888 against -589.023726);
`restarts = 3` does not move it, and dropping the smooth units finds
rellib-r5's optimum but brings back the NaN-SE runaway on seeds where
the units are what prevent it (seed 3: conv 1, theta -41.62). No cheap
change avoids it without moving the rest; recorded.

**m12.** `fit$autoscaled` (`R/fit.R`) records whether autoscale's plan
steered the optimizer; `diagnose()` keys "already standardized
internally" on it, not on `par_units`, which a smooth's units also set
(`R/confint.R`). Test in test-smooths.R (the reviewer's
`dev/fixes-rev2-diag.R` case), seen to fail on the punch-round-1 build
(`dev/fixes-log/p2-test-smooths-p1build.txt`: fail=2), lane pass.

### Functions edited in punch round 2

- `R/fit.R`: `check_nl_identified()` deleted with its call in `frm()`;
  new `nl_flat_message()` and `nl_flat_tol`; `fit_assembled()` (the
  call, gated on `announce_start`, and `autoscaled`); the `start` docs.
- `R/confint.R`: `diagnose()`'s print (two places).
- Tests: test-nl.R (the refusal blocks become warning blocks; the two
  B2 must-fit models), test-smooths.R (m12); frmtmb.ode's
  test-ode-nlf.R (allows the flat warning on its unidentified model).
- NEWS.

### Test runs of punch round 2

Every test file of core and the seven extensions on the final
build, one R process per file (`dev/fixes-run-all.sh`,
`dev/fixes-suite-p2/`):

```
# generated from dev/fixes-suite-p2/summary.txt
frmtmb           files 204 pass 16702 fail 0 err 0 skip 1 warn 0
frmtmb.coupling  files  11 pass   542 fail 0 err 0 skip 5 warn 0
frmtmb.eam       files  29 pass  1743 fail 0 err 0 skip 3 warn 0
frmtmb.latent    files  10 pass   359 fail 0 err 0 skip 2 warn 0
frmtmb.learn     files  15 pass   500 fail 0 err 0 skip 2 warn 0
frmtmb.ode       files  11 pass   549 fail 0 err 0 skip 1 warn 0
frmtmb.sample    files  46 pass  2558 fail 0 err 0 skip 1 warn 0
frmtmb.spline    files  15 pass   553 fail 0 err 0 skip 1 warn 0
lib: line names wt-fixes-lib's frmtmb in 341 of 341 logs; gated with FRMTMB_BRMS_FIT_TESTS and FRMTMB_DRMTMB_FIT_TESTS; skips are core test-fuzz.R (1) and each extension's test-scale.R
```

### R CMD check --as-cran, punch round 2

frmtmb, built from the final source (tarball 00:58, last `R/`
edit 00:48) in `dev/fixes-check/frmtmb/`: `Status: 1 NOTE`, the V8
math-rendering note on the HTML manual, as expected.

**Withdrawn on the page:** punch round 1's pre-fit refusal
(`check_nl_identified()`, below) is gone; its claims about what
is refused before the fit no longer hold.

## Punch round 1 (review `dev/reviews/2026-10-05-fixes.md`)

**B1. The nonlinear guard refused identified models.** `stats::D()`
differentiates `pnorm()` and `dnorm()` in their first argument only,
so the symbolic test read a mean and an sd as "identical"
(derivative 0). Replaced by `check_nl_identified()` (`R/fit.R`): the
Jacobian of every nonlinear predictor in the nonlinear parameters'
fixed-effect coefficients, taped by RTMB (exact), at three fixed
points (no random numbers). It refuses only when the Jacobian is rank
deficient at all three AND its null space is the same subspace of the
raw coefficients at all three: a line along which no body changes, at
any value the random effects add, so the likelihood is exactly flat.
That covers `a + b`, `a - b`, `a + 2 * b`, `(a + b) * x`, `log(exp(a)
* exp(b))` and partial overlaps with shared terms (naming only the
coefficients that move, `a_x, b_x`). It fails open, by design, on a
curved ridge (`a * b`, `exp(a)^k`, `plogis(lp) * plogis(lq)`), on a
body calling anything outside a list of elementary functions
(arithmetic, `exp`, `log`, `pnorm`, `dnorm`, `plogis`, ...), on a
coefficient in `betad`, on `mi()`/`me()` and on any tape error.

Why the narrowing (measured, not assumed): the first Jacobian form,
rank deficiency alone, refused four suite models in the reviewer's
survey rerun (`dev/fixes-nlwatch.sh`, first pass): the BCM latent
mixtures case study `plogis(lp) * plogis(lq)` (a scale ridge, which
the test fits for its probabilities), `c0 + exp(a)^k` of
frmtmb.sample's test-laplace-draws.R (`a * k` only), a `ps()` body
(the check dropped the penalized part and made a ridge that does not
exist) and an `frm_ode()` body. The first two are curved ridges, the
last two bodies whose tape the check now does not trust.

On the lane build (`dev/fixes-rev-nl-false.R`,
`dev/fixes-log/nl-false-p1.txt`): the reviewer's three models fit with
the base build's estimates and SEs: response preparation (SE 0.00889
0.00354 0.081), psychometric `pnorm(x, m0, exp(ls))` (0.0454 0.059),
Gaussian bump (0.017 0.0116 0.0102). Refused: `a + 2 * b`, `a - b`,
`b + a`, `(a + b) * x`, `a + b + 5`, `log(exp(a) * exp(b))`, shared
intercept with RE in one, the partial overlap. Fit as before (fail
open): `a * b` (SE 46.4 46.4) and `a + exp(b)` with intercept-only
`b` (SE 39.8, 124), both unidentified on base too.

One more false alarm was found after that, by drmTMB's gated
agreement test (`FRMTMB_DRMTMB_FIT_TESTS=true`, which the ordinary
gated tier skips): `b0 + exp(lsd) * zz` with `zz ~ 0 + (1 | g)` was
refused ("The coefficients lsd_w of the nonlinear parameter 'lsd'
cannot all be estimated"), because the check held the random effects
at zero, where `zz` is 0. It now evaluates each point with the random
effects at fixed nonzero values (one tape per point) and stands down on
a frame whose blocks expand their coefficients (`rr()`). The model is a
must-fit test in `test-nl.R`; the drmTMB file passes 131 of 131
gated (`dev/fixes-log/test-drmtmb-gated-p1.txt`). The failing run's
log was overwritten by the passing one; the message above is quoted
from that run.

False-alarm survey on the final check, the reviewer's 64 files rerun
with it wrapped (`dev/fixes-nlwatch.sh`, `dev/fixes-nlwatch-run.R`,
`dev/fixes-nlwatch-log/`, `dev/fixes-log/nlwatch-results.txt`): 64 of
64 RESULT lines, every one fail=0 err=0; 214 calls with a body;
refusals only in test-brms-suite-methods.R (row 955, 1) and test-nl.R's
own refusal tests (7).

Seen to fail: the three must-fit tests in `test-nl.R` run against the
pre-punch guard, injected into the lane build's namespace
(`dev/fixes-p0-nlcheck.R`, `dev/fixes-p0-run-test.R`,
`dev/fixes-log/test-nl-p0guard.txt`): pass=52 fail=8 err=3, the three
errors the three models' refusal (the failures are the new message
texts the old guard does not have); on the lane build pass=63 fail=0
(`dev/fixes-log/test-nl-p1.txt`). NEWS and `?frm`'s `start` section
now say what the check covers and what it does not (m2), and that a
prior then alone sets the split (m3); the "symbolic and exact proof"
claim is gone.

**m1. :955's starting values.** `fit_recovery_starts()` (`R/fit.R`)
offers, after a failed first attempt only, the prior-placed intercepts
with the nonlinear slopes at zero, when a prior placed a slope and the
objective is finite there. With fit2's brms priors
(`dev/fixes-u955-prior.R`, `dev/fixes-log/u955-prior-*.txt`):
rellib-r5 stops at its start ("NA/NaN gradient evaluation"); the lane
restarts from that start and converges (code 0) to a_Intercept 7.149,
b_Intercept 11.585, a_Age 1.647, b_Age -0.793 (SE 1.692 1.801 1.679
1.740), against brms's posterior means 6.883, 10.963, 1.603, -0.762;
`update(fit2, ..., prior = fit2's priors)` converges too. Moves no
other model: in the 64-file survey the recovery code was reached 4
times (test-diagnostics-ux.R 1, test-predfix.R 3) and offered this
start 0 times. **:955 stays a defect** and does not hold: the ported
fixture fits fit2 without its priors, so the update is the
unidentified model and is refused; and `is(up, "brmsfit")` cannot hold
under the class-name rule of 2026-09-17. Proposed reason: "the
fixture lacks fit2's priors; with them the update converges".

**m2, m3.** Covered by B1's narrowing and the docs above.

**m4. emmeans() on a variable from the formula environment.**
`frame_raw_vars()` (`R/frame.R`) also takes a variable found in the
formula environment with one value per data row. `poly(zz, 2) + f`
with `zz` outside the data (`dev/fixes-p1-misc.R`): rellib-r5
"undefined columns selected"; lane 0.8813863 0.9100011 1.123883
against glm's 0.8813865 0.9100009 1.123884, and at `zz = 1, 3` 0.8994122
1.008988 against 0.8994119 1.008988. Test in `test-emmeans.R`, seen to
fail (rellib-r5 err in that block).

**m5. vcov_cluster() on ordinal fits.** `cluster_scores_at()`
(`R/sandwich.R`) took its template from `obj$env$parameters`, where a
partly mapped component (disc held at 1) has the map's reduced shape;
it now takes `parList()`. `dev/fixes-rev-vcl.R`
(`dev/fixes-log/vcl-p1.txt`): cumulative and sratio, full = TRUE and
FALSE, now run. Test in `test-ordinal-thres-names.R`, seen to fail.

**m6.** The compat cells `disc`, `equidistant` and `sum_to_zero` x
`mixture` = "refused" must flip at consolidation if lane ordmix lands
ordinal mixtures; not changed here. With them, `ord_thres_linpred()`
needs brms's refusal "'incl_thres' is not supported for mixture
models."

**m7.** `confint(parm = "Intercept[2]")` on `cumulative()` now says the
threshold is held through a transform and points to `fixef()` and
`hypothesis(fit, "Intercept[2] = 0")` (`method = "profile"` refuses a
nonlinear combination, so the message does not name it). Test in
`test-ordinal-thres-names.R`.

**m8. Predict estimate__: a divergence everywhere, not only under
trunc().** `dev/fixes-ce-pred.R` (`dev/fixes-log/ce-pred.txt`, brms at
frmtmb's estimates, 4000 draws, seed 808): brms reports the median of
its predictive draws, frmtmb the expected response, on every family.
Poisson at x = -1, 0, 1: brms 1 2 4, frmtmb 1.46066 2.51970 4.34659;
lognormal: brms 1.22559 1.65879 2.24044, frmtmb 1.44327 1.95472
2.64742; gaussian equal within Monte Carlo error. The bands agree
(poisson 0 0 1 | 4 6 9 in both). Recorded as a divergence and not
matched: changing the predict estimate for every family is outside
this lane. The compatibility table has no conditional_effects()
feature to hold a cell, so this entry is the record for the ledger.

**m9. PENDING USER.** `posterior_linpred(incl_thres = TRUE)` on
`hurdle_cumulative()` is refused, a divergence: brms returns `hu` in
column 0 and `(1 - hu) * (tau_k - eta)` in columns 1 to K (the
reviewer's `dev/fixes-rev-thres2.R`: all.equal on 220 rows), which is
`posterior_epred()` at the identity link and no linear predictor.

**m10.** `test-nl.R` compares frm with lm in units of lm's own
standard errors (< 1e-3 SE); `test-incl-thres-draws.R` divides by the
draw's own `max(abs(ref))`.

**New item: the s() null-space scale (BREAKING).** brms builds `Xs`
and `Zs` with `mgcv::smoothCon(diagonal.penalty = TRUE)` (`data_sm()`,
`frame_basis_sm()` with `modCon = 3`); frmtmb did not. `R/frame.R` now
passes `diagonal.penalty = TRUE`. The predict path (`smooth_basis_at()`,
`PredictMat()` with mgcv's stored reparameterization) needed no
change. `dev/fixes-sx.R` (gamSim eg 6, n = 200, seed 1,
`dev/fixes-log/sx-*.txt`): the null-space columns equal brms's standata
`Xs` bitwise (max diff 0) for `s(x1) + s(x2)`, `s(x1, bs = "cr", k =
6)`, `s(x1, by = g)` (factor), `s(x1, by = z)` (numeric), `s(x1, x2)`
and `sigma ~ s(x0)`; on rellib-r5 every one differs (max diff 1.3 to
2.1). The wiggly `Zs` are now bitwise brms's too (base: Z Z' equal to
7.5e-16, columns not). `t2()` was already equal and is unchanged
(mgcv ignores the flag for several penalties). `te()` is refused by
frmtmb, as before. Newdata at the fitted rows reproduces the in-sample
fit to 1.2e-12 relative.

The fit is a reparameterization, but nlminb is not scale invariant.
brms's null-space columns have sd 0.13 to 0.16, so their coefficients
are large, and a smoothing SD that belongs at zero ran far down its
log scale (theta -14 to -112 where base stopped at -6 to -31) until
the outer Hessian was singular and the fixed-effect SEs NaN. Over
gamSim seeds 1 to 20 and seven smooth formulas (`dev/fixes-sx-conv3.R`,
140 fits per build): base conv != 0 in 1, non-finite fixef SE in 3;
the diagonal basis alone 2 and 17; `logLik` lane minus base from
-7.0e-05 to +0.638. The scale is the cause (`dev/fixes-sx-conv6.R`,
seed 3, `s(x1) + s(x2)`, the two columns multiplied by c: c = 0.04 or
0.2 theta -354.89 NaN SE; c = 1 -90.77 NaN; c = 5 -13.19 finite; c =
25 -8.43 finite; logLik equal to 1e-7). So the optimizer now steps the
null-space coefficients in unit-SD coordinates (`smooth_fx_units()`,
`R/autoscale.R`, through the existing `par_units`); the reported
coefficient stays brms's. With it: conv != 0 in 0 of 140, non-finite
fixef SE in 5 (a different set of boundary fits than base's 3),
`logLik` lane minus base from -6.2e-05 to +1.5e-05, at most 1.1e-07
relative. Fits whose SEs come back NaN have a smoothing SD at its
boundary; `fixef()` and `summary()` say nothing there and `vcov()`
warns, as on base (`dev/fixes-sx-conv5.R`), which is a pre-existing
gap filed below.

A prior with `coef = "sx1_1"` now penalizes brms's coefficient: the
column is brms's bitwise, so the parameter is the same one (an
identity). Tests: `test-smooths.R`, the design against brms's standata
bitwise and the fit against `mgcv::gam(method = "ML")`, seen to fail on
rellib-r5 (pass=42 fail=15) and pass on the lane (57);
`test-brms-agreement.R` now asserts the columns equal brms's;
`test-brms-likelihood.R` and the translator's comments say the map is
the identity.

The autoscale pre-fit's convergence test took the outer Hessian in
raw units with a 1e-4 condition cutoff on the coefficient block, and
brms's small null-space column put a smooth's block under it:
test-autoscale.R "autoscale leaves smooth and mo() columns untouched"
then warned that the pre-fit did not converge. The test now reads the
Hessian in the pre-fit's own `par_units`, as `autoscale_fit_sound()`
reads a fit (`autoscale_prefit_converged()`, `R/autoscale.R`); a model
with no smooth has no `par_units` there and is read exactly as before.

Filed, not fixed (`dev/test-backlog.md`): `fixef()` and `summary()`
return NaN standard errors without a warning when the outer Hessian is
singular (pre-existing, 3 of 140 gamSim smooth fits on base); only
`vcov()` warns. And m8's predict estimate divergence.

### Functions edited in the punch round, for the ordmix merge

- `R/fit.R`: `check_nl_identified()` (replaces
  `check_nl_sum_identified()`; the call in `frm()`), `fit_assembled()`
  (`par_units` from `smooth_fx_units()`), `fit_recovery_starts()` (the
  m1 start), the `start` docs.
- `R/frame.R`: the `smoothCon()` call of the smooth loop
  (`diagonal.penalty = TRUE`), `frame_raw_vars()` and its call.
- `R/autoscale.R`: new `smooth_fx_units()`, `autoscale_prefit_converged()`.
- `R/sandwich.R`: `cluster_scores_at()`.
- `R/confint.R`: `resolve_par_index()` (m7).
- Tests: `test-nl.R`, `test-smooths.R`, `test-emmeans.R`,
  `test-ordinal-thres-names.R`, `test-brms-agreement.R`,
  `test-brms-likelihood.R` (comment), `helper-brms.R` (comment);
  frmtmb.sample `test-incl-thres-draws.R`.
- NEWS of both packages; `dev/test-backlog.md`.

### Test runs of the punch round

Every test file of core and the seven extensions, one R process per
file, gated, on the lane build before the last check edit (random
effects at nonzero values); after that edit `test-nl.R`, the gated
drmTMB file and the 64-file survey were rerun, as above.

```
# generated from dev/fixes-suite-p1b/summary.txt
frmtmb           files 204 pass 16552 fail 0 err 0 skip 14 warn 0
frmtmb.coupling  files  11 pass   542 fail 0 err 0 skip 5 warn 0
frmtmb.eam       files  29 pass  1743 fail 0 err 0 skip 3 warn 0
frmtmb.latent    files  10 pass   359 fail 0 err 0 skip 2 warn 0
frmtmb.learn     files  15 pass   500 fail 0 err 0 skip 2 warn 0
frmtmb.ode       files  11 pass   547 fail 0 err 0 skip 1 warn 0
frmtmb.sample    files  46 pass  2558 fail 0 err 0 skip 1 warn 0
frmtmb.spline    files  15 pass   553 fail 0 err 0 skip 1 warn 0
every log's lib: line names wt-fixes-lib's frmtmb: 341 of 341
skips: core test-drmtmb-agreement.R 13 (its own gate, run gated: 131 pass) and test-fuzz.R 1; every extension skip is its test-scale.R
```

Seen to fail on rellib-r5 and pass on the lane, punch round: `test-nl.R`
(against the pre-punch guard, above), `test-smooths.R` (pass=42
fail=15 against 57), `test-emmeans.R` (pass=53 err=4 against 106),
`test-ordinal-thres-names.R` (pass=3 fail=12 err=3 against 34).

Ledger rows: none flips. `brmsfit-methods:955` stays a defect (m1).

### R CMD check --as-cran, punch round

Built and checked in `dev/fixes-check/<pkg>/` from the final source
(tarball 23:46, last `R/` edit 23:36):

- frmtmb: `Status: 1 NOTE`, the V8 math-rendering note on the HTML
  manual, as expected.
- frmtmb.sample: `Status: OK`.


## Round 0

## 1. `emmeans()` on a transformed predictor

**Cause.** `recover_data.frmtmb_fit()` handed emmeans the model frame,
which holds the transformed columns (`poly(z, 2)`, `scale(z)`) and not
`z`. emmeans's `recover_data.call()` selects `all.vars()` of the terms
(`z`, `f`) from that frame, so it stopped with "undefined columns
selected". Fixing only the data then exposed a second defect: the
design route rebuilt the grid's design from the predictor's terms
WITHOUT the fit's `predvars`, so `poly(z, 2)` at the grid's single `z`
stopped ("'degree' must be less than number of unique points") and
`scale(z)` gave `NA` or a grid-centered value (observed while
iterating; that intermediate build was not logged).

**Fix** (`R/interop.R`): `recover_data()` reads `ce_base_frame()`, the
model frame with the raw variables beside it (the frame
`conditional_effects()` already uses); `emm_terms()` and
`emm_basis_design()` patch the fit's frozen `predvars` onto the terms
(`patch_predvars()`, as `pred_design()` does); `emm_grid_vars()`
intersects with the base frame, so the grid route (`epred = TRUE`,
`re_formula = NULL`, nonlinear) holds `z` too.

**Against brms** (`dev/fixes-emm-brms.R`, brms at frmtmb's estimates
through `brms_fixed_fit()`, poisson, n = 120, seed 20260930,
`dev/fixes-log/emm-brms.txt`): every emmean agrees to the 9 digits
printed, and brms's grid holds `z` at 1.08863063 = mean(z), or at the
`at =` values:

    == scale(z) + f, at z = c(0, 4)
      frmtmb  1.31088209 1.98898608 1.62301737 2.30112135 1.65243614 2.33054012
      brms    1.31088209 1.98898608 1.62301737 2.30112135 1.65243614 2.33054012
      glm     1.31088213 1.98898621 1.62301729 2.30112136 1.65243600 2.33054007
    == poly(z, 2) + f, epred
      frmtmb  4.51773417 6.22605762 6.35935028
      brms    4.51773417 6.22605762 6.35935028

The other six cases (`poly(z, 2)`, at `z = c(-1, 2)`, `log(abs(z) +
1)`, `scale(z)`, at `z = 3`, `scale(z)` with an offset) agree the same
way. So brms freezes `scale()` at the training center and scale, as
the brief asked.

**Against glm()/lm()** (`dev/fixes-emm2.R`, `dev/fixes-log/emm2.txt`):
on rellib-r5 all 12 cases stop with "undefined columns selected"; on
the lane build the reference grids are `all.equal()` to emmeans's grids
for the glm/lm fit, and the largest relative emmean difference is
3.42e-07 (the two optimizers' gap). Poisson SEs agree to 1.46e-06. The
gaussian SEs differ by 0.0211, 0.0168 and 0.0253, which is the identity
sqrt(n / (n - p)) of ML against lm's residual df (120/115, 120/116,
120/114 give 0.0215, 0.0170, 0.0260 relative to the larger), not a
defect. The grid route (`dev/fixes-emm4.R`, `dev/fixes-log/emm4.txt`):
`epred = TRUE` against `regrid()` of the glm grid agrees to about 1e-6
relative, with and without `at = list(z = c(0, 3))`; a random
intercept at `re_formula = NULL` holds `z` at mean(z) = 1.0886306; a
nonlinear body reading `log(abs(z) + 1)` equals the hand computation
to the 8 digits printed.

**Tests** (`tests/testthat/test-emmeans.R`, three new blocks): the grid
holds mean(z) and the emmeans equal glm's within 100 times the
coefficient gap (SE within 10 times the SE gap) for four formulas, with
and without `at =`; `scale()` at `z = 3` and `z = c(0, 4)` equals the
training-centered value, also on the grid route; and a gated block
against brms at fixed parameters (`expect_exact_num()`). On rellib-r5:
pass=53 fail=0 err=3, every new block "undefined columns selected"
(`dev/fixes-log/test-emmeans-rellib-r5.txt`); lane: pass=102 fail=0
err=0 gated (`-wt-fixes-lib.txt`).

The gated block needed the brms translator to match brms's renamed
design columns (`polyz21` for `poly(z, 2)1`): `brms_fe_cols()` in
`tests/testthat/helper-brms.R`. Without it Stan's init failed with "no
more scalars to read".

## 2. Ledger row `brmsfit-methods:955`

**Withdrawn by punch round 1 (see the top of this file).** The
symbolic guard described below refused identified models (B1);
`check_nl_identified()` replaced it. The proposed verdict
"divergence" for `brmsfit-methods:955` is withdrawn: the row
stays a defect, and the start placement it names is fixed (m1).

brms's block (`tests.brmsfit-methods.R` lines 951 to 955) calls
`update(fit2, formula. = bf(count ~ a + b, nl = TRUE), testmode =
TRUE)` and asserts `is(up, "brmsfit")`; testmode never fits, and fit2
carries user priors `normal(2, 2)` on `a` and `normal(0, 3)` on `b`
(`brms:::brmsfit_example2$prior`), which brms's `update()` keeps.

**Why the default start fails** (`dev/fixes-u955.R` on rellib-r5,
`dev/fixes-log/u955.txt`). The update keeps fit2's `a ~ Age + (1 | ID1
| patient)` and `b ~ Age + (1 | ID1 | patient)`. At the zero start
`mu = a + b = 0` on `Gamma("identity")`, whose density needs `mu > 0`:
the objective is NaN and every gradient component NaN. So the cause is
the nonlinear starts on an identity link, with no bound.

**But no start fits it.** The body reads `a` and `b` only through
their sum, and both formulas hold an intercept and `Age`. From a
usable start (both intercepts at mean(count) / 2 = 9.8) rellib-r5
converges (code 0) to `a_Intercept = b_Intercept = 9.8991522`, `a_Age =
b_Age = 0.3886577`, every standard error NaN, and moving 1 or 3 units
from `b_Intercept` to `a_Intercept` changes the objective by 0, 0, 0
and -5.68e-14 (an identity, since the body is `a + b`; the residual is
rounding). The outer Hessian has eigenvalues 20.4, 3.72, 2.84, 2.42,
1.85e-06, 6.14e-07, 1.44e-15, -7.75e-16: two exact zeros (intercept,
Age) and two near zero (the shared random-effect block, where only the
sum is identified). The model is unidentified at any start.

**Fix: a clear refusal** (`check_nl_sum_identified()`, `R/fit.R`,
called from `frm()` after the frame is built, not under `dry_run`).
When `stats::D()` gives the same expression for two nonlinear
parameters (`a + b`, `exp(a + b)`, `(a + b) * x`), the body is a
function of their sum; if their fixed-effect designs then have a
joint rank below their column count, a direction of the coefficients
changes no row of the likelihood. The refusal names the aliased
coefficients and the two remedies. It is exempt when a prior moves
along every flat direction (the null space restricted to the
penalized coefficients keeps full rank), and it skips, failing open,
anything `D()` cannot differentiate and any parameter a second body
also reads. No starting value moves: the check reads the frame and
refuses or returns. The ledger fixture on the lane build
(`dev/rel067-update955.R`, `dev/fixes-log/update955-fix.txt`) now
stops with:

    The nonlinear parameters 'a' and 'b' enter `a + b` only through
    their sum, and their formulas share terms, so b_Intercept, b_Age
    cannot be told apart from the other parameters' coefficients: ...

Proposed verdict for `brmsfit-methods:955`: divergence (frmtmb refuses
an unidentified maximum-likelihood model that brms samples only under
fit2's priors; the row's `is(up, "brmsfit")` also cannot hold under
rule 2 of 2026-09-17). It does not flip to holding.

**Tests** (`tests/testthat/test-nl.R`, two blocks): the refusal on
`a + b` with shared `1 + x`, on `exp(a + b) * z`, under a prior on one
coefficient only (the Age direction stays open), and on the update of
a Gamma("identity") fit to `a + b`; and the absent case: disjoint
terms (equals `lm()` to 1e-6 of the largest coefficient),
`a * exp(b * x)`, a prior on every coefficient of `b`, and a parameter
a second body reads. On rellib-r5: pass=49 fail=6 err=0, all six the
refusal cases fitting instead (`dev/fixes-log/test-nl-rellib-r5.txt`);
lane: pass=55 fail=0.

**False alarms.** Every core and extension test file ran on the lane
build (section 6); no file failed with the refusal's text
(`grep "only through their sum"` over `dev/fixes-suite/`, result in
section 6).

**Found and not fixed.** With fit2's brms priors (`dev/fixes-u955.R`
arm 3, `dev/fixes-log/u955-fix.txt`) the model is identified and the
check lets it through, but the prior-located start places `a_Age` at
2 as well as `a_Intercept`, so `mu = 2 + 2 Age` runs from -0.78 to
7.10 (`dev/fixes-log/u955b.txt`) and the fit stops at its start. The
doctrine that a class-wide prior location places every coefficient of
the nonlinear parameter is frmtmb's own; changing it would move the
starts of every nonlinear model with such priors, so it is filed
rather than changed here.

## 3. `trunc()` and `se()` in `conditional_effects()`

**brms** (`dev/fixes-ce-brms.R` at frmtmb's estimates, 4000
Fixed_param draws, gaussian, seed 20261005, 137 rows,
`dev/fixes-log/ce-brms-base.txt` and `-fix.txt`): `prepare_conditions()`
holds every variable of `allvars` at its mean, the addition terms'
too, and the RESPONSE as well: `y = 1.1674001151`, `lo =
-0.7225620131`, `s = 0.6371016667`. An expression bound is evaluated
on the grid rows, so `trunc(lb = min(y) - 1)` is mean(y) - 1 = 0.1674,
not min(y) - 1 = -2.0702.

**frmtmb on rellib-r5** refused `trunc(lb = lo)`, the expression bound
and `se(s)` under `method = "posterior_predict"` ("cannot evaluate ...
Pin lo in conditions"); `trunc(ub = 10)` and the `se()` epred passed.

**Fix** (`ce_aterms()`, `R/conditional-effects.R`): the refusal of an
unpinned variable is gone. The grid already holds every base-frame
column at its reference value (`ce_build_nd()`), the response
included, so the terms evaluate where brms evaluates them. The
"could not evaluate" refusal stays for a term that does not evaluate.
`trials()` keeps its hold at 1, as in brms.

On the lane build every `estimate__` of `method = "posterior_epred"`
equals brms's to the 7 digits printed, for all five models
(truncated means 0.5365971 1.0615407 1.7199458 for `trunc(lb = lo)`
at x = -1, 0, 1; 1.051384 1.371119 1.799699 for the expression bound),
and the held values are identical. `method = "posterior_predict"`:
both packages' bands are Monte Carlo quantiles of the same predictive
distribution, and the new test checks frmtmb's against the exact
quantiles. brms's `estimate__` there is the median of its predictive
draws (0.4512 at x = -1 for `trunc(lb = lo)`), frmtmb's the truncated
mean (0.5366); that is pre-existing and not in this lane's scope.

**emmeans and frmtmb.sample** (`dev/fixes-ce-emm.R`,
`dev/fixes-log/ce-emm-rellib-r5.txt`): `emmeans(epred = TRUE)` had no
such refusal and already equals brms's for `trunc(lb = lo)`, the
expression bound and `se(s)` (brms's grid holds the response, frmtmb's
grid the bound variable; the means agree). frmtmb.sample's
`conditional_effects()` on draws never called `ce_aterms()` and
already held `lo` at its mean on rellib-r5.

**Tests**: `tests/testthat/test-ce-aterms.R` (new): the held values
for both methods; the truncated mean at the held bound for `lb = lo`,
the expression bound and `ub = 10` against the closed form; the
predictive band against the exact truncated-normal and normal
quantiles within 5 Monte Carlo standard errors at 4000 draws (`se(s)`
with sd mean(s), `se(s, sigma = TRUE)` with sqrt(sigma^2 + mean(s)^2));
and a gated block against brms. On rellib-r5: pass=0 err=4, each the
refusal (`dev/fixes-log/test-ce-aterms-rellib-r5.txt`); lane pass=37
gated. `test-simulate-ergonomics.R` asserted the refusal; it now
asserts the held value.

## 4. Ordinal leftovers

### a. `posterior_linpred(incl_thres = TRUE)`

brms (`dev/fixes-thres-brms.R`, `dev/fixes-thres-brms2.R`, logs
`thres-brms.txt`, `thres-brms2.txt`): an array draws x observations x
thresholds, dimnames `list(c("1", ...), NULL, c("1", ...))` (the draw
labels are 1..n under `draw_ids = c(2, 5)` too), layer k
`disc * (thres_k - mu)` for cumulative and sratio and
`disc * (mu - thres_k)` for cratio and acat, `mu` including the `cs()`
part at each threshold, `NA` past a `thres(gr = )` level's thresholds;
ignored with `dpar`, with `transform = TRUE`; `incl_thres = NA` an
error. `hurdle_cumulative` returns `hu` as a column "0" beside the
threshold predictors times `1 - hu` (0.215 in column 0 at the fixture),
which is no linear predictor; frmtmb refuses it and says why (a
divergence for the user to confirm). brms's sum_to_zero cumulative
logit program does not compile (`ordered_logistic_glm_lpmf` with an
`int` cutpoint, `dev/fixes-log/thres-brms.txt`), so that shape was run
as sratio and as cumulative probit.

**Fix**: core `ord_thres_linpred()` (`R/predict.R`, exported through
the sampling API in `R/sampling-api.R`), and frmtmb.sample's
`posterior_linpred.frmtmb_draws()` calls it per draw and stacks brms's
array. frmtmb's values at the same estimates equal brms's for the
cumulative and grouped shapes to every printed digit
(`dev/fixes-thres-probe.R`; newdata rows in `thres-brms2.txt`). Core
has no fit-level `posterior_linpred()` (brms's `fitted()` has no
`incl_thres`), so there was no fit-level analogue to change.

**Tests**: `extensions/frmtmb.sample/tests/testthat/test-incl-thres-
draws.R` (new): shape, dimnames, and two draws of each of cumulative,
sratio with `cs()`, cratio with `disc ~ 0 + z`, acat equidistant,
sratio sum_to_zero and grouped cumulative against brms's own
`dcumulative()`/`dsratio()`/`dcratio()`/`dacat()` at `link =
"identity"` on the draw's stored brms-named columns, to 1e3 epsilon;
the ignore rules; the hurdle refusal. rellib-r5: pass=0 err=4, each
the old refusal; lane: pass=43. `test-draws-spellings.R` asserted the
refusal on a gaussian; brms ignores the argument there, and the test
now says so.

### b. Per-threshold `Intercept` rows

brms (`dev/fixes-prior-brms.R`, `dev/fixes-log/prior-brms.txt`): under
flexible thresholds (cumulative, sratio, hurdle) `default_prior()`
lists `Intercept` and then `coef` "1" to "4"; grouped, the class row,
then per level the level row and its own count of coef rows (a: 1-4,
b: 1-3); under equidistant only the class row, and `coef = "1"` is
refused ("The following priors do not correspond to any model
parameter"); `prior(coef = "2")` gives `normal_lpdf(Intercept[2] |
...)` in the Stan code, `Intercept_1[2]` for level a.

**Fix** (`R/priors.R`): `default_prior()` lists the rows;
`ordinal_threshold_entry()` takes `coef`, and the resolver merges a
coef row into the vector's entry as a per-element density (dist kind
`"vec"` in `prior_base_logdens()`), so the ordered map's log-Jacobian
enters once. Refused, each by name: coef under equidistant (as brms),
under sum_to_zero (the class itself is refused there, unchanged from
0.67.0, see below), coef without `group` under grouped thresholds, a
number past the vector, and a coef that is not a number.

**Not done as asked: sum_to_zero rows.** brms lists the coef rows
under sum_to_zero too, but its `Intercept` there is the vector it
declares BEFORE centering, which frmtmb does not have; lane ordinal
refused the class for that reason, and listing rows the resolver then
refuses would advertise a slot nothing fills. Kept refused.

**Tests**: `tests/testthat/test-ordinal-thres-names.R` (new): the rows
against brms's lists; the MAP penalty at the estimates equals the
closed form (class + coef, coef alone, unordered, grouped) to 1e-10
relative; a tight coef prior moves the threshold more than ten times
closer to its location; the refusals. rellib-r5: pass=2 fail=11 err=2
(`dev/fixes-log/test-otn-rellib-r5.txt`); lane: pass=26.

### c. `confint()` names

rellib-r5 (`dev/fixes-confint-ord.R`, `confint-ord-base.txt`): rows
`tau_raw_1` .. `tau_raw_k` on every ordinal fit, and `parm =
"Intercept[1]"` or `"delta"` refused. The internal parameters are not
all the user-facing ones: `cumulative()` holds log increments, and
its equidistant distance as log(delta). Naming such a row
`Intercept[2]` would put a name on a value it does not have, so each
row is named by what it is (`outer_par_labels()`, `R/confint.R`):
`Intercept[k]`, `Intercept[a,k]`, `delta` where the internal value is
that parameter; `log(Intercept[2] - Intercept[1])`, `log(delta)` where
it is a transform; a `cs()` row its `fixef()` name. `vcov(full = TRUE)`
and `vcov_cluster(full = TRUE)` carry the same names, so confint's rows
still equal them. `parm =` takes the label, brms's `b_Intercept[1]` and
the old `tau_raw_k`; `parm = "delta"` on an ordered family addresses
`log(delta)` with a message, as the sd aliases do. The internal
`outer_par_names()` is unchanged, so nothing that matches by those
names moved (35 call sites; marginaleffects' `get_coef()` keeps
`tau_raw_k`). Lane output in `confint-ord-fix.txt`. Tests in
`test-ordinal-thres-names.R`.

### d. Compatibility rows

`R/compat.R`: features `disc`, `equidistant` and `sum_to_zero` (kind
"structure"), and rules: refused on every other family, works on the
ordinal group, works with `thres()`, refused in a mixture, works with
fitted, predict, simulate, residuals_osa, emmeans, REML, quadrature,
confint_profile and mvbf (each run by `dev/fixes-compat-probe.R`,
`dev/fixes-log/compat-probe.txt`; the ones the ordinal tests compare
cite them), prior works for disc and equidistant and is conditional
for sum_to_zero. A new block of `test-compat.R` checks the cells and
that the refusals they claim are the ones `frm()` makes. rellib-r5
fails it with "Unknown feature: 'disc'" (the weak form: the defect
was the missing rows). frmtmb.eam's tests, which read the table, ran in
section 6.

### e. A fixed `disc` in `variables()`

Measured (`dev/fixes-disc-vars.R`, `dev/fixes-log/disc-vars.txt`):
brms's program declares `real disc = 1;` in transformed parameters
when disc is not modeled, so `variables()` lists `disc` (seen on every
brms fixed fit, `thres-brms.txt`); with `disc ~ 0 + z` it is a local of
the model block and only `b_disc_z` is listed. frmtmb lists neither
`disc` on the fit nor on frmtmb.sample's draws. Not changed: a `disc`
variable would need a constant column in the draws (`as.matrix()`,
`summary()`, `hypothesis()`), a second rule beside the user's decision
of 2026-09-30 to hide the fixed disc, and it touches the draws columns
lane ordmix edits. Recorded. brms also lists `Intercept[k]` beside
`b_Intercept[k]` (the centered thresholds), which frmtmb does not.

## 5. Functions edited, for the ordmix merge

- `R/fit.R`: `frm()` (one call), new `check_nl_sum_identified()`, the
  `start` documentation.
- `R/interop.R`: `recover_data.frmtmb_fit()`, `emm_terms()`,
  `emm_basis_design()`, `emm_grid_vars()`; a doc line of
  `interop_coef_names()`.
- `R/conditional-effects.R`: `ce_aterms()` (signature now `(rspec, nd,
  n)`) and its two callers in `conditional_effects.frmtmb_fit()`; the
  `method` documentation.
- `R/predict.R`: new `ord_thres_linpred()` after
  `ord_linear_per_threshold()`. Nothing else in the ordinal code.
- `R/priors.R`: `ordinal_threshold_entry()` (inside the resolver), the
  resolver's ordinal branch, `prior_base_logdens()`, the default rows
  in `default_prior()`'s ordinal branch, the set_prior docs.
- `R/confint.R`: new `outer_par_labels()` and `ord_internal_labels()`,
  `resolve_par_index()`, `confint.frmtmb_fit()` (row names).
- `R/methods-fit.R`: `vcov.frmtmb_fit()` (dimnames, docs).
- `R/sandwich.R`: `vcov_cluster()` (dimnames).
- `R/compat.R`: `compat_features_build()` and the rule list.
- `R/sampling-api.R`: the export list and docs.
- frmtmb.sample `R/methods-draws.R`: `posterior_linpred.frmtmb_draws()`
  and its `incl_thres` doc. Note for ordmix: brms refuses `incl_thres`
  on an ordinal MIXTURE ("'incl_thres' is not supported for mixture
  models."); frmtmb has no ordinal mixtures today, so the branch only
  sees `type == "ordinal"` families. When ordmix adds them, that
  refusal belongs in `ord_thres_linpred()`.
- Not touched: `R/families.R`, `R/thres.R`.

## 6. Test runs

Every test file of core and the seven extensions, one R process per
file, 20 at a time, gated, on the lane library (`dev/fixes-run-all.sh`,
`dev/fixes-run-test.R`). The two core files that failed on the first
run asserted the old `confint()` names (`bcs2_k`, `o_tau_raw_1`); they
and the case-studies vignette now use the new names. `git diff` of
the vignette is the prose and the one index expression.

```
# generated from dev/fixes-log/core-results.txt and dev/fixes-suite-ext/summary.txt
core (gated: NOT_CRAN, FRMTMB_BRMS_FIT_TESTS), first run:
  files 204 pass 16502 fail 1 err 1 skip 14 warn 0
  not clean: RESULT test-case-studies.R pass=28 fail=0 err=1 skip=0 warn=0
  not clean: RESULT test-mv-gaps.R pass=57 fail=1 err=0 skip=0 warn=0
  re-run after the name updates:
    RESULT test-case-studies.R pass=30 fail=0 err=0 skip=0 warn=0
    RESULT test-mv-gaps.R pass=58 fail=0 err=0 skip=0 warn=0
    RESULT test-drmtmb-agreement.R pass=131 fail=0 err=0 skip=0 warn=0
  skips: test-drmtmb-agreement.R 13 (gate FRMTMB_DRMTMB_FIT_TESTS, run gated above), test-fuzz.R 1 (FRMTMB_FUZZ)
  files whose log holds the new refusal text: 0
extensions (gated), every log's lib: line names wt-fixes-lib's frmtmb:
  frmtmb.coupling  files  11 pass  542 fail 0 err 0 skip 5 warn 0
  frmtmb.eam       files  29 pass 1743 fail 0 err 0 skip 3 warn 0
  frmtmb.latent    files  10 pass  359 fail 0 err 0 skip 2 warn 0
  frmtmb.learn     files  15 pass  500 fail 0 err 0 skip 2 warn 0
  frmtmb.ode       files  11 pass  547 fail 0 err 0 skip 1 warn 0
  frmtmb.sample    files  46 pass 2558 fail 0 err 0 skip 1 warn 0
  frmtmb.spline    files  15 pass  553 fail 0 err 0 skip 1 warn 0
  every extension skip is its test-scale.R (FRMTMB_SCALE_TESTS)
```

## 7. Versions

Core: minor (a new refusal and new behavior, an export added to the
sampling API, and breaking name changes in `confint()` and
`vcov(full = TRUE)`). frmtmb.sample: minor (a new capability of
`posterior_linpred()`), its floor raised to the next frmtmb, for
`ord_thres_linpred()`.

## 8. Ledger rows

No ported-suite row flips. The gated brms-suite files of core and
frmtmb.sample pass as recorded (core `test-brms-suite-methods.R`
pass=162, `test-brms-suite-emmeans.R` pass=11, `test-brms-suite-priors.R`
pass=36, no stale "now HOLDS"). `brmsfit-methods:955` changes its
reason, not its state: it is a refusal of an unidentified model now,
and the proposed verdict is in section 2. No other ported assertion
changed state on the lane build.

## 9. R CMD check --as-cran

Built and checked in `dev/fixes-check/<pkg>/` (`dev/fixes-check.sh`),
R_LIBS = lane library, rellib-r5, user library:

- frmtmb: `Status: 1 NOTE`, the V8 math-rendering note on the HTML
  manual ("Skipping checking math rendering: package 'V8'
  unavailable"), which lane-rules.md lists as expected.
- frmtmb.sample: `Status: OK`.
