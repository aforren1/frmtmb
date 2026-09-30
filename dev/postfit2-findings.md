# Lane wt-postfit2: brms post-fit functions

Worktree `C:\Users\adf44\source\r\frmtmb-wt-postfit2`, branch
`wt-postfit2`, base `1f40800d` (frmtmb 0.65.0, frmtmb.sample 0.13.0).
Private library `C:/Users/adf44/source/r/wt-postfit2-lib`; the "before"
arm is the base build `C:/Users/adf44/source/r/rellib-r3`. brms 2.23.0's
source was read from its installed namespace
(`dev/postfit2-brmssrc.R`, `dev/postfit2-brmssrc2.R`, dumped into
`dev/postfit2-brmssrc/`), not recalled.

## 1. What changed

1. `conditional_smooths()`: a new generic (owner brms) and the fit
   method in frmtmb, `R/conditional-smooths.R`; the draws method in
   frmtmb.sample, `R/postfit-draws.R`.
2. `conditional_effects()` options: `surface`, `too_far`,
   `select_points` and `spaghetti` on the fit method, and the draws
   method takes the same four, which it refused. `make_conditions()`
   is in `R/brms-utilities.R`.
3. The nonlinear Wald band: the delta method through
   `frm_lp_basis()`, in frmtmb. The draws method never refused it.
4. `posterior_average()`: a new generic (owner brms) and the draws
   method, in frmtmb.sample.
5. `update_adterms()`: `R/brms-utilities.R`.

Supporting hunks, each small:

- `R/parse.R`: a smooth spec keeps its written term as the attribute
  `frm_term` (whitespace removed), which is brms's name for the term.
- `R/frame.R`: `sm_info` carries that as `term`, so the per-level bases
  of `s(z, by = f)` group back into one term.
- `R/predict.R`: the smooth-basis block of `pred_design()` moved into
  `smooth_basis_at()`, unchanged, so that `conditional_smooths()` draws
  exactly the basis prediction uses.
- `R/generic-owners.R`, `extensions/frmtmb.sample/R/generic-owners.R`:
  `conditional_smooths` and `posterior_average` are owned by brms, so
  they are active bindings like the other borrowed generics. The
  generic-collision tests of both packages restate the tables and were
  updated with them.
- `R/sampling-api.R`: the new engine functions join the extension API.
- `R/conditional-effects.R`: the options, the nonlinear band, the
  `surface`/`spaghetti` plotting, brms's `"points"` attribute, and the
  help page (the false `conditional_smooths()` sentence removed; the
  draws section said "mean and standard deviation" where the default
  is the median and MAD, fixed).

## 2. Seen to fail on the base build

`dev/postfit2-before.R`, log `dev/postfit2-log/before.txt` (seed 1,
n = 150, `y ~ f + s(x) + s(z, by = f)`):

    exists conditional_smooths: FALSE
    conditional_effects(fit, 'x') estimate__ minus the s(x) term: min 0.0213, max 0.0213
    its se__ over the term's own delta-method se: min 1.027, max 1.214
    nonlinear fit, default band: conditional_effects() cannot put a wald band on a nonlinear predictor ...
    spaghetti / select_points / too_far: conditional_effects() has no argument ...
    surface = TRUE: conditional_effects(surface = TRUE) is not implemented ...
    exists make_conditions: FALSE
    exists update_adterms: FALSE

So the help page's claim that `conditional_effects()` "also covers what
brms calls `conditional_smooths()`" was false as a behavior, not only
as a missing name: the display is the expected response, offset from
the term by the intercept and the held terms, with a wider band.

The new test files run on the base build (`dev/postfit2-runtest.R
<pkg> <file> --base`, logs `dev/postfit2-log/before-core-*.txt`):

```
RESULT frmtmb test-conditional-smooths.R: tests=0 failed=0 error=6 skipped=0 warning=0 passed=0
RESULT frmtmb test-brms-utilities.R: tests=0 failed=0 error=2 skipped=0 warning=0 passed=0
RESULT frmtmb test-ce-options.R: tests=2 failed=1 error=5 skipped=0 warning=0 passed=1
RESULT frmtmb test-ce-bands.R: tests=162 failed=0 error=2 skipped=0 warning=0 passed=162
RESULT frmtmb test-nlf.R: tests=69 failed=0 error=1 skipped=0 warning=0 passed=69
RESULT frmtmb.sample test-postfit-draws.R: tests=0 failed=0 error=4 skipped=0 warning=0 passed=0

first failures of each, base build:
[before-core-test-conditional-smooths.R.txt]
---- expectation_error in conditional_smooths() draws each term with brms's keys and grid 
Error in `conditional_smooths(fit)`: could not find function "conditional_smooths" 
---- expectation_error in the drawn term is the fitted term and its band the delta method 
Error in `conditional_smooths(fit, resolution = 25)`: could not find function "conditional_smooths" 
[before-core-test-brms-utilities.R.txt]
---- expectation_error in make_conditions() crosses levels and mean +/- sd 
Error in `make_conditions(d, c("a", "f"))`: could not find function "make_conditions" 
---- expectation_error in update_adterms() replaces, adds and removes addition terms 
Error in `update_adterms(form, ~trials(10))`: could not find function "update_adterms" 
[before-core-test-ce-options.R.txt]
---- expectation_error in surface = TRUE draws both predictors over their range 
<frmtmb_error/error/condition>
---- expectation_failure in select_points keeps the observations near the held values 
Expected `names(all_pts)` to be identical to `c("f", "resp__", "effect1__")`.
[before-core-test-ce-bands.R.txt]
---- expectation_error in a nonlinear fit finds the covariates of its nl body 
<frmtmb_error/error/condition>
---- expectation_error in surface = TRUE draws the two-variable effect over both ranges 
<frmtmb_error/error/condition>
[before-core-test-nlf.R.txt]
---- expectation_error in post-processing follows the nonlinear parameter a body names 
<frmtmb_error/error/condition>
[before-sample-test-postfit-draws.R.txt]
---- expectation_error in conditional_smooths() on draws summarizes the per-draw terms 
Error in `conditional_smooths(cs$d1, resolution = 15, ndraws = 12, spaghetti = TRUE)`: could not find function "conditional_smooths" 
---- expectation_error in conditional_effects() on draws takes spaghetti and surface 
<frmtmb_sample_error/frmtmb_error/error/condition>
```

## 3. Against brms at the same parameters

`dev/postfit2-brms-compare.R`, log `dev/postfit2-log/brms-compare.txt`.
brms is sampled with `algorithm = "fixed_param"`, one chain per
parameter vector, each initialized at a vector mapped from frmtmb, so
every brms draw IS that vector. A smooth's coefficients are mapped by
projecting frmtmb's term, evaluated at the data, on brms's columns for
the same mgcv basis. The projection residual is at rounding level
(3.9e-15 relative on the ML fit, 6.6e-16 on `t2()`, 5.2e-15 over ten
draws), so the two bases span the same space and the comparison means
something. Data seed 1, `frm_sample()` seed 20260929, posterior
average seed 5.

```
== A. make_conditions(), update_adterms(): identical() to brms
make_conditions(list(c("zBase", "zAge"))): identical TRUE, dim 9x3
make_conditions(list(c("zBase", "Trt"))): identical TRUE, dim 6x3
make_conditions(list(c("fz", "zAge"), digits = 1)): identical TRUE, dim 6x3
make_conditions(list("Trt", incl_vars = FALSE)): identical TRUE, dim 2x2
make_conditions(list(c("zAge", "zBase"), sep = " | ")): identical TRUE, dim 9x3
update_adterms(y | trials(size) ~ x, ~trials(10)) = y | trials(10) ~ x; identical TRUE
update_adterms(y | trials(size) ~ x, ~weights(w)) = y | trials(size) + weights(w) ~ x; identical TRUE
update_adterms(y | trials(size) ~ x, ~weights(w), replace) = y | weights(w) ~ x; identical TRUE
update_adterms(y ~ x, ~trials(10)) = y | trials(10) ~ x; identical TRUE
update_adterms(y | se(s, sigma = TRUE) + weights(w) ~ x + (1 | g), ~se(s2)) = y | weights(w) + se(s2) ~ x + (1 | g); identical TRUE
update_adterms(y | resp_se(s) ~ x, ~se(s3) + cens(c)) = y | se(s3) + cens(c) ~ x; identical TRUE
update_adterms(y | cens(c, y2) ~ x, ~weights(w), replace) = y | weights(w) ~ x; identical TRUE
update_adterms(y | weights(w) ~ x, ~1, replace) = y ~ x; identical TRUE
update_adterms(log(y) | trunc(lb = 0) ~ s(x), ~trunc(ub = 5)) = log(y) | trunc(ub = 5) ~ s(x); identical TRUE

== B1. conditional_smooths(): ML fit vs brms at the ML estimate
projection residual (relative): 3.86e-15
default: keys brms [mu: s(x), mu: s(z,by=f)] frmtmb [mu: s(x), mu: s(z,by=f)]
  mu: s(x): rows 100 / 100; columns identical TRUE; grid max diff 0
    estimate__ max |brms - frmtmb| = 3.77e-15 (scale 1.03)
  mu: s(z,by=f): rows 200 / 200; columns identical TRUE; grid max diff 0
    estimate__ max |brms - frmtmb| = 1.67e-16 (scale 0.589)
surface = FALSE, int_conditions: keys brms [mu: s(z,by=f)] frmtmb [mu: s(z,by=f)]
  mu: s(z,by=f): rows 4 / 4; columns identical TRUE; grid max diff 0
    estimate__ max |brms - frmtmb| = 1.67e-16 (scale 0.308)
brms s3: No valid smooth terms found in the model.
frmtmb s3: No valid smooth terms found in the model.

== B2. t2(x, z): surface and too_far
projection residual (relative): 6.57e-16
t2 surface, resolution 20, too_far 0: keys brms [mu: t2(x,z)] frmtmb [mu: t2(x,z)]
  mu: t2(x,z): rows 400 / 400; columns identical TRUE; grid max diff 0
    estimate__ max |brms - frmtmb| = 1.44e-15 (scale 1.41)
t2 surface, resolution 20, too_far 0.1: keys brms [mu: t2(x,z)] frmtmb [mu: t2(x,z)]
  mu: t2(x,z): rows 390 / 390; columns identical TRUE; grid max diff 0
    estimate__ max |brms - frmtmb| = 1.33e-15 (scale 1.41)
t2 surface = FALSE: keys brms [mu: t2(x,z)] frmtmb [mu: t2(x,z)]
  mu: t2(x,z): rows 60 / 60; columns identical TRUE; grid max diff 0
    estimate__ max |brms - frmtmb| = 1.33e-15 (scale 1.16)
conditional_effects x:z surface, too_far 0.2: keys brms [x:z] frmtmb [x:z]
  x:z: rows 225 / 225; columns identical TRUE; grid max diff 0.647 (in se__, lower__, upper__: brms 0 NA frmtmb 0.2605 0.2222)
    estimate__ max |brms - frmtmb| = 1.5e-15 (scale 1.57)

== B3. select_points: the observations each keeps
select_points 0.00: brms 150 rows, frmtmb 150 rows, same rows TRUE; frmtmb plus brms's response filter 150 rows, same as brms TRUE
select_points 0.05: brms 1 rows, frmtmb 1 rows, same rows TRUE; frmtmb plus brms's response filter 1 rows, same as brms TRUE
select_points 0.10: brms 3 rows, frmtmb 5 rows, same rows FALSE; frmtmb plus brms's response filter 3 rows, same as brms TRUE
select_points 0.30: brms 45 rows, frmtmb 58 rows, same rows FALSE; frmtmb plus brms's response filter 45 rows, same as brms TRUE
estimate__ of x at the ML estimate: max diff 2.22e-16

== B4. nonlinear predictor: estimate vs brms, band vs Monte Carlo
brms estimate__ vs frmtmb estimate__, z: max diff 0 (scale 4.36)
brms estimate__ vs frmtmb estimate__, z:f: max diff 0
R = 4000, seed 20260929: se__/MC sd over 100 grid points: min 0.9767, median 0.9812, max 1.0110
  |z| of (se__ - MC sd) against the MC error of an SD: max 2.08, points beyond 2: 21, beyond 3: 0
  band ends vs MC 2.5%/97.5% quantiles: max |diff| / MC width 0.0154

== C. draws at the same parameter vectors
projection residual over the 10 draws (relative): max 5.22e-15
brms draws: 10
conditional_smooths on draws: keys brms [mu: s(x), mu: s(z,by=f)] frmtmb [mu: s(x), mu: s(z,by=f)]
  mu: s(x): rows 100 / 100; columns identical TRUE; grid max diff 0
    estimate__ max |brms - frmtmb| = 3.55e-15 (scale 1.02)
  mu: s(z,by=f): rows 200 / 200; columns identical TRUE; grid max diff 0
    estimate__ max |brms - frmtmb| = 5.13e-16 (scale 0.641)
  mu: s(x) spaghetti: rows 1000 / 1000, names identical TRUE, estimate__ max diff 5.33e-15, sample__ identical TRUE
  mu: s(z,by=f) spaghetti: rows 2000 / 2000, names identical TRUE, estimate__ max diff 2e-15, sample__ identical TRUE
conditional_effects(x, spaghetti): rows 1000 / 1000, names brms [x,y,f,z,cond__,effect1__,estimate__,sample__] frmtmb [x,y,f,z,cond__,effect1__,estimate__,sample__]
  estimate__ max diff 5.36e-15; sample__ identical TRUE; band estimate__ max diff 4.27e-15
conditional_effects(x:f, spaghetti): rows 2000 / 2000, names brms [x,f,y,z,cond__,effect1__,effect2__,estimate__,sample__] frmtmb [x,f,y,z,cond__,effect1__,effect2__,estimate__,sample__]
  estimate__ max diff 5.44e-15; sample__ identical TRUE; band estimate__ max diff 4.27e-15
brms spaghetti + surface: Cannot use 'spaghetti' and 'surface' at the same time.
frmtmb spaghetti + surface: Cannot use 'spaghetti' and 'surface' at the same time.

== C2. posterior_average(): same draws, weights and seed
variables present: brms TRUE, frmtmb TRUE
weights 0.3/0.7: dim brms 10x3 frmtmb 10x3; names identical TRUE; max diff 1.39e-17; attr ndraws brms 3/7 frmtmb 3/7; weights identical TRUE
weights 1/0: dim brms 10x3 frmtmb 10x3; names identical TRUE; max diff 1.39e-17; attr ndraws brms 10/0 frmtmb 10/0; weights identical TRUE
weights 0.5/0.5: dim brms 10x3 frmtmb 10x3; names identical TRUE; max diff 1.39e-17; attr ndraws brms 5/5 frmtmb 5/5; weights identical TRUE
ndraws = 7: brms 7x3, frmtmb 7x3, max diff 1.39e-17
names(attr(, 'weights')): brms bK,bK | frmtmb dsK,dsK
unknown variable: brms 'Parameters 'nope' cannot be found in all of the models. Consider using argument 'missing'.' | frmtmb 'Parameters nope cannot be found in all of the models. Consider using argument 'missing'.'
missing = 0, a variable no model has: brms 'Parameters 'nope' cannot be found in any of the models.' | frmtmb 'Parameters nope cannot be found in any of the models.'
missing = 0 across two models: dims 10x2 / 10x2, max diff 0, ndraws 4/6 / 4/6
done
```

Reading it:

- `conditional_smooths()` on the ML fit: keys, rows, columns, `cond__`,
  the grid and the attributes `effects`, `response`, `surface` and
  `points` are identical to brms's; `estimate__` agrees to 3.8e-15 on
  `s(x)`, `s(z, by = f)`, `t2(x, z)` (surface, `too_far = 0.1`, which
  keeps 390 of 400 points in both) and `surface = FALSE`. `se__`
  differs by construction: brms's is the MAD of identical draws, 0;
  frmtmb's is the delta method.
- On ten draws: every column of `conditional_smooths()` equals brms's
  to 3.6e-15, and both spaghetti frames (1000 and 2000 rows) to
  5.3e-15 with identical `sample__`. `conditional_effects(spaghetti =
  TRUE)` equals brms's for `"x"` and `"x:f"`, rows, names and values.
- `posterior_average()`: identical dimensions, names, `ndraws` and
  `weights` attributes, values to 1.4e-17, for three weight vectors, a
  smaller `ndraws`, and `missing = 0` across two different models. A
  repeated model label stays repeated, as brms leaves it (the first
  run made it unique; brms's `bK,bK` against frmtmb's `dsK,dsK.1`
  showed the difference, and it was changed to match).
- `conditional_effects(surface = TRUE, too_far = 0.2)`: 225 rows in
  both, `estimate__` to 1.5e-15.
- `make_conditions()` on 5 inputs and `update_adterms()` on 9:
  `identical()` to brms's every time.

## 4. The nonlinear Wald band

`dev/postfit2-nlmc.R`, log `dev/postfit2-log/nlmc.txt`, the B4 model
`y2 ~ a * exp(b * z), a ~ 1 + f, b ~ 1`, grid `"z:f"`, 100 points:

```
coefficients: beta.a_(Intercept), beta.a_fb, beta.b_(Intercept)
se__ equals sqrt(diag(A V A')): max rel diff 0
1. Jacobian vs central differences (h = 1e-05): max |diff| / max |A| = 2.14e-11
2. V = 1, R = 20000, seed 20260929: se/MCsd min 0.9909 median 1.0014 max 1.0077; |z| max 1.81, mean z -0.02; MC error of one SD 0.0050
2. V = 1, R = 20000, seed 20260930: se/MCsd min 0.9969 median 1.0008 max 1.0033; |z| max 0.65, mean z 0.11; MC error of one SD 0.0050
3. V = 0.01, R = 20000, seed 20260929: se/MCsd min 0.9909 median 1.0014 max 1.0081; |z| max 1.81, mean z 0.01; MC error of one SD 0.0050
3. V = 0.01, R = 20000, seed 20260930: se/MCsd min 0.9972 median 1.0010 max 1.0032; |z| max 0.65, mean z 0.14; MC error of one SD 0.0050
```

1. The Jacobian `frm_lp_basis()` tapes agrees with central differences
   of `frm_linpred()` to 2.1e-11 relative, so the band differentiates
   the right function.
2. At the fit's covariance the median ratio of `se__` to the Monte
   Carlo SD is 1.0014 and 1.0008 for two seeds at R = 20000, whose
   Monte Carlo error on one SD is 0.005; no grid point is 2 errors
   away. Curvature is not visible at this spread.
3. At a hundredth of the covariance, where the delta method is exact,
   the ratios move by at most 0.0004 (the same seeds, so the same
   normal deviates), which says the same thing from the other side:
   there was no curvature to remove.

A dissolved signal, recorded rather than deleted: the first comparison
(`brms-compare.txt` B4, one seed, R = 4000, the grid `"z"`) gave a
median ratio 0.9812, 1.7 Monte Carlo errors low, with 21 of 100 points
beyond 2. The points share their draws, so those 21 are one deviation,
not 21. At R = 20000 and two seeds it is gone.

The estimate on the same fit equals brms's `conditional_effects()` at
the same parameters exactly (max difference 0), for `"z"` and `"z:f"`.

## 5. Decisions, with the construction

- **`spaghetti` on a fit.** brms draws one curve per posterior draw.
  A maximum-likelihood fit has draws only under `band = "boot"`, where
  a refit is one curve, so it takes `spaghetti` there and refuses it
  by name under `"wald"` and `"profile"`. `conditional_smooths()` has
  no `band`, so on a fit it refuses `spaghetti`, `ndraws` and
  `draw_ids` by name and names `frm_sample()`. Row
  `brmsfit-methods:273` therefore stays "cannot transfer", and so do
  167 and 169.
- **`surface = TRUE`** was refused on purpose before ("the display
  draws curves with bands, not a fitted surface"). `too_far` means
  nothing without it, and brms's `conditional_smooths()` defaults to
  it, so it is implemented: the grid, brms's numeric `effect2__`, and a
  `plot()` that draws an image with contours. `spaghetti` together
  with `surface` is refused with brms's own message, where brms raises
  it (first predictor numeric).
- **`select_points` measures predictors, not the response.** brms's
  `make_point_frame()` measures every variable its conditions hold,
  and its conditions hold the response at its mean, so brms also drops
  the observations whose RESPONSE is far from the mean response. That
  hides exactly the scatter the points are there to show, which is
  the tiebreaker's "clearly obvious". B3 of the compare log: at
  `select_points = 0.1` brms keeps 3 rows and frmtmb 5; applying
  brms's response filter to frmtmb's rows gives brms's 3 exactly, and
  so at 0.05 and 0.3 (1 and 45). Factors keep frmtmb's rule, measured
  only when a condition sets them (the documented divergence of
  0.62.0).
- **The nonlinear band** is refused where one Jacobian is not enough:
  an expected response through several predictors (zero-inflated,
  hurdle, `trials()`, truncation), a reported mixing weight, an
  ordinal display, and a new group under `re_formula = NULL`, where
  `frm_lp_basis()` itself refuses `allow_new_levels`. Each names
  `band = "boot"`. The zero-inflated case is in
  `tests/testthat/test-ce-options.R`, with its absent case
  (`dpar = "mu"` draws).
- **`update_adterms()` on frmtmb's formula objects.** brms returns a
  plain formula for a `brmsformula`; frmtmb returns the `bf()` object
  with its response formula updated and every other part kept, since
  dropping the distributional formulas would be a silent loss. An
  `mvbf()` is refused by name.
- **`posterior_average()` weights.** Numeric, `"stacking"`,
  `"pseudobma"`, `"loo"` and `"waic"` follow brms's `model_weights()`.
  `"kfold"` and `"bma"` are refused by name: `kfold()` and
  `post_prob()` on draws already refuse.
- **Class names.** The returned lists are `frmtmb_conditional_effects`,
  per rule 2 of 2026-09-17; row `brmsfit-methods:270` stays a
  divergence.

## 6. The ported rows this change flips

Do not edit `dev/brmsport-*` here (lane `defects` owns them). Run on
this lane's library with `FRMTMB_BRMS_FIT_TESTS=true`, each of these
now HOLDS where the ledger records it as not holding, so the tier
reports it as a stale verdict (logs
`dev/postfit2-log/gated-*-test-brms-suite-*.txt`; the base build's in
`gated-base-*`):

```
this lane's library:
core-test-brms-suite-brmsfit-helpers.R.txt: brms brmsfit-helpers:94 now HOLDS but is recorded as 'cannot transfer'
core-test-brms-suite-brmsfit-helpers.R.txt: brms brmsfit-helpers:95 now HOLDS but is recorded as 'cannot transfer'
core-test-brms-suite-brmsformula.R.txt: brms brmsformula:44 now HOLDS but is recorded as 'cannot transfer'
core-test-brms-suite-brmsformula.R.txt: brms brmsformula:48 now HOLDS but is recorded as 'cannot transfer'
core-test-brms-suite-brmsformula.R.txt: brms brmsformula:52 now HOLDS but is recorded as 'cannot transfer'
core-test-brms-suite-brmsformula.R.txt: brms brmsformula:56 now HOLDS but is recorded as 'cannot transfer'
core-test-brms-suite-methods.R.txt: brms brmsfit-methods:157 now HOLDS but is recorded as 'cannot transfer'
core-test-brms-suite-methods.R.txt: brms brmsfit-methods:170 now HOLDS but is recorded as 'cannot transfer'
core-test-brms-suite-methods.R.txt: brms brmsfit-methods:213 now HOLDS but is recorded as 'defect'
core-test-brms-suite-methods.R.txt: brms brmsfit-methods:215 now HOLDS but is recorded as 'defect'
core-test-brms-suite-methods.R.txt: brms brmsfit-methods:269 now HOLDS but is recorded as 'defect'
core-test-brms-suite-methods.R.txt: brms brmsfit-methods:275 now HOLDS but is recorded as 'defect'
core-test-brms-suite-methods.R.txt: brms brmsfit-methods:277 now HOLDS but is recorded as 'defect'
sample-test-brms-suite-brmsformula.R.txt: brms brmsformula:44 now HOLDS but is recorded as 'cannot transfer'
sample-test-brms-suite-brmsformula.R.txt: brms brmsformula:48 now HOLDS but is recorded as 'cannot transfer'
sample-test-brms-suite-brmsformula.R.txt: brms brmsformula:52 now HOLDS but is recorded as 'cannot transfer'
sample-test-brms-suite-brmsformula.R.txt: brms brmsformula:56 now HOLDS but is recorded as 'cannot transfer'
sample-test-brms-suite-methods.R.txt: brms brmsfit-methods:619 now HOLDS but is recorded as 'cannot transfer'
sample-test-brms-suite-methods.R.txt: brms brmsfit-methods:620 now HOLDS but is recorded as 'cannot transfer'
stale-verdict failures: 19, all of them rows above

the base build, same five files (every row holds its recorded verdict):
RESULT frmtmb test-brms-suite-brmsfit-helpers.R: tests=3 failed=0 error=0 skipped=0 warning=0 passed=3
RESULT frmtmb test-brms-suite-brmsformula.R: tests=16 failed=0 error=0 skipped=0 warning=0 passed=16
RESULT frmtmb test-brms-suite-methods.R: tests=162 failed=0 error=0 skipped=0 warning=0 passed=162
RESULT frmtmb.sample test-brms-suite-brmsformula.R: tests=16 failed=0 error=0 skipped=0 warning=0 passed=16
RESULT frmtmb.sample test-brms-suite-methods.R: tests=61 failed=0 error=0 skipped=0 warning=0 passed=61
```

Rows 627 and 628 of `posterior_average` are "not run" (the setup line
reaches `brms:::SW`) and stay so. 162 and 164 fail on `plot(plot =
FALSE)` and `stype`, not on `too_far`, and 167, 169 and 273 ask an ML
fit for spaghetti.

## 7. Tests run

This is the first round, before section 7b's fix. The final runs of
both suites, the gated tier and `R CMD check` are at the end of 7b.

```
whole ungated suites, one process per file (dev/postfit2-suite.sh):
jobs: 222
logs with RESULT: 222
failed=9 error=5 skipped=156 warning=0 passed=14290
files with failures or errors:
RESULT frmtmb.sample test-message-uniqueness.R: tests=7 failed=2 error=0 skipped=0 warning=0 passed=5
RESULT frmtmb test-ce-bands.R: tests=167 failed=2 error=0 skipped=0 warning=0 passed=165
RESULT frmtmb test-conditions.R: tests=145 failed=4 error=3 skipped=0 warning=0 passed=141
RESULT frmtmb test-data2.R: tests=29 failed=0 error=1 skipped=0 warning=0 passed=29
RESULT frmtmb test-id-kron.R: tests=59 failed=0 error=1 skipped=0 warning=0 passed=59
RESULT frmtmb test-nlf.R: tests=75 failed=1 error=0 skipped=0 warning=0 passed=74
files without a RESULT line:
(222 'Segmentation fault' lines from Rscript at process exit, after each RESULT line; every log has its RESULT)

after the fixes that run prompted, rerun one file per process:
RESULT frmtmb test-conditional-smooths.R: tests=39 failed=0 error=0 skipped=0 warning=0 passed=39
RESULT frmtmb test-brms-utilities.R: tests=28 failed=0 error=0 skipped=0 warning=0 passed=28
RESULT frmtmb test-ce-options.R: tests=29 failed=0 error=0 skipped=0 warning=0 passed=29
RESULT frmtmb test-ce-bands.R: tests=169 failed=0 error=0 skipped=0 warning=0 passed=169
RESULT frmtmb test-nlf.R: tests=77 failed=0 error=0 skipped=0 warning=0 passed=77
RESULT frmtmb test-conditions.R: tests=145 failed=3 error=3 skipped=0 warning=0 passed=142
RESULT frmtmb test-generic-collision.R: tests=56 failed=0 error=0 skipped=0 warning=0 passed=56
RESULT frmtmb test-adefects.R: tests=85 failed=0 error=0 skipped=0 warning=0 passed=85
RESULT frmtmb.sample test-postfit-draws.R: tests=28 failed=0 error=0 skipped=0 warning=0 passed=28
RESULT frmtmb.sample test-message-uniqueness.R: tests=6 failed=0 error=0 skipped=0 warning=0 passed=6
RESULT frmtmb.sample test-conditional-effects-draws.R: tests=61 failed=0 error=0 skipped=0 warning=0 passed=61
RESULT frmtmb.sample test-generic-collision.R: tests=58 failed=0 error=0 skipped=0 warning=0 passed=58
RESULT frmtmb.sample test-adefects.R: tests=20 failed=0 error=0 skipped=0 warning=0 passed=20

gated brms-suite tier, all files, this lane's library:
RESULT frmtmb test-brms-suite-brm.R: tests=23 failed=0 error=0 skipped=0 warning=0 passed=23
RESULT frmtmb test-brms-suite-brmsfit-helpers.R: tests=3 failed=2 error=0 skipped=0 warning=0 passed=1
RESULT frmtmb test-brms-suite-brmsformula.R: tests=16 failed=4 error=0 skipped=0 warning=0 passed=12
RESULT frmtmb test-brms-suite-brmsterms.R: tests=5 failed=0 error=0 skipped=0 warning=0 passed=5
RESULT frmtmb test-brms-suite-data-helpers.R: tests=6 failed=0 error=0 skipped=0 warning=0 passed=6
RESULT frmtmb test-brms-suite-emmeans.R: tests=11 failed=0 error=0 skipped=0 warning=0 passed=11
RESULT frmtmb test-brms-suite-families.R: tests=84 failed=0 error=0 skipped=0 warning=0 passed=84
RESULT frmtmb test-brms-suite-methods.R: tests=162 failed=7 error=0 skipped=0 warning=0 passed=155
RESULT frmtmb test-brms-suite-priors.R: tests=36 failed=0 error=0 skipped=0 warning=0 passed=36
RESULT frmtmb test-brms-suite-standata.R: tests=87 failed=0 error=0 skipped=0 warning=0 passed=87
RESULT frmtmb.sample test-brms-suite-brmsformula.R: tests=16 failed=4 error=0 skipped=0 warning=0 passed=12
RESULT frmtmb.sample test-brms-suite-brmsterms.R: tests=4 failed=0 error=0 skipped=0 warning=0 passed=4
RESULT frmtmb.sample test-brms-suite-families.R: tests=84 failed=0 error=0 skipped=0 warning=0 passed=84
RESULT frmtmb.sample test-brms-suite-helper-copy.R: tests=33 failed=0 error=0 skipped=0 warning=0 passed=33
RESULT frmtmb.sample test-brms-suite-methods.R: tests=61 failed=2 error=0 skipped=0 warning=0 passed=59

R CMD check --as-cran (dev/postfit2-check.ps1):
frmtmb Status: 1 NOTE
frmtmb.sample Status: OK
```

Reading it:

- The whole-suite run found three things of this lane's: two tests in
  `test-ce-bands.R` and one in `test-nlf.R` asserted the old refusals
  (the nonlinear Wald band, `surface = TRUE`), and frmtmb.sample's
  message-uniqueness test caught one duplicated refusal in
  `posterior_average()`. The tests now assert the new behavior; the
  duplicate was merged. All three files pass on the rerun, and the
  updated `test-ce-bands.R` and `test-nlf.R` fail on the base build
  (`before-core-test-ce-bands.R.txt`: 2 errors;
  `before-core-test-nlf.R.txt`: 1 error).
- `test-conditions.R`, `test-data2.R` and `test-id-kron.R` fail the
  same way on the base build under this runner
  (`before-core-test-conditions.R.txt` and siblings): they read
  internal functions that `testthat::test_file()`'s environment does
  not expose here. Under `R CMD check`, which runs the suite in the
  namespace, they pass.
- The 19 gated failures (15 rows; the four `brmsformula` rows fail in
  both packages) are the stale verdicts of section 6, and
  nothing else.
- frmtmb's NOTE is the expected V8 math-rendering note. frmtmb.sample
  was checked twice: the first run gave one NOTE, `arg_desc()` is not
  exported by frmtmb and `posterior_average()` called it. Replaced;
  the second run is `Status: OK`.

## 7b. A silent wrong answer from lane sampfix: an observed group

**WITHDRAWN IN PART: section 7c replaces the fix below.** Review
(`dev/reviews/2026-09-29-postfit2.md`) found that `ce_set_vars()` made
new silent wrong answers (an unseen level read as observed, a mixed
`conditions` column, a nested term, crossed terms under the bootstrap,
and `y ~ s(x)` under `re_formula = NULL`), and that the NEWS claim for
crossed terms was false. `ce_set_vars()` is removed. The diagnosis
below stands; the design and its numbers for the bootstrap do not.

Reported by the coordinator from lane sampfix
(`C:\Users\adf44\source\r\frmtmb-wt-sampfix\dev\sampfix-08-ce-level.R`,
read, not edited). `conditional_effects(re_formula = NULL, conditions =
data.frame(g = <an observed level>))` on draws drew a NEW level's
effects. The cause is where both methods build `na_vars`: every
grouping variable, whether or not `conditions` sets it. So
`ce_new_level_spec()` made a new-level spec for `g`, `ce_boot_grids()`
wrote the placeholder level (level one) over the grid's `g`, and every
draw (every bootstrap replicate, on the fit) drew fresh effects into
it.

The fix, the same line in both methods: `na_vars` is
`setdiff(ce_group_vars(x), ce_set_vars(conditions))`, where the new
core helper `ce_set_vars()` (exported on `?frmtmb-sampling-api`) names
the variables `conditions` sets to a value. A column that is all `NA`
is brms's "unset" and stays a new group. **The hunk in
`extensions/frmtmb.sample/R/conditional-effects-draws.R` is that one
statement** (`na_vars <- if (!pop_level) { setdiff(ce_group_vars(fit),
ce_set_vars(conditions)) } else { character(0) }`, with a two-line
comment), where lane sampfix also edits the file.

The fit method had the same defect under `band = "boot"`, and a second
one under it: its refits simulated responses with FRESH group effects
(`frm_bootstrap(re_formula = NA)`), so a refit's level was a different
group in every replicate and the band sat on the population curve.
`ce_boot_draws()` now simulates conditional on the fitted effects
(`re_formula = NULL`) when the display is at observed groups only (not
the population curve, no new-level spec), and `ce_boot_key()` records
which, so a population bootstrap cannot be reused for such a call.

`conditional_smooths()` has no `conditions` or `re_formula` and is not
reached. The spaghetti path on draws reads the same per-draw matrix and
was wrong the same way; it is fixed with it.

`dev/postfit2-celevel.R` (data seed 9, draws seed 1; the second fit
data seed 21, group SD 2), base build then this lane:

```
BASE
arm base: frmtmb 0.65.0, frmtmb.sample 0.13.0
draws, g = 1: curve moves by 0 0 0 when r_g[1,] moves by 10
draws, g = 2: curve moves by 0 0 0 when r_g[2,] moves by 10
draws: level 1 minus level 2: 0 0 0 (the median of r_g[1] - r_g[2]: 0.24946)
draws, g = 2: spaghetti: conditional_effects() on draws cannot honor `spaghetti`: this method does not implement it, so leave it at its default
fit, g = 2 (mode -0.2075): wald estimate -0.77695  0.60839  1.99374
fit, boot mean of the refit curves -0.58557  0.85491  2.29540; its offset from the estimate over the mode: 1.454
fit, boot band width / wald band width: 0.825
fit, boot band at g = 1 minus at g = 2: lower 0 0 0 (estimates differ by 0.2508 0.2508 0.2508)
fit (group SD 2, seed 21), g = 1, mode 1.983: boot mean minus estimate -1.699 -1.670 -1.642; over the mode: 0.857
  estimate inside the boot band at every point: TRUE; band width / wald width 3.117
LANE
arm lane: frmtmb 0.65.0, frmtmb.sample 0.13.0
draws, g = 1: curve moves by 10 10 10 when r_g[1,] moves by 10
draws, g = 2: curve moves by 10 10 10 when r_g[2,] moves by 10
draws: level 1 minus level 2: 0.311731 0.246607 0.266987 (the median of r_g[1] - r_g[2]: 0.24946)
draws, g = 2: spaghetti moves by 10 to 10
fit, g = 2 (mode -0.2075): wald estimate -0.77695  0.60839  1.99374
fit, boot mean of the refit curves -0.60763  0.75019  2.10802; its offset from the estimate over the mode: 0.816
fit, boot band width / wald band width: 0.682
fit, boot band at g = 1 minus at g = 2: lower 0.16473 0.15903 0.09726 (estimates differ by 0.2508 0.2508 0.2508)
fit (group SD 2, seed 21), g = 1, mode 1.983: boot mean minus estimate -0.17617 -0.13730 -0.09843; over the mode: 0.089
  estimate inside the boot band at every point: TRUE; band width / wald width 0.998
```

On the small-SD fit (mode -0.21) the refits' mean is still 0.82 modes
from the estimate after the fix. That is the ordinary double shrinkage
of a conditional parametric bootstrap at a heavily shrunk level (0.17
in absolute terms), not the defect: at a group SD of 2 it is 0.089.

Against brms at the same five draws (`dev/postfit2-celevel-brms.R`,
brms `fixed_param`, one chain per draw):

```
exp(theta) 0.220984; the VarCorr() sd:
           Estimate Est.Error       Q2.5    Q97.5
Intercept 0.2209842 0.1829048 -0.1375027 0.579471
brms draws 5; r_g[2,Intercept] brms vs frmtmb max diff 0
g = 1: estimate__ max |brms - frmtmb| 0, se__ 0, lower__ 0, upper__ 0; spaghetti rows 35 / 35, estimate__ 0
g = 2: estimate__ max |brms - frmtmb| 0, se__ 0, lower__ 0, upper__ 0; spaghetti rows 35 / 35, estimate__ 0
g = 3: estimate__ max |brms - frmtmb| 0, se__ 0, lower__ 0, upper__ 0; spaghetti rows 35 / 35, estimate__ 0
```

Tests: `tests/testthat/test-ce-options.R` "a boot band at an observed
group belongs to that group" and frmtmb.sample's
`test-postfit-draws.R` "a curve at an observed group reads that
group's draws" (hand-built draws, no sampler). Both fail on the base
build behaviorally: the curve does not move when the level's draws do,
and the two levels' boot bands are identical.

```
base build:
RESULT frmtmb test-ce-options.R: tests=5 failed=4 error=5 skipped=0 warning=0 passed=1   [before-core-test-ce-options.R.txt]
RESULT frmtmb.sample test-postfit-draws.R: tests=1 failed=1 error=5 skipped=0 warning=0 passed=0   [before-sample-test-postfit-draws.R.txt]
RESULT frmtmb test-brms-methods.R: tests=975 failed=0 error=1 skipped=0 warning=0 passed=975   [before-core-test-brms-methods.R.txt]
this lane, after the fix, one process per file:
RESULT frmtmb test-ce-options.R: tests=32 failed=0 error=0 skipped=0 warning=0 passed=32   [after3-core-test-ce-options.R.txt]
RESULT frmtmb test-ce-bands.R: tests=169 failed=0 error=0 skipped=0 warning=0 passed=169   [after3-core-test-ce-bands.R.txt]
RESULT frmtmb test-nlf.R: tests=77 failed=0 error=0 skipped=0 warning=0 passed=77   [after3-core-test-nlf.R.txt]
RESULT frmtmb.sample test-postfit-draws.R: tests=33 failed=0 error=0 skipped=0 warning=0 passed=33   [after3-sample-test-postfit-draws.R.txt]
RESULT frmtmb.sample test-conditional-effects-draws.R: tests=61 failed=0 error=0 skipped=0 warning=0 passed=61   [after3-sample-test-conditional-effects-draws.R.txt]
RESULT frmtmb.sample test-draws-methods.R: tests=148 failed=0 error=0 skipped=0 warning=0 passed=148   [after3-sample-test-draws-methods.R.txt]
RESULT frmtmb.sample test-re-formula-draws.R: tests=22 failed=0 error=0 skipped=0 warning=0 passed=22   [after3-sample-test-re-formula-draws.R.txt]
RESULT frmtmb.sample test-predfix-new-levels.R: tests=62 failed=0 error=0 skipped=0 warning=0 passed=62   [after3-sample-test-predfix-new-levels.R.txt]
the whole ungated suites again, after the fix (dev/postfit2-suite.sh):
jobs: 222
logs with RESULT: 222
failed=4 error=5 skipped=156 warning=0 passed=14306
files with failures or errors:
RESULT frmtmb test-conditions.R: tests=145 failed=4 error=3 skipped=0 warning=0 passed=141
RESULT frmtmb test-data2.R: tests=29 failed=0 error=1 skipped=0 warning=0 passed=29
RESULT frmtmb test-id-kron.R: tests=59 failed=0 error=1 skipped=0 warning=0 passed=59
files without a RESULT line:
the rest of the gated tier (dev/postfit2-gated.sh):
jobs: 23
logs with RESULT: 23
RESULT frmtmb.sample test-loo.R: tests=111 failed=0 error=0 skipped=0 warning=0 passed=111
RESULT frmtmb.sample test-sampling-ported.R: tests=213 failed=0 error=0 skipped=0 warning=0 passed=213
RESULT frmtmb test-bcm-bart.R: tests=12 failed=0 error=0 skipped=0 warning=0 passed=12
RESULT frmtmb test-bcm-binomial.R: tests=37 failed=0 error=0 skipped=0 warning=0 passed=37
RESULT frmtmb test-bcm-data-analysis.R: tests=46 failed=0 error=0 skipped=0 warning=0 passed=46
RESULT frmtmb test-bcm-esp.R: tests=17 failed=0 error=0 skipped=0 warning=0 passed=17
RESULT frmtmb test-bcm-gaussian.R: tests=23 failed=0 error=0 skipped=0 warning=0 passed=23
RESULT frmtmb test-bcm-gcm.R: tests=17 failed=0 error=0 skipped=0 warning=0 passed=17
RESULT frmtmb test-bcm-latent-mixtures.R: tests=52 failed=0 error=0 skipped=0 warning=0 passed=52
RESULT frmtmb test-bcm-model-selection.R: tests=29 failed=0 error=0 skipped=0 warning=0 passed=29
RESULT frmtmb test-bcm-mpt.R: tests=27 failed=0 error=0 skipped=0 warning=0 passed=27
RESULT frmtmb test-bcm-psychophysics.R: tests=28 failed=0 error=0 skipped=0 warning=0 passed=28
RESULT frmtmb test-bcm-retention.R: tests=32 failed=0 error=0 skipped=0 warning=0 passed=32
RESULT frmtmb test-bcm-signal-detection.R: tests=30 failed=0 error=0 skipped=0 warning=0 passed=30
RESULT frmtmb test-bcm-simple.R: tests=13 failed=0 error=0 skipped=0 warning=0 passed=13
RESULT frmtmb test-brms-agreement.R: tests=187 failed=0 error=0 skipped=0 warning=0 passed=187
RESULT frmtmb test-brms-likelihood.R: tests=440 failed=0 error=0 skipped=0 warning=0 passed=440
RESULT frmtmb test-brms-methods.R: tests=983 failed=0 error=0 skipped=0 warning=0 passed=983
RESULT frmtmb test-brms-port.R: tests=21 failed=0 error=0 skipped=0 warning=0 passed=21
RESULT frmtmb test-brms-priors.R: tests=103 failed=0 error=0 skipped=0 warning=0 passed=103
RESULT frmtmb test-drmtmb-agreement.R: tests=131 failed=0 error=0 skipped=0 warning=0 passed=131
RESULT frmtmb test-fuzz.R: tests=2 failed=0 error=0 skipped=0 warning=0 passed=2
RESULT frmtmb test-rl-example.R: tests=102 failed=0 error=0 skipped=0 warning=0 passed=102
files with a skip, a failure or an error:
  none
(test-loo.R, test-sampling-ported.R and test-brms-agreement.R are the rerun after the runner set TMP; their first logs errored on a Stan compile writing to C:\WINDOWS)
R CMD check --as-cran, final (dev/postfit2-check.ps1):
frmtmb Status: 1 NOTE
frmtmb.sample Status: OK
```

The two brms-suite files that reach `conditional_effects()` with
`re_formula = NULL` were rerun after the fix
(`dev/postfit2-log/gated2-*.txt`): core `test-brms-suite-methods.R`
162 tests, 7 failed, frmtmb.sample's 61 tests, 2 failed, and the
failures are the same rows as in section 6 (157, 170, 213, 215, 269,
275, 277; 619, 620). The ungated suites' 156 skips are the gated
blocks, which the gated tier above runs.

An instrument defect, recorded because it cost two runs: the first two
runs of `dev/postfit2-gated.sh` skipped every gated block (155
expectations; the first is kept in `dev/postfit2-gated-run1/`). The
switches were exported by the driver and did not reach R under
`xargs`, and the skip column said so while the failure column was
clean. The runner now takes `--gated` and sets them itself; the block
above is the third run, and its skip column is what to read.

Not fixed: a `conditions` data frame whose `g` column mixes levels and
`NA` rows. `g` counts as set, so its `NA` rows get no new-level draw;
before the fix its observed rows were overwritten instead. Neither is
right, and the case needs a per-row placeholder.

## 7c. Punch round 1: brms's per-row, per-term rule, and the bootstrap

Review `dev/reviews/2026-09-29-postfit2.md` found the fix of 7b wrong
on adjacent shapes (B1 to B7, M1), and the user decided on 2026-09-29
that no bootstrap redraws a smooth. The design was reworked rather
than patched.

**The rule, brms's.** brms predicts its display with
`allow_new_levels = TRUE`. Per grid row and per group-level term, the
term reads its fitted effects where every one of its grouping
variables is set in that row to a level the fit saw for that term, and
a NEW level everywhere else. `R/ce-levels.R` (new) implements it once
for both methods:

- `ce_level_plan()` reads, per row and per term, whether the row's
  level string is one of the term's fitted levels. A variable set in
  `conditions` or varied as an effect (M1) counts as set; an unset
  variable, `NA`, an unseen level (`"99"`, kept in the frame, m2) and a
  combination a nested term never saw are new. Rows with the same
  pattern of new levels are predicted together, at a placeholder
  level of each new term whose coefficients the draw overwrites. A
  placeholder moves only grouping columns that nothing else reads
  (`ce_locked_vars()`) and no observed term of those rows reads, so a
  nested `(1 | g / h)` row with `g` set and `h` unset keeps `g`'s level
  and draws a new `g:h` within it (B3). A pattern no placeholder can
  serve is refused by name.
- `ce_plan_eval()` evaluates a grid at one parameter vector, drawing
  each distinct new level once from that vector's covariance and
  sharing it between rows and panels.
- Both methods predict with `allow_new_levels = TRUE` where the
  display is not the population curve, so B1's and B2's "New levels
  ... NA" refusals are gone rather than reworded.
- The Wald band needed no change: `lp_extra_var()` already adds a new
  level's variance per row and per term.
- `ce_set_vars()` is removed. `ce_new_level_spec()`,
  `ce_boot_grids()` and `ce_draw_new_levels()` stay exported and
  unused. (Removed in punch round 2, with `ce_new_level_key()`.)
- "The Wald band needed no change" held for one-column groupings only;
  see the correction in section 9 on `mm()`.

**The fit's bootstrap.** A term some grid row reads at an observed
level is held at its fitted effects in the simulation; every other
group-level term is redrawn; a term read observed in some rows and new
in others is refused by name (B2, and the coordinator's allowance).
`ce_boot_draws()` passes that set to `frm_bootstrap()`'s body
(`boot_refits()`, new) as a `sim_re_plan()` list, since no single
`re_formula` spells a set read off the grid. Only group-level terms
enter the plan, so a model with none reads the same band under `NA`
and `NULL` (B5). The reuse key carries the held set (B7), and keys are
compared by value, because the unseen-level extension builds a
factor's levels with `union()` and their serialized bytes differ from
an equal factor read off the data (`dev/postfit2-p1-keydbg.R`; without
this the keydrop mutant was caught by the byte difference and the test
could not pin the field).

**B6.** The nonlinear guard asks `frm_lp_basis()` of the displayed
predictor alone and turns its classed new-levels refusal into the
band's own; `sigma ~ (1 | g)` beside a nonlinear `mu` no longer
refuses.

**frm_bootstrap() never redraws a smooth (user decision).** It now
reads `re_formula` exactly as `predict()` and `simulate()` do: a
smooth, `gp()` or `hsgp()` term, one indexed by a grouping factor
included, is held in every replicate and re-estimated by the refit;
`NA`, `~0` and `~1` redraw every true group-level effect. The
`redraw_smooths` argument of `sim_fit_draws()` and the `smooths`
argument of `sim_re_plan()` are removed; nothing else used them
(grep over core and the extensions). The CE bootstrap inherits the
rule. The test "frm_bootstrap() redraws a smooth by default, and NULL
holds it" is inverted: on `y ~ s(x)` `NA` and `NULL` give identical
replicates at the same seed, and beside `(1 | g)` `NA` still redraws
`g`. The `dev/simnewdata-boot.R` setup (`y ~ s(x)`, 200 rows, seed
21, 40 refits), by `dev/postfit2-p1-simboot.R`:

```
frmtmb from C:/Users/adf44/source/r/rellib-r3 
base: re_formula NA, estimate -0.2833 / 5.2011; bootstrap sd 0.0385 / 1.4376; converged 40 of 40
base: re_formula NULL, estimate -0.2833 / 5.2011; bootstrap sd 0.0432 / 0.3693; converged 40 of 40
base: intercept Wald se 0.0460
frmtmb from C:/Users/adf44/source/r/wt-postfit2-lib 
lane: re_formula NA, estimate -0.2833 / 5.2011; bootstrap sd 0.0432 / 0.3693; converged 40 of 40
lane: re_formula NULL, estimate -0.2833 / 5.2011; bootstrap sd 0.0432 / 0.3693; converged 40 of 40
lane: intercept Wald se 0.0460
```

**The round-0 lane build.** To show each new test failing on the
build the review ran, the tarballs of the round-0 `R CMD check`
(`dev/postfit2-check/*/`, copied before this round's check replaced
them) were installed into a second private library,
`C:/Users/adf44/source/r/wt-postfit2-r0lib`, used by
`dev/postfit2-runtest.R --r0` only.

Seen to fail, then pass (`dev/postfit2-log/p1-{r0,base,lane}-*`):

```
## core-test-ce-levels.R.txt
r0: RESULT frmtmb test-ce-levels.R: tests=25 failed=11 error=4 skipped=0 warning=0 passed=14   [p1-r0-core-test-ce-levels.R.txt]
  fails: a level the fit never saw is a new group, and keeps its name 
  fails: rows mixing a level and NA read one each, or boot says why not 
  fails: a nested term set at its parent draws a new level within it 
  fails: crossed terms, one set: the boot band sits on its estimate 
  fails: re_formula changes nothing on a model with no group term 
  fails: a group only another parameter reads does not refuse the band 
  fails: a population bootstrap is not reused for an observed group 
  fails: a grouping variable varied as an effect is an observed group 
  fails: the bootstrap band holds the smooths at their fit 
base: RESULT frmtmb test-ce-levels.R: tests=25 failed=9 error=4 skipped=0 warning=0 passed=16   [p1-base-core-test-ce-levels.R.txt]
  fails: a level the fit never saw is a new group, and keeps its name 
  fails: rows mixing a level and NA read one each, or boot says why not 
  fails: crossed terms, one set: the boot band sits on its estimate 
  fails: a group only another parameter reads does not refuse the band 
  fails: a population bootstrap is not reused for an observed group 
  fails: the bootstrap band holds the smooths at their fit 
lane: RESULT frmtmb test-ce-levels.R: tests=36 failed=0 error=0 skipped=0 warning=0 passed=36   [p1-lane-core-test-ce-levels.R.txt]
## core-test-simulate-newdata.R.txt
r0: RESULT frmtmb test-simulate-newdata.R: tests=57 failed=1 error=0 skipped=0 warning=0 passed=56   [p1-r0-core-test-simulate-newdata.R.txt]
  fails: frm_bootstrap() never redraws a smooth; NA still redraws a group 
base: RESULT frmtmb test-simulate-newdata.R: tests=57 failed=1 error=0 skipped=0 warning=0 passed=56   [p1-base-core-test-simulate-newdata.R.txt]
  fails: frm_bootstrap() never redraws a smooth; NA still redraws a group 
lane: RESULT frmtmb test-simulate-newdata.R: tests=57 failed=0 error=0 skipped=0 warning=0 passed=57   [p1-lane-core-test-simulate-newdata.R.txt]
## sample-test-postfit-draws.R.txt
r0: RESULT frmtmb.sample test-postfit-draws.R: tests=37 failed=3 error=1 skipped=0 warning=0 passed=34   [p1-r0-sample-test-postfit-draws.R.txt]
  fails: on draws, an unseen or unset row is a new group, per row 
  fails: on draws, a nested term set at its parent reads the parent 
  fails: on draws, a grouping variable varied as an effect is observed 
base: RESULT frmtmb.sample test-postfit-draws.R: tests=10 failed=6 error=5 skipped=0 warning=0 passed=4   [p1-base-sample-test-postfit-draws.R.txt]
  fails: conditional_smooths() on draws summarizes the per-draw terms 
  fails: conditional_effects() on draws takes spaghetti and surface 
  fails: posterior_average() takes draws in proportion to the weights 
  fails: posterior_average()'s stacking weights are loo's 
  fails: a curve at an observed group reads that group's draws 
  fails: on draws, an unseen or unset row is a new group, per row 
  fails: on draws, a nested term set at its parent reads the parent 
  fails: on draws, a grouping variable varied as an effect is observed 
lane: RESULT frmtmb.sample test-postfit-draws.R: tests=42 failed=0 error=0 skipped=0 warning=0 passed=42   [p1-lane-sample-test-postfit-draws.R.txt]
```

B7, the reuse key's held set, pinned by a mutant that drops it from
the key (`dev/postfit2-p1-mutant.R`):

```
RESULT keydrop-run none: failed=0 error=0 passed=21
[1] test   failed error 
<0 rows> (or 0-length row.names)
RESULT keydrop-run keydrop: failed=2 error=0 passed=19
                                                        test failed error
7 a population bootstrap is not reused for an observed group      2 FALSE
```

The review's own scripts, rerun on this build
(`dev/postfit2-log/p1-rev-*-lane.txt`), the lines the blockers are
about:

```
A1 list(g='99') anl=FALSE wald: est -0.5572 0.5889 1.7351 | lo -3.0134 -1.8166 -0.7105 | hi 1.8989 2.9944 4.1807
A1 list(g='99') anl=FALSE boot40 seed5: est -0.5572 0.5889 1.7351 | lo -2.6441 -1.4445 -0.3221 | hi 1.9172 2.7773 3.6961
A0: no conditions (a new group), for the width a new level has
A0 wald: est -0.5572 0.5889 1.7351 | lo -3.0134 -1.8166 -0.7105 | hi 1.8989 2.9944 4.1807
A0 boot40 seed5: est -0.5572 0.5889 1.7351 | lo -2.6441 -1.4445 -0.3221 | hi 1.9172 2.7773 3.6961
A2: rows mixing a level and NA
A2 wald cond lev2: est 0.5409 1.6871 2.8332 | lo -0.2093 1.0906 2.0655 | hi 1.2911 2.2835 3.6010
A2 wald cond unset: est -0.5572 0.5889 1.7351 | lo -3.0134 -1.8166 -0.7105 | hi 1.8989 2.9944 4.1807
A2 boot: frmtmb_error: band = "boot" cannot cover this call: the group-level term (1 | g) is read at an observed level in some grid rows and at a new level in others. The observed rows need refits that keep that term's fitted effects and th
A2 ref: lev2 alone wald: est 0.5409 1.6871 2.8332 | lo -0.2093 1.0906 2.0655 | hi 1.2911 2.2835 3.6010
B modes g: -1.2391 -0.7900 1.0980 -0.0818 -0.0565 -0.9680 1.9902 0.0472; level 7
B g set, h unset: wald: est 0.8967 2.5728 4.2489 | lo -0.2555 1.4874 3.0310 | hi 2.0490 3.6582 5.4668
B g set, h unset: boot60 seed5: est 0.8967 2.5728 4.2489 | lo -0.2962 1.2916 2.6669 | hi 2.5439 3.8464 5.5149
B boot mean minus wald estimate: -0.0476 -0.1364 -0.2251; over |mode| 0.113
B boot width / wald width: 1.2325 1.1769 1.1693
B both set: wald: est 0.7701 2.4462 4.1223 | lo 0.0185 1.8277 3.3131 | hi 1.5218 3.0647 4.9315
B both set: boot60 seed5: est 0.7701 2.4462 4.1223 | lo -0.0165 1.5842 3.0288 | hi 1.3178 2.9097 4.7356
B both set: boot mean minus wald estimate: -0.0947 -0.1282 -0.1616; over |mode| 0.081
B both set: boot width / wald width: 0.8876 1.0716 1.0546
D s(x) re_formula=NA boot60 seed5: est 0.0726 0.9459 -0.0702 -1.0271 -0.2963 | lo -0.0498 0.8490 -0.1561 -1.1492 -0.5204 | hi 0.3092 1.0868 0.0443 -0.8824 -0.0597
D s(x) re_formula=NULL boot60 seed5: est 0.0726 0.9459 -0.0702 -1.0271 -0.2963 | lo -0.0498 0.8490 -0.1561 -1.1492 -0.5204 | hi 0.3092 1.0868 0.0443 -0.8824 -0.0597
D width NULL / NA: 1.0000 1.0000 1.0000 1.0000 1.0000
D s(x) wald: est 0.0726 0.9459 -0.0702 -1.0271 -0.2963 | lo -0.1843 0.8348 -0.1906 -1.1565 -0.5436 | hi 0.3294 1.0570 0.0502 -0.8978 -0.0490
D width NULL-boot / wald: 0.6988 1.0704 0.8323 1.0314 0.9316; NA-boot / wald: 0.6988 1.0704 0.8323 1.0314 0.9316
E pop grid g column: 1,1,1; obs grid g column: 1,1,1
E pop boot into observed-level call: frmtmb_error: boot = was not produced by a conditional_effects(band = "boot") call on this grid: its draws are refits of a different simulation: this call reads observed groups, whose refits keep the fit
E observed boot into pop call (with conditions g=1): frmtmb_error: boot = was not produced by a conditional_effects(band = "boot") call on this grid: its draws are refits of a different simulation: those refits kept the fitted effects of ob
E observed boot into pop call (no conditions): frmtmb_error: boot = was not produced by a conditional_effects(band = "boot") call on this grid: its draws are refits of a different simulation: those refits kept the fitted effects of observed
E same-call reuse: ACCEPTED, identical band TRUE
E new-level boot into observed-level call: frmtmb_error: boot = was not produced by a conditional_effects(band = "boot") call on this grid: its draws are predictions for a different group: those draws carry a NEW group's effects (they came 
C-nested VarCorr g:h sd 0.9308
C-nested VarCorr g sd 0.9657
C-nested VarCorr residual__ sd 0.9581
C-nested sdr of theta: 0.1756 0.3229
C-nested g=2 set, g:h new: wald est 1.4097 2.5634 3.7172 lo -0.7405 0.4449 1.5549 hi 3.5600 4.6820 5.8794
C-nested nothing set: wald lo -3.1938 -2.0147 -0.8940 hi 2.3290 3.4573 4.6441
C-nested g=2 set, g:h new: boot60 seed5 lo -1.2822 -0.0460 1.0265 hi 2.7794 4.0819 5.4752; width/wald 0.9444 0.9742 1.0287
C-nested g and h both set (observed g:h): wald lo 0.3416 1.6098 2.6858 hi 2.1472 3.1865 4.4179
nl, group in nlpar a, observed g=3, re_formula NULL, wald: ok
nl, group in nlpar a, new group (no conditions), wald: REFUSED: conditional_effects() has no Wald band for a nonlinear predictor at a new group: New levels in grouping factor `g`: NA. A new level contributes its block's marginal variance, a
nl mu with NO group term; group only in sigma; re_formula NULL: ok
unseen g='99' anl=FALSE [cond 1]: est -0.5515 0.8041 2.2142 | lo -1.1382 0.3738 1.6329 | hi -0.0197 1.3512 2.8354 | width 1.1185 0.9774 1.2025
unseen g='99' anl=TRUE [cond 1]: est -0.5515 0.8041 2.2142 | lo -1.1382 0.3738 1.6329 | hi -0.0197 1.3512 2.8354 | width 1.1185 0.9774 1.2025
mixed rows g=2 / NA [cond lev2]: est -0.7834 0.6107 2.0002 | lo -1.0676 0.4426 1.6699 | hi -0.5096 0.7609 2.2813 | width 0.5580 0.3183 0.6114
mixed rows g=2 / NA [cond unset]: est -0.5515 0.8041 2.2142 | lo -1.1382 0.3738 1.6329 | hi -0.0197 1.3512 2.8354 | width 1.1185 0.9774 1.2025
mixed rows g=2 / NA, anl [cond lev2]: est -0.7834 0.6107 2.0002 | lo -1.0676 0.4426 1.6699 | hi -0.5096 0.7609 2.2813 | width 0.5580 0.3183 0.6114
mixed rows g=2 / NA, anl [cond unset]: est -0.5515 0.8041 2.2142 | lo -1.1382 0.3738 1.6329 | hi -0.0197 1.3512 2.8354 | width 1.1185 0.9774 1.2025
crossed g=7 set, h unset [cond 1]: est 0.8845 2.5466 4.2338 | lo -0.0762 1.6847 3.2728 | hi 1.9930 3.6401 5.3573 | width 2.0692 1.9554 2.0845
crossed g=7, h=1 both set [cond 1]: est 0.7736 2.4531 4.1349 | lo 0.5038 2.2803 3.7230 | hi 1.0506 2.6193 4.4637 | width 0.5469 0.3390 0.7407
crossed neither set [cond 1]: est -1.0578 0.6420 2.3602 | lo -3.2881 -1.5129 0.1511 | hi 1.0965 2.6666 4.3735 | width 4.3846 4.1795 4.2224
crossed g=7 set: curve moves by 10.0000 10.0000 10.0000 when r_g[7,] moves by 10
nested g=5 set [cond 1]: est 0.4336 1.5356 2.6783 | lo -1.4371 -0.3240 0.8725 | hi 2.7305 3.8988 4.9981 | width 4.1675 4.2227 4.1256
nested g=5, h=1 set [cond 1]: est 0.4336 1.5356 2.6783 | lo -1.4371 -0.3240 0.8725 | hi 2.7305 3.8988 4.9981 | width 4.1675 4.2227 4.1256
nested neither set [cond 1]: est -0.9212 0.2292 1.4110 | lo -3.6150 -2.3978 -1.1805 | hi 1.8249 2.9174 3.9175 | width 5.4399 5.3152 5.0980
1 crossed: brms r_g / r_h vs frmtmb: 0 0
1a both set g=7,h=4: rows 5/5; max |brms - frmtmb| estimate 0 se 0 lower 0 upper 2.22e-16 (scale 4.55)
1b g=7 set, h unset (h pinned at 0): rows 5/5; max |brms - frmtmb| estimate 1.6e-08 se 2.09e-08 lower 8.15e-09 upper 1.47e-08 (scale 4.36)
1c x:g effect, h unset (h pinned): rows 40/40; max |brms - frmtmb| estimate 1.6e-08 se 3.6e-08 lower 8.15e-09 upper 1.47e-08 (scale 4.36)
1d re_formula ~(1|g), g=7: rows 5/5; max |brms - frmtmb| estimate 0 se 3.33e-16 lower 0 upper 2.22e-16 (scale 4.36)
```

- B1: `"99"` now draws the new-group band, 4.56 4.22 4.02 against a
  Wald 4.91 4.81 4.89, where round 0 gave the population band.
- B2: the unset row of a mixed column gets the new-group band on the
  Wald band and on draws; the bootstrap refuses by name.
- B3: nested, `g` set and `h` unset, boot over Wald width 0.94 0.97
  1.03, where round 0 gave 0.25 0.22 0.29.
- B4: crossed, `g` set and `h` unset, the refits' mean is 0.113 modes
  from the estimate and the band contains it (round 0: 1.035 modes, and
  it missed at 2 of 3 points). Width 1.17 to 1.23 times the Wald band.
- B5: `y ~ s(x)`: `NA` and `NULL` identical, 0.70 to 1.07 times the
  Wald band.
- B6: the nonlinear `mu` beside `sigma ~ (1 | g)` answers.
- M1: `effects = "x:g"` equals brms at the same draws to 1.6e-8
  (`p1-rev-brms-lane.txt` 1c), where round 0 differed by 1.76.
- m1: a population boot reused for an observed-level call says
  "refits of a different simulation", not "a different grid".

Against brms at the same draws, with brms's `sample_new_levels =
"gaussian"` (the law frmtmb draws from); a new level is random in both,
so the new-level rows are compared by width (`dev/postfit2-p1-brms.R`,
400 draws; 8 repeats of each for the averaged line):

```
1: draws brms 400; r_g[2] max diff 2.78e-17
1 observed g=2 (no draw): estimate 0 lower 0 upper 0
1 population, for scale: width brms 0.5834 0.1938 0.5730 | frmtmb 0.5834 0.1938 0.5730 | ratio 1.0000 1.0000 1.0000
1 unset g (no conditions): width brms 1.0319 0.8798 1.0391 | frmtmb 1.1486 1.0414 1.1995 | ratio 1.1131 1.1837 1.1544
1 unseen g = '99': width brms 1.0725 0.8923 1.0012 | frmtmb 1.1486 1.0414 1.1995 | ratio 1.0709 1.1672 1.1982
1 mixed rows, lev2 (no draw): upper max diff 0
1 mixed rows, unset: width brms 1.0565 0.9090 1.0912 | frmtmb 1.1486 1.0414 1.1995 | ratio 1.0871 1.1457 1.0993
1 unset g, mean width over 8: brms 1.0053 0.8686 1.0351 (sd 0.0533 0.0360 0.0433) | frmtmb 1.0552 0.9191 1.0805 (sd 0.0535 0.0564 0.0720) | ratio 1.0497 1.0581 1.0438
2: brms groups in order: g, g:h
2: r_g:h[5_2] max diff 1.39e-17
2 g=5, h=2 both observed (no draw): estimate 5.55e-17 upper 0
2 g=5 set, h unset (new g:h within g = 5): width brms 4.0525 3.9013 3.9413 | frmtmb 4.1570 4.1412 4.2506 | ratio 1.0258 1.0615 1.0785
2 g=5 set, h unset: estimate brms -0.3751 0.7698 1.9053 | frmtmb -0.3780 0.7780 1.9663
done
```

Where nothing is drawn (an observed level, the lev2 row of a mixed
column, the nested term with both set) the values are equal to
rounding. The averaged new-level width is 1.04 to 1.06 times brms's,
with a standard error of that ratio of about 0.03 from the eight
repeats; about two standard errors, not resolved further. The nested
new `g:h` width is 1.03 to 1.08 times brms's in one realization.

CORRECTED in punch round 2: the 1.04 to 1.06 is noise from unmatched
seeds, not a difference. The review's `dev/postfit2-rev-p1-law.R`
calls both packages with the same seed (`set.seed(r)` before brms,
`seed = r` for frmtmb): the new-level effects agree to 8.9e-16 in z,
because frmtmb and brms draw from the same law with the same random
numbers, and the width ratio is 1.0000. Section 7d repeats the check
on `gr(g, by = f)` and `mm()`, where the bands agree to 8.9e-16.

The user decision's check on smooths (`dev/postfit2-p1-smooth.R`, 60
refits), boot width over Wald width, base build then this lane:

```
ARM base: frmtmb from C:/Users/adf44/source/r/rellib-r3/frmtmb
s(x), re_formula NA: boot/wald width 10.719 16.116 10.207 17.818 15.903; estimate inside the boot band TRUE
s(x), re_formula NULL: boot/wald width 10.719 16.116 10.207 17.818 15.903; estimate inside the boot band TRUE
s(x) + (1 | g), re_formula NA: boot/wald width 5.114 2.430 2.001 2.499 3.593; estimate inside the boot band FALSE
s(x) + (1 | g), re_formula NULL (new g): boot/wald width 1.539 1.301 0.952 1.137 1.386; estimate inside the boot band FALSE
s(x) + (1 | g), observed g = 3: boot/wald width 14.459 17.603 14.206 15.089 11.960; estimate inside the boot band FALSE
ARM lane: frmtmb from C:/Users/adf44/source/r/wt-postfit2-lib/frmtmb
s(x), re_formula NA: boot/wald width 0.840 0.941 0.913 1.029 0.730; estimate inside the boot band TRUE
s(x), re_formula NULL: boot/wald width 0.840 0.941 0.913 1.029 0.730; estimate inside the boot band TRUE
s(x) + (1 | g), re_formula NA: boot/wald width 0.910 0.961 0.936 0.994 0.820; estimate inside the boot band TRUE
s(x) + (1 | g), re_formula NULL (new g): boot/wald width 0.759 0.749 0.759 0.751 0.767; estimate inside the boot band TRUE
s(x) + (1 | g), observed g = 3: boot/wald width 0.810 0.812 0.871 1.032 0.953; estimate inside the boot band TRUE
```

All tiers after the rework, and the check:

```
core and frmtmb.sample, whole ungated suites (dev/postfit2-suite.sh):
jobs: 223
logs with RESULT: 223
failed=4 error=5 skipped=156 warning=0 passed=14352
files with failures or errors:
RESULT frmtmb test-conditions.R: tests=145 failed=4 error=3 skipped=0 warning=0 passed=141
RESULT frmtmb test-data2.R: tests=29 failed=0 error=1 skipped=0 warning=0 passed=29
RESULT frmtmb test-id-kron.R: tests=59 failed=0 error=1 skipped=0 warning=0 passed=59
files without a RESULT line:

gated tier (dev/postfit2-gated.sh):
jobs: 23
logs with RESULT: 23
RESULT frmtmb.sample test-loo.R: tests=111 failed=0 error=0 skipped=0 warning=0 passed=111
RESULT frmtmb.sample test-sampling-ported.R: tests=213 failed=0 error=0 skipped=0 warning=0 passed=213
RESULT frmtmb test-bcm-bart.R: tests=12 failed=0 error=0 skipped=0 warning=0 passed=12
RESULT frmtmb test-bcm-binomial.R: tests=37 failed=0 error=0 skipped=0 warning=0 passed=37
RESULT frmtmb test-bcm-data-analysis.R: tests=46 failed=0 error=0 skipped=0 warning=0 passed=46
RESULT frmtmb test-bcm-esp.R: tests=17 failed=0 error=0 skipped=0 warning=0 passed=17
RESULT frmtmb test-bcm-gaussian.R: tests=23 failed=0 error=0 skipped=0 warning=0 passed=23
RESULT frmtmb test-bcm-gcm.R: tests=17 failed=0 error=0 skipped=0 warning=0 passed=17
RESULT frmtmb test-bcm-latent-mixtures.R: tests=52 failed=0 error=0 skipped=0 warning=0 passed=52
RESULT frmtmb test-bcm-model-selection.R: tests=29 failed=0 error=0 skipped=0 warning=0 passed=29
RESULT frmtmb test-bcm-mpt.R: tests=27 failed=0 error=0 skipped=0 warning=0 passed=27
RESULT frmtmb test-bcm-psychophysics.R: tests=28 failed=0 error=0 skipped=0 warning=0 passed=28
RESULT frmtmb test-bcm-retention.R: tests=32 failed=0 error=0 skipped=0 warning=0 passed=32
RESULT frmtmb test-bcm-signal-detection.R: tests=30 failed=0 error=0 skipped=0 warning=0 passed=30
RESULT frmtmb test-bcm-simple.R: tests=13 failed=0 error=0 skipped=0 warning=0 passed=13
RESULT frmtmb test-brms-agreement.R: tests=187 failed=0 error=0 skipped=0 warning=0 passed=187
RESULT frmtmb test-brms-likelihood.R: tests=440 failed=0 error=0 skipped=0 warning=0 passed=440
RESULT frmtmb test-brms-methods.R: tests=983 failed=0 error=0 skipped=0 warning=0 passed=983
RESULT frmtmb test-brms-port.R: tests=21 failed=0 error=0 skipped=0 warning=0 passed=21
RESULT frmtmb test-brms-priors.R: tests=103 failed=0 error=0 skipped=0 warning=0 passed=103
RESULT frmtmb test-drmtmb-agreement.R: tests=131 failed=0 error=0 skipped=0 warning=0 passed=131
RESULT frmtmb test-fuzz.R: tests=2 failed=0 error=0 skipped=0 warning=0 passed=2
RESULT frmtmb test-rl-example.R: tests=102 failed=0 error=0 skipped=0 warning=0 passed=102
files with a skip:

frmtmb.sample, whole suite gated (dev/postfit2-sample-suite.sh):
jobs: 37
logs with RESULT: 37
failed=6 error=0 skipped=1 warning=0 passed=2195
files with a failure, an error or a skip:
RESULT frmtmb.sample test-brms-suite-brmsformula.R: tests=16 failed=4 error=0 skipped=0 warning=0 passed=12
RESULT frmtmb.sample test-brms-suite-methods.R: tests=61 failed=2 error=0 skipped=0 warning=0 passed=59
RESULT frmtmb.sample test-scale.R: tests=1 failed=0 error=0 skipped=1 warning=0 passed=0
files without a RESULT line:

core brms-suite tier, gated:
RESULT frmtmb test-brms-suite-brm.R: tests=23 failed=0 error=0 skipped=0 warning=0 passed=23   [p1-gated-core-test-brms-suite-brm.R.txt]
RESULT frmtmb test-brms-suite-brmsfit-helpers.R: tests=3 failed=2 error=0 skipped=0 warning=0 passed=1   [p1-gated-core-test-brms-suite-brmsfit-helpers.R.txt]
RESULT frmtmb test-brms-suite-brmsformula.R: tests=16 failed=4 error=0 skipped=0 warning=0 passed=12   [p1-gated-core-test-brms-suite-brmsformula.R.txt]
RESULT frmtmb test-brms-suite-brmsterms.R: tests=5 failed=0 error=0 skipped=0 warning=0 passed=5   [p1-gated-core-test-brms-suite-brmsterms.R.txt]
RESULT frmtmb test-brms-suite-data-helpers.R: tests=6 failed=0 error=0 skipped=0 warning=0 passed=6   [p1-gated-core-test-brms-suite-data-helpers.R.txt]
RESULT frmtmb test-brms-suite-emmeans.R: tests=11 failed=0 error=0 skipped=0 warning=0 passed=11   [p1-gated-core-test-brms-suite-emmeans.R.txt]
RESULT frmtmb test-brms-suite-families.R: tests=84 failed=0 error=0 skipped=0 warning=0 passed=84   [p1-gated-core-test-brms-suite-families.R.txt]
RESULT frmtmb test-brms-suite-methods.R: tests=162 failed=7 error=0 skipped=0 warning=0 passed=155   [p1-gated-core-test-brms-suite-methods.R.txt]
RESULT frmtmb test-brms-suite-priors.R: tests=36 failed=0 error=0 skipped=0 warning=0 passed=36   [p1-gated-core-test-brms-suite-priors.R.txt]
RESULT frmtmb test-brms-suite-standata.R: tests=87 failed=0 error=0 skipped=0 warning=0 passed=87   [p1-gated-core-test-brms-suite-standata.R.txt]

extension files calling frm_bootstrap(), simulate() or refit(),
this lane's core (dev/postfit2-ext-suite.sh lane):
jobs: 38
logs with RESULT: 38
failed=0 error=0 skipped=10 warning=4 passed=2304
files with a failure, an error or a skip:
RESULT frmtmb.coupling test-scale.R: tests=5 failed=0 error=0 skipped=5 warning=0 passed=0
RESULT frmtmb.eam test-scale.R: tests=3 failed=0 error=0 skipped=3 warning=0 passed=0
RESULT frmtmb.learn test-scale.R: tests=2 failed=0 error=0 skipped=2 warning=0 passed=0
files without a RESULT line:

the same files on the base build (dev/postfit2-ext-suite.sh base):
jobs: 38
logs with RESULT: 38
failed=0 error=0 skipped=10 warning=4 passed=2304
files with a failure, an error or a skip:
RESULT frmtmb.coupling test-scale.R: tests=5 failed=0 error=0 skipped=5 warning=0 passed=0
RESULT frmtmb.eam test-scale.R: tests=3 failed=0 error=0 skipped=3 warning=0 passed=0
RESULT frmtmb.learn test-scale.R: tests=2 failed=0 error=0 skipped=2 warning=0 passed=0
files without a RESULT line:

R CMD check --as-cran (dev/postfit2-check.ps1):
frmtmb Status: 1 NOTE
frmtmb.sample Status: OK
```

## 7d. Punch round 2

The re-check of round 1 found one more blocker on a term shape
(P1-B1), an untested lock (P1-B2), three majors and minors. Each fix
has a test that fails on the round-1 build. That build (the packages
installed at the end of round 1) was copied, before round 2 was
installed, into `C:/Users/adf44/source/r/wt-postfit2-r1lib`, which
`dev/postfit2-runtest.R --r1` and `dev/postfit2-p2-checks.R r1` use.

**P1-B1, `gr(g, by = f)`.** Such a term is one block per level of
`f`, and round 1 read every block in every row. A row with `f = "a"`
and `g = "3"` then read the `f = "b"` block at level `"3"`, which it
does not have, as a NEW level whose placeholder would move `g`; that
was refused, so `band = "boot"` and draws refused every
`re_formula = NULL` display, at observed levels too. `ce_by_reads()`
(new) evaluates the block's `by` expression on the grid, and a row
reads only the block of its own level; a block no row reads is
neither held nor drawn. With `g` unset or unseen, the new level is
drawn from its own `f` level's covariance. `ce_locked_vars()` now
also locks the `by` variables, so no placeholder moves a row into
another by-level's block.

brms at the same draws (`dev/postfit2-p2-brms.R`, 200 draws, seeds
1 to 3): at observed levels the values are equal, and at `g = "99"`
the estimate and band equal brms's `sample_new_levels = "gaussian"`
at the same seed to 8.9e-16, for `f = "a"` and `f = "b"`. With `g`
unset brms stops ("The following variables are missing in the draws
object: {'sd_g__Intercept:f'}"): its new level `NA` has no by-level.
frmtmb answers as for `"99"`, which is what brms does for an unseen
level of that `f`.

**P1-B2, the lock.** A test pins `ce_locked_vars()` on the review's
`y ~ x + trt + (1 | trt:subj)` (data seed 41, `trt = "b"`, `subj`
unset): the refits' mean must be nearer the `trt = "b"` estimate than
the `trt = "a"` one, and the frame must keep `trt = "b"`. The
review's `p1_nolock` mutant fails it (below). The lane reads boot
mean 1.0917 2.1500 3.2083 against a Wald estimate of 1.0120 2.0831
3.1542 (`dev/postfit2-p2-checks.R`).

**P1-M1, `mm(g1, g2)`.** brms's reading is tractable, so it is
implemented rather than refused:

- `ce_group_vars()` names the member columns, so `re_formula = NULL`
  blanks the unset ones. Round 1 held them at level `"1"`, an
  observed group, in every band.
- `ce_level_plan()` reads each member against the term's levels. A
  row whose members are all observed reads the term observed, and the
  bootstrap holds it. Each distinct new member value is one new level,
  drawn once and shared, at its own placeholder level, which no
  observed member of the row names; two unset members are one new
  level with their weights added. A row with one member observed and
  one new reads the term both ways, which one bootstrap cannot
  simulate: refused by name on `band = "boot"`, answered on draws.
- An `mm()` term with a `by` variable has one by column per member.
  Its observed rows are held; a new member is refused by name on
  `band = "boot"` and on draws, since each member would need a
  placeholder in its own by-level's block.

brms at the same draws: with nothing set and with `g1 = "2"`, `g2`
unset, the draws band equals brms's to 2.2e-16 at seeds 1 to 3. With
two DIFFERENT unseen members (`g1 = "99"`, `g2 = "98"`) brms returns
exactly its nothing-set band, one shared draw with full weight;
frmtmb draws the two independently, as brms does for two unseen
levels of `(1 | g)`, and so does frmtmb's Wald band (width 2.78
against 3.75 for one new level). This is a deliberate difference,
named on `?conditional_effects`.

The Wald band: the section 9 claim of round 1 that it "carries its
variance" was false, since the grid read level `"1"`. It now reads a
new level: width 3.75 at nothing set, against 0.72 in round 1.

**P1-M3, the reuse key.** `ce_boot_key()` carries `re_formula` (a
formula as its text), and a mismatch says so: "predictions under a
different re_formula (re_formula = NULL there, ~(1 | g) here)". On
the review's `(1 + x | g)`, `g = "3"`, round 1 returned the `NULL`
band for the `~ (1 | g)` estimate.

**P1-M2, crossed `(1 | g) + (1 | h) + (1 | g:h)` at an unseen `g:h`.**
It does not fall out of the rework: the placeholder of the new `g:h`
would move `g` or `h`, which the observed terms read. The Wald band
answers (width 3.46, against 0.91 at an observed `g:h`); boot and
draws refuse by name. The limitation is named on
`?conditional_effects`.

**Minors.**

- A factor-smooth or a smooth with a `by` factor, set to a level the
  fit never saw, now stops in `ce_grids_build()` and names the
  variable, the level and the smooth, under any `re_formula`. Round 1
  returned an `NA` estimate under `NULL` and, under `NA`, a message
  that called the level `NA`. A grid that leaves such a factor unset
  (a `(1 | g)` beside the smooth, under `NULL`) stops too.
- A test holds one new level shared across panels: `effects =
  c("x", "z")` with `int_conditions` putting both panels at the same
  point, on the bootstrap (core) and on draws (sample). The review's
  `p1_noshare` mutant fails both.
- `ce_new_level_spec()`, `ce_boot_grids()`, `ce_draw_new_levels()` and
  `ce_new_level_key()` are removed, and their exports and aliases with
  them.
- `?conditional_effects` has a new section, "New group levels": the
  `gr()` and `mm()` readings, the smooth refusal, how a placeholder is
  placed and when it is refused, and the crossed limitation.

The numbers, round 1 then this round (`dev/postfit2-p2-checks.R`):

```
ARM r1: frmtmb from C:/Users/adf44/source/r/wt-postfit2-r1lib/frmtmb
B1 f = a, g = 3 | wald: est -0.8359 0.7768 2.3894 | width 0.5491 0.4064 0.5797
B1 f = a, g = 3 | boot40: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = b]) at the grid rows with g = 3: to place one it would move a column that anot
B1 f = a, g = 3 | draws: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = b]) at the grid rows with g = 3: to place one it would move a column that anot
B1 f = b, g = 9 | wald: est 1.8634 3.4760 5.0886 | width 0.6030 0.4496 0.5899
B1 f = b, g = 9 | boot40: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = a]) at the grid rows with g = 9: to place one it would move a column that anot
B1 f = b, g = 9 | draws: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = a]) at the grid rows with g = 9: to place one it would move a column that anot
B1 f = a | wald: est -0.7242 0.8884 2.5011 | width 0.8871 0.7951 0.8858
B1 f = a | boot40: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = b]) at the grid rows with g = NA: to place one it would move a column that ano
B1 f = a | draws: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = b]) at the grid rows with g = NA: to place one it would move a column that ano
B1 f = b | wald: est -0.4065 1.2061 2.8187 | width 6.1676 6.1543 6.1660
B1 f = b | boot40: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = b]) at the grid rows with g = NA: to place one it would move a column that ano
B1 f = b | draws: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = b]) at the grid rows with g = NA: to place one it would move a column that ano
B1 f = a, g = 99 | wald: est -0.7242 0.8884 2.5011 | width 0.8871 0.7951 0.8858
B1 f = a, g = 99 | boot40: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = b]) at the grid rows with g = 99: to place one it would move a column that ano
B1 f = a, g = 99 | draws: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g [f = b]) at the grid rows with g = 99: to place one it would move a column that ano
B2 trt = b: wald est 1.0120 2.0831 3.1542 | boot mean of refits 1.0917 2.1500 3.2083
M1 VarCorr sd 0.9071
M1 population (re_formula NA) | wald: est -0.0884 1.0848 2.2579 | width 1.1862 1.1327 1.1897
M1 nothing set: frame g1 1 g2 1
M1 nothing set | wald: est -0.2411 0.9321 2.1053 | width 0.7229 0.5790 0.6362
M1 nothing set | boot60: est -0.2411 0.9321 2.1053 | width 3.4776 3.4733 3.4708
M1 nothing set: boot width / wald width 4.8109 5.9986 5.4559
M1 nothing set | draws: est -0.2491 0.9268 2.1045 | width 0.5679 0.2923 0.6094
M1 g1 = 2, g2 = 2 | wald: est -0.3167 0.8564 2.0296 | width 0.6731 0.5808 0.6915
M1 g1 = 2, g2 = 2 | boot60: est -0.3167 0.8564 2.0296 | width 3.4442 3.5050 3.5029
M1 g1 = 2, g2 = 2: boot width / wald width 5.1167 6.0353 5.0659
M1 g1 = 2, g2 = 2 | draws: est -0.3226 0.8595 2.0364 | width 0.5122 0.3179 0.6044
M1 g1 = 2 | wald: est -0.2789 0.8943 2.0674 | width 0.5584 0.4003 0.5151
M1 g1 = 2 | boot60: est -0.2789 0.8943 2.0674 | width 2.2684 2.3785 2.4887
M1 g1 = 2: boot width / wald width 4.0626 5.9423 4.8312
M1 g1 = 2 | draws: est -0.2760 0.8941 2.0803 | width 0.5388 0.2515 0.5968
M1 g1 = 99, g2 = 98 | wald: est -0.0884 1.0848 2.2579 | width 3.7484 3.7318 3.7495
M1 g1 = 99, g2 = 98 | boot60: est -0.0884 1.0848 2.2579 | width 1.0590 1.0675 1.1635
M1 g1 = 99, g2 = 98: boot width / wald width 0.2825 0.2860 0.3103
M1 g1 = 99, g2 = 98 | draws: est -0.0911 1.0824 2.2721 | width 0.4933 0.1847 0.5986
M3 re_formula NULL | boot20: est 1.0430 2.2202 3.3975 | width 1.1182 0.2911 1.0418
M3 ~ (1 | g) reusing the NULL bootstrap: est 0.6570 2.2078 3.7586 | width 1.1182 0.2911 1.0418
M3 ~ (1 | g) | wald: est 0.6570 2.2078 3.7586 | width 2.3579 0.4339 2.2113
M2 g = 1, h = 1 (never together) | wald: est 0.8020 1.9986 3.1951 | width 3.4602 3.4398 3.4848
M2 g = 1, h = 2 (observed together) | wald: est -0.3669 0.8297 2.0263 | width 0.9148 0.8463 1.0237
M2 g = 1, h = 1 | boot20: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g:h) at the grid rows with g = 1, h = 1: to place one it would move a column that ano
M2 g = 1, h = 1 | draws: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g:h) at the grid rows with g = 1, h = 1: to place one it would move a column that ano
fs g = 99, re_formula NULL | wald: est NA NA NA | width NA NA NA
fs g = 99, re_formula NA | wald: ERROR New levels in the factor-smooth term s(x,g): NA. The term has no curve for them. Use allow_new_levels = TRUE to predict them at the population level; re_formula
done
ARM lane: frmtmb from C:/Users/adf44/source/r/wt-postfit2-lib/frmtmb
B1 f = a, g = 3 | wald: est -0.8359 0.7768 2.3894 | width 0.5491 0.4064 0.5797
B1 f = a, g = 3 | boot40: est -0.8359 0.7768 2.3894 | width 0.4852 0.2901 0.3790
B1 f = a, g = 3 | draws: est -0.8299 0.7847 2.4086 | width 0.7680 0.2502 0.6813
B1 f = b, g = 9 | wald: est 1.8634 3.4760 5.0886 | width 0.6030 0.4496 0.5899
B1 f = b, g = 9 | boot40: est 1.8634 3.4760 5.0886 | width 0.4453 0.3505 0.4575
B1 f = b, g = 9 | draws: est 1.8459 3.4709 5.1038 | width 0.7310 0.3335 0.7118
B1 f = a | wald: est -0.7242 0.8884 2.5011 | width 0.8871 0.7951 0.8858
B1 f = a | boot40: est -0.7242 0.8884 2.5011 | width 0.5914 0.5550 0.7344
B1 f = a | draws: est -0.7034 0.8790 2.5078 | width 1.0721 0.8692 1.1232
B1 f = b | wald: est -0.4065 1.2061 2.8187 | width 6.1676 6.1543 6.1660
B1 f = b | boot40: est -0.4065 1.2061 2.8187 | width 5.9121 5.9938 6.0755
B1 f = b | draws: est -0.4544 1.1583 2.7923 | width 5.7590 5.7067 5.6078
B1 f = a, g = 99 | wald: est -0.7242 0.8884 2.5011 | width 0.8871 0.7951 0.8858
B1 f = a, g = 99 | boot40: est -0.7242 0.8884 2.5011 | width 0.5914 0.5550 0.7344
B1 f = a, g = 99 | draws: est -0.7034 0.8790 2.5078 | width 1.0721 0.8692 1.1232
B2 trt = b: wald est 1.0120 2.0831 3.1542 | boot mean of refits 1.0917 2.1500 3.2083
M1 VarCorr sd 0.9071
M1 population (re_formula NA) | wald: est -0.0884 1.0848 2.2579 | width 1.1862 1.1327 1.1897
M1 nothing set: frame g1 NA g2 NA
M1 nothing set | wald: est -0.0884 1.0848 2.2579 | width 3.7484 3.7318 3.7495
M1 nothing set | boot60: est -0.0884 1.0848 2.2579 | width 3.3490 3.2622 3.3824
M1 nothing set: boot width / wald width 0.8935 0.8742 0.9021
M1 nothing set | draws: est -0.1330 1.0410 2.2133 | width 3.6280 3.5238 3.6043
M1 g1 = 2, g2 = 2 | wald: est -0.3167 0.8564 2.0296 | width 0.6731 0.5808 0.6915
M1 g1 = 2, g2 = 2 | boot60: est -0.3167 0.8564 2.0296 | width 0.8310 0.5972 0.7246
M1 g1 = 2, g2 = 2: boot width / wald width 1.2346 1.0283 1.0480
M1 g1 = 2, g2 = 2 | draws: est -0.3226 0.8595 2.0364 | width 0.5122 0.3179 0.6044
M1 g1 = 2 | wald: est -0.2026 0.9706 2.1438 | width 1.9224 1.8910 1.9268
M1 g1 = 2 | boot60: ERROR band = "boot" cannot cover this call: the group-level term (1 | mm(g1, g2)) is read at an observed level in some grid rows and at a new level in others (or, for
M1 g1 = 2 | draws: est -0.1895 0.9473 2.1343 | width 1.9524 1.8566 1.9437
M1 g1 = 99, g2 = 98 | wald: est -0.0884 1.0848 2.2579 | width 2.7800 2.7576 2.7815
M1 g1 = 99, g2 = 98 | boot60: est -0.0884 1.0848 2.2579 | width 2.5196 2.5343 2.5928
M1 g1 = 99, g2 = 98: boot width / wald width 0.9063 0.9190 0.9322
M1 g1 = 99, g2 = 98 | draws: est -0.0737 1.0722 2.2819 | width 2.4881 2.3656 2.2443
M3 re_formula NULL | boot20: est 1.0430 2.2202 3.3975 | width 1.1182 0.2911 1.0418
M3 ~ (1 | g) reusing the NULL bootstrap: ERROR boot = was not produced by a conditional_effects(band = "boot") call on this grid: its draws are predictions under a different re_formula (re_formula ...
M3 ~ (1 | g) | wald: est 0.6570 2.2078 3.7586 | width 2.3579 0.4339 2.2113
M2 g = 1, h = 1 (never together) | wald: est 0.8020 1.9986 3.1951 | width 3.4602 3.4398 3.4848
M2 g = 1, h = 2 (observed together) | wald: est -0.3669 0.8297 2.0263 | width 0.9148 0.8463 1.0237
M2 g = 1, h = 1 | boot20: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g:h) at the grid rows with g = 1, h = 1: to place one it would move a column that ano
M2 g = 1, h = 1 | draws: ERROR conditional_effects() cannot draw a new level of the group-level term (1 | g:h) at the grid rows with g = 1, h = 1: to place one it would move a column that ano
fs g = 99, re_formula NULL | wald: ERROR conditions sets g to "99", a level the fit never saw, and the smooth s(x,g) reads g: it has no curve for that level. Set g to one of 1, 2, 3, 4, 5
fs g = 99, re_formula NA | wald: ERROR conditions sets g to "99", a level the fit never saw, and the smooth s(x,g) reads g: it has no curve for that level. Set g to one of 1, 2, 3, 4, 5
done
```

brms at the same draws (`dev/postfit2-p2-brms.R`):

```
1: draws brms 200; r_g[3] max diff 2.78e-17; r_g[9] max diff 0
1: sd_g__Intercept:fa max rel diff 0
1 observed f = a, g = 3 (no draw): max diff estimate 0 lower 0 upper 1.11e-16
1 observed f = b, g = 9 (no draw): max diff estimate 0 lower 0 upper 0
1 f = a, g unset, seed 1: brms ERROR The following variables are missing in the draws object: {'sd_g__Intercept:f'} | frmtmb answers
1 f = a, g unset, seed 2: brms ERROR The following variables are missing in the draws object: {'sd_g__Intercept:f'} | frmtmb answers
1 f = a, g unset, seed 3: brms ERROR The following variables are missing in the draws object: {'sd_g__Intercept:f'} | frmtmb answers
1 f = b, g unset, seed 1: brms ERROR The following variables are missing in the draws object: {'sd_g__Intercept:f'} | frmtmb answers
1 f = b, g unset, seed 2: brms ERROR The following variables are missing in the draws object: {'sd_g__Intercept:f'} | frmtmb answers
1 f = b, g unset, seed 3: brms ERROR The following variables are missing in the draws object: {'sd_g__Intercept:f'} | frmtmb answers
1 f = a, g = 99, seed 1: estimate brms -0.7034 0.8790 2.5078 | frmtmb -0.7034 0.8790 2.5078
1 f = a, g = 99, seed 1: lower brms -1.2963 0.4884 1.9873 | frmtmb -1.2963 0.4884 1.9873; upper brms -0.2242 1.3576 3.1105 | frmtmb -0.2242 1.3576 3.1105
1 f = a, g = 99, seed 1: max diff estimate 0 lower 0 upper 0
1 f = a, g = 99, seed 2: estimate brms -0.7236 0.8950 2.5093 | frmtmb -0.7236 0.8950 2.5093
1 f = a, g = 99, seed 2: lower brms -1.2471 0.5022 2.0584 | frmtmb -1.2471 0.5022 2.0584; upper brms -0.2138 1.2661 2.9426 | frmtmb -0.2138 1.2661 2.9426
1 f = a, g = 99, seed 2: max diff estimate 0 lower 0 upper 0
1 f = a, g = 99, seed 3: estimate brms -0.7479 0.8832 2.5194 | frmtmb -0.7479 0.8832 2.5194
1 f = a, g = 99, seed 3: lower brms -1.1631 0.5489 2.0311 | frmtmb -1.1631 0.5489 2.0311; upper brms -0.2561 1.2703 2.9252 | frmtmb -0.2561 1.2703 2.9252
1 f = a, g = 99, seed 3: max diff estimate 1.11e-16 lower 0 upper 1.11e-16
1 f = b, g = 99, seed 1: estimate brms -0.4544 1.1583 2.7923 | frmtmb -0.4544 1.1583 2.7923
1 f = b, g = 99, seed 1: lower brms -3.0370 -1.3432 0.3700 | frmtmb -3.0370 -1.3432 0.3700; upper brms 2.7220 4.3635 5.9777 | frmtmb 2.7220 4.3635 5.9777
1 f = b, g = 99, seed 1: max diff estimate 5.55e-17 lower 4.44e-16 upper 8.88e-16
1 f = b, g = 99, seed 2: estimate brms -0.4501 1.1322 2.7146 | frmtmb -0.4501 1.1322 2.7146
1 f = b, g = 99, seed 2: lower brms -3.4184 -1.8413 -0.1976 | frmtmb -3.4184 -1.8413 -0.1976; upper brms 2.5701 4.1510 5.7647 | frmtmb 2.5701 4.1510 5.7647
1 f = b, g = 99, seed 2: max diff estimate 5.55e-17 lower 0 upper 0
1 f = b, g = 99, seed 3: estimate brms -0.3531 1.3227 2.8305 | frmtmb -0.3531 1.3227 2.8305
1 f = b, g = 99, seed 3: lower brms -3.1014 -1.4308 0.2151 | frmtmb -3.1014 -1.4308 0.2151; upper brms 2.2201 3.8560 5.4693 | frmtmb 2.2201 3.8560 5.4693
1 f = b, g = 99, seed 3: max diff estimate 2.22e-16 lower 0 upper 0
2: draws brms 200; r_mmg1g2[2] max diff 2.78e-17
2 g1 = g2 = 2 (no draw): max diff estimate 0 lower 0 upper 0
2 nothing set, seed 1: estimate brms -0.1330 1.0410 2.2133 | frmtmb -0.1330 1.0410 2.2133
2 nothing set, seed 1: lower brms -1.7580 -0.4790 0.5952 | frmtmb -1.7580 -0.4790 0.5952; upper brms 1.8700 3.0448 4.1995 | frmtmb 1.8700 3.0448 4.1995
2 nothing set, seed 1: max diff estimate 0 lower 0 upper 0
2 nothing set, seed 2: estimate brms -0.1029 1.0548 2.1800 | frmtmb -0.1029 1.0548 2.1800
2 nothing set, seed 2: lower brms -1.9422 -0.8118 0.3652 | frmtmb -1.9422 -0.8118 0.3652; upper brms 1.7091 2.9281 4.0795 | frmtmb 1.7091 2.9281 4.0795
2 nothing set, seed 2: max diff estimate 1.11e-16 lower 0 upper 0
2 nothing set, seed 3: estimate brms -0.0511 1.1835 2.2303 | frmtmb -0.0511 1.1835 2.2303
2 nothing set, seed 3: lower brms -1.7656 -0.6008 0.5901 | frmtmb -1.7656 -0.6008 0.5901; upper brms 1.6821 2.8351 4.0199 | frmtmb 1.6821 2.8351 4.0199
2 nothing set, seed 3: max diff estimate 0 lower 0 upper 0
2 g1 = 2, g2 unset, seed 1: estimate brms -0.1895 0.9473 2.1343 | frmtmb -0.1895 0.9473 2.1343
2 g1 = 2, g2 unset, seed 1: lower brms -1.1493 0.1610 1.2624 | frmtmb -1.1493 0.1610 1.2624; upper brms 0.8031 2.0176 3.2062 | frmtmb 0.8031 2.0176 3.2062
2 g1 = 2, g2 unset, seed 1: max diff estimate 0 lower 0 upper 0
2 g1 = 2, g2 unset, seed 2: estimate brms -0.2221 0.9781 2.1073 | frmtmb -0.2221 0.9781 2.1073
2 g1 = 2, g2 unset, seed 2: lower brms -1.1691 -0.0258 1.2476 | frmtmb -1.1691 -0.0258 1.2476; upper brms 0.7472 1.8945 3.0881 | frmtmb 0.7472 1.8945 3.0881
2 g1 = 2, g2 unset, seed 2: max diff estimate 0 lower 0 upper 1.11e-16
2 g1 = 2, g2 unset, seed 3: estimate brms -0.2195 1.0201 2.1627 | frmtmb -0.2195 1.0201 2.1627
2 g1 = 2, g2 unset, seed 3: lower brms -1.0380 0.1118 1.2106 | frmtmb -1.0380 0.1118 1.2106; upper brms 0.6831 1.8836 3.0334 | frmtmb 0.6831 1.8836 3.0334
2 g1 = 2, g2 unset, seed 3: max diff estimate 2.22e-16 lower 0 upper 0
2 g1 = 99, g2 = 98, seed 1: estimate brms -0.1330 1.0410 2.2133 | frmtmb -0.0737 1.0722 2.2819
2 g1 = 99, g2 = 98, seed 1: lower brms -1.7580 -0.4790 0.5952 | frmtmb -1.3567 -0.0930 1.1723; upper brms 1.8700 3.0448 4.1995 | frmtmb 1.1314 2.2726 3.4166
2 g1 = 99, g2 = 98, seed 1: max diff estimate 0.0685 lower 0.577 upper 0.783
2 g1 = 99, g2 = 98, seed 2: estimate brms -0.1029 1.0548 2.1800 | frmtmb 0.0510 1.2476 2.4004
2 g1 = 99, g2 = 98, seed 2: lower brms -1.9422 -0.8118 0.3652 | frmtmb -1.3026 -0.1216 0.9875; upper brms 1.7091 2.9281 4.0795 | frmtmb 1.2895 2.4503 3.6620
2 g1 = 99, g2 = 98, seed 2: max diff estimate 0.22 lower 0.69 upper 0.478
2 g1 = 99, g2 = 98, seed 3: estimate brms -0.0511 1.1835 2.2303 | frmtmb -0.0454 1.1036 2.2943
2 g1 = 99, g2 = 98, seed 3: lower brms -1.7656 -0.6008 0.5901 | frmtmb -1.2333 -0.0402 1.1201; upper brms 1.6821 2.8351 4.0199 | frmtmb 1.0179 2.1860 3.3773
2 g1 = 99, g2 = 98, seed 3: max diff estimate 0.0799 lower 0.561 upper 0.664
done
```

Seen to fail on the round-1 build, then pass:

```
## frmtmb ce-levels
r1: RESULT frmtmb test-ce-levels.R: tests=55 failed=10 error=2 skipped=0 warning=0 passed=45   [p2-ce-levels-r1.txt]
  fails: a gr(g, by = f) row reads only the term of its own f level 
  fails: an mm() term's members are grouping variables 
  fails: an mm() term with a by variable holds its observed members 
  fails: a bootstrap is not reused under another re_formula 
  fails: a smooth's factor at a level it has no curve for is refused 
lane: RESULT frmtmb test-ce-levels.R: tests=65 failed=0 error=0 skipped=0 warning=0 passed=65   [p2-ce-levels-lane.txt]
## frmtmb.sample postfit-draws
r1: RESULT frmtmb.sample test-postfit-draws.R: tests=49 failed=2 error=1 skipped=0 warning=0 passed=47   [p2-postfit-draws-r1.txt]
  fails: on draws, a gr(g, by = f) row reads its own f level's term 
  fails: on draws, an mm() term reads a level per member 
lane: RESULT frmtmb.sample test-postfit-draws.R: tests=51 failed=0 error=0 skipped=0 warning=0 passed=51   [p2-postfit-draws-lane.txt]
```

Mutants, each undoing one fix in the process
(`dev/postfit2-p2-mutant.R`, and the review's
`dev/postfit2-rev-runtest.R` for `p1_nolock` and `p1_noshare`):

```
RESULT frmtmb test-ce-levels.R [byall]: failed=0 error=1 passed=58   [p2-mut2-core-byall.txt]
  fails: a gr(g, by = f) row reads only the term of its own f level      0  TRUE
RESULT frmtmb test-ce-levels.R [mmskip]: failed=5 error=0 passed=60   [p2-mut2-core-mmskip.txt]
  fails:              an mm() term's members are grouping variables      3 FALSE
  fails: an mm() term with a by variable holds its observed members      2 FALSE
RESULT frmtmb test-ce-levels.R [mmgv]: failed=5 error=0 passed=60   [p2-mut2-core-mmgv.txt]
  fails:              an mm() term's members are grouping variables      4 FALSE
  fails: an mm() term with a by variable holds its observed members      1 FALSE
RESULT frmtmb test-ce-levels.R [rekey]: failed=1 error=0 passed=64   [p2-mut2-core-rekey.txt]
  fails: a bootstrap is not reused under another re_formula      1 FALSE
RESULT frmtmb test-ce-levels.R [nosmooth]: failed=1 error=1 passed=61   [p2-mut2-core-nosmooth.txt]
  fails: a smooth's factor at a level it has no curve for is refused      1  TRUE
RESULT frmtmb.sample test-postfit-draws.R [byall]: failed=0 error=1 passed=49   [p2-mut2-sample-byall.txt]
  fails: on draws, a gr(g, by = f) row reads its own f level's term      0  TRUE
RESULT frmtmb.sample test-postfit-draws.R [mmskip]: failed=1 error=0 passed=50   [p2-mut2-sample-mmskip.txt]
  fails: on draws, an mm() term reads a level per member      1 FALSE
RESULT frmtmb.sample test-postfit-draws.R [mmgv]: failed=2 error=0 passed=49   [p2-mut2-sample-mmgv.txt]
  fails: on draws, an mm() term reads a level per member      2 FALSE
RESULT lane frmtmb test-ce-levels.R [p1_nolock]: tests=65 failed=1 error=0 skipped=0 warning=0 passed=64   [p2-mut-p1_nolock.txt]
  fails: a new level's placeholder never moves a fixed effect's column 
RESULT lane frmtmb test-ce-levels.R [p1_noshare]: tests=65 failed=2 error=0 skipped=0 warning=0 passed=63   [p2-mut-p1_noshare.txt]
  fails: a new level is drawn once and shared by every panel 
RESULT lane frmtmb.sample test-postfit-draws.R [p1_noshare]: tests=51 failed=3 error=0 skipped=0 warning=0 passed=48   [p2-mut-sample-p1_noshare.txt]
  fails: on draws, a new level is drawn once and shared by every panel 
```

All tiers (`dev/postfit2-p2-all.sh`) and the check:

```
tier target, jobs: 19
logs with RESULT: 19
failed=0 error=0 skipped=0 warning=0 passed=2485
files with a failure, an error, a warning or a skip:
  none
files without a RESULT line:
tier core, jobs: 186
logs with RESULT: 186
failed=3 error=5 skipped=152 warning=0 passed=12394
files with a failure, an error, a warning or a skip:
RESULT frmtmb test-bcm-bart.R: tests=9 failed=0 error=0 skipped=1 warning=0 passed=8
RESULT frmtmb test-bcm-binomial.R: tests=23 failed=0 error=0 skipped=6 warning=0 passed=17
RESULT frmtmb test-bcm-data-analysis.R: tests=26 failed=0 error=0 skipped=6 warning=0 passed=20
RESULT frmtmb test-bcm-esp.R: tests=7 failed=0 error=0 skipped=2 warning=0 passed=5
RESULT frmtmb test-bcm-gaussian.R: tests=14 failed=0 error=0 skipped=3 warning=0 passed=11
RESULT frmtmb test-bcm-gcm.R: tests=11 failed=0 error=0 skipped=2 warning=0 passed=9
RESULT frmtmb test-bcm-latent-mixtures.R: tests=35 failed=0 error=0 skipped=5 warning=0 passed=30
RESULT frmtmb test-bcm-model-selection.R: tests=16 failed=0 error=0 skipped=3 warning=0 passed=13
RESULT frmtmb test-bcm-mpt.R: tests=17 failed=0 error=0 skipped=2 warning=0 passed=15
RESULT frmtmb test-bcm-psychophysics.R: tests=9 failed=0 error=0 skipped=2 warning=0 passed=7
RESULT frmtmb test-bcm-retention.R: tests=18 failed=0 error=0 skipped=3 warning=0 passed=15
RESULT frmtmb test-bcm-signal-detection.R: tests=9 failed=0 error=0 skipped=3 warning=0 passed=6
RESULT frmtmb test-bcm-simple.R: tests=9 failed=0 error=0 skipped=1 warning=0 passed=8
RESULT frmtmb test-brms-agreement.R: tests=170 failed=0 error=0 skipped=2 warning=0 passed=168
RESULT frmtmb test-brms-likelihood.R: tests=54 failed=0 error=0 skipped=37 warning=0 passed=17
RESULT frmtmb test-brms-methods.R: tests=62 failed=0 error=0 skipped=46 warning=0 passed=16
RESULT frmtmb test-brms-priors.R: tests=12 failed=0 error=0 skipped=12 warning=0 passed=0
RESULT frmtmb test-conditions.R: tests=145 failed=3 error=3 skipped=0 warning=0 passed=142
RESULT frmtmb test-data2.R: tests=29 failed=0 error=1 skipped=0 warning=0 passed=29
RESULT frmtmb test-drmtmb-agreement.R: tests=13 failed=0 error=0 skipped=13 warning=0 passed=0
RESULT frmtmb test-fuzz.R: tests=1 failed=0 error=0 skipped=1 warning=0 passed=0
RESULT frmtmb test-id-kron.R: tests=59 failed=0 error=1 skipped=0 warning=0 passed=59
RESULT frmtmb test-rl-example.R: tests=98 failed=0 error=0 skipped=2 warning=0 passed=96
files without a RESULT line:
tier gated, jobs: 31
logs with RESULT: 31
failed=13 error=0 skipped=0 warning=0 passed=2752
files with a failure, an error, a warning or a skip:
RESULT frmtmb test-brms-suite-brmsfit-helpers.R: tests=3 failed=2 error=0 skipped=0 warning=0 passed=1
RESULT frmtmb test-brms-suite-brmsformula.R: tests=16 failed=4 error=0 skipped=0 warning=0 passed=12
RESULT frmtmb test-brms-suite-methods.R: tests=162 failed=7 error=0 skipped=0 warning=0 passed=155
files without a RESULT line:
tier sample, jobs: 37
logs with RESULT: 37
failed=6 error=0 skipped=1 warning=0 passed=2204
files with a failure, an error, a warning or a skip:
RESULT frmtmb.sample test-brms-suite-brmsformula.R: tests=16 failed=4 error=0 skipped=0 warning=0 passed=12
RESULT frmtmb.sample test-brms-suite-methods.R: tests=61 failed=2 error=0 skipped=0 warning=0 passed=59
RESULT frmtmb.sample test-scale.R: tests=1 failed=0 error=0 skipped=1 warning=0 passed=0
files without a RESULT line:
tier ext, jobs: 91
logs with RESULT: 91
failed=0 error=0 skipped=14 warning=4 passed=4234
files with a failure, an error, a warning or a skip:
RESULT frmtmb.coupling test-scale.R: tests=5 failed=0 error=0 skipped=5 warning=0 passed=0
RESULT frmtmb.eam test-scale.R: tests=3 failed=0 error=0 skipped=3 warning=0 passed=0
RESULT frmtmb.latent test-scale.R: tests=2 failed=0 error=0 skipped=2 warning=0 passed=0
RESULT frmtmb.learn test-scale.R: tests=2 failed=0 error=0 skipped=2 warning=0 passed=0
RESULT frmtmb.learn test-stan-identity.R: tests=75 failed=0 error=0 skipped=0 warning=4 passed=71
RESULT frmtmb.ode test-scale.R: tests=1 failed=0 error=0 skipped=1 warning=0 passed=0
RESULT frmtmb.spline test-scale.R: tests=1 failed=0 error=0 skipped=1 warning=0 passed=0
files without a RESULT line:

the core tier's three failing files with frmtmb attached, as tests/testthat.R attaches it (dev/postfit2-rev-runtest.R):
RESULT lane frmtmb test-conditions.R: tests=150 failed=0 error=0 skipped=0 warning=0 passed=150   [p2-testenv-lane-conditions.txt]
RESULT base frmtmb test-conditions.R: tests=150 failed=0 error=0 skipped=0 warning=0 passed=150   [p2-testenv-base-conditions.txt]
RESULT lane frmtmb test-data2.R: tests=31 failed=0 error=0 skipped=0 warning=0 passed=31   [p2-testenv-lane-data2.txt]
RESULT base frmtmb test-data2.R: tests=31 failed=0 error=0 skipped=0 warning=0 passed=31   [p2-testenv-base-data2.txt]
RESULT lane frmtmb test-id-kron.R: tests=60 failed=0 error=0 skipped=0 warning=0 passed=60   [p2-testenv-lane-id-kron.txt]
RESULT base frmtmb test-id-kron.R: tests=60 failed=0 error=0 skipped=0 warning=0 passed=60   [p2-testenv-base-id-kron.txt]

R CMD check --as-cran (dev/postfit2-check.ps1):
frmtmb Status: 1 NOTE
frmtmb.sample Status: OK
```

## 7e. Punch round 3 (the final check)

**P2-B1, the `mm()` rule "members with the same new value read one
new level, and their weights add" was not pinned.** The review's
`p2_mmsplit` mutant (each new member drawn on its own) passed every
file, yet made the nothing-set band 25 to 38% narrower on boot and
draws (`dev/postfit2-rev-p2-mmsplit.R`, data seed 43). Two tests now
compare the nothing-set band with the Wald band as a ratio, bounded
by 0.8 and 1.25:

- core `test-ce-levels.R`, boot over Wald: 0.87 to 0.90 on the lane,
  0.67 to 0.69 under the mutant;
- sample `test-postfit-draws.R`, 200 hand-built draws over Wald: 0.94
  to 0.97 on the lane, 0.60 to 0.66 under the mutant.

**P2-M1, `mm(g1, g2)` beside `(1 | g1)` with `g1` unset or `"99"`.**
Boot and draws refused, because `ce_plan_part()` placed the plain
term first and then refused the `mm()` member it had moved. The fix is
contained: a new `mm()` member that a plain new term already placed
keeps that placeholder, since the two terms' effects are separate
coefficient slots. It is still refused where the members sharing one
new value would need two different placeholders, or where the placed
value is an observed member's level. `?conditional_effects` says that
two new terms can share a placeholder column.

Against brms at the same 100 draws (`dev/postfit2-p3-brms.R`): the two
packages consume random numbers in different orders with two new
terms, so matched seeds do not give equal bands (the review's
`dev/postfit2-rev-p2-mm.R` shows both answering, with estimates
0.06 to 0.26 apart). The check is of the law: the per-draw new-level
offset over its model SD has variance 1.

```
R = 30:
brms sd_mmg1g2 / sd_g1 max rel diff 0 0
population spaghetti, brms vs frmtmb: max diff 2.22e-16
nothing set: 3000 offsets per package; var(offset / frmtmb's model SD): brms 0.9956, frmtmb 0.9959 (SE about 0.0258); brms's own model: 0.9956
g1 = 99, g2 unset: 3000 offsets per package; var(offset / frmtmb's model SD): brms 1.2272, frmtmb 1.0561 (SE about 0.0258); brms's own model: 0.9956
done
R = 120:
brms sd_mmg1g2 / sd_g1 max rel diff 0 0
population spaghetti, brms vs frmtmb: max diff 2.22e-16
nothing set: 12000 offsets per package; var(offset / frmtmb's model SD): brms 0.9994, frmtmb 1.0085 (SE about 0.0129); brms's own model: 0.9994
g1 = 99, g2 unset: 12000 offsets per package; var(offset / frmtmb's model SD): brms 1.2324, frmtmb 1.0157 (SE about 0.0129); brms's own model: 0.9994
done
```

With nothing set the two agree. With `g1 = "99"` and `g2` unset brms
follows its own shared draw (the `data_gr_local()` artifact above):
its offsets are the nothing-set offsets, bit for bit.

Tests, on the round-2 build (installed before this fix), then the
lane; and the review's `p2_mmsplit` mutant on the lane:

```
RESULT frmtmb test-ce-levels.R: tests=65 failed=0 error=1 skipped=0 warning=0 passed=65   [p3-ce-levels-r2.txt]
  fails: mm() beside (1 | g1) shares g1's placeholder when both are new 
RESULT frmtmb.sample test-postfit-draws.R: tests=53 failed=0 error=1 skipped=0 warning=0 passed=53   [p3-postfit-draws-r2.txt]
  fails: on draws, mm() beside (1 | g1) draws both new levels 
RESULT frmtmb test-ce-levels.R: tests=69 failed=0 error=0 skipped=0 warning=0 passed=69   [p3-lane-test-ce-levels.R.txt]
RESULT frmtmb test-ce-options.R: tests=32 failed=0 error=0 skipped=0 warning=0 passed=32   [p3-lane-test-ce-options.R.txt]
RESULT frmtmb.sample test-postfit-draws.R: tests=57 failed=0 error=0 skipped=0 warning=0 passed=57   [p3-lane-test-postfit-draws.R.txt]
RESULT lane frmtmb test-ce-levels.R [p2_mmsplit]: tests=69 failed=1 error=0 skipped=0 warning=0 passed=68   [p3-mmsplit-test-ce-levels.R.txt]
  fails: an mm() term's members are grouping variables 
RESULT lane frmtmb test-ce-options.R [p2_mmsplit]: tests=32 failed=0 error=0 skipped=0 warning=0 passed=32   [p3-mmsplit-test-ce-options.R.txt]
RESULT lane frmtmb.sample test-postfit-draws.R [p2_mmsplit]: tests=57 failed=2 error=0 skipped=0 warning=0 passed=55   [p3-mmsplit-test-postfit-draws.R.txt]
  fails: on draws, an mm() term reads a level per member 
  fails: on draws, mm() beside (1 | g1) draws both new levels 
```

R CMD check --as-cran, after the fix:

```
frmtmb 19:04 Status: 1 NOTE
frmtmb.sample 18:51 Status: OK
```

## 8. Not done, and defects found but not fixed

- `dev/brms-api-diff.md` rows 99 and 398 to 402 still list
  `spaghetti`, `select_points` and `too_far` as gaps. It is a record
  outside `dev/brmsport-*`; left for the consolidation so the record
  keeps one author per round.
- `plot()` of a conditional-effects object still refuses brms's
  `plot`, `rug` and `stype` arguments (rows 154, 162, 164, 169). Not
  in this lane's list.
- `posterior_average()`'s refusals quote the variable without brms's
  quotes; rule 1 allows it.
- `pp_average()` and `model_weights()` are not exported; brms has
  both. `posterior_average()` computes the weights internally.

## 9. Version

Both packages gain exported functions, so both bumps are minor.
frmtmb.sample calls `cs_build()`, `cs_coef()`, `cs_frame()`,
`cs_finalize()`, `cs_probs()`, `ce_spaghetti()`,
`ce_spaghetti_check()`, `ce_check_distance()` and the new arguments of
`ce_grids_build()` and `ce_finalize()`, so its frmtmb floor must rise
to the frmtmb version that carries this lane.

Punch round 1 adds to both lists:

- Not done: frmtmb draws a new level from the fitted Gaussian law;
  brms's default `sample_new_levels = "uncertainty"` picks an observed
  level's draw instead and gives narrower new-level bands (the review's
  "for the user" note). The comparison above is against brms's
  `"gaussian"`, which is the same law.
- Not done: a new level of a term whose grouping is not one column per
  variable (`mm()` membership, and similar composite groupings) gets no
  drawn level on the bootstrap or on draws, as before; the Wald band
  carries its variance. CORRECTED in punch round 2: this was false.
  With nothing set, the grid held both members at level `"1"`, an
  observed group, so no band carried a new level's variance: the Wald
  band was level 1's (width 0.58 to 0.72), and the boot band was 4.8
  to 6.0 times it because the term was redrawn. Section 7d fixes it.
- Found, not this lane's: `frmtmb.sample`'s `test-scale.R` and the
  extensions' `test-scale.R` skip behind `FRMTMB_SCALE_TESTS`, a gate
  no tier sets.
- Version: `frm_bootstrap()`'s change is BREAKING, so the frmtmb bump
  is minor at least, as already needed. frmtmb.sample now also calls
  `ce_level_plan()`, `ce_plan_eval()` and `ce_plan_has_new()`.
- Pre-existing, both arms: `frmtmb.learn`'s `test-stan-identity.R`
  lets 4 warnings escape on this lane's core and on the base build
  alike (`dev/postfit2-p1-ext-{lane,base}/`).

Punch round 2 adds:

- Not done (P1-M2): crossed `(1 | g) + (1 | h) + (1 | g:h)` with `g`
  and `h` at observed levels never seen together. The Wald band
  answers; `band = "boot"` and draws refuse by name, and brms answers.
  Drawing it would need the new level's contribution added outside the
  design (a placeholder cannot move `g` or `h`). Named on
  `?conditional_effects`.
- Not done: an `mm()` term with a `by` variable draws no new level on
  `band = "boot"` or on draws (refused by name); its observed rows
  work, and the Wald band answers.
- Deliberate difference: two DIFFERENT unseen `mm()` members are two
  new levels here, drawn independently; brms 2.23.0 gives them one
  shared draw (its band equals its nothing-set band exactly). Named on
  `?conditional_effects`. Round 3: the shared draw is an artifact of
  brms's `data_gr_local()`, which numbers unseen values per member, so
  the first unseen value of each member gets the same index whatever
  the value is (the reviewer's reading). It is not a design choice.
- brms 2.23.0 stops on `gr(g, by = f)` with `g` unset and
  `sample_new_levels = "gaussian"` ("missing sd_g__Intercept:f");
  frmtmb answers. A brms defect, not reported upstream by this lane.
- The ported brms-suite rows fail as in round 1: core
  brmsfit-helpers 2, brmsformula 4, methods 7; frmtmb.sample
  brmsformula 4, methods 2; 19 in all, all stale expectations of
  earlier lanes.
- Runner artifact: `dev/postfit2-p2-all.sh core` runs core files
  without attaching frmtmb, so `test-conditions.R`, `test-data2.R` and
  `test-id-kron.R` cannot find `bf()` or internal helpers there. With
  frmtmb attached they pass on this lane and on the base build
  (section 7d).
