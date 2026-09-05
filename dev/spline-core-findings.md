# SPLINE-CORE lane: findings

Worktree `C:/Users/adf44/source/r/frmtmb-wt-spline-core`, branch
`wt-spline-core`, base 2210aa1 (frmtmb 0.51.0, frmtmb.spline 0.1.0).
Private library `scratchpad/sc-lib`: frmtmb 0.51.0 (worktree),
frmtmb.spline 0.1.0, flexsurv 2.3.2, gratia 0.11.2, brokenstick 2.7.0.
R 4.6.1. Main checkout untouched.

Contract: `dev/spline-seam-proposal.md`, Parts 1a, 1b, 1c, 3, 2, then the
extension adoption.

Baseline counted before any edit: 109 files in `tests/testthat/test-*.R`.

---

## Order of work

1. Part 3, the `all.vars()` parser defect (smallest, and Part 2 needs it).
2. Part 1a, the joint-covariance accessor.
3. Part 1b, `frm_lp_basis()`.
4. Part 1c, the `lccdf` slot.
5. Part 2, `ps()`.
6. Extension adoption, 0.2.0.
7. Verification.

## Status

(filled in as each part lands)

---

## Part 3. The `all.vars()` parser defect. LANDED.

Reproduced on 2210aa1 before the fix, in one process:

```
all.vars(quote(f(x)(y)))           -> y            (x gone)
all.vars(quote(a * curry(tv)(zv))) -> a zv         (tv gone)
all.names(quote(f(x)(y)))          -> f x y

frm(bf(y ~ a * curry2(tv, zv), a ~ 1, nl = TRUE))       -> a = 2.005982
frm(bf(y ~ a * curry(tv)(zv),  a ~ 1, nl = TRUE))       -> Error: the
      nonlinear formula body could not be evaluated: object 'tv' not found
frm(bf(y ~ a * curry(tv)(zv) + 0 * tv, a ~ 1, nl = TRUE)) -> a = 2.005982
```

Fix: `nl_body_vars()` at **R/parse.R:1300** replaces `all.vars()` at the
three sites that read a nonlinear body, **R/parse.R:1285** (`nl_dpar()`),
**R/parse.R:1445** (the `nl = TRUE`-by-implication test) and
**R/parse.R:1506** (the "not used in the model formula" refusal).

The walker descends into a call in function position and keeps
`all.vars()`'s own conventions everywhere else: `$` and `@` walk only
their first argument, `::` and `:::` collect nothing, and a `function`
literal is handed back to `all.vars()` because its formals are a
pairlist rather than a call. Without those four cases the walker would
have collected `b` from `d$b` and `stats` from `stats::rnorm`, which is
a regression in the opposite direction.

After the fix the curried spelling fits to **a = 2.005982**, the same
number to seven figures as both working spellings.

## Part 1a. `frm_joint_cov()`. LANDED.

**R/predict.R:2785**. A thin exported wrapper on the memoized
`get_joint_cov()`, plus a `labels` element the internal object did not
carry: one label per row (`beta.<coefficient>`, `b.<block>.<level>.<name>`,
`theta.<i>`), built by `joint_coef_labels()` (**R/predict.R:2801**) and
`b_coef_labels()` (**R/predict.R:2820**).

It is the route for the autoscaled case because it reads the cache
`sdr_of()` filled through `autoscale_sdreport()`; nothing about the
accessor has to know that, which is the point.

## Part 1b. `frm_lp_basis()`. LANDED.

**R/predict.R:2914**, returning
`list(eta, A, coef_pos, V, coef_names, extra_var, nonest)`.

`predict(se.fit = TRUE)` is now written in terms of it
(**R/predict.R:1230-1239**): the seam's own test, and it passes.

Measured, one process, `y ~ s(x, k = 8) + (1 | g)`, n = 150, 12-point
grid:

| route | eta vs `predict()` | se vs `predict(se.fit = TRUE)` |
|---|---:|---:|
| `re.form = NA` (p = 8) | 0.000e+00 | 2.776e-17 |
| `re.form = NULL` (p = 29) | 0.000e+00 | 2.776e-17 |
| in sample | 0.000e+00 | 0.000e+00 |

The nonlinear branch, `y ~ ult * (1 - exp(-exp(lrc) * t))` with
`ult ~ 1 + (1 | id)`, 20 subjects: `A` is a 7 x 2 Jacobian, and against
a central difference of `predict()` in each coefficient it agrees to
**4.48e-10**, which is the finite difference's own error. `eta` matches
`predict()` to 0.000e+00. `predict(se.fit = TRUE)` stays refused for a
nonlinear predictor, as the proposal says it must; this is the route.

## Part 1c. The `lccdf` slot. LANDED, at five sites plus a sixth.

The five the proposal enumerated:

1. `frmtmb_family(lccdf = )` and its validation, **R/families.R:238** and
   **R/families.R:277**;
2. `frmtmb_ad_overload()` wrapping beside `lcdf`'s, **R/families.R:277**;
3. `fam_lccdf()`, the arity shim beside `fam_lcdf()`,
   **R/families.R:4106**, plus `has_lccdf()`;
4. the use site, **R/objective.R:100-124**;
5. the frame guard, **R/frame.R:1300**, which now admits `cens()` on a
   family with `lccdf` and no `lcdf` when the codes are right-censoring
   only.

The sixth is the one the proposal did not enumerate and the review
asked for by another route: **R/fit-end.R**, a new file, called once
from **R/fit.R:999**. A family may declare `post$fit_check(fit, resp)`
and it runs when the fit finishes. That is what closes the "logLik() and
AIC() cannot be gated from an extension" limit the proposal recorded:
a family cannot make `logLik()` refuse, but it can now say something at
the moment the fit is handed over.

### The censored term, before and after

One process, `row_lpdf()` called directly with the family's `lccdf`
present and removed. `log(1 - F)` is the "before" column and is the
arithmetic on 2210aa1.

Gaussian, right censoring, `log S` at `z` standard deviations:

| z | true | `log(1 - F)` | `lccdf` |
|---:|---:|---:|---:|
| 1.0 | -1.84102 | -1.84102 | -1.84102 |
| 5.0 | -15.06500 | -15.06500 | -15.06500 |
| 8.0 | -35.01344 | **-34.94504** | -35.01344 |
| 8.3 | -37.49422 | **-Inf** | -37.49422 |
| 20 | -203.91716 | -Inf | -203.91716 |
| 37 | -689.03059 | -Inf | -689.03059 |
| 100 | -5005.52421 | -Inf | -5005.52421 |
| 500 | -125007.13355 | -Inf | -125007.13355 |

Weibull (`shape = 3`), `log S = -(q/scale)^shape`: exact to -364.58 at
q = 8, where `log(1 - F)` is `-Inf` from q = 4 (`log S = -45.57`).

The GRADIENT is the failure, not the value. `d/d eta` of the censored
term, gaussian:

| z | true | `log(1 - F)` | `lccdf` |
|---:|---:|---:|---:|
| 5 | 5.18650e+00 | 5.18650e+00 | 5.18650e+00 |
| 8 | 8.12137e+00 | **7.58447e+00** | 8.12137e+00 |
| 10 | 1.00981e+01 | **Inf** | 1.00981e+01 |
| 20 | 2.00498e+01 | **Inf** | 2.00498e+01 |
| 40 | 4.00250e+01 | **NaN** | 4.00250e+01 |

The old form is already 6.6 percent low at z = 8 and unusable above it.
Even at z = 5, where it looks exact, it has lost eight digits
(-15.0649983938340 against -15.0649983939887); the test asserts that
too, so the claim is not "closer" but "exact where the other is not
representable".

### The five families wired, and the two not

`gaussian()`, `lognormal()`, `exponential()`, `weibull()`, `cox()`. The
last three have `log S` in closed form; the first two use
`RTMB::pnorm(..., lower.tail = FALSE, log.p = TRUE)`, which the
measurement above shows is exact to z = 500 with the right derivative.

The contract asked for six. There are five, and the sixth is a measured
refusal rather than an omission:

* `poisson()` is discrete and `cens()` is refused for discrete
  families, so the slot would be unreachable.
* `inverse.gaussian()` gains NOTHING.
  `RTMBdist::pinvgauss(lower.tail = FALSE, log.p = TRUE)` is computed on
  the probability scale and reaches `-Inf` at the same place
  `log(1 - F)` does. Measured at `mean = 1, shape = 2`: both give
  -33.79 where the truth is -33.75 at q = 30, and both give `-Inf` from
  q = 40 (truth -44.17). Wiring it would have been a slot that changed
  nothing.

### What it does not close

Left censoring is still `log(F(y) - Flb)`, interval censoring is still a
difference of CDFs, and the truncation normalizer at **R/objective.R:136**
is still `log(Fub - Flb)`. A LEFT-TRUNCATED survival model meets the
identical problem from the other side. Right censoring UNDER `trunc(ub)`
is exact, because the windowed difference is taken as a log difference,
`lS(y) + log1p(-exp(lS(ub) - lS(y)))`, and never round-trips through a
probability.

## The bug the recovery study found

Not in the contract, in the path of it, and it is core's.

`b ~ 0 + (1 | g)` - a nonlinear parameter that is purely a random
effect, which is how the paper's simulation designs write their
amplitude and phase - produced a design matrix with ZERO columns, and
the two coefficient-naming lines then produced the single name `"b_"`,
because `paste()` recycles to its longest argument and
`paste("b", character(0), sep = "_")` is `"b_"`, not `character(0)`.

The consequence was not cosmetic. The parameter template carried an
entry no linear predictor indexed (`idx` stayed `integer(0)`), it
entered no likelihood, and the outer Hessian was singular in exactly
that direction. Measured on the sine design: two such phantom entries,
both estimated at exactly 0, and the outer Hessian's two smallest
eigenvalues 1.2e-15 and -1.3e-13, spanning precisely those two
coordinates. `vcov()` and every standard error came back `NaN` with
"the model is probably overparameterized", so the recovery study could
not run at all: 0 of 25 replications produced a curve.

Reproduced on the main-built frmtmb 0.50.0 in the user library, so it
predates this lane. Fixed at **R/frame.R:2035**. After the fix the same
fit reaches the same optimum to every printed digit (0.8109, 16.4915,
-0.4987, ...) and the outer Hessian's smallest eigenvalue is 5.9e-03.
25 of 25 replications now converge.

## Part 2. `ps()`. LANDED.

**R/ps.R**, new file, plus the frame hook at **R/frame.R:1476-1517**,
the block resolution at **R/frame.R:2311**, the objective's evaluation
frame at **R/objective.R:290**, `eval_dpars()` at **R/predict.R:541**,
`predict(newdata =)` at **R/predict.R:1215**, and the refusals and the
fit-end coverage report through **R/fit.R:768** and **R/fit.R:999**.

### The basis, re-measured

Against `splines::splineDesign()`, cubic unless stated, one process.
The proposal measured 6.82e-14 at 8 interior knots; the figure grows
with `k`, which is why `k > 50` is refused.

| configuration | max abs | value with x as AD | taped d/dx vs `derivs = 1` |
|---|---:|---:|---:|
| k = 12 on [0, 1.5] (the proposal's) | 3.126e-13 | 3.126e-13 | 4.547e-13 |
| k = 15 on [-17, 131] (SMOCC) | 1.251e-12 | 1.251e-12 | 1.954e-14 |
| k = 15 on [0, 1] | 4.547e-13 | 4.547e-13 | 1.563e-12 |
| k = 20 on [0, 1] | 2.046e-12 | 2.046e-12 | 6.821e-12 |
| k = 30 on [0, 1] | 5.912e-12 | 5.912e-12 | 1.819e-11 |
| k = 40 on [0, 1] | 1.819e-11 | 1.819e-11 | 4.547e-11 |
| k = 15, quadratic | 7.105e-14 | 7.105e-14 | 1.137e-13 |

The AD column is bit-identical to the non-AD one at every size. The
relative error is larger (1e-8 to 1e-11) only because basis values near
a boundary knot are themselves near zero, which is the proposal's own
caveat. Partition of unity holds to the same absolute figure.

The summation ORDER matters and the implementation commits to one:
`ps_value()` sums per BASIS FUNCTION and then over coefficients, not per
knot. The two are algebraically identical; collapsing over knots first
performs the divided difference's cancellation across the whole knot
vector instead of inside each basis function's `ord + 1` window, and the
measured accuracy above is the per-basis one.

### The identity with the paper's likelihood

The SMOCC height model, `brokenstick::smocc_200`, 1906 rows after
dropping 36 missing heights, 200 children. Age in weeks
(`52.1775 * age`), gestational age centered at 40.

```r
frm(bf(hgt ~ int + exp(amp) * ps(age + shift, k = 15, pad = 0.25),
       int ~ sex + (1 | id), amp ~ 0 + sex,
       shift ~ 0 + ga + (1 | id), nl = TRUE), smocc, gaussian())
```

The reference is the paper's joint density written out in plain R, with
the basis from `splines::splineDesign()` rather than from the
divided-difference construction the objective tapes. Nothing was read
from github.com/matteodales/snmmTMB.

| point | frmtmb | reference | abs | rel |
|---|---:|---:|---:|---:|
| at the optimum | 3653.3851349042 | 3653.3851349043 | **9.413e-11** | 2.577e-14 |
| at a perturbed vector | 3917.4767377491 | 3917.4767377485 | 5.957e-10 | 1.521e-13 |

The eigensplit's defining properties, checked independently of the sign
convention `eigen()` happens to pick: `S U0` is 3.98e-15, `Us' S Us - I`
is 2.14e-13, and the sum-to-zero constraint holds to 6.07e-16.

### Against the paper's Table 1

| parameter | paper (95% Wald) | frmtmb |
|---|---|---:|
| intercept beta0 | 68.2 (67.6-68.7) | 68.754 |
| sex intercept beta1 | 1.80 (0.92-2.68) | 1.806 |
| sex scale beta2 | 0.00 (-0.01-0.01) | 0.000 |
| GA shift beta3 | 1.00 (0.73-1.27) | 1.050 |
| sd_b1 | 2.86 (2.56-3.17) | 2.868 |
| sd_b2 | 3.28 (2.89-3.68) | 3.324 |
| sigma | 1.05 (1.01-1.09) | 1.057 |

Six of seven inside the published interval and most of them at the
published precision. The intercept is 0.554 higher, and that is an
identifiability convention.

The attribution has to be careful, and the first draft of this document
was not. `ps()` constrains the spline's COEFFICIENTS to sum to zero. The
paper states only that it applies a sum-to-zero constraint "to ensure
identifiability of an intercept term" and never says on what, so what
can be claimed is a measurement rather than a reading:

```
beta0 (ps sum-to-zero on the coefficients) : 68.754333
mean(fitted curve at the data)             : -0.55609575
beta0 + mean(curve)                        : 68.198237
the paper's beta0                          : 68.2
difference                                 : -0.0018
```

Re-centering the curve to its own fitted values puts the intercept
within 0.002 cm of the published one, so the two describe the same
fitted surface. The centering is specific: on an even 0 to 130 week grid
the same arithmetic gives 76.61, because the data are not evenly spread
over age. Every shape-dependent quantity agrees regardless of the
convention.

Wall time 21.4 s for the whole 200 children at k = 15.

### Parameter recovery on the paper's designs

25 replications per setting, one process, `(n, m) = (20, 20)` so 21
observations on each of 20 subjects. Converged 25 of 25 in all three.

sine, `sigma^2 = 0.4`, D = diag(0.25, 0.16, 0.04):

| quantity | truth | mean | sd | bias |
|---|---:|---:|---:|---:|
| sigma | 0.6325 | 0.6344 | 0.0208 | +0.0019 |
| sd(b1) | 0.3162 | 0.3193 | 0.0642 | +0.0031 |
| sd(b2) | 0.2530 | 0.2361 | 0.0615 | -0.0169 |
| sd(b3) | 0.1265 | 0.0897 | 0.0632 | -0.0368 |

curve RMSE 0.1051, pointwise 95 percent coverage 0.971.

sine, `sigma^2 = 0.4`, D = diag(1, 0.25, 0.16):

| quantity | truth | mean | sd | bias |
|---|---:|---:|---:|---:|
| sigma | 0.6325 | 0.6342 | 0.0220 | +0.0017 |
| sd(b1) | 0.6325 | 0.6561 | 0.1107 | +0.0237 |
| sd(b2) | 0.3162 | 0.3005 | 0.0613 | -0.0158 |
| sd(b3) | 0.2530 | 0.2294 | 0.0739 | -0.0236 |

curve RMSE 0.1656, coverage 0.971.

bell, `sigma^2 = 0.2`, D = diag(2, 2):

| quantity | truth | mean | sd | bias |
|---|---:|---:|---:|---:|
| sigma | 0.4472 | 0.4507 | 0.0166 | +0.0035 |
| sd(b1) | 0.6325 | 0.6049 | 0.1385 | -0.0275 |
| sd(b2) | 0.6325 | 0.6663 | 0.3909 | +0.0338 |

curve RMSE 0.1749, coverage 0.925.

sigma is recovered to under one percent everywhere. The variance
components are recovered with a downward bias that grows as the
component shrinks: the sine design's phase-shift SD, truth 0.1265, comes
back at 0.0897 with a replication SD of 0.0632, so 20 subjects with 21
observations each barely identify it. That is a property of the design
rather than of the implementation, and the paper reports coverage rather
than variance-component bias for the same reason.

Coverage is pointwise and slightly conservative (0.925 to 0.971 against
a nominal 0.95); the paper reports SIMULTANEOUS coverage, which is a
different quantity and needs `frm_curve()`.

## The extension adoption. frmtmb.spline 0.2.0.

`Depends: frmtmb (>= 0.52.0)`.

### What the accessor replaced

`sp_joint_cov()`, `sp_curve_design()`, `sp_check_linear()`, `sp_probe()`
and `sp_coef_pos()` are **deleted**, not rewired: 200 lines of
`extensions/frmtmb.spline/R/curve-cov.R` become one call to
`frmtmb::frm_lp_basis()` inside `sp_curve_parts()`. With them go

* the read of `fit$cache$Vjoint`, which was the package's only
  unsanctioned dependency;
* the sparse `sdreport()` fallback and the autoscale refusal beside it,
  which existed only because there was no exported route to the
  autoscaled joint precision;
* the argument-ordering trick that forced the check call to run before
  the covariance, which existed only to keep the fallback unreachable;
* the unit-perturbation design rebuild and the linearity probe that
  guarded it.

`?frm_curve`'s section changes from "The one internal this reaches into"
to "The route to the covariance", and its Cost section drops the design
rebuild it used to quantify.

`n_predict` is **1** now, on every model, and it is the check call.
It was 11 on `s(x, k = 10)` and 32 on the 110-coefficient factor-smooth
model; both are pinned at 1 in `tests/testthat/test-curve.R`, and the
second is the one that matters, because the old count grew with the
number of contributing coefficients and this one does not.

### What the accessor bought

| row | 0.1.0 | 0.2.0 | measured |
|---|---|---|---|
| `frm_curve x rr` | conditional | **works** | 0.000e+00 at both `re.form` settings; 9 theta columns at `re.form = NULL`, which the perturbation could not see and which put its standard errors 27 percent out |
| `frm_curve x gp` | untested | **works** | 0.000e+00, `extra_var` 7.3e-07 to 8.5e-07 |
| `frm_curve x nl` | untested (refused in fact) | **works** | the SMOCC warped curve, drawn in `vignette("case-studies")` section 13 |

The nonlinear case is the one with no second route: `predict(se.fit =
TRUE)` refuses a nonlinear predictor, so `cov_rel_error` comes back `NA`
and `print()` says the check did not run rather than reporting a passing
one that never did.

One honesty note about the check. For a LINEAR predictor the two sides
now come from the same core helper, so their agreement is 0.000e+00 by
construction and is no longer independent evidence about core. What it
still verifies is this package's own indexing and assembly: the
transpose, the `extra_var` addition, the grid orientation. The
independent evidence is in core's own suite, where the pre-rewrite
comparison measured 2.776e-17.

### The censored term, in the extension

`royston_parmar()` supplies `lccdf`, in closed form on all three scales:
`-exp(eta)` (hazard, capped at `eta = 700` where the double runs out,
not at the log density's 30, which is a different question), 
`-log1p(exp(eta))` (odds), and
`pnorm(eta, lower.tail = FALSE, log.p = TRUE)` (normal).

On the 600-subject design the review built, one group-A subject censored
at t = 50 far beyond every event time:

| | frmtmb 0.51.0 (the review) | frmtmb 0.52.0 |
|---|---:|---:|
| reported `logLik` | -365.5803 | -575.5379429 |
| exact `logLik` at those coefficients | -22023.06 | -575.5379429 |
| difference | **2.166e+04** | **1.93e-12** |
| `H` on the censored row | 21 690, scored as -35.1274 | 55.73, scored as -55.73 |
| optimizer | converged, no warning | reports false convergence (8) and a max gradient of 0.00123, and WARNS |

The last row is the honest part. With the floor gone the likelihood is
no longer flat in that direction, so the optimizer has to fight the
outlier instead of ignoring it, and it says so. A fit that says it did
not converge is better than one that converges to a floor.

`rp_floored()`'s censored count is a diagnostic now and never refuses;
the monotonicity floor still does, because no core seam addresses it.
`royston_parmar()` also declares `post$fit_check`, so a non-monotone fit
warns as it is returned.

### The extension's own suite

258 tests at 0.1.0, **265 at 0.2.0**, 0 failed, 0 error, 0 warning, 0
skipped, in one process.

## Two more pre-existing defects found in the path

Neither is in the contract; both blocked work that is.

### `eval_dpars(fit, b = NULL)` recursed on every model with a block

`?frmtmb-extension-api` documents `b` as "Random-effect vector to
evaluate at, defaulting to the fit's own estimates. `NULL` drops the
random-effect contribution." It did not drop it: the code still formed
`lp$Z %*% b` with `b` at `NULL`, and Matrix's S4 dispatch on a sparse
matrix times `NULL` recurses until R reports "evaluation nested too
deeply: infinite recursion".

Reproduced on the main-built frmtmb 0.50.0 in the user library with
`y ~ x + (1 | g)`, so it affected every model with a random-effect
block and not only the ones this lane fits. Fixed at
**R/predict.R:571**, one condition. `ps()`'s closures take the matching
branch: a penalized block IS a random effect, so dropping it leaves the
curve's null space, which is what `ps_env()` (**R/ps.R:251**) now
returns.

### `ps(k = length(t))` was silently accepted as `k = 1`

Mine, caught by my own test. `ps_spec()` evaluated `k`, `degree`, `pad`
and `center` in `baseenv()` to insist they were constants. `baseenv()`
holds `t`, the transpose function, so `length(t)` evaluates there to 1
and a data-dependent basis size came back as a length-one number that
passed every subsequent check. They must be LITERALS, and are
(**R/ps.R:311**).


### The zero-column fix has a third witness, and it is a behavior change

`tests/testthat/test-tmb-examples.R` fits
`z ~ 0 + ar1(tim + 0 | g)` inside a nonlinear body and passed
`start = list(beta = c(0, 2, 3))`. The leading `0` was the phantom
entry's starting value. With the phantom gone the vector is
`c(2, 3)`, and both tests were updated.

That is the user-visible edge of the fix and it is worth stating plainly:
a model with a purely-random nonlinear parameter loses one entry from
`beta`, so a hand-written `start$beta` for such a model changes length.
The fit does not move; the entry never entered the likelihood.

### `nl_body_vars()` and an empty argument

Caught by `test-nl-lexical.R` under as-cran and by nothing before it:
`m[, 1]` in a nonlinear body holds the missing-argument symbol, which
can be EXTRACTED without complaint and cannot be passed on. The error
fires at the callee, so the `tryCatch` around the extraction, copied
from `test-message-uniqueness.R`'s walker, caught nothing. The value has
to be FORCED inside the handler. Fixed at **R/parse.R:1329** and at the
same pattern in `ps_extract()`, **R/ps.R:291**, with four regression
assertions and a fitted matrix-column model in
`tests/testthat/test-nl-body-vars.R`.

## Touched files

Core, all in the lane's list except where noted:

* `R/parse.R` - `nl_body_vars()`, the `ps_extract()` hook in `nl_dpar()`,
  three `all.vars()` call sites.
* `R/ps.R` - NEW. The whole `ps()` term.
* `R/fit-end.R` - NEW. `fit_end_checks()`.
* `R/frame.R` - the `ps()` frame hook, the block resolution, the
  `cens()`/`lccdf` guard, the zero-column `paste()` fix.
* `R/objective.R` - the `lccdf` branch of the censored term, the `ps()`
  closures in `ev`.
* `R/predict.R` - `frm_joint_cov()`, `frm_lp_basis()` and its nonlinear
  branch, `predict(se.fit = TRUE)` rewritten as a consumer, the three
  shared argument checks, the `ps()` closures in `eval_dpars()` and in
  `predict(newdata =)`, the `b = NULL` fix.
* `R/families.R` - the `frmtmb_family()` constructor's `lccdf` argument
  and its documentation, `fam_lccdf()`, `has_lccdf()`, and the slot on
  five built-in families. The mixture family's response scale, which
  wt-ce owns, was not touched.
* `R/fit.R` - THREE self-contained hunks, all outside the prior
  application: `require_frmtmb_fit()` beside `require_fitted()`, one
  `check_ps_fit()` call beside `check_structure_fit()`, and one
  `fit_end_checks()` call after `check_convergence()`.
* `R/compat.R` - **outside the lane's explicit file list**, and named
  here for that reason. Part 1b of the contract asks for a
  `frm_lp_basis` method feature and Part 2 asks that "the compat table
  say what it does with it", and the vocabulary and the rules both live
  in this file. No sibling lane owns it.
* `NEWS.md`, `DESCRIPTION` (0.52.0; `brokenstick` and `frmtmb.spline`
  added to Suggests for the new case-studies section).
* `vignettes/case-studies.Rmd` - one new section, 13.
* `man/` - by roxygen, idempotent on a second pass.
* `tests/testthat/test-nl-body-vars.R`,
  `tests/testthat/test-predict-lp-basis.R`,
  `tests/testthat/test-cens-lccdf.R`, `tests/testthat/test-ps.R` - NEW.
* `tests/testthat/test-compat.R` - the specials vocabulary gains `ps`.
* `dev/spline-seam-proposal.md`, `dev/spline-core-findings.md`.

Extension, all under `extensions/frmtmb.spline/`:

* `DESCRIPTION` (0.2.0, `frmtmb (>= 0.52.0)`, `Matrix` dropped from
  Imports because the sparse fallback it was there for is gone),
  `NEWS.md`.
* `R/curve-cov.R` - `sp_curve_parts()` rewritten onto
  `frm_lp_basis()`; `sp_joint_cov()`, `sp_curve_design()`,
  `sp_check_linear()`, `sp_probe()` and `sp_coef_pos()` deleted.
* `R/curve.R` - the two documentation sections and the `print()`
  branch for an unchecked covariance.
* `R/royston-parmar.R` - the `lccdf` slot and `post$fit_check`.
* `R/rp-check.R` - the censored count demoted to a diagnostic.
* `R/zzz.R` - six compat rows.
* `tests/testthat/test-curve.R`, `test-rp-floored.R`, `test-surface.R`.
* `man/` - by roxygen.

Nothing under `extensions/frmtmb.ode`, `extensions/frmtmb.ddm`,
`extensions/frmtmb.latent` or `extensions/frmtmb.sample` was edited, and
no file a sibling lane owns (`R/priors.R`, `R/sampling-api.R`,
`R/conditional-effects.R`, `R/methods-fit.R`, `R/covstruct.R`) was
edited.

## Verification

One process per item unless stated. `NOT_CRAN = true` for the suites;
`_R_CHECK_CRAN_INCOMING_ = false` and pandoc on PATH for the checks.

### The tests this lane added

| file | assertions | result |
|---|---:|---|
| `tests/testthat/test-nl-body-vars.R` | 23 | 0 failed, 0 error |
| `tests/testthat/test-predict-lp-basis.R` | 30 | 0 failed, 0 error |
| `tests/testthat/test-cens-lccdf.R` | 27 | 0 failed, 0 error |
| `tests/testthat/test-ps.R` | 43 | 0 failed, 0 error |

### The files the contract audits by name

| file | assertions | result |
|---|---:|---|
| `test-nl.R` | 14 | clean |
| `test-nlf.R` | 74 | clean |
| `test-nl-lexical.R` | 31 | clean (was the ONE failure the first as-cran found) |
| `test-nl-rtmb-scope.R` | 103 | clean |
| `test-nl-body-vars.R` | 23 | clean |
| `test-predict-newdata.R` | see the full table | clean |
| `test-predict-lp-basis.R` | 30 | clean |
| `test-cens-trunc.R` | 35 | clean |
| `test-cens-lccdf.R` | 27 | clean |
| `test-message-uniqueness.R` | 6 | clean, after three duplicate templates were factored into shared helpers |
| `test-bracket-access.R` | 8 | clean |
| `test-compat.R` | 265 | clean, after the specials vocabulary gained `ps` |
| `test-ps.R` | 43 | clean |

`test-message-uniqueness.R` failed first, on three templates
`frm_lp_basis()` had copied from `predict()`: the `re.form` refusal, the
unknown-response refusal and the missing-newdata-variable refusal. They
are one template each now, in `check_re_form()`,
`stop_unknown_response()` and `stop_newdata_missing()`
(**R/predict.R:34-62**), which is what the test's own note prescribes
for a shared validation helper.

### The extension suites

| package | failed | error | warning | skipped | passed |
|---|---:|---:|---:|---:|---:|
| frmtmb.spline 0.2.0 (one process) | 0 | 0 | 0 | 0 | **265** |
| frmtmb.ddm 0.2.0 (one process) | 0 | 0 | 0 | 1 | 923 |
| frmtmb.ode 0.1.0 (one process) | 0 | 0 | 0 | 0 | 211 |

frmtmb.spline was 258 at 0.1.0. Neither sibling extension was edited;
both use nonlinear bodies and censoring, which is why the contract
names them, and both are clean against the new core.

### as-cran

| package | status |
|---|---|
| frmtmb 0.52.0 | **OK** (0 ERROR, 0 WARNING, 0 NOTE; suite inside the check 375 s, OK) |
| frmtmb.spline 0.2.0 | **OK** |

Three things the checks found that the per-file runs did not, all fixed:

1. `nl_body_vars()` died on an empty argument. `m[, 1]` in a nonlinear
   body holds the missing-argument symbol, which can be extracted and
   not passed on, so the error fires at the callee and the `tryCatch`
   around the extraction caught nothing. Under the FIRST as-cran the
   whole core suite was FAIL 1, PASS 4012, SKIP 187, and this was the
   one.
2. The `\tabular` block in the new `lccdf` documentation was malformed:
   the shell heredoc that wrote it ate a backslash, so `\tab` became a
   literal tab character. Every `R/` and test file this lane touched was
   then scanned for stray tabs; there are none.
3. `frmtmb.spline` still declared `Matrix` in Imports, which the sparse
   Schur-complement fallback used and which went with it, and
   `print.frmtmb_curve()` read a zero-length `cov_rel_error` on a
   SUBSET of a curve, where the attributes are gone.


## What was left out, and why

**Consumer (ii), shape-constrained smooths through exponentiated
coefficients.** Not implemented, and `ps()` has no `shape` argument, so
nothing is half done. Consumer (i) needed the parse pass, the frame
block, the tape closure and a nonlinear delta method; the coefficient
map this consumer adds is a fifth piece that touches the same four. The
proposal's design still reads correctly against what landed: the block
is the same block, the map belongs where `ps_coefs()` builds the
coefficient vector, and the Jacobian this consumer would need for
standard errors is already free, because `frm_lp_basis()` tapes the body
rather than assuming it is linear. That is the piece that was expensive
and it is done.

**Consumer (iii), spline-valued inputs for `frmtmb.ode`.** Deferred by
the contract. Nothing under `extensions/frmtmb.ode` was touched, and the
note for that package's owner is in the proposal. What changed for them
is that the capability exists: `ps()` evaluates its basis at an
arbitrary tape-valued argument, measured against `splineDesign()` with
an AD input.

**`ps(by = )` and `ps(id = )`, the factor-smooth surface.** Refused by
name rather than accepted, because a `by =` factor smooth is several
blocks with a shared variance and none of that was measured.
`ps(penalty = )` likewise: the penalty is the second difference and
saying so is better than an argument with one legal value.

**REML, quadrature, profiling and `mvbf()` with `ps()`.** Refused by
name. The first three integrate out something a `ps()` block has already
put inside the Laplace approximation. `mvbf()` is untested rather than
known to be wrong, and is refused for that reason.

**The importance correction with `ps()`.** Refused, but not by anything
this lane wrote: core already refuses EVERY nonlinear predictor, because
a body mixes parameter values with raw data columns and the corrected
objective evaluates it once per draw. The proposal asked for this row to
be marked `works` only after it was run; it cannot be run, so it is
`refused` with the inherited reason named.

**`emmeans::emm_basis()` rewritten in terms of `frm_lp_basis()`.** Not
done. The proposal lists it as an eventual ask rather than an immediate
one, and it is a separate behavior change: it would make a smooth term
visible to `emmeans` on a frmtmb fit, which is a feature and not a
seam.

**Left truncation.** `lccdf` fixes right censoring and the proposal says
so. Left censoring, interval censoring and the truncation normalizer are
unchanged, so a left-truncated survival model meets the same
representability problem from the other side. Closing it needs a
windowed log-difference slot.

**Simultaneous coverage in the recovery study.** The paper reports
simultaneous coverage of the population curve; the table above reports
POINTWISE coverage, which is a different quantity. The simultaneous
critical value needs `frmtmb.spline::frm_curve()`, which now works on a
`ps()` model (`vignette("case-studies")` section 13 reports 2.5683 on
the SMOCC fit), but running 300 replications of it at the paper's 16
settings is a study rather than a check.

**The refusal wording the proposal asks about at the end of Part 3.**
Not changed. With the walker fixed, "object 'tv' not found" is no longer
reached by the case that prompted the suggestion, and rewording a
message on the strength of a case that can no longer produce it is how a
message stops describing what happened.

### The full core suite, one file per process, audited by name

113 test files at the end of this lane, 109 at 2210aa1 plus the four it
adds. Every file was run in its own process, and the file list was
compared against `ls test-*.R` rather than counted:

```
expected 113   ran 113   missing 0   extra 0   duplicates 0
files 113 | failed 0 | error 0 | warning 0 | skipped 87 | passed 6255
```

The audited names:

| file | passed | skipped | result |
|---|---:|---:|---|
| `test-nl.R` | 14 | 0 | clean |
| `test-nlf.R` | 74 | 0 | clean |
| `test-nl-lexical.R` | 31 | 0 | clean |
| `test-nl-rtmb-scope.R` | 103 | 0 | clean |
| `test-nl-body-vars.R` (new) | 23 | 0 | clean |
| `test-predict-newdata.R` | 12 | 0 | clean |
| `test-predict-lp-basis.R` (new) | 30 | 0 | clean |
| `test-cens-trunc.R` | 35 | 0 | clean |
| `test-cens-lccdf.R` (new) | 27 | 0 | clean |
| `test-message-uniqueness.R` | 6 | 0 | clean |
| `test-bracket-access.R` | 8 | 0 | clean |
| `test-compat.R` | 265 | 0 | clean |
| `test-ps.R` (new) | 43 | 0 | clean |
| `test-tmb-examples.R` | 28 | 0 | clean |

One process-hygiene note, because it produced a wrong number before it
produced a right one. An earlier run of this loop was stopped and its
shell survived the stop; a second run then appended to the same output
file and the table came back with 87 of 113 files and a contiguous
alphabetical hole in the middle. That was two writers, not a test
result. The surviving processes were killed by matching `sc-coresuite.sh`
and `sc-onefile.R` on the command line, which is this lane's own prefix,
and the suite was re-run alone to a fresh file. The table above is that
run.

---

# Punch round, 2026-09-05

Against `dev/review-spline-core.md` (GO WITH FIXES, ten items). The
reviewer's own two edits, `R/ps.R:96-100` and `man/ps.Rd:101-105`, were
kept: the roxygen source they corrected is the source `ps.Rd` is
generated from, so re-running roxygen reproduces their text rather than
reverting it. Verified after regeneration at **man/ps.Rd:116-120**, from the roxygen
at **R/ps.R:111-115**.

## 1. `?ps` claimed the importance correction applies. FIXED BY THE REVIEWER.

Nothing to do but keep it and prove roxygen does not undo it. It does
not: `?ps` and `R/compat.R:1123` now say the same thing, which is that
the correction is refused for every nonlinear predictor and not by
anything `ps()` declares.

## 2. `predict(newdata =)` past the knot span was silent. FIXED.

`ps_span_warning()` at **R/ps.R:553**, wired at **R/predict.R:1256-1260**.

The check goes inside the closure rather than beside it, because the
closure is the only place the evaluated argument exists: the expression
may name nonlinear parameters, so its value is not known until the body
has been walked. `is.numeric()` keeps it off the tape, and it is armed
only when `newdata` is not `NULL`, because in sample the fit-end report
has already said it and repeating it on every `fitted()` call would be
noise.

Measured, on a `pad = 0.3` fit whose `knot_range` is [-0.271, 1.271]:

| grid | before | after |
|---|---|---|
| inside the span | silent | silent |
| 2 of 3 rows outside | silent | warns, naming both ends of the span and the count |
| far past the last knot | silent, returns exactly `lev` | warns; still returns exactly `lev`, which is the honest value |

The arithmetic did not change and was never wrong. What changed is that
a curve which bends smoothly to its intercept now says why.

Pinned in `tests/testthat/test-ps.R`, "predict() at newdata past the
knot span says so": the warning fires outside and not inside, the
message carries both ends of the span and the row count, the far-field
prediction equals `fixef(fit)$lev` to 1e-8, and `predict(fit)` and
`fitted(fit)` stay silent.

## 3. `fit_end_checks()` ran a family's hook with no guard. FIXED.

**R/fit-end.R:18-44**. Both calls are wrapped: the family's
`post$fit_check` and core's own `ps_coverage_warning()`, which
evaluates a nonlinear body and is arbitrary code for the same reason.

A failing check degrades to a warning naming the family and the hook,
and the fit is returned. The reviewer confirmed the old behavior
destroyed a completed fit; the new behavior is pinned four ways in the
new file `tests/testthat/test-fit-end-hook.R`:

* a hook that runs sees a FINISHED fit (`fit$opt` present), which is the
  only thing the slot exists for;
* a hook that warns reaches the user and the fit survives;
* a hook that THROWS leaves a `frmtmb_fit` whose `logLik()` is finite
  and identical to the same model fitted with no hook at all
  (tolerance 1e-10), with a warning naming the family (`hooked`), the
  hook (`post$fit_check`), the underlying message and "complete and
  unaffected";
* the guard is exercised directly on a real fit with a throwing hook
  planted on its family, and `fit_end_checks()` returns rather than
  propagating.

`post$fit_check` is now documented in `?frmtmb_family`'s `post`
argument, including the instruction to warn rather than stop.

## 4. NEWS filed the `start$beta` change as a plain FIX. FIXED.

**NEWS.md:93-111**. Relabelled `BEHAVIOR CHANGE, and a fix`, and the
sentence the reviewer objected to ("Affected models fitted to the same
optimum; only the inference was lost", which describes the old bug
rather than the new break) is replaced by a bullet that says what
breaks, what to do, and what does not move: `beta` loses one entry per
purely-random nonlinear parameter, a hand-written `start$beta` for such
a model now errors, drop the element that stood for the parameter with
no fixed design, `par_template(fit)` names what is left, and the
optimum, the log-likelihood and every coefficient are unchanged.

## 5. The beta0 explanation attributed a convention the paper does not state. FIXED.

**vignettes/case-studies.Rmd** and this document. The paper says only
that it applies a sum-to-zero constraint "to ensure identifiability of
an intercept term" and never says on what, so the attribution is gone
and a measurement is in its place. Re-derived here rather than taken
from the review:

```
beta0 (ps sum-to-zero on the coefficients) : 68.754333
mean(fitted curve at the data)             : -0.55609575
beta0 + mean(curve)                        : 68.198237
the paper's beta0                          : 68.2
difference                                 : -0.0018
```

The vignette now runs that arithmetic in a chunk of its own
(`smocc-beta0`), so the claim is computed on the page rather than
asserted. One detail worth keeping that the review does not mention: the
centering is specific to the data's own age distribution. On an even 0
to 130 week grid the same arithmetic gives 76.61, because the children
are not evenly spread over age, so "re-centred to the fitted values" is
the statement and "re-centred" alone is not.

## 6. `ps()`'s knot rule is not the paper's, and the text implied it was. FIXED.

**R/ps.R:71-84** and **vignettes/case-studies.Rmd**. `?ps` now says
plainly that this is not a reimplementation of the paper's section
2.3.1, what the paper does (rescale the argument to [0, 1] between
bounds that move with the current variance estimates, so its knots
follow the fit), what `ps()` does instead, and the two reasons: moving
knots make the penalty matrix and its eigendecomposition functions of
the parameters, so the fixed-random split would have to be redone inside
the objective; and a basis rebuilt between fitting and prediction is a
different basis, which is `poly()`'s rule and the reason `pad =` exists.

## 7. The derived compat note was wrong for `ps()`. FIXED.

**R/compat.R:704**. A `kind:family` x `ps()` rule now overrides the
`kind:special` default, whose sentence ("Predictor specials build design
columns before the family sees them") is false of the one special that
builds no design columns.

Measured after the change: **36 of 36** `family x ps()` rows carry the
new note and **0** carry the old one, and `gaussian x s()` is untouched.
The status is `works` on every row, as before; only the reason moved.
The review says 68 rows; the vocabulary has 36 family features, so 36 is
the count here.

## 8. `use_re` was duplicated verbatim. FIXED.

`re_form_keeps()` at **R/predict.R:62**, replacing the identical three
lines at what were `R/predict.R:1085-1087` and `R/predict.R:2952-2954`.
It sits beside the three message helpers factored out for the same
reason, and the reason is worth stating: section 3 of the review
measures `predict()` and `frm_lp_basis()` returning the same numbers,
which they can do only while they read `re.form` the same way.

## 9. A third `n_predict` assertion described deleted machinery. FIXED.

**extensions/frmtmb.spline/tests/testthat/test-curve.R**. The test named
"a grouping block costs one probe, not one call per level" asserted
`expect_lt(n_predict, 30L)` under a six-line comment about chunks of 24
and a straddled boundary. There are no chunks; `n_predict` is the
constant `1L`, so the assertion passed vacuously and the comment
documented machinery this lane deleted.

Rewritten as "a grouping block costs nothing, because there is no
probe": `expect_equal(n_predict, 1L)` on a 40-level grouping block plus
a `k = 9` smooth, with the comment saying what the test used to pin and
why the count no longer responds to the block at all. The extension
suite goes from 265 to **266**.

## 10. `Suggests` out of alphabetical order. FIXED.

**DESCRIPTION:41**. `brokenstick` now sits after `bayesplot`.

## Verification after the punch round

Every file this round touched, plus the set the contract audits by name,
one process each.

| file | passed | result |
|---|---:|---|
| `test-nl.R` | 14 | clean |
| `test-nlf.R` | 74 | clean |
| `test-nl-lexical.R` | 31 | clean |
| `test-nl-rtmb-scope.R` | 103 | clean |
| `test-nl-body-vars.R` | 23 | clean |
| `test-predict-newdata.R` | 12 | clean |
| `test-predict-lp-basis.R` | 30 | clean |
| `test-message-uniqueness.R` | 6 | clean |
| `test-compat.R` | 265 | clean |
| `test-ps.R` | 51 | clean (43 before this round; the knot-span pin adds 8) |
| `test-fit-end-hook.R` (NEW) | 18 | clean |
| `test-cens-lccdf.R` | 27 | clean |
| `test-cens-trunc.R` | 35 | clean |
| `test-bracket-access.R` | 8 | clean |

**0 failed, 0 error, 697 passed** over the fourteen.

`test-message-uniqueness.R` passing matters here specifically: this round
adds three warning templates (the knot-span warning, the coverage-report
guard and the hook guard), and they are distinct from each other and from
the fit-end coverage warning they sit beside.

frmtmb.spline, one process: **266** passing, 0 failed, 0 error, 0
warning, 0 skipped. It was 265 before item 9 rewrote the stale test.

as-cran, both packages, `_R_CHECK_CRAN_INCOMING_=false`, pandoc on PATH:

| package | status |
|---|---|
| frmtmb 0.52.0 | **OK** (0 ERROR, 0 WARNING, 0 NOTE; suite inside the check OK) |
| frmtmb.spline 0.2.0 | **OK** (same) |

Roxygen is idempotent after the round, and the reviewer's `?ps` text is
what it regenerates.

The vignette's new `smocc-beta0` chunk was knitted before it was
committed to the page: it prints 68.75433, 68.19824, 68.20000. No
knot-span warning fires anywhere in the section, which is the right
outcome, because the 0 to 130 week grid is inside the padded span
[-34.93, 174.64].
### The full core suite after the punch round, audited by name

114 test files now, 113 before this round plus `test-fit-end-hook.R`.
Every file in its own process, the list compared against `ls test-*.R`
rather than counted:

```
expected 114   ran 114   missing 0   extra 0   duplicates 0
files 114 | failed 0 | error 0 | warning 0 | skipped 87 | passed 6281
```

6255 before the round over 113 files; the 26 added are the 18 of
`test-fit-end-hook.R` and the 8 the knot-span pin adds to `test-ps.R`.

No stray `sc-` loop was running when this started, checked first,
because the collision that produced a wrong count last round came from
one that was.
