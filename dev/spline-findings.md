# frmtmb.spline: findings

Worktree `C:\Users\adf44\source\r\frmtmb-wt-spline`, branch `wt-spline`,
base commit 5dfdd84 (core 0.50.0). Nothing here is committed.

## 1. The covariance seam, read before anything was written

The task named two exported routes to check. Both were read end to end.

### `vcov(fit, full = TRUE)` (`R/methods-fit.R:377`)

Returns the covariance of the OUTER parameter vector. Two branches:

* ML, no profiling: `sdr_of(object)$cov.fixed`, renamed to
  `outer_par_names(object)`. The outer vector is what `nlminb` optimized,
  which under the Laplace approximation is `beta`, `betad`, `theta`,
  `thetaac`, `thetar`. The random-effect vector `b` is NOT in it: it is
  integrated out, so it never was an outer parameter.
* REML or `control$profile`: the blocks come from
  `sdr_of(object)$jointPrecision`, and the code then explicitly drops
  `b`: `comps <- setdiff(names(par_template), c("b", "miss", "beta"))`.

So `full = TRUE` never returns a `b` row under either branch. It cannot:
its documented invariant (`test-methods-audit.R`) is that its row names
are exactly `confint()`'s.

### `predict(fit, se.fit = TRUE)` (`R/predict.R:1019`)

This one DOES use the joint covariance, and does it correctly:

    jc <- get_joint_cov(object)                       # R/predict.R:12
    da <- lp_delta_A(object, lp, ed, newdata, ...)
    V  <- jc$V[da$coef_pos, da$coef_pos, drop = FALSE]
    var_eta <- pmax(rowSums((A %*% V) * A), 0)

`get_joint_cov()` solves `sdr_of(fit)$jointPrecision` (falling back to
`autoscale_sdreport(fit, jp = TRUE)` and then to `cov.fixed`), so `V`
covers `beta` AND `b` and their cross-block. `lp_delta_A()` builds the
per-row design row `A` over `coef_pos`, which for a smooth spans the
null-space columns in `X` and the penalized columns in `b`.

What comes back is `sqrt(diag(A V A'))`, one number per row. The
off-diagonal of `A V A'` is formed and thrown away inside `rowSums()`.

### Every other exported candidate, checked

| route | what it gives | why it is not enough |
|---|---|---|
| `vcov(fit)` | `beta`, `betad` block | no `b` |
| `vcov(fit, full = TRUE)` | outer parameters | no `b`, by construction |
| `predict(se.fit = TRUE)` | pointwise SE of `eta` | diagonal only |
| `ranef(condVar = TRUE)` | `sqrt(diag.cov.random)` | diagonal of the `b` block only, no cross-block with `beta` |
| `emmeans::emm_basis` | `X`, `bhat`, `V = vcov(fit)[idx, idx]` | fixed block only; a smooth's wiggly columns are absent |
| `marginaleffects::get_vcov` | `vcov(fit)` | same |
| `insight::get_varcov` | `vcov(fit)` | same |
| `lme4::getME(fit, c("X","Z","Zt","b"))` | designs and modes IN SAMPLE | no covariance, and no newdata design |
| `hypothesis()` | delta method over `beta, betad, theta, thetaac, thetar` (`hyp_par_cov()`, `R/confint.R:1660`) | `b` excluded there too |
| `marginaleffects::set_coef` | sets `beta` and `betad` | cannot set `b`, so the design cannot be recovered by unit perturbation |
| `frm_bootstrap(FUN =)` | nsim refits, any numeric functional | DOES reach the curve, by refitting; see below |

**Conclusion.** No exported route yields the joint covariance of a grid
prediction. Exactly two things are missing and both are one line of core:

1. the joint covariance itself (`get_joint_cov()` is `@noRd`), and
2. the grid design over the same coefficient positions
   (`lp_eta_design()` + `lp_delta_A()` are `@noRd`).

`predict(se.fit = TRUE)` already computes both and discards the second.
The seam goes in `dev/spline-seam-proposal.md`.

`frm_bootstrap()` is the one exported route that reaches a full curve
covariance, by paying for nsim refits.

## 2. The workaround, and why it is safe

`predict()` is linear in the coefficient vector `c = (beta, b)`:
`eta(grid) = C c`. So column `j` of `C` is exactly

    predict(fit_with_c_j_plus_1, newdata = grid) - predict(fit, newdata = grid)

with no truncation error, because the difference of a linear function is
its own derivative. `fit$estimates` is the slot the documented
extension-API example already reads
(`dpar_linpred(fit$frame, fit$estimates, ...)`), and `fit$obj` is the
slot `frmtmb.sample::frm_sample()` already reads (`sample.R:1271`), so
neither is a `:::` reach.

`V` comes from `RTMB::sdreport(fit$obj, getJointPrecision = TRUE)`,
which is what `autoscale_sdreport()` itself calls when `par_units` is
absent (`R/autoscale.R:137`).

**The construction proves itself.** `predict(se.fit = TRUE)` is the one
exported route that overlaps: it returns `sqrt(diag(C V C'))`. So the
package computes its own `diag(C V C')` and compares. Measured on
`y ~ s(x, k = 10)`, gaussian, n = 300, 60-point grid:

    max abs difference  5.90e-16
    max rel difference  1.02e-14

That is machine precision. `C` and `V` are the same `C` and `V` core
uses. The check runs on every `frm_curve()` call and refuses rather than
returning a covariance it cannot vouch for, so an autoscaled or
profiled fit whose joint precision this route reconstructs wrongly is
caught rather than believed.

Cost of the workaround: one `predict()` call per candidate coefficient,
plus one probe per chunk of 24 that contributes nothing. On this model
that is 11 calls in about 0.04 s, and `RTMB::sdreport()` a further
0.04 s. On the vignette's 20-subject factor-smooth model, which carries
110 random coefficients, it is 32.

## 3. Measured against mgcv and gratia

frmtmb reproduces the mgcv ML fit (the core case-studies recipe):

    logLik(frm)      -169.5139
    -gam$gcv.ubre    -169.5139
    max |eta_frm - eta_mgcv| on the grid   5.72e-07

The curve standard errors do NOT match mgcv's `Vp` exactly, and should
not:

    frmtmb vs mgcv Vp                      max rel  2.93e-02
    frmtmb vs mgcv vcov(unconditional=TRUE) max rel 1.37e-02

TMB's joint precision is over `(b, beta, betad, theta)` together, so
inverting it and taking the `(beta, b)` block already carries the
smoothing-parameter uncertainty. mgcv's `Vp` does not; its
`unconditional = TRUE` correction (Wood 2016) adds a first-order
version of the same thing, and lands 1.4 percent away. frmtmb is
nearer the unconditional matrix than the conditional one, which is the
expected ordering.

### The simultaneous critical value, algorithm identity with gratia

gratia's `confint.gam(type = "simultaneous")` standardizes the
SMOOTH-ONLY deviation `Cg %*% (b* - b)` by `.se` from
`smooth_estimates()`, and that `.se` is the FULL linear predictor's
standard error, intercept column included (verified: gratia `.se` equals
`sqrt(rowSums((Xp %*% Vp) * Xp))` to the printed digit, and differs from
the smooth-only `sqrt(rowSums((Cg %*% Vp[cs,cs]) * Cg))` by up to 13.5
percent). The band is still an exact 95 percent simultaneous band,
because the same `.se` divides in the calibration and multiplies in the
band, but the critical VALUE is not comparable across two packages
unless the divisor is the same.

Feeding the same `Cg`, the same `Vp` block and gratia's own divisor to
this package's max-deviation simulation reproduces gratia:

    nsim     mine     mcse      gratia    difference
    10000    2.74068  0.01344   2.73337   0.54 mcse
    100000   2.73762  0.00448   2.74590   1.85 mcse
    500000   2.74001  0.00198   2.73683   1.60 mcse

Every difference is inside two Monte Carlo standard errors. The Monte
Carlo SE is the quantile SE, `sqrt(p(1-p)/n) / f(q)` with `f` a kernel
density estimate at the quantile.

On the self-consistent quantity, the whole linear predictor on the grid
standardized by its own standard error, the critical values are

    mgcv Vp                 2.96054  mcse 0.00334
    mgcv unconditional      2.97196  mcse 0.00329
    frmtmb joint precision  2.96874  mcse 0.00340

at nsim = 200000, so the covariance flavor moves the critical value by
about 0.011, roughly 3 Monte Carlo standard errors, and the two
frequentist-with-theta-uncertainty matrices agree to under one.

## 4. Central differences: the measurement that sets `eps`

gratia's `derivatives()` uses a fixed `eps = 1e-7`. Swept against the
exact derivative of `f(t) = 2 sin(pi t) + 0.6 t` on 41 points of
`[0.05, 0.95]`:

    eps     max |d1 error|   max |d2 error|
    1e-02   1.021e-03        1.623e-03
    1e-03   1.021e-05        1.624e-05
    1e-04   1.021e-07        2.234e-07
    1e-05   1.040e-09        9.235e-06
    1e-06   4.033e-10        8.219e-04
    1e-07   5.660e-09        8.165e-02
    1e-08   4.652e-08        1.007e+01
    1e-10   2.925e-06        8.884e+04

First order: truncation falls as `eps^2` and cancellation rises as
`u/eps`, meeting near `eps = 1e-6` at 4e-10. gratia's 1e-7 is 14 times
worse than the optimum and still negligible.

Second order: the second difference divides by `eps^2`, so cancellation
rises as `u/eps^2`. The optimum is near `eps = 1e-4` at 2.2e-07, and
gratia's 1e-7 is wrong by **8.2e-02**, eight orders of magnitude worse.
A fixed `eps` cannot serve both orders. This package therefore scales
`eps` with the covariate range and sets it per order: `1e-6 * range` for
the first derivative and `1e-4 * range` for the second.

## 5. Derivatives, measured against gratia and against mgcv

Same model, same 25-point grid, `frm_curve_deriv(order = 1)` against
`gratia::derivatives(type = "central")` on the mgcv ML fit:

    estimate, max abs difference                      7.07e-06
    estimate vs mgcv's own same-stencil derivative    7.10e-06
    standard error, max rel difference                7.11e-02

The estimate difference is the FIT difference amplified by
differentiation: the two fitted curves differ by 5.7e-07 and dividing a
difference of that size by `2 eps` is what 7e-06 is. Forcing gratia onto
this package's step size moves it by 3e-08, so the step size is not the
cause. The standard errors differ by 7.1 percent for the same reason the
curve's do: mgcv's `Vp` is conditional on the smoothing parameter and
the joint precision is not, and differentiation amplifies a 2.9 percent
difference into 7.1.

Order 2 is where the fixed step matters:

    frmtmb (eps = 1e-4 * range) vs mgcv same-stencil   8.62e-05
    gratia (eps = 1e-7) vs frmtmb                      1.76e+00

At the first grid point gratia reports -5.1244 and this package -3.3677,
a 52 percent difference on a second derivative of order 5. The mgcv
same-stencil calculation agrees with this package, so the disagreement
is gratia's cancellation, not this package's. That is the whole argument
for scaling `eps` by the covariate range and setting it per order.

## 6. Royston-Parmar: the log-likelihood identity with flexsurv

`flexsurv`'s `bc` data (686 rows, 299 events, 387 right censored). For
each `scale` and each knot count, flexsurv is fitted first, its knots
are handed to `royston_parmar()`, the frmtmb objective is built with
`frm(dry_run = "objective")` and evaluated at flexsurv's OWN fitted
coefficients. The two log likelihoods are then the same number:

    scale    k   flexsurv         frmtmb at its par   rel difference
    hazard   0   -811.9419254350  -811.9419254350     7.00e-16
    hazard   1   -792.8637978132  -792.8637978132     8.60e-16
    hazard   3   -790.0060434492  -790.0060434492     1.15e-15
    odds     0   -799.0522855601  -799.0522855601     7.11e-16
    odds     1   -788.9819017983  -788.9819017983     4.32e-16
    odds     3   -786.7919335307  -786.7919335307     1.59e-15
    normal   0   -790.7650654392  -790.7650654392     1.44e-16
    normal   1   -785.6232568284  -785.6232568284     1.30e-15
    normal   3   -783.4023677785  -783.4023677785     7.26e-16

Every one is inside ten units in the last place of a number near 800.
This is an identity, not an agreement: the two packages evaluate the
same function.

The two OPTIMA then agree as well, and frmtmb's is never the worse of
the two:

    scale    k   frmtmb optimum minus flexsurv's   max abs coef difference
    hazard   0   +1.89e-05                          6.65e-04
    hazard   1   +6.48e-09                          1.55e-05
    hazard   3   +1.12e-06                          3.32e-04
    odds     1   +8.38e-11                          2.54e-06
    normal   3   +5.60e-06                          8.83e-04

The coefficient differences are the optimizer's, not the model's: they
are the flat directions of a log likelihood the two packages agree about
to 1e-15.

## 7. The floor that was wrong in double precision

The obvious smooth positive part,

    sp_floor_pos(u) = 0.5 * (u + sqrt(u * u + eps2))

is right on paper and WRONG in floating point on exactly the side it
exists for. At `u = -35.75` and `eps2 = 1e-14`, `sqrt(u^2 + eps2)`
rounds to `|u|` because `1e-14` is lost against `1278`, and the sum
cancels to exactly zero, so `log()` returns `-Inf` after all.

Measured before the fix, on flexsurv's `bc` data at a parameter point
the optimizer visits: **647 of 686 rows** returned `-Inf` from the log
density, and 150 of 400 random parameter vectors gave a `NaN` objective
or gradient. Two of the package's own fits emitted "NA/NaN function
evaluation" while converging.

The first repair attempt was also wrong, and for a related reason. Using
the identity `(u + s)(s - u) = eps2` to get a cancellation-free
expression on each side and then COMBINING them,

    0.25 * ((w + r) + sign(u) * (w - r))     with w = |u| + s, r = eps2/w

fails because adding `r = 1.4e-16` to `w = 71.5` loses `r` entirely, so
both halves round to `w` and the subtraction cancels again. The measured
count went from 647 to 647.

The form that works SELECTS rather than combines:

    0.25 * ((1 + sign(u)) * w + (1 - sign(u)) * (eps2 / w))

`sign()` is available on an RTMB advector and differentiates to zero,
which is the right derivative here because the two branches agree at
`u = 0`, where both are `sqrt(eps2)`, so the function is continuous.

After the change: 0 rows of `-Inf`, 0 `NaN`, both optimizer warnings
gone, and the flexsurv identity unchanged at 1.4e-16 to 1.6e-15
relative. The measurement is what found this; the paper formula looked
correct at every reading.

A second floor was needed for the same class of reason. frmtmb forms a
right-censored contribution as `log(Fub - F(y))`, which is `log(1 - F)`,
so a family that returns `F == 1` returns `-Inf`. Measured: a slope of
30 in log time puts 572 of 686 rows at `F == 1` exactly. `F` is now
squeezed into the open interval with the same smooth positive part at
each end, exact in double precision wherever the distance to the
boundary exceeds 6.7e-08. That is a workaround, and the core fix, an
`lccdf` slot, is in the seam proposal.

## 8. The branch-free B-spline basis, measured

For the seam proposal, and run rather than asserted. A B-spline is a
divided difference of truncated powers, and the truncated power is
`0.5 * (e + abs(e))`; `abs()` works on an advector and `pmax()` does not
("Comparison is generally unsafe for AD types"), and RTMB exports no
`CondExp`, so this is the only construction available.

Cubic (`ord = 4`), 8 interior knots on [0, 1.5], 37 evaluation points,
against `splines::splineDesign`:

    value, plain numeric        max abs error   6.82e-14
    value, plain numeric        max rel error   3.64e-11
    value, x an AD parameter    max abs error   4.64e-14
    d/dx by reverse mode, vs
      splineDesign(derivs = 1)  max abs error   1.52e-13

The relative figure is the larger only because basis values near a
boundary knot are near zero. No tape comparison is needed anywhere and
the reverse-mode derivative IS `splineDesign`'s derivative to 1.5e-13.
The one caveat the measurement also shows: divided differences need
DISTINCT knots, so the repeated boundary knots `splineDesign` normally
takes have to be replaced by distinct knots outside the data range.

## 9. The parser defect, reproduced on 5dfdd84

`all.vars()` drops the arguments of a call in FUNCTION POSITION:

    all.vars(quote(f(x)(y)))            # "y"           -- x is gone
    all.vars(quote(a * curry(tv)(zv)))  # "a" "zv"      -- tv is gone
    all.names(quote(f(x)(y)))           # "f" "x" "y"   -- all.names sees it

Core collects a nonlinear body's data variables with `all.vars()`, so a
variable appearing only inside a function-position call is never
collected. Run on 5dfdd84 with `curry <- function(u) function(v) u * v`:

    a * curry2(tv, zv)           -> fitted, a = 2.006
    a * curry(tv)(zv)            -> REFUSED: object 'tv' not found
    a * curry(tv)(zv) + 0 * tv   -> fitted, a = 2.006

The third line is the diagnosis: naming `tv` where `all.vars()` can see
it makes the identical model fit to the same coefficient. The fix and
the message wording are in `dev/spline-seam-proposal.md`, Part 3. Note
that the refusal reads "object 'tv' not found" rather than the "not used
in the model formula" the brief anticipated; the mechanism is the one
described, the wording is not.

## 10. Ownership, confirmed rather than asserted

`git -C . diff --name-only 5dfdd84` returns exactly two files:

    README.md          (one table row plus one word, four/five)
    dev/build-docs.R   (one subsite entry plus the count in its comment)

Untracked: everything under `extensions/frmtmb.spline/**`,
`.github/workflows/check-frmtmb-spline.yaml`, `dev/spline-findings.md`
and `dev/spline-seam-proposal.md`. **Zero files under the root `R/`**,
by

    (git diff --name-only 5dfdd84; git ls-files --others --exclude-standard) | grep -E "^R/" | wc -l
    0

## 11. The CI workflow, compared line by line

`.github/workflows/check-frmtmb-spline.yaml` against
`check-frmtmb-ddm.yaml`, with the package name substituted out on both
sides so the comparison is of structure and not of spelling: the two
differ in the header comment ONLY. Every path filter, job name, step,
action version and argument is identical, including the
`R CMD INSTALL --no-multiarch --with-keep.source .` of the core from the
checkout and the DESCRIPTION-reading dependency step that skips
anything already installed. The comment differs because the dependency
list does (flexsurv, gratia, mgcv, survival rather than RWiener,
numDeriv, brms) and because it now says why the two heavy suggests are
installed rather than skipped: a run that skipped them would pass while
asserting nothing about the two comparisons this package exists to make.

## 12. Verification, counted by name

### The package suite, one process, `NOT_CRAN = true`

    file                        expectations   tests
    test-curve.R                          39      12
    test-deriv.R                          37      10
    test-gratia.R                         12       5
    test-message-uniqueness.R              4       1
    test-royston-parmar.R                 49       9
    test-surface.R                        46       9
    TOTAL                                187      46

    pass 187   fail 0   error 0   skip 0

Zero skipped: flexsurv, gratia, mgcv, survival and tinyplot are all
installed in the private library, so every `skip_if_not_installed()`
falls through and every comparison actually runs. The 42 tests by name
are in `sp-logs/named.log` under the session scratchpad; the ones
carrying the headline claims are

* `the assembled covariance reproduces predict(se.fit) exactly` (4)
* `the max-deviation simulation is gratia's, on gratia's own inputs` (3)
* `the objective at flexsurv's coefficients IS flexsurv's log likelihood`
  (18, being three scales times three knot counts times two assertions)
* `the second derivative is the one gratia's fixed eps gets wrong` (2)
* `the monotonicity floor keeps the log density finite` (3)
* `a reduced-rank block is caught by the check, not by the probe` (2)
* `a curve on a dpar other than mu finds its coefficients` (5)

Note that two of those carry `skip_on_cran()`, so the count under
`R CMD check --as-cran` is lower by design: the flexsurv identity sweep
fits nine models and the gratia comparison simulates 100000 draws, and
neither belongs in a CRAN run. They are what the developer suite is for.

## 13. The `rr` row was wrong, and the check is what said so

The compatibility table first claimed `frm_curve x rr` was `refused`,
with the reason that a reduced-rank block makes the predictor nonlinear
in its coefficients and the linearity probe catches it. Both halves were
wrong, and running it is what found that out.

A reduced-rank block's loadings live in `theta`, not in `b`
(`expand_b()`, `R/covstruct.R:1648`: `cvec[c_idx] <- as.vector(L %*% Fm)`
with `L` built from `theta`). So `eta` IS linear in `b` at fixed
`theta`, the linearity probe passes, and the construction goes ahead.
What it cannot see is the derivative with respect to the loadings, which
core's own delta method carries separately as `rr_jacobians()`
(defined at `R/predict.R:587`, used at 1230 and 1549).

Measured on `y ~ s(x, k = 6) + rr(0 + sp | site, d = 2)`, 5 species and
30 sites:

    re.form = NA     works, cov_rel_error  1.18e-14
    re.form = NULL   REFUSED, standard errors 27 percent away

`re.form = NA` works because the reduced-rank block drops out of a
population prediction entirely, so there is nothing to miss.
`re.form = NULL` is refused by the COVARIANCE CHECK, not by the probe.

The row is now `conditional` with that note, and there is a test. This
is the case the check was written for, and it is worth recording that it
fired on a case its author had reasoned about incorrectly.

The refusal message on the linearity probe was corrected at the same
time: it named `rr` as the thing it catches, and what it actually
catches is a nonlinear (`nl = TRUE`) body.

## 14. A second bug the check caught, in this package's own design

Reading a curve off a dpar other than `mu` refused, with the covariance
disagreeing by a relative 1. The cause was in this package: every dpar
except the location one keeps its coefficients in `betad`, not in
`beta`, and the design loop covered `beta` and `b` only. So the fixed
part of a `sigma` or `gamma1` curve was missing from `C` entirely.

Measured before and after, on
`bf(recyrs | cens(censored) ~ group, gamma1 ~ s(age, k = 5))` with
`royston_parmar(df = 2)`:

    before   REFUSED, covariance disagrees by 1.00 relative
    after    works,   cov_rel_error 2.89e-15

Two details the fix had to get right, both of which would have been
silent errors rather than refusals:

* A `betad` entry held FIXED (a dpar pinned to a constant, as
  `wiener(bias = 0.5)` pins one) moves the prediction when perturbed and
  has NO row in the joint precision, so including it would build an `A`
  wider than `V`. Fixed entries are excluded before the probe runs, and
  there is a test with `sigma = 0.5`.
* The joint precision numbers `betad` rows over the ESTIMATED entries
  only, while the design indexes the full `betad` vector. The two differ
  by however many are fixed, so the index has to be mapped through
  `setdiff(seq_along(betad), betad_fixed_idx)`. Every other component is
  one to one.

Adding `betad` to the loop costs one extra `predict()` call on a model
whose only `betad` entry is a `sigma` intercept that cannot move the
`mu` curve: the column comes back all zero and is dropped, and 10 live
columns now cost 11 calls. That is the price of not knowing in advance.

Both of this session's real defects, this one and the floor in section
7, were found by a check or a measurement rather than by reading the
code, and in both cases the code read correctly.

## 15. Final state

`R CMD check --as-cran --no-manual` with `_R_CHECK_CRAN_INCOMING_=false`
against core installed in the private library:

    Status: OK

Zero errors, zero warnings, **zero notes**. The one note an earlier run
carried, `Non-standard file/directory found at top level: '_pkgdown.yml'`,
was a missing `.Rbuildignore`: every sibling extension has one and this
package did not, because `ls -R` does not list dotfiles and the omission
was invisible until the check said so. ddm's file was copied unchanged.

A second undeclared dependency went the other way. `utils::head()` in one
print method meant `importFrom(utils, head)` in NAMESPACE with `utils`
absent from DESCRIPTION. R CMD check did not flag it here, and CRAN's
incoming checks might. The `head()` call was replaced with base
subsetting and the import dropped, which removes a dependency rather
than declaring one.

Both vignettes rebuild inside the check (12 s) and the suite runs inside
it (16 s, with the two `skip_on_cran()` tests skipped by design).

Nothing was committed. `git log --oneline -1` is still 5dfdd84.

---

# Punch round, 2026-09-05

Against `dev/review-spline.md` (789 lines, read end to end). Verdict
taken as given: curve half GO WITH FIXES, `royston_parmar()` NO-GO. The
review's punch table is the contract; every item below names the
file:line it landed at and what was measured.

Base is still 5dfdd84 and nothing was committed. Main has moved to
a57657f during the round; this branch was NOT rebased.

## 1. BLOCKER: the log S floor. Refusal by name, never a floor.

**New file `extensions/frmtmb.spline/R/rp-check.R`.** Exported
`rp_floored(object, action = c("error", "report"), max_nlogS = 19.2)`,
which refuses by default.

**Reproduced first.** The review's 600-subject design, one group-A
subject censored at t = 50, `df = 3` hazard scale:

    converged 0, no warning
    reported logLik  -468.82      AIC 947.63
    mu.grpB          -3.53230
    max H on the censored row     2424.30
    true log S there              -2424.30
    scored as                     -35.127363

Same class as the review's numbers (theirs: H = 21 690, mu.grpB -4.68);
the difference is the seed. The reported log likelihood is short by
about 2.4e+03 on one row.

**One threshold, derived and then measured, not three.** The error in
the scored `log S` is `eps / S`, because core forms `log(1 - F)` on the
probability scale and `1 - F` carries absolute error about
`.Machine$double.eps`. That depends on `S` ALONE, so it is the same on
every scale. Measured by forming `F` for a given `-log S` on each scale
and reading `log(1 - F)` back:

    -log S    computed        abs error    eps / S
    10        -10             1.3e-13      4.9e-12
    15        -15             9.0e-11      7.3e-10
    19.2      -19.2           2.4e-10      4.8e-08
    25        -25             3.8e-06      1.6e-05
    30        -29.99983       1.7e-04      2.4e-03
    36        -34.94504       1.05         9.6e-01
    40        -35.12736       4.87         5.2e+01

The three scales agree to **every printed digit at every row**. So the
review's two thresholds are one threshold said twice: `eta = 6` on the
normal scale is `-log S = 20.737`, and `-log S = 19.2` is `eta = 5.745`.
The single check is therefore slightly STRICTER on the normal scale than
the review asked for. `-log S` is `exp(eta)` on `"hazard"` (the
cumulative hazard `H`), `log1p(exp(eta))` on `"odds"` and
`-log(Phi(-eta))` on `"normal"`.

**Where it refuses from.** The review is right that this cannot go in
`family_finalize()`. It cannot go at fit end either: `logLik()` reads
`object$opt$objective` directly (`R/methods-fit.R:233-240`) and the
family protocol has exactly three post-fit slots, `post$mean_fn`,
`post$var_fn` and `post$dev_fn`, none of which core calls when a fit
finishes. **There is no hook.** So `logLik()` and `AIC()` cannot be
gated by any extension, and the refusal message says so in those words.

What IS gated, and it is every entry this package owns:

* `R/curve.R:131`, `R/curve-deriv.R:90`, `R/curve-feature.R:91` -
  `sp_rp_gate()` at the top of all three exported curve functions. The
  gate sits at the ENTRY rather than inside `sp_curve_parts()` because
  `frm_curve_feature()` returns early when its grid brackets no root and
  would otherwise never reach the covariance. The first draft had it in
  `sp_curve_parts()` and the test caught the hole.
* `post$mean_fn` already refuses, so `fitted()`,
  `predict(type = "response")` and `residuals()` were never at risk.

**The message** names the row count, the maximum, the threshold, the
reason and the remedy, as required. Rendered in the vignette:

    rp_floored(): this fit's reported likelihood is not the model's.
    1 of 1 censored rows are scored past the accurate region: -log S,
    which on the hazard scale is the cumulative hazard H, reaches 2424.3
    where this family stays accurate only to 19.2. frmtmb forms a
    right-censored term as log(1 - F) on the probability scale and core
    has no complementary log-CDF (lccdf) slot a family could use
    instead, so the scored log S is floored at -35.127363 and its
    gradient is exactly zero past -log S of 30. [...] The remedy for the
    censored term is the lccdf slot proposed in
    dev/spline-seam-proposal.md [...] Note that this check is POST-FIT:
    logLik() reads the optimizer's own value and no family hook runs
    when a fit finishes, so nothing here could have refused earlier.

**No check inside the objective.** There is no cheap tape-safe way. RTMB
refuses a comparison on an AD type and has no `CondExp`, and the
R-level density is called only at TAPE time, not during optimization, so
a check written there would run once at the starting values and never
again. Documented instead, in the message above and at
`R/royston-parmar.R:58-118`.

**Pinned** by `tests/testthat/test-rp-floored.R` (new, 6 tests, 33
expectations) on the review's own 600-subject design, including that the
fit converges without a warning, that all three curve functions refuse
it, and that `bc` passes with `max_nlogS < 5`.

## 2. `cens()` and `trunc()` downgraded

`R/zzz.R:48` and `R/zzz.R:50`. Both are now `conditional`, not `works`.

The `cens()` note carries the accuracy table, the measured 2.4e+03
failure, the fact that the identity test reaches only `-log S = 2.01`
and so cannot see it, and the condition for going back to `works`.
`trunc()` goes through the same term and additionally through
`R/objective.R:112` (`ll - log(Fub - Flb)`), so it is downgraded for
both reasons; its note records that an `lccdf` slot alone would NOT
close left truncation.

Three new rows register the checker itself: `rp_floored x cens()`,
`rp_floored x royston_parmar` and `rp_floored x frm_curve`, all `works`,
under a new `rp_floored = "method"` feature declared at `R/zzz.R:31`.

## 3. The monotonicity floor

`rp_floored()` counts rows whose fitted `d(eta)/d(log t)` is
non-positive on an observed row, and refuses on them, in the same call
and the same message template as item 1. I could NOT show the floor
never activates at a converged optimum, so it refuses.

Documented at `R/royston-parmar.R:58-85`, which used to say the floor
"is inert wherever the model is sane" and now says what it does where
the model is not: it makes `logLik()` and `AIC()` a pseudo-likelihood,
measured by the review at 6 rows and 3952 units.

## 4 and 5. The scaling wall and the fresh sdreport, one fix

`R/curve-cov.R:141-205`. Both items had the same answer and it is better
than either asked for: **read core's cached covariance instead of
building one.**

`predict(se.fit = TRUE)` runs `get_joint_cov()`, which memoizes the
inverted joint precision in `fit$cache$Vjoint`; `fit$cache` is a
`new.env(parent = emptyenv())` (`R/fit.R:938`), so the memo survives the
copy an extension makes. `sp_curve_parts()` makes that call anyway for
the check, so the object is free. The call was reordered to run first
(`R/curve-cov.R:255`).

This removes the second `sdreport()`, removes the Schur complement
entirely, and removes `as.matrix(Q)` with it, so item 4's line is gone
rather than fixed. The fallback that remains, for a fit where core never
formed a joint covariance, keeps `Q` sparse and uses `Matrix::solve()`
as item 4 asked. `Matrix` added to `Imports`.

Measured, `s(t, k = 8) + (1 + t | subject)`, 20-point grid,
`re.form = NA`, cold process each:

    random coefs   frm_curve before    after     peak RSS before   after
    506            (not run)           0.94 s    -                 288 MB
    2006           (not run)           2.74 s    -                 294 MB
    4006           3.18 s (review)     3.94 s    215 MB            476 MB
    8006           118.47 s (review)   8.24 s    2108 MB           1209 MB

**14x faster and 1.7x smaller at 8000 coefficients**, and the growth is
no longer in a dense matrix.

Item 5's correctness half is fixed rather than documented: the cached
object is `autoscale_sdreport()`'s, so an autoscaled fit now WORKS where
the old code would have been refused by its own check. Measured on
`y ~ s(x, k = 10) + z` with `z` at scale 1e6, `par_units` spanning
9.7e-07 to 1: `cov_rel_error` 6.66e-15, and agreement with
`predict(se.fit = TRUE)` to 3.68e-16 absolute. Registered as
`frm_curve x autoscale = works` in `R/zzz.R` and pinned in
`tests/testthat/test-curve.R`.

One more thing fell out. A fit with NO random-effect block used to be
refused, and that refusal was an artifact of this package recomputing a
joint precision that does not exist. Core's `get_joint_cov()` falls back
to `cov.fixed`, so such a fit now works and agrees with
`predict(se.fit = TRUE)` to **1.39e-17**. The test that pinned the
refusal now pins the agreement.

## 6 and 7. The Cost section

`R/curve.R:43-75` rewritten and `man/frm_curve.Rd` regenerated from it.
It now leads with the term that scales. Measured per stage,
`re.form = NA`, 20-point grid, one process each:

    model                                   coefs  design        predict(se.fit)  cov subset
    s(x, k = 10)                                8  0.01 s (11)   0.29 s           0.00 s
    s(t,k=8) + (1+t|subject), 1000 subjects  2006  0.07 s (101)  0.98 s           0.00 s
    the same, 4000 subjects                  8006  0.28 s (351)  6.87 s           0.00 s

The design rebuild is a tenth of the cost at every size and the
covariance subset is free. The section now says to size a job from the
number of coefficients in the FIT, not from the grid or the call count.

Item 7: the "110 coefficients, 32 calls" figure is MINE and it
reproduces, on the model the vignette fits,
`v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)` over 20 subjects:
110 random coefficients, 32 calls, `rel` 1.11e-15, and the population
curve is a real bell (`sd(.estimate) = 0.356`). The review's
non-reproduction is a DIFFERENT model, an `fs` term with no population
smooth, where at `re.form = NA` the term contributes nothing and the
curve is a constant. Both are now pinned in
`tests/testthat/test-curve.R`, which is what the review actually asked
for: no `bs = "fs"` model appeared anywhere in `tests/` before.

## 8. The identity test's reach

`tests/testthat/test-royston-parmar.R:52-63`. The identity sweep now
fits each configuration and asserts `rp_floored(action = "report")` on
it: `max_nlogS < 5` and `n_censored_floored == 0`. Measured on `bc`, the
largest `-log S` on any of the 387 censored rows is about 2.01, an order
of magnitude below the 19.2 where accuracy starts to go, so a reader
sees the test's reach beside its result.

## 9. CI

`.github/workflows/check-frmtmb-spline.yaml`. `NOT_CRAN: "true"` was
already set at the job level, inherited from ddm's workflow; the review
inferred the gap from a local `R CMD check` without it. It is now
repeated on the `check-r-package` step with a comment giving the reason
and the MEASURED counts, so it cannot be lost to a change in how the
action passes the environment down.

Measured rather than asserted: the suite runs **254** expectations with
`NOT_CRAN` and **207 with 4 skips** without it. Inside `R CMD check`,
where `test-message-uniqueness.R` skips itself because an installed
package has no sources to scan, that is **250 against 203**. The current
check log shows `FAIL 0 | WARN 0 | SKIP 1 | PASS 250`, and the one skip
is that self-guard, so both headline comparisons now run in CI.

## 10. The log S limit, documented

`R/royston-parmar.R:86-118`, the section now titled
`@section Censoring, truncation, and the accuracy limit on log S:`,
regenerated into `man/royston_parmar.Rd`. It carries the accuracy table,
the floor value -35.127363, the zero gradient past 30, the measured
converged-and-wrong fit, and pointers to `rp_floored()` and to the core
seam. `vignette("royston-parmar")` gains a section that runs
`rp_floored()` on a good fit and shows the refusal on the bad one.

## 11. The reviewer's edit

Kept. `R/rp-basis.R:49` still reads `sqrt(eps2) / 2`, which is correct:
both branches return 5e-08 at `u = 0`, and `sqrt(1e-14)` is 1e-07.

## 12 and 13. The seam proposal

`dev/spline-seam-proposal.md`, all four corrections folded in plus the
fifth the review gave in passing:

* `rowSums` is `R/predict.R:1235`; 1019 is the head of
  `predict.frmtmb_fit`. Both are cited now.
* A new section, "Core already caches both halves", giving `sdr_of()`
  (`R/fit.R:1209-1217`) and `get_joint_cov()` (`R/predict.R:12`), what
  reading them fixed, and the measured before and after. The cost
  argument for the seam is explicitly downgraded: it is worth exporting
  for reach and correctness, not for speed.
* `lccdf` is **five sites**, enumerated: the `frmtmb_family()` argument
  and its validation, the `frmtmb_ad_overload()` wrapping at
  `R/families.R:259`, an arity shim beside `fam_lcdf()`
  (`R/families.R:3898-3903`), the use site at `R/objective.R:100`, and
  the frame guard at `R/frame.R:1249`. And it does **not** close left
  truncation, which divides by the window at `R/objective.R:112`.
* The `ps()` closure moves from `nl_env` to `ev`, in the frame table and
  in a new paragraph quoting the evaluation at `R/objective.R:286-293`
  and the parse-time freeze at `R/parse.R:1357`. The parse gap is stated
  as the harder half: `nl_dpar()` learns the body from `all.vars(body)`
  at `R/parse.R:1247` and would never see a `ps()` call at all.
* The nl `se.fit` refusal is TWO sites, and the link-scale one at
  `R/predict.R:1179` is the one a nonlinear `frm_lp_basis()` must lift.

Item 13: Part 3 is retitled "A parser defect found on the way (a CORE
item)" and opens by saying it waits on nothing here and should be filed
and fixed on its own. All nine core citations in the document were
re-checked against the tree this round; every one resolves.

## Verification after the round

    file                        expectations   tests
    test-curve.R                          55      15
    test-deriv.R                          37      10
    test-gratia.R                         12       5
    test-message-uniqueness.R              4       1
    test-royston-parmar.R                 67       9
    test-rp-floored.R                     33       6
    test-surface.R                        46       9
    TOTAL                                254      55

    pass 254   fail 0   error 0   skip 0

Up from 187 over 46 tests. `test-message-uniqueness.R` passes unchanged
in shape (1 test, 4 expectations) and `test-surface.R` at 9 tests and 46
expectations.

`R CMD check --as-cran --no-manual`, `_R_CHECK_CRAN_INCOMING_=false`,
`NOT_CRAN=true`, core installed in the private library:

    Status: OK

Zero errors, zero warnings, zero notes. Both vignettes rebuild inside
it.

---

# Residuals, 2026-09-05

Against the re-check's "Residual items" (`dev/review-spline.md`, now
1132 lines). Verdict taken as given: GO for 0.1.0. R2, R3, R4 and R5
closed; R1 and R6 recorded rather than fixed, as instructed.

Base is still 5dfdd84, nothing committed, main left at a57657f and not
rebased onto.

## R2. The sparse fallback's ordering is now an assertion

`R/curve-cov.R:198-227`.

The reviewer's finding: the fallback bypasses `autoscale_sdreport()` and
so would hand an autoscaled fit an unscaled covariance, and it is
unreachable only because the `predict(se.fit = TRUE)` check at
`R/curve-cov.R:301` runs first and warms the cache. An ordering that
load-bearing should not rest on a comment.

**The ordering is now structural rather than asserted.** `sp_joint_cov()`
takes a third argument, `ref`, which IS the check call's return value:

    sp_joint_cov <- function(fit, des, ref) {
      if (missing(ref) || is.null(ref) || is.null(ref$se.fit)) {
        stop("frm_curve(): the joint covariance was asked for before the ",
             "predict(se.fit = TRUE) check that fills core's cache. That ",
             "ordering is not cosmetic: without the cache this falls back to ",
             "a fresh sdreport(), which goes round autoscale_sdreport() and ",
             "would return an unscaled covariance on an autoscaled fit. ",
             "Call the check first and pass its result", call. = FALSE)
      }

A caller cannot reach the covariance before the call it depends on,
because it needs that call's answer to get in. Moving the comment, or
moving the call, now fails loudly instead of silently taking the wrong
route.

**And the underlying correctness bug is closed too**, which the review
offered as the alternative fix. When the cache is empty AFTER the check
has run, the fallback refuses outright on an autoscaled fit
(`R/curve-cov.R:217-225`) rather than building an unscaled covariance:

    frm_curve(): this fit is autoscaled (par_units is not all 1) and core
    formed no cached joint covariance for it, so the only route left here
    is a fresh sdreport() on the UNSCALED objective, which would not be
    this fit's covariance. frmtmb exports no accessor for the autoscaled
    joint precision; see dev/spline-seam-proposal.md Part 1a

So the fallback is now correct on every fit it will answer for, rather
than correct only because nothing reaches it.

Pinned by `tests/testthat/test-curve.R:317-336`: `sp_joint_cov()` with
the argument missing, with `NULL`, and with an object that is not a
`predict(se.fit = TRUE)` result all refuse by name, and with the real
result it returns a matrix.

## R3. The internal read is named in the package's own documentation

Two places, as asked, neither of them `dev/`.

**`R/curve.R:43-83`**, a new `@section The one internal this reaches
into:` on `frm_curve()`, regenerated into `man/frm_curve.Rd`. It says
outright that `fit$cache$Vjoint` is a read into frmtmb's internals and
not a sanctioned seam, that `get_joint_cov()` is unexported and
undocumented, that neither the slot nor its `list(V =, names =)` shape
appears in `frmtmb-extension-api`, and that unlike `fit$obj` and
`fit$estimates` this one has no precedent to point at. It then gives the
trade (the alternative was slower AND wrong on an autoscaled fit) and
the three bounded failure modes the reviewer established:

* name or shape changes: reads `NULL`, sparse fallback runs, same answer
  about 16 times slower. No wrong number.
* meaning of `V` changes without the name changing: the covariance check
  refuses. No wrong number.
* absent: the ordinary case on every FRESH fit, because it is a memo
  rather than a slot. Not a dependency on a warm cache; the check call
  warms it, which is why it is ordered first and why R2's argument now
  enforces that.

**`R/frmtmb.spline-package.R`**, a new `@section What this package reads
that frmtmb does not promise:` on the package doc, listing all three
reads and marking which two have precedent and which one does not, with
a pointer to the `frm_curve()` section for the failure modes and to
Part 1a for the fix.

The user-visible claim is now the honest one: this reach cannot produce
a wrong answer, only a slow one or a refusal.

## R4. The proposal asks for both, and 1b's cost argument is gone

`dev/spline-seam-proposal.md`, Part 1 restructured into three asks:

* **Part 1a, an accessor for the cached joint covariance** (new). Export
  `get_joint_cov()` or a `frm_joint_cov(fit)` wrapper returning
  `list(V, names)`; five lines, since the function exists and already
  memoizes. The section says why this one is urgent: it is **the only
  unsanctioned dependency in the package**, and exporting it removes the
  reach outright. It also records that there is no exported route to the
  AUTOSCALED joint precision at all, which is why R2's fallback now
  refuses that combination rather than guessing.
* **Part 1b, `frm_lp_basis()`** (was the whole of the ask). **The cost
  argument is withdrawn in the document, in those words.** The measured
  figure is the reviewer's: the design rebuild this would replace is
  **0.28 s of a 5.91 s call** at 8006 random coefficients, so removing it
  entirely saves under five percent. What remains, and what nothing else
  buys, is REACH: `A` for a nonlinear body, where it is a Jacobian
  rather than a design, and for a reduced-rank block at
  `re.form = NULL`, where the loadings live in `theta` and unit
  perturbation cannot see them.
* **Part 1c, the `lccdf` slot** (renamed from "a second, one-line seam").
  Extended with the alternative core seam that would close the same
  class: a POST-FIT FAMILY HOOK. Core has none - `post$mean_fn`,
  `post$var_fn` and `post$dev_fn` are the only post-fit slots and core
  calls none of them at fit end - so `logLik()` and `AIC()` cannot be
  gated by any extension. Either seam closes R1, and the hook would
  close it for every family rather than for survival alone.

## R5. The seed-specific figure is gone

Four sites carried "a reported log likelihood wrong by 2.4e+03", which
is my seed; the reviewer's is 2.166e+04 on the same design. All four now
give the MECHANISM and both numbers as a range:

* `R/zzz.R:49` (the `cens()` compat note)
* `R/royston-parmar.R:114-127` (the accuracy-limit section)
* `R/rp-check.R:107-117` (the `rp_floored()` censored-floor section)
* `NEWS.md`

The wording they share: a floored censored row contributes -35.127363
instead of its own `-log S`, so the reported log likelihood is short by
about `-log S - 35` per floored row, which is thousands to tens of
thousands as soon as one censored time sits well past the event times;
two runs of the same 600-subject design differing only in seed give
2.4e+03 and 2.166e+04. Both converged without a warning and both put the
treatment coefficient out by tens of percent.

That is a property of the construction rather than of a seed, which is
what the review asked for.

## R1 and R6, recorded rather than fixed

As instructed, both are limits rather than work items.

**R1** was already documented in three places and is now in a fourth,
the package doc (`R/frmtmb.spline-package.R`, `@section Two limits that
are core's to fix:`), so it appears in `?frmtmb.spline` rather than only
on the family and the checker. It is also now in the proposal's Part 1c
with the two core seams that would close it, and in the "Two limits
recorded rather than proposed" section.

**R6** was NOT recorded anywhere before this round and now is, in the
same two places: the package doc and the proposal. Re-measured rather
than copied: `frmtmb_control()` takes no `map` argument, so a mapped `b`
block has no supported route; the nearest reachable analogue is a
distributional parameter held fixed, which sets `betad_fixed_idx` and
takes the same index-remapping path through `sp_coef_pos()`. On
`bf(y ~ s(x, k = 8), sigma = 0.5)` that is 8 `predict()` calls at a
`cov_rel_error` of **1.55e-15** (the reviewer measured 3.22e-15 on their
own build; both are machine precision, and the cited figure is now the
one this tree produces).

## Verification

    file                        expectations   tests
    test-curve.R                          59      16
    test-deriv.R                          37      10
    test-gratia.R                         12       5
    test-message-uniqueness.R              4       1
    test-royston-parmar.R                 67       9
    test-rp-floored.R                     33       6
    test-surface.R                        46       9
    TOTAL                                258      56

    pass 258   fail 0   error 0   skip 0

Up from 254 over 55. `test-message-uniqueness.R` passes at 1 test and 4
expectations, so the two new refusal templates added for R2 are unique
against the rest of `R/`. `test-surface.R` unchanged at 9 tests and 46
expectations.

roxygen regenerates with zero unresolved links; `man/frm_curve.Rd` and
`man/frmtmb.spline-package.Rd` carry the new sections.

`R CMD check --as-cran --no-manual`, `_R_CHECK_CRAN_INCOMING_=false`,
`NOT_CRAN=true`, core installed in the private library:

    Status: OK

Zero errors, zero warnings, zero notes. Inside the check the suite is
`FAIL 0 | WARN 0 | SKIP 1 | PASS 254`; the single skip is
`test-message-uniqueness.R` declining to scan an installed package that
has no sources, which is its own guard, and both headline comparisons
run. Both vignettes rebuild.

Scope, re-audited: `git -C . diff --name-only 5dfdd84` is still
`README.md` and `dev/build-docs.R`, zero root `R/` files are touched by
the diff or by the untracked set, and `git log --oneline -1` is still
5dfdd84.
