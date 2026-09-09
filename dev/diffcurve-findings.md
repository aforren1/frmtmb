# A difference curve in frmtmb.spline

Phase 1, item 1.3 of `dev/extension-gaps-plan.md`. Branch
`wt-diffcurve`, worktree `frmtmb-wt-diffcurve`, base 780dec1.

`frm_curve(object, newdata, contrast = newdata2)` returns `A1 - A2`
with covariance `(A1 - A2) V (A1 - A2)'`, pointwise and simultaneous.
`frm_curve_deriv()` and `frm_curve_feature()` carry the second grid, so
`frm_curve_feature(type = "crossing", at = 0)` on a difference locates
where the two curves meet.

## After the review

`dev/reviews/2026-09-08-diffcurve.md`. Verdict: not mergeable as it
arrived, mergeable after five fixes the reviewer made in this worktree.
Everything below is what changed, so that this file reads as the state
of the item rather than as its history.

**The reviewer's HIGH finding, and it is the one I should have found.**
`sp_spec()` resolved the contrast from `is.null(newdata)`:

    nd <- if (is.null(newdata)) s$newdata else newdata
    ct <- if (is.null(newdata)) s$contrast else contrast

so `frm_curve_deriv(dif, var = "x", newdata = g2)` dropped the second
grid and returned the FIRST curve's derivative, `identical()` to it,
with `what` missing ", differenced" and `spec$contrast` NULL. That is
the same silent wrong answer I reverted `sp_spec()` to demonstrate,
reached through a door I did not test, and it is what that function's
own roxygen says must not happen. Fixed by the reviewer: a new
`newdata` on a difference object now requires a new `contrast` beside
it, refused by name otherwise. Seen failing first at PASS=46 FAIL=2.

The lesson worth recording: I wrote the reversion test for the branch I
had just written and stopped. The block has two inputs, `newdata` and
`contrast`, and I tested one combination of them.

**Fixed by the reviewer, not redone here:** that refusal and its Rd
sections; `print.frmtmb_feature()`'s missing NA guard, which claimed a
covariance check on a nonlinear fit where none had run; the dead `D2`;
the `group_means = TRUE` sentence in the help page; and two of the
three absolute tolerances.

**Corrected by me, on the reviewer's measurements:** the gratia
standard-error story, in `NEWS.md`, the `s()` compat row and this file;
the `sp_one_basis()` comment, which claimed arithmetic independence
worth more than the one ulp it is; the API note, which claimed a silent
positional reinterpretation that does not happen; and the unreachable
"a grouping level the fit never saw" in the `extra_var` refusal, its Rd
and the `curve-cov.R` header.

**Closed by me:** the `gp()` refusal, narrowed to what is actually
unavailable, with the core seam filed; the third absolute tolerance.

`expect_gt(gap_se, 0.005)` guarded `2 * gap_se` from being a vacuous
bound. It is now a ratio to what the run measures about the same two
packages: the estimate gap relative to the curve's own scale, 6.5e-09,
against the standard-error gap, 0.031. The assertion is
`gap_se > 1000 * rel_est`, and it says something true and load-bearing
rather than merely non-zero: the two packages agree about the CURVE and
disagree about its standard error by orders of magnitude more, because
they are reading two covariances of one fit rather than two fits.


## The plan's claim, checked before it was built on

The plan says "everything it needs is already assembled in
`curve-cov.R`". Mostly true, and three things were missing. None is
large; two of them are the honest limits of the path and had to become
refusals rather than code.

1. **The two designs must sit at the same rows of `V`.**
   `frm_lp_basis()` returns `coef_pos`, so the question is answerable,
   but nothing asked it. Measured on the by-factor fixture:
   `identical(lbA$coef_pos, lbB$coef_pos)` is `TRUE`, 16 columns each.
   It is now compared on every difference call and a mismatch refuses
   by name rather than pairing columns that belong to different
   parameters.

2. **`extra_var` has no cross term, and that is a SEAM limit rather
   than a mathematical one.** The seam returns the variance that is not
   coefficient uncertainty as a per-ROW number and returns no
   covariance between two grids. The only source reachable through
   `frm_curve()` is an exact `gp()` off the observed positions: the new
   grouping-level source needs `allow_new_levels`, which this function
   does not offer and core refuses first (review finding 3, measured).

   It is not needed in the ordinary case. Where the two grids sit at
   the same `gp()` position they load the SAME kriging residual, it
   cancels exactly, and the answer is the matrix the code already
   forms. That case is now computed rather than refused; see "The
   `gp()` narrowing" below for the criterion and the identity. What
   remains refused is a difference whose grids sit at different `gp()`
   positions, and that is the core seam filed at the end of this file.

3. **The built-in check has nothing to check the difference against.**
   `predict(se.fit = TRUE)` returns a marginal standard error per row
   and never the covariance BETWEEN two grids, which is the whole
   content of a difference. So the check runs on each half,
   `cov_rel_error` is the worse of the two, `n_predict` is 2, and what
   the number licenses has changed from "this answer was verified" to
   "both designs were read correctly". `print()` says so in its own
   sentence rather than reusing the one-grid line. That is the reason
   `test-difference.R` carries three identities: the check cannot be
   one of them.

Everything else was already there and was reused unchanged:
`sp_catch_span()`, `sp_span_stop()`, `sp_grid_pinned()`,
`sp_sim_crit()` (which takes `S` and a divisor and does not care where
they came from), and `sp_assemble()`.

## The gratia identity

`y ~ fac + s(x, by = fac, k = 8)`, n = 400, seed 3, grid of 25 points
on [0.02, 0.98]. frmtmb ML against `mgcv::gam(method = "ML")`:
log-likelihood -167.6401 against `-gm$gcv.ubre` -167.6401.

`gratia::difference_smooths(select = "s(x)", group_means = TRUE)` is
the comparable call. The DEFAULT (`group_means = FALSE`) zeroes the
intercept and the parametric `fac` columns and reports the smooth-only
difference, which on this model is 0.522 away from the difference of
the two linear predictors; `group_means = TRUE` keeps them and is the
same quantity `contrast = ` reports.

| quantity | one curve | the difference |
|---|---|---|
| max abs estimate gap vs mgcv/gratia | 1.292e-08 | 1.954e-08 |
| max abs `se` ratio minus 1 | 0.0311 | 0.0156 |

The estimate gap is 1.51 times the one-curve gap, and the attribution
is provable rather than argued. The reviewer measured
`max |gapD - (gapA - gapB)| = 5.55e-16`: the difference's disagreement
with gratia IS the difference of the two curves' own disagreements, so
the difference arithmetic contributes nothing to it. It is the fit
difference.

**The standard-error sentence I first recorded was wrong, and the
correction matters.** I wrote that the difference's gap is half either
curve's because the common part of the disagreement cancels. The
reviewer measured all three, and 3.1 percent is curve B's alone:

| | conditional (`Vp`) | `unconditional = TRUE` |
|---|---|---|
| curve A | [0.9961, 1.0149] | [0.9885, 0.9955] |
| curve B | [0.9962, 1.0311] | [0.9842, 0.9952] |
| difference | [0.9961, 1.0156] | [0.9868, 0.9954] |

The difference's worst deviation, 1.56 percent, sits BETWEEN curve A's
1.49 and curve B's 3.11 rather than under both, and does the same under
`unconditional = TRUE` (1.32 between 1.15 and 1.58). Nothing the
measurement can see cancels; the difference tracks the better curve.
The reason for the gap at all still stands: mgcv's `Vp` is conditional
on the smoothing parameter and the joint precision this package inverts
is not, and neither of gratia's two covariances is this package's.
My error was taking the worse of the two curves and reading the
comparison as an improvement. `NEWS.md`, the `s()` compat row and the
sentence above are corrected.

`test-difference.R` asserts both against the run's own measurement
(`4 * gap_est`, `2 * gap_se`) rather than against a constant. `gap_se`
is the max over both curves, so `2 * gap_se` is a loose bound and is
now labelled as one.

## Two identities the check cannot give

**Core's own `vcov()`, bit for bit.** On `y ~ fac * x` with no random
effect anywhere, the difference of two predictions is a linear contrast
of `beta` and `vcov(fit)` is its covariance. `frm_curve(contrast = )`'s
`.se` is `identical()` to the one formed from `vcov()` and the same
contrast matrix. This is the one case where core reaches the number by
a second route, and it is exact rather than close.

**A grid differenced with itself.** `A1 - A2` is the zero matrix, not a
small one, so `.estimate` and `.se` are `identical()` to `rep(0, n)`.
Checked with `identical()` rather than a printed zero. The same call
with `simultaneous = TRUE` refuses, because a deviation process with no
uncertainty anywhere has no maximum to take a quantile of.

That second one is also the demonstration that a difference is not two
curves compared afterwards: adding the two curves' standard errors in
quadrature, which is what a reader doing it by hand would do, returns
`sqrt(2)` times a positive number for a quantity that is exactly zero.
On the by-factor fixture the two happen to agree to within 0.02 percent
(range of quadrature over the difference: 0.99983 to 1.00016), because
a `by=` smooth puts the two curves in almost disjoint coefficient
blocks. That near-agreement is a property of THAT model and not of the
method, which is why the self-difference is the case the test pins.

## The `gp()` narrowing, and the core seam it stops at

The first version of this lane refused every difference on a predictor
carrying variance that is not coefficient uncertainty. The reviewer
showed that on the lane's own fixture the refusal blocks an answer that
is exactly right, and left the narrowing to me. It is done.

### What is now computed

An exact `gp()` off the observed positions contributes a kriging
residual `g(x) - Xr g(obs)` per row. Its variance comes back from the
seam as `extra_var`; the covariance BETWEEN two grids does not. But
when the two grids sit at the same `gp()` position they load the same
residual, `var(g1 - g2)` is exactly zero, and the difference is
`(A1 - A2) V (A1 - A2)'` with nothing left over.

Measured on `y ~ fac + gp(x)`, n = 90, grid at the 89 midpoints between
observed positions, contrasting `fac`:

| quantity | value |
|---|---|
| difference standard error, every row | 0.0618114534972 |
| `sqrt(vcov(fit)["facB", "facB"])` | 0.0618114534973 |
| relative gap | 1.54e-12, about 6900 ulps |
| `sd(.se) / mean(.se)` | 0 exactly |
| spread of `.estimate` | 3.331e-16, about 6 ulps of 0.4516 |

The gp cancels completely and what is left is the `facB` contrast,
which core's `vcov()` reaches by a genuinely different route: it reads
`sdreport()`'s OUTER Hessian inverse, `cov.fixed`, a 3 x 3 matrix,
where this package inverts the 95 x 95 joint precision and takes a
block of it. The two are equal only through the Schur complement, and
they are NOT bitwise equal. That is why the agreement is evidence: the
same arithmetic gets `identical()` in this file's fixed-effect-only
test, and this gets 6934 ulps, which is what a 95 x 95 solve against a
3 x 3 inverse costs. The
standard error is flat BIT FOR BIT because the two designs' gp columns
are bit-identical, so `A1 - A2` is exactly zero there; the estimate is
flat only to a few ulps, because `eta1 - eta2` subtracts two numbers
that each carry the gp contribution.

### How sameness is established, and why not on the numbers

`sp_same_latent()` compares the DESIGN, not the variances: the two
grids must agree bit for bit on every column of `A` whose coefficient
is not a fixed effect. The component of each column comes from
`frm_joint_cov(fit)$names[coef_pos]`, which is a documented return
value, so there is no name parsing and no reach into core.

The coordinator's warning was right and the measurement bears it out
twice.

* **Two levels of one grouping block** have identical marginal
  variances by construction and are different draws. Measured on
  `y ~ fac + s(x, k = 6) + (1 | g)` at `re.form = NULL`, level 1
  against level 2: `identical(extra_var1, extra_var2)` is TRUE and the
  non-fixed design columns are NOT identical. A test on the numbers
  passes it; this one refuses. It is not reachable through
  `frm_curve()`, because an in-sample level carries no extra variance
  at all, so it is evidence about the predicate rather than a live
  hazard.
* **A mirrored `gp()` grid** is the reachable version of the same trap.
  On an exactly symmetric observed design,
  `seq(-5, 5, length.out = 41)`, the grids `x` and `-x` have kriging
  variances agreeing to 1.11e-16, a relative 1e-10 on a variance of
  1e-06. They are different draws: their residuals are correlated below
  one, so the difference does carry a variance the seam cannot supply.
  Any tolerant comparison of the two variances accepts it and
  understates the standard error. The weights are not mirrored, so the
  design test refuses it.

**Seen failing first, and the counts here are the second set.** The
first pair I recorded, 57/2 and 57/0, was stale in my own favour. I had
patched the CALL SITE, leaving `sp_same_latent()` in the namespace for
`test-difference.R:424` to call directly, so the one assertion that
reaches the predicate rather than the door could not move. The faithful
reversion replaces the predicate BODY, which is what an alternative
implementation would look like. Both arms below were rebuilt that way,
installed and rerun from the shipped files:

| predicate body | PASS | FAIL | failing assertions |
|---|---|---|---|
| `isTRUE(all.equal(extra_var1, extra_var2))` | 56 | 3 | 372, 399, 424 |
| `identical(extra_var1, extra_var2)` | 58 | 1 | 424 |
| the shipped design test | 59 | 0 | none |

All three total 59, so it is the same file. The call-site patch gives
57/2 and 57/0 for the reason above; the body patch is the number to
record.

That makes the conclusion STRONGER than I first wrote it, not weaker. I
said the design test was "only luckily different from the bitwise one",
because I thought the bitwise gate passed the suite. It does not: it
fails at line 424, on an assertion I wrote deliberately. The honest
statement is that the design test and the bitwise numeric test agree on
everything reachable through `frm_curve()` today, and the suite
separates them anyway, because the grouping-block case reaches the
predicate directly. The design test is strictly better than both
numeric tests and this file proves it.

Three lanes this round reported a count that was stale, capped or an
artifact rather than a measurement, and mine was one of them. The two
numbers above were re-derived from the shipped files rather than
carried forward.

**And the instrument had the same defect.** My per-file runner summed
`failed` and not `error`. A `test_that()` block that THROWS aborts the
rest of its own assertions, so an unmatched `expect_error()` message
showed up as PASS=51 FAIL=0 on a file that asserts 59: eight assertions
silently did not run and the line read clean. It happened once in this
lane, on the refusal text I rewrote for R6, and it is the same failure
mode as the stale counts. The runner now prints `ERROR` and `BAD` and
the suite below is reported from it. Every count in this file was taken
after that change except the two reversion arms, whose failures were
all Failures rather than Errors and which the reviewer measured
independently to the same numbers.

### Sensitivity, and what it assumes

On the `y ~ fac + gp(x)` fixture the predicate returns TRUE for a
contrast across `fac` at one `x`, and FALSE for a second grid whose `x`
moved by 0.05, by 1e-06, by 1e-10, or was mirrored about the range
midpoint.

The assumption is that the kriging weights fix the position uniquely.
They do unless every observed position is equidistant from the two
grids' positions: in one dimension that needs every observation at one
point, which makes the kernel matrix singular; in two it needs them
collinear with the two positions mirrored across that line. On such a
design the extra term would be understated rather than refused. That is
recorded in `sp_same_latent()`'s own block. Everything else fails
closed.

### CORE SEAM: the conditional cross-covariance of an exact `gp()`

Filed for core, precisely enough to act on. Not folded into
`dev/extension-gaps-plan.md` from this lane, because a sibling lane
owns that file this round.

**What core computes today.** `R/predict.R:383-405`, inside the `gps`
loop of `lp_eta_design()`, for a grid with any unseen position:

    th  <- fit$estimates[["theta"]][bk[["theta_idx"]]]
    K   <- unname(covstruct_registry[["gp"]]$vcov(th, bk))  # obs x obs
    Ks  <- gp_cross_cov(th, bk, Xc, pos)                    # n x obs
    Xr  <- t(solve(K, t(Ks)))                               # n x obs
    kss <- exp(2 * th[1]) * (1 + 1e-6)
    extra_var <- pmax(kss - rowSums(Xr * Ks), 0)             # length n

So `K`, `Ks` and `Xr` all exist, and `gp_cross_cov()` already computes
a kernel between two arbitrary position sets. What is thrown away is
the off-diagonal: the conditional covariance between prediction rows
`i` and `j` is `k(x_i, x_j) - Xr_i K Xr_j'`, and every factor of it is
in scope at that line. `lp_extra_var_vec()` then reduces the whole
thing to a length-`n` vector before `frm_lp_basis()` returns it
(`man/frm_lp_basis.Rd`, `extra_var`).

**What core would have to expose.** Either of:

* `frm_lp_basis(newdata, extra_cov = TRUE)` returning the `n x n`
  conditional covariance in place of, or beside, `extra_var`. Its
  diagonal is `extra_var`, so nothing existing changes meaning.
* or `frm_extra_cov(object, newdata, newdata2, ...)`, an accessor for
  the cross-covariance between two grids.

**The first is sufficient for a difference, and that is the point the
first draft of this filing missed.** A difference gets its cross term
from ONE call on the stacked grid `rbind(newdata, contrast)`, reading
the off-diagonal `n x n` block of the `2n x 2n` result. So the `n x n`
shape is not "one grid only"; the second accessor is not needed for
this item at all. It is also the smaller change and it closes a second
defect at the same time: a SINGLE grid's simultaneous band over an
exact `gp()` off the observed positions is too narrow today, by at
least 17 percent where the grid extrapolates. That one is shipped in
0.4.0, is not a difference-curve defect, and now has its own `NEWS.md`
bullet and an entry in `dev/feature-gaps.md`.

**The matrix must carry the new-level source too.** This filing named
only `gp()`, and `extra_var` also carries a new grouping level's
marginal variance. `extra_var_blocks()` (`R/predict.R:1660`) already
decides which rows share a draw, through its `new_key` grouping, so the
information exists; a seam that returned a matrix for `gp()` and a
diagonal for new levels would be an object whose meaning depends on the
model, which is worse than either.

**Which extension items wait on it.**

* **1.3, this item.** A difference curve whose two grids sit at
  different exact `gp()` positions. Refused today, by name, with the
  reason given.
* **4.1, RP predictions.** The plan routes `rp_predict()`'s bands
  "through the same covariance path `frm_curve()` uses". A Royston and
  Parmar fit with a `gp()` term on any spline coefficient meets the
  same limit for any contrast between two covariate profiles.
* **The predicate's own over-strictness.** `sp_same_latent()` asks
  every non-fixed column to match, where only the block carrying
  `extra_var` has to, so `y ~ fac + s(x, by = fac) + gp(x)` is refused
  with its `gp()` columns bit-identical. A matrix-valued `extra_cov`
  removes the predicate and the restriction with it. See "Found, not
  fixed" below.
* Anything later that wants a simultaneous band over an exact `gp()`
  evaluated off the observed positions, for the reason in the first
  bullet above.

The size looks small: the pieces are already local to one loop, and the
work is returning a matrix instead of its diagonal plus one more
`gp_cross_cov()` call. No extension can do it. The plan forbids `:::`,
and `K`, `gp_cross_cov()` and `covstruct_registry` are all internal.


## Simultaneous-band coverage, as a rate over seeds

200 seeds. Per seed: n = 300, two groups, `y ~ fac + s(x, by = fac,
k = 6)`, truth `2 sin(pi x)` against `1.4 sin(pi x + 0.6)`, noise 0.35,
grid of 30 points on [0.05, 0.95], `nsim = 4000`. Whole-curve coverage
is "the true difference is inside the band at every one of the 30
points".

| quantity | measured |
|---|---|
| simultaneous, whole curve | 0.970 |
| pointwise, whole curve | 0.685 |
| pointwise, per point | 0.956 |
| mean simultaneous critical value | 2.818 |
| binomial mcse of a rate at 0.95, 200 seeds | 0.0154 |

27.4 s for the 200 seeds in a fresh process through the public
`frm_curve()`, which is two `predict(se.fit = TRUE)` checks per
replicate.

The band is conservative rather than short, which is the direction a
simultaneous band is allowed to err. The three assertions in
`test-difference.R` are all multiples of a spread the run measures, the
way `test-gratia.R` bounds its critical-value comparison by
`3 * mine$mcse`:

* `mean(hit_sim) > 0.95 - 3 * mcse`, that is 0.970 against 0.904.
  Under-coverage is the failure mode.
* `mean(hit_pt) < mean(hit_sim) - 6 * mcse`, that is a gap of 0.285
  against 0.092. This is the control, and it is the whole reason the
  simultaneous band exists.
* per-point pointwise coverage within `5 * se` of 0.95, where `se` is
  the spread ACROSS replicates divided by sqrt(200), not a binomial
  standard error: 30 points on one curve are not 30 draws. Measured
  |0.956 - 0.95| = 0.006 against 0.0274.

No absolute numeric threshold is written anywhere in the file.

## What the 0.55.0 span gate does to a difference curve

It does not fire, and that had to be constructed rather than reasoned
about.

The gate asks whether a column other than the search variable varies
DOWN a grid. A difference is built from two grids that differ from ONE
ANOTHER by construction, and those are different questions. Each grid
of a normal difference curve is pinned in every column but `var`
(`fac` is `"A"` all the way down `newdata` and `"B"` all the way down
`contrast`), so `sp_grid_pinned()` is `TRUE` for both and no extra
prediction is made. Measured directly on the fixture below: both
`TRUE`.

The hazard is an implementation that STACKS the two grids and hands the
stack to the gate. There `fac` varies down the stack, the gate fires,
and on a model with a `ps()` term the extra prediction can refuse a
difference curve that is entirely inside every span. This
implementation never stacks: `sp_curve_parts()` reads the seam twice
and subtracts.

**Fixture, and why it is not the existing two-`ps()` one.** The
existing `sp_span_fit2()` is additive, so a difference across any of
its columns is constant in `t`, the scan finds no sign change, and the
second gate is never reached (it sits after the zero-root early
return). `test-span.R` gains `sp_span_fit3()`, the same two `ps()`
terms with an interaction:

```
y ~ lev + ps(t, k = 10, pad = 0.3) + b1 * w + b2 * w * t
      + ps(z, k = 8, pad = 0.02),  nl = TRUE
```

The difference across `w` is `b1 + b2 t`, which crosses zero. Fitted at
n = 400: `b1` = -0.539 (0.048), `b2` = 1.436 (0.084), crossing at
0.3753 (se 0.0188) against the truth 0.4. `t` span [-0.2970, 1.2981],
`z` span [-0.0197, 1.0192].

Three cases, all in `test-span.R`:

| grid | 0.55.0 gate, extended | result |
|---|---|---|
| `nd` z = 0.5 w = 1, `ct` z = 0.5 w = 0 | both pinned | ACCEPTED, root found |
| `nd$z[10]` out of span | `nd` not pinned | refused, quotes `ps(z, k = 8` |
| `ct$z[10]` out of span | `ct` not pinned | refused, tags the contrast |

The gate is not weakened for the case it was built for: it is the same
predicate, asked once per grid, and the one-grid path is byte for byte
what 0.55.0 shipped. Every span message now carries which grid raised
it, because a reader told "a value is outside the frozen knot span" of
a two-grid call otherwise has to guess.

### Seen failing first

Two reversions, each installed and run against the new tests.

**The gate on `newdata` only** (drop the `|| (!is.null(ct) &&
!sp_grid_pinned(ct, var))` clause). `test-span.R` PASS=68 FAIL=1:

```
1. Failure ('test-span.R:441:3'): an out-of-span row refuses from
   either grid of a difference
Expected `frm_curve_feature(...)` to throw a error.
```

That is the 0.55.0 defect reproduced from the other side: the same
out-of-span row is refused or ignored depending on which of the two
grids holds it, and the search never predicts at it either way.

**`sp_spec()` dropping the contrast** (return `contrast = NULL` from
the `frmtmb_curve` branch). `test-difference.R` PASS=32 FAIL=5,
including:

```
4. Failure: Expected max(abs(d1$.estimate - a1$.estimate)) >
   0.1 * max(abs(a1$.estimate)).  Actual comparison: 0.00 <= 0.69
5. Failure: nrow(frm_curve_feature(dif, ...)) to equal 1L.
   actual: 0
```

`0.00` is the point: the derivative of the difference came back
bit-identical to the derivative of the FIRST curve, with no complaint,
on an object whose every printed row is a difference. That is the
silent wrong answer the plan's rule 3 puts first, and it is why the
contrast is threaded through `sp_spec()` rather than accepted only as
an argument.

With both reversions undone: `test-span.R` PASS=73 FAIL=0,
`test-difference.R` PASS=42 FAIL=0.

## The recorded e2 defect, and whether it is this item's

`dev/debts-findings.md`, item 3, "Found, not fixed": for a stationary
point the plus and minus `e2` stencil points enter `.se` through `f2`,
so a maximum located within `e2` of a knot has its curvature computed
across the partition-of-unity cliff. For a crossing they were said to
enter nothing reported. A difference curve searches for a crossing, so
the question is whether that is now this lane's problem.

It is not, and the reason is measured rather than read. `eps` moves
`e1` and `e2` together through the public argument, so the reported
block was replicated with `e1` held at 1e-6 of the grid range and `e2`
moved by a factor of 100 (1e-4 to 1e-2 of the range), on the by-factor
fixture at the located crossing:

| reported | identical | value |
|---|---|---|
| `.se` | TRUE | 0.023344393861040272 |
| `.value` | TRUE | 4.4408920985006262e-16 |
| `.value_se` | TRUE | 0.075175120356382269 |
| `slope` attribute (`f1`) | TRUE | 3.2202643942639648 |
| `curvature` (`f2`) | FALSE | 2.1959008926409398 then 2.1936734673 |

Every number a crossing reports is bit-identical under a hundredfold
change in `e2`; only the `curvature` attribute moves. The recorded
defect is confined to `type` in `c("maximum", "minimum", "extremum")`,
which this item does not use, and nothing here changes it. The script
is `diffcurve-probe7.R` in the lane scratchpad.

## The check works only because its arithmetic is redundant

Worth recording, because it was nearly optimized away in this lane.

A difference discards both halves' full covariances and keeps only the
difference's, so `sp_one_basis()` looks like it should form the
diagonal directly, `rowSums((A %*% V) * A) + extra_var`, rather than
`diag(A V A')`. That is also exactly how core writes it: the
`frm_lp_basis()` help says "`var(eta)` is
`rowSums((A %*% V) * A) + extra_var`, and `predict(se.fit = TRUE)` is
written that way".

Which is the reason not to. Written core's way the check compares two
runs of one formula, `cov_rel_error` comes back 0 exactly on this
package's own fixtures, and `test-curve.R`'s "a tolerance no covariance
could meet refuses rather than returning" stops refusing at `tol = 0`.
It was built, installed and run: `test-curve.R` PASS=61 FAIL=1, the one
failure being that `frm_curve(tol = 0)` no longer throws.

**I first recorded the reason as "an independent arithmetic route", and
that overstates it.** The reviewer measured what the independence is
worth. Both spellings read the SAME `A`, `V` and `extra_var` out of the
same `frm_lp_basis()` call, so the only thing separating them is
summation order, one ulp: with `V`'s first two rows and columns
swapped, BOTH spellings refuse at 0.101 relative. A real misreading of
the seam is caught either way. What the triple product actually buys is
that `tol` can be set to zero and still mean something, and one ulp of
arithmetic redundancy on top. That is a smaller claim and it is the
true one; the comment on `sp_one_basis()` now makes it.

The decision does not change. `sp_one_basis()` keeps `diag(A V A')`,
and a difference pays for two half covariances it discards. The
documented cost of this call is dominated by the single
`predict(se.fit = TRUE)` joint-precision solve (`@section Cost:` of
`frm_curve()`), so the trade spends work that is not the bottleneck.
The measurement is in the comment, so the next reader who spots the
waste finds it before deleting the line.


## Found, not fixed

* **`D2` in `curve-feature.R` was dead; the reviewer removed it.**
  `D2 <- (blk(5L) - 2 * blk(3L) + blk(1L)) / e2^2` was assigned and
  never read, an `nr x p` product formed and discarded on every feature
  call. I recorded it and left it as out of scope. The reviewer
  confirmed it was dead, confirmed it was pre-existing at 780dec1
  rather than this lane's, and removed it, with no behavior change and
  no test moved. Recorded here because leaving it was my call and it
  was the wrong one: it is one line, it is provably unread, and "no
  test fails before the fix" is not a reason to keep waste in a file I
  had already rewritten.

* **The degenerate-band refusal names the wrong usual cause for a
  self-difference.** `frm_curve(newdata = g, contrast = g,
  simultaneous = TRUE)` refuses with "every point on this grid has a
  standard error of exactly zero", which is true, and then offers "past
  a `ps()` term's knot span" as the usual way to arrive there, which
  for a self-difference it is not. The first sentence is correct and
  the second is offered as a hint, so the message is not wrong, only
  incomplete. Changing it would touch a template `test-span.R` matches
  on; left for whoever next edits that refusal.

* **`gratia::difference_smooths()`'s default is a different quantity,
  and the help page now says so.** `group_means = FALSE`, the default,
  reports the SMOOTH-only difference and drops the parametric group
  offset; on the fixture here that is 0.522 away from the difference of
  the two linear predictors. I recorded it and asked whether the help
  page should say it. The reviewer answered that it must, because the
  sentence "on the same mgcv fit the two agree" is false at gratia's
  defaults, and wrote it. Agreed, and the same holds for
  `dev/diffcurve-findings.md`'s own gratia section, which names
  `group_means = TRUE` throughout.

* **A single grid's simultaneous band over an exact `gp()` is too
  narrow, and I under-reported it.** `sp_sim_crit()` draws from
  `Sigma = A V A'` and standardizes by
  `se = sqrt(diag(A V A') + extra_var)`, which off the observed
  positions are not the same object. I wrote that it could not be
  measured without the missing seam. It can, and the reviewer did:
  `extra_var` comes from `frm_lp_basis()` and the divisor from
  `predict(se.fit = TRUE)`, both public, so rerunning the same
  simulation standardized by `sqrt(diag(Sigma))` bounds the correction
  directly. On `y ~ fac + gp(x)`, 60 observations on [0, 6], noise 0.2,
  `nsim = 20000`: 1.00009 of the correct critical value inside the
  data, 1.00007 at the edge, **1.17485 extrapolating**, and that is a
  LOWER bound because rescaling keeps `A V A'`'s correlation. So a
  simultaneous band is at least 17 percent too narrow exactly where it
  should widen. Pre-existing, confirmed on a 780dec1 build with
  identical numbers, and a ONE-GRID defect with nothing to do with a
  difference. Not fixed: the correct fix is the core seam above. It now
  has its own `NEWS.md` bullet under the development heading and an
  entry in `dev/feature-gaps.md`, because a reader of a seam filing
  would not find it.

* **`sp_same_latent()` is stricter than the mathematics needs.** It
  asks every non-fixed column of `A` to match, where only the block
  carrying `extra_var` has to for the residual to cancel. Measured by
  the reviewer on `y ~ fac + s(x, by = fac, k = 6) + gp(x)`, n = 120,
  contrasting `fac`: `identical(extra_var1, extra_var2)` TRUE, the
  `gp()` columns of `A` bit-identical, all non-fixed columns NOT
  identical, and the call refuses. That is one step from the model this
  whole item is built on. It fails closed, so it is a missing feature
  rather than a wrong answer, and the refusal message, the help page,
  the `gp` compat row and `NEWS.md` now all say what the test actually
  keys on rather than blaming the `gp()` positions. Narrowing it needs
  the block identity, whose only public route is the
  `b.<block>.<level>` shape of `frm_joint_cov()$labels`, which is the
  name parsing this package deliberately avoids; the core seam removes
  the predicate and the restriction together.

## Decided against

* **Recycling a one-row `contrast`.** `A1 - A2` is row by row, and a
  one-row contrast silently recycled against a fifty-row grid is a
  different quantity from the one the argument name promises. Equal
  height is the contract and the refusal tells the user to repeat the
  row themselves, which says so in their own code.

* **A "fixed reference profile" mode in the derivative and the
  feature.** Both move `var` in each frame, so a `contrast` whose `var`
  column disagrees with `newdata`'s has had it overwritten. The
  alternative was a mode switch ("a contrast pinned in `var` is a fixed
  reference; one that tracks the grid moves with it"), which is two
  rules where one will do and a silent change of meaning between them.
  `frm_curve()` itself accepts such a contrast, because a curve minus a
  reference profile is well defined; only the two functions that MOVE
  `var` refuse it.

* **Differencing `extra_var`.** See point 2 above. Where both grids
  hold the same `gp()` covariate values the contribution does cancel,
  but nothing in the seam's output distinguishes that from the case
  where it does not, so guessing would be a silent wrong answer in the
  second case.

## API note for consolidation

`contrast` is inserted as the THIRD positional argument of
`frm_curve()` (before `dpar`), the fifth of `frm_curve_deriv()` and the
sixth of `frm_curve_feature()`, which is where the plan's spelling
`frm_curve(object, newdata, contrast = newdata2)` puts it. Every call
site in the repository passes `dpar`, `resp` and the rest by name, so
nothing in tree breaks; `R CMD check` rebuilt both vignettes clean.

**I wrote that an out-of-tree positional caller would silently
reinterpret its third argument. It does not.** The reviewer measured
all three signatures: the old third positional was `dpar`, a character,
and `sp_check_contrast()` refuses a character by name, so
`frm_curve(fit, gA, "sigma")` and both siblings raise "`contrast` must
be a data frame with at least one row". A silent reinterpretation would
need a data frame passed positionally where `dpar` was expected, which
is not a call anyone writes. The real side effect is smaller and worth
knowing: partial matching now resolves `c = `, `co = ` and `con = ` to
`contrast`, where they used to be unused-argument errors.

The reviewer also settled the placement, and I agree: keep it third.
None of the three signatures has a `...`, so "after `...`, so it must
be named" is not available without adding one, and a `...` added to
force naming would start swallowing misspelled arguments in exchange
for a risk the measurement shows is not there.

## Verification

Installed into the lane's private library. Test files run one per R
process, `NOT_CRAN=true`:

| file | pass | fail |
|---|---|---|
| test-bracket-access.R | 1 | 0 |
| test-curve.R | 62 | 0 |
| test-deriv.R | 38 | 0 |
| test-difference.R | 59 | 0 |
| test-gratia.R | 13 | 0 |
| test-message-uniqueness.R | 4 | 0 |
| test-royston-parmar.R | 67 | 0 |
| test-rp-floored.R | 42 | 0 |
| test-span.R | 73 | 0 |
| test-surface.R | 48 | 0 |
| total | 407 | 0 |

390 before the review, 396 after the reviewer's fixes, 407 now: the
extra assertions are the two `gp()` blocks that replace the single
refusal test, and the reworked `gap_se` guard. Counted with a runner
that sums errors as well as failures, for the reason above.

`test-difference.R` is the slow one, 29 s to 57 s of wall clock
depending on what else was on the machine, almost all of it the
200-seed coverage loop. It is `skip_on_cran()`. The reviewer notes it
is about 60 percent of the package's test wall clock and should move
behind a gate rather than get shorter if it grows; agreed.

`R CMD check --as-cran` on the built tarball, with pandoc and TinyTeX
on PATH: **1 WARNING, 1 NOTE**.

* The NOTE is the expected one: "checking HTML version of manual ...
  Skipping checking math rendering: package 'V8' unavailable".
* The WARNING is pre-existing and untouched by this lane: CRAN incoming
  feasibility reports "New submission", "Strong dependencies not in the
  CRAN or BioC software repositories: frmtmb", and a 301 on the
  `DESCRIPTION` URL `https://aforren1.github.io/frmtmb/frmtmb.spline`
  (missing trailing slash). `DESCRIPTION` was not modified.
* `checking tests ... [86s] OK`, `checking examples ... OK`,
  `checking Rd files ... OK`, `checking Rd line widths ... OK`,
  `checking for code/documentation mismatches ... OK`, `checking
  re-building of vignette outputs ... [27s] OK`.
* One check in this lane came back with an extra NOTE, "Examples with
  CPU or elapsed time over 5s: frm_curve", and the reviewer's did too.
  Neither is a regression: no example in the diff changed, and both
  runs had a probe or an install on the same machine. Alone,
  `checking examples ... OK`. Recorded because it is the same
  instrument-before-measurement trap as the runner defect above, and
  because a check run under contention should not be reported.
* The reviewer built the same check at 780dec1 and on the lane as
  received, and the three `00check.log` bodies are identical line for
  line apart from timings and paths, so the WARNING being pre-existing
  is measured rather than asserted.

Nothing was committed.

## Files

* `extensions/frmtmb.spline/R/curve-cov.R`: `sp_one_basis()`,
  `sp_cov_check()`, `sp_span_both()`, `sp_grid_span()`,
  `sp_same_latent()`, and `sp_curve_parts(contrast = )`. The unused
  `se_ref` element was dropped from the parts list; nothing in `R/` or
  `tests/` read it.
* `extensions/frmtmb.spline/R/curve.R`: `frm_curve(contrast = )`, the
  `transform` refusal, `sp_check_contrast()`, the `spec$contrast`
  slot, and `print()`'s difference branch.
* `extensions/frmtmb.spline/R/curve-deriv.R`: `sp_spec()` carries the
  contrast, `sp_check_contrast_var()`, the stacked contrast stencil.
* `extensions/frmtmb.spline/R/curve-feature.R`: the two-frame
  `eta_at()`, the two-frame stencil, the second span gate asked of each
  grid, and `print()`'s difference branch.
* `extensions/frmtmb.spline/R/zzz.R`: `frm_curve_contrast` registered
  as a method, with ten compat rows.
* `extensions/frmtmb.spline/NEWS.md`: a
  `# frmtmb.spline (development version)` section.
* `extensions/frmtmb.spline/man/frm_curve.Rd`,
  `frm_curve_deriv.Rd`, `frm_curve_feature.Rd`: roxygenised.
* `extensions/frmtmb.spline/tests/testthat/test-difference.R`: new.
* `extensions/frmtmb.spline/tests/testthat/test-span.R`:
  `sp_span_fit3()` and the two difference-curve span tests.
