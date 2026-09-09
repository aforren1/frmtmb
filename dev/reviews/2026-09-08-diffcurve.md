# Review: the difference curve in frmtmb.spline

Item 1.3 of `dev/extension-gaps-plan.md`. Lane worktree
`frmtmb-wt-diffcurve`, branch `wt-diffcurve`, base 780dec1. Reviewer did
not write the code. Nothing was committed in either tree.

## Verdict

**Not mergeable as it arrived. Mergeable now**, with findings 1, 6,
7, 8 and 9 fixed in the lane worktree and every affected test rerun.

One HIGH finding: a second door into the exact silent wrong answer the
lane built `test-difference.R` to close. Passing a new `newdata` to
`frm_curve_deriv()` or `frm_curve_feature()` on a difference curve
dropped the contrast and returned the FIRST curve's answer, without a
word, on an object whose every row is a difference. Measured on the
lane build as received: `max|d1 - a1| = 0.000e+00`, the same "0.00" the
lane records as the mark of that defect. Fixed, with the test seen
failing first.

The rest is documentation accuracy, and two of those claims are wrong
in a way that will mislead a reader: the standard-error comparison
against gratia (F4) and the justification for the `gp()` refusal (F2).

The arithmetic is right. Every identity the lane records reproduces,
including both `identical()` claims, and the gratia attribution that the
task called load-bearing is not just plausible but provable: the
difference's disagreement with gratia IS the difference of the two
curves' own disagreements, to 5.55e-16.

## What I ran

Private libraries only: `rvdc-lib` (lane), `rvdc-main-lib` (780dec1),
and six single-purpose arms (`rvdc-lib-corespell`, `-permv`,
`-corespell-permv`, `-nogate`, `-nospec`, `-e2`). One test file per R
process, `NOT_CRAN=true`.

`frmtmb.spline`, after my fixes:

| file | pass | fail |
|---|---|---|
| test-bracket-access.R | 1 | 0 |
| test-curve.R | 62 | 0 |
| test-deriv.R | 38 | 0 |
| test-difference.R | 48 | 0 |
| test-gratia.R | 13 | 0 |
| test-message-uniqueness.R | 4 | 0 |
| test-royston-parmar.R | 67 | 0 |
| test-rp-floored.R | 42 | 0 |
| test-span.R | 73 | 0 |
| test-surface.R | 48 | 0 |
| total | 396 | 0 |

The lane's own tally, 390 / 0, reproduced exactly before my changes.

Core, four files that the seam this item reads could plausibly reach:
`test-predict-lp-basis.R` 30/0, `test-compat.R` 573/0,
`test-compat-register.R` 103/0, `test-predict-newdata.R` 12/0. Core is
untouched by the lane.

## Findings

### 1. HIGH (fixed). A new grid silently returned the first curve

`sp_spec()` resolved the contrast from whether `newdata` was absent:

```r
nd <- if (is.null(newdata)) s$newdata else newdata
ct <- if (is.null(newdata)) s$contrast else contrast
```

So a caller who reuses a difference curve WITH a new grid, and does not
happen to repeat the contrast, gets the second grid dropped. Measured on
the lane build as received, on the `y ~ fac + s(x, by = fac, k = 8)`
fixture:

* `frm_curve_deriv(dif, var = "x", order = 1, newdata = g2)` returns an
  estimate `identical()` to `frm_curve_deriv(fit, newdata = g2)`;
  `max|d_new - a_new| = 0.000e+00`.
* `attr(., "what")` is `"derivative of order 1 in x"`, with no
  `", differenced"`, and `attr(., "spec")$contrast` is `NULL`.
* `frm_curve_feature(dif, var = "x", type = "crossing", at = 0,
  newdata = g2)` reports 0 roots where the difference has 1.

That is the plan's rule 3 case, and it is what `sp_spec()`'s own
roxygen block says must not happen: "Dropping it would be the worst
kind of answer: the derivative or the feature of the FIRST curve,
returned without complaint."

**Fix.** `sp_spec()` refuses when the object carries a contrast,
`newdata` is given and `contrast` is not. Both grids together still
work and still difference. An ordinary curve is untouched.

**Seen failing first.** With the new test block against the unfixed
build: `test-difference.R` PASS=46 FAIL=2, both failures being the two
`expect_error()` calls. After the fix: PASS=48 FAIL=0.

Rd sections for `frm_curve_deriv()` and `frm_curve_feature()` say so,
roxygenised with roxygen2 8.1.0, which is the version `DESCRIPTION`
pins. Only the two files I touched were rewritten.

### 2. MEDIUM. The `extra_var` refusal is over-broad, for a wrong reason

The lane reports this as a structural limit. It is a seam limit, and
narrower than recorded.

On the exact fixture `test-difference.R` uses, `y ~ fac + gp(x)` at
n = 90 with the grid offset off the observed positions:

* the two grids' designs differ in ONE column, `facB`; all 90
  `gp(x)` columns of `A` are bit-identical between them;
* `identical(lbA$extra_var, lbB$extra_var)` is `TRUE`, range
  [6.819147e-07, 8.911105e-07].

The kriging weights `Xr = K* K^-1` ARE those columns of `A`. Equal
weights mean the same conditional residual random variable, so
`var(g1 - g2) = 0` and the correct difference variance is exactly
`(A1 - A2) V (A1 - A2)'`, which the code has already formed and then
throws away. Computed: a flat 0.0618115 across the grid. The refusal
blocks a case whose answer is exactly right, and the sentence "the
difference has no variance to report" is false for it.

The general case is derivable too, just not from the public seam. Core
already has every piece at `R/predict.R:383-405`: `K`,
`gp_cross_cov(th, bk, Xc, pos)` and `Xr`. The conditional
cross-covariance between two prediction points is
`k(x1, x2) - Xr1 K Xr2'`. What is missing is that `frm_lp_basis()`
collapses it to a length-`n` vector (`man/frm_lp_basis.Rd`,
`extra_var`), and the plan forbids an extension reaching core with
`:::`.

**Not fixed.** Narrowing needs either a core seam (`frm_lp_basis()`
returning the cross-covariance, or an `extra_cov(nd1, nd2)` accessor)
or a name test on `coef_names` to find the block columns. Both are
larger than a review edit. This belongs in the plan's "Core seams"
list with item 1.3 named as the waiting consumer; the findings file
records it as a closed question instead.

**Blast radius, measured.** I suspected the refusal would also fire on
an ordinary random-effect difference and it does not:
`y ~ fac + s(x, by = fac, k = 6) + (1 | g)` at `re.form = NULL` with an
in-sample level gives `extra_var` all zero and the difference goes
through, with the same standard error as `re.form = NA` because the
group effect cancels. Exact `gp()` at unseen positions is the only
reachable source.

### 3. MEDIUM. The refusal names a cause that cannot occur

The message, `man/frm_curve.Rd`, and `curve-cov.R`'s header comment all
say the variance is "an exact `gp()` kriging variance or a grouping
level the fit never saw". The second cannot happen.
`frm_curve()` has no `allow_new_levels` argument and `sp_one_basis()`
calls `frm_lp_basis()` without one, so a new level errors in core
first. Measured: `frm_curve(f2, newdata = ndnew, re.form = NULL)` gives
"New levels in grouping factor `g`: NEW. Use allow_new_levels = TRUE to
predict them at the population level".

Not fixed: it is the same prose in four places, including roxygen the
author should own.

### 4. MEDIUM. The standard-error claim is an artifact of taking the worse curve

Recorded, in `dev/diffcurve-findings.md`, in `NEWS.md`, and in the
`s()` compat row: the difference's standard errors agree with gratia to
1.6 percent "against 3.1 percent for either curve alone", and the
improvement is "the part of that disagreement common to both curves
cancelling".

Reproduced construction (n = 400, seed 3, `k = 8`, 25 points on
[0.02, 0.98], mgcv ML; log-likelihood -167.640105 both ways).
`frmtmb .se / gratia .se`:

| | conditional (`Vp`) | `unconditional = TRUE` |
|---|---|---|
| curve A | [0.9961, 1.0149] | [0.9885, 0.9955] |
| curve B | [0.9962, 1.0311] | [0.9842, 0.9952] |
| difference | [0.9961, 1.0156] | [0.9868, 0.9954] |

The difference's worst deviation, 1.56 percent, sits BETWEEN the two
curves' 1.49 and 3.11, not at half of them, and it does the same under
`unconditional = TRUE` (1.32 between 1.15 and 1.58). "Either curve
alone agrees to 3.1 percent" is false of curve A, which agrees to 1.49.
Nothing the measurement can see cancels; the difference tracks the
better curve.

The test itself is not affected: it bounds the difference's gap by
`2 * gap_se` where `gap_se` is the max over both curves, so it is a
loose bound and it holds. What needs correcting is the prose and the
compat note.

### 5. LOW, the subtlest one. The redundancy buys one ulp

The lane's reversion reproduces exactly. Four arms, each installed into
its own library, same fixture:

| arm | one-grid `cov_rel_error` | `tol = 0` | `V` rows 1 and 2 swapped |
|---|---|---|---|
| `diag(A V A')` | 2.22044604925031308e-16 | REFUSES | 0.101, refused |
| `rowSums((A %*% V) * A)` | exactly 0 | RETURNS | 0.101, refused |

`test-curve.R` under the core spelling: PASS=61 FAIL=1, the one failure
being "a tolerance no covariance could meet refuses rather than
returning". Exactly as recorded.

But the recorded reason overstates it. `sp_one_basis()`'s comment and
the findings file say the triple product is "an independent arithmetic
route to the same number, and that independence is what makes the
agreement evidence rather than a tautology". Both routes read the SAME
`A`, `V` and `extra_var` out of the same `frm_lp_basis()` call, so the
only independence between them is summation order, and it is worth
exactly `.Machine$double.eps`. The permuted-`V` arm is the proof: a
real misreading of the seam is caught at 0.101 relative by BOTH
spellings. What core's spelling actually costs is the ability of
`tol = 0` to mean anything, which is a real cost and is why the
reversion was right to undo, but it is not the power to catch a wrong
read.

The general hazard the lane names is real and worth keeping in the
record. The sentence that states it should say what it buys.

Not fixed: rewording a comment the author reasoned their way to is
their call, and the decision it defends is the right one either way.

### 6. LOW (two of three fixed). Absolute numeric tolerances in new tests

`dev/diffcurve-findings.md`: "No absolute numeric threshold is written
anywhere in the file." Measured, three:

* `test-difference.R:80`, `expect_gt(gap_se, 0.005)`.
* `test-difference.R:307`, `expect_lt(attr(dif, "check")$cov_rel_error,
  1e-10)`.
* `test-span.R`, `expect_equal(ft$.estimate, 0.4, tolerance = 0.1)`.

**Fixed two.** The `cov_rel_error` bound is now
`10 * max(<one-grid cov_rel_error>, .Machine$double.eps)`, a ratio to
the same check run on one grid of the same fit with a machine floor, so
it stays portable across BLAS implementations where a bare `1e-10` is
arbitrary in one direction and `8 * eps` would be brittle in the other.
The crossing is now judged in units of the standard error the run
reports for it, `|.estimate - 0.4| < 6 * ft$.se`: measured 0.0247
against 0.113.

**Left.** `expect_gt(gap_se, 0.005)` guards `2 * gap_se` from being a
vacuous bound. What it should be a ratio to is a judgment about what
the guard is for, and I am not the author of that intent.

### 7. LOW (fixed). `print.frmtmb_feature()` claimed a check that never ran

On the nonlinear fixture, `print()` of a feature said

```
  covariance checked against predict(se.fit = TRUE) to NA relative
```

and the lane's new difference branch said

```
  each grid checked against predict(se.fit = TRUE) to NA relative; the
  difference itself has no second route
```

`predict(se.fit = TRUE)` is refused for a nonlinear predictor, so
nothing was checked. `print.frmtmb_curve()` has the NA guard and says
so correctly; `print.frmtmb_feature()` did not. Pre-existing for one
grid, duplicated into the difference branch. The `nl` compat row's
"print() says so" was true of `frm_curve()` and false of
`frm_curve_feature()`.

Fixed with the same NA guard, so the row is now true of both.
`test-message-uniqueness.R` collects `stop`, `warning`, `message` and
the two condition constructors, not `cat`, so the shared wording does
not trip it; it still passes 4/0.

### 8. LOW (fixed). `D2` at `curve-feature.R:305`

Confirmed dead. `D1` feeds the stationary branch through `qf(D1)`, `C0`
feeds the crossing branch and `.value_se`, and `D2` appears nowhere
else in the file. It is an `nr x p` product formed and discarded on
every feature call. Pre-existing at 780dec1, not the lane's. Removed;
no behavior changes and no test moves.

### 9. LOW (fixed). The help page stated a non-default agreement

`@section A difference curve:` said "It is the quantity
`gratia::difference_smooths()` reports for a factor-by smooth, and on
the same mgcv fit the two agree." Measured on the fixture, gratia's
default `group_means = FALSE` is 0.5219 away, because it zeroes the
intercept and the parametric `fac` columns. A reader following the help
page against gratia's defaults sees a large disagreement that is not an
error in either package. The lane recorded this under "Found, not
fixed" and asked whether the help page should say it. It should: the
sentence as written is false at the defaults. Now names
`group_means = TRUE` and says what the default reports.

### 10. LOW. The API note overstates the compatibility cost

`contrast` is the third positional argument of `frm_curve()`. The
findings say an out-of-tree caller passing a third positional argument
"would silently reinterpret it". Measured, it does not:

```
frm_curve(fit, gA, "sigma")                       -> ERROR
frm_curve_deriv(fit, "x", 1L, gA, "sigma")        -> ERROR
frm_curve_feature(fit, "x", "crossing", 0, gA, "sigma") -> ERROR
```

all three raising "`contrast` must be a data frame with at least one
row: it is the second grid the curve is differenced against". The old
third positional was `dpar`, a character, and `sp_check_contrast()`
refuses a character by name. A silent reinterpretation would need a
data frame passed positionally where `dpar` was expected, which is not
a call anyone writes.

**On the placement itself:** keep it third. None of the three
signatures has a `...`, so "after `...`, so it must be named" is not
available without adding one, and a `...` added to force naming would
start swallowing misspelled arguments in exchange for a compatibility
risk the measurement shows is not there. `contrast` is the second most
important thing a caller passes and it reads correctly next to
`newdata`. Every call site in the repository passes by name; I checked
`R/`, `tests/`, both vignettes and `docs/`.

One side effect worth knowing: partial matching now resolves `c = `,
`co = ` and `con = ` to `contrast`, where they used to be
unused-argument errors.

### 11. LOW. The degenerate-band refusal, as recorded

Reproduced verbatim. `frm_curve(newdata = g, contrast = g,
simultaneous = TRUE)` refuses with "every point on this grid has a
standard error of exactly zero ... Past a `ps()` term's knot span the
basis is exactly zero and so is the derivative design, which is the
usual way to arrive here". First sentence true, second is a hint that
does not fit a self-difference. `test-span.R` matches the template.
Left alone, as the lane did.

## What I verified and found correct

### The span gate, by construction

* `sp_grid_pinned(nd, "t")` TRUE and `sp_grid_pinned(ct, "t")` TRUE on
  the fixture, so no extra prediction is made and the gate cannot fire
  for the reason a difference exists.
* The hazard is real: `sp_grid_pinned(rbind(nd, ct), "t")` is FALSE. An
  implementation that stacked the two grids would refuse a difference
  curve entirely inside every span.
* Nothing in tree is that implementation. The only `rbind()` of grids in
  `extensions/frmtmb.spline/R/` is `curve-deriv.R:126` and `:130`, and
  each stacks a grid with ITSELF, three copies with `var` shifted, never
  `nd` with `ct`. Even that stack stays pinned in every column but
  `var`, so a gate asked there would not fire either.
  `sp_grid_pinned()` is called at exactly two places,
  `curve-feature.R:282-283`, once per grid.
* Core has no analogue of the predicate and no `rbind()` in
  `R/predict.R` or `R/conditional-effects.R`.

So the gate is one refactor away from firing on correct models, exactly
as the lane says, and no such refactor exists yet.

The out-of-span refusal fires from either side, correctly tagged;
`test-span.R:441` is the assertion, and it fails when the gate is asked
of `newdata` only (below).

### The `sp_span_fit3()` fixture, and its crossing, by another route

Reproduced independently in a fresh process: t span [-0.2970, 1.2981],
z span [-0.0197, 1.0192], `b1` = -0.5390 (0.0484), `b2` = 1.4363
(0.0839).

The recorded crossing is 0.3753 (se 0.0188). The difference across `w`
is `b1 + b2 t`, so its root and the delta-method standard error can be
computed from `coef()` and `vcov()` alone, with no part of
`frm_curve_feature()` involved:

| route | root | se |
|---|---|---|
| `frm_curve_feature()` | 0.375296 | 0.018779 |
| closed form from `vcov()` | 0.375296 | 0.018779 |

Same to six decimals. The fixture is sound and the number is right.

### The three failing-first results

All three reproduce exactly, each on its own patched install:

| reversion | file | result |
|---|---|---|
| gate on `newdata` only | test-span.R | PASS=68 FAIL=1, line 441 |
| `sp_spec()` drops the contrast | test-difference.R | PASS=32 FAIL=5 |
| core's spelling of the diagonal | test-curve.R | PASS=61 FAIL=1 |

and the `sp_spec()` arm prints the recorded line: "Expected
`max(abs(d1$.estimate - a1$.estimate))` > `0.1 * max(abs(a1$.estimate))`.
Actual comparison: 0.00 <= 0.69".

### The gratia identity, with the attribution PROVEN rather than argued

| quantity | measured |
|---|---|
| `logLik(fit)` against `-gm$gcv.ubre` | -167.640105 both |
| max abs estimate gap, curve A | 7.892798e-09 |
| max abs estimate gap, curve B | 1.291884e-08 |
| max abs estimate gap, the difference | 1.954472e-08 |
| ratio | 1.5129 |

The lane attributes the 1.51 to the fit difference rather than the
difference arithmetic. That is not just plausible, it is exact:

* `max |gapD - (gapA - gapB)| = 5.551115e-16`, which is 2.84e-08 of the
  difference gap. The difference's disagreement with gratia IS the
  difference of the two curves' own disagreements.
* `max |dif$.estimate - (a$.estimate - b$.estimate)| = 0` exactly, on
  the frmtmb side.
* `max |gd$.diff - (pa$fit - pb$fit)| = 4.440892e-16` on gratia's.

So the difference arithmetic contributes nothing to the gap in either
package. The attribution stands.

The standard-error half of the same paragraph does not: see finding 4.

### Coverage over 200 seeds

Re-measured in a fresh process, same construction, 28.3 s (lane
recorded 27.4 s):

| quantity | lane | mine |
|---|---|---|
| simultaneous, whole curve | 0.970 | 0.9700 |
| pointwise, whole curve | 0.685 | 0.6850 |
| pointwise, per point | 0.956 | 0.9560 |
| mean simultaneous critical value | 2.818 | 2.8177 |
| binomial mcse at 0.95 | 0.0154 | 0.0154 |

Assertion margins: 0.9700 against a bound of 0.9038; 0.6850 against
0.8775; |0.0060| against 0.0274.

**0.970 against 0.95 in mcse units:** 1.30 at the nominal rate, 1.65 at
the measured rate (0.0121). Over-coverage that 200 seeds cannot
separate from nominal. That does not weaken the test, because the
assertion is one-sided and asks only that under-coverage is absent, and
the direction it errs in is the direction a simultaneous band is
allowed to err.

**Are the assertions ratios to a run-measured spread?** The third is:
`sd(per_point) / sqrt(200)`, measured across replicates, correctly not
binomial. The first two use `sqrt(lev * (1 - lev) / n_seed)`, computed
from the nominal level and the seed count. That is the textbook
quantity and not a magic constant, but "all multiples of a spread the
run measures" is not exactly what they are. No absolute threshold in
any of the three.

**Is 200 enough?** The loop is fully deterministic: `set.seed(s)` per
replicate and `seed = s` into the band, so all three numbers reproduce
bit for bit and the mcse describes the estimator, not run-to-run risk.
The exposure is optimizer or RNG drift across R versions, and the
margins (0.066, 0.19, 0.021) are wide enough to absorb a good deal of
it.

**Cost.** 28 to 60 s depending on machine load, the slowest file in the
package by about 4x, `skip_on_cran()`. Acceptable now; it is about 60
percent of the package's test wall clock, so if it grows it should move
behind a gate rather than get shorter.

### The two identities the check cannot give

Both reproduce, both `identical()`:

* `y ~ fac * x`, no random effect. `cv$.se` is
  0.13190741863121563, 0.09318196316020555, 0.14082170654557177 and the
  `vcov()` route gives the same three doubles;
  `identical(cv$.se, se)` TRUE.
* A grid differenced with itself: `identical(z$.estimate, rep(0, 25))`
  and `identical(z$.se, rep(0, 25))` both TRUE, and both print as
  0.00000000000000000e+00 at 17 digits, not as a rounded zero.

Quadrature over the difference on the by-factor fixture: my first run
at n = 250 gave [0.99978, 1.00073], not the recorded [0.99983,
1.00016]. The construction is the n = 400 fixture, where I get
[0.99983, 1.00016] exactly. Record found, not a false number.

### The e2 stencil defect

Replicated with `e2` scaled by an option so `e1` stays put, which is
what the public `eps` cannot do. On the difference's crossing, a
hundredfold change in `e2`:

| reported | identical | value |
|---|---|---|
| `.se` | TRUE | 0.023344393861040272 |
| `.value` | TRUE | 4.4408920985006262e-16 |
| `.value_se` | TRUE | 0.075175120356382269 |
| `slope` | TRUE | 3.2202643942639648 |
| `curvature` | FALSE | 2.1959008926409398 then 2.1936734673028222 |

Every figure the lane records, to the last digit.

I also measured the other side, which the lane did not: on a STATIONARY
point of the same fit the defect is live. `.estimate` moves
0.53736407036424594 to 0.537364070351386 and `.se` moves
0.023069714372094628 to 0.023081239769917671, about 5e-4 relative under
the same hundredfold change. So the confinement claim holds from both
directions, the recorded defect in `dev/debts-findings.md` stays open,
and it is not this item's.

### The compat rows

Ten rows registered under `frm_curve_contrast` and all ten render:
`frm_compat("frm_curve_contrast")` returns 112 rows, of which the ten
carry the lane's notes and the other 102 are the auto-filled defaults.
Statuses: `smooth` works, `gp` refused, `s()` works, `t2()` untested,
`ps()` works, `rr` untested, `mvbf` untested, `nl` works, `predict`
works, `frm_lp_basis` works.

No row promises something a refusal contradicts. Three carry claims the
measurement does not support, all covered above: the `s()` row's
standard-error sentence (finding 4), the `gp` row's "there is no
covariance between the two grids to difference it against" (finding 2),
and the `nl` row's "print() says so", which was half true until finding
7 was fixed.

The `nl` row is otherwise solid: on the two-`ps()` nonlinear body the
crossing of the difference agrees with the closed-form delta method to
six decimals, and `cov_rel_error` is `NA` with `n_predict` 0, for the
same reason it is on one curve.

`frm_curve_contrast` is registered as a feature but does not appear in
`frm_compat_features()`; neither does `frm_curve`, so that is how
package-registered features behave, not something this lane broke.

### `R CMD check --as-cran`, verified against the base rather than accepted

Both arms built from a tarball, each in its own check library, with
pandoc and TinyTeX on PATH.

| arm | status |
|---|---|
| 780dec1, built into `rvdc-main-lib` | 1 WARNING, 1 NOTE |
| the lane, as received | 1 WARNING, 1 NOTE |
| the lane, after my fixes | 1 WARNING, 1 NOTE |

The three `00check.log` bodies are identical line for line apart from
timings and the check directory path. So the WARNING is pre-existing,
as the lane says, and now that is measured: CRAN incoming feasibility
reporting "New submission", "Strong dependencies not in the CRAN or
BioC software repositories: frmtmb", and a 301 on the `DESCRIPTION`
URL, which wants a trailing slash. The NOTE is the expected V8 one.
`checking for code/documentation mismatches ... OK`, `checking Rd line
widths ... OK`, `checking examples ... OK`, vignettes rebuild clean.

### House style

Every added line in `R/`, `NEWS.md` and both test files is within 80
columns. The only over-80 lines are the compat note strings in `zzz.R`,
which is the existing convention there (24 such lines at 780dec1). No
em dash, no en dash, no emoji, no spaced hyphen standing in for an em
dash; every ` - ` in the added lines is R subtraction. Comments say why.
No `@noRd` block is followed by loose text and no two roxygen blocks
run together, which `R CMD check` confirms by building the manual
clean.

## Suspicions I chased and disproved

* **That `sp_span_both()` could swallow one of its two messages**,
  because it ends in `unique()`. It cannot: the contrast's message is
  tagged with "In `contrast`: " BEFORE the `unique()`, so two identical
  span messages stay two entries. Measured directly on the helper.
* **That the `extra_var` refusal would fire on an ordinary
  random-effect difference at `re.form = NULL`**, which would have made
  it very intrusive. It does not. `lp_extra_var()` only produces an
  entry for a row whose level index is NA, so an in-sample level gives
  `extra_var` all zero. Measured on
  `y ~ fac + s(x, by = fac, k = 6) + (1 | g)`: the difference goes
  through at both `re.form = NA` and `re.form = NULL`, with the same
  maximum standard error 0.11142, because the shared group effect
  cancels in the difference.
* **That core's spelling of the diagonal would let a real seam misread
  through**, which is what "the check becomes vacuous" reads like. It
  does not. With `V`'s first two rows and columns swapped, both
  spellings refuse at 0.101 relative. Only `tol = 0` is lost. That is
  finding 5.
* **That the quadrature range [0.99983, 1.00016] was wrong**, because
  my first attempt gave [0.99978, 1.00073]. I had used the n = 250
  fixture; the record is the n = 400 one, and it reproduces exactly.

## Files I changed in the lane worktree

* `extensions/frmtmb.spline/R/curve-deriv.R`: the `sp_spec()` refusal
  (finding 1) and the Rd section that documents it.
* `extensions/frmtmb.spline/R/curve-feature.R`: the same Rd section,
  the `print.frmtmb_feature()` NA guard (finding 7), and the dead `D2`
  removed (finding 8).
* `extensions/frmtmb.spline/R/curve.R`: the `group_means = TRUE`
  sentence in `@section A difference curve:` (finding 9).
* `extensions/frmtmb.spline/man/frm_curve.Rd`,
  `frm_curve_deriv.Rd`, `frm_curve_feature.Rd`: roxygenised, 8.1.0.
* `extensions/frmtmb.spline/tests/testthat/test-difference.R`: the new
  test block for finding 1, and the `cov_rel_error` bound made a ratio
  (finding 6).
* `extensions/frmtmb.spline/tests/testthat/test-span.R`: the crossing
  judged in units of its own standard error (finding 6).

Nothing else. `NEWS.md`, `zzz.R` and `curve-cov.R` are as the lane left
them, and findings 2, 3, 4, 5, 6's third item, 10 and 11 are left for
the author.

## What consolidation should do before this ships

1. Correct the standard-error sentence in `NEWS.md`, the `s()` compat
   row and `dev/diffcurve-findings.md` (finding 4). The number 3.1
   percent is curve B's; curve A agrees to 1.49 and the difference to
   1.56.
2. Drop the unreachable "grouping level the fit never saw" from the
   `extra_var` refusal, its Rd and the `curve-cov.R` comment
   (finding 3).
3. Rewrite the `gp` refusal's justification to say what is true: the
   public seam returns no cross-covariance, so an extension cannot form
   the term, and where the two grids carry the same `gp()` positions it
   is exactly zero (finding 2). File the core seam it waits on.
4. Say what the triple product buys, which is one ulp and a `tol` that
   can be set to zero, rather than arithmetic independence
   (finding 5).
5. Decide what `expect_gt(gap_se, 0.005)` should be a ratio to
   (finding 6).

---

# Re-check, 2026-09-08 (second pass): the `gp()` narrowing

The lane closed review finding 2 with new code rather than a wording
change. This section reviews that code, and the two other things the
coordinator flagged. Same discipline: private libraries only, one test
file per R process, nothing committed.

## Verdict for the re-check

**Mergeable.**

The narrowing is sound, the predicate is the right one, and the design
test is better than the lane's own evidence claims. Two of the lane's
recorded reversion counts are stale in a direction that UNDERSELLS it,
and both are corrected below.

Two things I would still change, neither blocking:

* **R6, MEDIUM.** The predicate is stricter than the shipped claim. On
  `y ~ fac + s(x, by = fac, k = 6) + gp(x)` the `gp()` columns are
  bit-identical and the residual demonstrably cancels, but the
  by-factor smooth's columns differ, so the call is refused. That is
  one step from the fixture this whole item is built on, and `NEWS.md`
  says it works. It fails closed, so it is a missing feature and not a
  wrong answer.
* **R9, HIGH for the package but not for this lane.** The
  `sp_sim_crit()` defect is real, shipped in 0.4.0, and at least 17
  percent on a grid that extrapolates past the observed `gp()`
  positions. It needs a NEWS bullet and a backlog entry.

## R1. The narrowing is sound, and the identity reproduces

`frm_joint_cov(fit)$names` is a documented return value, and documented
as exactly what this predicate needs: "the PARAMETER COMPONENT each row
belongs to (`beta`, `betad`, `b`, `theta`, ...), which is what
`frm_lp_basis()$coef_pos` indexes" (`man/frm_joint_cov.Rd`). No name
parsing, no reach into core. The claim holds.

The flat-difference identity, reproduced in a fresh process on
`y ~ fac + gp(x)`, n = 90, grid at the 89 observed midpoints:

| quantity | measured |
|---|---|
| `sd(.se) / mean(.se)` | 0 exactly |
| `identical(.se, rep(.se[1], n))` | TRUE |
| `.se[1]` | 0.0618114534972 |
| `sqrt(vcov(fit)["facB", "facB"])` | 0.0618114534973 |
| relative gap | 1.540e-12, 6934 ulps |
| spread of `.estimate` | 3.331e-16 on -0.4516 |

The findings record the estimate spread as 2.2e-16, about 4 ulps; I
measure 3.331e-16, about 6. The assertion is a ulp ratio with a 16-ulp
bound, so both pass and neither number is load-bearing.

`sp_same_latent()` also fails closed everywhere it cannot answer: a
component vector that does not line up with the design, or no non-fixed
column at all, returns `FALSE`. `theta` columns count as non-fixed and
so must match too, which is stricter than needed and safe.

## R2. `vcov()` is an independent route, not a restatement

This was worth checking, and it survives.

| quantity | measured |
|---|---|
| `dim(frm_joint_cov(fit)$V)` | 95 x 95: b 90, beta 2, betad 1, theta 2 |
| `dim(vcov(fit))` | 3 x 3 |
| `identical(jc$V[facB, facB], vcov(fit)["facB", "facB"])` | FALSE |

`vcov.frmtmb_fit()` on a non-REML, non-profile fit reads
`sdr_of(object)$cov.fixed` (`R/methods-fit.R:448`), the inverse of
`sdreport()`'s OUTER Hessian. `frm_joint_cov()` inverts the JOINT
precision and the extension takes a block of it. They are two different
matrices from two different solves, equal only through the Schur
complement identity, and they are NOT bitwise equal: they differ by
1.54e-12.

That is the signature of an independent route. Where the two routes
really are the same arithmetic, this package's own fixed-effect-only
test gets `identical()`; here it gets 6934 ulps, which is what a 95 x 95
solve against a 3 x 3 Hessian inverse costs. So the identity is
evidence, and the ulp-count bound in the test is the right shape for it.

## R3. The mirrored `gp()` construction is real, and reproduces exactly

`seq(-5, 5, length.out = 41)`, grids `x` and `-x` at
`c(0.125, 1.125, 2.125)`:

| quantity | measured |
|---|---|
| `max abs(extra_var1 - extra_var2)` | 1.1102e-16 |
| `max abs(extra_var1 / extra_var2 - 1)` | 1.2041e-10 |
| the variance being compared | 9.22e-07 to 9.26e-07 |
| `isTRUE(all.equal(extra_var1, extra_var2))` | TRUE |
| `identical(extra_var1, extra_var2)` | FALSE |
| non-fixed columns identical | FALSE |
| `max abs(A1 - A2)` over them | 3.7314e-01 |
| `frm_curve(contrast = )` | refuses |

Every figure the lane records. The case is not a curiosity: computing
the kriging cross-covariance from the fit's own kernel gives residual
correlations of 0.2117, -0.0171 and 0.0049, so the two grids really do
load different draws, and the omitted `var(g1 - g2)` is 1.45e-06 to
1.88e-06, an omitted standard error of about 1.2e-03 to 1.4e-03. A
tolerant gate does not just accept a case it should refuse; it drops a
term that exists.

## R4. Both recorded reversion counts are stale, and both undersell the lane

Each reversion built as the findings describe it, installed into its own
private library, `test-difference.R` run in its own process:

| gate | recorded | measured |
|---|---|---|
| `all.equal(extra_var1, extra_var2)` gate | PASS=57 FAIL=2 | PASS=56 FAIL=3 |
| `identical(extra_var1, extra_var2)` gate | PASS=57 FAIL=0 | PASS=58 FAIL=1 |
| shipped design test | PASS=59 FAIL=0 | PASS=59 FAIL=0 |

All three total 59, which is what the shipped file asserts, so the file
is the same one; the counts are not.

The extra failure is the same assertion in both arms:
`expect_false(sp_same_latent(...))` at `test-difference.R:424`, the
two-levels-of-one-grouping-block case. The tolerant arm also lets the
mirrored grid and the 1e-10 shift through, exactly as recorded.

This matters for the lane's own conclusion. The findings say "the
design test is provably better than the tolerant one and only luckily
different from the bitwise one", on the basis that the bitwise gate
passes the suite. It does not: it fails, on an assertion the lane wrote
deliberately. The honest statement is stronger. The design test and the
bitwise numeric test agree on everything reachable through
`frm_curve()` today, and the suite separates them anyway, because the
grouping-block case reaches the predicate directly. The design test is
strictly better than both numeric tests and the file proves it.

## R5. The grouping-block case is unreachable, verified from both sides

Not accepted, constructed.

* Forward: `lp_extra_var()` (`R/predict.R:1605`) creates an entry only
  for rows whose level index is NA (`if (!length(nas)) next`), that is,
  a level the fit never saw. A new level needs `allow_new_levels = TRUE`
  in `frm_lp_basis()`. `formals(frm_curve)` is `object, newdata,
  contrast, dpar, resp, re.form, level, simultaneous, nsim, transform,
  seed, tol` and `sp_one_basis()` passes no such argument, so core
  refuses first: measured, `frm_curve(fit, newdata = <new level>,
  re.form = NULL)` gives "New levels in grouping factor `g`: NEW. Use
  allow_new_levels = TRUE to predict them at the population level".
* Backward: an in-sample level carries no extra variance at all.
  Measured on `y ~ fac + s(x, k = 6) + (1 | g)` at `re.form = NULL`,
  `all(extra_var == 0)` is TRUE, so the guard never runs and the
  difference across two levels goes through.

"Evidence about the predicate rather than a live hazard" is right.

## R6. MEDIUM. The predicate is stricter than the claim, on the item's own model

`sp_same_latent()` requires EVERY non-fixed column of `A` to match, not
only the columns of the block that can carry `extra_var`. Measured on
`y ~ fac + s(x, by = fac, k = 6) + gp(x)`, n = 120, grid at the
observed midpoints, contrasting `fac`:

| quantity | measured |
|---|---|
| `identical(extra_var1, extra_var2)` | TRUE |
| `gp()` columns of `A` identical | TRUE |
| all non-fixed columns identical | FALSE |
| non-fixed components present | `b`, 128 columns |
| `frm_curve(contrast = )` | REFUSES |

The kriging residual cancels there, provably: the `gp()` columns are
bit-identical, so `A1 - A2` is exactly zero on them. The answer is
available and is not returned, because the by-factor smooth's own `b`
columns differ, which is what a by-factor smooth is for.

`NEWS.md` says "Works when the two grids sit at the SAME `gp()`
positions, which is the ordinary case". That is not what the code does:
it also requires every other latent column to agree. A user who adds a
`gp()` nugget to the `s(x, by = fac)` model this item is built around
gets a refusal whose message is about `gp()` positions that are in fact
identical.

Not fixed. Narrowing the kept columns to the block that carries
`extra_var` needs the block identity, and the only public route to it is
the `b.<block>.<level>` shape of `frm_joint_cov()$labels`, which is the
name parsing the lane deliberately avoided. The real fix is the same
core seam: a matrix-valued `extra_cov` makes the predicate unnecessary.
What I would do now is one sentence in `NEWS.md` and in the Rd naming
the restriction, and this case added to the seam's waiting list.

## R7. The two-dimensional escape: real, and not reachable with a healthy fit

The stated assumption is that the kriging weights fix the position
unless every observed position is equidistant from both grid positions,
which in two dimensions needs the observations collinear and the grids
mirrored across that line. I built exactly that: `gp(u, v)`, 60
observations all at `v = 0`, grids at `v = +0.75` and `v = -0.75`,
`u` in `c(2.3, 4.7, 7.1)`.

The escape is real and it is large:

| quantity | measured |
|---|---|
| `identical(extra_var1, extra_var2)` | TRUE |
| non-fixed columns identical | TRUE |
| `max abs(A1 - A2)` over them | 0 exactly |
| `var(g)` | 2.5136e-01 |
| `cov(g1, g2)` | -1.4322e-01 |
| omitted `var(g1 - g2)` | 7.8917e-01 |
| omitted standard error | 0.888 |

So `sp_same_latent()` does return TRUE on that design. But the design
cannot be fitted: frmtmb's `gp()` carries one length scale per
dimension, so observations on a line do not identify the scale across
it, and the fit came back with "The joint precision matrix is singular
... theta_3 - zero gradient and an empty Hessian row". `frm_curve()`
refused, at the covariance check rather than at the latent gate.

And the escape does not survive contamination. Put three observations
off the line: the fit is healthy (`all(is.finite(frm_joint_cov(f2)$V))`
TRUE), and `identical(extra_var)` is FALSE, `max abs(A1 - A2)` over the
non-fixed columns is 3.5244, and the difference is refused correctly.

The two conditions are therefore in tension: the escape needs every
observation equidistant from both grid points, which in two dimensions
needs them all collinear, which is what stops the kernel being
identified. That is a second line of defence rather than a designed
one, and it is not airtight: a fit that pins the second length scale by
a prior or a bound would be non-singular and would take the escape.
Worth one sentence in `sp_same_latent()`'s block, since the block
currently says the term "would be understated rather than refused" and
the measurement says the fit refuses first.

## R8. The core seam filing is accurate, and its preferred shape is right

Checked against `R/predict.R:383-405`. The transcription is exact: `th`
from `fit$estimates[["theta"]][bk[["theta_idx"]]]`, `K` from the `gp`
registry `vcov()`, `Ks <- gp_cross_cov(th, bk, Xc, pos)`,
`Xr <- t(solve(K, t(Ks)))`, `kss <- exp(2 * th[1]) * (1 + 1e-6)`,
`extra_var <- pmax(kss - rowSums(Xr * Ks), 0)`.

The missing quantity is stated correctly.
`Cov(g_i, g_j | obs) = k(x_i, x_j) - Ks_i K^-1 Ks_j'`, and
`Xr_i K Xr_j' = Ks_i K^-1 K K^-1 Ks_j' = Ks_i K^-1 Ks_j'`, so
`k(x_i, x_j) - Xr_i K Xr_j'` is right, and its diagonal is
`kss - rowSums(Xr * Ks)`, exactly what core returns.

The preferred shape, `frm_lp_basis(newdata, extra_cov = TRUE)`, is the
right one, for a reason the filing does not give: a DIFFERENCE gets its
cross term from ONE call on the stacked grid `rbind(newdata, contrast)`,
reading the off-diagonal block of the `2n x 2n` result. Without that
sentence a reader may conclude the `n x n` shape serves only one grid
and that the second option is needed for a difference. It is not; say
so, because it is what makes the preferred shape sufficient rather than
merely convenient.

One gap. The filing names only the `gp()` source. `extra_var` also
carries a new grouping level's marginal variance, and
`extra_var_blocks()` (`R/predict.R:1660`) already decides which rows
share a draw through its `new_key` grouping. A seam that returns a
matrix should carry that too, or `predict()`-side consumers get an
object that is a covariance for `gp()` and a diagonal for new levels.
One sentence.

## R9. HIGH for the package. The `sp_sim_crit()` defect is real and shipped

**The reading is confirmed.** `curve.R:276` calls
`sp_sim_crit(Sigma, se, nsim, level, seed)` with `Sigma = A V A'` and
`se = sqrt(diag(A V A') + extra_var)`. The drawn deviation process has
marginal standard deviation `sqrt(1 - r_j)` with
`r_j = extra_var_j / se_j^2`, not 1, so its maximum and the quantile
read off it are too small and the simultaneous band is too narrow.

**It is pre-existing.** Identical numbers from a build of 780dec1 in
its own library. `frmtmb.spline 0.4.0` in both arms.

**It can be bounded without the missing seam.** `extra_var` comes from
`frm_lp_basis()` and `se` from `predict(se.fit = TRUE)`, both public,
so `r_j` is computable. The scale part of the correction lies between
`1 / sqrt(1 - min r)` and `1 / sqrt(1 - max r)`, and running the same
simulation re-standardized by `sqrt(diag(Sigma))` gives it directly.
On `y ~ fac + gp(x)`, 60 observations on [0, 6], noise 0.2,
`nsim = 20000`, seed 1:

| grid | `r` range | crit shipped | crit rescaled | ratio |
|---|---|---|---|---|
| inside the data, [1, 5] | [0.00008, 0.00025] | 2.76566 | 2.76591 | 1.00009 |
| at the edge, [5.5, 6.5] | [0.00010, 0.00288] | 2.37258 | 2.37275 | 1.00007 |
| extrapolating, [7, 12] | [0.01649, 0.70908] | 2.01622 | 2.36877 | 1.17485 |

So the band is under a hundredth of a percent too narrow on a grid
inside the data, and at least 17 percent too narrow on one that
extrapolates past the observed `gp()` positions, which is exactly where
a reader expects a band to widen rather than shrink. The ratio is a
LOWER bound on the correction: it fixes the marginal scale but keeps
`A V A'`'s correlation, and a kriging residual that decorrelates faster
than the mean function pushes the critical value higher still.

**Not fixed.** The correct fix needs the matrix the seam does not
return. Two things are available now and both are the author's call: a
partial mitigation that standardizes the draw by `sqrt(diag(Sigma))`
instead of `se`, which is right in scale and still misses the
correlation, so it is strictly closer than what ships but changes a
released number; or a refusal or a warning on `simultaneous = TRUE`
when `max(r)` is large, with `r` measured from public quantities as
above.

The coordinator is right that this needs its own `NEWS.md` bullet and a
backlog entry rather than a paragraph inside a core-seam filing. It is
shipped in 0.4.0, it is a ONE-GRID defect with nothing to do with a
difference, and a reader of the seam filing will not find it.

## R10. LOW, fixed. Stale `@noRd` text

`sp_curve_parts()`'s block still listed, as one of the three things a
difference must establish, that `extra_var` "cannot be differenced and
is refused". The narrowing made that false. Rewritten to say the cross
term is exactly zero where the two grids load the same draw, which
`sp_same_latent()` decides, and refused otherwise. `@noRd`, so no Rd
regeneration and no roxygen run.

## R11. What the coordinator asked me to confirm

* `NEWS.md` and the `s()` compat row now say what I measured: the
  difference's 1.56 percent "sits BETWEEN curve A's 1.49 and curve B's
  3.11 rather than below both", the same under `unconditional = TRUE`
  (1.32 between 1.15 and 1.58), and "nothing measurable cancels".
* The `gp` compat row is `conditional` and its note matches the code,
  including the seam filing and the `1e-10` sensitivity.
* The third absolute threshold is gone: `expect_gt(gap_se, 0.005)` is
  now `expect_gt(gap_se, 1000 * rel_est)`, a ratio between two things
  the run measures, and it says something true and load-bearing.
* `D2` is gone; the `group_means = TRUE` sentence is in; my `sp_spec()`
  refusal and the `print.frmtmb_feature()` NA guard are intact.

## Verification for the re-check

`frmtmb.spline`, shipped lane, one file per process, `NOT_CRAN=true`:
test-bracket-access.R 1/0, test-curve.R 62/0, test-deriv.R 38/0,
test-difference.R 59/0, test-gratia.R 13/0, test-message-uniqueness.R
4/0, test-royston-parmar.R 67/0, test-rp-floored.R 42/0, test-span.R
73/0, test-surface.R 48/0. Total 407 / 0, which is the number the
coordinator reports.

After my `@noRd` edit, test-curve.R 62/0, test-difference.R 59/0,
test-span.R 73/0.

`R CMD check --as-cran` on the built tarball: 1 WARNING, 1 NOTE, the
same pre-existing pair verified against 780dec1 in the first pass.

## A suspicion I chased and disproved, on myself

My first `R CMD check` of this round came back 1 WARNING and 2 NOTEs,
the new NOTE being "Examples with CPU (user + system) or elapsed time
over 5s: frm_curve 4.86 0.77 6.89". I nearly filed it as a regression.
It is not one. The lane changed no example: the `@examples` block of
`frm_curve()` is untouched in the diff, and so is the Rd's example
section. I had a probe and a second package install running against the
same machine while the check ran. Re-run with nothing else going on,
`checking examples ... OK` and the status is 1 WARNING, 1 NOTE. The
instrument, not the package. It is the trap the round rules put first,
and this time it caught the reviewer.

## Files I changed in the re-check

* `extensions/frmtmb.spline/R/curve-cov.R`: the stale `@noRd` sentence
  in `sp_curve_parts()` (R10). Nothing else.
