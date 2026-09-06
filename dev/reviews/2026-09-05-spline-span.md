# Review: wt-spline-span

Reviewer pass over branch `wt-spline-span` (worktree
`C:\Users\adf44\source\r\frmtmb-wt-spline-span`), branched from 68d6782,
all work uncommitted. Verdict and findings below; measurements are mine
unless marked as reproduced from the lane's own script.

Library: private `rvss-lib`, `frmtmb 0.52.0` and `frmtmb.spline 0.2.0`
both resolved from it (checked with `packageVersion()` and
`system.file()`).

## (a) Diff hygiene

23 paths, matching the claim: 21 modified, 2 untracked
(`dev/spline-span-findings.md`,
`extensions/frmtmb.spline/tests/testthat/test-span.R`).

- No `dev/*.log`.
- `git diff 68d6782 -- NAMESPACE extensions/frmtmb.spline/NAMESPACE
  DESCRIPTION extensions/frmtmb.spline/DESCRIPTION` is EMPTY. No
  NAMESPACE change, and no version bump on either package. See (g) for
  why the missing floor matters.
- `extensions/frmtmb.spline/R/curve-cov.R` is modified but was not named
  in the lane's report. It is where `sp_catch_span()` and
  `sp_span_stop()` live, so the omission is in the report, not the work.

(work in progress)

## (b) The collector design

`R/ps.R:286` passes `check` (FALSE / TRUE / environment) into the ps
closure; `R/ps.R:572` `ps_span_warning(x, pt, sink)` collects or warns;
`R/ps.R:598` `ps_span_flush()`; `R/ps.R:617` `ps_span_signal()` raises a
`frmtmb_ps_span_warning`. `R/predict.R:3226-3243` arms the collector,
`:3283-3285` runs the off-tape pass, `:3292` flushes.

### The dedup is load-bearing, not incidental

I counted closure entries with `trace("ps_span_warning", where =
asNamespace("frmtmb"))`. The body IS evaluated more than once, so the
lane's stated reason for collecting is the real one:

| call | closure entries | warnings raised |
| --- | --- | --- |
| fit A `frm_lp_basis(re.form = NA)`, 1 ps term | 2 | 1 |
| fit A `frm_lp_basis(re.form = NULL)`, 1 ps term | 2 | 1 |
| fit A `predict()` | 1 | 1 |
| fit B `frm_lp_basis()`, 2 ps terms | 4 | 2 |
| fit B `predict()` | 2 | 2 |

Fit A is `y ~ lev + ps(t + shift, k = 8, pad = 0.3)` with
`shift ~ 0 + (1 | id)`; fit B is two plain-data `ps()` terms. On the
basis route every term is entered exactly twice (off-tape pass, then the
tape build) and emerges as one warning. Warning in place would have
double-fired; the collector is doing real work.

### The tape guard

Confirmed by observation, not by reading. In fit A with
`re.form = NULL` the trace records the two passes as
`x=numeric` then `x=advector`; with `re.form = NA` both are numeric
(the dropped random effect leaves no estimated coefficient in the
`ps()` argument). So the advector case is real and reachable, and
`inherits(x, "advector")` at `R/ps.R:573` is what catches it. The old
`is.numeric()` guard would have let it through, and `<` on an advector
raises. The lane's account of the bug is accurate.

### The off-tape grid is the right grid

The off-tape call is `eta_fun(chat)` (`R/predict.R:3284`), the same
closure the tape is built from and closed over the same `dl <-
lp_basis_nl_data(object, newdata)`, so it is on the newdata grid by
construction. Checked anyway: the traced `x` values on the basis route
and the predict route are the same vector (`3 4 5` for both terms of fit
B), and `lb$eta` equals `predict()`'s value at `tolerance = 0` under
both `re.form = NA` and `re.form = NULL`.

### Cost of the off-tape pass, measured

Fit C, one `ps(t, k = 10)` term, `ncol(A) = 10`. Whole-call on/off
timing in isolated processes was inside the run-to-run noise, so I
bounded the pass directly: it is one nonlinear body evaluation, and one
`predict(type = "link")` on the same grid is an upper bound on it.

| grid N | `frm_lp_basis()` | one body eval | share |
| --- | --- | --- | --- |
| 200 | 37.2 ms | 4.4 ms | 11.7% |
| 1000 | 565.8 ms | 3.9 ms | 0.7% |
| 5000 | 7304.1 ms | 4.6 ms | 0.1% |

The pass is flat in N (a vectorized body evaluation) while the call is
dominated by the taped Jacobian, which is superlinear. The overhead is
acceptable and shrinks exactly where it would matter.

## (c) Text identity

Byte-identical on both routes, `identical()` on `conditionMessage()`,
with `re.form = NA` and `re.form = NULL`:

- fit A, 3-row newdata past the span: basis 1 warning, predict 1
  warning, `identical(msgs) = TRUE`, both routes carry class
  `frmtmb_ps_span_warning/warning/condition`.
- fit B with `t` past the span and `u` inside it: exactly 1 warning,
  naming `ps(t, ...)` and not `ps(u, ...)`, `identical()` to
  `predict()`'s.
- fit B with both past: 2 warnings, `identical()` to `predict()`'s.
- silent where it should be: `frm_lp_basis(fit)` in-sample 0 warnings,
  `frm_lp_basis(newdata = )` inside the span 0 warnings.

`ps_span_signal()` builds the condition with `warningCondition(..., class
= )`, whose `call` defaults to NULL, so the old `call. = FALSE` printing
is preserved. Verified: no `In ... :` prefix on either route.

### One correction to the lane's own account

The lane's `dev/spline-span-findings.md:73` says a tape build "evaluates
the R body ONCE" and calls the collector protection against a property
RTMB does not state. That is right about RTMB, but it understates the
case: the second entry I measured is the lane's OWN off-tape pass
(`R/predict.R:3284`). The collector is not defensive, it is required by
the design as written. The comment at `R/ps.R:263-266` and
`R/predict.R:3227-3232` both explain the collector as being about RTMB's
re-entry, which is the hypothetical reason rather than the actual one. A
reader debugging a double-fire later will look in the wrong place.

## (d) Extension

Fit C, `y ~ lev + ps(t, k = 10, pad = 0.1)`, knot span
`[-0.09366693, 1.09810595]`.

Inside the span all three are silent, past it the counts are as claimed,
`suppressWarnings()` still suppresses both curve functions (0 warnings),
and the refusal names the term and the range:

```
frm_curve_feature(): the search bracket leaves a ps() term's knot span,
... Narrow newdata to the span. 32 of 50 predicted values of t lie
outside the frozen knot span of ps(t, k = 10, pad = 0.1) [-0.09367, 1.098]. ...
```

The refusal is raised at `curve-feature.R:155`, after the grid scan and
before the `for (i in cross)` Newton loop at `:169`, so no root is
refined against the decaying partial sum. Confirmed by reading and by
the error arriving with zero warnings (the old 11 are gone).

`conditionCall()` is NULL on both core routes, so the old
`call. = FALSE` printing survives the move to `warningCondition()`.

## (e) rp_floored()

The lane left no script for this, so I wrote one from the
`sp_far_censored()` fixture (`test-rp-floored.R:18`). Every number
reproduces exactly:

| quantity | lane | mine |
| --- | --- | --- |
| `convergence` | 1 | 1 |
| `logLik` | -575.537942891 | -575.537942891 |
| `n_censored_deep` | 1 | 1 |
| `max_nlogS` | 55.7302136 | 55.7302136 |
| the deep row's exact term, `-H` | -55.73021360299 | -55.73021360299 |
| `S = exp(-H)` | 6.26146e-25 | 6.261462e-25 |
| `n_nonmonotone` | 0 | 0 |
| `rp_floored(action = "error")` | did not refuse | did not refuse |

`fitted()` and `predict(type = "response")` both still refuse, with the
message unchanged: "royston_parmar: a survival time has no mean on the
response scale here. ...".

The rename IS a user-detectable behavior change: `rp_floored()` is
exported, and on the returned list `r[["n_censored_floored"]]` is now
NULL. `$` partial matching does not rescue it, because "floored" is not
a prefix of "deep", so old code reading that field gets NULL silently
and fails one line later. `extensions/frmtmb.spline/NEWS.md` carries it
under an explicit "BEHAVIOR CHANGE:" heading with the count, the
`-log S` and the old floor value, which is the right treatment. No
deprecation shim, which is defensible at 0.2.0 but is a choice, not an
oversight to leave unremarked.

## (f) Docs

All four claims verified by running, not by reading.

- `lccdf`: present on `royston_parmar(df = 3)`, a function.
- `post$fit_check`: `names(fam$post)` is `mean_fn, fit_check`, and core
  reads it at `R/fit-end.R:32`.
- `frm_lp_basis()`: exported and documented; no `fit$cache$Vjoint` reach
  survives in `extensions/frmtmb.spline/R/` outside comments that
  describe the reach as historical.
- `frm_curve()` on an rp fit: draws, returns `frmtmb_curve`.

The vignette's `{r floored-bad}` chunk lost its `error = TRUE`, which is
correct now that `rp_floored()` does not refuse; run verbatim it prints
`n_censored_deep: int 1` and `max_nlogS: num 55.7`, as the vignette
text claims.

Vignette execution: pandoc is not installed here, so
`rmarkdown::render()` cannot finish; `knitr::knit()` with
`opts_chunk$set(error = FALSE)` runs every chunk, which is the
substantive check, and all of them ran. The only `Error:` in the output
is the deliberate `{r refuse, error = TRUE}` chunk at line 198
demonstrating the `fitted()` refusal. I did NOT verify the pandoc stage.

## (g) Tests

Fresh process per file, from the worktree sources against the rvss-lib
installs.

| suite | pass | fail | warn | skip | error |
| --- | --- | --- | --- | --- | --- |
| core `test-ps` | 71 | 0 | 0 | 0 | 0 |
| core `test-predict-lp-basis` | 30 | 0 | 0 | 0 | 0 |
| core `test-predict-newdata` | 12 | 0 | 0 | 0 | 0 |
| core `test-message-uniqueness` | 6 | 0 | 0 | 0 | 0 |
| `frmtmb.spline` full suite | 302 | 0 | 0 | 0 | 0 |

The spline total matches the lane's claim of 302. Per file the suite
reports curve 62, deriv 38, gratia 12, message-uniqueness 4,
royston-parmar 67, rp-floored 42, span 29, surface 48.

The `test-message-uniqueness.R` extension to pool `errorCondition()` and
`warningCondition()` with the plain forms is the right change: without
it the new `warningCondition()` template would have collected as the
empty string and dropped out of the pool unasserted. It skips `class =`,
`call =` and friends by NAME, which is correct, and descends one
`paste0()` deeper for the constructor's single message argument.

### Assertion quality

Mostly strong: `test-span.R` asserts counts (`expect_length(ws, 1L)`),
exact fragments (`"15 of 15"`, `"45 of 45"`, `"three rows per grid
point"`) and the span numbers themselves. Two weak spots, both noted as
findings below: one bare `expect_error()` with no regexp, and the
systematic use of `simultaneous = FALSE` on every past-span call.

## (h) Cross-lane: the predict.R merge with wt-car-jacobian

Every hunk this lane makes to `R/predict.R` is inside `lp_basis_nl()`.
New-file line ranges:

| hunk | new lines | what |
| --- | --- | --- |
| 1 | 3226-3243 | `has_ps` scan and the `span` collector |
| 2 | 3257 | `ps_env(lpk, ..., check = span)` in `ev_one()` |
| 3 | 3277 | `ps_env(lp, ..., check = span)` in the top-level body |
| 4 | 3280-3288 | the off-tape pass |
| 5 | 3292 | `ps_span_flush(span)` |

wt-car-jacobian owns `rr_jacobians()` at `R/predict.R:637` and
`lp_delta_A()` at `R/predict.R:1485`. There is NO textual overlap: no
added or removed line in this lane mentions either name, and the nearest
hunk is ~1740 lines away. `R/ps.R` hunks are 260, 286, 553-569, 572-575
and 580-628, all in `ps_env()` / `ps_span_*()`, which wt-car-jacobian
does not own.

The one coupling to watch is semantic rather than textual:
`lp_basis_nl()` CALLS `lp_delta_A()` (in its `build()` closure, around
`:3190`) and `rr_jacobians()`. If wt-car-jacobian changes that
signature, the call site is in neither lane's hunks, so git will merge
cleanly and leave the call site as whichever side wrote it. Worth one
build after the merge rather than a conflict resolution. The collector's
gate (`has_ps`, `is.null(newdata)`) reads nothing either function
produces, so it is independent of their result shapes.

---

# Findings, ranked

## 1. MEDIUM-HIGH. frm_curve_deriv() past the span warns and then dies on default arguments, and every test dodges it

`frm_curve_deriv()` defaults to `simultaneous = TRUE`
(`R/curve-deriv.R:81`). Past the knot span the basis is exactly zero, so
the derivative's standard error is exactly zero, and the simultaneous
critical value divides by it:

```r
# extensions/frmtmb.spline/R/curve-cov.R:167
mx <- apply(abs((L %*% z) / div), 2L, max)   # div == 0 -> NaN
crit <- unname(stats::quantile(mx, level, type = 8))   # :168, errors
```

Reproduction (fit C, `y ~ lev + ps(t, k = 10, pad = 0.1)`, span
`[-0.0937, 1.0981]`):

```r
nout <- data.frame(t = seq(0.5, 1.0981 + 1, length.out = 25))
frm_curve_deriv(fitC, "t", newdata = nout)
#> Warning: frm_curve_deriv(): this grid leaves a ps() term knot span, ...
#> Error in quantile.default(mx, level, type = 8) :
#>   missing values and NaNs not allowed if na.rm is FALSE
frm_curve_deriv(fitC, "t", newdata = nout, simultaneous = FALSE)   # ok
```

`sp_sim_crit()` is untouched by this lane, so the crash is PRE-EXISTING
and the lane did not cause it. What the lane owns is that it shipped
around it:

* `man/frm_curve_deriv.Rd` gains a "Past a ps() knot span" section
  promising "Warned once per `ps()` term per call". On default arguments
  the user gets the warning and then an error from `quantile.default`,
  which names neither the span nor the function.
* `tests/testthat/test-span.R` passes `simultaneous = FALSE` on EVERY
  past-span call (the `frm_curve()`, `frm_curve_deriv()` and no-ps
  control tests alike). The default path is never exercised, so the
  suite is green on a surface that does not work.

`frm_curve()` survives the same grid because its own standard error
keeps the intercept contribution and never reaches zero.

## 2. MEDIUM. Gridding exactly to the knot span, which is what the message tells the user to do, trips a spurious warning and a spurious REFUSAL

The span is checked on the difference stencil, not on the user grid.
`frm_curve_deriv()` builds `c(x - e, x, x + e)` (`R/curve-deriv.R:103`)
and the stationary-point `gfun` of `frm_curve_feature()` evaluates at
`c(tv - e1, tv + e1)` (`R/curve-feature.R:127`), with
`e = diff(range(x)) * 1e-6` for order 1 and `* 1e-4` for order 2
(`R/curve-deriv.R:194`). A grid whose endpoint sits on the knot is
therefore outside by `e`, and the refusal at `R/curve-feature.R:155`
fires on a bracket that is entirely inside the span.

The refusal message says "Narrow newdata to the span." and prints the
span. Following that literally is the failing input:

```
knot span [-0.09366692726, 1.09810595175]
grid = exactly the knot span     feature(max) = REFUSED   deriv = 1 warning
grid pulled in by 1e-07 * span   feature(max) = REFUSED   deriv = 1 warning
grid pulled in by 1e-06 * span   feature(max) = ok        deriv = 0 warnings
grid = the DATA range (usual)    feature(max) = ok        deriv = 0 warnings
```

The common case is safe, because `pad =` puts the knots outside the data
range. The user who acts on the message is the one who is bitten, and
for `frm_curve_feature()` the result is a hard error. `type =
"crossing"` is unaffected, because its scan does not widen the grid.

## 3. MEDIUM. The extension frmtmb floor is not raised, and the degradation against a released 0.52.0 is silent

`extensions/frmtmb.spline/DESCRIPTION:34` is still
`frmtmb (>= 0.52.0)`; core `DESCRIPTION:3` is still `Version: 0.52.0`.
The `frmtmb_ps_span_warning` class is new in the UNRELEASED development
state of 0.52.0, so the declared floor cannot express the dependency.

Against a released 0.52.0 the whole seam degrades with no error:
`sp_catch_span()` matches nothing, `frm_curve()` and
`frm_curve_deriv()` go silent again, and `frm_curve_feature()` reverts
to the eleven raw warnings with no refusal. Nothing tells the user. Core
needs a version bump and the extension the matching floor at
consolidation. The lane report flags the floor as consolidation work,
but neither NEWS file states the requirement.

## 4. LOW. The deriv warning loop shadows the grid row count

`R/curve-deriv.R:123` writes `for (m in parts$span)`, reusing `m`, which
is `nrow(nd)` from `:109` and is what `lo`, `mid` and `hi` are built
from at `:110-112`. Nothing reads `m` after the loop today, so it is
harmless, but the next edit that adds a line below the loop inherits a
character vector where the row count should be. Rename the loop
variable.

## 5. LOW. The collector last-write-wins can undercount, and the test bakes that in

`R/ps.R:583` keys on the term label and overwrites:
`assign(pt[["label"]], list(out, length(x), pt), envir = sink)`.
`tests/testthat/test-ps.R` asserts the semantics directly: collect
`c(5, 6)`, then `c(7, 8, 9)`, and the flush reports "3 of 3". The first
pass count is discarded rather than merged.

Today no user-visible undercount is reachable: I measured the two passes
and they always see the same `x`. But the documented invariant is "one
warning per term", and the implementation delivers "the last pass count
for that term", which are different statements. Two `ps()` terms sharing
a label would also collide on the key. Worth a comment that says which
one is guaranteed, or a merge instead of an overwrite.

## 6. LOW. The rationale comments name a cause that is not the operative one

`R/ps.R:263-266` and `R/predict.R:3227-3232` both explain the collector
as protection against RTMB re-entering the closure an unknown number of
times. Measured, RTMB enters once; the second entry is the lane OWN
off-tape pass at `R/predict.R:3284`. The collector is required by this
design rather than defensive against RTMB, and the comments send a
future reader to the wrong place. The lane
`dev/spline-span-findings.md:73` states the RTMB fact correctly and then
draws the same wrong conclusion.

## 7. LOW. expect_error() with no regexp

In `tests/testthat/test-span.R`, the test "the refusal replaces a
warning per Newton step, not adds to it" wraps `expect_error(...)`
around the `withCallingHandlers()` call with no `regexp`. It would pass
on any error at all, including one raised before the span check. The
`expect_length(ws, 0L)` beside it carries the real assertion; the error
assertion should name the refusal, as the neighboring test does.

## 8. INFO, no action. Things I checked that are correct

* Text identity between the two routes is exact, class is on both, one
  warning per term, the right term named when only one is past the span,
  under both `re.form = NA` and `re.form = NULL`.
* `conditionCall()` is NULL on both routes, so `warningCondition()`
  preserved the old `call. = FALSE` printing.
* The advector guard fix is real and the case is reachable (measured:
  `x=advector` on the tape pass of fit A at `re.form = NULL`).
* The off-tape pass is on exactly the newdata grid and costs about
  4 ms flat: 11.7% of the call at N=200, 0.1% at N=5000.
* All 600-subject rp numbers reproduce to every printed digit;
  `fitted()` and `predict(type = "response")` refuse unchanged.
* No NAMESPACE change, no `dev/*.log`, no file owned by another lane
  (`families.R`, `compat.R`, the RL vignette, `lp_delta_A()` and
  `rr_jacobians()` are all untouched).

---

# Verdict: PUNCH

The core work is sound and I could not break it: the collector fires
once per term by construction, the text is byte-identical to the one
`predict()` raises, the class is on both routes, the tape guard fix is
real and reachable, and the off-tape pass is cheap and on the right
grid. Findings 5 and 6 are comments and a latent invariant, not defects.

It is a PUNCH on the extension side, for three things in scope for
exactly what this lane was asked to deliver.

## Punch list

1. **frm_curve_deriv() on default arguments.** Past the span it warns
   and then errors out of `quantile.default`. Either guard the zero
   divisor in `sp_sim_crit()` (`R/curve-cov.R:167`), or have
   `frm_curve_deriv()` refuse the way `frm_curve_feature()` does, or at
   minimum say in `man/frm_curve_deriv.Rd` that the simultaneous band is
   not available past the span. Whichever is chosen, ADD A TEST ON THE
   DEFAULT ARGUMENTS: every past-span call in `test-span.R` currently
   sets `simultaneous = FALSE`, so the default path is untested.

2. **The stencil boundary.** Check the span on the user grid rather than
   on the widened stencil, or give the check a tolerance of the stencil
   width, so that a grid ending exactly on a knot is neither warned
   about nor refused. As it stands, following the advice the refusal
   itself gives reproduces the refusal. Add the exact-boundary grid as a
   test for all three functions.

3. **The version floor.** Bump the core `Version:` and raise
   `extensions/frmtmb.spline/DESCRIPTION` to match, so that
   `sp_catch_span()` cannot silently match nothing against a released
   0.52.0. State the requirement in the extension NEWS.

Optional, cheap, and I would take them in the same pass: rename the
`for (m in ...)` loop variable (`R/curve-deriv.R:123`); give the bare
`expect_error()` in `test-span.R` a regexp; and correct the two comments
that attribute the second closure entry to RTMB rather than to the lane
own off-tape pass.

## Counts I measured

* 23 paths, 21 modified plus 2 untracked, matching the claim.
* core `test-ps` 71 pass, `test-predict-lp-basis` 30, `test-predict-newdata`
  12, `test-message-uniqueness` 6; spline suite 302 pass. Zero failures,
  zero errors, zero skips anywhere.
* Closure entries per call: 2 per term on the basis route, 1 per term on
  the predict route. Warnings raised: 1 per term on both.
* Off-tape pass: about 4 ms, 11.7% of the call at N=200 and 0.1% at
  N=5000.
* rp 600-subject: logLik -575.537942891, `n_censored_deep` 1, `max_nlogS`
  55.7302136, deep row term -55.73021360299, `n_nonmonotone` 0.

---

# Re-check, 2026-09-06

Second pass over the same worktree after the lane addressed the punch
list. Still 23 paths, no new files, both DESCRIPTIONs and both
NAMESPACEs still untouched. Growth since the first pass: `curve-cov.R`
54 to 109 added lines, `curve-deriv.R` 15 to 33, `curve-feature.R` 31 to
48, `R/ps.R` 99 to 103, `R/predict.R` 32 to 34, spline NEWS 43 to 73.
Both packages reinstalled into `rvss-lib` and re-verified to load from
it. Everything below I ran; nothing is taken from the lane's report.

## Punch item 1: the zero divisor. FIXED, and the argument is exact

The guard is at `extensions/frmtmb.spline/R/curve-cov.R:207`
(`keep <- is.finite(div) & div > 0`), the all-degenerate refusal at
`:209`. Line numbers as claimed.

**Bit-identity.** I rebuilt the pre-guard body from the shipped function
by deleting the refusal block and undoing the row subsetting, installed
it with `assignInNamespace()`, and compared at a fixed seed on a grid
with no zero-se row:

```
guarded   curve=2.9022059772043827 deriv=2.8832279842129589
unguarded curve=2.9022059772043827 deriv=2.8832279842129589
identical(): TRUE   bitwise equal: TRUE
```

Identical, not merely close. The reason is structural and worth
recording: `z <- matrix(rnorm(m * nsim), m, nsim)` is still drawn at
full `m` before any subsetting, so the guard cannot move the RNG stream.
Had the lane subset `S` before the draw, this would have failed.

**Can a zero-se row be the argmax?** No, and by more than the PSD
argument. On the real past-span derivative grid (25 rows, 8 with se
exactly zero):

```
rows with se exactly 0 : 8 of 25
max |L[zero rows, ]|   : 0        <- exactly zero, not 1e-18
max |L[kept rows, ]|   : 3.304
max |numerator| at zero rows : 0
including them gives   : NaN
```

The zero rows of `L` are exactly zero in floating point, so the
deviation there is exactly deterministic; the standardized value is
`0 / 0`, which is the `NaN` that used to reach `quantile.default()`.
The lane's PSD argument is sound and, here, survives finite precision.

**On the defaults.** All three functions now work past the span with
`simultaneous = TRUE`:

```
frm_curve(default)          warns=1  ok
frm_curve_deriv(default)    warns=1  ok
frm_curve_deriv(order = 2)  warns=1  ok
frm_curve_feature(max)      warns=0  REFUSED by name
  band returned: 25 rows, 8 zero-se rows, .crit_sim 2.705058, all CI finite
```

An all-past grid refuses by name with a message that says what to do
("Use simultaneous = FALSE, or move the grid inside the span"). The two
new tests (`test-span.R:96` and `:126`) exercise the default path and
assert the band is finite, that the zero-se rows collapse to the
estimate, and that `.crit_sim > .crit` so it is still a simultaneous
band. That is a stronger assertion set than I asked for.

## Punch item 2: the stencil boundary. FIXED

`sp_span_on_grid()` at `curve-cov.R:92`, consulted at
`curve-deriv.R:137`, `curve-feature.R:166` and `:229`.

A grid laid exactly on the knot span is now silent everywhere, my
original reproduction included:

```
frm_curve(edge)                warns=0  ok
frm_curve_deriv(edge, o=1)     warns=0  ok
frm_curve_deriv(edge, o=2)     warns=0  ok
frm_curve_feature(edge, max)   warns=0  ok
frm_curve_feature(edge, min)   warns=0  ok
frm_curve_feature(edge, cross) warns=0  ok
```

And the gate does not let a real excursion through. A grid pushed past
the upper knot by 0.1, 1 and 5 percent of the span still warns
(derivative) and still refuses (feature) at every one of the three.

The derivative warning now counts grid rows: "16 of 25" on a 25-row
grid, where the stencil reading would have said 48 of 75.

**Cost.** I traced `sp_span_on_grid()`. It fires exactly once, and only
when the widened call already reported something:

```
deriv, inside span   grid re-checks=0  warnings=0
deriv, past span     grid re-checks=1  warnings=1
feature, inside      grid re-checks=0  warnings=0
feature, past        grid re-checks=1  warnings=0
```

A clean grid pays nothing, and no path double-warns.

## The monotonicity caveat, tested

The comment at `curve-cov.R:80-84` says the feature-scan gate holds
"whenever the spline argument is monotone in `var`". I fitted the case
it does not cover, `y ~ lev + ps(abs(t - 0.5), k = 8, pad = 0.1)`, whose
argument falls to zero at `t = 0.5` and rises again, knot span
`[-0.04930, 0.54966]` in the argument:

```
feature(max), all past       REFUSED
feature(cross), all past     REFUSED
deriv, all past              warns=1
feature(max), straddling (7 of 21 grid points past hi)   REFUSED
deriv, straddling            warns=1
feature(max), fully inside, straddling the kink          ok, silent
deriv, fully inside          ok, silent
feature(max), edge on the non-monotone argument          ok, silent
```

No genuine past-span point is missed. I also swept 401 grids whose upper
end walks through the knot and counted the cases where the grid is
outside the span but the scan stencil is not: **zero**.

The reason is stronger than monotonicity. The only miss would need the
argument to have a local maximum just outside the span with
`arg(x) > hi` while `arg(x +/- e1) < hi`, which for
`e1 = 1e-6 * range(x)` demands a curvature around `1e12`. Continuity
plus a stencil a millionth of the range wide is what actually carries
it. The comment's condition is therefore sufficient but stronger than
necessary, and it errs on the safe side, so I am not asking for a
change. Worth knowing that the guarantee does not in fact fall over on a
non-monotone argument.

## Punch item 3: the version floor. Done as instructed

Both DESCRIPTIONs untouched, per the coordinator's instruction. The
requirement is stated in `extensions/frmtmb.spline/NEWS.md:12-18`, and
it names the right file: the dependency is under `Depends:`, not
`Imports:`, and the NEWS text says `Depends:`. (My first-pass note
called it Imports; the line number was right and the section was not.)

## The smaller items

* Loop variable renamed to `msg` at `curve-deriv.R:141` and
  `curve.R:171`. Verified on disk.
* The bare `expect_error()` now carries
  `"the search bracket leaves", fixed = TRUE`, at `test-span.R:216-223`,
  with a comment saying why a bare one was wrong.
* Both rationale comments corrected: `R/ps.R:270-271` and the
  `R/predict.R` span block now say the closure is entered twice, once
  for the off-tape pass and once for the tape build, and both add that
  RTMB enters once so a reader chasing a double-fire looks at the
  off-tape pass. That is the correction I asked for.
* Finding 5 (collector last-write-wins) left deliberately. I accept the
  disposition: I measured both passes seeing the same `x`, and
  identical `ps()` calls in one body share a `deparse1()` label, so
  they are one term. The invariant is narrower than "one per term" but
  no user-visible undercount is reachable.
* `man/` is in sync with the roxygen sources: both new `@section`
  blocks are present in `frm_curve_deriv.Rd` and
  `frm_curve_feature.Rd`.

## One new NOTE, not a blocker

The second refusal in `frm_curve_feature()` (`curve-feature.R:229`)
calls `sp_span_on_grid(sp$fit, nd, ...)`, the same arguments as the
first at `:166`. Newton is clamped inside the grid, so the located roots
are inside `nd`; if the first check did not refuse then `nd` is clean,
and the second call asks the identical question and gets the identical
empty answer. It can never refuse when the first did not, so as a
refusal it is dead code, and on the boundary case it costs one extra
`predict()`. Harmless and defensible as belt-and-braces, but if a later
reader wonders whether it does anything, the answer is no. Either check
the ROOT stencil there (which is what its comment describes) or drop it.

## Counts, re-measured

| suite | pass | fail | warn | skip | error |
| --- | --- | --- | --- | --- | --- |
| core `test-ps` | 71 | 0 | 0 | 0 | 0 |
| core `test-predict-lp-basis` | 30 | 0 | 0 | 0 | 0 |
| core `test-predict-newdata` | 12 | 0 | 0 | 0 | 0 |
| core `test-message-uniqueness` | 6 | 0 | 0 | 0 | 0 |
| `frmtmb.spline` full suite | 321 | 0 | 0 | 0 | 0 |

Spline by file: curve 62, deriv 38, gratia 12, message-uniqueness 4,
royston-parmar 67, rp-floored 42, span 48, surface 48. Sums to 321.
`span` was 29 on the first pass, so the lane added 19 assertions, and
the claim of 321 is exact.

---

# Final verdict: MERGE

All three punch items are fixed, and I verified each by running it
rather than by reading the diff. The two I could have been wrong about,
I attacked directly: the guard is bit-identical to the unguarded code on
a clean grid at a fixed seed (`identical()` TRUE, and the RNG stream is
structurally protected because `z` is drawn before the subsetting), and
the dropped rows have `L` exactly zero so they cannot be the argmax in
floating point any more than in theory. The stencil-boundary fix does
not open a hole: a grid 0.1 percent past the knot is still caught, and a
401-grid sweep plus a deliberately non-monotone `ps(abs(t - 0.5))` fit
found no missed past-span point.

The new tests are better than the ones they replace: they exercise the
default arguments that the first round dodged, and they assert the
band's content rather than only its absence of error.

Remaining, and none of it blocking: the dead second refusal in
`frm_curve_feature()` noted above, the collector invariant left as-is by
agreement, and the version floor, which is consolidation work by
instruction and is written down in NEWS where consolidation will find
it.
