# Review of lane splinecurve

Round of 2026-09-29, reviewed 2026-09-29. The change lives in the
release integration tree `C:\Users\adf44\source\r\frmtmb-wt-release`,
under `extensions/frmtmb.spline/**`. Arms:

- `after`: the worker's build, `C:/Users/adf44/source/r/wt-splinecurve-lib`.
- `before`: the release library, `C:/Users/adf44/source/r/rellib-r3`,
  frmtmb.spline 0.8.1 without this change.

Every reviewer script is `dev/splinecurve-rev-*.R`; logs are in
`dev/splinecurve-rev-log/`. Nothing else in the tree was changed. No git
operation was run beyond read-only `git diff`, `git status` and
`git show`. No library was installed into.

## Verdict: NOT MERGEABLE

Blockers, most severe first:

1. **The flag opens three wrong-number routes at an unseen level of a
   bar term, and the lane guards only the fourth (the contrast).** With
   `allow_new_levels = TRUE` and a `re_formula` that keeps `(1 | g)` or
   `(1 + t | g)`, `frm_curve()`'s simultaneous band has a critical value
   BELOW 1.96 and covers 35 percent (intercept) or 15 percent (slope)
   instead of 95; `frm_curve_deriv()`'s standard error is 0.40 to 0.50
   of the right one with a random slope; `frm_curve_feature()`'s
   `.value_se` and crossing `.se` are about a quarter of the right ones.
   The new `?frm_curve` section invites exactly this call. Before the
   lane every one of these calls refused. Details under "Finding A".
2. **The stored-curve OR rule overrides an explicit
   `allow_new_levels = FALSE`.** Claim 2 below. It combines with
   blocker 1: a stored curve, a new grid at an unseen level, and the
   default flag return a `.value_se` of 0.0146 where `frm_curve()` at
   the same point says 0.0432.

Not blockers, to fix in the same pass:

3. Record defect: `dev/splinecurve-findings.md` says "611 passed". Its
   own table sums to 531, and so does the rerun.
4. Two stale doc passages, confirmed (claim 6).

## Integrity of the arms

`dev/splinecurve-rev-00-integrity.R`, logs
`00-integrity-wt-splinecurve-lib.txt`, `00-integrity-rellib-r3.txt`.

- frmtmb 0.65.0 is the same build in both libraries: md5 over every
  deparsed namespace function `443c188d7b1120d082e18270078198f2` in both.
- `after`: all 55 frmtmb.spline functions deparse identically to the
  tree's current `R/` source.
- `before`: 12 of 55 differ, and they are exactly the functions the diff
  touches (`frm_curve`, `frm_curve_deriv`, `frm_curve_feature`,
  `sp_assemble`, `sp_cov_check`, `sp_curve_parts`, `sp_grid_span`,
  `sp_linkinv`, `sp_one_basis`, `sp_predict_eta`, `sp_span_on_grid`,
  `sp_spec`). `before` has no `allow_new_levels` formal.

## Claim 1: the band at an unseen fs level. HOLDS for fs

See Finding A for an unseen bar-term level.

`dev/splinecurve-rev-01-band.R`, log `01-band.txt`. Data: the
vignette's simulation, `set.seed(4)`, 7200 rows, 20 subjects.

A. `v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)`, 80-point grid
at the unseen level `"population"`:

| Quantity | Reviewer | Worker |
|---|---|---|
| fs columns of `A` all exactly 0 | TRUE (100 of 112) | TRUE |
| `extra_var` all exactly 0 | TRUE | TRUE |
| `.se` `identical()` to hand restriction to 12 columns | TRUE | TRUE |
| max abs estimate diff to mgcv ML `exclude =` | 7.889e-07 | 7.9e-07 |
| se ratio frm / mgcv, range | 0.99999 to 1.00357 | 1.000 to 1.004 |
| simultaneous crit, `nsim = 20000, seed = 1` | 2.710666008 | 2.710666 |
| crit mcse | 0.0112596 | 0.01125959 |

The hand-restriction identity is an IDENTITY (exact zero columns), as
the worker says; `identical()` is its numerical check. The mgcv ratio is
a measurement.

`re_formula = NULL` against `NA` at the unseen fs level: the two
`frm_curve()` results are `identical()` once the stored `re_formula` is
aligned. Median `.se` 0.0277 at the unseen level, 0.0073 at subject 3.

B. `s(t) + s(t, subject, bs = "fs") + (1 | subject)` at the unseen
subject. `extra_var` is exactly `sd(subject)^2` at `re_formula = NULL`
(relative error 0) and exactly 0 at `NA`, and
`.se^2 = diag(A V A') + extra_var` to 2.8e-16 relative. So at `NULL`
the band IS population smooth plus between-subject intercept variance.
On this design it does not matter: the fs term absorbs the per-subject
intercept and `sd(subject)` fits to 6.99e-05, so `extra_var` is 4.9e-09.
It is not documented beyond the sentence that a `(1 | g)` term's unseen
level adds its variance. Finding A applies to it whenever that sd is not
negligible.

C. `s(t) + s(subject, bs = "re")` and D.
`t2(t, subject, bs = c("cr", "re"), k = c(5, 20))` at the unseen level:
core refuses under `TRUE` with "allow_new_levels = TRUE cannot predict
the new level(s) population of the smooth term ... whose basis has an
`re` factor in it", at `NA` and `NULL`, in `frm_curve()`,
`frm_curve_deriv()` and `frm_curve_feature()`. `FALSE` refuses with
class `frmtmb_new_levels`. A SEEN level under `TRUE` answers. No zero
row, no wrong answer.

E. `v ~ f + s(t, subject, bs = "fs", by = f, k = 5)`: every combination
of seen/unseen, `TRUE`/`FALSE`, `NA`/`NULL` is refused. Under `TRUE`,
and under `FALSE` at a seen level, it is core's "predict(newdata = ) is
not supported for the smooth s(t,subject):fb" refusal; `FALSE` at an
unseen level hits the new-level refusal first. `frm_linpred()` refuses
the same way. The flag opens no path around the refusal.

## Finding A: an unseen bar-term level. BLOCKER.

`dev/splinecurve-rev-02-extravar.R`, log `02-extravar.txt`. Same data.
60-point grid at an unseen subject, `re_formula = NULL`,
`allow_new_levels = TRUE`. `extra_var` equals `diag(Z S Z')` exactly,
with `Z` the new level's design rows and `S` the block covariance from
`VarCorr()`, so the full covariance of the grid is
`C V C' + Z S Z'`: ONE draw shared by every row. frm_lp_basis() returns
only its diagonal.

| | `(1 \| subject)` | `(1 + t \| subject)` |
|---|---|---|
| pointwise `.se` vs `sqrt(diag(full))` | 0 rel. | 0 rel. |
| share of `.se^2` that is `extra_var`, median | 0.940 | 0.947 |
| `.crit_sim` reported (`nsim = 20000, seed = 1`) | **0.646** | **0.615** |
| crit from the full covariance | 2.146 | 2.472 |
| coverage of the reported band under the full covariance (20000 draws, seed 99) | **0.352** | **0.148** |
| `frm_curve_deriv()` `.se` / right se | 1 | **0.40 to 0.50** |
| feature `.value_se` at the two roots of `at = 0.5` | 0.0091, 0.0092 | 0.0112, 0.0116 |
| `frm_curve()` `.se` at the same `t` | 0.0376, 0.0376 | 0.0488, 0.0502 |
| crossing `.se` reported | 0.0025, 0.0025 | 0.0031, 0.0032 |
| crossing `.se` with the new level's variance | 0.0102, 0.0104 | 0.0132, 0.0139 |

The causes, all in code the lane routes the flag through:

- `sp_sim_crit(Sigma, se, ...)` draws from `Sigma = C V C'`, which has
  no `extra_var`, and divides by `se`, which has it. The standardized
  process has variance about 0.06, so the maximum sits below 1.96.
- `frm_curve_deriv()` forms `D V D'` and drops `extra_var` entirely.
  Right for a random intercept (its derivative is 0), wrong for a slope.
- `frm_curve_feature()` forms `.value_se` and the crossing `.se` from
  `qf(C0) = sqrt(diag(C0 V C0'))`, again without `extra_var`, so a
  feature and a band on one fit now DISAGREE, which the file's own
  comment says cannot happen.
- The `"Sigma"` attribute, documented as "the grid covariance", has
  `diag(Sigma)` 94 percent short of `.se^2`.

Reachability, by construction (`dev/splinecurve-rev-03-or.R`,
`dev/splinecurve-rev-07-na-level.R`, both arms): on `before` the same
calls refuse with `frmtmb_new_levels`, at an unseen label, at an `NA`
label, and with the grouping column absent. So these wrong numbers are
NEW in this lane. Under `TRUE` core also fills an absent `(x | g)`
column with `NA`, so a grid with no `subject` column reaches the same
route on a bar-term model.

The new `?frm_curve` section ends: "An unseen level of a `(1 | g)` term
is different: it adds its block's marginal variance, so a curve read
there at `re_formula = NULL` is a new subject's curve". It presents the
case as supported.

Pre-existing and NOT this lane's: an exact `gp()` off the observed
positions reaches the same three code paths with its kriging variance.
The recorded magnitude there is 7e-07 (`R/zzz.R`), so the effect was
negligible; the same fix covers it. For the backlog.

What would clear it: extend the lane's own refusal. When
`allow_new_levels = TRUE` and any row carries nonzero `extra_var`,
refuse in `frm_curve_deriv()`, `frm_curve_feature()`, and
`frm_curve(simultaneous = TRUE)`; the pointwise `frm_curve()` answer is
right and may stay. Rewrite the `(1 | g)` sentence to say what is
refused and why. Add a test per route that FAILS on the current build
with the numbers above (a crit below `qnorm(0.975)` is a run-measured
ratio, not an absolute tolerance).

## Claim 2: the stored-curve OR rule. FALSIFIED as a design.

`dev/splinecurve-rev-03-or.R`, logs `03-or-wt-splinecurve-lib.txt`,
`03-or-rellib-r3.txt`. Fixture: `test-new-levels.R`'s own
(`set.seed(4)`, 12 subjects, 1200 rows).

- A curve built with `TRUE` on a grid of SEEN levels gives numbers
  `identical()` to the default. The worker's "changes no number on a
  seen row" holds.
- `frm_curve_deriv(stored_TRUE, newdata = unseen)` with the default:
  ANSWERED, 30 rows. The same call on the FIT refuses.
- The same with an EXPLICIT `allow_new_levels = FALSE`: ANSWERED. The
  caller's explicit value is silently overridden. The `@param` text says
  so ("unless this call passes TRUE"), but a documented override of an
  explicit argument is still one.
- On the `(1 | subject)` fit, a curve stored with `TRUE` at a seen level
  and `re_formula = NULL`, then `frm_curve_feature()` on an unseen grid
  with the default: answered with `.value_se` 0.0146 and 0.0143 where
  `frm_curve()` at the same `t` gives 0.0432 and 0.0431 (Finding A).

Reusing the stored value when the caller says nothing is consistent with
how `re_formula` is reused (the caller's `re_formula = NULL` on a stored
`NA` curve is ignored, measured in the same log). OR is not: it cannot
express "refuse". The rule should be "the caller's value when given,
the stored value otherwise", with `missing(allow_new_levels)` passed
from the two exported functions to `sp_spec()`. That keeps the worker's
wanted case (stored `FALSE`, new unseen grid, caller `TRUE`) and honors
an explicit `FALSE`. A test needs the explicit-`FALSE` case.

## Claim 3: the new contrast refusal. HOLDS

Its false alarms are the documented ones.

`dev/splinecurve-rev-04-contrast.R`, log `04-contrast.txt`. Vignette
data plus a within-subject factor `cond` (`set.seed(11)`, +0.05 on
`"b"`). The unguarded arm is a copy of the lane source with lines
400 to 411 of `curve-cov.R` deleted, loaded with pkgload from
`%TEMP%\splinecurve-rev-noguard`.

| Case | Guarded | Unguarded |
|---|---|---|
| 1. two unseen `(1 \| g)` levels A minus B, `NULL` | refused | est 0, `.se` exactly 0 on 5 rows |
| 2. ONE unseen level, `cond` b minus a, `NULL` | refused | est 0.048891, se 0.0025604 |
| 2'. same at `re_formula = NA` | est 0.048891, se 0.0025604 | same |
| 3. two unseen fs levels, pointwise | est 0, se 0 | same |
| 3'. same, `simultaneous = TRUE` | refused: every se is 0 | same |
| 4. seen fs level 3 minus unseen fs level | answered | same |

- Case 1 reproduces the worker's defect: `.se` exactly 0 where
  `sqrt(2) * sd(subject)` is 0.0517 on this fit (the worker's 0.0515
  is the same quantity on the data without `cond`).
- Case 2 is the false alarm the worker named. The draw really is shared
  and the unguarded answer is right; the guarded call refuses, and the
  advice in the message (`re_formula = NA`) returns the same numbers to
  the last printed digit. Acceptable because it fails closed and says
  how to get the answer.
- Case 3: two unseen fs levels give a difference of exactly 0 with
  standard error 0. That is correct under the documented meaning of an
  unseen fs level (the population curve, twice), and it matches what
  `(1 | g)` at `re_formula = NA` gives, which the worker's test asserts.
  It is NOT the difference between two new subjects, and the docs do
  not say so. A reader who asks "how different are two new subjects"
  gets "not at all, with certainty". Worth one sentence.
- Case 4: `.se` equals `(A3 - Anew) V (A3 - Anew)'` by hand with
  relative error 0, and the estimate equals the difference of two
  `frm_linpred()` calls exactly.

The guard fires on a correct request in two cases: case 2 above and an
exact `gp()` contrast under `TRUE` with no unseen level. Both are in the
refusal list and the message.

## Claim 4: everything else unchanged. HOLDS.

`dev/splinecurve-rev-05-bitwise.R` run once in each library,
`dev/splinecurve-rev-05b-compare.R`, log `05-compare.txt`. 22 calls
modelled on the suite's fixtures: `s()` with simultaneous band, first
and second derivative, stored-curve derivative, maximum and crossing
features; Poisson `transform = TRUE`; `(1 | g)` at `NA`, at `NULL` on a
seen level, as a contrast of two seen levels, and at an unseen level
under the default (an error in both arms); a `sigma` dpar; a
`fac + s(x, by = fac)` difference with its derivative and crossing;
`gp()` curve and difference; a nonlinear `ps()` body inside and past the
knot span and its derivative; an fs model at a seen subject.

- 22 of 22 are `identical()`, value and warnings, once
  `attr(, "spec")$allow_new_levels` is removed. The `"fit"` attribute
  was removed in both arms before saving because it holds a TMB pointer.
- Positive control on the instrument: 16 of 22 are NOT identical before
  the spec field is removed (the 5 feature objects and the error carry
  no spec), so the comparison sees a one-field difference.

## Claim 5: tests and check. HOLDS, with a count error in the record.

One test file per `Rscript` process, `NOT_CRAN=true`, runner
`dev/splinecurve-run-tests.R` (read, not changed), logs in
`dev/splinecurve-rev-log/suite/`.

| File | Library | pass | fail | err | skip |
|---|---|---|---|---|---|
| test-new-levels.R | before | 6 | 0 | 6 | 0 |
| test-bracket-access.R | after | 1 | 0 | 0 | 0 |
| test-conditions-census.R | after | 11 | 0 | 0 | 0 |
| test-conditions.R | after | 5 | 0 | 0 | 0 |
| test-curve.R | after | 65 | 0 | 0 | 0 |
| test-deriv.R | after | 38 | 0 | 0 | 0 |
| test-difference.R | after | 59 | 0 | 0 | 0 |
| test-frailty.R | after | 30 | 0 | 0 | 0 |
| test-gratia.R | after | 13 | 0 | 0 | 0 |
| test-message-uniqueness.R | after | 4 | 0 | 0 | 0 |
| test-new-levels.R | after | 37 | 0 | 0 | 0 |
| test-royston-parmar.R | after | 67 | 0 | 0 | 0 |
| test-rp-floored.R | after | 81 | 0 | 0 | 0 |
| test-scale.R | after | 0 | 0 | 0 | 1 |
| test-span.R | after | 73 | 0 | 0 | 0 |
| test-surface.R | after | 47 | 0 | 0 | 0 |

15 of 15 files on `after`: **531** passed, 0 failed, 0 errors, 1 skip
(the scale tier, gated by `FRMTMB_SCALE_TESTS`). Every per-file count
matches the worker's table; the worker's total of 611 is an addition
error. The `before` run of `test-new-levels.R` is the weak form of
seen-to-fail ("unused argument"); the noguard run (claim 3) and the
refusals on `before` (claim 2 log) are the behavioral form.

`test-new-levels.R` has no absolute numeric tolerance: `expect_equal()`
defaults, `identical()`, exact `== 0`, a ratio to a run-measured gap
(`1e-8` of the seen-subject gap) and a bound of one grid step. It has no
test for Finding A or for an explicit `FALSE` on a stored curve.

`dev/splinecurve-check-after.log` (UTF-16; read after conversion):
`R CMD check --as-cran` on frmtmb.spline 0.8.1, `Status: 1 NOTE`. The
NOTE is "Skipping checking math rendering: package 'V8' unavailable",
which is environmental. Tests OK in 39 s, vignette re-building OK, PDF
manual OK. The `NativeCommandError` block in the log is PowerShell
wrapping R's stderr, not a check result.

## Claim 6: vignette and docs. HOLDS; two stale passages confirmed

`dev/splinecurve-rev-06-rd.R`, log `06-rd.txt`: the three Rd files
render with no dropped characters. The new section and the `@param`
read correctly for fs. The `(1 | g)` sentence is the problem of
Finding A.

Vignette: `dev/splinecurve-vignette-064.txt` and
`dev/splinecurve-vignette-after.txt` differ, after whitespace collapse,
only in the new `subject` column and where the table wraps. The 0.64.0
side was checked against `git show 9b4bb650:docs/frmtmb.spline/articles/
curve-inference.md`, which carries 2.710666, 0.01125959, 0.05442608
and 0.513161. Independently, `rev-01` reproduces crit 2.710666008 and
mcse 0.0112596 on `after`. One cosmetic nit: the reflowed paragraph
after the `peak` chunk leaves "of the height" alone on a source line; it
renders as one paragraph.

Stale passage 1, `?frm_curve` "What the covariance is, and how it is
checked": says frmtmb exports no route to the joint covariance, that the
function rebuilds the design by unit perturbation, and that agreement is
"at the tenth significant figure or better". The next section says the
opposite and is right. It should say: the design `A`, the joint
covariance `V` at `A`'s columns, and `extra_var` come from
`frmtmb::frm_lp_basis()`; `Sigma` is `A V A'`; the pointwise standard
error adds `extra_var`; every call compares it with
`frm_linpred(se.fit = TRUE)` and refuses past `tol`; agreement is at
machine precision (2.2e-16 on the vignette's model). The simplest fix
is to delete the section and keep "The route to the covariance".

Stale passage 2, `R/zzz.R` compat row `frm_curve` / `predict`: "the
design is the difference between predictions one coefficient apart".
It should say that the design and covariance come from
`frmtmb::frm_lp_basis()` in one call, and that the answer is checked
against `frm_linpred(se.fit = TRUE)`. The `gp` row's "adds it to the
diagonal where it belongs" is also incomplete: the simultaneous band
and the derivative and feature standard errors do not add it (Finding
A, pre-existing part).

## What was run

| Script | Library | Log |
|---|---|---|
| rev-00-integrity.R | both | 00-integrity-*.txt |
| rev-01-band.R | after | 01-band.txt |
| rev-02-extravar.R | after | 02-extravar.txt |
| rev-03-or.R | both | 03-or-*.txt |
| rev-04-contrast.R | after, and after with the guard deleted | 04-contrast.txt |
| rev-05-bitwise.R, rev-05b-compare.R | both | 05-*.txt |
| rev-06-rd.R | none | 06-rd.txt |
| rev-07-na-level.R | both | 07-na-*.txt |
| splinecurve-run-tests.R, 16 processes | see table | suite/ |

## MERGEABLE or not

NOT MERGEABLE. Blockers: Finding A (wrong simultaneous band, derivative
and feature standard errors at an unseen bar-term level, newly reachable
and documented as supported), then the OR rule (an explicit `FALSE` is
ignored). The fs population curve itself, which is the lane's purpose,
is right, and claims 1, 3, 4 and 6 hold as stated for it.

## Re-check after punch round 1, 2026-09-29

Build under review: the worker's rebuilt library, frmtmb.spline 0.9.0
(the consolidating session's version). `rev-00` on it: all 57 functions
deparse identically to the tree's `R/`, and core is the same md5 as
before. Logs are `dev/splinecurve-rev-log/r2-*`.

### Blocker 1 (Finding A). CLEARED.

`dev/splinecurve-rev-08-recheck.R`, log `r2-08-recheck.txt`. Same data
and grid as `rev-02`. The four routes that gave wrong numbers, plus
`transform = TRUE` and `order = 2`, on `(1 | subject)` and
`(1 + t | subject)`:

| Route | unseen, `NULL` | column absent, `NULL` | unseen, `NA` |
|---|---|---|---|
| `frm_curve(simultaneous = FALSE)` | answers | answers | answers |
| `frm_curve(simultaneous = TRUE)` | refuses | refuses | answers |
| same with `transform = TRUE` | refuses | refuses | answers |
| `frm_curve_deriv()`, order 1 and 2 | refuses | refuses | answers |
| `frm_curve_feature()`, maximum and crossing 0.5 | refuses | refuses | answers |
| crossing 0.5 on a rootless window, t in [0.9, 1] | 0 rows | 0 rows | 0 rows |

- The pointwise `.se` that still answers is right: relative error 0
  against `sqrt(diag(A V A') + extra_var)` and 2.2e-16 against
  `frm_linpred(se.fit = TRUE)`, on both models.
- False-alarm side: at `re_formula = NA` with `TRUE`, `extra_var` is
  exactly 0 and nothing refuses. The pointwise curve, simultaneous
  curve, first derivative and crossing are `identical()` in `.estimate`
  and `.se` to the same calls at a SEEN subject with the default.
- The vignette's fs population curve answers all four routes, with
  crit 2.710666008, mcse 0.01125959, peak 0.513161, se 0.003104161,
  unchanged.
- Rootless window: the refusal sits after the root search, so a window
  with no root returns 0 rows. It carries no number, so nothing wrong
  gets through. A window WITH roots at the unseen level refuses
  (crossing 0.5 above). The search itself reads only the estimate, which
  does not depend on the new level's variance. The stencil is row 1 of
  the grid repeated, so a grid whose row 1 is seen and later rows unseen
  answers for row 1 only, which is correct for row 1.

New false alarm, not a blocker: with `TRUE` and an exact `gp()`
off the observed positions, the kriging variance is nonzero
(8.2e-07 to 1.1e-06), so `frm_curve(simultaneous = TRUE)` and
`frm_curve_deriv()` now refuse with no unseen level anywhere. `FALSE`
answers. This is the rule as the new section states it, and it fails
closed. The message, however (`rev-10`, log `r2-10-gpmsg.txt`), says
"a new level's marginal variance" and advises `re_formula = NA` or a
seen level. Following that advice still refuses, because `NA` keeps
`gp()`. The advice should also say "or pass allow_new_levels = FALSE
when no row is at an unseen level", as the contrast refusal already
does.

### Blocker 2 (the OR rule). CLEARED.

`dev/splinecurve-rev-09-stored.R`, log `r2-09-stored.txt`, and `rev-03`
rerun, log `r2-03-or.txt`. Fixture of `test-new-levels.R`. Stored value
by caller's value, on a new grid at the unseen fs level:

| Stored | Caller | `frm_curve_deriv()` | `frm_curve_feature()` |
|---|---|---|---|
| TRUE | missing | answers | answers |
| TRUE | FALSE | refuses, `frmtmb_new_levels` | refuses |
| TRUE | TRUE | answers | answers |
| FALSE | missing | refuses | refuses |
| FALSE | FALSE | refuses | refuses |
| FALSE | TRUE | answers | answers |

Every answer from deriv is `identical()` to deriv on the fit with
`TRUE`. On the stored grid itself (seen levels), an explicit `FALSE`
answers. A curve stored at the unseen level refuses when re-read with
an explicit `FALSE` and answers when the argument is missing. The
`(1 | subject)` chain that returned `.value_se` 0.0146 in round 1 now
refuses.

### Tests, check, Rd

- Suite, one file per process on the rebuilt library, logs
  `dev/splinecurve-rev-log/suite2/`: all 15 files ran, 0 fail, 0 error,
  1 skip (`test-scale.R`, gated). `test-new-levels.R` passes 59. Every
  other count is unchanged from round 1. My own sum of the pass column
  is **553**, which matches the worker.
- `dev/splinecurve-check-r2.log`: frmtmb.spline 0.9.0,
  `Status: 1 NOTE`. The note is the environmental V8 one. Tests OK in
  34 s. The checked tarball `final-r2/frmtmb.spline_0.9.0.tar.gz`
  matches the current tree in `R/`, `man/`, `vignettes/` and `NEWS.md`.
  `test-new-levels.R` was edited 100 s after the tarball was built; the
  only difference is one `test_that()` title, shortened. It does not
  affect the result.
- Rd, rendered once (`rev-06`, log `r2-06-rd.txt`):
  - The stale covariance section is gone, merged into "The route to
    the covariance". That section now states `extra_var` and the
    2.2e-16 agreement.
  - The new section "An unseen level of a '(1 | g)' term" says what
    answers, what refuses and why, with the round-1 numbers. It also
    says the `"Sigma"` attribute does not carry `extra_var`.
  - The fs section now says two unseen fs levels contrast to 0 with se
    0, and that this is not the difference between two new subjects.
  - The fs plus `(1 | g)` case is documented.
  - No dropped characters.
- `R/zzz.R` compat rows, read in the diff:
  - The `frm_curve` / `predict` row now says the design and covariance
    come from one `frm_lp_basis()` call, checked against
    `frm_linpred(se.fit = TRUE)` at 2.2e-16.
  - The `gp` row now says the simultaneous band, the derivative and the
    feature standard errors do not carry the kriging variance, and gives
    its size. That is true under the default `FALSE`; under `TRUE`
    those routes now refuse.
  - Both are correct.
- Not reproduced: the worker's "12 failures on the round-1 build". That
  build no longer exists. `rev-02`'s round-1 log is the behavioral
  evidence that the new refusal assertions would have failed there.

### Verdict after re-check: MERGEABLE

Both blockers are cleared by construction on the rebuilt library. One
non-blocking fix is recommended before or after merge: the refusal text
for the `gp()`-only case should name `allow_new_levels = FALSE` as the
way out.
