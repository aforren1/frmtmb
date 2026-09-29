# Lane splinecurve: `allow_new_levels` in frmtmb.spline

Released as frmtmb.spline 0.9.0, not 0.8.1: the consolidating session
raised the bump to minor because `allow_new_levels` is a new argument on
three exported functions. "0.8.1" below is the version the tree carried
while the lane worked.

Lane of 2026-09-29, worked inside the release integration tree
`frmtmb-wt-release`, frmtmb 0.65.0 core with frmtmb.spline 0.8.1. Only
`extensions/frmtmb.spline/**` and `dev/splinecurve-*` were touched. No
git operation was run by the lane beyond one read-only `git diff --stat`.

## The defect

frmtmb 0.65.0 keeps every smooth under `re_formula = NA`, an
`s(t, subject, bs = "fs")` term included. The population curve of such a
model is read at an UNSEEN level with `allow_new_levels = TRUE`. The
three curve functions took no `allow_new_levels`, so:

- `R CMD build` of frmtmb.spline failed at
  `Quitting from curve-inference.Rmd:86-90 [curve]` with core's
  "needs the grouping column `subject`" refusal.
  Log: `dev/splinecurve-build-before.log` (release library `rellib-r3`).
- No call of `frm_curve()`, `frm_curve_deriv()` or `frm_curve_feature()`
  on 0.8.1 could return that population curve. At an unseen level each
  refused with `frmtmb_new_levels`, at `re_formula = NA` and `NULL`
  alike, while core's `frm_linpred(allow_new_levels = TRUE)` answered.
  Script `dev/splinecurve-before.R`, log `dev/splinecurve-before.log`.

## What changed

- `frm_curve()`, `frm_curve_deriv()`, `frm_curve_feature()`: new argument
  `allow_new_levels = FALSE`, after `re_formula`, validated with
  `sp_check_flag()`.
- It is passed to EVERY core call that reads a grid: `frm_lp_basis()` in
  `sp_one_basis()`, `frm_linpred()` in `sp_predict_eta()` (the feature
  root search and the `ps()` span re-check go through it),
  `frm_linpred(se.fit = TRUE)` in `sp_cov_check()`, and both
  `frm_linpred()` calls in `sp_linkinv()` for `transform = TRUE`. The
  internal helpers take it with NO default, so no route can fall back
  to `FALSE` silently.
- A `frmtmb_curve` stores it in its `"spec"` attribute.
  `frm_curve_deriv()` and `frm_curve_feature()` on a stored curve use
  the caller's value when the argument is SUPPLIED (`FALSE` included)
  and the stored value only when it is missing. Round 1 used an OR of
  the two, which overrode an explicit `FALSE`; the review falsified it
  (round 2, blocker 2 below).
- NEW REFUSAL, needed by the new argument: `frm_curve(contrast = ,
  allow_new_levels = TRUE)` refuses when either grid has nonzero
  `extra_var`. See "The difference guard" below.
- NEW REFUSALS, round 2: with `allow_new_levels = TRUE` and nonzero
  `extra_var` on any grid row, `frm_curve(simultaneous = TRUE)`,
  `frm_curve_deriv()` and `frm_curve_feature()` refuse, naming the
  grouping factor. The pointwise `frm_curve()` band is answered. See
  "Round 2" below.
- Roxygen: `@param allow_new_levels` on `frm_curve()` (inherited by the
  other two), a new section "The population curve of a factor-smooth
  model", a new item in the difference-curve refusal list, and a note on
  `object` in `frm_curve_deriv()`. Rendered with `tools::Rd2txt` and
  read: `dev/splinecurve-rd.R`, `dev/splinecurve-rd.log`.
- Vignette `curve-inference.Rmd`: see below.
- `NEWS.md`: three bullets added under the existing entry in round 1,
  two more in round 2. The heading and `DESCRIPTION` read 0.9.0, set by
  the consolidating session and kept. `frmtmb (>= 0.65.0)` unchanged.
- New test file `tests/testthat/test-new-levels.R`.

Exported functions that forward `re_formula` to core, found by grep of
`R/` for `re_formula`, `frm_linpred`, `frm_lp_basis`: only the three
above. NOT changed:

- `rp_floored()`: calls `frm_linpred(object, type = "link", dpar = p)`
  with no `newdata`, on the training rows, so every level is seen and
  the flag has nothing to do. It takes no `re_formula` either.
- `royston_parmar()`: a family constructor. It makes no prediction.

## The band at the unseen fs level

Script `dev/splinecurve-band.R`, log `dev/splinecurve-band.log`. Data:
the vignette's simulation, `set.seed(4)`, 7200 rows over 20 subjects.
Model `v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)`, gaussian.
Grid: 80 points on [0, 1] at `subject = "population"`, a level the fit
did not see. `frm_curve(..., nsim = 20000, seed = 1)`.

What core returns, `frm_lp_basis(re_formula = NA, allow_new_levels =
TRUE)`, and identically at `re_formula = NULL`:

- `A` has 112 columns: `beta` 1, `beta.s(t)` 1, `b.s(t)` 10,
  `b.s(t,subject)` 100.
- `extra_var`: exactly 0 on every row (`all(extra_var == 0)` TRUE).
- All 100 `fs` columns of `A` are exactly 0.

So core adds NO between-subject variance for an unseen `fs` level. The
reason is structural: `pred_design()` skips `smooth` blocks when it
builds `re_parts`, and the new-level variance is formed only from
`re_parts`. The band is therefore the uncertainty of the intercept and
the `s(t)` coefficients alone, which is the band the population curve
needs. The vignette does not have to qualify it.

Comparisons:

| Quantity | Result |
|---|---|
| `.estimate` vs `frm_linpred(allow_new_levels = TRUE)` | `identical()` TRUE |
| `.se` vs `frm_linpred(se.fit = TRUE)` | max relative 2.2e-16 |
| self-check `cov_rel_error` at the unseen level | 2.220446e-16, passes |
| (b) `.se` vs the seam restricted by hand to the 12 intercept and `s(t)` columns | `identical()` TRUE, max relative 0 |
| (a) `.estimate` vs mgcv `predict.gam(exclude = "s(t,subject)")`, ML fit | max abs diff 7.9e-07, against a curve range of 1.006 |
| (a) `.se` ratio frm / mgcv | 1.000 to 1.004, median 1.002 |

The (b) identity is an IDENTITY, not a measurement: with the `fs`
columns exactly zero, `A V A'` over 112 columns and over the 12 kept
columns are the same sum plus exact zeros. `identical()` is the numerical
check of it on this BLAS (R's reference BLAS). The (a) comparison is a
measurement against an independent fit.

For contrast, a BAR term's unseen level does carry variance: on
`v ~ s(t, k = 12) + (1 | subject)` at `re_formula = NULL`, `extra_var` is
0.001328057 on every row, which is `sd(subject)^2 = 0.03644^2`.

Every printed number of the vignette's population curve is the same as
0.64.0's rendered article (`docs/frmtmb.spline/articles/curve-inference.md`
at commit 9b4bb650), where `re_formula = NA` dropped the `fs` term:
simultaneous critical value 2.710666, mcse 0.011, first-row `.se`
0.05442608. The "about 2.7" and "2.72 from someone else's 2.70" in the
prose stand. See "Vignette" for the full comparison.

## The difference guard

The new argument opened a silent wrong answer, found while wiring it.
An unseen level of a bar term loads NO column of `A` and carries its
marginal variance in `extra_var`. `sp_same_latent()` decides whether two
grids of a difference load the same latent draw by comparing their
non-fixed design columns and their `extra_var`. Two DIFFERENT unseen
levels match on both, so the predicate returns TRUE and the variance
cancels.

Measured (`dev/splinecurve-band.R` section 6 and
`dev/splinecurve-noguard.R`, same data, `(1 | subject)`, grids at unseen
levels "A" and "B", `re_formula = NULL`):

- `sp_same_latent()` returns TRUE.
- With the guard deleted, `frm_curve(contrast = )` returns `.se` exactly
  0 on all 5 rows. The right value is at least
  `sqrt(2 * extra_var) = 0.0515375`.

The guard refuses any contrast with `allow_new_levels = TRUE` and nonzero
`extra_var` in either grid. It fails closed. What it costs:

- A contrast of the SAME unseen level across a fixed effect is refused
  though its draw would cancel. `re_formula = NA` gives the same
  difference there, and the message says so.
- An exact `gp()` contrast with `allow_new_levels = TRUE` and no unseen
  level is refused. `allow_new_levels = FALSE` answers it, and the
  message says so.

An unseen `fs` level carries no `extra_var`, so the vignette's use and
any `fs` contrast are not affected.

## Vignette

`vignettes/curve-inference.Rmd`:

- Section 1: the grid names `subject = pop`, a level the fit did not
  see, and passes `allow_new_levels = TRUE`. The prose that said
  `re_formula = NA` drops the per-subject curves is replaced by what
  `NA` does now and how to get the population curve. A paragraph says
  what the band means (no between-subject variance).
- The "check" paragraph said frmtmb exports no route to the grid
  covariance and that the two routes agree "to about twelve significant
  figures". Both were stale since frmtmb 0.52.0 (`frm_lp_basis()`), and
  0.64.0's own render printed 2.2e-16. Now: the covariance comes from
  `frm_lp_basis()` and agrees with `frm_linpred(se.fit = TRUE)` to
  machine precision.
- The "ncall" paragraph described the old perturbation rebuild
  ("skipped in blocks"). Now: one call, the check, on 110 random
  coefficients.
- Sections 2 and 3 (derivative, `eps`, peak, onset): every chunk passed
  `g2` without `subject`; `g2` now carries `subject = pop` and each call
  passes `allow_new_levels = TRUE`.
- "Per-subject curves": said "Drop `re_formula = NA`" to get a subject's
  curve, which is no longer the switch. Now: name a seen subject and
  leave out `allow_new_levels`; the chunk no longer passes
  `re_formula = NULL`.
- "What to report": one bullet added.
- Teaching points kept: pointwise against simultaneous band, the Monte
  Carlo standard error, the self-check.

No number in the prose moves. Every `#>` output line of the built
vignette (54 lines, from the tarball's `inst/doc/curve-inference.html`)
was diffed against 0.64.0's render
(`docs/frmtmb.spline/articles/curve-inference.md` at 9b4bb650), with
runs of whitespace collapsed. The ONLY difference is the new `subject`
column in the printed curve, which shifts where the table wraps; every
value is the same, including the critical value 2.710666, its mcse
0.01125959, the check 2.220446e-16, the rising and falling intervals,
the peak 0.513161 (se 0.003104161), both onset crossings and subject 3's
peak. Files: `dev/splinecurve-vignette-064.txt`,
`dev/splinecurve-vignette-after.txt`.

## Tests

`tests/testthat/test-new-levels.R`, 6 blocks, 37 expectations, on a
smaller copy of the vignette design (`set.seed(4)`, 12 subjects, 1200
rows):

1. At an unseen `fs` level `frm_curve(allow_new_levels = TRUE)` equals
   `frm_linpred(allow_new_levels = TRUE, se.fit = TRUE)` in estimate and
   standard error, and the self-check ran and passed.
2. The population band carries no between-subject variance:
   `extra_var == 0`, and the gap between the unseen-level `.se` and the
   hand restriction to the intercept and `s(t)` columns is under 1e-8
   times the gap a SEEN subject shows. That scale is measured by the
   run.
3. The default refuses the unseen level by name, class
   `frmtmb_new_levels`, in all three functions; a grid without `subject`
   is still refused with `TRUE`.
4. The flag is validated (`NA`, `"yes"`, `c(TRUE, TRUE)`, `1`) in all
   three functions.
5. The derivative and the feature at the unseen level; a stored curve
   carries the flag; a stored `FALSE` curve with a new unseen grid and a
   caller's `TRUE` is honored.
6. The difference guard, with `sp_same_latent()` asserted TRUE on the
   two unseen levels so the line flips if the predicate ever learns to
   tell them apart.

No absolute numeric tolerance: `expect_equal()` defaults,
`identical()`, and one ratio to a run-measured gap.

Seen to fail:

- Against `rellib-r3` (0.8.1 unfixed): pass=6 fail=0 err=6 skip=0, every
  block an "unused argument" error. That is the weak form. The
  behavioral form is the vignette's `R CMD build` failure and the five
  refusals in `dev/splinecurve-before.log`. Log:
  `dev/splinecurve-test-new-levels-before.log`.
- Block 6 without the guard: pass=36 fail=1, "Expected `frm_curve(...)`
  to throw a error", and the unguarded `.se` is 0 on every row. Log:
  `dev/splinecurve-noguard.log`.
- Against the lane library: pass=37 fail=0 err=0 skip=0. Log:
  `dev/splinecurve-test-new-levels-after.log`.

## Verification

Whole frmtmb.spline suite, one file per R process, against
`wt-splinecurve-lib` (frmtmb 0.65.0 installed from this tree, spline from
this lane). Driver `dev/splinecurve-suite.ps1`, runner
`dev/splinecurve-run-tests.R`, summary `dev/splinecurve-suite-after.log`,
per-file logs in `dev/splinecurve-suite-after/`. 15 of 15 files ran:

| File | pass | fail | err | skip |
|---|---|---|---|---|
| test-bracket-access.R | 1 | 0 | 0 | 0 |
| test-conditions-census.R | 11 | 0 | 0 | 0 |
| test-conditions.R | 5 | 0 | 0 | 0 |
| test-curve.R | 65 | 0 | 0 | 0 |
| test-deriv.R | 38 | 0 | 0 | 0 |
| test-difference.R | 59 | 0 | 0 | 0 |
| test-frailty.R | 30 | 0 | 0 | 0 |
| test-gratia.R | 13 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 |
| test-new-levels.R | 37 | 0 | 0 | 0 |
| test-royston-parmar.R | 67 | 0 | 0 | 0 |
| test-rp-floored.R | 81 | 0 | 0 | 0 |
| test-scale.R | 0 | 0 | 0 | 1 |
| test-span.R | 73 | 0 | 0 | 0 |
| test-surface.R | 47 | 0 | 0 | 0 |

531 passed. The one skip is the scale tier, which needs
`FRMTMB_SCALE_TESTS=true` by design.

CORRECTION (round 2). Round 1 of this file said "611 passed". The table
above sums to 531, and so does the reviewer's rerun. The 611 was typed
into the sentence by hand from the per-file lines, not generated from
the log, which is the failure `dev/lane-rules.md` names ("Generate
counts into the document; do not type them"). The round-2 totals below
were summed by a script over the RESULT lines of the log and pasted:
the same script over the round-1 log gives `files 15 pass 531 fail 0
err 0 skip 1`.

`R CMD build` then `R CMD check --as-cran`, once, in
`dev/splinecurve-check/final/`, driver `dev/splinecurve-check.ps1`:
build OK (vignettes rebuilt, log `dev/splinecurve-build-after.log`),
check `Status: 1 NOTE` (log `dev/splinecurve-check-after.log`). The NOTE
is the expected one: "Skipping checking math rendering: package 'V8'
unavailable". Tests OK in 39 s, vignette re-building OK, PDF manual OK.

## Round 2: the review's blockers

Review: `dev/reviews/2026-09-29-splinecurve.md`, verdict NOT MERGEABLE,
with the fs population curve itself holding on every check.

### Blocker 1: an unseen level of a bar term

With `allow_new_levels = TRUE` and a `re_formula` that keeps `(1 | g)`
or `(1 + t | g)`, the three non-pointwise routes answered with wrong
numbers where 0.8.1 refused (review Finding A, `rev-02`: simultaneous
critical value 0.646 and 0.615 against 2.146 and 2.472 from the full
covariance, coverage 0.352 and 0.148; derivative se 0.40 to 0.50 of the
right one; feature `.value_se` 0.0091 against 0.0376). Cause:
`frm_lp_basis()` returns a new level's variance one number per row, as
`extra_var`, while the grid's covariance is `A V A' + Z S Z'`, one draw
shared by every row. The simulation, the derivative and the feature use
`A V A'` and so omit it. My round-1 claim that the flag "changes no
number" was true only of seen rows, and round 1 did not look past the
contrast route. The finding is the reviewer's.

Fix, as instructed: REFUSE. New `sp_new_level_stop()` in
`R/curve-cov.R`, called from `frm_curve()` when `simultaneous = TRUE`,
from `frm_curve_deriv()` and from `frm_curve_feature()` (after the root
search, so a window with no root still returns its zero-row answer). It
fires when `allow_new_levels = TRUE` and any row has nonzero
`extra_var`, and names the grouping factor and its unseen levels
through the public `ngrps()` and `ranef()` (an absent column is named as
absent; a term it cannot match by name, or a `gp()`, gets a generic
phrase). The pointwise `frm_curve()` band is still answered: it adds
`extra_var` row by row and the reviewer measured it right to 0
relative. `?frm_curve` gains the section "An unseen level of a `(1 | g)`
term", which says what is answered, what refuses and why, and that the
`"Sigma"` attribute is `A V A'` without the new level's variance. The
round-1 sentence presenting "a new subject's curve at `re_formula =
NULL`" as supported is gone.

Seen to fail, on the round-1 build before the round-2 code was
installed (`dev/splinecurve-r2-before.log`, the new test file;
`dev/splinecurve-r2-before-numbers.log`, script
`dev/splinecurve-r2-before.R`, fixture `set.seed(4)`, 12 subjects, 1200
rows, 30-point grid at an unseen level, `re_formula = NULL`):

- test file: pass=47 fail=12. EIGHT of the failures are blocker 1,
  behavioral: "Expected `frm_curve(...)` / `frm_curve_deriv(...)` /
  `frm_curve_feature(...)` to throw a error", four per model, because
  each call ANSWERED. Two are blocker 2 (next section). The other two
  were a defect in my new test, not in the package: the pointwise
  identity compared a named vector with an unnamed one. Fixed with
  `unname()` before the round-2 run; the eight behavioral failures do
  not depend on it.
- numbers, same build: simultaneous critical value 0.91152 on
  `(1 | subject)` and 0.88297 on `(1 + t | subject)`, both BELOW
  `qnorm(0.975)` = 1.96; feature `.value_se` 0.014337 against
  `frm_curve()`'s 0.043118 at the same `t`, and 0.014367 against
  0.04466.
- after (`dev/splinecurve-r2-after-numbers.log`, same script on the
  round-2 build): all six calls REFUSED with the new message.

### Blocker 2: the stored-curve rule

Round 1 combined a stored curve's flag with the caller's by OR, so an
EXPLICIT `allow_new_levels = FALSE` was overridden. Now the caller's
value wins when the argument is supplied: `frm_curve_deriv()` and
`frm_curve_feature()` pass `!missing(allow_new_levels)` to `sp_spec()`,
and the stored value applies only when it is missing. The round-1 case
that motivated OR (stored `FALSE`, new unseen grid, caller's `TRUE`)
still works and is still tested.

Seen to fail on the round-1 build: the two explicit-`FALSE`
expectations failed ("Expected `frm_curve_deriv(...)` to throw a error
with class <frmtmb_new_levels>", and the same for
`frm_curve_feature()`), and the numbers script shows both calls
ANSWERED, the feature with `.value_se` 0.03724812. After: both refuse
with core's `New levels in the factor-smooth term s(t,subject): new`.

### Record and docs

- The suite total: see the correction under Verification.
- `?frm_curve`: the stale section "What the covariance is, and how it
  is checked" is deleted and its still-true content is merged into
  "The route to the covariance": `A`, `V` and `extra_var` come from one
  `frm_lp_basis()` call, `Sigma` is `A V A'`, the pointwise standard
  error is `sqrt(diag(Sigma) + extra_var)`, and the check agrees with
  `frm_linpred(se.fit = TRUE)` to machine precision (2.2e-16 on the
  vignette's model).
- `R/zzz.R`: the `frm_curve` / `predict` row now says the same; the
  `frm_curve` / `gp` row now says the pointwise band carries the kriging
  variance and the other three routes do not.
- Two sentences the review asked for: two unseen `fs` levels contrast to
  exactly 0 with se 0, which is the population curve twice and NOT the
  difference between two new subjects; and a model with both an `fs`
  term and `(1 | g)` read at an unseen level with `re_formula = NULL`
  gets population plus the `(1 | g)` variance in its pointwise band,
  while the other routes refuse.
- Rendered with `Rd2txt` and read: `dev/splinecurve-rd-r2.log`.

### Round 2 verification

`test-new-levels.R` now has 8 blocks. Against the round-2 build: pass=59
fail=0 err=0 skip=0 (`dev/splinecurve-r2-test-new-levels-after.log`).

Whole suite, one file per R process, driver
`dev/splinecurve-suite.ps1 -Tag r2`, summary
`dev/splinecurve-suite-r2.log`, per-file logs
`dev/splinecurve-suite-r2/`. Total from `dev/splinecurve-sum.sh`:
`files 15 pass 553 fail 0 err 0 skip 1`.

| File | pass | fail | err | skip |
|---|---|---|---|---|
| test-bracket-access.R | 1 | 0 | 0 | 0 |
| test-conditions-census.R | 11 | 0 | 0 | 0 |
| test-conditions.R | 5 | 0 | 0 | 0 |
| test-curve.R | 65 | 0 | 0 | 0 |
| test-deriv.R | 38 | 0 | 0 | 0 |
| test-difference.R | 59 | 0 | 0 | 0 |
| test-frailty.R | 30 | 0 | 0 | 0 |
| test-gratia.R | 13 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 |
| test-new-levels.R | 59 | 0 | 0 | 0 |
| test-royston-parmar.R | 67 | 0 | 0 | 0 |
| test-rp-floored.R | 81 | 0 | 0 | 0 |
| test-scale.R | 0 | 0 | 0 | 1 |
| test-span.R | 73 | 0 | 0 | 0 |
| test-surface.R | 47 | 0 | 0 | 0 |

`R CMD build` then `R CMD check --as-cran`, once to completion, on
frmtmb.spline 0.9.0 in `dev/splinecurve-check/final-r2/`, driver
`dev/splinecurve-check-r2.ps1`: build OK (`dev/splinecurve-build-r2.log`),
check `Status: 1 NOTE` (`dev/splinecurve-check-r2.log`), the same V8
math-rendering NOTE. Tests OK in 33 s, vignette re-building OK. The 54
`#>` lines of the rebuilt vignette are identical to round 1's
(`dev/splinecurve-vignette-r2.txt`).

The test title of block 7 was shortened to fit 80 columns after the
suite and the check ran; nothing else changed after them. That file was
rerun alone afterwards: pass=59 fail=0 err=0 skip=0
(`dev/splinecurve-r2-test-new-levels-final.log`).

One aborted start, recorded so the log is not misread: the first
round-2 driver was a sed copy of the round-1 one whose directory
substitution silently failed, so it pointed at `final/` and globbed
`frmtmb.spline_*.tar.gz`, which would have picked round 1's 0.8.1
tarball. I stopped it during the build, before any check ran and before
any tarball was written; no R process of the lane was left. The driver
now names the tarball from `DESCRIPTION`'s version.

## Not done, and why

- Core is unchanged. Nothing in core needs to change for this lane: an
  unseen `fs` level adds no variance, which is the right answer for the
  population curve.
- Round 1 reported two stale doc passages as outside the lane. The
  review handed them to the lane and round 2 fixed both (below). The
  NEWS sentence round 1 flagged was already reworded by the
  consolidating session.
- The new-level variance is REFUSED on three routes, not added. Adding
  it to the simulation needs the covariance between rows, `Z S Z'`,
  which `frm_lp_basis()` does not return; the coordinator ruled it out
  of this round.

### For the consolidating session to file

The lane may not edit `dev/test-backlog.md`. Please file:

- **Exact `gp()` kriging variance on the three non-pointwise routes.**
  Under `allow_new_levels = FALSE` an exact `gp()` evaluated off the
  observed positions still reaches `frm_curve(simultaneous = TRUE)`,
  `frm_curve_deriv()` and `frm_curve_feature()` with a nonzero
  `extra_var` that each route omits: the simulation draws from `A V A'`
  and divides by a standard error that includes it, and the derivative
  and feature standard errors are built from `A V A'` alone. It is the
  same mechanism as round 2's blocker 1, older than this lane, and
  small where measured: `extra_var` 7.3e-07 to 8.5e-07 on `y ~ gp(x)`
  over 90 points (`R/zzz.R`, `frm_curve` / `gp` row). The fix that
  makes every route right is a core seam returning the between-row
  block (`k(x1, x2) - Xr1 K Xr2'` for `gp()`, `Z S Z'` for a new level)
  rather than its diagonal. Under `TRUE` the lane's refusal already
  covers `gp()` too.

## Found by the release scale tier, after the review

`tests/testthat/test-scale.R`, which runs only under
`FRMTMB_SCALE_TESTS=true` and so skipped in this lane's suite, drew
the population curve of `s(t) + s(t, subject, bs = "fs")` on a grid
with no `subject` column: the same defect as the vignette. The
consolidating session gave the grid an unseen subject and passed
`allow_new_levels = TRUE` to its five curve calls. Before: err=1 with
core's "needs the grouping column" error; after: pass=2
(`dev/release/scale-spline.log`). A gated file is the one a lane
suite never reaches, which is what the release scale tier is for.
