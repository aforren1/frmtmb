# Lane splitsurv: curve inference into core, frmtmb.survival, frmtmb.spline retired

Base 19b162b9 (frmtmb 0.73.1). Worktree `frmtmb-wt-splitsurv`, branch
`wt-splitsurv`, library `C:/Users/adf44/source/r/wt-splitsurv-lib`.
Nothing is committed. This file is the record; the logs it cites are
local (`dev/splitsurv-*-out/`, gitignored).

## What changed, in one paragraph

`frm_curve()`, `frm_curve_deriv()`, `frm_curve_feature()`, their print
methods, their compat rows and the curve-inference article moved from
frmtmb.spline into core. The one Royston-Parmar branch in them, the
`rp_floored()` gate, became a family hook, `post$curve_check`. A new
extension, `extensions/frmtmb.survival`, holds `cox()` and
`cox_baseline()` (from core) and `royston_parmar()` and `rp_floored()`
(from frmtmb.spline). Core's cox-specific branches (the baseline's
chart, its face search, its draw names on the sampler, its `sbhaz`
prior class and default row) became one seam,
`frmtmb_family(simplex_pars =)`, and the baseline's data became
`frmtmb_family(aterm_data =)`. `extension_family_homes` gained
`cox = "frmtmb.survival"`. Nothing is left in frmtmb.spline, and it is
retired in this worktree (section "Retirement").

## Guardrail 3, measured before any move

Script `dev/splitsurv-measure.R`, output `dev/splitsurv-out/measure.txt`.
Lines are from `git show 19b162b9:<path>`; "code" is a line that is not
blank and not a comment. Test time is the sum of per-test `real`
seconds from the before run on rellib-r11 (`dev/splitsurv-par.sh`, 8
processes at once on a shared machine, `NOT_CRAN=true`, logs in
`dev/splitsurv-before-out/`), so it is an order of magnitude, not a
benchmark. Commits are `git log --oneline 4f5ea39f..19b162b9 -- <paths>`
(4f5ea39f is 0.68.1, so the range covers the 0.69 to 0.73 rounds, 37
commits in all).

| | (a) curve inference | (b) Royston-Parmar | (c) `cox()` |
|---|---|---|---|
| lines in `R/` (all / code) | 1846 / 805 in `R/curve*.R`, plus 43 / 40 of compat rules | 1139 / 411 in `royston-parmar.R`, `rp-basis.R`, `rp-check.R`, plus 37 / 36 of compat rules | 324 / 174 in `R/families.R` (family, bases, `cox_baseline()`, docs, constructor), 12 lines of compat rules, and branches in `objective.R` (10 lines name cox or bhaz), `fit.R` (18), `priors.R` (30) |
| test time, share of the source suite | 62.7 s over 73 tests, 72.9% of frmtmb.spline's 86.0 s | 22.1 s over 36 tests, 25.7% of frmtmb.spline's 86.0 s | 1.2 s over 13 tests, 0.09% of core's 1428.7 s (2742 tests); and frmtmb.sample's `test-cox-draws.R`, 6.0 s, 1.00% of that suite |
| commits in the last five rounds | 1 to the code (f5770ca0, accurate SEs on ill-conditioned designs, which also added `test-illcond.R`) | 0 | 1 substantive (58b97a2c, the cox chart, 0.70.0 round); 3 more touched `cox` in a family list only (2edf7ab5, 12b92070, 0a8d9a02, from `git log -G"cox|bhaz"`) |
| dependencies only it needs | gratia (tests; was a frmtmb.spline Suggests, now a core Suggests). mgcv is already a core Import | flexsurv, rstpm2 (tests and vignette); survival (tests) | splines2 (tests only; removed from core Suggests); survival is used by other core tests too |
| ported brms rows that depend on it | 0 | 0 (not a brms feature) | 0. brms's own suite tests cox in `tests.stancode.R:1752` and `tests.standata.R:1017`; neither is ported (`dev/brmsport-verdicts.tsv` has no cox, bhaz or sbhaz row) |

Reading the numbers against the rule. (a) is large, stable and needs one
test-only dependency; it moves INTO core, so the measure is what core
gains: 805 code lines and about a minute of tests. (b) is stable (no
commit in five rounds) and is the bulk of frmtmb.survival. (c) is small
(174 code lines in the family) and stable, needs no dependency of its
own beyond a test package, and its tests are 0.09% of core's time. By
`dev/core-boundary.md`'s own words "a feature that is small, stable,
and needs no extra dependency is a weak candidate", so the case for (c)
is the layout, not a measured cost: it is the one family that made core
carry family-specific branches in four files, and the seam that
replaces them (`simplex_pars`) is the measurable gain. The user decides.

## The seams

Every one is public and documented, and frmtmb.survival uses only these
and other exported functions (guardrail 2).

| seam | where | what it replaces | useful to a third party because |
|---|---|---|---|
| `frmtmb_family(simplex_pars = c(<brms name> = "<extra par>"))` | `R/families.R`, read by `family_simplex_terms()` in `R/objective.R` | the `Zbhaz` branches of `mo_chart_frame()`, `simplex_frame_terms()`, `mo_search_terms()`, `cox_face_leave()` (now `simplex_face_leave()`), the `sbhaz` default prior row and `sbhaz_targets()` (now `family_simplex_targets()`) | any family with a simplex of weights gets the mo() chart, the pole and face search, the sampler's softmax chart and brms draw names, and a Dirichlet prior class of its own name. Shown by `tests/testthat/test-simplex-pars.R`, a location mixture over fixed centers that is no cox() at all |
| `frmtmb_family(aterm_data =)` | `R/families.R` | the `fam[["aterm_data"]] <-` slot cox() set by hand after construction | data built once from the observed response (a basis over observed times); the slot existed, now it is an argument with a contract |
| `post$curve_check(fit, resp)` | documented in `?frmtmb_family`, run by `curve_family_check()` in `R/curve.R` | `sp_rp_gate()`, the call to `rp_floored()` by name inside the curve functions | a family whose likelihood can be a floor artifact refuses a curve off such a fit. Shown by `tests/testthat/test-curve-check.R` on a family of the test's own |

New exports in core: `frm_curve()`, `frm_curve_deriv()`,
`frm_curve_feature()` (moved), and the S3 methods `print.frmtmb_curve`,
`print.frmtmb_feature`. Removed from core: `cox()`, `cox_baseline()`.
No `frmtmb:::` export was added for frmtmb.survival.

`set_prior()` and a family's simplex name. `set_prior()` reads a class
that is not frmtmb's own as a dpar's name and stores it as that dpar's
`Intercept`. A declared simplex name (`"tw"` in the test) arrives that
way, so `prior_simplex_check()` gives it its class back at the fit, and
a `dirichlet()` on a name that is no simplex is refused there instead
of at `set_prior()`. This is the one behavior change outside the moves
(BREAKING in NEWS): `set_prior("dirichlet(2)", class = "sigma")` used
to refuse when written and now refuses in `frm()`; the message adds
"or a simplex a family declares" to the classes it names. The print
methods of a curve and a feature now call `frm_check_dots()`, as every
core print method does (`test-arg-refusal.R` requires it); in
frmtmb.spline they ignored an unknown argument. `"sbhaz"` keeps its `set_prior()`-time checks, as brms's
class. `check_simplex_pars()` refuses a simplex name that is one of
frmtmb's own classes or one of the family's dpars.

What did NOT become a seam, and why:

- `fam_lcdf()` / `fam_lccdf()`, the arity shim that passes the extra
  parameters to a 4-argument `lcdf` / `lccdf`. It was never cox-specific
  (it reads `formals()`), so it stays as it is.
- `check_mv_extra_family()` refuses any family with extra parameters in
  a multivariate model except ordinal and logistic_normal. cox() is one
  of the refused; the guard is generic and stays.
- `links-brms.R`, generated from brms, lost its two `cox` lines
  (`dev/famlink-gen-links.R` no longer lists it). frmtmb.survival
  carries brms's set as `cox_links` and `test-cox-links.R` checks it
  against the installed brms.
- `R/predict.R` keeps one doc sentence and one comment naming cox, both
  still true; it is the mvnlatent lane's area, so it was left alone.

## Evidence

### Guardrail 1, both arms

- Core, `tests/testthat/test-cox-home.R`: without frmtmb.survival,
  `frm(family = "cox")`, `brmsfamily("cox")` and
  `brmsfamily("Cox", link = "identity")` stop naming frmtmb.survival, and
  core's registry has no cox; with it, `family_ctor("cox")` is
  frmtmb.survival's `cox`, and `brmsfamily("cox")` builds it. The
  session-dependent arm skips when the other is the one in force, as
  `test-families-extreme.R` does for wiener.
- frmtmb.survival, `tests/testthat/test-cox-home.R`: both arms in child
  processes, so neither depends on the session. Without the package:
  "The cox family is in the extension package frmtmb.survival, which is
  not loaded. Install it and run library(frmtmb.survival)" for the
  string and for `brmsfamily()`, and the child confirms it never loaded
  the package (the absent case). With it: `family = "cox"`,
  `family = brmsfamily("cox")` and `family = cox()` give the same fit,
  `identical(logLik())`.

### frm_compat() before and after

Script `dev/splitsurv-compat-dump.R` (before: core + frmtmb.spline on
rellib-r11; after: core + frmtmb.survival on the lane build), compared by
`dev/splitsurv-compare.R`, keyed by the unordered pair:

<!-- splitsurv:compat -->
```
== frm_compat()
pairs before 7955 after 7955 ; same pair set: TRUE 
cox                  pairs 86; status differs 0; note differs 0
royston_parmar       pairs 86; status differs 0; note differs 1
rp_floored           pairs 135; status differs 0; note differs 0
frm_curve            pairs 135; status differs 0; note differs 4
frm_curve_contrast   pairs 135; status differs 0; note differs 1
all pairs: status differs 0 ; note differs 6 
                      key status.b status.a
    autoscale x frm_curve    works    works
  fitted x royston_parmar  refused  refused
           frm_curve x nl    works    works
           frm_curve x rr    works    works
          frm_curve x s()    works    works
 frm_curve_contrast x s()    works    works
features identical (name, kind): TRUE 
unresolved before: 0 after: 0 
```
<!-- /splitsurv:compat -->

Six notes differ, and no status, on all 7955 pairs. Five are core's
own frm_curve rows, reworded where they said "this package" (now
"frm_curve()"). The sixth is `royston_parmar x fitted`: its note said
"core::cox() refuses for the same reason", which the move made false,
and it says "cox(), in this package too, refuses for the same reason"
now (review finding 6). The 86 `cox` pairs are equal in status and
note, and the 86 `royston_parmar` pairs in status, and in note but
for that one. Three `cox` rows (trials(), vint(),
vreal()) are derived from the family's `accepts_aterms`; frmtmb.survival
reproduces them with the exported `compat_aterm_rules()`.

### frm_curve() and companions, base build against lane build

Script `dev/splitsurv-curve-dump.R` (cases from the moved tests: a
`s(x)` gaussian curve with a simultaneous band, the transform, first
and second derivatives, a maximum and a crossing, a difference curve
and its derivative; on Royston-Parmar fits the `mu` curve and
derivative, `rp_floored()` reports, the gate's refusal from all three
functions on a non-monotone fit, and the flexsurv `bc` data). Compared
by `dev/splitsurv-compare.R`:

<!-- splitsurv:curve -->
```
== curve outputs
curve              columns identical: TRUE; attributes identical: TRUE; bitwise
curve_transform    columns identical: TRUE; attributes identical: TRUE; bitwise
deriv1             columns identical: TRUE; attributes identical: TRUE; bitwise
deriv2             columns identical: TRUE; attributes identical: TRUE; bitwise
feat_max           columns identical: TRUE; attributes identical: TRUE; bitwise
feat_cross         columns identical: TRUE; attributes identical: TRUE; bitwise
contrast           columns identical: TRUE; attributes identical: TRUE; bitwise
contrast_deriv     columns identical: TRUE; attributes identical: TRUE; bitwise
rp_curve_mu        columns identical: TRUE; attributes identical: TRUE; bitwise
rp_deriv_mu        columns identical: TRUE; attributes identical: TRUE; bitwise
rp_floored_report  identical: TRUE
rp_bad_report      identical: TRUE
rp_bad_curve       error message identical: TRUE
rp_bad_deriv       error message identical: TRUE
rp_bad_feature     error message identical: TRUE
bc_curve_gamma1    columns identical: TRUE; attributes identical: TRUE; bitwise
bc_curve_mu        columns identical: TRUE; attributes identical: TRUE; bitwise
bc_deriv_mu        error message identical: TRUE
bc_deriv_mu_pw     columns identical: TRUE; attributes identical: TRUE; bitwise
```
<!-- /splitsurv:curve -->

### cox() and royston_parmar() fits, base build against lane build

Script `dev/splitsurv-fit-dump.R`, compared by
`dev/splitsurv-fit-compare.R`: log likelihood, `opt$par`, objective,
convergence, standard errors, `fixef_by_dpar()`, `cox_baseline()`,
the mo() search record and the warnings.

<!-- splitsurv:fits -->
```
cox_cens           identical: TRUE; warnings identical: TRUE
cox_face40         identical: TRUE; warnings identical: TRUE
cox_frailty        identical: TRUE; warnings identical: TRUE
cox_df6_identity   identical: TRUE; warnings identical: TRUE
cox_trunc          identical: TRUE; warnings identical: TRUE
cox_prior          identical: TRUE; warnings identical: TRUE
cox_string         identical: TRUE; warnings identical: TRUE
cox_brmsfamily     identical: TRUE; warnings identical: TRUE
cox_get_prior      identical: TRUE
rp_bc              identical: TRUE; warnings identical: TRUE
rp_bc_odds         identical: TRUE; warnings identical: TRUE
```
<!-- /splitsurv:fits -->

### No reach into frmtmb from frmtmb.survival

<!-- splitsurv:grep -->
```
$ grep -rn "frmtmb:::\|getFromNamespace\|asNamespace(\"frmtmb\")\|asNamespace('frmtmb')" extensions/frmtmb.survival/
(no output; grep exit status 1)
```
<!-- /splitsurv:grep -->

`extensions/frmtmb.survival/tests/testthat/test-no-reach.R` asserts the
same on every source, test and vignette file, with the inverse case
(each pattern matches the reach it describes).

### The floor

`bash dev/splitsurv-floor.sh <scratch>` installs frmtmb.survival over
rellib-r11's frmtmb 0.73.0 from a copy with the floor lowered. It does
not install: `frmtmb_register_compat() called the rules = builder ...
unused arguments (aterm_data = ...)`. So frmtmb.survival needs this
round's core; its DESCRIPTION says `frmtmb (>= 0.73.1)` as a
placeholder, and the number must be the core release that carries this
change.

<!-- splitsurv:tests -->
### Every test file, one per process, lane build

`dev/splitsurv-par2.sh` over `dev/splitsurv-after.lst` (every test file
of core and of the seven extensions, frmtmb.spline retired),
`NOT_CRAN=true`, library `wt-splitsurv-lib` first, run whole after
fix round 1, which changed paths every suite reaches (frame
assembly, the objective, the standard errors). Every log's `lib:`
line names the lane's frmtmb.

Files: 432 listed, 432 with a RESULT line.

| package | files | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 273 | 19194 | 0 | 0 | 203 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 32 | 1838 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 11 | 368 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 430 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 65 | 3026 | 0 | 0 | 9 | 0 |
| frmtmb.survival | 14 | 376 | 0 | 0 | 0 | 0 |

Without a RESULT line (every block gated off, as on the base):  

Files with a failure, an error or an escaped warning: none 

The base run (rellib-r11, base tree, same driver): 407 files with a RESULT line; failing: test-dirichlet-multinomial.R (a 0.73.1 test of the `fitted()` clamp, which rellib-r11's 0.73.0 lacks). The other 15 listed files are the gated brms-suite files, on which that first runner (`dev/splitsurv-run-tests.R`) stopped at an empty result; `dev/splitsurv-run-tests2.R` reports them.

Every new or moved test file, with its counts:

| package | file | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | test-arg-refusal.R | 122 | 0 | 0 | 0 | 0 |
| frmtmb | test-brms-families.R | 1257 | 0 | 0 | 0 | 0 |
| frmtmb | test-ci-siblings.R | 16 | 0 | 0 | 0 | 0 |
| frmtmb | test-compat-extensions.R | 2 | 0 | 0 | 0 | 0 |
| frmtmb | test-compat.R | 767 | 0 | 0 | 0 | 0 |
| frmtmb | test-cox-home.R | 8 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-check.R | 8 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-deriv.R | 48 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-difference.R | 71 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-gratia.R | 13 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-illcond.R | 4 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-new-levels.R | 57 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-scale.R | 0 | 0 | 0 | 1 | 0 |
| frmtmb | test-curve-span.R | 73 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve.R | 85 | 0 | 0 | 0 | 0 |
| frmtmb | test-extension-home.R | 23 | 0 | 0 | 0 | 0 |
| frmtmb | test-famgaps.R | 64 | 0 | 0 | 0 | 0 |
| frmtmb | test-hurdle-nb-zoib.R | 86 | 0 | 0 | 0 | 0 |
| frmtmb | test-mv-gaps.R | 58 | 0 | 0 | 0 | 0 |
| frmtmb | test-simplex-pars.R | 37 | 0 | 0 | 0 | 0 |
| frmtmb | test-simulate-density.R | 620 | 0 | 0 | 0 | 0 |
| frmtmb | test-update-pool.R | 21 | 0 | 0 | 0 | 0 |
| frmtmb.eam | test-simulate-density.R | 69 | 0 | 0 | 0 | 0 |
| frmtmb.ode | test-compat.R | 15 | 0 | 0 | 0 | 0 |
| frmtmb.sample | test-cox-draws.R | 12 | 0 | 0 | 0 | 0 |
| frmtmb.sample | test-simplex-family-notes.R | 4 | 0 | 0 | 0 | 0 |
| frmtmb.sample | test-simulators.R | 54 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-bracket-access.R | 1 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-conditions-census.R | 11 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-conditions.R | 5 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-cox-chart.R | 15 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-cox-home.R | 12 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-cox-links.R | 31 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-cox.R | 45 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-frailty.R | 30 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-message-uniqueness.R | 4 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-no-reach.R | 3 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-royston-parmar.R | 67 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-rp-floored.R | 81 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-rp-surface.R | 44 | 0 | 0 | 0 | 0 |
| frmtmb.survival | test-saved-fits.R | 27 | 0 | 0 | 0 | 0 |

### Tests, accounted

`dev/splitsurv-account.R`, matching tests by name between the base
run and the lane run. A change in a core file's count is the moved
assertions; the list after the table says where each went.

| before file | tests | pass before | after file | tests | pass after |
|---|---|---|---|---|---|
| frmtmb.spline `test-curve.R` | 19 | 82 | frmtmb `test-curve.R` | 20 | 85 |
| frmtmb.spline `test-deriv.R` | 12 | 48 | frmtmb `test-curve-deriv.R` | 12 | 48 |
| frmtmb.spline `test-difference.R` | 13 | 71 | frmtmb `test-curve-difference.R` | 13 | 71 |
| frmtmb.spline `test-gratia.R` | 5 | 13 | frmtmb `test-curve-gratia.R` | 5 | 13 |
| frmtmb.spline `test-illcond.R` | 1 | 4 | frmtmb `test-curve-illcond.R` | 1 | 4 |
| frmtmb.spline `test-new-levels.R` | 8 | 57 | frmtmb `test-curve-new-levels.R` | 8 | 57 |
| frmtmb.spline `test-span.R` | 13 | 73 | frmtmb `test-curve-span.R` | 13 | 73 |
| frmtmb.spline `test-scale.R` | 1 | 0 | frmtmb `test-curve-scale.R` | 1 | 0 |
| frmtmb.spline `test-royston-parmar.R` | 9 | 67 | frmtmb.survival `test-royston-parmar.R` | 9 | 67 |
| frmtmb.spline `test-rp-floored.R` | 15 | 81 | frmtmb.survival `test-rp-floored.R` | 15 | 81 |
| frmtmb.spline `test-frailty.R` | 4 | 30 | frmtmb.survival `test-frailty.R` | 4 | 30 |
| frmtmb.spline `test-surface.R` | 9 | 47 | frmtmb.survival `test-rp-surface.R` | 8 | 44 |
| frmtmb.spline `test-bracket-access.R` | 1 | 1 | frmtmb.survival `test-bracket-access.R` | 1 | 1 |
| frmtmb.spline `test-conditions-census.R` | 2 | 11 | frmtmb.survival `test-conditions-census.R` | 2 | 11 |
| frmtmb.spline `test-conditions.R` | 2 | 5 | frmtmb.survival `test-conditions.R` | 2 | 5 |
| frmtmb.spline `test-message-uniqueness.R` | 1 | 4 | frmtmb.survival `test-message-uniqueness.R` | 1 | 4 |
| frmtmb `test-cox-chart.R` | 3 | 17 | frmtmb.survival `test-cox-chart.R` | 2 | 15 |
| frmtmb `test-famgaps.R` | 21 | 100 | frmtmb `test-famgaps.R` | 13 | 64 |
| frmtmb `test-brms-families.R` | 23 | 1278 | frmtmb `test-brms-families.R` | 23 | 1257 |
| frmtmb `test-hurdle-nb-zoib.R` | 16 | 87 | frmtmb `test-hurdle-nb-zoib.R` | 16 | 86 |
| frmtmb `test-mv-gaps.R` | 9 | 58 | frmtmb `test-mv-gaps.R` | 9 | 58 |
| frmtmb `test-simulate-density.R` | 61 | 621 | frmtmb `test-simulate-density.R` | 60 | 620 |
| frmtmb `test-update-pool.R` | 5 | 22 | frmtmb `test-update-pool.R` | 5 | 21 |
| frmtmb.sample `test-cox-draws.R` | 3 | 12 | frmtmb.sample `test-cox-draws.R` | 3 | 12 |
| frmtmb.sample `test-simulators.R` | 10 | 54 | frmtmb.sample `test-simulators.R` | 10 | 54 |

Tests of a before file whose name is not in its after file, and where
they are now:

- frmtmb.spline `test-surface.R` "print methods say what was checked" (3 passed): frmtmb test-curve.R (3 passed)
- frmtmb.spline `test-bracket-access.R` "no frmtmb.spline function reads a hazard container" (1 passed): NOT FOUND BY NAME
- frmtmb.spline `test-conditions-census.R` "frmtmb.spline raises no bare stop(), warning() or message()" (4 passed): NOT FOUND BY NAME
- frmtmb.spline `test-conditions.R` "a frmtmb.spline refusal is caught as a frmtmb_error" (3 passed): NOT FOUND BY NAME
- frmtmb `test-cox-chart.R` "a face the likelihood rises from is left; the maximum is not" (2 passed): frmtmb test-simplex-pars.R (5 passed)
- frmtmb `test-famgaps.R` "the M-spline and I-spline bases match splines2" (4 passed): frmtmb.survival test-cox.R (4 passed)
- frmtmb `test-famgaps.R` "the I-spline basis is the integral of the M-spline basis" (6 passed): frmtmb.survival test-cox.R (6 passed)
- frmtmb `test-famgaps.R` "cox() matches a hand-rolled M-spline PH likelihood exactly" (5 passed): frmtmb.survival test-cox.R (5 passed)
- frmtmb `test-famgaps.R` "cox() agrees with survival::coxph on the hazard ratio" (2 passed): frmtmb.survival test-cox.R (2 passed)
- frmtmb `test-famgaps.R` "cox() with no censoring is a plain flexible parametric PH fit" (3 passed): frmtmb.survival test-cox.R (3 passed)
- frmtmb `test-famgaps.R` "cox() frailty models come out of the Laplace approximation" (7 passed): frmtmb.survival test-cox.R (7 passed)
- frmtmb `test-famgaps.R` "cox() handles left and interval censoring, and refuses a mean" (6 passed): frmtmb.survival test-cox.R (6 passed)
- frmtmb `test-famgaps.R` "cox() validation" (3 passed): frmtmb.survival test-cox.R (3 passed)
- frmtmb `test-simulate-density.R` "cox() states why it has no simulator rather than only that it has none" (1 passed): NOT FOUND BY NAME

Tests in an after file whose name is not in its before file (new or
rewritten):

- frmtmb `test-curve.R` "print methods say what was checked" (3 passed)
- frmtmb.survival `test-bracket-access.R` "no frmtmb.survival function reads a hazard container" (1 passed)
- frmtmb.survival `test-conditions-census.R` "frmtmb.survival raises no bare stop(), warning() or message()" (4 passed)
- frmtmb.survival `test-conditions.R` "a frmtmb.survival refusal is caught as a frmtmb_error" (3 passed)
- frmtmb `test-simplex-pars.R`: 6 tests, 37 passed (new file)
- frmtmb `test-curve-check.R`: 2 tests, 8 passed (new file)
- frmtmb `test-cox-home.R`: 2 tests, 8 passed (new file)
- frmtmb `test-extension-home.R`: 3 tests, 23 passed (new file)
- frmtmb.survival `test-saved-fits.R`: 1 tests, 27 passed (new file)
- frmtmb.sample `test-simplex-family-notes.R`: 1 tests, 4 passed (new file)
- frmtmb.survival `test-cox.R`: 13 tests, 45 passed (new file)
- frmtmb.survival `test-cox-home.R`: 3 tests, 12 passed (new file)
- frmtmb.survival `test-cox-links.R`: 4 tests, 31 passed (new file)
- frmtmb.survival `test-no-reach.R`: 1 tests, 3 passed (new file)


The 21 assertions core's `test-brms-families.R` lost are cox()'s
rows of three registry loops: 4 links built, 13 roster links
refused (the roster has 17 and cox() takes 4), and 4 of "link and
link_<dpar> agree with links on every registry family". The first
record said 4 + 17 and missed the third test (review finding 1).
frmtmb.survival's `test-cox-links.R` runs all three for cox().

### One session per suite, three orders

`bash dev/release/run-sessions.sh splitsurv 8` with `FRMTMB_RELLIB` the
lane library and `NOT_CRAN=true` (log `dev/release/splitsurv.log`,
copied to `dev/splitsurv-out/sessions-all.log`). LEAK lines without
`[with load of ...]`: none.

```
SESSION frmtmb forward files=263/273 pass=19194 fail=0 err=0 skip=203 warn=0
SESSION frmtmb reverse files=263/273 pass=19294 fail=0 err=0 skip=203 warn=0
SESSION frmtmb preload files=263/273 pass=19387 fail=0 err=0 skip=205 warn=0
SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward files=11/11 pass=368 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent reverse files=11/11 pass=368 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent preload files=11/11 pass=368 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample forward files=60/65 pass=3026 fail=0 err=0 skip=9 warn=0
SESSION frmtmb.sample reverse files=60/65 pass=3026 fail=0 err=0 skip=9 warn=0
SESSION frmtmb.sample preload files=60/65 pass=3026 fail=0 err=0 skip=9 warn=0
SESSION frmtmb.survival forward files=14/14 pass=376 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=14/14 pass=376 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=14/14 pass=376 fail=0 err=0 skip=0 warn=0
SPLITSURV ran 24 of 24
```

Tagged LEAK lines that name a package this lane added:

```
LEAK test-cox-home.R compat features: 2 -> 3 entries [with load of frmtmb.survival]
LEAK test-cox-home.R compat rules: 1 -> 2 entries [with load of frmtmb.survival]
LEAK test-cox-home.R compat features: 3 -> 4 entries [with load of frmtmb.survival]
LEAK test-cox-home.R compat rules: 2 -> 3 entries [with load of frmtmb.survival]
LEAK test-cox-draws.R compat features: 6 -> 7 entries [with load of frmtmb.survival]
LEAK test-cox-draws.R compat rules: 6 -> 7 entries [with load of frmtmb.survival]
LEAK test-simulators.R compat features: 1 -> 4 entries [with load of frmtmb.latent frmtmb.survival]
LEAK test-simulators.R compat rules: 1 -> 4 entries [with load of frmtmb.latent frmtmb.survival]
```

`files=262/272` and `59/64` are the 10 and 5 ported brms-suite
files, gated off without `FRMTMB_BRMS_FIT_TESTS`, which report no
RESULT line.

### OpenBLAS

`dev/ciharden-openblas.sh 0.3.26 lapack`, then the new numerical
tests and the moved curve and cox() tests through that R
(`dev/splitsurv-par3.sh`, `OPENBLAS_NUM_THREADS=4`):

```
frmtmb RESULT test-curve-check.R pass=8 fail=0 err=0 skip=0 warn=0 elapsed=1.26
frmtmb RESULT test-curve-deriv.R pass=48 fail=0 err=0 skip=0 warn=0 elapsed=4.02
frmtmb RESULT test-curve.R pass=85 fail=0 err=0 skip=0 warn=0 elapsed=5.50
frmtmb RESULT test-simplex-pars.R pass=33 fail=0 err=0 skip=0 warn=0 elapsed=1.94
frmtmb.survival RESULT test-cox-chart.R pass=15 fail=0 err=0 skip=0 warn=0 elapsed=1.38
frmtmb.survival RESULT test-cox.R pass=45 fail=0 err=0 skip=0 warn=0 elapsed=1.74
ran 6 of 6
```
<!-- /splitsurv:tests -->

The four "NOT FOUND BY NAME" tests are renamed, not dropped. The three
per-package ones carry the package name in their title and are
frmtmb.survival's `test-bracket-access.R`, `test-conditions-census.R`
and `test-conditions.R` (the last regenerated by
`dev/conditions-ext-tests.R` with a survival refusal, `cox(df = 0)`,
since `frm_curve()` now raises core's class); core's own census,
uniqueness and bracket tests scan the curve code in core's `R/`.
"cox() states why it has no simulator rather than only that it has
none" is frmtmb.survival's "cox() states why it has no simulator, at
both entry points" in `test-cox.R`, which adds the two entry-point
assertions core's rewritten no-simulator test no longer makes for cox.

### R CMD check --as-cran

`bash dev/splitsurv-check.sh <pkg>`, built and checked inside
`dev/splitsurv-check/<pkg>/` with `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`,
pandoc and TinyTeX on PATH, and `R_LIBS` the lane library, rellib-r11,
the user library.

<!-- splitsurv:checks -->
```
== frmtmb (fix round 1)
* checking HTML version of manual ... [22s] NOTE
Status: 1 NOTE
Skipping checking math rendering: package 'V8' unavailable
== frmtmb.survival (fix round 1)
* checking package subdirectories ... NOTE
* checking HTML version of manual ... NOTE
Status: 2 NOTEs
Problems with news in 'NEWS.md':
No news entries found.
Skipping checking math rendering: package 'V8' unavailable
== frmtmb.sample (fix round 1)
Status: OK
== frmtmb.learn (first round; not changed since)
Status: OK
```
<!-- /splitsurv:checks -->

The V8 NOTE is the expected one (dev/lane-rules.md). frmtmb.survival's
NEWS NOTE comes from the heading lane-env.md prescribes,
"# frmtmb.survival (development version)": R reads no version from it,
and the package has no other section yet. It goes away when the release
gives the section its number.

## Retirement of frmtmb.spline

A separate decision for the user. Nothing is left in the package after
the moves; every file has a home:

| frmtmb.spline file | new home |
|---|---|
| `R/curve.R`, `R/curve-cov.R`, `R/curve-deriv.R`, `R/curve-feature.R` | core `R/`, same names; `sp_` helpers renamed `curve_`, `sp_rp_gate()` replaced by `curve_family_check()` |
| `R/royston-parmar.R`, `R/rp-basis.R`, `R/rp-check.R` | frmtmb.survival `R/`, `sp_` renamed `sv_`; `sp_rp_gate()` is now the family's `post$curve_check` |
| `R/zzz.R` | curve rows into core `R/compat.R` (text reworded where it said "this package"); RP and `rp_floored` rows into frmtmb.survival `R/zzz.R`, verbatim |
| `R/frmtmb.spline-package.R` | frmtmb.survival's package page; its history sections stay in git history |
| `vignettes/curve-inference.Rmd` | core `vignettes/curve-inference.Rmd` |
| `vignettes/royston-parmar.Rmd` | frmtmb.survival `vignettes/royston-parmar.Rmd` |
| `dev/scale-row.md` | `dev/scale-row-curves.md` |
| `NEWS.md` | git history; frmtmb.survival's and core's NEWS say what moved |
| `tests/testthat/*` | see "Tests, accounted" |

Files the retirement touches (separable from the moves only where noted):

- deleted: `extensions/frmtmb.spline/` (41 files, listed in
  `dev/splitsurv-out/spline-files.txt`), `.github/workflows/check-frmtmb-spline.yaml`;
- `.github/workflows/test-coverage.yaml`: the `frmtmb.spline` matrix entry
  is replaced by `frmtmb.survival` (same two lines);
- `.github/workflows/scale-helper-identical.yaml`: counts core's copy
  (the curves row moved to core) in place of frmtmb.spline's, still 7;
- `.github/workflows/pkgdown.yaml`, `dev/release/build-docs.R`,
  `dev/release/run-session.R`, `dev/release/run-sessions.sh`,
  `dev/release/run-check.ps1`, `dev/release/run-suite.ps1`: the name
  in the package list is replaced (one token each);
- `dev/release/run-scale.ps1`: the spline row runs from core's
  `tests/testthat/test-curve-scale.R`;
- `_pkgdown.yml`: the menu entry is removed and 11 redirects send the
  old subsite's pages to their new homes
  (`dev/splitsurv-retire-pkgdown.R`, separate from
  `dev/splitsurv-edit-pkgdown.R`, which does the moves' part);
- `dev/conditions-census-gen.R`, `dev/conditions-ext-tests.R`: the
  package entry is frmtmb.survival's;
- `README.md`: the package table row and install line.

Not touched: the historical `dev/release/install-r*.R` scripts and the
many `dev/*-findings.md` and `dev/reviews/*` records that name
frmtmb.spline.

## Overlaps with the other lanes

- `R/families.R` `extension_family_homes`: restructured to one entry per
  line with `cox = "frmtmb.survival"`; mvnlatent adds an entry to the
  same vector, so the hunk conflicts textually and resolves by keeping
  both lines.
- `R/compat.R`: the method vocabulary line gains `frm_curve` and
  `frm_curve_contrast` (line 870), four lines below mvnlatent's
  `mixture_mvn` structure line; the cox rows and group entries are cut.
- `_pkgdown.yml`: `cox_baseline` removed one line from mvnlatent's
  `mixture_mvn`; the menu gains frmtmb.survival where spatial likely
  adds frmtmb.spatial.
- `.github/workflows/test-coverage.yaml`, `dev/release/run-sessions.sh`,
  `dev/release/run-session.R`: package lists, as spatial's.
- `NEWS.md`: a new top section.
- `tests/testthat/test-compat-extensions.R`: the extension list (spatial
  adds frmtmb.spatial).
- `README.md`: the package table and the families paragraph (mvnlatent's
  `mixture_mvn` sentence is two lines above the cox one).
- frmtmb.sample: `DESCRIPTION` Suggests gains frmtmb.survival,
  `check-frmtmb-sample.yaml` installs it, and `test-cox-draws.R` and
  `test-simulators.R` call `frmtmb.survival::cox()` and skip without it.
  Its R code is unchanged: it reads `simplex_frame_terms()`, whose
  `kind` is still `"cox"` for a cox() baseline.
- frmtmb.learn: two texts that named `core::cox()` and
  `frmtmb.spline::royston_parmar()` (a comment and one compat note in
  `R/zzz.R`) now name frmtmb.survival. That changes one learn compat
  note.

## Not done, and why

- The ported brms ledger is not regenerated (consolidation does it). No
  row depends on cox, RP or the curve functions, so none flips.
- `dev/rel073-floor.sh` is the last release's script for frmtmb.sample;
  the next release's floor script needs a block for frmtmb.survival, and
  `dev/splitsurv-floor.sh` is that block.
- The history sections of `frmtmb.spline-package.R` ("What this package
  reads that frmtmb does not promise", "One limit that is core's to
  fix") were not carried: they describe reaches and limits of a package
  that no longer exists. They are in git history.

## Bumps and floors

- frmtmb: minor (two exports moved out, three moved in, two new
  `frmtmb_family()` arguments and a family hook, one prior check moved
  from `set_prior()` to `frm()`).
- frmtmb.survival: new, 0.1.0 in DESCRIPTION. Floor: the core release
  that carries this change (measured above: it does not install on
  0.73.0). DESCRIPTION says `frmtmb (>= 0.73.1)` until that number
  exists.
- frmtmb.sample: patch (a Suggests and two test files; no R code). Its
  floor on core is unchanged.
- frmtmb.learn: patch (one compat note and one comment). Floor
  unchanged.
- frmtmb.spline: retired. If the user keeps it, it is an empty package
  and needs a major bump that says so.
- coupling, eam, latent, ode: unchanged; their suites pass on the lane's
  core from rellib-r11's builds.

## Fix round 1 (review dev/reviews/2026-10-23-splitsurv.md)

- **F1.** The 21 assertions are 4 + 13 + 4 (see "Tests, accounted").
  frmtmb.survival's `test-cox-links.R` gains "link and linkinv agree
  with links on cox()", the 4 assertions of the third loop.
- **F2.** One helper, `require_model_packages()` in core's new
  `R/extension-home.R`, built on `require_extension(pkg, what)`. Both
  are internal and documented there; another lane that moves a model
  part (mixture_mvn, car, spde) calls `require_extension()` with its
  own phrase, or `require_model_packages()` for a family. The family's
  package is the namespace of the function that called
  `frmtmb_family()` (`frm_family_package()`, documented now in
  `?frmtmb_family`); a `mixture()` records its components' packages in
  the attribute `frmtmb_component_packages`. The check runs where core
  calls a model's own functions, found by tracing each failing method
  of the review's case (`dev/splitsurv-saved-trace.R`):
  `assemble_frame()` (frm, update, refits, get_prior), `build_objective()`
  (every retape), `sdr_of()` (summary, print, vcov, confint),
  `response_mean()` (fitted, residuals, predict), `sim_draw()`
  (simulate) and `curve_family_check()` (the curve functions). The
  refusal has class `frmtmb_extension_missing`. Measured with
  `dev/splitsurv-saved.R` (saved fits read with frmtmb alone in a
  scratch library ahead of rellib-r11): print, summary, vcov, fitted,
  residuals, update, confint, simulate (RP) and frm_curve (RP) all name
  frmtmb.survival; logLik, fixef, frm_linpred(type = "link") and
  family() still answer; `simulate()` of cox() gives cox()'s own stated
  refusal, which needs none of its functions. Tests: core
  `test-extension-home.R` (a family a test package "owns", whose
  namespace no session has, at the fit, in a mixture, and in four
  methods, with the guard-absent case), and frmtmb.survival
  `test-saved-fits.R` (real cox() and RP fits saved and read in a child
  process whose library holds a copy of frmtmb and no frmtmb.survival).
- **F3.** frmtmb.sample's default-prior notes announce every simplex a
  family declares, by its brms stem and family name, and no longer read
  `kind == "cox"`. New sample test `test-simplex-family-notes.R` with
  the third-party family of core's `test-simplex-pars.R`; the cox()
  note still reads "sbhaz: dirichlet(1), brms's default".
- **F4.** `check_simplex_extra()` refuses, at frame assembly, a
  `simplex_pars` value that names no parameter `extra_pars()` returns,
  or one with no coordinate, and names it. Not at construction:
  `extra_pars` is a function of the data, so its names exist only once
  the response is in hand. Tested in `test-simplex-pars.R`, with the
  correct declaration beside it.
- **F6.** The `royston_parmar x fitted` note reworded; the compat
  comparison now differs by six notes, not five (section "frm_compat()
  before and after").
- **F9.** `?frm` (R/fit.R) and an internal note in R/confint.R say
  frmtmb.survival.
- **F10.** Redirects for the old subsite's `articles/index.html` (to
  core's articles index, which lists curve-inference) and
  `authors.html` (to frmtmb.survival's). frmtmb.spline's NEWS is carried
  as an article of frmtmb.survival, `vignettes/articles/spline-history.Rmd`
  (pkgdown builds it; `.Rbuildignore` keeps it out of the tarball), with
  every heading two levels down, and `frmtmb.spline/news/index.html`
  redirects there. It is not a section of frmtmb.survival's NEWS.md: R's
  NEWS reader takes a version number from a heading at any level,
  whatever package it names, so frmtmb.spline's 0.1.0 to 0.10.2 would
  have been read as this package's versions (measured with
  `tools:::.build_news_db_from_package_NEWS_md()`: 22 entries under
  frmtmb.spline's version numbers). The comment above `redirects:` and
  `pkgdown.yaml`'s install-order comment are corrected.

### For consolidation (F5, not this round)

- frmtmb.survival's floor: `frmtmb (>= 0.73.1)` is satisfied by the
  released 0.73.1, which lacks `aterm_data =` and `simplex_pars =`.
  Write the release's number.
- `dev/rel073-install.sh` (or its successor) roxygenises and installs
  `frmtmb.spline` (lines 36, 38 and 42) and needs `frmtmb.survival` in
  its place, installed before frmtmb.sample, which suggests it. The
  floor script needs a block for frmtmb.survival; `dev/splitsurv-floor.sh`
  is one.

### Defects found, not fixed (for dev/test-backlog.md)

- brms's `bhaz()` addition term (`time | bhaz(df = 6) ~ x`, with
  `knots =` and `Boundary.knots =`) is refused only generically
  ("Addition term `bhaz()` is not supported"), on the base and on this
  build alike. frmtmb.survival could take it into `cox()`, and refuse
  `gr =` by name, since frmtmb has no stratified baseline.
- `brms::cox(bhaz = list(df = 6))` passed to `frm()` stops with
  "unused argument (bhaz = list(df = 6))", which names neither the
  argument's frmtmb spelling (`cox(df = 6)`) nor a reason.
