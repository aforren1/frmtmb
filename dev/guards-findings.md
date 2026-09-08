# guards lane, 0.55 round

Items 4, 5 and 6 of "Follow-ups carried out of the 0.54.0 round" in
`dev/feature-gaps.md`. All three are fixed. Items 4 and 5 turned out to
be one guard, and fixing 5 without 4 would have broken a model that
works today, so they are one change.

Everything below was measured on this worktree, R 4.6.1, with the
package installed into the lane's private library.

## What changed

| File | Change |
| --- | --- |
| `R/families.R` | `frmtmb_family(se_dpar =)`, its validator `check_se_dpar()`, the reader `family_se_dpar()`, and the declaration on gaussian and student |
| `R/parse.R` | the `se()` scale guard: it reads the family's declaration and refuses a FORMULA for the replaced dpar, by any route |
| `R/predict.R` | `se_unused_sigma()` reports the dpar the FAMILY says `se()` replaces, not the name `sigma` |
| `R/zzz.R` | `frmtmb_notice_state` and `notify_once()`: one-shot session notices |
| `R/frame.R` | the discrete censoring notice, last of the censoring guards |
| `R/compat.R` | the `se()` registry note now names `se_dpar` and the formula rule |
| `vignettes/brms-migration.Rmd` | `se_dpar`, the formula rule, and the notice |
| `tests/testthat/setup.R` | new: the suite runs with `frmtmb.notices` off |
| `tests/testthat/test-custom-family.R` | 7 new tests for items 4 and 5 |
| `tests/testthat/test-cens-trunc.R` | 5 new tests for item 6 |
| `NEWS.md` | one bullet per item under a development heading |
| `man/frmtmb_family.Rd` | roxygenised |

## Item 4: what separates the two families

The refused case and the accepted case have the SAME shape from the
core's side: a family that declares `se()`, has no dpar named `sigma`,
and carries one more free dpar. What differs is inside the density.

* second SCALE: the density reads `aterms[["se"]]` and never reads the
  dpar. The direction is flat, and every standard error is `NaN`.
* genuine SHAPE (a skew, a tail index): the density reads BOTH. The
  model is identified and the dpar must stay free.

Nothing the core can see tells them apart. There is no per-dpar
metadata to infer from: `dpars` is a character vector, and the link is
no help, because a positive shape and a scale take the same links. The
lpdf is an opaque closure. So the family has to say, and the new
argument is how:

* `se_dpar = "tau"` names the dpar the known standard error replaces.
  The core maps THAT one out, exactly as it maps out `sigma`.
* `se_dpar = NA` says it replaces none, because the known standard
  error is the whole scale. Every remaining dpar is then the family's
  business, which is the trust student's `nu` has always had.
* omitted keeps the convention: `sigma` if the family has one, and the
  refusal otherwise.

### The two families, measured

Both are in `tests/testthat/test-custom-family.R`. `scale_tau` is a
gaussian whose scale dpar is called `tau` and whose density ignores
`tau` when `se()` is present. `skew_se` is a skew normal whose scale IS
the known standard error and whose `alpha` the density reads
(400 rows, `alpha` drawn at 3).

| model | before | after |
| --- | --- | --- |
| `scale_tau`, `y \| se(sev) ~ x` | refused | refused (message now names `se_dpar`) |
| `scale_tau`, `bf(y \| se(sev) ~ x, tau ~ 1)` | FITTED, all 3 SEs `NaN` | refused |
| `scale_tau` with `se_dpar = "tau"` | not expressible | fits, logLik -62.4977718, SEs 0.0338019 / 0.0371666 |
| `skew_se`, `y \| se(sev) ~ x` | refused | fits with `se_dpar = NA` |
| `skew_se` undeclared | refused | refused |

The accepted case in full: logLik -285.4588975, `alpha` 2.955556 with
a standard error of 0.351892, which is 0.126 standard errors from the
3 the data were drawn at, and every standard error finite. The same
family with `alpha` pinned at 0 (the symmetric special case) reaches
-365.5154179, so the shape buys 80.06 log units: it is a parameter
with data behind it, not a flat direction.

The refused case in full, taken through the item-5 hole so the fit
could be seen at all (`bf(y | se(sev) ~ x, tau ~ 1)`): logLik
-62.49777177 and coefficients 1.009836 and 0.692258, which are the
reference model's to every digit printed, and standard errors
`NaN NaN NaN`. `diagnose()`'s flat-direction detector names
`tau_(Intercept)`: "zero gradient and an empty Hessian row".

`se_dpar = "tau"` reaches exactly the model the "name the scale
`sigma`" escape reaches: same logLik to 1e-10, same coefficients to
1e-8, and `predict(dpar = "tau")` reports 0, as `predict(dpar =
"sigma")` does for the built-ins.

## Item 5: the dpar formula, and what it did on the built-in gaussian

Verified first, as asked. On 120 rows, `bf(y | se(sev) ~ x)` against
`bf(y | se(sev) ~ x, sigma ~ 1)` with the BUILT-IN gaussian:

| | reference | with `sigma ~ 1` |
| --- | --- | --- |
| logLik | -62.49777177 | -62.49777177 |
| (Intercept) | 1.009836 | 1.009836 |
| x | 0.692258 | 0.692258 |
| standard errors | 0.0338019, 0.0371666 | `NaN NaN NaN` |
| `frm()` itself | silent | silent |
| `vcov()` / `summary()` | silent | warns |

So the claim in `feature-gaps.md` holds: this is pre-existing and not a
0.54.0 regression. It is also not entirely silent. The fit says
nothing, but asking for the covariance raises the flat-direction
warning, which names `sigma_(Intercept)` and says "moving it does not
change the likelihood at all". The extra parameter changes no fitted
value: the density never reads it, because `resid_sd()` returns
`aterms[["se"]]` alone.

`+ lf(sigma ~ 1)` produced the identical fit, so the hole was not
specific to `bf()`.

### The fix

The guard used to count a dpar as "spoken for" if it had a formula OR a
constant. A formula is not a pin: it estimates. The guard now counts
only a CONSTANT, and reads the formulas from every route that reaches
the frame: `bf()`, `lf()`, `nlf()` bodies, and a family's own
`default_forms` (merged further down `parse_one_response()`, so the
guard reads them where they are declared).

The three documented escapes were rerun after the change:

| escape | result |
| --- | --- |
| no formula (the core maps the scale out) | fits, SEs finite, `sigma(fit)` 0 |
| pin it: `bf(y \| se(sev) ~ x, sigma = 1)` | fits, logLik equal to the reference to 1e-10 |
| `se(sev, sigma = TRUE)` | fits, `sigma` estimated, and it takes a `sigma ~ 1` formula, as it always did |
| name the scale `sigma` (custom family) | fits, SEs finite |

The `default_forms` route was measured directly, because nothing in
core ships one for a scale: the same family with and without
`fam[["default_forms"]] <- list(tau = ~1)` fits at -62.497772 without
it and is refused with it, and pinning the dpar (`tau = 1`) beats the
default form and fits at -62.497772 again. There is a test.

The refusal is per response, measured on a bivariate model: with
`mvbf(bf(y1 | se(s1) ~ x), bf(y2 ~ x, sigma ~ x))` the second
response's own `sigma ~ x` is untouched and the model fits
(-217.870747), while moving the formula onto the `se()` response
refuses. An ordinary distributional model alongside a meta-analytic
one is not caught.

Back compatibility, measured: a family object with the `se_dpar` field
REMOVED, which is the 0.54.0 shape and the shape of any family inside
a fit saved before this change, reads the convention and behaves as it
did. Same logLik (-62.49777177), `sigma()` 0, `predict(dpar =
"sigma")` 0, standard errors finite, and the new formula refusal still
fires on it.

### The false-alarm question

The new refusal fires on a model whose extra parameter cannot enter the
likelihood by construction: with `se()` and no `sigma = TRUE` the
density is handed `aterms[["se"]]` and the replaced dpar is unread. The
measurement above is the evidence, not the argument.

Population checks:

* 16 core test files write a `se()` addition term (plus one helper).
  None writes it together with a dpar formula for the scale, so the new
  refusal fires in none of them. All 16 were rerun one per process;
  `R CMD check` reruns all 132 test files.
* No family in this repository or in the five extensions declares
  `se()` without also having `sigma`: only gaussian and student declare
  it at all. The breaking change therefore hits zero shipped families.

It is NOT a zero false-alarm change for an out-of-tree family. An
undeclared family whose extra dpar is a shape, given a formula
(`bf(y | se(sev) ~ x, alpha ~ 1)`), fitted correctly before and is
refused now: measured, that model produced logLik -285.4588975 with
finite standard errors, the same model `se_dpar = NA` now reaches by
declaration. The trade is deliberate and matches the rule 0.54.0 set
for the term itself: "did not say" cannot mean yes when the difference
between a correct model and a silently unidentified one is invisible to
the core. The refusal names the one-argument fix.

## Item 6: the once-per-session notice

There was no one-shot machinery to follow. `R/zzz.R` had only
`.onLoad()`; no `getOption("frmtmb...")`, no session-state environment
for messages, no `packageStartupMessage`. Searched for all three. So
the mechanism is new, and it is the smallest one that meets the brief:
a package-level environment keyed by tag, plus `notify_once()`.

The divergence, measured again on this lane rather than quoted: 200
poisson draws at `lambda = 4`, right censored at 6, 44 censored rows.
At the inclusive fit's `mu` of 3.9839046 the two readings are
-366.2873978 (frmtmb, `P(Y >= k)`) and -395.6308736 (brms,
`P(Y > k)`), a difference of 29.34 log units, and each reading's own
maximum puts `mu` at 3.9839046 against 4.1850586, a ratio of 1.0505.
The 0.54.0 NEWS reports 20.8 log units on its own draw; the size
depends on how many rows are censored, the sign and the cause do not.

Where it fires: `assemble_frame()`, after every censoring guard, when
the family type is `discrete` AND some `cens()` code is 1 (right) or 2
(interval). Placed last so a call that is about to be refused does not
spend the session's one notice; there is a test for that.

Measured, one fresh session per row:

| model | notices |
| --- | --- |
| poisson, right censored | 1 |
| the same fit again, same session | 0 |
| poisson, interval censored | 1 |
| poisson, LEFT censored | 0 |
| gaussian, right censored | 0 |
| poisson, `trunc(lb = 2)` alone | 0 |
| after those three silent fits, poisson right censored | 1 |

The last row matters: the silent cases do not spend the notice, so a
model that would diverge still gets it later in the session.

It fires from frame assembly, so `get_prior()` on such a model gets it
too and the fit that follows is then silent (measured: 1 message then
0). That is the intended reading of "once per session": the user is
told the first time the convention applies to something they asked
for, whichever entry point that was.

Suppression: it is a `message()` and nothing else, so
`suppressMessages()` silences one call, and
`options(frmtmb.notices = FALSE)` silences every notice for the
session. Both measured. A suppressed notice is spent, not saved.

`frmtmb.notices` is the package's first user-facing option: `grep` for
`getOption(` in `R/` returned only the marginaleffects class list
before this change. It is documented in the notice text itself and in
the migration vignette. There is no options help topic to add it to.

Test-suite discipline: `tests/testthat/setup.R` sets the option off for
the package's own suite, so no assertion depends on which test file
censors a count first. `testthat::source_test_setup()` is confirmed to
run it (option `unset` before, `FALSE` after). The tests that exercise
the notice turn the option back on and clear the state environment,
with a `withr::defer()` to clear it again. `test-message-uniqueness.R`
parses `R/` rather than running fits, so the notice cannot reach it;
its 6 assertions pass, and the new message templates are unique.

## Tests, seen failing first

Run against the UNFIXED 0.54.0 build installed in the lane library,
with the new test files from this worktree.

| file | before | after |
| --- | --- | --- |
| `test-custom-family.R` | 94 pass, 6 fail, 6 error | 126 pass, 0 fail |
| `test-cens-trunc.R` | 53 pass, 0 fail, 5 error | 69 pass, 0 fail |

The 5 errors in `test-cens-trunc.R` are the 5 new notice tests: nothing
answers to `frmtmb:::frmtmb_notice_state` in the unfixed build. In
`test-custom-family.R`, "a formula for the replaced scale is refused,
every route" failed 3 assertions before the fix, one per route on the
BUILT-IN gaussian, because no error was raised at all; "the same family
without the declaration is still refused" failed 3, because the old
message names three ways out and not `se_dpar`. The remaining errors
are the new `se_dpar` argument, which the old build rejects as unused.

## Test files rerun, one per R process

`NOT_CRAN=true`, the package installed in the lane's private library.
20 files, 0 failures, 0 errors.

| file | pass | skip |
| --- | --- | --- |
| test-message-uniqueness.R | 6 | 0 |
| test-custom-family.R | 126 | 0 |
| test-cens-trunc.R | 69 | 0 |
| test-cens-lccdf.R | 27 | 0 |
| test-parse.R | 43 | 0 |
| test-compat.R | 573 | 0 |
| test-families.R | 216 | 0 |
| test-autocor.R | 190 | 0 |
| test-diagnostics-ux.R | 116 | 0 |
| test-simulate-density.R | 435 | 0 |
| test-brms-agreement.R | 168 | 2 |
| test-brms-likelihood.R | 20 | 32 |
| test-brms-methods.R | 18 | 45 |
| test-bcm-data-analysis.R | 20 | 6 |
| test-bcm-esp.R | 5 | 2 |
| test-case-studies.R | 30 | 0 |
| test-gp-multidim.R | 42 | 0 |
| test-review-fixes.R | 17 | 0 |
| test-review-v25.R | 60 | 0 |
| test-v14.R | 81 | 0 |

Those 20 include every core test file that writes a `se()` addition
term (16 of them), the two censoring files, `test-message-uniqueness.R`
and `test-families.R`. None writes `se()` together with a formula for
the scale, which is why the new refusal fires in none.

## R CMD check

`R CMD build` then `R CMD check --as-cran` on the tarball, with
TinyTeX and pandoc on PATH so the manual stages actually run.

Run twice. `Status: 1 NOTE` both times, 58 stages, manual PDF 577074
bytes. The NOTE is `checking HTML version of manual ... NOTE /
Skipping checking math rendering: package 'V8' unavailable`, which is
the environmental one this project expects. Every other stage is OK,
examples included (35s, and 49s with `--run-donttest`).

The suite inside the second check, which is the tree as delivered,
reported `FAIL 0 | WARN 1 | SKIP 250 | PASS 5881`. The first check ran
before the `default_forms` test was written and reported PASS 5878,
which is the same run plus that test's three assertions. The one
warning and the 250 skips are the profile the previous lane's checks
recorded on this same tree (`FAIL 0 | WARN 1 | SKIP 250 | PASS 5398`,
and 5377 before that), so the warning is pre-existing and not this
lane's.

## What I decided NOT to do

* **`sigma()` still means the dpar named `sigma`.** A family that
  declares `se_dpar = "tau"` gets `sigma(fit) == 1`, the value any
  family without a `sigma` dpar gets. Generalizing `sigma()` to report
  a declared scale is a separate question about what `sigma()` means,
  and it would change the answer for families that have nothing to do
  with `se()`. `predict(dpar = "tau")` does report the mapped-out 0,
  because that path is about `se()` and nothing else.
* **No numerical flat-direction test at parse time.** The package
  already detects the flat direction after the fit and names the
  parameter, so a fit-time probe (evaluate the objective at two values
  of the suspect dpar) would work and needs no family API at all. It is
  the wrong trade here: it moves a call-site refusal to post-fit, it
  costs objective evaluations on every `se()` model, and a shape whose
  contribution vanishes at its starting value would be a false alarm.
  Recorded because it is the fallback if `se_dpar` proves too much to
  ask of family authors.
* **No inference of scale versus shape from links or names.** A
  positive shape and a scale admit the same links, so any rule would
  guess. Refused to guess.

## Defect found and not fixed

`sigma.frmtmb_fit()` returns 1 for a family with no dpar called
`sigma`, including one whose scale is genuinely known and mapped out
through `se_dpar`. For such a fit the residual scale is per-observation
and 0 is the answer `se()` models get elsewhere (`sigma(fit)` is 0 for
gaussian with `se()`). It is a reporting inconsistency between two
families that fit the same model under different dpar names, not a
likelihood defect, and it predates this lane for every family without a
`sigma`. Left for whoever settles what `sigma()` reports for a family
that does not have one.
