# Lane drop1raw: refits read the data the fit holds

Base: 47cf659d (frmtmb 0.71.0), measured against `rellib-r9`. Lane
library `C:/Users/adf44/source/r/wt-drop1raw-lib`. Every log named
here is in `dev/drop1raw-out/` (local, gitignored); every number below
is pasted from it.

## The defect

`drop1()` built each refit from `fit$frame$data_frame`, the MODEL
frame. It holds `log(x2)` and not `x2`, so a refit that keeps a term
which reads a raw variable stopped with "The model uses `x2`, which is
not a column of `data`". The same class of defect, a post-fit path that
reads the model frame or the workspace where it needs the raw data,
reached seven more paths (item 1).

## 1. Survey of post-fit paths on rellib-r9

Scripts `dev/drop1raw-survey.R` (seed 5, n = 200; logs
`survey-base.txt`, `survey-lane.txt`), `dev/drop1raw-fresh.R` (a fit
saved by the survey and read in a new `Rscript` process with no `d`;
`fresh-<arm>.txt`), `dev/drop1raw-battery.R` (21 model types, seed 11;
`battery-<arm>.txt`), `dev/drop1raw-probe-args.R`,
`dev/drop1raw-probe-more.R`, `dev/drop1raw-probe-dt.R`,
`dev/drop1raw-probe-subset.R`.

Failing on rellib-r9 (all pass on the lane):

| path | case | rellib-r9 |
|---|---|---|
| `drop1()` | `log()`, `scale()`, `I()`, `ns()`, `poly() + x`, `x * log(x2)`, transformed response, `weights()` + `log()`, poisson, `(1 \| g)` + `log()`, `scope = ~ x` | "`x2` (or `x`, `z`) is not a column of `data`" |
| `drop1()` | `z \| cens(cc, z2) ~ x + x2`, no transform at all | "cens(<code>, z2) reads `z2`, which is not a column of `data`" |
| `drop1()` | NA rows in another variable | same error; once refits run, the bases must be frozen to match `lm` (item 3) |
| `drop1()` | `family = fam`, a local of a function whose formula came from outside | "object 'fam' not found" |
| `drop1()` | a prior on the dropped coefficient | "Prior target not found (class=b, coef=x3)" |
| `step()` | any step that removes a term, even `z ~ x + x3` | "object 'change' not found" (update() left step()'s text delta `"~ . - x3"` unevaluated) |
| `step()` | transformed term | the `drop1()` error |
| `update()` without data | data removed, or fresh session after `readRDS()` | "object 'd' not found" |
| `update()` without data | `x2` changed in `d` after the fit | logLik -3597.2830163650 for a fit of -3536.9112541110 |
| `update()` without data | 10 rows deleted from `d` | nobs 1990 for a fit of 2000 |
| `influence()` default data | `log(x2)` + `(1 \| g)` | "all 8 deletion refits failed" |
| `influence(data = d)` | `poly(x2, 2)`, `scale(x)` | runs, but each deletion recomputes the basis: differs from lme4 by 1.30 and 0.326 of the largest displacement |
| `predict(re_formula = )` | the fit's own term `(1 + log(x2) \| g)` | "matches no group-level term of this fit" |
| `make_conditions(fit, "x2")` | `log(x2)` in the model | "`x2` is not a column of the data. It has: z, x, log(x2), g" (brms reads its stored data, which holds `x2`) |
| `frmtmb.latent::hmm_starts()` | data removed | "the model's own call could not be re-evaluated from here" |
| `frm()` on a tibble | any transformed term | warning "Setting row names on a tibble is deprecated." |
| `frm()` on a named list | any transformed term | `raw_vars` empty, so every refit fails |

Working on rellib-r9 and unchanged: `add1()` with data present,
`anova(refit = TRUE)` (refits the stored frame), `confint(method =
"boot")` and `frm_bootstrap()` (`refit()` replaces the response in the
stored frame), `confint(method = "profile")`, `frm_allfit()`,
`emmeans()`, `marginaleffects::avg_slopes()`, `conditional_effects()`
(these three read `frame$raw_vars` through `ce_base_frame()`),
`insight::get_data()`, `frm_multiple()` + `anova()` (each fit is given
its data), `predict()`, `simulate()`.

Not a refit, by construction: `frm_simulate()` takes a formula and data,
not a fit. `frmtmb.sample::frm_sample(fit)` retapes the stored frame;
its `kfold()`, `reloo()` and `update()` for draws refuse. frmtmb.eam,
.spline, .ode, .coupling and .learn have no path that evaluates a fit's
call or formula again (grep of `eval(`, `$call`, `update(`, `frm(` in
their `R/`).

The first version of this lane claimed `step()` fixed; the review
(`dev/reviews/2026-10-19-drop1raw.md`, finding 1) showed it moved
through the workspace data. Fixed in round 1 (section 5).

Left as they are (record):

- `add1()` with a new variable reads the call's data, as for `lm()`,
  and only while those data are still the fit's (round 1, finding 3).
  With the data removed or changed it refuses with brms's words: the
  fit cannot hold a variable it never read.
- `vcov_cluster(cluster = ~ g)` with `g` outside the model: the fit does
  not hold `g`, and it falls back to the call's data. Removed data: it
  fails, as on base.
- `update()` takes the formula, family, control, start and the flags
  REML, quadrature, importance and drop_unused_levels from the fit
  (round 1, finding 6). It still evaluates the call's other symbols
  (`na.action`, `se`, `verbose`) where `update()` is called.

## 2. The fix and the choice between keeping and recovering the data

Option (a), keep the raw data on the fit, was mostly in place: 0.71.0
already stores `frame$raw_vars`, the raw variables the formula reads
through a transform on the fit's rows, for `conditional_effects()` and
`emmeans()`. The model frame plus `raw_vars` is the refit data
(`fit_refit_data()`, new `R/refit-data.R`); it is built on call and
shares the model frame's vectors, so it costs no stored bytes.

Measured on 1e5 rows with 50 unused columns, `z ~ x + log(x2) +
splines::ns(x3, 3) + (1 | g)` (`dev/drop1raw-memory.R`, seed 3, logs
`memory-<arm>.txt`):

```
object                              in memory      rds raw     rds gzip
the data, all 55 columns              41.60 MB     41.58 MB     39.35 MB
the data, the 5 model columns          3.45 MB      3.44 MB      2.73 MB
fit$frame$data_frame                   4.98 MB      4.96 MB      4.51 MB
fit$frame$raw_vars                     1.53 MB      1.53 MB      1.12 MB   (lane; 7.63 MB in memory on base)
whole fit                             32.19 MB     43.97 MB     30.75 MB   (lane; 38.30 MB in memory on base)
whole fit without raw_vars            30.67 MB     42.45 MB     29.64 MB
fit_refit_data(fit), built on call     6.51 MB      6.49 MB      5.63 MB
refit data shares the model frame's z: TRUE
```

So the raw variables cost 1.53 MB of a 32.19 MB fit in memory (4.8%)
and 1.12 of 30.75 MB gzipped (3.6%); the 50 unused columns (38 MB) are
not kept. The lane's `raw_vars` is 6.1 MB smaller in memory than base's
because base set 1e5 character row names on it; the lane copies the
model frame's compact row names.

Option (b), recover the data from the call (`update()`'s route on
0.71.0), is wrong after the workspace changes: data removed and a
fresh session give "object 'd' not found"; a changed `x2` gives
logLik -3597.2830163650 for a fit of -3536.9112541110; 10 deleted rows
give nobs 1990 for 2000 (`memory-base.txt`). The lane gives the fit's
-3536.9112541110 and 2000 in both (`memory-lane.txt`). Only (a) is
correct after `saveRDS()`/`readRDS()` and after the workspace changes,
and it is what brms does (`update.brmsfit()` uses `object$data` when
`newdata` is missing; read from brms 2.23.0's source).

What changed:

- `frame_raw_vars()` (R/frame.R) also keeps the variables of an
  interval `cens(, y2)` bound (`raw_rhs()`), which the frame leaves out
  because they may be `NA`; handles a named list and builds its result
  column by column (no row names set on a tibble; a data.table read by
  `[[`); returns a zero-column frame, not `NULL`, when the model frame
  already holds everything, so `NULL` means "could not keep them".
- `subset_cut_data()` keeps the rows a univariate `subset()` cut, so
  `fit_refit_data(rows = "data")` gives brms's rows and `nobs()` of an
  `update()` is the fit's (160, not the 107 fitted rows;
  `probe-subset-<arm>.txt`).
- A channel into `frm()` (`refit_channel()`): the refit data carry an
  attribute `frmtmb_refit` with the call the refit records, the fit's
  `predvar_map` and the fit's `na_action`. `frm()` records that call
  (so the data, formula and family read as the original call spelled
  them), passes `predvar_map` to `assemble_frame()`, which patches the
  combined frame's `predvars` as `influence()` does, and stores the
  fit's NA rows as the refit's `na_action`. Round 1 replaced the first
  version's formula rewriting (`freeze_bases()`) and its fix-up of the
  call after `eval()`.
- `drop1()` refits with `do.call(frm, ...)` on the channelled refit
  data, the fit's family (on the `bf()` object), its prior (marked as
  `update()` marks it, so a bound on the dropped term is dropped with a
  message naming `drop1()`), control, data2, quadrature, importance and
  `drop_unused_levels`, and `se = FALSE`. The kept terms keep the
  fit's bases and names, so a submodel is nested in the fit as in
  `stats::drop1()`. It refuses a fit whose prior adds a density
  (round 1, finding 5), and a refit whose df is `NA`.
- `update()` takes from the fit, by value, what it does not replace:
  the data (`fit_refit_data(rows = "data")`, channelled), the formula
  (`fit$bform`; `frm()` keeps the unresolved family of a
  `categorical()`-type model as attribute `frmtmb_written` for
  `newdata`), the family, the control, the start and the flags. The
  call it evaluates is `frmtmb::frm(...)` with those values;
  `evaluate = FALSE` returns it, which is what `stats::step()` and
  `add1()` evaluate. A formula that reads a variable the fit does not
  hold (`formula_vars()`, less parameter names) reads the call's data
  only when `refit_data_matches()` finds the fit's rows and values in
  it; otherwise brms's refusal. A character `formula.` (step()'s) is
  read with `as.formula()`; an argument given as `NULL` is removed.
- `influence()` defaults to the refit data and passes the fit's
  `predvar_map` to `assemble_frame()` (new argument), which patches the
  combined frame's `predvars`, so the deletion refits keep the bases
  and the coefficient names.
- `re_term_cnms()` (re_formula) and `make_conditions()` read the refit
  data.
- `prior_drop_unmatched()` names the caller (`drop_who` attribute).
- frmtmb.latent `hmm_starts()` refits through `stats::update()`,
  evaluated from the caller's frame.

Rows: the refit data are the rows after `na.action` (and the subset cut
re-applied), with the fit's weights column, so a refit never gains the
rows `na.action` removed. Under `na.exclude` the refit takes the fit's
NA rows, so its residuals have the fit's 200 rows, as lm's update has
(round 1, finding 8).

## 3. Comparison with stats::drop1() and lme4

`dev/drop1raw-compare.R` (seed 5, n = 200; NA case: `x3` NA at rows 3,
17, 40), log `compare-lane.txt`. LRT for `lm` is `n log(RSS_i/RSS_0)`
from `drop1.lm`'s table, for `glm` the deviance difference, for `lmer`
lme4's `LRT` column; "max rel diff" is `max |a - b| / max(|b|, 1)`.

```
lm log                             rows TRUE | df TRUE | LRT max rel diff 3.32e-15 | LRT 74.815925 55.721306
lm scale                           rows TRUE | df TRUE | LRT max rel diff 7.14e-15 | LRT 67.636826 44.456944
lm I                               rows TRUE | df TRUE | LRT max rel diff 1.09e-13 | LRT 0.190332 32.245302
lm ns                              rows TRUE | df TRUE | LRT max rel diff 5.70e-15 | LRT 73.536988 54.870890
lm poly_x                          rows TRUE | df TRUE | LRT max rel diff 1.08e-14 | LRT 51.751735 72.184250
lm inter                           rows TRUE | df TRUE | LRT max rel diff 7.39e-14 | LRT 3.196594
lm resp                            rows TRUE | df TRUE | LRT max rel diff 1.00e-15 | LRT 70.927276 49.545978
lm scope = ~ x                     rows TRUE | df TRUE | LRT max rel diff 1.03e-14 | LRT 72.819930
lm NA in x3, ns()                  rows TRUE | df TRUE | LRT max rel diff 2.01e-13 | LRT 36.711310 2.032690
lm NA in x3, scale()               rows TRUE | df TRUE | LRT max rel diff 6.47e-13 | LRT 54.827215 0.060402
lm NA in x3, poly()                rows TRUE | df TRUE | LRT max rel diff 5.15e-14 | LRT 35.042424 1.939604
glm log                            rows TRUE | df TRUE | LRT max rel diff 1.00e-11 | LRT 48.639791 34.826189
glm scale_ns                       rows TRUE | df TRUE | LRT max rel diff 2.70e-12 | LRT 46.154433 37.278694
glm inter                          rows TRUE | df TRUE | LRT max rel diff 1.01e-10 | LRT 0.002227
glm NA in x3, ns()                 rows TRUE | df TRUE | LRT max rel diff 2.92e-11 | LRT 39.783550 1.740210
lmer log                           rows TRUE | df TRUE | LRT max rel diff 6.95e-12 | LRT 119.136068 125.840711
lmer scale_ns                      rows TRUE | df TRUE | LRT max rel diff 1.22e-10 | LRT 116.948201 121.199094
lmer inter                         rows TRUE | df TRUE | LRT max rel diff 1.24e-11 | LRT 17.425976
lmer poly_x                        rows TRUE | df TRUE | LRT max rel diff 4.09e-11 | LRT 106.647559 110.683551
lmer NA in x3, ns()                REF ERROR number of rows in use has changed: remove missing values?
brms has a drop1 method: FALSE
```

- The NA rows are where freezing matters: the survey's `z ~
  splines::ns(x2, 3) + x3` with 3 NA in `x3` gives AIC `<none>`
  752.8936170441, `ns` 804.4436392253, `x3` 751.0199384715 on the lane
  and the same 10 decimals from `drop1(lm)` plus the constant
  (`survey-lane.txt`).
- lme4's `drop1()` refits through `update()`, which re-evaluates the
  call on the full data, so dropping `x3` gains the 3 rows and lme4
  stops ("number of rows in use has changed"). frmtmb gives the `lm`
  answer on the fit's rows.
- lme4's own fragility, met while writing the comparison: called with a
  function argument named `data`, lme4's refits resolved `data` in the
  formula's environment to `utils::data()` and failed with "first
  argument must be a named list"; the script names the data by a
  global symbol.
- brms 2.23.0 has no `drop1()` method (checked in the namespace and by
  `getS3method()`), so there is nothing to match there. Its
  `update()`, which the lane follows, reads `object$data`.
- On base the comparison script cannot reach the transformed-term
  error: `ft()` passes `family = fam`, a local, and every base refit
  stops first with "object 'fam' not found" (`compare-base.txt`), the
  defect of the table above.

influence() against lme4's `influence.merMod()` (`dev/drop1raw-influence.R`,
logs `influence-<arm>.txt`), max |frmtmb - lme4| over the
deleted-group fixed effects, as a ratio to the largest displacement:

```
base  z ~ poly(x2, 2) + x + (1 | g)   default data  ERROR all 8 deletion refits failed
base  z ~ poly(x2, 2) + x + (1 | g)   data = d      ratio 1.30e+00
base  z ~ scale(x) + log(x2) + (1 | g) data = d      ratio 3.26e-01
lane  z ~ poly(x2, 2) + x + (1 | g)   both          ratio 4.62e-06
lane  z ~ scale(x) + log(x2) + (1 | g) both          ratio 8.34e-06
```

## 4. Tests (first round; section 5 supersedes the counts)

`tests/testthat/test-refit-data.R` (new) and
`extensions/frmtmb.latent/tests/testthat/test-hmm-starts.R` (one test
replaced, one added). Run against both builds with
`dev/drop1raw-runtest.R`.

`test-refit-data.R`, 16 tests (logs `t-refit-<arm>.txt`):

```
lane      RESULT test-refit-data.R pass=71 fail=0 error=0 skip=0 warn=0 tests=16
rellib-r9 RESULT test-refit-data.R pass=3 fail=5 error=15 skip=0 warn=1 tests=16
```

Every one of the 16 tests fails on rellib-r9, each for its own defect:

| test | rellib-r9 |
|---|---|
| drop1() on `log()`, `scale()`, `I()`, `ns()`, `poly() + x`, `x * log(x2)`, transformed response, against lm() | "`x2` is not a column of `data`" |
| drop1(scope = ) | same |
| drop1() poisson against glm() | same |
| drop1() mixed against lme4 | same |
| drop1() with NA rows, frozen `ns()`, `scale()`, `poly()` against lm() | same |
| drop1() with interval `cens(cc, z2)` | "reads `z2`, which is not a column of `data`" |
| drop1() with a local `family = fam` and a prior on the dropped coef | "object 'fam' not found" |
| a fit whose data were removed (`local()`): drop1() and update() | "`x2` is not a column of `data`" |
| update() after `x2` changed and 10 rows deleted | nobs 190 for 200, logLik -371.2 for -388.8 |
| update() of a univariate `subset()` fit with `d` shadowed by NULL | "`data` is NULL" |
| saveRDS(), readRDS() in a new Rscript process | the child stops on the drop1() error |
| step() removing a term | "object 'change' not found" |
| influence() against lme4, `data = d` and default | ratio 1.3009 >= 1e-4; "all 8 deletion refits failed" |
| `re_formula = ~ (1 + log(x2) \| g)` and `make_conditions(fit, "x2")` | "matches no group-level term" |
| a named list | "`x2` is not a column of `data`" |
| a tibble | warning "Setting row names on a tibble is deprecated." |

The `subset()` test also guards this lane's own change: the lane's
first version, before `subset_cut_data()`, gave `nobs(update(fit))`
107 for the fit's 160 (`dev/drop1raw-probe-subset.R`, seen in a run of
the probe, not of the test), and the test asserts 160. That version was
not rerun against the test.

frmtmb.latent `test-hmm-starts.R` (logs `t-hmmstarts-<arm>.txt`):

```
lane      RESULT test-hmm-starts.R pass=88 fail=0 error=0 skip=0 warn=0 tests=19
rellib-r9 RESULT test-hmm-starts.R pass=85 fail=0 error=1 skip=0 warn=0 tests=19
```

The rellib-r9 error is the new "a fit whose data the caller cannot see
refits from its own" ("could not be re-evaluated from here"). The
replaced test pinned that refusal; its replacement keeps the refusal for
a call symbol other than the data (`control = ctl_gone`), and passes on
both builds, as a guard must.

Every suite, one file per process on the lane library (driver
`dev/drop1raw-par.sh`, runner `dev/drop1raw-runtest.R`, silent reporter
so no count is capped; log `all2.log`, the run after the last code
change; every per-file log's `lib:` line is the lane library):

```
ran 383 of 383
frmtmb files=237 pass=17511 fail=0 error=0 skip=185 warn=0
frmtmb.coupling files=11 pass=542 fail=0 error=0 skip=5 warn=0
frmtmb.eam files=30 pass=1808 fail=0 error=0 skip=3 warn=0
frmtmb.latent files=10 pass=364 fail=0 error=0 skip=2 warn=0
frmtmb.learn files=15 pass=430 fail=0 error=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 error=0 skip=1 warn=0
frmtmb.sample files=53 pass=2697 fail=0 error=0 skip=4 warn=0
frmtmb.spline files=16 pass=594 fail=0 error=0 skip=1 warn=0
```

Gated (`FRMTMB_BRMS_FIT_TESTS=true`), the two gated files that call
`update()`, `drop1()`, `influence()`, `make_conditions()` or `step()`
(logs `gated-*-<arm>.txt`), the same on both builds, so no ported-suite
row flips there:

```
test-brms-suite-methods.R         lane pass=162 fail=0 error=0 skip=0 | rellib-r9 pass=162 fail=0 error=0 skip=0
test-brms-suite-brmsfit-helpers.R lane pass=3 fail=0 error=0 skip=0   | rellib-r9 pass=3 fail=0 error=0 skip=0
```

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
`OPENBLAS_NUM_THREADS=4`; logs `ob-*.txt`): `test-refit-data.R` pass=71
fail=0 error=0 skip=0; `test-hmm-starts.R` pass=88 fail=0 error=0;
`dev/drop1raw-compare.R` all 19 reference rows match in rows and df,
LRT max rel diff 2.0e-16 to 1.22e-10 (the lm rows move in the last
digits, 4.34e-15 against 3.32e-15 for `lm log`, so the BLAS did
change). The tests' tolerances are ratios to the statistic itself
(1e-6 for drop1, 1e-4 of the largest displacement for influence, 1e-8
across processes); the bitwise `identical()` checks compare two fits in
one process.

Core's whole suite in one session (`dev/drop1raw-onesession.R`, log
`onesession.txt`):

```
ONESESSION files=227 pass=17506 fail=9 error=0 skip=185 warn=0
FAILED: test-compat-register.R :: expects = accepts a name another package already supplies
FAILED: test-compat-register.R :: a declared but unresolved rule side is reported, not dropped
FAILED: test-compat.R :: every rule wins at least one pair
FAILED: test-compat.R :: rule patterns name real kinds, groups, and features
FAILED: test-compat.R :: declared post-fit methods are callable
```

All nine failed expectations are in `test-compat.R` and
`test-compat-register.R`, the one-session failures known on 47cf659d
(fixed by lane testiso, not in this base). Both files pass in their own
processes above. No test of `test-refit-data.R` failed in the one
session.

`R CMD check --as-cran`, built and checked in `dev/drop1raw-check/`
(`R_LIBS` = lane, rellib-r9, user library; pandoc and TinyTeX on PATH;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`):

```
frmtmb         Status: 1 NOTE   (HTML manual: math rendering skipped, V8 unavailable; tests 13m OK)
frmtmb.latent  Status: 1 NOTE   (the same V8 note; tests 14s OK)
```

## 5. Fix round 1 (review `dev/reviews/2026-10-19-drop1raw.md`)

Every new test was run against the lane build as installed before the
round (logs `r1-t-refit-oldlane.txt`, `r1-t-hmmstarts-oldlane.txt`) and
against rellib-r9 (`r1-t-*-base.txt`), then against the fixed build
(`r1-t-*-lane.txt`). The reviewer's own scripts `dev/drop1raw-rev-*.R`
were rerun on the fixed build, read-only, with their output in my logs
`rev-r1-<script>-lane.txt`.

```
test-refit-data.R, 24 tests
  rellib-r9          pass=9   fail=10 error=21 skip=0 warn=1   (all 24 tests fail)
  lane before round  pass=82  fail=15 error=3  skip=0 warn=1   (9 tests fail)
  lane after round   pass=111 fail=0  error=0  skip=0 warn=0
test-hmm-starts.R (frmtmb.latent), 20 tests
  rellib-r9          pass=87 fail=0 error=1   (the data-gone test)
  lane before round  pass=88 fail=0 error=1   (the start test, finding 7)
  lane after round   pass=90 fail=0 error=0
```

The 9 tests that fail on the lane as it was before the round, with the
observed failure:

| finding | test | lane before round |
|---|---|---|
| 1 | step() moves on the fit's data | data changed: logLik -451.5 for -427.5; data removed: "object 'd' not found" |
| 1 | a saved fit refits in a new session (now with step()) | child: "object 'd' not found" in `eval.parent` |
| 2 | update() and drop1() keep the fit's bases when rows are dropped | `update(fit)` 1.76e-5 relative from the fit; drop1 AIC 1.76e-5 from `AIC(update(fit, ~ . - x3))`; `poly()` and `scale()` coefficients 0.0086 and 0.0118 apart |
| 4 | drop1() keeps a kept term's name, so its prior stays | the bound on `b_scalex` dropped from both refits |
| 3 | update() refits on the fit's data after the object changes | no refusal; the changed `d` became the model |
| 3 | a new variable is read from data that still hold the fit's | `update()` and `add1()` refused nothing after one value changed |
| 6 | update() takes the formula, family and control from the fit | "object 'form' not found" |
| 5 | drop1() and anova() refuse a fit whose prior adds a density | no refusal from `drop1()`, `step()` or `anova()` |
| 8 | update() on an na.exclude fit pads as the fit does | 197 rows for lm's 200 |

The poisson `step()` test (`stats::step()` on `glm()`, the same path
and deviance differences to 1e-6) passed before the round too; it is
the parity check the coordinator asked for, not a pin.

Per finding, on the fixed build:

1. **step()** (HIGH). `update(evaluate = FALSE)` returns the call it
   evaluates: `frmtmb::frm()` with the fit's data (channelled, with the
   original call as the spelling `frm()` records), formula, family and
   control by value. The reviewer's `dev/drop1raw-rev-step2.R`
   (seed 21, n = 300): data in place `z ~ x + log(x2)`, logLik
   -413.02051460; data gone `z ~ x + log(x2)`; after `d$x2` changed
   `z ~ x + log(x2)`, logLik -413.02051460 (was -435.26018233); `add1`
   with the data gone: ok. `dev/drop1raw-rev-step.R`: lm path and AIC
   differences unchanged (-1.14e-13, 1.02e-12), glm path the same
   (4.55e-13, 8.50e-10, 1.45e-09). With NA rows in `x4`, lm's
   `step()` stops ("number of rows in use has changed") and frmtmb's
   stays on the fit's 297 rows. My survey (`survey-lane.txt`,
   `fresh-lane.txt`): `step()` and `add1(~ . + I(x3^2))` now also run
   with the data removed and in a new session.
2. **Bases on update()** (MEDIUM). The refit's `assemble_frame()` gets
   the fit's `predvar_map` through the channel, in `update()` and in
   `drop1()`; `freeze_bases()` is gone. `dev/drop1raw-rev-probe.R`
   section 2: `ns(x2, 3) + x3` fit -361.02064484, `update(fit)`
   -361.02064484 (diff 0.00e+00; was 1.49e-2); drop1 AIC(-x3)
   733.73080917 and `AIC(update(fit, ~ . - x3))` 733.73080917;
   `update(family = student())` -361.02064489, as `frm()` gives. The
   same holds for `poly()`, `scale()` and `ns() * g`. brms 2.23.0
   computes the basis again here (the reviewer's `brms2.txt`) and is
   wrong, so this does not follow it. The comparison with lm, glm and
   lme4 (`compare-lane.txt`) is unchanged to the last printed digit.
3. **New variable** (MEDIUM). brms 2.23.0, run
   (`dev/drop1raw-brms-newvar.R`, `brms-newvar.txt`): "New variables
   found: 'w'\nPlease supply your data again via argument 'newdata'."
   (two variables: "'w', 'v'"). frmtmb gives that text, then the
   reason, when the call's data are gone or no longer hold the fit's
   rows and values (`refit_data_matches()`); otherwise it reads them,
   so `add1()` works. The call's data are looked up where `update()` is
   called and then in the formula's environment, where `add1()`
   evaluates. Probe section 1: after `d` changed, `update(fit, ~ . +
   w)` and the global `wv` are refused (before: nobs 188 and logLik
   -284.1072416771, and an error naming the wrong cause).
4. **Coefficient names** (MEDIUM). Follows from 2: a bound on
   `b_scalex` is dropped only from the refit without `scale(x)`.
5. **MAP fits** (MEDIUM, the coordinator's decision). `drop1()`, and
   so `step()`, and `anova()` refuse a fit whose resolved prior entries
   add a density (`fit_prior_density()`). A bound-only prior has no
   entries and is not refused: its logLik equals the ML fit's to 1e-8,
   and `drop1()` and `anova()` run. Before the round the reviewer's
   probe printed an "LRT" of -1.562 with p = 1. No existing test in any
   suite called either on a MAP fit (all suites pass).
6. **Formula, family, control** (MEDIUM). `dev/drop1raw-rev-locals.R`:
   `update(fit)` -271.13497122 (was "object 'form' not found");
   `update(fit, ~ . - x3)` -271.27068236; with a global `form = z ~ x`
   and `fam`, `update(fit)` keeps `z ~ x + log(x2) + x3` and the fit's
   -271.13497122 (was `z ~ x`, -336.49411568). The start and the flags
   REML, quadrature, importance and drop_unused_levels come from the
   fit too. A `categorical()`-type family is kept unresolved on the fit
   (attribute `frmtmb_written` of `fit$bform`) for `newdata` and
   new-variable refits, so it reads its categories again from those
   data; a refit on the fit's own rows keeps the resolved one.
7. **hmm_starts() start** (LOW). The probe passes `start = NULL`;
   `update()` removes an argument given as `NULL` (assigning `NULL` to
   a name a call lacks is an out-of-bounds error). The reviewer's
   `dev/drop1raw-rev-hmmstart.R`: "ok, status
   original,converged,converged". The refusal test now uses
   `na.action = <gone local>`, which `update()` does not take from the
   fit.
8. **na.exclude** (LOW, done). The channel carries the fit's
   `na_action`, classed by the refit's own `na.action`; probe: residual
   lengths 800 and 800 (4 columns of 200), lm 200.
9. **NEWS and findings**: rewritten; each defect gives the release it
   dates from (found with `git log -S` on the function or line).
10. **Style**: "center" in every new or changed line (the remaining
    "centre" in `NEWS.md`, `R/confint.R` and other files predate the
    lane); the long `update()` roxygen sentence is now three paragraphs
    of short sentences.
11. **catmulti**: section 6.

Other changes this round: `prior_drop_unmatched()` labels a bounds-only
specification "bounds" (it printed an empty string); `drop1()` stops
with `df_why` when a refit's df is `NA` (for the catmulti merge).

Every suite, one file per process on the lane library after the round
(`all3.log`; every per-file `lib:` line is the lane library):

```
ran 383 of 383
frmtmb files=237 pass=17551 fail=0 error=0 skip=185 warn=0
frmtmb.coupling files=11 pass=542 fail=0 error=0 skip=5 warn=0
frmtmb.eam files=30 pass=1808 fail=0 error=0 skip=3 warn=0
frmtmb.latent files=10 pass=366 fail=0 error=0 skip=2 warn=0
frmtmb.learn files=15 pass=430 fail=0 error=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 error=0 skip=1 warn=0
frmtmb.sample files=53 pass=2697 fail=0 error=0 skip=4 warn=0
frmtmb.spline files=16 pass=594 fail=0 error=0 skip=1 warn=0
```

Core in one session (`onesession.txt`): `ONESESSION files=227
pass=17546 fail=9 error=0 skip=185 warn=0`; the 9 failed expectations
are the 5 known tests of `test-compat.R` and `test-compat-register.R`
on 47cf659d, nothing else.

Gated (`FRMTMB_BRMS_FIT_TESTS=true`, `r1-gated-*-lane.txt`):
`test-brms-suite-methods.R` pass=162 fail=0 error=0 skip=0;
`test-brms-suite-brmsfit-helpers.R` pass=3. These are the counts
rellib-r9 gave in the first round, so no ported-suite row flips there.

OpenBLAS 0.3.26 with LAPACK (`OPENBLAS_NUM_THREADS=4`, logs
`r1-ob-*.txt`): `test-refit-data.R` pass=111, `test-update-pool.R`
pass=22, `test-confint-anova.R` pass=35, `test-methods-audit.R`
pass=58, `test-hmm-starts.R` pass=90, all fail=0 error=0; the
comparison script matches in all 19 rows, largest LRT rel diff
1.22e-10. The OpenBLAS R copy was deleted afterwards.

`R CMD check --as-cran` (`dev/drop1raw-check/<pkg>/check.txt`):
frmtmb `Status: 1 NOTE`, frmtmb.latent `Status: 1 NOTE`; both NOTEs
are the HTML manual's skipped math rendering (V8 unavailable).

## 5b. Fix round 2 (re-check appended to the review)

The re-check confirmed findings 1 to 11 fixed and added N1 to N3.

- **N1, the channel outlived its call (MEDIUM).** `refit_channel()`
  now stores the formula and family values it was made for and
  `NROW(data)`; `refit_channel_for()` in `frm()` uses the channel only
  when all three are `identical()`, and otherwise strips the
  attribute, records `frm()`'s own `match.call()` and takes no
  `predvar_map` or `na_action`. `update()` evaluates a formula or a
  given family that is still an expression to its value first, so the
  values can be compared; `drop1()` makes one channel per refit. The
  reviewer's `dev/drop1raw-rev-r1.R` section A on the fixed build
  (`rev-r2-r1-lane.txt`): `frm(bf(z ~ x), data = dd)` records
  `frm(formula = bf(z ~ x), data = dd)` and has no NA rows (before:
  the update's call and rows 3, 17, 40). `dev/drop1raw-rev-r1e.R`: the
  edited call fits `z ~ x` and records its own call (formula by value,
  as `match.call()` of that call holds it), and on other data gives
  -335.3828985832, the value `frm()` gives on that data. Before the
  fix the recorded call named `z ~ splines::ns(x2, 3) + x + x3`.
- **N2, a missing column passed the data check (MEDIUM).** The frame
  now records `data_vars`, the variables of the model the data
  supplied by name (taken before the `intercept` fill), and
  `refit_data_matches()` returns `FALSE` when one of them is not in
  the candidate data. This is stricter than `intersect(names(rd),
  formula_vars(fit$bform))`: a variable the fit read from the formula
  environment is not required in the data. Section B: "x2 column
  missing" and "x2 column missing, global x2 vector" are now refused
  with brms's text (before: a wrong-cause error and a silent logLik
  -406.3372068809 for -380.7151487902). The other shadowing cases are
  unchanged.
- **N3, filed (see "To file").**
- **Docs.** `?update.frmtmb_fit` says that the `evaluate = FALSE` call
  holds the data by value (printing it prints them; 34600 characters
  deparsed at n = 200 in `rev-r2-r1-lane.txt`), and that an edited
  call records its own call and takes nothing from the old fit.
- **Console message.** The one message `test-refit-data.R` let through
  (`drop1(fb)`, bounds on `b_x3`) is now asserted with
  `expect_message()`.

Tests (`r2-t-refit-*.txt`), `test-refit-data.R`, now 27 tests:

```
lane before round 2  pass=116 fail=6 error=0 skip=0 warn=0 tests=27
lane after round 2   pass=122 fail=0 error=0 skip=0 warn=0 tests=27
OpenBLAS 0.3.26      pass=122 fail=0 error=0 skip=0 warn=0 tests=27
```

The three new tests fail on the round-1 build: "an edited update()
call records the model it fits" (recorded call identical to the old
fit's; `na.action()` the 3 old rows), "data taken from an update()
call fit as plain data" (logLik -420.8280 for -420.8205 from the old
frozen knots; the old NA rows; the call names `d`, not `dd`), and
"data that lack a variable the fit read are not its data" (no
refusal).

Every suite after round 2, one file per process on the lane library
(`all4.log`; this covers every file that calls `update()`, `drop1()`,
`step()` or `anova()`, since `frm()` itself changed):

```
ran 383 of 383
frmtmb files=237 pass=17562 fail=0 error=0 skip=185 warn=0
frmtmb.coupling files=11 pass=542 fail=0 error=0 skip=5 warn=0
frmtmb.eam files=30 pass=1808 fail=0 error=0 skip=3 warn=0
frmtmb.latent files=10 pass=366 fail=0 error=0 skip=2 warn=0
frmtmb.learn files=15 pass=430 fail=0 error=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 error=0 skip=1 warn=0
frmtmb.sample files=53 pass=2697 fail=0 error=0 skip=4 warn=0
frmtmb.spline files=16 pass=594 fail=0 error=0 skip=1 warn=0
```

Core in one session (`onesession.txt`): `ONESESSION files=227
pass=17557 fail=9 error=0 skip=185 warn=0`, the 9 being the 5 known
`test-compat.R` and `test-compat-register.R` tests of 47cf659d.
`R CMD check --as-cran` of frmtmb (`dev/drop1raw-check/core`):
`Status: 1 NOTE`, the V8 math-rendering note. frmtmb.latent did not
change in round 2; its round-1 check stands.

## To file

- **N3 (LOW, not a regression).** After a delta update,
  `getCall(update(fit, ~ . - x3))` holds the deparsed `bf()` object
  (`frm(formula = structure(list(formula = z ~ x + (1 | g), ...`, 2474
  characters), on rellib-r9 and on the lane alike. A `family =` update
  records a readable call. The recorded call still names the user's
  data, and `update()` of an updated fit reproduces it. Fix: record the
  delta through `bform_call()`, as the complete-formula branch does.

## Files touched

- `R/refit-data.R` (new): `fit_refit_data()`, `formula_vars()`,
  `refit_channel()`, `refit_channel_for()`, `fit_prior_density()`,
  `prior_lrt_refusal()`, `refit_data_matches()`, `data_label()`.
- `R/fit.R`: `frm()` reads the refit channel (call, `predvar_map`,
  `na_action`) and keeps the unresolved family (`frmtmb_written`).
- `R/frame.R`: `frame_raw_vars()`, `raw_rhs()`, `subset_cut_data()`,
  `refit_na_rows()`, `assemble_frame(predvar_map = , refit_na_action
  = )`, the frame's `subset_cut` and `data_vars` elements.
- `R/confint.R`: `anova.frmtmb_fit()` (MAP refusal, documentation),
  `drop1.frmtmb_fit()` (refit, MAP refusal, df check, documentation),
  `update.frmtmb_fit()` (rewritten after the formula branches, and its
  documentation).
- `R/influence.R`, `R/re-formula.R` (`re_term_cnms()`),
  `R/brms-utilities.R` (`make_conditions()`), `R/priors.R`
  (`prior_drop_unmatched()` message).
- `man/anova.frmtmb_fit.Rd`, `man/drop1.frmtmb_fit.Rd`,
  `man/update.frmtmb_fit.Rd`, `man/influence.frmtmb_fit.Rd`,
  `man/make_conditions.Rd` (roxygen).
- `tests/testthat/test-refit-data.R` (new), `NEWS.md`.
- `extensions/frmtmb.latent/R/hmm-starts.R`,
  `extensions/frmtmb.latent/tests/testthat/test-hmm-starts.R`,
  `extensions/frmtmb.latent/NEWS.md`.
- `dev/drop1raw-*.R`, `dev/drop1raw-*.sh`, this file.

Named overlaps with other lanes: `R/fit.R`, `R/frame.R`,
`R/confint.R`, `R/priors.R`, `NEWS.md`.

## 6. Overlap with lane catmulti

The reviewer repeated the three-way `git merge-file` after converting
all three versions to LF. A merge of the CRLF working copies against
git's LF base reports 6 false conflicts in `R/priors.R`, so normalize
the line endings first. The real conflicts:

- `R/confint.R`: 1, the `drop1()` roxygen; both lanes add a paragraph
  after "not part of the deletion scope." Keep both.
- `man/drop1.frmtmb_fit.Rd`: 1, generated; roxygenise after the merge.
- `NEWS.md`: 1, two `# frmtmb (development version)` headings; join
  them.
- `R/priors.R` merges cleanly. The `drop1()` code merges cleanly:
  catmulti's `dp_main` scope and `df0` NA refusal sit above this lane's
  refit loop and stay intact.
- `R/fit.R`: since round 1 both lanes edit `frm()`. This lane adds the
  channel read after `match.call()`, the `frmtmb_written` lines after
  `resolve_deferred_families()`, and two arguments to the
  `assemble_frame()` call; catmulti adds one statement after
  `carry_finalized_responses()`. They are separate hunks; check that
  part of the merge by hand.

How they compose: each `drop1()` refit is `frm(formula = <bform minus
the term>, data = <channelled fit_refit_data()>)`, and catmulti's
`refcat_first_outer(spec, bform, data, ...)` in `frm()` assembles the
reference-category twin from that `data`, which holds `x2`, so
`log(x2)` evaluates. For the merge:

- This lane's refit loop already stops on a refit whose df is `NA`,
  with `attr(ll_i, "df_why")`; until catmulti's `loglik_df()` lands,
  that branch cannot fire. Add a merged test that reaches it.
- catmulti's comments in `R/fit.R` (at the `refcat_first_outer()`
  call) and in `loglik_df()` say the fit keeps only its model frame,
  where `log(x)` has no `x`. After the merge that is false: change
  them, or compute the twin lazily from `fit_refit_data()`. A twin
  built from the fit's rows should also take the fit's `predvar_map`
  (`assemble_frame(predvar_map = )`), or its spline bases differ when
  `na.action` removed rows.

After the merge, run `test-refit-data.R` with catmulti's tests, and
`drop1()` on a `categorical(refcat = NA)` fit with `log(x2)`, which
neither lane alone runs.

## Defects found and not fixed

- `vcov_cluster()` on a variable outside the model reads the call's
  data and fails when they are gone.
- `update()` still evaluates `na.action`, `se` and `verbose` of the
  stored call where it is called.
- A smooth (`s(x)`) in an `influence()` deletion refit is rebuilt from
  the subset's rows (mgcv's knots move); `predvars` do not cover it.
- `AIC()` of a MAP fit stays the penalized quantity, as documented.
- An optional cleanup the reviewer named, not done: `ce_base_frame()`
  could call `fit_refit_data()`.

## Bumps and floors

- frmtmb: minor bump. Behavior changes marked BREAKING: `update()`
  takes the data, formula, family and control from the fit; a new
  variable is read only from unchanged data; `drop1()`, `step()` and
  `anova()` refuse a fit whose prior adds a density; `drop1()` and
  `influence()` keep the fit's bases.
- frmtmb.latent: patch bump, and its `frmtmb (>=)` floor must move to
  the frmtmb release carrying this change: the data-gone
  `hmm_starts()` test fails against 0.71.0's `update()`.
- No other extension changes; all seven suites were run against the
  lane's core.
