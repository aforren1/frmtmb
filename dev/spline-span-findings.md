# SPLINE-SPAN lane findings

Worktree `C:/Users/adf44/source/r/frmtmb-wt-spline-span`, branch
`wt-spline-span`, base 68d6782 (frmtmb 0.52.0, frmtmb.spline 0.2.0).
Private library
`scratchpad/ss-lib`: frmtmb 0.52.0 (worktree), frmtmb.spline 0.2.0,
flexsurv 2.3.2, gratia 0.11.2, brokenstick 2.7.0. R 4.6.1.

The task is the residual the spline-core review left
(`dev/reviews/2026-09-05-spline-core.md:1454`): `predict(newdata =)`
warns past a `ps()` knot span and `frm_lp_basis()` does not, so the
extension's three curve functions, the doors a user actually draws a
curve through, are silent.

## 1. The fire count before any change

Model `bf(y ~ lev + ps(t, k = 10, pad = 0.3), lev ~ 1, nl = TRUE)`,
n = 200, knot span [-0.297464, 1.291776]; a 12-row grid wholly past the
span.

| door | span warnings |
|---|---:|
| `predict(newdata = past-span)` | **1** |
| `predict(newdata = inside)` | 0 |
| `frm_lp_basis(newdata = past-span)` | **0** |
| `frm_lp_basis(newdata = inside)` | 0 |
| `frm_lp_basis(newdata = NULL)` | 0 |
| `frm_curve(newdata = past-span)` | **0** |
| `frm_curve_deriv(newdata = past-span)` | **0** |
| `frm_curve_feature(newdata = past-span)` | **11** |
| `frm_curve_feature(newdata = inside)` | 0 |

Two separate defects, not one.

`frm_curve()` and `frm_curve_deriv()` are silent because their only
route to the curve is `frm_lp_basis()`, which does not arm the check.

`frm_curve_feature()` is the opposite: it already fires **11 times on
one call**, because its root scan and Newton refinement reach the curve
through `sp_predict_eta()`, which is `predict(newdata =)`, and that door
IS armed. One warning per grid scan and one per Newton iteration per
root. So the extension today is both silent and noisy, in different
functions, about the same fact.

## 2. What naive arming gives, measured

Arming means `check = TRUE` on `lp_basis_nl()`'s two `ps_env()` calls
(`R/predict.R:3239`, `:3259` at the base). Measured by giving `ps_env()`'s `check`
argument a default of `TRUE` in the installed namespace:

| shape | door | armed fires |
|---|---|---:|
| one `ps()` term | `frm_lp_basis(newdata =)` | 1 |
| **two** `ps()` terms in one body | `frm_lp_basis(newdata =)` | 2 (one per term, both unique) |
| a nonlinear parameter with its own linear predictor | `frm_lp_basis(newdata =)` | 1 |
| one `ps()` term | `frm_curve()` | 1 |
| one `ps()` term | `frm_curve_deriv(order = 1)` | 1, counting 36 of 36 (the 3-point stencil over a 12-row grid) |
| one `ps()` term | `frm_curve_feature()` | **12** (the 11 above plus the basis) |

So naive arming does not multiply per row or per coefficient: `A` is
12 x 10 and the count stays at one per term. The reason is that
`lp_basis_nl()` tapes the body with `RTMB::MakeTape()`, which evaluates
the R body ONCE; `tp(chat)` and `tp$jacobian(chat)` replay the tape and
never re-enter the closure. That is an accident of RTMB's
implementation rather than a property the code states, which is why the
fix below collects rather than warns from inside the closure.

## 3. The defect naive arming exposes: `is.numeric()` is not a tape guard

`ps_span_warning()` (`R/ps.R:554`) guards with
`if (!is.numeric(x) || !length(x))`, and the review records the claim
that this "is what keeps it off the tape". It does not.

```
MakeTape(function(x) ...) :  is.numeric TRUE   is.double FALSE
                             class "advector"  inherits(x, "advector") TRUE
```

`is.numeric()` is **TRUE** for an RTMB advector. The guard has never
been exercised on a tape because the only two armed call sites were
`predict()`'s nonlinear branch, which evaluates off the tape. Arming
`frm_lp_basis()` puts it on one, and on the model whose `ps()` argument
names a nonlinear parameter,

```
bf(y ~ lev + ps(t + sh, k = 8, pad = 0.2), lev ~ 1, sh ~ 1, nl = TRUE)
```

naive arming turns `frm_lp_basis(newdata =)` into

```
Error: Comparison is generally unsafe for AD types
```

`predict(newdata =)` fires once on the same fit and the same grid, so
the two doors would disagree by an error rather than by a warning. The
guard has to be `inherits(x, "advector")`, which is the predicate core
already uses at `R/importance.R:665`.

## 4. What landed

### Core

* **`R/ps.R:286`** - the closure hands `chk` to the check instead of
  testing it as a flag: `if (!isFALSE(chk)) ps_span_warning(x, term, chk)`.
* **`R/ps.R:572`** - `ps_span_warning(x, pt, sink = TRUE)`. The guard is
  now `inherits(x, "advector") || !is.numeric(x) || !length(x)`, which is
  the defect in section 3. The message is built once into `msg`; an
  environment `sink` collects it keyed on the term's label, anything else
  warns where it stands.
* **`R/ps.R:598`** - `ps_span_flush(sink)`, one warning per collected term.
* **`R/ps.R:617`** - `ps_span_signal(out, n, pt)`, the single raise site. The
  warning now carries the class **`frmtmb_ps_span_warning`**, so a
  consumer catches this one warning without matching on its text.
  `frmtmb.spline` does exactly that.
* **`R/predict.R:3239`** - `lp_basis_nl()` creates the collector,
  `span <- if (is.null(newdata)) FALSE else new.env(parent = emptyenv())`,
  matching `predict()`'s `check = !is.null(newdata)` arming.
* **`R/predict.R:3257`, `:3277`** - both `ps_env()` calls take it.
* **`R/predict.R:3287`** - one evaluation of the body OFF the tape before
  `MakeTape()`, wrapped in `tryCatch(suppressWarnings(...))`, so a `ps()`
  argument that names a nonlinear parameter is checked on numbers rather
  than on advectors. Errors are swallowed: a diagnostic must not take the
  basis with it, which is the discipline `fit_end_checks()` already
  adopted.
* **`R/predict.R:3292`** - `ps_span_flush(span)` after the Jacobian.

### The extension

* **`R/curve-cov.R:39`** - `sp_catch_span(expr)`, catching by CLASS.
* **`R/curve-cov.R:60`** - `sp_span_stop(span)`, one refusal template.
* **`R/curve-cov.R:104`** - `sp_curve_parts()` catches core's warning and
  returns it as `parts$span` on both of its two exits.
* **`R/curve.R:173`** and **`R/curve-deriv.R:126`** - one warning per
  `ps()` term per call, under the name of the function the user called.
* **`R/curve-feature.R:155`** and **`:213`** - refuses, at the grid scan
  before any root is refined and again on the difference stencil at the
  located roots. The Newton loop between them is muffled, because it is
  clamped inside a bracket the scan cleared.

### The tests that assert the property

Both packages' `test-message-uniqueness.R` now pool `errorCondition()`
with `stop()` and `warningCondition()` with `warning()`, skipping the
`class =` and `call =` arguments by name so a class string is not read as
message text. Without that the new idiom would have left its templates
unasserted, which is the coverage the old walker had by accident.

## 5. The fire count after

| door | before | after |
|---|---:|---:|
| `predict(newdata = past-span)` | 1 | 1 |
| `frm_lp_basis(newdata = past-span)` | **0** | **1** |
| `frm_lp_basis(newdata = inside)` | 0 | 0 |
| `frm_lp_basis(newdata = NULL)` | 0 | 0 |
| `frm_lp_basis`, two `ps()` terms | 0 | **2**, one per term |
| `frm_lp_basis`, `ps(t + sh)` (the tape case) | 0 | **1** (naive arming: ERROR) |
| `fitted()`, `simulate()` | 0 | 0 |
| `frm_curve(past-span)` | **0** | **1** |
| `frm_curve_deriv(past-span)` | **0** | **1** |
| `frm_curve_feature(past-span)` | **11 warnings** | **refuses, 0 warnings** |
| `frm_curve*(inside)` | 0 | 0 |

`predict()` and `frm_lp_basis()` now produce **byte-identical** warning
text on the same fit and grid, checked at `re.form = NA` and `NULL`, on
an in-span grid that the nonlinear shift pushes out (2 of 9) and on a
wholly out-of-span one (9 of 9).

Values did not move: `max |frm_lp_basis()$eta - predict()|` is **0** on
both the plain and the taped-argument model, and on a non-`ps()` smooth
`frm_curve()`'s covariance check is unchanged at 2.22e-16 relative with
`n_predict` still 1.

## 6. The extension's three doors, after

| door | grid inside the span | grid past it |
|---|---|---|
| `frm_curve()` | silent | **1 warning**, naming `frm_curve()`, carrying the span and counting the grid (15 of 15) |
| `frm_curve_deriv()` | silent | **1 warning**, naming `frm_curve_deriv()`, and saying that its count (45 of 45) is over the three-point difference stencil |
| `frm_curve_feature()` | silent, and still finds the peak of `sin(2 pi t)` at 0.25 | **REFUSES**, by name, with **0 warnings** |

The feature refusal replaces the eleven warnings rather than adding a
twelfth: measured at 0 raised warnings on the refusing call. It fires at
the grid scan, before any root is refined, and again on the difference
stencil at the located roots; both go through one message template, so
the package's message-uniqueness property holds. A grid that STRADDLES
the span is refused too.

A model with no `ps()` term is untouched, including on a grid far
outside the DATA range, which for an ordinary `s()` smooth is
extrapolation and not a partition-of-unity cliff: `cov_rel_error` stays
below 1e-10 and `n_predict` stays 1.

## 7. The residual documentation (brief item 3)

The brief names "the vignette's 'What a built-in family would add'
paragraph in frmtmb.spline". **There is no such paragraph in
frmtmb.spline.** That heading exists only once in the repository, at
`vignettes/reinforcement-learning.Rmd:453`, which is core's and the
protocol lane's. Grepped for it and for its wording across
`extensions/frmtmb.spline/`: absent.

What IS there is the same staleness under different headings, and that
is what was corrected:

* **`R/frmtmb.spline-package.R`**, section "What this package reads that
  frmtmb does not promise", said in the PRESENT tense that
  `fit$cache$Vjoint` "has no precedent and is an internal" that
  `frm_curve()` reads. That has been false since 0.2.0, which the
  package's own NEWS says in its first bullet. Replaced.
* The same file's "Two limits that are core's to fix" named the fix as
  "an `lccdf` slot, or a post-fit family hook". Both landed in 0.52.0
  and this package declares both. Replaced with what core has since
  supplied and, in the protocol lane's names, what is still missing: a
  **per-row log-likelihood slot** (`loo()`, `waic()`) and a
  **per-group** one (`frm(importance =)`).
* **`vignettes/royston-parmar.Rmd`**, the "Why it matters" paragraph,
  said core "gives a family no complementary log-CDF slot to hand back
  `log S` directly" and described the floor as current. Rewritten to the
  past tense with what replaced it.
* The same vignette's "The refusal is POST-FIT and cannot be otherwise"
  paragraph said "the family protocol has no hook that runs when a fit
  finishes" and that `frm_curve()` "will not draw a curve off a fit like
  that one". Both false: `post$fit_check` exists and this family
  declares it, and `frm_curve()` draws on that fit (measured, and
  already pinned by `test-rp-floored.R`'s "the curve functions draw a
  deeply censored fit now"). Rewritten, and the `floored-bad` chunk
  loses its `error = TRUE` because nothing errors there any more.

The RL vignette was left alone, as instructed.

## 8. `rp_floored()` on the review's 600-subject fit (brief item 4)

```
convergence   1   ("false convergence (8)", as the review recorded)
logLik        -575.537942891
n_censored_deep    1
max_nlogS          55.7302136
threshold          19.2
n_nonmonotone      0
scale              "hazard"
n_obs              600
rp_floored(action = "error")   did NOT refuse
fitted()                       REFUSED, "no mean on the response scale"
predict(type = "response")     REFUSED, same reason
frm_curve()                    drew a curve
```

`n_nonmonotone` is 0 and the monotonicity refusal does not fire, which
is the brief's expectation. `fitted()` refuses by design, and so does
`predict(type = "response")`; that refusal is about a survival time
having no mean on the response scale and has nothing to do with the
floors, so an exact likelihood does not change it.

**One row survives the count, and here is why.** The field counts
censored rows whose fitted `-log S` passes 19.2. There is one, the
deliberately extreme censored subject at `t = 50`. It is NOT floored:

```
-log S (= H)                 55.7302136
its exact term, -H           -55.73021360299
S = exp(-H)                  6.26146e-25
the old form, log(1 - (1-S)) -Inf
the 0.51.0 floor             -35.127363
```

The row is scored exactly, and the whole fit reproduces an independent
log-likelihood written from the hazard-scale definition
(`eta + log(deta/dlog t) - log t - exp(eta)` for an event, `-exp(eta)`
for a censor) to **4.206e-12**, which is the review's own figure to
every digit it printed.

So the count is right and the NAME was wrong. `n_censored_floored`
reported a row that the `lccdf` slot had made exact, which is the same
defect class as the help page that told users the importance correction
applies. Renamed to **`n_censored_deep`**, with the docs, three test
assertions and the vignette moved with it. Nothing about what is
counted, or about the threshold, changed.

## 9. What the off-tape pass costs

`frm_lp_basis()` gains one evaluation of the nonlinear body off the
tape, and only on a body that HAS a `ps()` term at non-NULL `newdata`;
`has_ps` gates it, so every other nonlinear model is byte-for-byte the
old path with `check = FALSE`.

Measured on `bf(y ~ lev + ps(t, k = 12, pad = 0.3), lev ~ 1 + (1 | g))`,
n = 600, 40 groups, five blocks of 60 calls, minimum reported. (The
first two attempts at this measurement were wrong and are worth
recording: `force(promise)` caches after the first call, so a timing
loop over a promise times one call and 59 no-ops. The figures below use
a zero-argument function.)

| call | ms |
|---|---:|
| `frm_lp_basis(newdata = NULL)`, no pass | 85.17 |
| `frm_lp_basis(newdata = )`, 600 rows | 93.00 |
| difference | 7.83, and it is mostly not the pass: it also carries the newdata column extraction and the design rebuild |
| `predict(newdata = )`, which IS one body pass plus a link inverse plus argument checks | **1.00** |

So the pass is bounded above by 1.00 ms against a 93 ms call, about
**1%**. The call is dominated by the joint-precision solve, as
`?frm_curve`'s cost section already says. Six sibling lanes were
running R processes throughout, so treat the absolute figures as an
upper bound and the ratio as the result.

## 10. Runs

### The frmtmb.spline suite, one file per process, audited by name

```
expected 8   ran 8   unique 8   missing 0   extra 0
```

| file | passed | failed | error | warning | skipped |
|---|---:|---:|---:|---:|---:|
| test-curve.R | 62 | 0 | 0 | 0 | 0 |
| test-deriv.R | 38 | 0 | 0 | 0 | 0 |
| test-gratia.R | 12 | 0 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 | 0 |
| test-royston-parmar.R | 67 | 0 | 0 | 0 | 0 |
| test-rp-floored.R | 42 | 0 | 0 | 0 | 0 |
| **test-span.R (new)** | **29** | 0 | 0 | 0 | 0 |
| test-surface.R | 48 | 0 | 0 | 0 | 0 |
| **total** | **302** | **0** | **0** | **0** | **0** |

302 against the 267 the brief records at 0.2.0. The 35 new assertions
are the 29 in `test-span.R` and 6 in `test-rp-floored.R`'s new
`fitted()` / field-name test; nothing was removed.

### The core files the contract audits by name

| file | passed | failed | error |
|---|---:|---:|---:|
| test-ps.R | **71** (51 before) | 0 | 0 |
| test-predict-lp-basis.R | 30 | 0 | 0 |
| test-predict-newdata.R | 12 | 0 | 0 |
| test-message-uniqueness.R | 6 | 0 | 0 |

### Message-template coverage, measured on the source

The walker in both `test-message-uniqueness.R` files now pools
`errorCondition()` with `stop()` and `warningCondition()` with
`warning()`, skips `class =` and `call =` by name, and descends one
`paste0()` into a constructor's message argument. Run over the sources:

| package | error templates | warning templates | message templates | duplicates |
|---|---:|---:|---:|---:|
| frmtmb | 646 | 28 | 14 | **0** |
| frmtmb.spline | 29 | 3 | 0 | **0** |

Core's warning pool is 28 where the unextended walker found 27: the
`ps()` span template is the one it could not see, and it is asserted
now. The extension's warning pool is 3 where the old walker found 1
(`royston_parmar()`'s fit-end notice); the two new curve notices are
the difference. Without those two extensions to the walker the new
idiom would have been unasserted, which is the coverage the old walker
had only by accident.

## 11. What was left out, and why

* **The extension's version floor was not raised.** `frmtmb.spline`
  declares `Depends: frmtmb (>= 0.52.0)`, and the condition class its
  three curve functions now catch is new in core's "(development
  version)". Against core 0.52.0 exactly, `sp_catch_span()` catches
  nothing: the curve functions go quiet again and
  `frm_curve_feature()` stops refusing. That is a degradation, not a
  break, and the numbers stay right. The floor has to rise to whatever
  version the consolidation gives core, and this lane cannot know that
  number: the spline-core review records that lane as the only one
  bumping `DESCRIPTION`, so guessing here would collide.
  **A consolidation action, not an omission.**

* **`frm_curve_deriv()`'s row count is over the difference stencil, not
  the grid.** Rebuilding the span check on the grid alone would cost a
  second `frm_lp_basis()` call, which is the expensive part of the
  whole function. The message says which it is counting instead.

* **`frm_curve_feature()` refuses on a grid that only PARTLY leaves the
  span**, even when the peak it would find is well inside. That is
  deliberate and it is the brief's instruction, but it is a judgment
  call worth recording: the alternative, refusing only when a LOCATED
  root lies outside, would let the grid scan run against the decaying
  partial sum first, and a sign change manufactured out there can put a
  spurious root inside the span. The message names the remedy ("Narrow
  newdata to the span").

* **`simulate(newdata = )` is still silently ignored.** The
  spline-core review filed this as out of scope and pre-existing
  (`simulate.frmtmb_fit` has no `newdata` argument, so it lands in
  `...`), and it is not this lane's.

* **`frm_lp_basis()`'s LINEAR branch is not armed.** It does not need
  to be: `ps()` is refused outside a nonlinear body, so a `ps()` term
  can only ever reach `lp_basis_nl()`.

* **The core reinforcement-learning vignette was not touched**, as
  instructed. Its "What a built-in family would add" section is the
  protocol lane's.

## 12. A harness failure of mine, and what it cost

The first attempt at the verification runs produced two "failures" that
were not failures, and the diagnosis is worth recording because the
symptom looked exactly like a real regression.

I ran three heavy jobs at once: the core `R CMD check`, the
`frmtmb.spline` `R CMD check`, and the 117-file core sweep. All three
used `--library=<ss-lib>` or read from it. `R CMD check --library=L`
INSTALLS the package under test into `L`, so the core check reinstalled
`frmtmb` into the very library the other two were reading from. During
that window a reader fell through to the user library:

```
Error: package 'frmtmb' 0.50.0 was found, but >= 0.52.0 is
required by 'frmtmb.spline'
```

That one line explains all of it:

| symptom | truth |
|---|---|
| `frmtmb.spline` check: 3 ERRORs, 3 WARNINGs, 1 NOTE, including four "Missing link(s)" for `frmtmb::ps()` and `frmtmb::frm_lp_basis()` | it was resolving links and running examples against frmtmb **0.50.0**, which has neither |
| sweep: `test-custom-family.R` failed=5 error=3 | passes standalone: **41 passed, 0 failed, 0 error** |
| sweep: `test-importance.R` NOTRUN, then failed=3 error=16 | passes standalone: **185 passed, 0 failed, 0 error** |

A second, independent harness bug inflated the same sweep: the aborted
first sweep's `sh` loop survived the stop and kept appending to the
output path I had recreated, so the tally read 127 result lines for 117
files with 18 duplicates and the last 8 files missing. The rewritten
sweep writes one result FILE per test file, which cannot race with
itself, and the check script no longer passes `--library` at all so a
check installs into its own `.Rcheck` directory.

Sibling lanes were not involved and were not harmed: `p2`, `sm` and
`cj` each run their own `R CMD check` against their own `p2-lib`,
`sm-lib` and `cj-lib`, and all five sibling lanes were still running
after the cleanup.

**The core `R CMD check` itself was not affected** and its result
stands: it installs `frmtmb` into `ss-lib` at the start and then reads
what it wrote, and no other process wrote `frmtmb` there. Status **OK**,
zero ERROR, WARNING or NOTE lines anywhere in the log; examples [65s],
`--run-donttest` [119s], vignette rebuild [551s]; inner suite
`FAIL 0 | WARN 1 | SKIP 205 | PASS 4224`.

### The base, installed separately, settles it

`git archive 68d6782` into `scratchpad/ss-base-src`, installed into
`scratchpad/ss-base-lib`, base sources AND base tests, brokenstick
reachable so nothing skips for want of a Suggests package:

| file | base 68d6782 | this lane |
|---|---|---|
| test-prior-compat.R | 0 failed, 0 error, **1 warning**, 161 passed | identical: 0/0/**1**/161 |
| test-importance.R | 0 failed, 0 error, 185 passed | identical |
| test-custom-family.R | 0 failed, 0 error, 41 passed | identical |
| test-message-uniqueness.R | 6 passed | 6 passed, with two more template pools asserted |
| test-ps.R | **51** passed | **71** passed |
| test-predict-lp-basis.R | 30 passed | 30 passed |

The `WARN 1` inside the core `R CMD check` is `test-prior-compat.R`'s,
in "coef and group narrow the classes that read them". It is present on
the untouched base, it is in `R/priors.R` territory (the wt-priors2
lane's), and nothing this lane touched can reach it. Not mine, and not
new.

test-ps.R's 51 is the spline-core review's own post-punch figure to the
digit, which also says the base install is the tree the review measured.

## 13. Merge notes for the consolidator

* **`R/predict.R`**: all five hunks are inside `lp_basis_nl()`, lines
  3226 to 3287. `lp_delta_A()`, which wt-car-jacobian owns, begins at
  line 1485, and the ranef condSD path is elsewhere again. Closest
  approach is 1679 lines.
* **`R/ps.R`**: this lane's alone.
* **`tests/testthat/test-message-uniqueness.R`** is the one core file
  outside `R/` that another lane might also want, because it is the
  file every lane's new messages have to satisfy. This lane widened its
  walker (two new template pools plus a `paste0()` descent) rather than
  changing what it asserts, so a sibling that only ADDS templates merges
  cleanly; a sibling that also rewrites `collect()` will conflict there
  and the resolution is to keep both widenings.
* **No sibling-owned file is touched.** `R/structure.R`, `R/families.R`,
  `R/objective.R`, `R/importance.R`, `R/compat.R`, `R/frame.R`,
  `R/priors.R` and `vignettes/reinforcement-learning.Rmd` are all
  absent from the changed list.
* **`DESCRIPTION` is not touched** in either package. See section 11 on
  the extension's version floor, which the consolidation has to set.

## 14. Touched files

Modified (21): `NEWS.md`, `R/predict.R`, `R/ps.R`,
`tests/testthat/test-message-uniqueness.R`, `tests/testthat/test-ps.R`,
`extensions/frmtmb.spline/NEWS.md`,
`extensions/frmtmb.spline/R/{curve.R, curve-cov.R, curve-deriv.R,
curve-feature.R, frmtmb.spline-package.R, rp-check.R}`,
`extensions/frmtmb.spline/man/{frm_curve.Rd, frm_curve_deriv.Rd,
frm_curve_feature.Rd, frmtmb.spline-package.Rd, rp_floored.Rd}`,
`extensions/frmtmb.spline/tests/testthat/{test-message-uniqueness.R,
test-royston-parmar.R, test-rp-floored.R}`,
`extensions/frmtmb.spline/vignettes/royston-parmar.Rmd`.

Added (2): `dev/spline-span-findings.md`,
`extensions/frmtmb.spline/tests/testthat/test-span.R`.

Nothing committed. `man/` is roxygen output only, and roxygenise is
idempotent on both packages; core's `man/` and both `NAMESPACE` files
are unchanged, because every core roxygen edit was to an `@noRd` block.

## 15. Final runs, after the machine restart

The machine restarted at about 22:35 on 2026-09-05 and killed every
process. `scratchpad/ss-lib` survived and was re-validated before it was
trusted: frmtmb 0.52.0 and frmtmb.spline 0.2.0 load, `ps_span_flush()`
and `ps_span_signal()` are present, and every fire count in section 5
reproduces (`predict(out)` 1, `frm_lp_basis(out)` 1,
`frm_lp_basis(in)` 0, `frm_curve()` 1, `frm_curve_deriv()` 1,
`frm_curve_feature()` refuses).

The core sweep had 95 of 117 results on disk. One,
`test-smooth-population.R`, had been killed mid-run and carried no
RESULT line; it was deleted and re-run rather than trusted. The
remaining 23 ran one process at a time.

### The full core suite, one file per process, audited by name

```
expected 117   ran 117   unique 117   missing 0   extra 0   duplicates 0
failed 0 | error 0 | warning 1 | skipped 91 | passed 6565
```

The single warning is `test-prior-compat.R`'s, in "coef and group narrow
the classes that read them". Section "The base, installed separately"
shows it firing identically on the untouched 68d6782.

| file the contract audits | passed | failed | error |
|---|---:|---:|---:|
| test-ps.R | **71** (base: 51) | 0 | 0 |
| test-predict-lp-basis.R | 30 (base: 30) | 0 | 0 |
| test-predict-newdata.R | 12 | 0 | 0 |
| test-message-uniqueness.R | 6 | 0 | 0 |
| test-compat.R | 266 | 0 | 0 |
| test-cens-lccdf.R | 27 | 0 | 0 |
| test-nl-body-vars.R | 23 | 0 | 0 |
| test-tmb-examples.R | 28 | 0 | 0 |
| test-bracket-access.R | 8 | 0 | 0 |

### as-cran

`R CMD build` then `R CMD check --as-cran --no-manual`, with
`_R_CHECK_CRAN_INCOMING_=false` and pandoc from
`RSTUDIO_PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"`
on PATH. No `--library`, so each check installs into its own `.Rcheck`.

| package | Status | ERROR | WARNING | NOTE | inner suite |
|---|---|---:|---:|---:|---|
| frmtmb 0.52.0 | **OK** | 0 | 0 | 0 | FAIL 0, WARN 1, SKIP 205, PASS 4224 |
| frmtmb.spline 0.2.0 | **OK** | 0 | 0 | 0 | FAIL 0, WARN 0, SKIP 11, PASS 222 |

The strings `ERROR`, `WARNING` and `NOTE` appear **zero** times in the
`frmtmb.spline` log. Both checks were verified to have built from
sources byte-identical to the worktree's `R/` and `tests/testthat/`
(`diff -rq`, only the untracked `_snaps` directory differs), so both OK
statuses certify the tree as delivered.

Stages that cost something: core examples [65s], `--run-donttest`
[119s], vignette rebuild [551s]; spline examples OK, vignette rebuild
[32s]. The corrected `royston-parmar` vignette rebuilds and its live
chunk prints `n_censored_deep: int 1` and `max_nlogS: num 55.7`, so the
numbers in its prose are computed on the page rather than typed.

### One thing that moved underneath this lane

Main advanced past 68d6782 during the session, to b131fe1 (`612cde1`
CI, `de9d639` DESCRIPTION Suggests, `b131fe1` docs). This worktree is
still based at 68d6782 and was not rebased. `de9d639` touches
`DESCRIPTION`, which this lane does not, and the other two are CI and
`docs/`; none is a file this lane changed, so the merge is unaffected.

---

# Punch round, against dev/reviews/2026-09-05-spline-span.md

Verdict PUNCH, on the extension side; the core work verified clean and
the reviewer could not break it. Both behavior items reproduced before
they were fixed.

## Item 1. frm_curve_deriv() died on its own defaults past the span

Reproduced on `y ~ lev + ps(t, k = 10, pad = 0.1)`, span
[-0.09499071741, 1.09801176277], grid `seq(0.5, span[2] + 1, len = 25)`:

```
frm_curve_deriv(default simultaneous = TRUE)  warn 1  ERROR:
    missing values and NaNs not allowed if na.rm is FALSE
frm_curve_deriv(simultaneous = FALSE)         warn 1  ok
frm_curve(default simultaneous = TRUE)        warn 1  ok
```

with **8 of 25** derivative standard errors exactly zero against **0 of
25** for the curve, which is why only the derivative died.

**Guarded the divisor**, at `R/curve-cov.R:207`. The argument is not a
patch: `S` is positive semi-definite, so a zero diagonal entry forces
that whole row of `S` to zero, hence that row of `L` to zero, hence
`(L z)` exactly zero. The point's deviation is DETERMINISTIC, it is
covered with probability one, and it cannot be the argmax. It leaves the
maximization rather than contributing a zero over a zero:

```r
keep <- is.finite(div) & div > 0
mx <- apply(abs((L[keep, , drop = FALSE] %*% z) / div[keep]), 2L, max)
```

Dropping nothing is bit-identical to the old arithmetic, and the suite
says so: `test-curve.R` 62 and `test-deriv.R` 38 both pin simultaneous
critical values at fixed seeds and neither moved.

A grid on which EVERY row is deterministic has no band at all, so
`R/curve-cov.R:209` refuses by name instead of returning one with no
content ("every point on this grid has a standard error of exactly
zero ... Use simultaneous = FALSE, or move the grid inside the span").

Measured after the fix, on the same grid: the call returns, the
simultaneous critical value is 2.73784 against a pointwise 1.95996, the
5 zero-standard-error rows get a band of width exactly 0, and a grid
`seq(span[2] + 5, span[2] + 9)` with 12 of 12 rows deterministic
refuses.

**Tests on the DEFAULT arguments**, which the reviewer correctly said
every past-span test had dodged: "frm_curve_deriv() past the span works
on its own defaults" and "a grid with no uncertainty anywhere refuses by
name" in `test-span.R`.

## Item 2. A grid ending exactly on a knot

Reproduced, and the reviewer's sharpest point reproduced with it: the
refusal says "Narrow newdata to the span", and doing exactly that
reproduced the refusal.

| grid | before: deriv | before: feature(max) |
|---|---|---|
| exactly the knot span | 1 warning | **REFUSED** |
| pulled in by 1e-7 of the span | 1 warning | **REFUSED** |
| pulled in by 1e-6 of the span | 0 | ok |

The span was checked on the widened stencil. `frm_curve_deriv()` hands
core `c(x - e, x, x + e)` and the stationary-point scan of
`frm_curve_feature()` hands it `c(x - e1, x + e1)`, with `e` a millionth
of the grid's range (measured: 1.193e-06 at order 1, 1.193e-04 at order
2), so a grid laid on the knot is outside by `e`.

**Checked the user grid instead**, via `sp_span_on_grid()`
(`R/curve-cov.R:90`), one `predict()` on the grid that was passed. It is
consulted at `R/curve-deriv.R:137` and at `R/curve-feature.R:166` and
`:229`, and only when the widened call already reported something, so an
ordinary call pays nothing for it. That gate is exact for the
derivative, whose stencil point set CONTAINS the grid, so a clean
stencil implies a clean grid; for the feature scan it holds whenever the
spline argument is monotone in `var`, which is every shape this package
documents. Said so in the helper's comment rather than left implicit.

After the fix all four calls are clean on the exact-boundary grid and on
grids pulled in by 1e-7, 1e-6 and 1e-5 of the span, and a grid pushed
1 percent of the span OUTSIDE still warns and still refuses, so the
check was moved rather than loosened.

The `frm_curve_deriv()` warning now counts GRID rows, so the message
drops its "three rows per grid point" clause, and the test that pinned
"45 of 45" pins "15 of 15" and asserts "45 of 45" is absent.

**Test** "a grid ending exactly on a knot is inside the span" covers all
three functions plus `type = "crossing"` and the order-2 default, where
the stencil is widest.

## Item 3. The version floor

Left both `DESCRIPTION` files untouched, as instructed; core goes to
0.53.0 at consolidation and the extension floor rises with it. Added the
requirement to `extensions/frmtmb.spline/NEWS.md` as its own bullet: the
surfacing is driven by the `frmtmb_ps_span_warning` class, and against a
released frmtmb 0.52.0 the class is never raised, so the two curve
functions fall silent and `frm_curve_feature()` stops refusing, with
nothing to say so.

## The cheap items

* `for (m in ...)` renamed to `for (msg in ...)` at
  `R/curve-deriv.R:141`. The reviewer flagged only that one, where `m`
  shadows `nrow(nd)`; `R/curve.R:171` had the same spelling with nothing
  to shadow, and was renamed too rather than left as the next trap.
* The bare `expect_error()` in `test-span.R` now carries
  "the search bracket leaves", with a comment saying why a bare one was
  wrong: it would pass on an error raised before the span check ran.
* **The two comments crediting RTMB are corrected**, and the reviewer is
  right on the substance. They measured 2 closure entries per term on
  the basis route against 1 on the predict route, and the second entry
  is this lane's OWN off-tape pass, not RTMB re-entry. The collector is
  therefore REQUIRED by this design rather than defensive against an
  unstated RTMB property. `R/ps.R` (the `ps_env()` block) and
  `R/predict.R` (the `span <-` block) now say so, and point a reader
  chasing a double-fire at the off-tape pass rather than at the tape.
  Section 2 of this document overstated the same thing and the
  correction belongs there too: RTMB enters once, and the collector
  earns its place on the second entry that this lane itself added.

## Items I did not act on, and why

* **Finding 5, the collector last-write-wins.** Left, deliberately, and
  the reviewer marked it "no user-visible undercount is reachable". I
  measured the same thing: the two passes see the same `x` by
  construction, because the off-tape pass is `eta_fun(chat)` and the
  tape build is `MakeTape(eta_fun, chat)`, the same closure over the
  same `dl`. A merge would be strictly more code for a case that cannot
  arise while both passes are handed the same grid, and the collision
  the reviewer raises (two `ps()` terms sharing a label) cannot happen
  either: `pt[["label"]]` is `deparse1()` of the call, so two identical
  calls in one body are one term. Worth revisiting only if a third
  caller ever arms the collector on a different grid.
* Core `R CMD check` was not re-run. The only core edits this round are
  two comment blocks; `roxygenise()` produces no change under core
  `man/`, both `NAMESPACE` files are unchanged, and the four named core
  files pass unchanged (71 / 30 / 12 / 6).

## Punch-round runs

frmtmb.spline suite, one process per file, audited by name: 8 expected,
8 ran, 8 unique, 0 missing, 0 extra. failed 0, error 0, warning 0,
skipped 0, **passed 321**.

test-curve 62, test-deriv 38, test-gratia 12, test-message-uniqueness 4,
test-royston-parmar 67, test-rp-floored 42, **test-span 48** (was 29),
test-surface 48. 321 against 302 before the punch: +19, all in
`test-span.R`, and no count anywhere else moved.

| core file | passed |
|---|---:|
| test-ps.R | 71 |
| test-predict-lp-basis.R | 30 |
| test-predict-newdata.R | 12 |
| test-message-uniqueness.R | 6 |

`roxygenise()` idempotent on both packages, both `NAMESPACE` files
unchanged. `frmtmb.spline` as-cran **Status: OK**, zero occurrences of
`ERROR`, `WARNING` or `NOTE` anywhere in the log, vignette rebuild
[28s] OK, inner suite FAIL 0 / WARN 0 / SKIP 14 / PASS 222, built from
sources verified byte-identical to the worktree `R/` and
`tests/testthat/`.
