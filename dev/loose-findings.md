# Loose-ends lane: recorded debt, verified and cleared

Worktree `frmtmb-wt-loose`, branch `wt-loose`, at `c18253e` (two commits
past `main`'s 54d4d92; the extra pair is the `wt-api` merge). Nothing is
committed from this lane.

Each entry states whether the recorded item was still real, what the
measurement was, and what changed.

## 1. frmtmb.spline: the second span refusal (WRONG FIRST, corrected)

**I got this one wrong, a review caught it, and the change is
reverted.** Recording the error rather than quietly fixing it, because
the mistake is instructive.

### What I claimed, and why it was false

I argued the second refusal was dead: both refusals end in the same
call, `sp_span_stop(sp_span_on_grid(sp$fit, nd, ...))`, so a non-empty
result would always have fired at the first refusal and control could
never reach the second with anything to say. I checked the argument
against the suite, found no test reaching the second refusal, and read
that as confirmation.

The argument has a hole. The two calls are textually identical but ask
about **different data**. The scan and the five-point stencil are built
from `row1` replicated (`curve-feature.R:120`, `:125-129`, `:222-223`),
so every point they evaluate holds every column except `var` pinned at
row 1's value. `sp_span_on_grid()` (`R/curve-cov.R:90-92`) instead runs
one `predict()` on the **whole** `nd`. My step "a grid point outside
the span puts the scan outside it too" is true for the `var` column and
false for every other column. With a second `ps()` term the grid can
leave that term's span in a row the scan never touches, so `scan$span`
is empty, the first guard is false, and the second refusal is the only
one left.

### Measured, not argued

The reviewer's construction (`scratchpad/rvlo-reach2.R`), on default
`eps` with no exotic argument: a two-`ps()` fit over `t` and `z`, a
grid whose `t` stays inside `t`'s span but whose row 10 puts `z` past
`z`'s span, and an `at` placing the root within `e2` of the grid end so
the stencil crosses `t`'s outer knot. Run on my own builds:

| build | result |
|---|---|
| branch, block deleted | **no error, no warning**, returns a root at 1.298035 |
| branch, block restored | **ERROR**, "the search bracket leaves a ps() term's knot span" |
| control (all `z` inside), either build | no error, root returned |

So the deletion turned an error into a silent answer. The sibling
settles it independently: `frm_curve_deriv()` keeps a structurally
identical gated re-ask at `R/curve-deriv.R:136-141`, and on the same
grid it warns while the feature path had gone quiet. If my reasoning
made one dead it would make the other dead, and it demonstrably is not.

My **empirical** claim was correct and the reviewer confirmed it: over
`test-span.R` the first refusal fires 7 times, the second 0, and the
stencil guard is never even true. That is a statement about test
coverage, not reachability, and I read it as the latter. My "routine"
cost claim was unsupported too: across the whole spline suite the guard
was true zero times.

### Done

- **Restored** the block at
  `extensions/frmtmb.spline/R/curve-feature.R:235-237`. The executable
  code is now byte-identical to main's (`git diff` on the file, with
  comment lines filtered, is empty); only the comment differs, and it
  now records why the block is reachable so the next reader does not
  repeat the deletion.
- **Added the missing coverage**, `test-span.R:277-309`, a
  `sp_span_fit2()` helper plus a test that pins both halves: the
  control, where the stencil guard is true but the grid is clean and
  the re-ask must decline to refuse; and the two-`ps()` case, where the
  refusal must fire. Verified as a real gate rather than a passing
  assertion, by building a copy with the block removed
  (`scratchpad/lo-noblock-lib`) and re-running: **51 passed, 1 failed**
  at `test-span.R:306` ("Expected to throw an error") without the
  block, **52 passed, 0 failed** with it.

### Left undone, and recorded

The reviewer notes a genuine defect underneath that both main and my
branch share: in the control case the stencil really did leave the
span, `parts$span` is non-empty, and neither build says anything. Main
also refuses on evidence unrelated to the search, since the offending
row is one the search never evaluates. That incoherence is real, but
fixing it is a behavior change to a refusal that other lanes may rely
on, not debt-clearing, so it is recorded in `dev/test-backlog.md`
rather than done here.

## 2. Unqualified roxygen links that resolve to another package

**Real, and narrower than feared.** `man/frm_lp_basis.Rd:44` and `:78`
carried `\code{\link[brms:gp]{gp()}}`. The source is `R/predict.R:3085`
and `:3103`, which wrote `[gp()]`; roxygen resolved the unqualified name
against installed brms, because brms is in core's `Suggests`. frmtmb
implements `gp()` itself, as a formula special: parsed at
`R/parse.R:960`, registered at `R/covstruct.R:1587`. It is not exported
(no `export(gp)` in `NAMESPACE`) and has no Rd topic (no
`\alias{gp}` anywhere in `man/`).

**Fix: qualified by removing the link**, not by giving `gp()` a topic.
The reason is the package's own convention. After this change, all 14
roxygen mentions of `gp()` in core, spread over `R/bf.R`,
`R/conditional-effects.R`, `R/confint.R`, `R/predict.R`, `R/priors.R`
and `R/sandwich.R`, write it as plain `` `gp()` `` code with no link,
which is what 12 of them already did. `gp()`
is documented in prose in `vignettes/frmtmb.Rmd`, `inputs.Rmd`,
`compatibility.Rmd` and `brms-migration.Rmd`. The two linked mentions
were the outliers. `ps()` has a topic because `ps()` is an exported
function a nonlinear body calls; `gp()` is a formula special that is
never called directly, so a topic for it is a documentation feature and
not debt, and it is recorded as such in `dev/test-backlog.md`.

`R/predict.R:3085` and `:3103` now read `` `gp()` ``.

### Cross-package link audit, all seven packages

Two measurements, both over every `man/*.Rd` in
`frmtmb`, `frmtmb.eam`, `frmtmb.sample`, `frmtmb.latent`, `frmtmb.ode`,
`frmtmb.spline` and `frmtmb.learn`.

**Undeclared package anchors: zero.** Every `\link[pkg:...]` anchor in
all seven packages names a package in that package's own
`Depends`/`Imports`/`Suggests`, or a standard package. Anchors seen:
core uses TMB, RTMB, stats, splines, mgcv, Matrix, graphics (all
Imports or standard), lme4, DHARMa, brms (all Suggests); eam uses
frmtmb; sample uses frmtmb, tmbstan, mgcv; latent uses frmtmb; ode uses
frmtmb and deSolve; spline uses frmtmb; learn uses frmtmb and stats.
So this class produces no `R CMD check` note anywhere, and `brms:gp`
did not either. The defect was semantic, not a check failure.

**Package-anchored links whose visible text does not name that package**
(the shape that lets a link masquerade as local): 39 total, of which
exactly 2 are defects.

| package | count | verdict |
|---|---|---|
| frmtmb | 2 | **the `brms:gp` pair, fixed** |
| frmtmb.sample | 29 | all on `reexports.Rd:42`, roxygen's own generated form for `@importFrom`-re-exported objects, which are literally frmtmb's functions. Not a defect. |
| frmtmb.latent | 8 | `frmtmb:frm_allfit`, `set_prior`, `mixture_probs`, `mvbf`, `mixture_mvn`, in `hmm.Rd` and `lca.Rd`. frmtmb is a hard `Depends`, so those names are on the user's search path and the link lands on the function the reader will actually call. Not a defect. |
| frmtmb.eam, .ode, .spline, .learn | 0 | |

A third shape, `\link[=topic]` with no local alias of that name, is 69
instances across core (50), sample (10), ode (5), latent (2) and learn
(2). Every one resolves to a `stats` or `base` generic that the package
declares and implements a method for: `simulate`, `predict`, `fitted`,
`residuals`, `confint`, `vcov`, `AIC`, `BIC`, `anova`, `sigma`,
`influence`, `cooks.distance`, `poly`, `message`, `as.numeric`,
`set.seed`, `gaussian`, `poisson`, `inverse.gaussian`, `summary`,
`logLik`. Landing on the generic's page is the right destination for a
reader, so this is not the same defect and nothing was changed.

## 3. frmtmb.sample: sampling.Rmd shows code and no output

**Real, deliberate, and for a reason nobody had written down where a
reader could see it.**

All three articles in the package are unevaluated, not just this one:
`posterior-diagnostics.Rmd` and `brms-posterior.Rmd` set
`eval = FALSE` globally in their `include = FALSE` setup chunk, and
`sampling.Rmd` set it on all nine of its visible chunks. That
consistency is what settles deliberate against accidental. The two
siblings even record a reason, but in a code comment inside a hidden
chunk, which no reader of the built page ever sees.

One of those recorded reasons is wrong. `posterior-diagnostics.Rmd`
said the chunks are "minutes of compute". Measured on this machine
against the model those pages actually use (80 rows, `y ~ x + (1 | g)`,
gaussian):

| step | elapsed |
|---|---|
| `library(frmtmb)` + `library(frmtmb.sample)` | 0.77 s |
| `frm()`, the ML fit | 3.66 s |
| `get_prior()`, either route | 0.03-0.05 s |
| `frm_sample(fit, chains = 4)`, default iter | 14.6 s |
| `frm_sample(fit, chains = 2, iter = 500)` | 2.3 s |

Evaluating the whole of `sampling.Rmd` would cost roughly 90 seconds of
sampling plus the post-processing calls. That is not a reason to
withhold output. The real reasons are that tmbstan and rstan are
`Suggests`, so a build machine need not have them, and that the page's
own last section documents a tmbstan built against StanHeaders 2.39
which samples a standard normal in silence. Baking numbers from a
possibly-wrong Stan into the site would be worse than baking none. The
package's examples already follow that policy: every Stan-touching
example in `R/sample.R`, `R/loo.R` and `R/methods-draws.R` is wrapped in
`if (requireNamespace("tmbstan", quietly = TRUE) && ...)`.

**Done.**

- `vignettes/sampling.Rmd`, after the opening paragraphs: a note saying
  the chunks are shown and not run, that speed is not the reason, and
  what the reason is.
- The opening chunk was split. Library loading and the data frame now
  evaluate; `frm()` and `frm_sample()` stay `eval = FALSE`.
- The two `get_prior()` calls under **The default priors** now
  evaluate. That is the cheapest demonstration per second in the file:
  it needs no sampler (the route is answered by the hook registered in
  `R/zzz.R:20`, and the tmbstan gate is at `R/sample.R:1251`, inside
  `frm_sample()` only), it costs about 0.8 s in total including the
  library loads, and it prints the two prior tables the surrounding
  prose is entirely about.
- `vignettes/posterior-diagnostics.Rmd`: the hidden comment's cost
  claim corrected to the measured 14.6 s, and a reader-visible sentence
  added saying the page is unevaluated on purpose.
- `vignettes/brms-posterior.Rmd`: the same reader-visible sentence.

## 4. local_mocked_bindings() without .package

**Real, still reproduces, and the recorded diagnosis was slightly
wrong.** The note guessed a pkgload upgrade might have fixed it.
pkgload 1.5.3 and testthat 3.3.2 are installed and it fails anyway,
because the missing piece was never pkgload itself. `local_mocked_bindings()`
resolves `.package %||% dev_package()`, and `dev_package()` aborts with
"No packages loaded with pkgload" whenever no package is loaded THROUGH
pkgload, which is exactly the case when a file is run against an
installed package. `test_check()` hides it.

Measured before, one process per file, bare `test_file()` against the
lo-lib installs:

- `tests/testthat/test-influence-plot.R`: 13 passed, **2 errors**, both
  "No packages loaded with pkgload", at lines 68 and 130.
- `extensions/frmtmb.sample/tests/testthat/test-parallel-chains.R`:
  0 passed, **1 error** - but at `dev_namespaces_of` on the line BEFORE
  the mock, an unqualified internal, so the file never reached the
  recorded defect. Fixing only the mock would have left the file broken.

Measured after:

- `test-influence-plot.R`: **21 passed, 0 errors**.
- `test-parallel-chains.R`: **4 passed, 0 errors**.

**Done.**

- `tests/testthat/test-influence-plot.R:78` and `:139`: `.package =
  "frmtmb"` added, with one comment saying why it is named rather than
  inferred. The two calls in the same file that already had
  `.package = "graphics"` were correct and untouched.
- `extensions/frmtmb.sample/tests/testthat/test-parallel-chains.R:53`
  and `:58`: the bare `dev_namespaces_of()` qualified to
  `frmtmb.sample:::dev_namespaces_of()`, and `.package =
  "frmtmb.sample"` added to the mock. Both were needed; the first test
  in that file skips under `load_all()`, so the file is designed to run
  against the installed package and could not.
- `extensions/frmtmb.sample/tests/testthat/test-loo.R:496` and `:507`:
  the same defect, found by the audit rather than recorded. Fixed.
  The file's bare-run error count stays at 7, but the mock error is
  gone and the remaining seven are a different class (the file calls
  `mixture`, `mvbf`, `cumulative`, `rescor_matrix`, `expose_functions`
  and `loo_matrix` without attaching frmtmb).

Audited every `local_mocked_bindings()` call in all seven packages: 14
sites. The partition, re-measured site by site after a review caught
the first count wrong (it said 3/7/2/2 and named three line numbers
under a "2"):

| group | count | sites |
|---|---|---|
| already correct | 4 | test-ce-bands.R:268, test-ce-facets.R:150, test-influence-plot.R:24 and :122 |
| fixed here | 5 | test-influence-plot.R:71 and :134, test-parallel-chains.R:59, test-loo.R:494 and :505 |
| left in frmtmb.ode | 2 | test-ode.R:128 and :132 |
| left in core | 3 | test-compat-register.R:347, test-structure.R:304 and :326 |

The ode pair is recorded in the backlog: verifying a change there means
a full check of a package this lane does not otherwise touch. The three
core sites never reach their mock at all, because those files die
earlier on unqualified internals, so fixing the mock alone would not
make them run; that shape is recorded in the backlog too and is not
this item.

## 5. dev/test-backlog.md triage

The file held 17 entries under an `Open` heading and about 58 more under
`Done`, `Fixed`, `Already worked`, `Mitigated` and `Verified immune`. I
re-measured the 17 by running them, not by reading for them; the archive
sections were left alone, because that is what they are for.

**Result: 6 done, 1 moot, 1 halved, 9 remain, 2 of the 9 restated.**

Closed as done:

| item | measurement |
|---|---|
| Interval censoring (cens code 2 + y2) | `cens_code_map` has `interval = 2` at R/frame.R:531, NA y2 allowed off interval rows (:1086) and refused on them (:1383); a gaussian fit over all four codes converges |
| Discrete truncation needs F(lb - 1) | poisson `trunc(lb = 2)` logLik matches a hand-rolled `ppois(lb - 1)` reference to 1.1e-13. Was already recorded under "Verified immune", so it was a duplicate |
| `\|ID\|` new-level prediction from the joint block | new level under `(1 \| ID \| g)` across mu and sigma returns fit 1.5254, se 0.4976 |
| dpar-formula reordering must not change logLik | `student()` with `sigma` and `nu` swapped: difference 0.0e+00 |
| `(1\|a:b)` level ordering | tests/testthat/test-aliased-grouping.R:144 pins it against `levels(droplevels(a:b))`; `(1\|a*b)` is under "Verified immune" |
| dpar/nlpar names with dots or underscores | refused at R/bf.R:157 with the message the entry asks for. Already listed under "Addressed in v0.6", so a duplicate |
| zi predict grid | `type = "response"` equals `(1 - zi) * conditional` to 0.0e+00 |

Closed as moot: "zprob on a non-zi model returns 0 not garbage". It now
refuses by name, `Unknown dpar: 'zi' for response 'y'. Available: mu`,
which is a better answer than the one the entry asked for.

Halved: the offset cluster. `y ~ 0 + offset(o)` fits, and offset inside
a dpar formula was already done; only the `frm(..., offset = )` argument
form remains, and it is still an "unused argument" error.

Restated, because what they predict is not what happens:

- Profile CIs on boundary parameters "must warn rather than hang". They
  do not hang, and they fail two ways, not one. Over 17 no-signal
  gaussian fits, `confint(parm = "theta_1", method = "profile")` either
  dies inside `approx()` with "need at least two non-NA values to
  interpolate" (10 of 17), or returns `lwr = NA` with a finite `upr`
  and warns only "NA/NaN function evaluation", naming no component
  (7 of 17). A review's sweep split the same two modes 6 and 11, so the
  ratio is seed-dependent and both modes are always present. The
  half-answer is the worse of the two and my first pass recorded only
  the error; the backlog entry now carries both. One fix covers both:
  catch it and return NA bounds with a named warning, as the Wald path
  already does.
- Constant-weight invariance for gaussian was filed as high priority.
  The behavior is already correct: a constant `weights(w)` of 2 leaves
  the coefficients invariant to 5.9e-06 and exactly doubles the logLik
  (-347.5174 against -173.7587). Only the test and the doc line are
  missing, so it stays open but is no longer a correctness risk.

Confirmed still real, unchanged: `gam`-style `exclude=`
(`predict(exclude = "s(z)")` warns "ignoring unknown arguments"); no
`link_zi`/`link_hu` (`zero_inflated_poisson(link_zi = "probit")` is an
unused-argument error); `cens()` refused for discrete families (the
message says so); `dharma_residuals()` on a censored fit returns
silently, neither warning nor refusing;
`conditional_effects(method = "predict")` drops an unevaluable
`vint()`/`vreal()` payload, because `ce_aterms()`'s `strict` set at
R/conditional-effects.R:14 is trials, se, trunc_lb and trunc_ub and
neither payload is in it. Two entries are standing design limits rather
than work: OSA under `cens()` covering the uncensored rows only, and the
`quadrature` breakdown on hard likelihoods, whose fix is an integrator
to build.

The file gained a dated triage header stating those counts, two
"Closed at the 2026-09-07 triage" subsections holding the closures with
their reasons, and a "Recorded by the loose-ends lane" section for three
things this lane found and did not do: a topic for `gp()`, the two
remaining `local_mocked_bindings()` sites in frmtmb.ode, and the test
files that only run under `pkgload::load_all()`.

## Verification

Private library only, `scratchpad/lo-lib`, holding the worktree builds
of frmtmb, frmtmb.sample and frmtmb.spline. Nothing was installed into
the shared user library, which was read last through `R_LIBS`.

### Suites, one process each

| package | files | passed | failed | errors | warnings | skipped |
|---|---|---|---|---|---|---|
| frmtmb | 131 | 7440 | 0 | 0 | 1 | 133 |
| frmtmb.sample | 12 | 982 | 0 | 0 | 0 | 2 |
| frmtmb.spline | 8 | 326 | 0 | 0 | 0 | 0 |

The one core warning is `test-prior-compat.R:505`, "Optimizer did not
report convergence: singular convergence (7)". It is in a file this
lane did not touch, on the priors and compatibility path this lane did
not touch.

frmtmb.spline was re-run after the item-1 revert and the new test:
**326 passed**, up from 322, the four added by `test-span.R`'s
two-`ps()` case. `test-span.R` alone goes 48 to 52 assertions. The two
item-4 files were re-run under a bare `test_file()` against the
installed packages and still pass clean: `test-influence-plot.R` 21
passed / 0 errors, `test-parallel-chains.R` 4 passed / 0 errors.

Every skip is an opt-in tier, not a silent hole: 91 core and 2 sample
on `FRMTMB_BRMS_FIT_TESTS=true` for brms fit tests, 41 core on the same
variable for the Stan tier, and 1 core on `FRMTMB_FUZZ=true` for the
grammar fuzzer.

### Roxygen

Idempotent on all three, each checked by `git status` being
byte-identical before and after a run. On frmtmb the first run rewrote
exactly one file, `man/frm_lp_basis.Rd`, which is the `gp()` change and
nothing else, and a second run wrote nothing. On frmtmb.spline and
frmtmb.sample a run after the edits wrote nothing at all.

Two roxygen link warnings were raised and left alone, because both live
in files another lane owns: `families.R:1563` cannot resolve
`fam_binomial`, and `links.R:243` cannot resolve `log_inv_logit`.
Neither is a `R CMD check` finding.

### R CMD check --as-cran, with the manual

Run on a tarball from `R CMD build`, with pandoc 3.8.3 and TinyTeX on
PATH. A first attempt raised a `checking top-level files` NOTE,
"Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc'
being installed". That was the harness, not the package: `R CMD check`
probes `Sys.which("pandoc")` and `RSTUDIO_PANDOC` alone reaches only
rmarkdown, which is why the vignettes built while the check still
called pandoc missing. Putting the same directory on PATH cleared it,
and both packages were re-checked from scratch.

**frmtmb.spline: 1 WARNING, 1 NOTE.** `checking PDF version of manual
... OK`, and `frmtmb.spline-manual.pdf` was produced.

- NOTE, `checking HTML version of manual`: "Skipping checking math
  rendering: package 'V8' unavailable". The one expected NOTE.
- WARNING, `checking CRAN incoming feasibility` [196s], three
  sub-items, all pre-existing: "New submission"; "Strong dependencies
  not in the CRAN or BioC software repositories: frmtmb", because core
  is not published; and a 301 on the URL
  `https://aforren1.github.io/frmtmb/frmtmb.spline` named by
  DESCRIPTION and `man/frmtmb.spline-package.Rd`, which wants a
  trailing slash.

**frmtmb.sample: 1 WARNING, no NOTE.** `checking PDF version of manual
... [21s] OK`, `checking HTML version of manual ... OK` (no V8 NOTE
here, because this manual has no math for V8 to render), and
`checking re-building of vignette outputs` passed, which is the stage
that exercises the newly evaluated chunks.

- WARNING, `checking CRAN incoming feasibility` [166s], all
  pre-existing: "New submission"; "Strong dependencies not in the CRAN
  or BioC software repositories: frmtmb"; and "Suggests or Enhances not
  in mainstream repositories: frmtmb.latent". Both are unpublished
  packages of this monorepo.

The evaluated chunks were checked in the built artifact, not just
assumed: the tarball's `inst/doc/sampling.html` carries 20 output
lines, among them both `get_prior()` tables, each naming its route on
its first line, which is the thing the surrounding prose describes.

**frmtmb: 3 NOTEs, no WARNING, no ERROR.** `checking PDF version of
manual ... [15s] OK`, and `frmtmb-manual.pdf` was produced. The
in-check suite ran `Running 'testthat.R' [30m]` and reported OK.

- NOTE, `checking CRAN incoming feasibility` [327s]: "New submission",
  and nothing else. Core is a WARNING lighter than the extensions here,
  because it depends only on CRAN packages and has no redirecting URL.
- NOTE, `checking examples` [107s]: four examples over the 5s advisory,
  `profile.frmtmb_fit` (19.87s elapsed), `residuals.frmtmb_fit` (9.30s),
  `pp_check` (9.63s) and `dharma_residuals` (6.64s). All pre-existing,
  in topics this lane did not touch, and the elapsed figures are
  inflated by three R processes sharing this machine.
- NOTE, `checking HTML version of manual` [49s]: "Skipping checking
  math rendering: package 'V8' unavailable". The expected NOTE.

### Summary

| package | ERRORs | WARNINGs | NOTEs |
|---|---|---|---|
| frmtmb | 0 | 0 | 3 (CRAN new submission; slow examples; V8) |
| frmtmb.sample | 0 | 1 (CRAN new submission, unpublished deps) | 0 |
| frmtmb.spline | 0 | 1 (CRAN new submission, unpublished dep, URL 301) | 1 (V8) |

No ERROR anywhere. Every WARNING and NOTE is either the expected V8 one
or a pre-existing consequence of an unpublished monorepo, and none was
introduced by this lane.

### Re-verification after the review round

The item-1 revert changed `frmtmb.spline`'s R code and tests, and the
item-3 correction changed two `frmtmb.sample` vignettes, so both
packages were checked again from scratch rather than left on a stale
result. Both statuses are **identical to the pre-revert run**:

| package | before | after the round |
|---|---|---|
| frmtmb.spline | 1 WARNING, 1 NOTE | 1 WARNING, 1 NOTE |
| frmtmb.sample | 1 WARNING, 0 NOTE | 1 WARNING, 0 NOTE |

Same causes as before: the CRAN-incoming WARNING and, on spline, the
V8 NOTE. Manuals built on both (`checking PDF version of manual ...
OK`). The in-check suites ran `[92s]` and `[299s]` clean. The rebuilt
`inst/doc/sampling.html` still carries its 20 output lines and both
`get_prior()` route tables, so removing the second-count from the prose
did not disturb the evaluated chunks. `frmtmb` core was not re-checked:
nothing under it changed in this round.

Roxygen re-confirmed idempotent on frmtmb.spline after the revert and
the new test: a run wrote nothing.

### Worktree state

Nothing committed. Eleven files modified and one added:

- R/predict.R, man/frm_lp_basis.Rd, tests/testthat/test-influence-plot.R
- extensions/frmtmb.spline/R/curve-feature.R
- extensions/frmtmb.spline/tests/testthat/test-span.R
- extensions/frmtmb.sample/tests/testthat/test-loo.R
- extensions/frmtmb.sample/tests/testthat/test-parallel-chains.R
- extensions/frmtmb.sample/vignettes/sampling.Rmd
- extensions/frmtmb.sample/vignettes/posterior-diagnostics.Rmd
- extensions/frmtmb.sample/vignettes/brms-posterior.Rmd
- dev/test-backlog.md
- dev/loose-findings.md (new)

`dev/reviews/2026-09-08-loose.md` is also untracked but is the
reviewer's, not this lane's.

No dev/*.log and no scratch script was left in the worktree; every
build, check and probe artifact lives under the scratchpad. Core's
DESCRIPTION was not touched and still lists no extension in Imports or
Suggests. None of the five reserved areas was edited: R/links.R,
R/families.R's ordinal link switch, R/frame.R's se() and cens() gates,
the addition-term allow-list, R/predict.R's band code (only two roxygen
comment lines in the frm_lp_basis block changed) and
extensions/frmtmb.learn.
