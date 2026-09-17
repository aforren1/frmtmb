# Classed conditions, item 2.6e (lane wt-conditions)

Base: e031e8c, frmtmb 0.59.0 and frmtmb.sample 0.7.0. No Version field
changed.

## The contract

brms 2.23.0 raises its refusals through `stop2()`, which calls
`rlang::abort(.subclass = c(.subclass, "brms_error"))`. The class of such
an error is `c("brms_error", "rlang_error", "error", "condition")`, so
`tryCatch(brms_error = )` catches any refusal. Before this lane, no
condition that frmtmb raised had a class of its own, except
`frmtmb_fit_error` and two warning classes.

## Design decisions

Evidence for decisions 1 to 3 is in `dev/conditions-log/design.txt`,
made by `dev/conditions-design-evidence.R` from the brms 2.23.0 CRAN
tarball (`dev/brms-suite/brms`) and the installed libraries.

### 1. Base R conditions, not rlang

The class vectors are:

| kind | class vector |
|---|---|
| error | `c("frmtmb_error", "error", "condition")` |
| warning | `c("frmtmb_warning", "warning", "condition")` |
| message | `c("frmtmb_message", "message", "condition")` |

Reasons, from the evidence file:

- brms's own R sources name `brms_error` on 2 lines, both in the
  definition and documentation of `stop2()`. They name `rlang_error`,
  `catch_cnd`, `cnd_` and `trace_back` on 0 lines. brms never catches its
  own class.
- brms's test suite names `brms_error` on 5 lines, all in
  `tests.stop2.R`, which tests `stop2()` itself. It names `rlang_error`
  on 0 lines. The only thing brms-shaped caller code relies on is the
  package class, which base R conditions give.
- rlang is not in frmtmb's dependency closure. The closure of the 13
  Imports has 39 packages and rlang is not one of them. Adding rlang
  would add a 1.9 MB dependency (rlang 1.3.0 as installed) for the
  second class name and for rlang's backtrace printing, which nothing
  measured uses.
- `rlang_error` changes how an uncaught error prints (`Error in f():`
  becomes rlang's bullet format), which would change what a user sees.
  Base conditions print exactly as `stop()` did.

The `simple*` classes are dropped, which matches brms: its class vector
has no `simpleError` either. One test asserted `class = "simpleError"`
(`tests/testthat/test-compat-register.R`); it now asserts
`frmtmb_error`.

`call.` is kept per site. The helpers take `call.` with the default of
the base functions, and the rewrite changed only the function name, so
each site keeps what it had. On the base, `dev/conditions-scan2.R` found
`call.` passed at 1,209 `stop()` sites and 44 `warning()` sites.
brms's `stop2()` always records
the calling function; frmtmb does not copy that, because it would change
what every refusal prints.

### 2. One helper family in core, with a subclass per extension

`frm_stop()`, `frm_warning()` and `frm_message()` are in `R/conditions.R`,
exported and documented on `?frmtmb-conditions`, which
`?frmtmb-extension-api` links. They take the arguments of the base
functions, so a call site is a rename. They build the message with
`.makeMessage()`, as the base functions do, and record `sys.call(-1)`
(`frm_stop`, `frm_warning`) or `sys.call()` (`frm_message`), which is
the call each base function records. `tests/testthat/test-conditions.R`
asserts the message, the call and the text `try()` prints against the
base function, including a call forced inside a promise.

An extension raises its own subclass as well:
`c("frmtmb_eam_error", "frmtmb_error", "error", "condition")`. The
subclass is found from `topenv()` of the calling frame, so no call site
passes it and the rewrite stays a rename. The name is the package name
with dots made underscores, which is the naming the two existing
extension classes already used (`frmtmb_eam_units_warning`). The
subclass names the package whose FUNCTION raised the condition: a core
validator called by an extension raises a plain `frmtmb_error`. A
function defined outside a namespace gets frmtmb's classes only.

`class =` puts more classes in front, as brms's `.subclass` does. It
carries the three classes that already existed:
`frmtmb_fit_error`, `frmtmb_ps_span_warning` and
`frmtmb_eam_units_warning`.

### 3. Warnings and messages are classed too, which brms does not do

brms has `warning2()`, which is `warning(..., call. = FALSE)` and adds no
class, and no `message2()`; its 42 `message()` calls are unclassed. So
matching brms literally would leave frmtmb's warnings and messages
unclassed. They are classed anyway, as the task specifies, because the
class is additive: every handler for `warning` or `message`, including
brms-shaped caller code, still matches, and
`suppressWarnings(classes = "frmtmb_warning")` becomes possible. This is
the one place the lane goes beyond brms rather than matching it, and it
is listed for the user below.

### 4. Message text

No message text changed. `dev/conditions-rewrite.R verify` proves it
structurally: for every rewritten file, the parse tree of the base file
put through the rename deparses identically to the parse tree of the new
file, so no string literal and no argument changed. The lines that the
longer name pushed past 80 columns (listed as `OVER80` in
`dev/conditions-log/rewrite.txt`) were reflowed by hand between
arguments, and the proof was run again after the reflow and after the
second pass; its result is in the rewrite block below. The suites that
assert message text, and every `test-message-uniqueness.R`, pass
(below).

### 5. Conditions frmtmb does not raise

Errors from RTMB, TMB, Matrix, base R and the optimizers keep their
class. The one place frmtmb catches and rethrows them,
`fit_error_context()`, already built its own condition; it now calls
`frm_stop(msg, call. = FALSE, class = "frmtmb_fit_error")`, so an
optimizer failure is a `frmtmb_error`.

## The rewrite

`dev/conditions-rewrite.R rewrite <snapshot>` copies each R/ tree and
rewrites it from `getParseData()`: it renames the head token of each
`stop`, `warning` and `message` call and moves the continuation lines of
a multi-line call right by the same four columns. Its counts
(`dev/conditions-log/rewrite.txt`):

<!-- BEGIN GENERATED: rewrite counts -->
```
frmtmb exempt stop(<condition>) rethrow                          3
frmtmb message                                                  15
frmtmb stop                                                    770
frmtmb stop(structure(frmtmb_fit_error)) special                 1
frmtmb warning                                                  27
frmtmb warning(warningCondition()) special                       1
frmtmb.coupling stop                                            26
frmtmb.coupling warning                                          1
frmtmb.eam stop                                                111
frmtmb.eam warning(warningCondition()) special                   1
frmtmb.latent stop                                              46
frmtmb.latent warning                                            1
frmtmb.learn stop                                               22
frmtmb.ode stop                                                 92
frmtmb.ode warning                                               3
frmtmb.sample message                                            9
frmtmb.sample stop                                             105
frmtmb.sample warning                                           11
frmtmb.spline stop                                              38
frmtmb.spline warning                                            1
frmtmb.spline warning(warningCondition()) special                2
frmtmb-inst stop                                                42

converted per package (every non-exempt row above summed):
frmtmb             814
frmtmb-inst         42
frmtmb.coupling     27
frmtmb.eam         112
frmtmb.latent       47
frmtmb.learn        22
frmtmb.ode          95
frmtmb.sample      125
frmtmb.spline       41
ALL               1325
exempt               3

proof over every file, after the hand reflow and the second
pass (rewrite-verify.txt):
files compared 113, mismatched 0, lines over 80 columns 0
proof after the second pass over inst/ (rewrite-inst.txt):
files compared 113, mismatched 0, lines over 80 columns 0
```
<!-- END GENERATED: rewrite counts -->

### Sites decided one by one

Exempt, and named in the census:

- `R/fit.R`, `fit_assembled()`: `stop(e)`. It rethrows the optimizer's
  own error after every restart has failed; `fit_error_context()` then
  wraps it as a `frmtmb_fit_error`.
- `R/fit.R`, `optimizer_from_best()`: `stop(e)`. The same error when
  there is no better point to restart from, wrapped the same way.
- `R/fit.R`, `fit_error_context()`: `stop(e)`. It rethrows a
  `frmtmb_fit_error` from the inner autoscale fit, which is already
  classed and must not be wrapped twice.
- `R/conditions.R`: `stop(cnd)`, `warning(cnd)` and `message(cnd)`, the
  helpers themselves.

Converted by the script, with a shape of their own:

- `R/fit.R`, `stop(structure(class = c("frmtmb_fit_error",
  "simpleError", ...), list(message = msg, call = NULL)))` is now
  `frm_stop(msg, call. = FALSE, class = "frmtmb_fit_error")`.
- `warning(warningCondition(..., class = K))` in `R/ps.R`, frmtmb.eam
  `R/ddm-shared.R`, and frmtmb.spline `R/curve.R` and
  `R/curve-deriv.R` is now `frm_warning(..., class = K, call. = FALSE)`,
  because `warningCondition()` records no call.
- frmtmb.ode `ode_warn_once()`: `warning(..., call. = FALSE)`, a
  wrapper whose warnings are frmtmb.ode's.
- `R/zzz.R`, `notify_once()`: `message(...)`, the once-per-session
  notices.
- `tryCatch()` handlers that call `stop(conditionMessage(e), ...)`:
  they build a new message, so the condition is frmtmb's.
- 42 `stop()` calls in `inst/bcm/*.R` and `inst/rl/rw-delta.R`, in a
  second pass. These are the worked-example families that the vignettes
  and tests load with `source()`; the first sweep found three of their
  refusals unclassed. They run from the global environment, so they get
  `frmtmb_error` with no subclass, and they need frmtmb attached, which
  they already did for `frmtmb_family()`.

Not converted:

- 22 `stopifnot()` calls in core and 1 in frmtmb.sample (below).
- `simpleError()`, `simpleWarning()`, `signalCondition()` and
  `errorCondition()`: `dev/conditions-scan.R` found no call of any of
  them.

`stopifnot()` is not converted. brms has 318 `stopifnot()` calls and
leaves them unclassed, and converting them would change their text
(`length(pars) == 2 is not TRUE`), which this lane must not do. Most are
internal assertions, but some are reachable from user input. The ones in
`R/priors.R` that check a prior string are, measured by
`dev/conditions-prior-probe.R` (`dev/conditions-log/prior-probe.txt`):
`set_prior("normal(0)")` raises a `simpleError` with the text
`length(pars) == 2 is not TRUE`, and so do `student_t(3, 0)`,
`exponential(-1)`, `gamma(1)` and `lkj()` with their own conditions.
That is a refusal a user meets, with a message that does not say what
is wrong, and it is not classed. It is a follow-up for the user (below),
not a silent gap.

A defect found on the way and not fixed: `frmtmb.latent::hmm_starts(1)`
fails with base R's `$ operator is invalid for atomic vectors` instead of
a refusal, because `hmm_rspec()` reads `fit$spec` before it checks that
`fit` is a fit.

## Guards

### (a) The census

`dev/conditions-census-gen.R` writes
`tests/testthat/test-conditions-census.R` into frmtmb and each of the
seven extensions from one template. It parses every file in `R/` and
fails on a bare `stop()`, `warning()` or `message()`, on one called as
`base::stop()`, and on one passed as a value (`error = stop`,
`do.call("message", ...)`, a formal default). An exemption names a site
as `file | enclosing top-level name | call`, and an exemption that
matches no site fails as stale. The test also parses a planted snippet
and asserts it finds exactly the five planted shapes and not
`e$message`, `list(message = )` or `frm_stop()`.

Seen failing, with `dev/conditions-census-run.R` on scratch copies:

<!-- BEGIN GENERATED: census -->
```
lane tree (census-lane.txt):
RESULT frmtmb           pass 4 fail 0 error 0 skip 0
RESULT frmtmb.coupling  pass 4 fail 0 error 0 skip 0
RESULT frmtmb.eam       pass 4 fail 0 error 0 skip 0
RESULT frmtmb.latent    pass 4 fail 0 error 0 skip 0
RESULT frmtmb.learn     pass 4 fail 0 error 0 skip 0
RESULT frmtmb.ode       pass 4 fail 0 error 0 skip 0
RESULT frmtmb.sample    pass 4 fail 0 error 0 skip 0
RESULT frmtmb.spline    pass 4 fail 0 error 0 skip 0

base sources with the census added (census-absent-base.txt):
RESULT frmtmb           pass 2 fail 2 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / actu
    Expected `stale` to be identical to `character(0)`. / Differences: /    
RESULT frmtmb.coupling  pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / actu
RESULT frmtmb.eam       pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / actu
RESULT frmtmb.latent    pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / actu
RESULT frmtmb.learn     pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / actu
RESULT frmtmb.ode       pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / actu
RESULT frmtmb.sample    pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / actu
RESULT frmtmb.spline    pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / actu

lane plus a bare stop() in core R/utils.R and `error = stop`
in frmtmb.eam R/ddm-shared.R (census-absent-plant.txt):
RESULT frmtmb           pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: /     
RESULT frmtmb.coupling  pass 4 fail 0 error 0 skip 0
RESULT frmtmb.eam       pass 3 fail 1 error 0 skip 0
    Expected `bare` to be identical to `character(0)`. / Differences: / `act
RESULT frmtmb.latent    pass 4 fail 0 error 0 skip 0
RESULT frmtmb.learn     pass 4 fail 0 error 0 skip 0
RESULT frmtmb.ode       pass 4 fail 0 error 0 skip 0
RESULT frmtmb.sample    pass 4 fail 0 error 0 skip 0
RESULT frmtmb.spline    pass 4 fail 0 error 0 skip 0

lane with fit_error_context()'s stop(e) renamed, so its
exemption is stale (census-absent-stale.txt):
RESULT frmtmb           pass 3 fail 1 error 0 skip 0
    Expected `stale` to be identical to `character(0)`. / Differences: / `ac
RESULT frmtmb.coupling  pass 4 fail 0 error 0 skip 0
RESULT frmtmb.eam       pass 4 fail 0 error 0 skip 0
RESULT frmtmb.latent    pass 4 fail 0 error 0 skip 0
RESULT frmtmb.learn     pass 4 fail 0 error 0 skip 0
RESULT frmtmb.ode       pass 4 fail 0 error 0 skip 0
RESULT frmtmb.sample    pass 4 fail 0 error 0 skip 0
RESULT frmtmb.spline    pass 4 fail 0 error 0 skip 0
```
<!-- END GENERATED: census -->

### (b) Runtime catch, and the sweep

`tests/testthat/test-conditions.R` in core, and a generated
`test-conditions.R` in each extension (`dev/conditions-ext-tests.R`),
catch one real refusal per package with `tryCatch(frmtmb_error = )` and
assert the full class vector, and in each extension catch it again by
its subclass alone. On the base build every one of the eight fails:
the refusal is not caught and the test errors with the refusal's own
message (`dev/conditions-log/ext-conditions-tests.txt`; for core,
`dev/conditions-log/core-conditions-base.txt`, where the frm() block
fails behaviorally and the helper blocks fail because the helpers do
not exist there).

The sweep runs every test file of every suite, one per process, with a
`trace()` on testthat's `expect_error()` (both the namespace and the
attached binding, which `dev/conditions-traceprobe.R` showed are both
needed) that records the class of every condition an `expect_error()`
caught. `dev/conditions-sweep-summary.R` counts them and matches each
non-`frmtmb_error` message against the string literals of every
`stop()`, `frm_stop()` and `stopifnot()` call in frmtmb's sources, to
separate an unclassed error frmtmb raised from one it did not.

<!-- BEGIN GENERATED: sweep -->
```
lane build (sweep-lane.txt):
== sweep: conditions the suites' own expect_error() caught ==
package          caught  frmtmb_error  other
frmtmb             1419          1397     22
frmtmb.coupling      33            32      1
frmtmb.eam          193           192      1
frmtmb.latent        55            55      0
frmtmb.learn         56            56      0
frmtmb.ode           78            78      0
frmtmb.sample       142           138      4
frmtmb.spline        56            55      1
ALL                2032          2003     29

frmtmb_error class vectors seen:

frmtmb_coupling_error/frmtmb_error/error/condition 
                                                28 
     frmtmb_eam_error/frmtmb_error/error/condition 
                                               170 
                      frmtmb_error/error/condition 
                                              1488 
     frmtmb_fit_error/frmtmb_error/error/condition 
                                                 3 
  frmtmb_latent_error/frmtmb_error/error/condition 
                                                42 
   frmtmb_learn_error/frmtmb_error/error/condition 
                                                34 
     frmtmb_ode_error/frmtmb_error/error/condition 
                                                76 
  frmtmb_sample_error/frmtmb_error/error/condition 
                                               108 
  frmtmb_spline_error/frmtmb_error/error/condition 
                                                54 

not frmtmb_error, by class and whether the message carries a literal
of a stop()/frm_stop()/stopifnot() call in frmtmb's sources:
                                        frmtmb_literal
class                                    FALSE TRUE
  brms_error/rlang_error/error/condition     3    1
  simpleError/error/condition               24    1

every distinct non-frmtmb_error message (class | package | text):
[foreign] brms_error/rlang_error/error/condition | frmtmb | Invalid addition
[foreign] brms_error/rlang_error/error/condition | frmtmb | The following va
[foreign] brms_error/rlang_error/error/condition | frmtmb | New factor level
[foreign] simpleError/error/condition | frmtmb | 'arg' should be one of \res
[foreign] simpleError/error/condition | frmtmb | brms parameterizes a bym2 f
[foreign] simpleError/error/condition | frmtmb | 'arg' should be one of \epr
[foreign] simpleError/error/condition | frmtmb | object 'z' not found
[foreign] simpleError/error/condition | frmtmb | factor f has new levels zz
[foreign] simpleError/error/condition | frmtmb | bad x: 3
[foreign] simpleError/error/condition | frmtmb | 'arg' should be one of \fit
[foreign] simpleError/error/condition | frmtmb | Stan parameter simo_sigma_1
[foreign] simpleError/error/condition | frmtmb | invalid type (list) for var
[foreign] simpleError/error/condition | frmtmb | variable lengths differ (fo
[foreign] simpleError/error/condition | frmtmb | object 'no_such_thing' not 
[foreign] simpleError/error/condition | frmtmb | 'softit' link not recognise
[foreign] simpleError/error/condition | frmtmb | unused argument (priors = p
[foreign] simpleError/error/condition | frmtmb | unused argument (lower = c(
[foreign] simpleError/error/condition | frmtmb | unused argument (upper = c(
[foreign] simpleError/error/condition | frmtmb | invalid type (list) for var
[foreign] simpleError/error/condition | frmtmb.coupling | 'arg' should be on
[foreign] simpleError/error/condition | frmtmb.eam | 'arg' should be one of 
[foreign] simpleError/error/condition | frmtmb.sample | The following variab
[foreign] simpleError/error/condition | frmtmb.sample | The following variab
[foreign] simpleError/error/condition | frmtmb.sample | The following variab
[foreign] simpleError/error/condition | frmtmb.sample | object '(Intercept)'
[foreign] simpleError/error/condition | frmtmb.spline | 'arg' should be one 
[LITERAL] brms_error/rlang_error/error/condition | frmtmb | The model does n
[LITERAL] simpleError/error/condition | frmtmb | error in evaluating the arg

not frmtmb_error, by origin (patterns below, first match wins):
    7  R evaluation: an argument or object that does not exist
    6  base match.arg() inside a frmtmb function
    4  brms itself (the agreement tests)
    4  stats::model.frame() and its relatives
    3  test code: a test helper or a deliberate absent case
    3  the posterior package
    1  R's S4 dispatch rewraps a frmtmb_error as simpleError
    1  stats::make.link()
```
<!-- END GENERATED: sweep -->

How to read it. The first sweep, before the second rewrite pass, found
three refusals from `inst/rl/rw-delta.R` and `inst/bcm/marginal.R` that
were unclassed; that is what brought `inst/` into the rewrite and the
census, and the sweep above is the run after it. Of the conditions that
are not a `frmtmb_error`, the only one frmtmb raised is the S4 case.
`getME(fit, "Zt")` on a multivariate fit evaluates
`Matrix::t(mu_Z())` (`R/interop.R` line 421), `mu_Z()` raises the
`frmtmb_error` "Multiple responses have dpar 'mu'", and R's S4 method
dispatch catches it while it selects a method and raises a new
`simpleError` with the text `error in evaluating the argument 'x' in
selecting a method for function 't': ...`. Computing `mu_Z()` before
the call to `t()` would keep the class, but it would also remove that
prefix from the message, so it is left for the user rather than done
in a lane that must not change message text. The rest are raised by
base R, stats, posterior, brms itself or
the test code. The `match.arg()` refusals are frmtmb functions refusing
an argument, raised by base R; brms uses `match.arg()` the same way and
they are left alone under rule 5.

The sweep was not run on the base build. On the base, the count of
`frmtmb_error` is 0 by construction, since no code there raises that
class (`frmtmb_fit_error` did not inherit from it), so the run would
measure an identity. The runtime tests above are the base-build
evidence instead: all eight fail there.

### (c) Message text

The structural proof in decision 4, plus the suites below, which assert
message text by regular expression, plus the eight
`test-message-uniqueness.R` files. Those tests found their templates by
the call name `stop`, `warning` or `message`, so each now pools the base
name with its `frm_` name; without that change they would have found no
templates at all.

Other tests that keyed on the name `stop` were updated the same way:
`tests/testthat/test-arg-refusal.R` (`ar_refusers` and the two
unconditional refusal helpers, whose body is now `frm_stop()`).

## Suites and checks

<!-- BEGIN GENERATED: suites -->
```
lane build (sweep-lane.txt):
== suites, one file per R process ==
runner dev/conditions-sweep-onefile.R, sums failures AND errors
frmtmb           files 142 (no result line 0) blocks 1551 PASS 11733 FAIL 0 
   .start 06:50:45 
   .end 06:59:05 
frmtmb.coupling  files   9 (no result line 0) blocks   65 PASS   438 FAIL 0 
frmtmb.eam       files  26 (no result line 0) blocks  275 PASS  1640 FAIL 0 
frmtmb.latent    files   9 (no result line 0) blocks   69 PASS   339 FAIL 0 
frmtmb.learn     files  14 (no result line 0) blocks   98 PASS   467 FAIL 0 
frmtmb.ode       files  10 (no result line 0) blocks  117 PASS   333 FAIL 0 
frmtmb.sample    files  23 (no result line 0) blocks  231 PASS  1702 FAIL 0 
frmtmb.spline    files  14 (no result line 0) blocks   92 PASS   441 FAIL 0 
   .start 06:39:24 
   .end 06:47:13 

```
<!-- END GENERATED: suites -->

The frmtmb.sample logs record StanHeaders 2.32.10 in all 23 files, so
the sampler tests ran against the pinned library.

`R CMD check --as-cran`, once each, without `--no-manual`
(`dev/conditions-check.ps1`; logs `dev/conditions-log/check-core.txt`
and `check-sample.txt`):

```
core:   * checking HTML version of manual ... [25s] NOTE
core:   Skipping checking math rendering: package 'V8' unavailable
core:   Status: 1 NOTE
sample: Status: OK
```

The one core NOTE is the expected V8 math-rendering note.

## Needs the user

- Warnings and messages are classed although brms classes neither
  (decision 3). Removing it is one change to `frm_condition_class()`.
- Every extension's `frmtmb (>= 0.59.0)` floor must move to the frmtmb
  release that exports `frm_stop()`. No Version field was changed.
- `stopifnot()` sites reachable from user input stay `simpleError`, and
  their text does not say what is wrong: `set_prior("normal(0)")` gives
  `length(pars) == 2 is not TRUE`. Converting them to `frm_stop()`
  changes their text, which this lane was not to do.
- `frmtmb.latent::hmm_starts(1)` fails with a base `$` error instead of
  a refusal. Found, not fixed.
- `match.arg()` refusals inside frmtmb functions stay `simpleError`.
- `getME(mv_fit, "Zt")` loses the class through S4 dispatch (sweep
  section). The fix is one line in `R/interop.R` and drops R's prefix
  from the message.

## Punch round 1

Review: `dev/reviews/20260917-conditions.md`. Scripts are
`dev/conditions-punch-*.R`; logs are in `dev/conditions-punch-log/`.
The build is `conditions-lib`, installed from this tree. The user
decided A (warnings and messages stay classed), B (every user-reachable
`match.arg()` and `stopifnot()` refusal becomes a named `frmtmb_error`)
and C (the `getME()` S4 fix). This section is assembled by
`dev/conditions-punch-findings.R` from the logs named in each block.

### What changed

- **Helpers** (`R/conditions.R`, review MINOR 4). Given one condition
  object, `frm_stop()`, `frm_warning()` and `frm_message()` signal that
  object with the frmtmb classes added in front, as `stop(cond)` does.
  `frm_warning(immediate. = TRUE)` holds `options(warn = 1)` for the
  one signal when `warn` is 0, which is what `immediate.` does in base.
  The subclass lookup skips frames of base R closures, so
  `lapply(x, frm_stop)` gets the subclass of the function that passed
  it, never `base_error`. `package =` overrides the lookup. The default
  `call. = TRUE` stays, and its extra `Calls:` line is documented; the
  one site that used the default (`cr_adjust()`) now passes
  `call. = FALSE`. `noBreaks.` stays ignored and is documented so.
- **`frm_match_arg()`** replaces `match.arg()` at the 34 user sites of
  the review's list: 18 in core and 16 in the extensions. It takes the
  choices from the caller's formals and keeps `NULL`, identical
  choices, partial matching and `several.ok`. It refuses with
  `` `type` must be one of "a", "b", not character "x" `` and records
  the call of its caller. The two internal sites (`autoscale_map()`,
  `ln_recurse()`) keep `match.arg()`.
- **`stopifnot()`** at every user site of the review's list. The prior
  string checks are one arity table, `prior_dist_params`, checked
  before the `switch()`, so all nine densities refuse the same way; the
  exponential rate has its own message. The five internal sites stay.
- **`getME(mv, "Zt")`**: `mu_Z()` is forced before `Matrix::t()`.
- **Other unclassed errors.** `hmm_rspec()` refuses a `fit` that is not
  a list with `spec`. The test is not "is a `frmtmb_fit`", because the
  simulator passes an unclassed shim. `frm_task_simulate()` refuses a
  `family` that is not a list. `check_newdata_frame()` refuses a
  missing `newdata` variable and a new factor level before each
  `model.frame()` of the `newdata` path (`pred_design()`, its
  random-effect part and `mm_member_designs()`). It uses the lookup of
  `model.frame()`, so it never refuses a variable that `model.frame()`
  would find. `check_frame_variables()` does the same in
  `assemble_frame()` for an unknown variable, a list used as a whole
  term, and an environment object of the wrong length. A `data` that is
  not a data frame, an environment or a named list is refused before
  that. `hyp_parse_all()` refuses a name in backticks. frmtmb.sample
  refuses an exact variable name that the draws lack
  (`draws_subset_variable()`) before posterior sees it, which covers
  `posterior_interval()`, `as.matrix()`, `rhat()` and `neff_ratio()`.
- **Subclass** (review MINOR 1). `frmtmb_family()` and
  `frmtmb_structure()` record the package of their caller as the
  attribute `frmtmb_package`, and `frm_family_package()` reads it. An
  attribute rather than a list element keeps the names of the family,
  and the partial matches of `$` on them, as they were.
  `structure_gate()` and 40 core refusals about a declaration of a
  family pass `package = frm_family_package(<family>)`. The sites and
  the exclusions are in `family-tag-plan.txt`, made from parse data.
  frmtmb.eam's non-decision-time seam records the calling namespace in
  `ndt_bound()` and `ndt_bound_pending()`, as an attribute of the bound
  record. `ndt_apply()` looks the namespace up on its refusal branch
  only, so the per-trial path does no extra work.
- **insight** (review MINOR 2). `model.matrix.frmtmb_fit()` adds
  `simpleError` after `frmtmb_error` to its refusal.
  `dev/conditions-rev-interop2.R lane` now differs from the base log
  only in the printed class vector (`interop2-lane.txt`).
- **Census** (review MINOR 3). The template is a file of its own,
  `dev/conditions-census-template.R`. It counts exemptions
  (`census_compare()`), scans `R/` and `inst/` recursively for
  `.R .r .S .s .q`, and reads the `what` of `do.call()` by name in any
  position. It flags `get`, `get0`, `mget`, `match.fun`,
  `getExportedValue` and `dynGet` of a condition function's name;
  `call`, `as.call`, `as.name` and `as.symbol` of one; text given to
  `str2lang`, `str2expression` or `parse` that calls one; and the
  raisers `.Deprecated`, `.Defunct`, `.signalSimpleWarning`,
  `.dfltStop` and `.dfltWarn`. A source tree is recognized by its
  DESCRIPTION and by having no `Meta/`. The anchor file is then
  asserted, so a rename fails. The real trees needed no new exemption.

### Message text that changed

Every refusal above has new text, listed in `NEWS.md` of frmtmb and of
each extension. No other text changed. The review's printed-output
comparison of 61 sites included none of these sites.

### Seen failing before the fix

Every new or changed test was run on the build from before this round
(`before-*.txt`) and on the new build (`after-*.txt`). Before:
`test-conditions.R` failed in 7 of its 14 blocks (104 failures, 2
errors); `test-sugar.R` (insight), `test-interop.R` (getME),
`test-nl-lexical.R` (2 blocks), `test-tabular-inputs.R` and
`test-arg-refusal.R` errored; the punch block of each extension's
`test-conditions.R` failed (coupling 1, eam 4, latent 2, learn 4, ode
2, sample 2, spline 1 failures); and the changed expectation errored in
`test-cross-spectrum.R`, `test-family.R`, `test-gddm-family.R`, both
`test-surface.R`, `test-hmm.R`, `test-rlddm-ndt.R`,
`test-royston-parmar.R`, `test-draws-methods.R` (2 blocks) and
`test-sample-direct.R`. After: every one of these files has 0 failures
and 0 errors.

The immediate-warning test runs R in a child process, because
testthat muffles every warning in a test. On the old build the child
printed the deferred `Warning message:` form. For the census, the old
census passed or skipped on each plant below
(`dev/conditions-rev-log/census-rerun.txt`, `census2.txt`).

### User-reachable refusals, probed

`dev/conditions-punch-probe.R`, seed 20260917, one probe per site of the
review's lists, including all nine prior densities and both
exponential cases:

<!-- BEGIN GENERATED: punch probe -->
```
before: probes 65  frmtmb_error 1  not frmtmb_error 64  expected subclass 1 
after:  probes 65  frmtmb_error 65  not frmtmb_error 0  expected subclass 65 
```
<!-- END GENERATED: punch probe -->

"expected subclass" is the extension's subclass for an extension site,
and no subclass for a core site.

### Census plants

`dev/conditions-punch-census.R`:

<!-- BEGIN GENERATED: punch census -->
```
core_control                       pass 10 fail 0 error 0 skip 0  -> PASSES as required
eam_control                        pass 10 fail 0 error 0 skip 0  -> PASSES as required
core_exempt_second_site_fit        pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_exempt_second_site_helper     pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_exempt_owner_redefined        pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_lowercase_r                   pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_subdir_windows                pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_subdir_unix                   pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_ext_S                         pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_ext_q                         pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_ext_lower_s                   pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_inst_subdir_S                 pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_docall_what_first             pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_docall_what_late              pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_get                           pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_match_fun                     pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_call_string                   pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_str2lang                      pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_parse_text                    pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_deprecated                    pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_signal_simple_warning         pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_sentinel_renamed              pass  9 fail 1 error 0 skip 0  -> FAILS  as required
eam_sentinel_renamed               pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_bare                          pass  9 fail 1 error 0 skip 0  -> FAILS  as required
eam_base_message                   pass  9 fail 1 error 0 skip 0  -> FAILS  as required
core_value_map                     pass  9 fail 1 error 0 skip 0  -> FAILS  as required
```
<!-- END GENERATED: punch census -->

### The sweep, and the subclass before and after

Both use the review's runner, `dev/conditions-rev-sweep-onefile.R`
through `dev/conditions-rev-runset.sh`, and one summarizer,
`dev/conditions-punch-summary.R`. "Before" is the review's sweep of the
lane build (`dev/conditions-rev-log/sweep2-lane`); the summarizer gives
the review's 606, 504 and 102 from it. "After" is
`dev/conditions-punch-log/sweep-lane`.

<!-- BEGIN GENERATED: punch sweep -->
```
before (summary-before.txt):
== extension subclass, in each extension's own suite ==
(the construction of dev/conditions-rev-subclass-sweep.R)
                caught with_own without pct_with
frmtmb.coupling     32       28       4     87.5
frmtmb.eam         192      162      30     84.4
frmtmb.latent       55       42      13     76.4
frmtmb.learn        56       34      22     60.7
frmtmb.ode          78       76       2     97.4
frmtmb.sample      138      108      30     78.3
frmtmb.spline       55       54       1     98.2
ALL                606      504     102     83.2

after (summary-after.txt):
== suites, one file per R process ==
test files on disk 247  result lines 247 
frmtmb           files 142 of 142 (no result line 0) blocks 1559 PASS 11860 FAIL 0 ERROR 0 SKIP 1
frmtmb.coupling  files   9 of   9 (no result line 0) blocks   66 PASS   446 FAIL 0 ERROR 0 SKIP 5
frmtmb.eam       files  26 of  26 (no result line 0) blocks  276 PASS  1651 FAIL 0 ERROR 0 SKIP 3
frmtmb.latent    files   9 of   9 (no result line 0) blocks   70 PASS   348 FAIL 0 ERROR 0 SKIP 2
frmtmb.learn     files  14 of  14 (no result line 0) blocks   99 PASS   478 FAIL 0 ERROR 0 SKIP 2
frmtmb.ode       files  10 of  10 (no result line 0) blocks  118 PASS   342 FAIL 0 ERROR 0 SKIP 1
frmtmb.sample    files  23 of  23 (no result line 0) blocks  232 PASS  1711 FAIL 0 ERROR 0 SKIP 1
frmtmb.spline    files  14 of  14 (no result line 0) blocks   93 PASS   449 FAIL 0 ERROR 0 SKIP 1
start 11:20:13  end 11:38:00 
frmtmb.sample logs 23  lines naming compileCode 0 

== expect_error() conditions ==
                caught frmtmb_error other
frmtmb            1444         1432    12
frmtmb.coupling     33           33     0
frmtmb.eam         193          193     0
frmtmb.latent       55           55     0
frmtmb.learn        56           56     0
frmtmb.ode          78           78     0
frmtmb.sample      142          142     0
frmtmb.spline       56           56     0
ALL               2057         2045    12
```
<!-- END GENERATED: punch sweep -->

The 45 extension-suite refusals without their own subclass, each by the
rule on `?frmtmb-conditions` (the messages are in
`summary-after.txt`):

- frmtmb.eam, 5: a `vint()` value that is not an integer and a
  `newdata` without an addition-term column, which are core checks of
  data for any family; a constant outside the range of a link (two); a
  prior on `bs`. None is about a declaration of the family.
- frmtmb.latent, 1: `quadrature = TRUE` with a block that is not a
  scalar intercept, a core limit of the quadrature.
- frmtmb.learn, 6: `simulate(newdata =)`, a core argument; a constant
  `ndt` outside its link; and four that frmtmb.eam functions raise
  for rules of frmtmb.eam, so they are `frmtmb_eam_error`:
  `frmtmb.eam::ndt_time()` called by the test (two), the `ndt_group()`
  coercion that frmtmb.eam registers, and the frame check of
  frmtmb.eam that no family read `ndt_group()`.
- frmtmb.ode, 2: one message, core's wrapper of an R error in a
  nonlinear body ("The nonlinear formula body could not be evaluated").
- frmtmb.sample, 31: refusals of core functions that the sample suite
  calls, such as the maximum-likelihood `loo()`, `waic()` and
  `bayes_R2()` stubs, argument spellings, an object built from a
  formula only, and unknown parameters; and two `hmm()` refusals that
  are `frmtmb_latent_error`. frmtmb.sample raises none of the 31.

### Unclassed errors that remain

The `expect_error()` conditions of the sweep that are not a
`frmtmb_error` went from 29 to 12 (listed at the end of
`summary-after.txt`). frmtmb raises none of them:

- 4 from brms itself, in the agreement tests.
- 2 from `stan_pars_from_fit()`, a test helper in
  `tests/testthat/helper-brms.R`, and 1 deliberate absent case in
  `test-conditions.R`.
- 1 from `stats::binomial(link = "softit")`.
- 4 "unused argument" (`priors =`, `lower =`, `upper =`), left to R as
  decided.

Not in the sweep, and left: `log_lik()` on a maximum-likelihood fit
still stops at "no applicable method". Its generic belongs to
frmtmb.sample, and core has no `log_lik()` refusal near its `loo()` and
`waic()` stubs to extend. `frm(data = d)` without a formula stays R's
"argument is missing".

### NIT 3: the sweep counts

The lane's runner caught 2,054 `expect_error()` conditions: the
`CAUGHT` lines of `dev/conditions-log/lane-core` and `lane-ext` sum to
1,441 and 613. The findings printed 2,032 because
`dev/conditions-sweep-summary.R` read the `.tsv` files back with
`read.delim(quote = "\"")`. `write.table()` had written the quotes
inside messages of `test-prior-compat.R` as `\"`, which `read.delim()`
does not undo, so rows merged. `dev/conditions-punch-nit3.R`
(`nit3.txt`) finds that file the only one that differs, 23 rows against
45. The review's rows also number 2,054; its 2,053 is its join with the
base build, which drops one row found only on the lane. So 2,032 and
2,003 were an undercount of the same construction. The counts in this
section come from the review's `.rds` rows, which are never parsed back
from text.

### Checks

<!-- BEGIN GENERATED: punch checks -->
```
check-core.txt:
* checking HTML version of manual ... [21s] NOTE
Status: 1 NOTE

check-sample.txt:
Status: OK

```
<!-- END GENERATED: punch checks -->

### Needs the user

- The floors `frmtmb (>= 0.59.0)` of the seven extensions must move to
  the release that exports `frm_match_arg()` and
  `frm_family_package()`, as they must for `frm_stop()`.
- The subclass rule is now "the package that wrote the refusal", with
  the exceptions above. frmtmb.sample's 31 are core refusals and stay
  unclassed by design. If the user wants every refusal met through an
  extension's object to carry that extension's subclass, the mechanism
  would have to be dynamic (the families of the model in hand) rather
  than per site.
- `log_lik(<frmtmb_fit>)` could get a refusal in frmtmb.sample like the
  `loo()` stub; not done, as instructed.

## Punch round 2

Review: the "Recheck, round 1" section of
`dev/reviews/20260917-conditions.md`. Logs are
`dev/conditions-punch-log/p2-*`. Assembled by
`dev/conditions-punch2-findings.R`.

### R1. A variable in the parent of an environment `data`: fixed

`check_frame_variables()` now searches an environment `data` with
inheritance and does not consult the formula environment, which is
what `model.frame()` does when it evaluates in that environment. A data
frame or list `data` is unchanged. The test "the frame check searches
an environment data as model.frame()" in `test-conditions.R` fits the
reviewer's construction and compares its log-likelihood with the same
model on a data frame, bitwise. It errored on the build before the fix
(`p2-before-test-conditions.txt`).

### R2. The frame check named the wrong problem: fixed

- `sar()` and `fcor()`, brms autocorrelation terms that frmtmb does not
  implement, are refused by name when the formula is parsed
  (`R/parse.R`, beside the `mm()` and `mo()` term refusals), before the
  frame check sees their `data2` argument. frmtmb had no refusal for
  them to reuse: base reached `model.frame()` and "could not find
  function".
- A `.` in a formula is refused as unsupported, with the advice to
  write the predictors out. frmtmb does not expand `.`; base refused it
  through `model.frame()` with "'.' in formula and no 'data' argument".

`dev/conditions-punch2-probe.R`, before and after:

<!-- BEGIN GENERATED: punch2 probe -->
```
before:
env_parent  frmtmb_error/error/condition | The model uses `xq`, which is not a column of `data` and not an object that R finds from the formula
sar_data2   frmtmb_error/error/condition | The model uses `Wm`, which is not a column of `data` and not an object that R finds from the formula
fcor_data2  frmtmb_error/error/condition | The model uses `Vm`, which is not a column of `data` and not an object that R finds from the formula
dot         frmtmb_error/error/condition | The model uses `.`, which is not a column of `data` and not an object that R finds from the formula.
after:
env_parent  fits
sar_data2   frmtmb_error/error/condition | sar() is a brms autocorrelation term that frmtmb does not support: sar(Wm). The supported residual c
fcor_data2  frmtmb_error/error/condition | fcor() is a brms autocorrelation term that frmtmb does not support: fcor(Vm). The supported residual
dot         frmtmb_error/error/condition | A formula with `.` is not supported: frmtmb does not expand `.` into the columns of `data`. Write th
```
<!-- END GENERATED: punch2 probe -->

### R3. Census over-reach: narrowed, one shape recorded

In `dev/conditions-census-template.R`:

- Names a top-level expression binds for itself (assignment targets,
  formals, `for` variables) are local objects there. Such a name is
  not flagged as a value, and a local function named like a lookup or
  a builder is not treated as one. A call of `stop()`, `warning()` or
  `message()` is flagged whatever is bound; the planted
  `function(message) stop(message)` checks that.
- The left side of an assignment is not a value.
- `as.name()` and `as.symbol()` are no longer builders. A symbol raises
  nothing; they are flagged where they become a call, inside
  `as.call()` (planted: `as.call(list(as.name("stop"), "p"))`).
- Code text that parses to formulas only is not flagged.
- Recorded as a known false positive, not narrowed:
  `get("warning", envir = e)`. With the default `inherits = TRUE` it can
  reach base's function, so the rule cannot tell the two apart from the
  call. It fails loudly and takes an exemption.

The reviewer's `dev/conditions-rev-census-overreach.R` on the new
census:

<!-- BEGIN GENERATED: punch2 overreach -->
```
local_var_message      clean
list_field_message     clean
local_get_fn           clean
get_column_env         FLAGGED: <get> warning
identical_as_name      clean
call_inspection        clean
switch_labels          clean
local_fn_named_call    clean
parse_formula_text     clean
mget_names             clean
getOption_warn         clean
formals_default_msg    clean
R6_like                clean
dollar_call            clean
```
<!-- END GENERATED: punch2 overreach -->

The census test itself now asserts the six legitimate shapes are
clean. All eight real trees pass with no new exemption, and the 26
plants of `dev/conditions-punch-census.R` give the same result as in
round 1:

<!-- BEGIN GENERATED: punch2 census -->
```
core_control                       pass 11 fail 0 error 0 skip 0  -> PASSES as required
eam_control                        pass 11 fail 0 error 0 skip 0  -> PASSES as required
core_exempt_second_site_fit        pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_exempt_second_site_helper     pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_exempt_owner_redefined        pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_lowercase_r                   pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_subdir_windows                pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_subdir_unix                   pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_ext_S                         pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_ext_q                         pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_ext_lower_s                   pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_inst_subdir_S                 pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_docall_what_first             pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_docall_what_late              pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_get                           pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_match_fun                     pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_call_string                   pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_str2lang                      pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_parse_text                    pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_deprecated                    pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_signal_simple_warning         pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_sentinel_renamed              pass 10 fail 1 error 0 skip 0  -> FAILS  as required
eam_sentinel_renamed               pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_bare                          pass 10 fail 1 error 0 skip 0  -> FAILS  as required
eam_base_message                   pass 10 fail 1 error 0 skip 0  -> FAILS  as required
core_value_map                     pass 10 fail 1 error 0 skip 0  -> FAILS  as required
```
<!-- END GENERATED: punch2 census -->

### Filed, not fixed: newdata without a grouping column

`predict(fit, newdata = <no grouping column>, allow_new_levels = TRUE)`
stops at base R's "object 'g' not found" on the base build and on this
lane, for a random intercept, a random slope and a binomial mixed model
(`dev/conditions-rev-log/newdata-brms-cmp.txt`). brms's
`validate_newdata()` accepts it. `check_newdata_frame()` guards the
fixed-effect `model.frame()`; the grouping lookup in `pred_design()`
evaluates the bar's right side with `eval()` and is not guarded. This is
a user-reachable unclassed error, and with `allow_new_levels = TRUE` it
is arguably a missing feature too (brms predicts at the population
level). It belongs to a later item.

### Suites

Every test file of all eight packages, one R process each, on the
build with these changes (`p2-sweep`, `p2-summary.txt`):

<!-- BEGIN GENERATED: punch2 suites -->
```
== suites, one file per R process ==
test files on disk 247  result lines 247 
frmtmb           files 142 of 142 (no result line 0) blocks 1561 PASS 11865 FAIL 0 ERROR 0 SKIP 1
frmtmb.coupling  files   9 of   9 (no result line 0) blocks   66 PASS   447 FAIL 0 ERROR 0 SKIP 5
frmtmb.eam       files  26 of  26 (no result line 0) blocks  276 PASS  1652 FAIL 0 ERROR 0 SKIP 3
frmtmb.latent    files   9 of   9 (no result line 0) blocks   70 PASS   349 FAIL 0 ERROR 0 SKIP 2
frmtmb.learn     files  14 of  14 (no result line 0) blocks   99 PASS   479 FAIL 0 ERROR 0 SKIP 2
frmtmb.ode       files  10 of  10 (no result line 0) blocks  118 PASS   343 FAIL 0 ERROR 0 SKIP 1
frmtmb.sample    files  23 of  23 (no result line 0) blocks  232 PASS  1712 FAIL 0 ERROR 0 SKIP 1
frmtmb.spline    files  14 of  14 (no result line 0) blocks   93 PASS   450 FAIL 0 ERROR 0 SKIP 1
start 12:21:05  end 12:35:09 
frmtmb.sample logs 23  lines naming compileCode 0 
```
<!-- END GENERATED: punch2 suites -->

No export and no Rd changed, so `R CMD check` was not rerun. The
extensions were roxygenised in a fresh process after core was
installed; their NAMESPACE files did not change.
