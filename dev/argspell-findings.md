# Lane argspell: plan item 2.5e

The argument-name audit against brms, the spelling changes that came out
of it, and the refusal that stops the class rather than the instance.

Worktree `frmtmb-wt-argspell`, branch `wt-argspell`, base `64237e2`.
Private library `C:/Users/adf44/source/r/argspell-lib`; the base to
measure against is the shared reference library
`C:/Users/adf44/source/r/rellib-r3`, which holds frmtmb 0.57.0 and
frmtmb.sample 0.5.0.

## 1. The defect, reproduced

`dev/argspell-defect.R`, run against the reference library, seeds 2505
for the design and 25051 for the response draw. One gaussian fit,
`y ~ x + (1 | g)`, 120 rows over 8 groups, `sd_g__Intercept = 1.2`
against `sigma = 0.4`, so the conditional and population answers are far
apart on purpose.

    frmtmb 0.57.0 from C:/Users/adf44/source/r/rellib-r3

    fitted(fit)[1]                   = 1.291651
    fitted(fit, re_formula = NA)[1]  = 1.291651
    predict(re.form = NA)[1]         = 1.954262
    identical(fitted, swallowed)     = TRUE
    max |fitted - population|        = 3.065438
    misspelling swallowed silently   = TRUE

`identical()` rather than a printed difference: `%.6f` renders 1.8e-11
as `0.000000`, and what is claimed here is that the two are the same
object, not that they are close.

The task statement's own figures are 1.31023 and 1.74231 on a different
fixture. Those are not reproduced here and are not contradicted; this is
a separate construction of the same defect, and its seeds and script are
recorded so it can be re-derived.

## 2. The audit

`dev/argspell-audit.R` reads the registered method table with
`parseNamespaceFile()` and the real `formals()` of each method, for
frmtmb and frmtmb.sample, and pairs each with the `brmsfit` method of
the same generic. Full output in `dev/argspell-formals.csv`;
`dev/argspell-dots.csv` carries the dots analysis. A grep over NAMESPACE
would have missed the multi-name `export(...)` block in
`R/sampling-api.R`, which is where the internals shared with extensions
are declared, so nothing here is grepped.

### 2.1 What the diff showed, by class

**One method was wrong, and it was the one the plan named.**
`fitted.frmtmb_fit(object, ...)` against
`fitted.brmsfit(object, newdata, re_formula, scale, resp, dpar, nlpar,
ndraws, draw_ids, sort, summary, robust, probs, ...)`: twelve
brms arguments, none of them present, all of them swallowed.

**The seam the plan described is real and was in one signature.**
`predict.frmtmb_fit` paired `re.form`, lme4's spelling, with
`allow_new_levels`, brms's, and accepted `allow.new.levels` out of its
dots as well. `conditional_effects()` already took `re_formula`, and
refused `re.form` with a message that explained the two-dialect rule.
`frmtmb.sample`'s `posterior_predict()` took both, as brms does.

**Most remaining diffs are one class and are not defects.** For the
great majority of the methods with a brms counterpart, every brms-only
argument is drawn from `ndraws`, `draw_ids`, `sort`, `summary`,
`robust`, `probs` and `pars`: brms summarizes draws and a maximum
likelihood fit has one estimate. `dev/argspell-formals.csv` has the
per-method diff; the count is not restated here because it would be a
typed number with a file behind it that already says the same thing.
Those arguments are refused by name with the reason rather than
adopted; see 3.3 and 3.4.

**Twelve methods have no brms counterpart at all**, and are not
expected to. From `dev/argspell-audit.R`'s own output:

    find_formula.frmtmb_fit      find_random.frmtmb_fit
    find_statistic.frmtmb_fit    get_coef.frmtmb_fit
    get_parameters.frmtmb_fit    get_predict.frmtmb_fit
    get_varcov.frmtmb_fit        get_vcov.frmtmb_fit
    getME.frmtmb_fit             link_function.frmtmb_fit
    link_inverse.frmtmb_fit      set_coef.frmtmb_fit

Eleven belong to insight and marginaleffects, and `getME()` is lme4's.
Their names come from the package that owns the generic, which is the
same rule applied. `emm_basis` and `recover_data` are NOT in this list,
because brms defines `emm_basis.brmsfit` and `recover_data.brmsfit`
itself; they are exempt from the dots refusal for a different reason
(section 4), not for want of a counterpart.

### 2.2 The argument the audit settled by looking rather than by memory

Where brms carries BOTH `re_formula` and `re.form` is a fact about brms,
not a judgment, and it decides where the alias may stay. Measured from
`brms` 2.23.0's own formals:

    posterior_epred.brmsfit      re_formula, re.form
    posterior_linpred.brmsfit    re_formula, re.form
    posterior_predict.brmsfit    re_formula, re.form
    predictive_error.brmsfit     re_formula, re.form
    predictive_interval.brmsfit  neither (prob, ...)
    pp_check.brmsfit             neither; forwards to
                                 prepare_predictions(), whose formals
                                 carry re_formula alone

So the alias stays on four methods and was removed from two,
`predictive_interval()` and `pp_check()`, where frmtmb.sample had it
without brms precedent. `tests/testthat/test-draws-spellings.R` now
COMPARES against `formals(<generic>.brmsfit)` rather than restating a
list, so the two cannot drift apart.

## 3. What changed

### 3.1 Spellings: brms wins, the lme4 name is dropped

Was, then is. Every "now" name is brms's.

* `predict.frmtmb_fit()`: `re.form`, and `allow.new.levels` accepted
  out of the dots, become `re_formula` and `allow_new_levels`. Both
  lme4 names are refused by name.
* `simulate.frmtmb_fit()`: `re.form` becomes `re_formula`.
* `fitted.frmtmb_fit()`: `...` alone becomes `newdata`, `re_formula`,
  `scale`, `resp`, `dpar`.
* `frm_bootstrap()`, `dharma_residuals()`: `re.form` becomes
  `re_formula`.
* `pp_check.frmtmb_fit()`, `pp_check.frmtmb_draws()`,
  `predictive_interval.frmtmb_draws()`: carried `re_formula` AND
  `re.form`; the alias is dropped, because brms's counterparts carry
  neither.
* `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`,
  `predictive_error()` on draws: unchanged, both spellings, because
  brms carries both.

The rename reached five other extension packages, which call
`predict(re.form = )` into the core and spell their own arguments the
same way: `frmtmb.coupling` (`frm_coherence()`, `frm_phase()`),
`frmtmb.spline` (`frm_curve()` and the derivative and feature
surfaces), `frmtmb.latent` (two refusal messages), `frmtmb.eam` (one
doc line). They are installed, not test-run, in this lane; see 6.

`fitted()` also changed a MESSAGE, because it delegates now: a family
with no mean used to answer "fitted() is not defined for family
'multinomial'" and now answers the family's own "family 'multinomial'
declares no mean: it has no dpar named mu and no post$mean_fn ... Ask
for type = "link" or a dpar by name", which says what to do instead.
`test-multinomial.R` pinned the old substring and was updated;
`test-mean-fn.R` already allowed either.

### 3.2 `fitted()` is now one call to `predict()`

`fitted.frmtmb_fit()` had its own copy of the structured-family,
ordinal and categorical branches, and the documented identity
`fitted() == predict(type = "response")` was two bodies that agreed
rather than one body. It is now:

    predict(object, newdata = newdata,
            type = if (scale == "response") "response" else "link",
            dpar = dpar, resp = resp, re_formula = re_formula)

with one guard kept in front of it: `fitted()` on a multivariate fit
with no `resp` still refuses, where `predict()` defaults to the first
response.

### 3.3 The refusal, which is the part that prevents the class

`frm_check_dots(...)` in `R/utils.R`, exported to siblings through the
`@rawNamespace` block in `R/sampling-api.R`. It reads the CALLER's name
and formals out of the calling frame, so the call site is the same one
line everywhere and cannot drift out of step with a renamed argument. It
reports the generic rather than the method, because the caller typed
`fitted()` and not `fitted.frmtmb_fit()`, and it offers the nearest
formal when the edit distance is at most a quarter of the name's length.

Three refusal modes:

* an unknown name: `fitted() has no argument 'ndraws'` plus the list of
  what it does take;
* an argument with NO name, which a method with only `object` and `...`
  has nowhere to put;
* a name in `.unsupported`: a REAL brms argument, or a retired lme4
  spelling, refused with the reason and with where the argument does
  work.

`.allow` exists for the two methods that forward their dots on to
somebody else: `residuals(type = "osa")` checks them against
`formals(TMB::oneStepPredict)`, and `pp_check()` forwards everything to
bayesplot's `ppc_*` but still refuses `re.form` itself first, because
bayesplot would have accepted that name and done nothing with it.

### 3.4 Arguments accepted-but-not-implemented are now refused with the reason

`fitted()`: `ndraws`, `draw_ids`, `sort`, `summary`, `robust`, `probs`,
`nlpar`. `predict()`: the same seven plus `transform`, and the two
retired lme4 spellings. `residuals()`: `newdata`, `re_formula`,
`method`, `resp`, and the five draws arguments.

`residuals(newdata = )` and `residuals(re_formula = )` are the two that
are a missing FEATURE rather than a refusal of principle, and the
message says so: a residual needs an observed response, and this package
evaluates residuals on the fitted rows only. They are filed here rather
than implemented, because implementing them is a second change and this
lane's job was to stop the silence.

## 4. What I did NOT change, and why

* **`predict(type = )` stays, and its default stays `"link"`.** brms's
  `predict.brmsfit()` has no `type` and always returns the response
  scale; `frmtmb`'s default is the link scale and is documented as a
  deliberate divergence in `?frmtmb-scales`. Changing the default return
  SCALE of `predict()` is a different decision from changing an argument
  NAME, it is not what item 2.5e asks for, and it would move numbers in
  every downstream extension. `transform` is refused with a message
  pointing at `type`.
* **`residuals()` does not gain `newdata`/`re_formula`.** Feature, not
  spelling; refused explicitly instead. See 3.4.
* **The draws methods do not gain `dpar`, `nlpar`, `draw_ids` or
  `sort`.** Same reasoning. `posterior_epred.frmtmb_draws` lacks `dpar`
  where brms has it; that is recorded here as an open gap rather than
  closed, because adding it is a behavior change that wants its own
  numbers.
* **Nine interop methods keep swallowing dots.** This was thirteen and
  the count was both too wide and misstated: the prose said thirteen
  where the generated block said twelve, because
  `recover_data.frmtmb_fit` READS its dots, so its entry did nothing.
  Section 10.3 has the re-measurement, the three that were closed, and
  the two that had to be REOPENED after a run. The first pass's
  evidence, kept because it is what the emmeans entries still rest on,
  was `dev/argspell-exempt-evidence.R` walking the installed namespaces
  for call sites that forward dots into each generic:

      emmeans: 3 call sites forwarding dots into emm_basis(
      emmeans: 14 call sites forwarding dots into recover_data(
      insight: 7 call sites forwarding dots into find_formula(
      insight: 12 call sites forwarding dots into get_parameters(
      marginaleffects: 1 call sites forwarding dots into get_predict(

  Two of them name the argument that would break:
  `insight:::.n_parameters_component` calls
  `get_parameters(x, component = component, ...)` and
  `insight:::.n_parameters_effects` calls
  `get_parameters(x, effects = effects, ...)`, and
  `get_parameters.frmtmb_fit(x, ...)` has neither formal. A refusal
  there breaks insight, not a user's typo. The thirteen are listed by
  name in `tests/testthat/test-arg-refusal.R` and in
  `dev/argspell-insert.R`.
* **A method whose whole body is an unconditional refusal keeps its
  dots, and is exempt by SHAPE rather than by name.** Guarding it is
  worse than not guarding it: every call already errors, and the guard
  only replaces a message that says what to do instead with one that
  says the argument is unknown. On one method it made the message
  WRONG. `conditional_effects(frm_multiple_result, "x")` passes the
  effect positionally, and with the guard in front the answer became
  "conditional_effects() was given 1 argument with no name, which it
  has nowhere to put. It takes: x" instead of "no pooled version for a
  frm_multiple() result". That was caught by the suite, not by review;
  `dev/argspell-unstop.R` removed the guard from the 38 methods of that
  shape, and `ar_refuses_always()` in the test recognises the shape so
  the invariant stays total without a name list. The near miss is
  asserted too: a refusal reached only on a BRANCH is still a
  swallower.
* **`getME(name = )` keeps lme4's spelling.** `getME()` IS lme4's
  function; brms has no counterpart, so there is no disagreement for
  brms to break.
* **No version number is chosen.** NEWS carries
  `# frmtmb (development version)`, and the code comments say "no longer
  accepted" rather than naming a release.

## 5. Evidence

### 5.1 The test was seen failing against 0.57.0

`tests/testthat/test-arg-refusal.R` against the reference library, with
no part of the change installed (`dev/argspell-run1.R . arg-refusal base`):

    frmtmb from: C:/Users/adf44/source/r/rellib-r3
    version    : 0.57.0
    ARGSPELL arg-refusal pass 7 fail 22 error 2 skip 0 in 2 s

and against the change installed (the file grew by three assertions
after the first run, which is why the pass count is 47 and not 44):

    frmtmb from: C:/Users/adf44/source/r/argspell-lib
    ARGSPELL arg-refusal pass 47 fail 0 error 0 skip 0 in 2 s

The 22 failures are behavioural, not "the symbol does not exist": the
first is `fitted(fit)` and `fitted(fit, re_formula = NA)` compared and
found equal.

### 5.2 The structural guard, constructed in its absent case first

The last block of `test-arg-refusal.R` asserts that NO registered S3
method has a `...` its body never mentions. The detector is exercised on
a synthetic swallower and a synthetic forwarder in the test above it,
because a guard whose positive branch has never been observed to fire is
not a guard:

    expect_true(ar_swallows(function(object, ...) object))
    expect_false(ar_swallows(function(object, ...) list(...)))
    expect_false(ar_swallows(function(object) object))
    expect_false(ar_swallows(function(object, ...) stop("no")))
    expect_true(ar_swallows(function(object, ...) {
      if (is.null(object)) stop("no")
      object
    }))

The last pair is the one that matters after the unconditional-refusal
exemption went in: a refusal on a BRANCH leaves the other branches
swallowing, so it must still count as a swallower.

Against 0.57.0 the detector reports 145 methods across the two
packages.

### 5.3 Two runner defects, neither of them the change

**The test environment.** The first whole-suite run reported errors like
`could not find function "autocor_natural"` in nine files. That was the
RUNNER: `testthat::test_dir()` without `package =` does not put the
package namespace behind the test environment, so every test that calls
an unexported helper by its bare name errors. Fixed in
`dev/argspell-run1.R`; `test-autocor.R` went from 9 errors to
`pass 191 fail 0 error 0`, which is its baseline count exactly. Recorded
because a runner defect that only shows up on SOME files is exactly the
shape that gets read as a regression.

**PowerShell and a native command's stderr.** `dev/argspell-check.ps1`
ran with `$ErrorActionPreference = "Stop"`. `R CMD build` writes its
progress to stderr, PowerShell 5.1 wraps a native command's stderr in an
ErrorRecord, and the script therefore aborted after a SUCCESSFUL build
and never ran the check. The tarball was on disk and the log ended at
`* building 'frmtmb_0.57.0.tar.gz'`, which reads as a build failure and
is not one. Changed to `"Continue"`, with the reason in the script.

## 6. The extensions, and an ordering constraint for consolidation

The five extension packages other than `frmtmb.sample` had `re.form`
renamed in R sources, tests and vignettes, because they call
`predict(re.form = )` into the core and spell their own arguments the
same way. All seven extensions were roxygenised and INSTALLED into the
lane library, so the rename is known to build and document; their own
suites were not run, and their own S3 methods did NOT get the dots
refusal. Both are follow-on work.

**The ordering constraint is real and it cost a green run here.**
`frmtmb.sample`'s `test-simulators.R` asserts the text of a refusal
that lives in `frmtmb.latent`. Against the reference library's
`frmtmb.latent` 0.3.0, which still says `re.form`, that block errored:

    test-simulators.R  pass 48 fail 0 error 1 skip 0

With `frmtmb.latent` rebuilt from this worktree it is its baseline
count:

    ARGSPELL ^simulators$ pass 50 fail 0 error 0 skip 0 in 3 s

So a consolidating session that installs core and `frmtmb.sample` and
stops will see one failure that is not a defect. Every extension in
this worktree has to be reinstalled together.

## 7. Open, filed rather than fixed

* `posterior_epred.frmtmb_draws()` has no `dpar`, where
  `posterior_epred.brmsfit()` does. `posterior_linpred()` here already
  has it, so the gap is one method wide.
* `residuals()` has no `newdata` or `re_formula`. Refused with the
  reason; see 3.4.
* `anova.frmtmb_fit()` and `anova.frmtmb_multiple()` select their extra
  models with `Filter(inherits, list(...))`, so a stray NAMED argument
  is dropped rather than refused. That is brms's and stats's own
  behavior for a variadic model list, and refusing it needs a rule for
  telling a model apart from a typo. Left as is, recorded here.
* The five extensions other than `frmtmb.sample` still have methods
  that swallow their dots.

## 8. Generated counts

Emitted by the scripts named beside them and pasted verbatim. None of
these numbers is typed.

### 8.1 `Rscript dev/argspell-report.R`, against the installed change

Regenerated after review round 2. The detector reads the parse tree, so
a string literal cannot inflate the swallower count, and `..N` is
matched by pattern rather than enumerated.

    ---- GENERATED: method counts ----
    frmtmb
      registered S3 methods (deduped) : 123
      taking `...`                    : 114
      with `...` never read           : 9
      of those, exempt by contract    : 9
      with a brms counterpart         : 84
    frmtmb.sample
      registered S3 methods (deduped) : 68
      taking `...`                    : 65
      with `...` never read           : 0
      of those, exempt by contract    : 0
      with a brms counterpart         : 68
    
    METHODS STILL SWALLOWING THEIR DOTS: 0
    frm_check_dots(...) call sites: frmtmb 72, frmtmb.sample 36
    ---- END GENERATED: method counts ----

    ---- GENERATED: re.form / re_formula surface ----
    predict.frmtmb_fit               re_formula TRUE re.form FALSE
    fitted.frmtmb_fit                re_formula TRUE re.form FALSE
    simulate.frmtmb_fit              re_formula TRUE re.form FALSE
    residuals.frmtmb_fit             re_formula FALSE re.form FALSE
    conditional_effects.frmtmb_fit   re_formula TRUE re.form FALSE
    pp_check.frmtmb_fit              re_formula TRUE re.form FALSE
    frm_bootstrap                    re_formula TRUE re.form FALSE
    dharma_residuals                 re_formula TRUE re.form FALSE
    posterior_epred.frmtmb_draws     ours TRUE  brms declares TRUE  accepts TRUE
    posterior_linpred.frmtmb_draws   ours TRUE  brms declares TRUE  accepts TRUE
    posterior_predict.frmtmb_draws   ours TRUE  brms declares TRUE  accepts TRUE
    predictive_error.frmtmb_draws    ours TRUE  brms declares TRUE  accepts TRUE
    predictive_interval.frmtmb_draws ours TRUE  brms declares FALSE accepts TRUE
    pp_check.frmtmb_draws            ours FALSE brms declares FALSE accepts TRUE
    ---- END GENERATED: re.form / re_formula surface ----

`pp_check.frmtmb_draws` is the one row where `ours` and `brms accepts`
disagree, and 11.1 is why: brms honors the alias there while warning
that it ignored it, and a warn-then-honor path is not a decision to
copy. The column is printed rather than hidden.

The formals diff for the prediction family is in
`dev/argspell-report-log.txt` in full. What remains as a brms-only
argument is the draws class, `ndraws`, `draw_ids`, `sort`, `summary`,
`robust`, `probs`, `pars`, plus `nlpar`, and every method that can be
handed one refuses it by NAME with the reason (10.5).


### 8.1b The exemption evidence, `Rscript dev/argspell-exempt2.R`

    ---- GENERATED: exemption evidence, re-measured ----
    A call site COUNTS when it passes a name the method has no formal for.
    
    get_parameters.frmtmb_fit    sites  29  reads_dots FALSE
        e.g.  insight::.boot_em_df passes summary 
    get_varcov.frmtmb_fit        sites  17  reads_dots FALSE
        e.g.  insight::.get_simulation_from_zi passes component 
    recover_data.frmtmb_fit      sites   8  reads_dots TRUE 
        e.g.  emmeans::recover_data.aovlist passes na.action 
    get_vcov.frmtmb_fit          sites   4  reads_dots FALSE
        e.g.  marginaleffects::comparisons passes vcov, type 
    emm_basis.frmtmb_fit         sites   1  reads_dots FALSE
        e.g.  emmeans::emm_basis.glmgee passes vcov. 
    find_formula.frmtmb_fit      sites   1  reads_dots FALSE
        e.g.  insight::.recover_data_from_environment passes skip_dot_formula 
    find_random.frmtmb_fit       sites   0  reads_dots FALSE
    find_statistic.frmtmb_fit    sites   0  reads_dots FALSE
    get_coef.frmtmb_fit          sites   0  reads_dots FALSE
    get_predict.frmtmb_fit       sites   0  reads_dots FALSE
    link_function.frmtmb_fit     sites   0  reads_dots FALSE
    link_inverse.frmtmb_fit      sites   0  reads_dots FALSE
    set_coef.frmtmb_fit          sites   0  reads_dots FALSE
    
    exemptions with measured need : 6
    exemptions inert (reads dots) : 1
    exemptions with NO evidence   : 7
    ---- END GENERATED: exemption evidence ----

This is the WRITTEN-call-site count. It is necessary and not
sufficient: `get_coef` and `set_coef` score 0 here and still needed
their exemption, which only running marginaleffects showed. See 10.3.


### 8.2 The core suite, ONE TEST FILE PER PROCESS

`Rscript dev/argspell-suite.R .`, full log in
`dev/argspell-core-suite-log.txt`. Run after review round 2:

    ---- GENERATED COUNTS (paste verbatim) ----
    package        : frmtmb
    files run      : 137
    crashed        : 0
    pass 8688  fail 0  error 0  skip 133  in 1023 s
    files below their 0.57.0 baseline count: 2
                     file pass pass_base skip skip_base
     test-api-spellings.R   33        36    0         0
      test-brms-methods.R   16        18   45        45
    baseline files with no run: 0
    ---- END GENERATED COUNTS ----

Both files below baseline are assertion restructuring, not lost
coverage, and the arithmetic is per block:

* `test-api-spellings.R` 33 against 36, -3. "carries both spellings"
  3 to 2; "giving both spellings is refused" 4 to 3; "identical
  pp_check output" 1 to 0. Every assertion removed tested an alias that
  no longer exists.
* `test-brms-methods.R` 16 against 18, -2, split 3 to 2 on the gaussian
  case and 2 to 1 on the ordinal one. Each case loses its
  `expect_length(w, 1L)`, which counted a warning VECTOR that cannot
  exist once the warning is an error.

The three runs of this suite across the rounds: 8614 with one error
(`hypothesis(method = "wald", ytol = 8)`, a guard that refused on every
branch when only one forwards), then 8670 clean, then 8688 clean. The
growth is `test-arg-refusal.R`, which went 47, 100, 118 as each round
added the assertions its fixes needed.


### 8.3 The frmtmb.sample suite

`Rscript dev/argspell-suite.R extensions/frmtmb.sample`, full log in
`dev/argspell-sample-suite-log.txt`. All seven extensions installed
together, so the ordering constraint of section 6 is satisfied:

    ---- GENERATED COUNTS (paste verbatim) ----
    package        : frmtmb.sample
    files run      : 18
    crashed        : 0
    pass 1204  fail 0  error 0  skip 3  in 244 s
    files below their 0.57.0 baseline count: 0
    baseline files with no run: 0
    ---- END GENERATED COUNTS ----

Nothing below baseline. `test-draws-spellings.R` is 52 against its
baseline 48: restoring the `predictive_interval()` alias put its
equivalence block back, and the detector that derives "brms accepts it"
from brms's own code added assertions of its own.


### 8.4 `R CMD check --as-cran`, once, no `--no-manual`

`dev/argspell-check.ps1`; full log copied to `dev/argspell-00check.log`.
Run after review round 2:

    * checking examples ... [38s] OK
    * checking examples with --run-donttest ... [39s] OK
    * checking tests ...
      Running 'testthat.R' [178s]
     [178s] OK
    * checking re-building of vignette outputs ... [202s] OK
    * checking PDF version of manual ... OK
    * checking HTML version of manual ... NOTE
    Skipping checking math rendering: package 'V8' unavailable
    * DONE
    Status: 1 NOTE

The one NOTE is the expected V8 math-rendering one. The in-check test
run, which is the whole suite in ONE process rather than one file per
process, reported:

    [ FAIL 0 | WARN 1 | SKIP 265 | PASS 6088 ]

identical across all three rounds, so nothing in the review work moved
it. The counts differ from 8688 because `R CMD check` leaves `NOT_CRAN`
unset, so every `skip_on_cran()` file skips, `test-arg-refusal.R` among
them. Its evidence is the per-file run in 5.1.

The one WARN is identified in 10.8 and is not this lane's.


## 9. The version bump, stated but not chosen

`DESCRIPTION` is untouched in both packages: frmtmb stays 0.57.0 and
frmtmb.sample stays 0.5.0 in this worktree, because choosing the number
is the consolidating session's call.

What the number has to carry: this is a BREAKING change to the exported
signature of both packages, not a patch. Calls that worked at 0.57.0
now error rather than returning a different answer, which is the point.
Concretely, a call breaks if it uses `re.form` or `allow.new.levels`
anywhere on the fit surface, `re.form` on `pp_check()` or
`predictive_interval()`, or any argument a method does not have. NEWS
in both packages is written under `(development version)` and needs
only the heading replaced.

## 10. Review round 1: what was falsified, and what fixing it turned up

`dev/reviews/20260916-argspell.md`. Two BLOCKERs, five named fixes,
four nits. The interop result and the `fitted()` identity both held.

### 10.1 BLOCKER: the swallow detector was fooled by a string

`any(grepl("...", deparse(body(fn)), fixed = TRUE))` matches three dots
inside a STRING LITERAL. `print.frmtmb_par_template()` prints
`cat("... ", length(v) - n, " more\n")` and never reads its dots, so it
was scored as a dots user and went on swallowing while the generated
line said `METHODS STILL SWALLOWING THEIR DOTS: 0`. The count that
certified the whole class fix was false, in the direction that matters.

Both detectors now read the PARSE TREE:

    if (any(all.names(body(fn)) %in%
            c("...", "..1", "..2", "..3",
              "...length", "...names", "...elt"))) return(FALSE)

and `print.frmtmb_par_template()` is guarded. My synthetic cases had
never included a literal containing dots, which is exactly why it
survived; the detector test now covers `message("a ... b")` and
`cat("... ", 1, " more\n")` as swallowers, and `unclass(x)[[i, ...]]`,
`match.call(expand.dots = FALSE)$...` and `...length()` as users, so
the rewrite is exercised in both directions.

### 10.2 BLOCKER: over-refusal on `stats::step()`, and a second hop

`stats::step()`, `add1.default()`, `drop1.default()` and
`sigma.default()` call `nobs(object, use.fallback = TRUE)`.
`use.fallback` is `stats::nobs.default`'s own formal, and
`nobs.frmtmb_fit` refused it, so `step(fit)` errored where it had
worked. Fixed with `.allow = "use.fallback"`.

Fixing it exposed a SECOND hop the reviewer's run could not see,
because `step()` died at the first: `stats::step()` then calls
`drop1(object, scope, scale = , trace = , k = )`, and
`drop1.frmtmb_fit` refused `scale` and `trace`. Found by the new
behavioural test, not by reading. `.allow = c("scale", "trace")` there,
and `test-arg-refusal.R` now runs `step(fit, trace = 0)` end to end and
separately checks that a name stats does NOT pass is still refused.

### 10.3 The exemption list, and a measurement that was not enough

`dev/argspell-exempt2.R` counts, per generic, the WRITTEN call sites in
emmeans, insight and marginaleffects that pass a name our method has no
formal for, which is what a guard would actually refuse. That is a
stronger test than the first pass's "forwards dots anywhere", which
counted sites that can never reach a `frmtmb_fit`.

    get_parameters.frmtmb_fit    sites  29
    get_varcov.frmtmb_fit        sites  17
    recover_data.frmtmb_fit      sites   8   reads_dots TRUE
    get_vcov.frmtmb_fit          sites   4
    emm_basis.frmtmb_fit         sites   1
    find_formula.frmtmb_fit      sites   1
    find_random, find_statistic, get_coef, get_predict,
    link_function, link_inverse, set_coef    sites 0

`recover_data` left the list because it reads its dots, so exempting it
was inert. I then guarded the five the reviewer named as unevidenced.

**That was wrong for two of them, and only RUNNING the third-party
packages showed it.** Rerunning the reviewer's own interop harness
against the lane library (`dev/argspell-interop.R`,
`dev/argspell-interop-log.txt`) broke four marginaleffects entry
points:

    FAIL  avg_slopes(fit)      get_coef() has no argument `variables`
                               (and 2 more: numderiv, internal_call)
    FAIL  slopes(fit)          same
    FAIL  avg_comparisons(fit) same
    FAIL  hypotheses(...)      same

marginaleffects assembles that call rather than writing it out, so
BOTH static passes, the reviewer's and mine, scored `get_coef` at zero
and both were wrong. `get_coef` and `set_coef` are exempt again, on
runtime evidence, and the comment at the exemption says so. The three
that survive the harness (`find_statistic`, `link_function`,
`link_inverse`) are guarded.

The list is thirteen -> nine. After the restore the harness is back to
the reviewer's numbers: emmeans 6 of 6, insight 11 of 11, and
marginaleffects 6 of 6, with the only failures the three insight calls
that fail identically at 0.57.0 and three packages that are not
installed.

**The rule this earned:** a count of zero WRITTEN call sites is not
evidence that no call site exists. An exemption may be removed only
against a run.

### 10.4 `predictive_interval()`: declares is not accepts

The reviewer falsified the formals-deep reading.
`predictive_interval.brmsfit` declares neither spelling and its whole
body is `posterior_predict(object, ...)`, and
`posterior_predict.brmsfit` HAS a `re.form` formal, so brms accepts and
honors `predictive_interval(x, re.form = NA)`. The alias is restored;
the dual set is five, not four. `dev/argspell-brms-accepts.R` measures
declared-versus-accepted for all six, and the test now derives
"brms accepts" structurally instead of comparing declared formals.

**And the same test points the other way for `pp_check()`, which I have
NOT acted on.** The instruction was that the `pp_check()` half stands
because `prepare_predictions.brmsfit` carries `re_formula` alone. That
premise does not hold for the yrep path. Measured:

    method <- "posterior_epred"
    method <- "posterior_predict"
    pred_args <- list(object, newdata = , resp = , draw_ids = , ...)
    yrep <- do_call(method, pred_args)

`pp_check.brmsfit`'s dots reach `posterior_predict`, not
`prepare_predictions`, so by the very argument that restored
`predictive_interval`'s alias brms accepts `pp_check(x, re.form = NA)`
too. The generated surface block now prints both columns and shows the
divergence rather than hiding it:

    pp_check.frmtmb_draws   ours FALSE  brms declares FALSE  accepts TRUE

I left the behaviour alone, because the instruction was explicit and
narrow and the re-check is scoped to the two blockers and the
`predictive_interval` fix. The evidence is static only: I cannot
execute brms's `pp_check()` without a fitted `brmsfit`. **Filed for the
reviewer to settle**, and the test comment says the assertion states
what this package does rather than claiming brms settles it.

### 10.5 brms arguments refused with a reason rather than as unknown

The principle `fitted_no_draws` states was applied to three methods and
nothing else. Now: `print()` and `print(summary())` name `digits` and
`short`; `summary()` names `priors`, `prob`, `mc_se`; `coef()`,
`fixef()`, `ranef()`, `VarCorr()` and `summary()` share
`brms_draws_summary_args` (`summary`, `robust`, `probs`, `pars`);
`nobs()` and `family()` name `resp`; `vcov()` names `correlation` and
`pars`. Each message says what to use instead.

`digits` is refused rather than implemented, which is the narrower of
the two options the reviewer offered. brms's `print.brmsfit` has it and
brms is the tiebreaker, so **implementing `digits` on
`print.frmtmb_fit()` and `print.summary.frmtmb_fit()` is left as an
open item**, not a settled decision.

`calling_fun_name()` reported `print.summary()` for
`print.summary.frmtmb_fit`, which is not a function. The pattern takes
an optional `summary.` now and the test asserts the message says
`print()` and does not say `print.summary(`.

### 10.6 The migration vignette says what does not transfer

"Method conventions" advertised the name transfer and was silent about
the scale. It now says that `predict()` keeps the LINK default where
the brms call it replaces returns the response scale, gives the
measured 0.114223 against 1.121002 on a poisson fit, names
`type = "response"` and `fitted()` as the two ways to get brms's scale,
and records that the 0.57.0 warning on that exact call is gone. The
default itself is still left to the consolidating session, for the
reasons in section 4.

### 10.7 Nits closed

* `man/fitted.frmtmb_fit.Rd` said "lme4's `re_formula` is not accepted
  here" where it meant `re.form`, so the page denied the argument it
  documents. Fixed at the roxygen source.
* `frm_check_dots()` degrades to the HELPER's name and an empty formal
  list one frame down. A new test counts every mention of the helper in
  both packages' sources two ways, `all.names()` anywhere against a
  structured pass over top-level statements of dots-taking top-level
  functions, and fails if they differ. The helper's own definition is
  subtracted and the subtraction is asserted to be exactly one.
* `fit_no_draws()` and `multiple_no_draws()` carry 38 shape exemptions
  between them and nothing asserted they were unconditional. A test now
  asserts each body is a braced block of exactly one `stop()`, with the
  absent case (a body with a branch) asserted not to satisfy it.
* `expect_error(pp_check(ds, re.form = NA))` has its `"re.form"`
  pattern back.

### 10.8 Filed, not fixed

The check's one WARN is `test-prior-compat.R`, "Optimizer did not
report convergence: singular convergence (7)", from a two-correlated-
blocks fixture on a small design. It is pre-existing and not this
lane's. Recording the reviewer's point rather than the count: a
convergence warning living inside a PRIOR test means a REAL convergence
regression would be absorbed into a count that is already one. Either
the fixture should be strengthened or the warning expected explicitly.

## 11. Review round 2: the rule was wrong, not the patches

`dev/reviews/20260916-argspell.md`, appendix. One BLOCKER needing a
DESIGN change, one record fix, one nit.

### 11.1 `pp_check()`: the mechanism, the reason, and the refusal

The reviewer settled it by RUNNING brms, building a real `brmsfit`
from a cached program with
`rstan::sampling(algorithm = "Fixed_param")`:

    pp_check(ndraws = 5, re.form = NA)  ok, warns
                                        "unrecognized and ignored: re.form"
      re.form = NA changes the answer vs default : TRUE
      re.form = NA equals re_formula = NA        : TRUE
      a TYPO (re_frmula) changes the answer      : FALSE

So brms HONORS `re.form` on `pp_check()` while warning that it ignored
it. My mechanism was right: the value rides the dots into
`do_call("posterior_predict", pred_args)`. The
`prepare_predictions` premise, which I inherited and repeated, was
wrong, and filing it rather than acting on it was the right call.

**The refusal stands, and the REASON changed in three places.** brms
warn-then-honor is a leak, not a decision, and warn-then-honor is the
exact failure this item exists to kill: the argument changes the answer
and the message says it did not. Matching brms means matching what brms
decided. `predictive_interval()` is the contrast, and it is why the
distinction is real rather than convenient: there brms honors the alias
SILENTLY, with a maximum difference of 0 against `re_formula`, and this
package follows it.

Rewritten: the refusal message in `R/conditional-effects.R` and
`extensions/frmtmb.sample/R/methods-draws.R`, the `@param` text, the
*Argument spellings* section, the frmtmb.sample NEWS entry, and the
header and solo-case comment of `test-draws-spellings.R`.

### 11.2 BLOCKER: one table, because the RULE was wrong

Three rounds each uncovered the next hop, and the reviewer's scan found
a third and a fourth: `as.formula(fit)` (which is literally
`formula(object, env = baseenv())`), `terms(data =)`,
`model.frame(data =)`, `confint(trace =)`, `residuals(na.rm =)`,
`coef(complete =)`. The collision is structural, and the reviewer
states it exactly: the refusal's rule was "refuse what the METHOD does
not have", and R's S3 contract is "a method tolerates the arguments its
GENERIC carries". That is why every `stats` default method declares
`use.fallback`, `env`, `data`, `na.rm` or `trace` and ignores them.

`s3_contract_args` in `R/utils.R` is that contract, consulted by
`frm_check_dots()` for every method at once. The two ad-hoc `.allow`
entries it subsumes, `nobs`'s `use.fallback` and `drop1`'s
`scale`/`trace`, are gone, so there is one mechanism and not two. The
`.allow` argument survives only for its real purpose: a method that
FORWARDS its dots somewhere specific, which is `residuals(type = "osa")`
against `formals(TMB::oneStepPredict)`, `confint()` and `hypothesis()`
against whatever the chosen method forwards to, and `pp_check()`'s
catch-all.

The table is derived, not typed. `dev/argspell-contract.R` walks the
ASTs of R's own packages for real calls to every generic either package
registers a dots-refusing method on, and reports the named arguments
our method has no formal for. The block below is what the SCAN found;
the shipped `s3_contract_args` is slightly wider, because it also
carries the names the coordinator's own table named that R passes at
sites this scan filters (`model.frame`'s `subset` and
`drop.unused.levels`, `coef`'s `complete`, `residuals`'s `na.rm`,
`confint`'s `trace` and `test`, `drop1`/`add1`'s `k` and `test`) and
because `plot`'s row is the generated graphical set rather than only
the 37 names seen here:

    ---- GENERATED: S3 generic contract scan ----
    packages scanned: stats, base, utils, graphics, grDevices, methods 
    generics with a dots-refusing method: 65
    call sites, table EMPTY             : 62
    call sites, table APPLIED           : 0
    
    what the table carries, generic to the names R itself passes:
      as.data.frame  fix.empty.names, optional, row.names, stringsAsFactors
      as.matrix      rownames.force
      confint        type
      drop1          scale, trace
      formula        env
      influence      do.coef
      model.frame    data, na.action, xlev
      model.matrix   contrasts.arg, data
      nobs           use.fallback
      plot           add, angle, asp, axes, bg, border, cex, ci, col,
                     col.lab, density, font.lab, freq, horiz, labels,
                     leaflab, lend, log, lty, lwd, main, mgp, panel.first,
                     pch, sub, type, verticals, xaxs, xaxt, xlab, xlim,
                     xy.labels, xy.lines, y, yaxs, ylab, ylim
      predict        terms
      print          quote, steps, useS4
      summary        dispersion
      terms          allowDotAsName, data, specials
      weights        type
    ---- END GENERATED: S3 generic contract scan ----

62 to 0. `tests/testthat/test-arg-refusal.R` runs the same scan through
the same file, `dev/argspell-contract-scan.R`, so the table and the test
cannot drift; the NEXT hop fails the suite instead of reaching a user.
The test is non-vacuous in both directions: it asserts the scan reached
more than 40 generics, and that with the table EMPTIED it still reports
`formula`/`env` and `nobs`/`use.fallback`, so a scan that silently found
nothing would fail.

Two scoping decisions, both of which cost a run to learn:

* **R's own packages only.** The first scan included emmeans, insight,
  marginaleffects, lme4, nlme and mgcv and returned 224 sites, almost
  all of them a package calling a generic on its OWN objects:
  `mgcv::plot.gam` calling `plot(..., too.far =)` is not a claim about
  what our method must tolerate. The S3 contract is R's.
* **A method of the same generic is skipped.** `plot()` inside
  `plot.dendrogram()` dispatches on a dendrogram.

`plot`'s row is 37 names and they are all graphical parameters, which is
correct rather than a capitulation: a plot method that refuses `col` is
broken whatever it does with it. `s3_plot_args` is generated as the
union of `par()`, `plot.default()`'s formals and the names R's other
plot methods pass. The message does NOT list them back, because
printing 95 names would bury the answer; they are used for the
"did you mean" search only.

Verified by running the reviewer's own scripts against the lane library
(`dev/argspell-hop3-log.txt`, `dev/argspell-hop4-log.txt`):

    as.formula(fit)                  ok formula
    formula(fit, env = baseenv())    ok formula
    update.formula(fit, ~ . + w)     ok formula
    stats::add.scope(fit, ~ x+z+w)   ok character
    terms(fit, data = dd)            ok terms/formula
    model.frame(fit, data = dd)      ok data.frame
    confint(fit, trace = FALSE)      ok matrix/array
    residuals(fit, na.rm = TRUE)     ok numeric
    step(one), step(one, trace = 0)  ok character
    nobs(one, use.fallback = TRUE)   ok integer

and the guard keeps its value: `nobs(fit, nosucharg = 1)`,
`formula(fit, envv =)`, `terms(fit, re.form = NA)` and
`coef(fit, complet =)` still error, and an `.unsupported` reason still
wins, so `print(fit, digits = 3)` keeps its message.

Two calls the reviewer listed stay REFUSED, deliberately:
`df.residual(fit, na.rm = TRUE)` and `coef(fit, digits = 3)`. Neither
name is passed to that generic anywhere in R's own packages, so they are
outside the contract, and tolerating them would be tolerating a typo.
That is the rule the coordinator set: everything not in the table stays
refused.

Out of scope and left alone, as instructed: `step(fit)` on a model with
two droppable fixed terms fails with "object 'change' not found" at
BOTH versions. What the change did fix is the MASKING: `add1()` and
`influence.measures()` now fail at their own pre-existing cause rather
than at a refusal.

### 11.3 Nit: `..N` by pattern

`ar_dots_names` enumerated `..1`, `..2`, `..3`, so `..4` and `..11`
were legitimate dots reads scored as swallowers. Nothing in either
package uses one today, which is exactly why a run would not have found
it. Matched by `grepl("^[.][.][0-9]+$", x)` now, in the test and in
`dev/argspell-report.R`.
