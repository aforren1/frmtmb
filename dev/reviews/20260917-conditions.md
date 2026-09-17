# Review: lane wt-conditions, item 2.6e (classed conditions)

Reviewer, round 1 (restarted after a shutdown). Worktree
`frmtmb-wt-conditions` at base e031e8c. Lane build `conditions-lib`,
base build `rellib-r3`, `pinlib` ahead of the user library.

## Verdict

Mergeable after one punch round. No BLOCKER and no MAJOR. The core
claim holds: on 2,024 errors, 72 warnings and 39 messages that the
suites catch, every converted condition has the base message, the base
call (errors and warnings) and the base `try()` text, and the printed
output of 61 real sites is byte-identical to the base build in every
mode tested. Four MINOR findings and three NITs follow. The MINOR items
are one interop regression, a documented but large gap in the extension
subclass, census holes that no current file uses, and helper behavior
that no current site reaches.

## The instrument

The lane build matches the worktree: `dev/conditions-rev-verifylib.R`
compares 1,591 top-level functions in the 8 namespaces with the parsed
sources, 0 differ (`dev/conditions-rev-log/verifylib.txt`, re-run this
round). Its control against the base build finds 308 of 1,024 core
functions differ, so it can see a difference. The installed `inst/`
trees are identical to the worktree (`diff -rq`).

All scripts below are in `dev/conditions-rev-*`; logs are in
`dev/conditions-rev-log/`. The earlier reviewer's scripts were re-run,
not trusted.

## Findings

### MINOR 1. The extension subclass is absent on 17% of extension refusals

Construction: `dev/conditions-rev-subclass-sweep.R` on the lane sweep
(`sweep2-lane`, every test file of every suite, one per process).

In each extension's own suite, the caught `frmtmb_error` conditions
that lack that extension's subclass:

| package | caught | with own subclass | without |
|---|---|---|---|
| frmtmb.coupling | 32 | 28 | 4 |
| frmtmb.eam | 192 | 162 | 30 |
| frmtmb.latent | 55 | 42 | 13 |
| frmtmb.learn | 56 | 34 | 22 |
| frmtmb.ode | 78 | 76 | 2 |
| frmtmb.sample | 138 | 108 | 30 |
| frmtmb.spline | 55 | 54 | 1 |
| all | 606 | 504 | 102 |

This follows the stated rule (the raising function decides), but the
rule defeats the documented use. Examples:

- frmtmb.latent writes its refusal text as data
  (`refusals = list(reml = "REML = TRUE cannot be combined with hmm() ...")`,
  `extensions/frmtmb.latent/R/hmm.R` line 621) and core raises it:
  `frm(bf(y ~ 1), family = hmm(...), data = d, REML = TRUE)` is a plain
  `frmtmb_error`. `tryCatch(frmtmb_latent_error = )` does not catch it.
- `frm(bf(rt ~ 1, bias = 0.5), family = wiener(), data = dat)`: "wiener:
  the density needs one of `dec` or `vint1`" is a plain `frmtmb_error`.
- frmtmb.learn's `rlddm()` refusal "rlddm(): this model bounds the
  non-decision time ..." is raised by `frmtmb.eam` code and carries
  `frmtmb_eam_error`, not `frmtmb_learn_error`.

Base: no subclass exists. Lane: `?frmtmb-conditions` says
`tryCatch(frmtmb_eam_error = )` "catches only the refusals of
frmtmb.eam"; it catches 84% of them. The page does add that the package
"is not always the package the user called", so this is MINOR, not
MAJOR. A fix: when core raises a refusal a family declared, pass the
family's owning package as `class =`; or document the subclass as
best-effort and name these cases.

### MINOR 2. insight::get_predicted_ci() on a multivariate fit now errors

Construction: `dev/conditions-rev-interop2.R base|lane`, logs
`interop2-base.txt` and `interop2-lane.txt`.

insight 1.5.4 `.get_predicted_ci_modelmatrix()` catches the error of
`get_modelmatrix()` and tests `inherits(mm, "simpleError")`.
`model.matrix()` on an `mvbf()` fit refuses ("Multiple responses have
dpar 'mu'").

    mv <- frm(mvbf(bf(y ~ x), bf(z ~ x)), data = d)
    insight::get_predicted_ci(mv, predictions = rep(0, 60), data = d)

Base: a warning "Something went wrong with computing standard errors
..." and a data frame of NA. Lane: `Error: argument is of length zero`.
Reachable only where a frmtmb `model.matrix()` refuses; a nonlinear fit
returns `NULL` there and is unaffected (same script). emmeans,
marginaleffects and the rest of insight gave identical output on a
gaussian, a binomial and a mixed fit, including their error paths
(`dev/conditions-rev-interop.R`, `interop-base.txt` against
`interop-lane.txt`: the only differences are 4 environment addresses).
The fix is to add `"simpleError"` to the error class vector, or accept
this and say so; the first costs brms parity (brms has no
`simpleError` either).

### MINOR 3. The census has holes that no current file uses

Construction: `dev/conditions-rev-census.R` (re-run,
`census-rerun.txt`) and `dev/conditions-rev-census2.R`
(`census2.txt`). Each case plants one site in a scratch copy and runs
the package's own census test. The census DOES fail on all 8 base
trees (`dev/conditions-rev-census-base.R`, `census-base.txt`), on a bare
call in `R/`, on `base::`, on `Map(warning, ...)`, on
`do.call(what = "stop", list())` and on a bare call in an extension.

It PASSES with each of these plants:

- a second `stop(e)` in the same file and function as an exempt one
  (`setdiff()` merges identical site keys), in `R/fit.R` and in
  `R/conditions.R`;
- `R/windows/*.R` and `R/unix/*.R`, which `R CMD INSTALL` collates, and
  `R/*.S`, `R/*.q`;
- `get("warning")(...)`, `match.fun("stop")(...)`,
  `eval(call("stop", ...))`, `eval(parse(text = "stop()"))`,
  `eval(str2lang(...))`, `.signalSimpleWarning()`, `.Deprecated()`;
- `do.call(args = list("p"), what = "stop")`, so the claim "through
  `do.call`" holds only for a positional or leading `what`;
- `lapply(x, frm_stop)`, which is not bare but gets the wrong subclass
  (MINOR 4).

It SKIPS, and so passes silently, when the sentinel file it uses to
recognize a source tree is renamed (`R/objective.R` in core,
`R/<pkg>-package.R` in each extension). This is the house pattern of
`test-bracket-access.R` and `test-message-uniqueness.R`, and under
`R CMD check` the census skips for the same reason.

Current trees: 0 `R/windows` or `R/unix` directories, 0 non-`.R` files
in any `R/`, and no `get`, `match.fun`, `call("stop")` or
`.signalSimpleWarning` use (grep). The duplicate-key hole is the one
worth closing: count sites per key, not `setdiff()`.

### MINOR 4. Helper behavior that differs from base, reached by no current site

Construction: `dev/conditions-rev-subclass.R` (re-run,
`subclass-rerun.txt`) and `dev/conditions-rev-shapes.R` (re-run,
`shapes-rerun.txt`, 54 shapes, each a fresh `Rscript` per arm).

`?frmtmb-conditions` tells extension authors the helpers "take the same
arguments as the base functions and give the same message text and the
same call". That is false in these cases:

1. A condition object. `stop(e)` gives message `boom`;
   `frm_stop(e)` gives `Error in f(1): boom`. `warning(w)` gives
   `careful`; `frm_warning(w)` gives `simpleWarning: careful`. No site
   passes a condition (grep of every single-argument call: 0).
2. A helper passed as a value. `lapply("x", frm_stop)` from an
   extension function gives `base_error/frmtmb_error/...`;
   `Map(frm_warning, "x")` gives `base_warning/...`. The subclass names
   base R. No site does this (grep: 0).
3. `immediate. = TRUE` is ignored: base prints `Warning in f() : imm`
   before the next line of output; the lane defers it to the end. No
   site passes `immediate.` (grep: 0).
4. With the default `call. = TRUE`, an uncaught error under `Rscript`
   prints an extra line (`Calls: f -> frm_stop`), which is R's concise
   traceback with the helper frame. Of 1,325 converted sites, 1 error site uses
   the default (`cr_adjust()`, an unreachable `switch()` fallback in
   `R/sandwich.R`); all 48 warning sites and 1,252 error sites pass
   `call. = FALSE`, where the output is identical
   (`dev/conditions-rev-sitetable.R`).

Fix 1 and 2 in the helpers (pass a condition through; find the caller
with `sys.function(-1)` and refuse a non-closure), or state the limits
on the help page. 3 and 4 need a sentence each.

### NIT 1. frm_message() records its own call

`conditionCall()` of a message is `frm_message(...)` where base gives
`message(...)`: 39 of 39 caught messages in the sweep, 24 sites. Base
never prints a message's call, and the printed output is identical.

### NIT 2. testthat labels an uncaught frmtmb error

testthat 3.3.2 `as.expectation.error()` prefixes
`<frmtmb_error/error/condition>` to an error that is not
`simpleError`, so a test that errors reports one more line. Developer
output only.

### NIT 3. The sweep count is a different construction

This review's tracer (`dev/conditions-rev-sweep-onefile.R`, one row per
condition returned from each `expect_error()` call) counts 2,053 errors
caught, 2,024 `frmtmb_error`, 29 other. The findings record 2,032,
2,003 and 29. The 29 are the same set. The difference of 21 is in how
repeated `expect_error()` calls are counted; not a record defect unless
the worker's script says otherwise.

## Decided by the user, listed but not ranked

- `getME(mv_fit, "Zt")`: besides the class, the lane changes the
  printed call. `try()` prints `Error in h(simpleError(msg, call)) :`
  on base and `Error in (function (cond)  :` on the lane; the message is
  identical. It is the only printed difference in all 61 real sites
  and 4 modes. R's S4 dispatch calls its handler with a C-built call
  for a C-level error and not for `stop(<condition>)`. The one-line fix
  removes both. A grep for other S4 generics applied to a raising call
  found only this one (heuristic).
- `stopifnot()` and `match.arg()`: the full user-reachable list is
  below.

Also carried, not a finding: each extension imports `frm_stop` while
its DESCRIPTION still says `frmtmb (>= 0.59.0)`, which does not export
it. The floor must move at release.

## Every user-reachable unclassed refusal

Probed on the lane build: `dev/conditions-rev-unclassed-probe.R` (re-run,
26 probes, identical to the earlier log) and
`dev/conditions-rev-unclassed-probe2.R` (22 probes). Every probe below
raised a `simpleError` that `tryCatch(frmtmb_error = )` does not catch.

### match.arg() on a user argument (37 sites)

Core:

| call | site |
|---|---|
| `predict(fit, type = )` | `R/predict.R` 1217 |
| `fitted(fit, scale = )` | `R/predict.R` 1918 |
| `residuals(fit, type = )` | `R/predict.R` 2548 |
| `confint(fit, method = )` | `R/confint.R` 377 |
| `drop1(fit, test = )` | `R/confint.R` 1730 |
| `hypothesis(fit, method = )` | `R/confint.R` 2720 |
| `hypothesis(fit, scope = )` | `R/confint.R` 2693 |
| `conditional_effects(fit, band = )` | `R/conditional-effects.R` 1701 |
| `conditional_effects(fit, method = )` | `R/conditional-effects.R` 458 |
| `frmtmb_control(check_nlev_1 = )` | `R/fit.R` 1533 |
| `frmtmb_control(check_olre = )` | `R/fit.R` 1534 |
| `default_prior()` and `get_prior(route = )` | `R/priors.R` 1438 |
| `validate_prior(route = )` | `R/priors.R` 1517 |
| `vcov_cluster(type = )` | `R/sandwich.R` 463 |
| `frm_periodogram(taper = )` | `R/spectral.R` 161 |
| `frm_periodogram(detrend = )` | `R/spectral.R` 162 |
| `anova(<frmtmb_multiple>, method = )` | `R/multiple.R` 243 |
| `anova(<frmtmb_multiple>, use = )` | `R/multiple.R` 244 |

Extensions:

| call | site |
|---|---|
| `frm_cross_spectrum(window = )` | frmtmb.coupling `cross-spectrum.R` 229 |
| `gddm_control(tridiagonal = )` | frmtmb.eam `gddm.R` 761 |
| `gddm(lapse = )` | frmtmb.eam `gddm.R` 988 |
| `gddm_simulate(lapse = )` | frmtmb.eam `gddm.R` 1813 |
| `hmm(init = )` | frmtmb.latent `hmm.R` 408 |
| `bandit2arm_dual(split = )` | frmtmb.learn `bandit2arm-dual.R` 75 |
| `frm_task_design(task = )` | frmtmb.learn `task.R` 42 |
| `frm_lincmt(output = )` | frmtmb.ode `lincmt.R` 758 |
| `frm_ode(on_error = )` | frmtmb.ode `ode.R` 1652 |
| `loo_compare(<draws>, criterion = )` | frmtmb.sample `loo.R` 501 |
| `hypothesis(<draws>, scope = )` | frmtmb.sample `methods-draws.R` 283 |
| `pp_check(<draws>, prefix = )` | frmtmb.sample `methods-draws.R` 777 |
| `predictive_error(<draws>, method = )` | frmtmb.sample `methods-draws.R` 1315 |
| `frm_curve_feature(type = )` | frmtmb.spline `curve-feature.R` 125 |
| `royston_parmar(scale = )` | frmtmb.spline `royston-parmar.R` 290 |
| `rp_floored(action = )` | frmtmb.spline `rp-check.R` 164 |

Not user-reachable: `autoscale_map(to = )` (`R/autoscale.R` 69) and
`ln_recurse(mode = )` (frmtmb.learn `engine.R` 176), both internal with
constant arguments.

### stopifnot() on user input (19 sites)

| call | site | text |
|---|---|---|
| `set_prior(1)` | `R/priors.R` 961 | `is.character(prior) is not TRUE` |
| `set_prior("normal(0)")` | 982 | `length(pars) == 2 is not TRUE` |
| `set_prior("student_t(3, 0)")` | 986 | `length(pars) == 3 is not TRUE` |
| `set_prior("cauchy(0)")` | 990 | same shape |
| `set_prior("exponential(-1)")` | 994 | `pars[1] > 0 is not TRUE` |
| `set_prior("lkj()")` | 999 | same shape |
| `set_prior("logistic(0)")` | 1003 | same shape |
| `set_prior("gamma(1)")` | 1007 | `length(pars) == 2 is not TRUE` |
| `set_prior("inv_gamma(1)")` | 1011 | same shape |
| `set_prior("beta(1)")` | 1015 | same shape |
| `set_prior(...) + "x"` | 1041 | `inherits(e2, "frmtmb_priorlist") is not TRUE` |
| `frm(..., prior = list(1))` | 3220 | `!is.null(names(prior)) is not TRUE` |
| `diagnose(lm(...))` | `R/confint.R` 1187 | `inherits(fit, "frmtmb_fit") is not TRUE` |
| `frmtmb_family(family = 1, ...)` | `R/families.R` 463 | `is.character(family) is not TRUE` |
| `check_custom_family(1)` | `R/interop.R` 64 | `inherits(family, "frmtmb_family") ...` |
| `cluster_scores(lm(...), g)` | `R/sandwich.R` 296 | `inherits(object, "frmtmb_fit") ...` |
| `vcov_cluster(lm(...), g)` | `R/sandwich.R` 451 | same |
| `frmtmb_register_prior_defaults(1)` | `R/sampling-api.R` 410 | `is.function(provider) is not TRUE` |

The cauchy, lkj, logistic, inv_gamma and beta rows are the same
`switch()` in `parse_prior_dist()` as the probed rows; they were not
probed one by one.

Not user-reachable (internal invariants): `R/autoscale.R` 175,
`R/brms-names.R` 45, `R/parse.R` 2065, `R/predict.R` 2239, frmtmb.sample
`sample.R` 2153 (after sampling in `check_laplace()`).

### Other unclassed errors a user meets through a frmtmb function

From the lane sweep's 29 (`sweep2-others.txt`) and the probes:

| call | error | raised by |
|---|---|---|
| `lme4::getME(mv_fit, "Zt")` | S4 rewrap | decided fix |
| `frmtmb.latent::hmm_starts(1)` | `$ operator is invalid for atomic vectors` | frmtmb reads before checking |
| `predict(fit, newdata = <missing column>)` | `object 'z' not found` | stats; brms raises `brms_error` |
| `predict(fit, newdata = <new factor level>)` | `factor f has new levels zz` | stats; brms raises `brms_error` |
| `hypothesis(<draws>, "`(Intercept)` > 0")` | `object '(Intercept)' not found` | evaluation of the hypothesis |
| `frm(bf(y ~ b0 * no_such, b0 ~ 1, nl = TRUE), ...)` | `object 'no_such' not found` | model.frame |
| `frm(bf(y ~ b0 * short, ...), ...)` | `variable lengths differ` | model.frame |
| `frm(... ~ tab ..., data = <list column>)` | `invalid type (list) for variable` | model.frame |
| `frm(y ~ nope, data = d)` | `object 'nope' not found` | model.frame |
| `frm(y ~ x, data = "notdf")` | `'data' must be a data.frame ...` | model.frame |
| `frm(..., priors = p)`, `lower =`, `upper =` | `unused argument` | R argument matching |
| `frm(data = d)` | `argument "formula" is missing` | R |
| `frm_task_simulate(1, d)` | `subscript out of bounds` | frmtmb.learn reads before checking |
| `posterior_interval()`, `rhat()`, `neff_ratio()` on draws with an unknown variable | `The following variables are missing ...` | posterior |
| `log_lik(<frmtmb_fit>)` | `no applicable method` | UseMethod |

## What held

- **Class vectors.** Every classed condition seen has exactly
  `c(<ext>_kind?, "frmtmb_kind", kind, "condition")`, plus `class =`
  in front: 9 distinct vectors in the sweep.
- **Message, call and try() text.** Sweep join
  (`dev/conditions-rev-sweep-join.R`, `sweep2-join.txt`): on 2,024
  errors, 72 warnings and 39 messages, the message is identical base
  against lane; the call and the `try()` text are identical for every
  error and warning; 1,772 distinct messages. The one unclassed row
  that differs is getME.
- **Printed output on 61 real converted sites** (`dev/conditions-rev-printed.R`;
  44 errors, 10 warnings, 7 messages, across all 8 packages, including
  the `frmtmb_fit_error` rethrow from `tryCatch()`, warnings from
  `lapply()`, two warnings in one call, a message inside `tryCatch()`,
  and a warning under `options(warn = 2)` inside `local()`/`on.exit()`).
  Four modes, each a fresh `Rscript` per arm, stdout plus stderr and
  stderr alone: `try()` with deferred warnings, `options(warn = 1)`,
  `suppressWarnings(suppressMessages())`, and user `tryCatch(error =,
  warning =, message =)` plus `withCallingHandlers()` with muffle
  restarts. Line counts equal in all 8 pairs; the only differences are
  the getME call line and, by construction, the class name the catch
  mode prints. Messages go to stderr in both arms and vanish under
  `suppressMessages()` in both.
- **Shapes** (`shapes-rerun.txt`): the 54 synthetic shapes reproduce
  the earlier log; every difference is MINOR 4 item 3 or 4, NIT 1, or
  the call text of a raise written inside a promise argument, which
  contains the helper's own name by construction.
- **Subclass at run time.** Correct when the extension is loaded but not
  attached, through `lapply()` closures, `tryCatch()` handlers,
  `do.call()` by function and by name, `Reduce()`, `local()`, and
  `eval()` in a child environment.
- **Cost.** `frm_condition_class()` counted with a trace
  (`dev/conditions-rev-hotcount.R`, `hotcount.txt`): 0 raises in a
  gaussian, binomial and mixed fit, `predict(se.fit = TRUE)`, Wald
  `confint()`, `simulate(nsim = 20)`, `hypothesis()`,
  `frm_bootstrap(nsim = 20)` and a refit `anova()`; 1 in
  `conditional_effects()`. The lookup costs 5.9 us in core and 9.5 us
  in an extension per raise (160,000 calls per block, minimum of 7
  rounds, control arm 1.000). No raise is on a taped path: the one
  R callback during optimization, frmtmb.eam's `ADjoint` tridiagonal
  solver, raises nothing.
- **Rewrite proof.** `dev/conditions-rev-proof.sh` (re-run,
  `proof2.txt`) runs `verify` against the main checkout, whose `R/`
  and `inst/` are the base commit: 113 files, 0 mismatched. A one-letter
  change to a core message, a dropped `call. = FALSE`, and a changed
  extension message each give exactly 1 mismatch. 113 is every `.R`
  file in the 8 `R/` trees and `inst/` except `R/conditions.R`.
- **Counts.** From the parse tree (`dev/conditions-rev-sitetable.R`):
  814 core, 42 `inst/`, 469 in extensions (27, 112, 47, 22, 95, 125,
  41), 3 exempt rethrows; 22 core `stopifnot()` and 1 in frmtmb.sample.
- **Census on base.** Fails on all 8 base trees.
- **Suites.** Every test file of all 8 packages, one per process, on the
  lane build: 247 files, 0 with a failure or an error, 0 without a
  result line (`sweep2-lane/summary.txt`).
- **Interop.** No code in emmeans 2.0.4, marginaleffects 1.0.0, broom
  1.0.13, lme4, glmmTMB, loo, posterior, brms or bayesplot names a
  `simple*` class (`interop-grep.txt`); broom.mixed is not installed.
  insight is MINOR 2; testthat is NIT 2; future and evaluate build
  their own `simple*` conditions and do not inspect ours.

## Recheck, round 1

Scope: the punch round's changes and the regressions the organizer
named. The lane build was reinstalled; `dev/conditions-rev-verifylib.R`
finds 1,600 top-level functions in the 8 namespaces, 0 differ from the
worktree (`recheck1-verifylib.txt`), and the installed `inst/` trees
match. Scripts are `dev/conditions-rev-*`; logs are in
`dev/conditions-rev-log/`.

### Verdict

Mergeable after one small fix (R1 below), which can go into the final
punch round with the NITs. No BLOCKER and no MAJOR. One design that
fits on base is refused on the lane, and it is narrow.

### Remaining findings

#### MINOR R1. A variable in the parent of an environment `data` is refused

Construction: `dev/conditions-rev-falsealarm2.R base|lane`, case
`data_env_parent_var` (`falsealarm2-cmp.txt`).

    pe <- new.env(); assign("xq", rnorm(60), envir = pe)
    de <- new.env(parent = pe)
    assign("y", d$y, envir = de); assign("x", d$x, envir = de)
    frm(y ~ x + xq, data = de)

Base: fits. Lane: "The model uses `xq`, which is not a column of
`data` and not an object that R finds from the formula".

`check_frame_variables()` looks in `data` with `inherits = FALSE` and
then in the formula environment. For an environment `data`,
`model.frame()` evaluates with `envir = data`, which ignores the
formula environment and searches the parents of `data` instead. So the
claim that the check "never refuses a variable that `model.frame()`
would find" is false for this one container. `?frm` does not document
environment `data` (it names a data frame, tibble, data.table or named
list), so this is MINOR. Fix: for an environment `data`, use
`exists(v, envir = data)` with inheritance and do not consult `env`.

#### NIT R2. The new frame check can name the wrong problem

The refusal fires before `model.frame()`, so an unsupported term whose
argument lives in `data2` reports that argument as an unknown variable:

- `y ~ x + sar(Wm)` with `data2 = list(Wm = ...)`: base says
  `could not find function "sar"`; lane says "The model uses `Wm`,
  which is not a column of `data` ...". `fcor(Vm)` is the same.
- `y ~ .` (unsupported on base too): base says
  `'.' in formula and no 'data' argument`; lane says "The model uses
  `.`, which is not a column of `data`".

Both builds error; only the message is less accurate.

#### NIT R3. The census flags some legitimate shapes

Construction: `dev/conditions-rev-census-overreach.R`
(`recheck1-census-overreach.txt`), which runs `census_sites()` from the
generated core census on one-function snippets. Flagged:

- a local variable named `message` passed to a function
  (`message <- paste(...); frm_stop(message)`) and
  `list(message = message)` (a rule from the first round);
- `get("warning", envir = e)` for an object named `warning`, and a
  local `get <- function(k) ...; get("message")`;
- `identical(h, as.name("stop"))`, a common check on a call head;
- a local `call <- function(x) x; call("stop")`;
- `str2lang("y ~ message(x)")`.

Clean: `switch(k, message = )`, `%in% c("stop", "warning")`,
`obj$get("message")`, `mget()` of a vector,
`getOption("warning.length")` and a `c("stop", "warning")` default. The
real trees contain none of the flagged shapes, and each failure is
loud and exemptible, so this is a NIT.

#### Carried: an unclassed newdata error remains

`predict(fit, newdata = <no grouping column>, allow_new_levels = TRUE)`
still stops at base R's `object 'g' not found`, on base and on the
lane, for a random intercept, a random slope and a binomial mixed
model. brms's `validate_newdata()` accepts the same input
(`newdata-brms-cmp.txt`). `check_newdata_frame()` covers the
fixed-effect `model.frame()` path, not the grouping lookup. It is not a
regression; it belongs on the list of user-reachable unclassed
refusals.

### What held

- **No false alarm from `check_frame_variables()` except R1.**
  `dev/conditions-rev-falsealarm.R` (152 cases) and
  `dev/conditions-rev-falsealarm2.R` (32 cases), each run on both
  builds (`falsealarm-cmp.txt`, `falsealarm2-cmp.txt`): 144 OK on both,
  39 error on both, 1 regression (R1), 0 the other way. The cases:
  - formula-environment variables (global, function local, calling
    function, a formula built in another function), closures used as
    transforms, scalars inside `I()` and `poly()`;
  - `splines::ns`, `stats::poly`, `base::I`, `poly`, `ns`, `I`, `log`,
    `pi`; matrix, Date and character columns;
  - `offset()`, `weights()`, `trials()` (column, literal, environment
    vector), `cens()` right and interval, `trunc()`, `mi()` response
    and predictor, `se()`, `mo()`, `cs()`;
  - `s(x, by = f)`, `s(x, g, bs = "fs")`, `t2()`, `gp()`, `(1 | g)`,
    `(f | g)`, `(1 | g1:g2)`, `(x || g)`, `(1 | p | g)`,
    `gr(g, cov = A)` from `data2` and from the environment, `mm()` with
    and without weights, `car()` from `data2` and from the environment,
    `ar()` and `arma()` with `cov = TRUE`, `cosy()`;
  - `sigma ~`, `hu ~`, `zi ~`, `mixture()`, categorical with a
    character response, `cbind(k, n - k)`, `mvbf()`, `bf() + bf()`,
    `mvbind()`, transformed responses;
  - four nonlinear shapes: a covariate and a random effect in an nlpar,
    a local helper function in the body, an nlpar that shadows an
    environment object, and a data column shadowed by an environment
    scalar;
  - named list, list plus an environment variable, environment, tibble,
    data.table, `POSIXlt` in a list, columns named `get` and `c`, and a
    subset data frame.

  Where a scalar `trials()`, `weights()`, `se()` or `cens()` variable
  is refused, base refused it too ("variable lengths differ").
- **No false alarm from `check_newdata_frame()`.** 48 newdata cases in
  the first harness, and 180 in `dev/conditions-rev-newdata-brms.R`
  (5 models, 12 newdata shapes, and `re_formula = NULL`,
  `re_formula = NA` or `allow_new_levels = TRUE`): 0 cases OK on base
  and refused on the lane, 0 the other way.
  - `allow_new_levels = TRUE` predicts a new grouping level, a new
    `mm()` member and a new level in a multivariate fit.
  - `re_formula = NA` works without the grouping column.
  - A missing response, missing `weights`, `cens` and `trials` columns,
    character for factor, an `NA` level, `droplevels()` and an unused
    extra factor level all predict.
  - A new fixed-effect level is refused on both builds, and brms
    refuses it too. On a slope-only factor with `re_formula = NA`, all
    three accept it.
  - `conditional_effects()`, emmeans and marginaleffects
    (`avg_predictions`, `avg_slopes`, `avg_comparisons`) agree on both
    builds for gaussian, binomial, weighted, censored, `s(by)`, `poly`,
    mixed and random-slope fits.
- **Rd examples.** All 163 topics that both builds have, across the 8
  packages, run to completion on both, `\donttest` sections included
  (`dev/conditions-rev-rdex.R`, `rdex/cmp.tsv`). The lane's extra
  topic, `frmtmb-conditions`, also runs.
- **`frm_match_arg()` matches `match.arg()`.** 60 comparisons
  (`recheck1-matcharg.txt`): 20 inputs (missing, `NULL`, exact,
  partial, ambiguous, empty, wrong case, the full vector, a
  permutation, a vector with a bad element, zero length, numeric,
  `NA`, factor, list) through a default-vector formal, the same with
  `several.ok = TRUE`, and explicit choices. 0 differ.
  - It reads the right frame in an S3 method, after `NextMethod()`
    (the default method's formals, as `match.arg()` does), in a nested
    function with a shadowing formal, through `do.call()`, after the
    argument is reassigned, and with explicit choices in a helper.
  - Inside an anonymous `vapply()` function, both functions error.
  - Valid values pass through 14 real sites.
- **Priors.** 25 prior strings through `set_prior()`, and 4 fits with
  priors (with bounds, on `sd`, `sigma` and nlpars), give the same
  outcome on both builds. Whitespace, `1e3`, `.5` and every accepted
  density parse. `constant()`, `lkj_corr_cholesky`, `horseshoe`,
  `R2D2`, `uniform`, named arguments and `flat` are refused on both
  with the same text. `prior_dist_params` has the base arities.
- **Package tag.**
  - `frm_family_package()` gives `frmtmb.eam`, `frmtmb.latent` and
    `frmtmb` for `wiener()`, `hmm()` and `gaussian()`.
  - The tag survives `saveRDS()`/`readRDS()`, `$<-`, `modifyList()`,
    the finalized family in `fit$spec`, `saveRDS()` of the fit, and
    `update()`. A deviance refusal on a fitted wiener model is a
    `frmtmb_eam_error` before and after `readRDS()`.
  - The tag is lost under `[` subsetting and in `mixture()`, which
    report `frmtmb`.
  - A user family made at top level carries no attribute.
  - `identical()` of two `gaussian()` or two `wiener()` objects is
    `FALSE` because of their closures, and `all.equal()` is `TRUE`.
  - `ndt_apply()`'s non-refusal path is textually identical to base
    (diff). The lookup is an argument of `frm_stop()`, so it is
    evaluated only on the refusal.
- **Census.** The worker's 26 plants reproduce exactly
  (`recheck1-punch-census.txt` against `census.txt`). The census still
  fails on all 8 base trees (`dev/conditions-rev-census-base.R`). The
  anchor is asserted after the source-tree test, so a renamed
  `objective.R` or `<pkg>-package.R` fails. The source-tree test still
  skips an installed package, which is the house pattern.
- **Suites.** 133 test files on the lane build, one per process: all
  98 extension files, plus 35 core files that the round touches
  (census, conditions, interop, sugar, nl-lexical, tabular inputs,
  argument refusal, brms families and methods, priors, newdata,
  conditional effects, famgaps, custom family, structure, dates,
  autocor, ps). 1,376 blocks, 9,508 passes, 0 failures, 0 errors, 15
  skips and 0 `compileCode` lines (`recheck1-spot/summary.txt`). The
  pass count of every one of the 133 files equals the worker's
  `sweep-lane`.
