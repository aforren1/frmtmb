# autotest triage

This record covers the run of rOpenSci's
[`autotest`](https://docs.ropensci.org/autotest/) over `frmtmb` before
the rOpenSci statistical-software submission. It says what the tool
reported, how each report was triaged, what was fixed, and what was
suppressed with the reason.

The tool is `autotest` 0.1.1.010 with `typetracer` 0.2.5 on R 4.6.1,
run with `NOT_CRAN=true`. `dev/autotest-run.R` reproduces the run.

## What autotest does

`autotest` parses each `.Rd` file, runs the example code under
`typetracer` to learn the class, storage mode and length of every
argument each function actually receives, then calls the function again
with one argument mutated. The mutations are listed by
`autotest_types()`: a length-1 argument given length 2, a logical given
a character or an integer, a character given a different case, a
`data.frame` given as a `tibble` or a `data.table`, a vector given as a
list column, and so on. Each mutation carries an expectation. Some say
the result must not change; some say the call must refuse.

The interesting reports are therefore of two kinds. Either a mutation
that should have been harmless was not, or a mutation that should have
been refused was accepted.

## Read this first: autotest tests the INSTALLED package

The first run of this exercise produced 450 rows, 251 of which looked
like a single spectacular artifact: every `frm()` example came back as

```
:quote(structure(list(formula = y ~ x, pforms = list(),
  : unused argument (data2 = base::quote(list()))
```

The obvious reading was that `autotest` had failed to rebuild the call.
The real cause was worse, and is worth stating plainly for anyone else
running this tool.

`typetracer::pre_install()` installs the source package into a
temporary library and loads it from there, which is correct. At the end
of tracing, `typetracer::reload_pkg()` deletes that temporary library.
Every mutation `autotest` then performs resolves the package through
`.libPaths()`, which means **whatever version happens to be installed
on the machine**. Here that was `frmtmb` 0.28.0 from the user library,
eleven releases behind the source. The mutations were run against a
package whose functions do not all have a `data2` argument, so a third
of the rows were version skew rather than diagnostics.

Nothing in the output says which version was tested. The check that
catches it is one line:

```r
cat(as.character(utils::packageVersion("frmtmb")),
    exists("check_flag", asNamespace("frmtmb")), "\n")
```

after `autotest_package()` returns. The run recorded here pins a
library built from the exact source under test at the front of
`.libPaths()` and keeps `reload_pkg()` from deleting the traced
install, so both the before and after runs test the source they name.
The numbers below are from those corrected runs. The uncorrected run is
not reported, because none of it can be trusted.

## Upstream workarounds

None is a `frmtmb` defect. All five are recorded because the run cannot
be reproduced without them. Four are Windows bugs in `typetracer`; one
is a robustness gap in `autotest`. They live in `dev/autotest-run.R`
and are applied with `assignInNamespace()`.

**1. `typetracer::insert_counters_in_tests()`** pastes the trace
directory and the package directory into R string literals. On Windows
those are backslash paths, so `C:\Users\...` contains the escape
sequence `\U` and every test file of the traced copy fails to parse:

```
Error: '\U' used without hex digits in character string
  (test-aliased-grouping.R:12:27)
```

The run aborts there. The workaround writes both paths with forward
slashes. Worth knowing: this function rewrites test files in place, but
it works on the pre-installed copy in the `typetracer` cache, so the
worktree is untouched. That was verified with `git status` after the
first aborted run.

**2. `typetracer::reload_pkg()`** contains

```r
ifelse(grepl(tempdir(), lib_path), lib_path, tempdir())
```

which uses a Windows temp path as a REGULAR EXPRESSION. The backslashes
read as back references and the call dies with `invalid regular
expression ..., reason 'Invalid back reference'`, at the end of tracing
and before `autotest` has seen a single trace. The workaround compares
with `fixed = TRUE`.

**3. `typetracer::reload_pkg()` again**, for the reason in the previous
section: it deletes the temporary library holding the traced build. The
workaround keeps it.

**4. `typetracer::trace_package_tests()`** calls
`testthat::test_package()`, which throws when any test fails. Under the
rewriting in (1) some tests do fail, and one failure aborts the whole
`autotest` run rather than the test pass alone. The workaround degrades
to examples-only tracing.

This one has a cost worth stating. Under tracing the test pass takes
roughly fifty minutes, and `autotest` mutates only traces whose source
is `"examples"` (`autotest_package()` filters on that), so the test
pass buys parameter coverage and no mutations. The recorded run is
examples-only. The package's own suite is run separately and in full.

**5. `autotest:::pass_one_rect_as_other()`** calls

```r
do.call(data.frame, x$params[[x$i]], quote = TRUE)
```

which needs the parameter it classified as rectangular to be a list. A
matrix is rectangular and is not a list, so the call dies with `second
argument must be a list` and takes the whole run with it. The
workaround skips that one mutation.

## Results

| | Rows | Diagnostics |
| --- | --- | --- |
| Before | 455 | 106 |
| After | 412 | 60 |

Per class, before and after:

| Test class | Before | After | Delta |
| --- | --- | --- | --- |
| `single_par_as_length_2` | 38 | 13 | -25 |
| `subst_int_for_logical` | 16 | 5 | -11 |
| `subst_char_for_logical` | 14 | 4 | -10 |
| `random_char_string` | 8 | 5 | -3 |
| `par_matches_docs` | 2 | 0 | -2 |
| `single_char_case` | 15 | 18 | +3 |
| `negate_logical` | 221 | 226 | +5 |
| `par_is_demonstrated` | 69 | 69 | 0 |
| `par_is_documented` | 23 | 23 | 0 |
| `return_desc_includes_class` | 13 | 13 | 0 |
| `return_successful` | 11 | 11 | 0 |
| `int_range` | 1 | 1 | 0 |
| `vector_to_list_col` | 1 | 1 | 0 |

Two rises, both understood.

`single_char_case` gains exactly the three entries `random_char_string`
loses: `frmtmb_family(type =)`, `custom_family(type =)` and
`frm(dry_run =)`. Giving each a fixed set of permitted strings is what
stops an unrecognized value being accepted, and it also makes the
upper-cased value refuse, which is by definition case dependence. They
join the class already suppressed below for the reason every other
identifier in the package is.

`negate_logical` gains five rows of `frm(verbose = TRUE)` stage output.
The count of stage lines varies with how many restarts the optimizer
takes, so this number moves a little between runs.

Triage of the 455 baseline rows:

| Bucket | Rows |
| --- | --- |
| (a) Real defect, fixed | 46 |
| (b) Missing documentation, fixed | 2 |
| (c) False positive or not applicable | 278 |
| (d) Deferred | 129 |

By test class:

| Test class | Rows | (a) | (b) | (c) | (d) |
| --- | --- | --- | --- | --- | --- |
| `negate_logical` | 221 | 7 | 0 | 196 | 18 |
| `par_is_demonstrated` | 69 | 0 | 0 | 0 | 69 |
| `single_par_as_length_2` | 38 | 25 | 0 | 9 | 4 |
| unnamed | 23 | 0 | 0 | 5 | 18 |
| `par_is_documented` | 23 | 0 | 0 | 23 | 0 |
| `subst_int_for_logical` | 16 | 0 | 0 | 16 | 0 |
| `single_char_case` | 15 | 0 | 0 | 15 | 0 |
| `subst_char_for_logical` | 14 | 10 | 0 | 4 | 0 |
| `return_desc_includes_class` | 13 | 0 | 0 | 0 | 13 |
| `return_successful` | 11 | 0 | 0 | 5 | 6 |
| `random_char_string` | 8 | 4 | 0 | 4 | 0 |
| `par_matches_docs` | 2 | 0 | 2 | 0 | 0 |
| `int_range` | 1 | 0 | 0 | 0 | 1 |
| `vector_to_list_col` | 1 | 0 | 0 | 1 | 0 |

The largest bucket, 186 rows, is `negate_logical` flipping
`frm(verbose = FALSE)` to `TRUE`. A verbose fit then narrates its
stages, which is what the flag is for, and `autotest` files every
message as its own row.

## The one defect class

Everything in bucket (a) is the same fault wearing different names. A
scalar argument was read by an idiom that turns a wrong value into a
plausible one instead of into a refusal.

`isTRUE(x)` is the main offender. It answers `FALSE` for a length-2
logical, for a string, for an integer and for `NA`. A flag set by
mistake therefore selected the OTHER option, and the fit ran to
completion reporting a model the user did not ask for. Nineteen
arguments were read that way.

The second offender is recycling. A confidence level becomes a normal
quantile that is added to a vector of estimates, so a length-2 level
makes a length-2 quantile that recycles: one table came back with its
odd rows at 90% and its even rows at 95%, and no column recorded it.
That is a wrong answer, not a poor error message.

The third is `[[` on a list, which indexes RECURSIVELY when given a
character vector and by POSITION when given an integer. The link
registry is a list of lists, so `link = c("log", "name")` handed back
the string `"log"` as though it were a link object, and `link = 1L`
selected the identity link without a word.

The fix is one family of checks in `R/utils.R`, each raising a single
condition message, called at the function entries.

## (a) Real defects fixed

`autotest` reached 46 of these rows. The probes described at the end
reached more, on functions whose arguments the examples do not vary.

| Function | Argument | Old behavior | New refusal |
| --- | --- | --- | --- |
| `bf()` | `nl` | `nl = "yes"` built a LINEAR model in silence | `` `nl` must be TRUE or FALSE, not character "yes" `` |
| `mvbf()` | `rescor` | fitted without residual correlation | `` `rescor` must be TRUE or FALSE, not ... `` |
| `set_rescor()` | `rescor_value` | returned `rescor = FALSE` for any non-`TRUE` value | `` `rescor_value` must be TRUE or FALSE, not ... `` |
| `frmtmb_control()` | `profile`, `sparse_x`, `autoscale` | the option was silently off | `` `profile` must be TRUE or FALSE, not ... `` |
| `frmtmb_control()` | `restarts` | `"x"` reached `seq_len()`; `c(1, 2)` used only its first element | `` `restarts` must be a single whole number of at least 0, not ... `` |
| `frmtmb_control()` | `grad_tol` | a string was compared as a string at fit time | `` `grad_tol` must be a single finite positive number, not ... `` |
| `frmtmb_control()` | `optimizer` | a length-2 value reached `switch()` and reported `EXPR must be a length 1 vector` once per element | `` `optimizer` must be a single optimizer name or a function, not ... `` |
| `frmtmb_control()` | `optCtrl` | a non-list reached the optimizer | `` `optCtrl` must be a list of options for the optimizer ... `` |
| `frm()` | `REML`, `se`, `quadrature` | `the condition has length > 1`; `argument is not interpretable as logical` | `` `REML` must be TRUE or FALSE, not ... `` |
| `frm()` | `dry_run` | `"nonsense"`, `TRUE` and `c("spec", "frame")` all fitted the full model | `` `dry_run` must be one of "spec", "frame", "objective", not ... `` |
| `frm()` | `start` | `start = "a"` was iterated over an empty `names()` and ignored, so the fit ran from the DEFAULT start and reported it as the user's | `` `start` must be a named list, e.g. start = list(beta = c(0, 1)), not ... `` |
| `frm()` | `control` | `control = 5` gave `$ operator is invalid for atomic vectors`; `control = list()` gave `argument must be coercible to non-negative integer` from inside the optimizer | `` `control` must be a list from frmtmb_control(), not ... `` and `` `control` must come from frmtmb_control(); this list is missing ... `` |
| `frm()` | `na.action` | `attempt to apply non-function` | `` `na.action` must be a function such as stats::na.omit, or its name as a string, not ... `` |
| `frm()` | `priors` | `is.list(priors) is not TRUE` | `` `priors` must be a set_prior() specification or a named list of prior objects, not ... `` |
| `assemble_frame()` | `data` | `data = NULL` gave `object 'y' not found`, which reads as a formula typo | `` `data` is NULL: frm() needs the data frame holding the model variables ... `` |
| `get_link()` | `name` | `c("log", "name")` returned the STRING `"log"` as a link; `1L` selected the identity link | `A link must be named by a single string, e.g. link = "logit", not ...` |
| `confint()` | `level` | a length-2 level recycled and reported DIFFERENT coverages on different rows of one table; `1.5` gave a table of `NaN`; `-1` gave `Inf` | `` `level` must be a single number strictly between 0 and 1, not ... `` |
| `confint()` | `nsim` | reached the bootstrap as a length | `` `nsim` must be a single whole number of at least 1, not ... `` |
| `confint()` | `parm` | `parm = 5` gave `subscript out of bounds`, beside an excellent message for an unknown NAME | `` `parm` must be a character vector of parameter names, or NULL for all of them, not ... `` |
| `conditional_effects()` | `prob` | the same recycling, alternating the band coverage along the grid | `` `prob` must be a single number strictly between 0 and 1, not ... `` |
| `conditional_effects()` | `resolution`, `ndraws`, `profile_points` | a length-2 value used only its first element | `` `resolution` must be a single whole number of at least 1, not ... `` |
| `conditional_effects()` | `surface`, `conditions` | read through `isTRUE()` and an empty `names()` | flag and named-list refusals |
| `predict()` | `se.fit`, `allow_new_levels` | `invalid argument type` | `` `se.fit` must be TRUE or FALSE, not ... `` |
| `predict()` | `newdata` | a list gave `non-conformable arrays` from inside the design builder | `` `newdata` must be a data frame, or NULL to predict on the training data, not ... `` |
| `predict()` | `re.form` | `re.form = "x"` returned the POPULATION prediction without saying so | `` `re.form` must be NULL to keep every random effect, NA to drop them all, or a one-sided formula ... `` |
| `simulate()` | `nsim`, `censored` | `invalid 'length' argument`; `nsim = 2.5` silently truncated | count and flag refusals |
| `frm_simulate()` | `nsim` | as above | count refusal |
| `frm_bootstrap()` | `nsim` | as above, raised from inside the bootstrap loop | count refusal |
| `frm_multiple()` | `level` | recycled down the pooled coefficient table | coverage refusal |
| `diagnose()` | `quiet` | `invalid argument type`; `the condition has length > 1` | flag refusal |
| `vcov_cluster()` | `full` | `invalid argument type` | flag refusal |
| `ranef()` | `condVar` | `argument is not interpretable as logical` | flag refusal |
| `lca()` | `na.rm` | `na.rm = "yes"` changed which subjects were kept, in silence | flag refusal |
| `cox()` | `df`, `degree` | a length-2 value sized the I-spline baseline basis one way and recorded it another | count refusal |
| `cox()` | `intercept` | read through `isTRUE()` | flag refusal |
| `frmtmb_family()`, `custom_family()` | `type` | any string accepted, giving a family that behaved as none of the four kinds | `` `type` must be one of "continuous", "discrete", "ordinal", "categorical", not ... `` |
| `frmtmb_family()`, `custom_family()` | `drop_intercept` | read through `isTRUE()` | flag refusal |
| `prior_normal()`, `prior_t()` | `location` | a length-2 location built a prior that recycled against the parameter block | `` `location` must be a single finite number, not ... `` |
| `prior_normal()`, `prior_t()` | `scale`, `df` | `stopifnot()` reported `scale > 0 is not TRUE`, which names the test rather than the argument | `` `scale` must be a single finite positive number, not ... `` |
| `check_custom_family()` | `tol` | a length-2 tolerance compared elementwise | `` `tol` must be a single finite positive number, not ... `` |

No correct code path changed. Every entry is a guard at a function
entry that refuses an input which used to either fail somewhere
unhelpful or return a result the caller did not ask for.

### Two things the test suite caught, and both were right to

The new coverage check was first named `check_coverage()`, which is
already the name of an internal in `R/simulate-new.R` that checks
whether supplied parameters cover the model's slots. The new one is
`check_probability()`.

More interesting: `frmtmb_control(optimizer = )` was first given a full
value check against `c("nlminb", "optim")`. Three tests failed, and
they were asserting a deliberate design decision. An unknown optimizer
name is refused at FIT time, so the message can carry the family and
mode it was raised from:

```
Unknown optimizer 'nope' (use "nlminb", "optim", or a function)
  (raised while fitting: gaussian, ML, nope)
```

`frmtmb_control()` cannot know any of that. The value check was
reverted and only the SHAPE of `optimizer` is checked in the
constructor, which still fixes the length-2 case and leaves the tested
behavior alone. This is the difference between validation and a
behavior change, and it is exactly the line the exercise was not
supposed to cross.

## (b) Missing documentation, fixed

`par_matches_docs` reported that `frm_simulate()` and `get_prior()`
describe their `data` argument without naming its class. Both said
"Model data"; both now say "A data frame of model data". `frm_sample()`
got the same correction. 2 rows.

Two further documentation corrections came from the probes:

* `?frm` now states what `data` accepts. A `tibble`, a `data.table` and
  a plain named list all work, because each reaches
  `stats::model.frame()` unchanged. A matrix column is a supported
  model variable. A list column is not one, though it may sit unused in
  the frame.
* The `G2.16` block in `R/fit.R` claimed that `Inf`, `-Inf` and `NaN`
  are each rejected with their own message. The check is
  `any(!is.finite(y) & !is.na(y))` and `is.na(NaN)` is `TRUE`, so `NaN`
  takes the `na.action` route and never reaches that message. The claim
  now describes the division the code actually makes, which is
  `stats::lm()`'s: between values that are undefined and values that
  are absent, not between the `Inf` and `NaN` spellings.

## (c) False positive or not applicable

There is no comment pragma for `autotest` in the way `srr` has
`@srrstatsNA`. The sanctioned mechanism is the `test_data` argument:
`autotest_package(test = FALSE)` returns one row per test it would run,
each carrying a `test` flag, and setting flags to `FALSE` and handing
the table back runs everything else. The exclusions and their reasons
are in `dev/autotest-run.R`.

### The verbose flag doing its job (186 rows)

`negate_logical` flips `frm(verbose = FALSE)` to `TRUE`. The fit then
prints one `message()` per stage, which is the documented purpose of
the flag, and `autotest` files each message as a row:

```
frmtmb: parse [0.00s]
frmtmb: frame [0.02s]: 60 obs, 2 linear predictors
frmtmb: tape [0.00s]: 3 outer, 0 inner parameters
```

`verbose` is also the one flag deliberately left unvalidated.
`verbose_level()` accepts `TRUE`, `FALSE` and integer levels and maps
anything else to silent, and `tests/testthat/test-verbose.R` asserts
that `verbose_level(list(verbose = "no")) == 0L`. A `check_flag()`
there would have been a behavior change, not validation.

### Identifiers are case-sensitive (15 rows)

`single_char_case` upper-cases a character argument and expects the
same result. The arguments reported are link names, prior classes,
`vcov(type =)`, `frm_ode(method =)` and `frm_ode(output =)`, dpar name
vectors and `mixture_mvn(model =)`. Each is an identifier, not free
text, and R matches identifiers by case: `stats::binomial(link =
"LOGIT")` errors for the same reason, and `match.arg()` is
case-sensitive throughout base R. Accepting `"LOGIT"` would be a new
feature.

The package already makes the one case exception its users expect:
`confint(method = "Wald")` is accepted beside `"wald"`, because brms
spells it with the capital.

### Combined roxygen parameter blocks (23 rows)

`par_is_documented` reported 23 undocumented parameters. Every one is
documented, in a combined block. `@param lb,ub`,
`@param location,scale,df`, `@param nsim,seed`, `@param atol,rtol`,
`@param feature_a,feature_b`, `@param formula,...` and
`@param data2,start,control,na.action,REML` each reach the `.Rd` as a
single `\item{lb, ub}{...}`. `autotest` looks for an `\item` naming one
parameter and does not split the combined form.

This was checked, not assumed: the `\arguments` block of every reported
`.Rd` was read and each reported parameter appears in one.

### Deliberate refusals autotest expects to succeed (16 rows)

`subst_int_for_logical` expects `profile = 1L` to keep working. The
package now refuses it. The reason is the fix itself: these flags were
read through `isTRUE()`, which calls `1L` FALSE, so the integer
spelling silently selected the opposite option. Refusing an integer is
a smaller surprise than accepting it and doing the reverse.

### Arguments autotest misclassified (21 rows)

* `set_prior(lb =, ub =)` are numeric bounds whose "unset" marker is
  `NA`. `autotest` reads the class off the default and calls them
  logical, so the logical mutations test a type they never had.
* `frm_bootstrap(re.form =)` defaults to `NA` for the same reason. It
  takes `NA`, `NULL` or a formula.
* `custom_family(dpars =, primary_dpars =)` are character VECTORS by
  design, one entry per distributional parameter. The examples happen
  to pass one name.
* `set_prior(coef =, dpar =, group =)` name a coefficient, a
  distributional parameter and a grouping factor of a model
  `set_prior()` has not seen. Any string is a well-formed name; whether
  it matches is checked when the prior is applied to a fit, which does
  refuse an unmatched name and lists the available ones.
* `frmtmb_family(sim_refusal =)` is the TEXT of a refusal message, so
  changing its case or length changing the output is correct.
* `check_custom_family(y =)` is the response vector a candidate density
  is evaluated at, passed to that density directly rather than through
  a model frame, so `vector_to_list_col` has nothing to test.

### Calls autotest could not build (17 rows)

Genuine reconstruction limits, now that version skew is out of the way.
`autotest` calls an S3 method by its full name
(`could not find function "vcov.frmtmb_fit"`, same for
`simulate.frmtmb_fit`), and it drops a required argument
(`argument "formula" is missing, with no default` for `frm()`,
`argument "hypothesis" is missing` for `hypothesis()`).

### Rectangular class replacement

Not reported by the corrected run: converting `data` to a `tibble` or a
`data.table`, and adding a class on top of `data.frame`, all work.
Probed separately, stripping `data.frame` from the class leaves an
object `as.data.frame()` refuses with `cannot coerce class '"myframe"'
to a data.frame`, which is correct: an object that no longer claims to
be rectangular cannot be a model frame.

## (d) Deferred

129 rows, and no code change for any of them.

### Examples do not exercise every parameter (69 rows)

`par_is_demonstrated` reports a parameter that appears in no example.
The reported ones are the optional arguments of `frm_ode()`,
`frm_sample()`, `custom_family()`, `frm()` and the rest of the surface
in ones and twos.

Every exported object already has a worked example, and the check that
matters for the submission, `fn_without_example`, reported nothing.
Demonstrating each remaining optional argument would mean adding a
fitted model per argument to example blocks that already run at
submission time. That trades worse documentation for a better score, so
it is declined. Arguments that are subtle rather than merely optional
are demonstrated in the vignettes.

### rstan sampler diagnostics (14 rows)

`frm_sample()` examples run short chains, and `rstan` reports low ESS,
one divergent transition and R-hat up to 1.07 through `warning()`.
These are the sampler's own diagnostics on a deliberately small
example, not `frmtmb` conditions.

### Correct refusals of degenerate calls (6 rows)

`return_successful` rows where the mutated call cannot succeed and the
package says so by name: `bf()` refusing `nl = TRUE` with no parameter
formula, `mvbf()` refusing a single response, `loo_compare()` refusing
one model, `check_custom_family()` reporting a tape failure for a
density assembled from mutated arguments. The same shows up in
`negate_logical`, where negating `quadrature` or `REML` hits the
documented incompatibilities with `hmm()` and the mixture families.
Each refusal is the right one.

### Return-class descriptions (13 rows)

`return_desc_includes_class` wants the `@return` block to name the
class string. The blocks describe the structure the object holds
instead, which is more useful to a reader. Logged rather than churned.

### Integer range probes (1 row)

`frm_ode(n_ss = 20)` with "does not respond appropriately for
specified/default input". The test looks for a documented permissible
range; `n_ss` is a steady-state iteration cap with no upper bound worth
documenting.

### `frm_bootstrap(seed =)`, `frm_simulate(seed =)`, `frm_compat()`

Reported by `single_par_as_length_2`. `seed` is handed to `set.seed()`,
which takes the first element itself; `frm_compat()` looks a feature
name up in a table and returns nothing for a vector. Neither can
produce a wrong number, and guarding them would add condition messages
for no gain.

## The frmtmb-specific probes

`autotest` mutates only what the examples demonstrate, so a set of
inputs was probed directly as well. These ran against the worktree
source under `pkgload::load_all()`.

### Rectangular input

| Input | Result | Verdict |
| --- | --- | --- |
| `tibble` as `data` | fits | correct |
| `data.table` as `data` | fits | correct |
| `data.frame` with an extra class | fits | correct |
| class replaced, no `data.frame` left | refused by `as.data.frame()` | correct |
| `list` as `data` | fits | correct, now documented |
| `matrix` as `data` | refused by `model.frame()` | correct |
| zero-row `data.frame` | `` `data` has no rows; nothing to fit `` | correct |
| `NULL` as `data` | was `object 'y' not found` | fixed, see (a) |

### List columns

A list column that is not a model variable is carried through and
ignored; the fit is unaffected. A list column named in the formula is
refused by `stats::model.frame()` with `invalid type (list) for
variable 'lc'`, which names the offending column. That refusal is
correct and specific, so no guard was added. What was missing is the
statement in `?frm` that a list column cannot be a model variable while
a matrix column can, which is now documented.

### Single-row data

`frm(bf(y ~ x) + gaussian(), data = d[1, ])` fits, and is not silent
about being degenerate. The rank check drops the aliased column with a
message and the optimizer reports non-convergence as a warning:

```
Fixed-effect design of 'y.mu' is rank deficient; dropping column(s): x
Optimizer did not report convergence: function evaluation limit reached
```

That is the documented `RE2.4` and `RE3.0` behavior working on the
smallest possible input, so nothing was changed.

### Factor columns given as character

A character grouping variable and a character fixed-effect predictor
both fit, matching `stats::lm()` and `lme4`. Correct.

### Integer and double responses

An integer response to `gaussian()` fits. A whole-number double
response to `poisson()` fits. A non-integer double response to
`poisson()` is refused with `poisson: response must be non-negative
integers`. A logical response to `bernoulli()` fits. All correct.

### Matrix columns

Matrix columns are a feature and stay one. `frm(bf(y ~ M) + gaussian())`
with `M` an n-by-3 matrix column fits, and `dry_run = "frame"` returns
the assembled design. No mutation of a matrix column was converted into
a refusal, and the one `autotest` mutation that would have crashed on a
matrix is the upstream bug in workaround (5).

The one refusal in this area is a two-column `cbind()` response to
`binomial()`, refused by the response check. That is the brms grammar,
where the trials count is an addition term (`y | trials(n)`) rather
than a second response column, so the refusal is correct.

### Missing values

`NA` in a predictor is handled by `na.action`: dropped with a message
under the default `na.omit`, refused under `na.fail`. `Inf` in the
response is refused separately with its own message. `NaN` in the
response is dropped by `na.action`, because `is.na(NaN)` is `TRUE`.

No `NA` path was double-guarded. The only change here is the `G2.16`
wording described under (b).

## Reproducing

```
Rscript dev/autotest-run.R
```

Results land in `$AUTOTEST_OUT`, or a directory under `tempdir()`. The
script prints how many candidate tests were suppressed, and the version
of `frmtmb` it actually tested. Read that version line before reading
anything else, for the reason in the second section.

The first run is slow, and the reason is worth knowing before you start
it. `autotest_package()` re-traces every example on each call, and the
script calls it twice: once with `test = FALSE` to learn which tests
exist, so the suppressions can be applied per function, and once with
`test = TRUE` to run them. The `test = FALSE` pass is much the slower
of the two on this package, because it emits a row for every candidate
test rather than only the failures, and `autotest` accumulates them
with `rbind()` and de-duplicates the whole frame once per trace. It ran
past eighty minutes here without finishing.

The plan from that pass is cached as `autotest-plan.rds`, so a repeat
run skips it and takes about ten minutes. Delete the cache after
changing an example or a signature.

If all you want is the reports, skip the script:
`autotest_package(PKG, test = TRUE)` with the five patches applied
takes ten minutes and is what produced the before and after numbers
above.

The numbers reported above come from unsuppressed runs
(`autotest_package(test = TRUE)` with no `test_data`), so that the
before and after totals count the same things. The suppression script
is for reading the surviving reports, not for producing those totals.
