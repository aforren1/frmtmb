# Lane wt-adefects: the silent wrong answers brms's suite found

Item 2.6f, the silent half. Worktree
`C:\Users\adf44\source\r\frmtmb-wt-adefects`, base `f8b45ef` (frmtmb
0.60.0, frmtmb.sample 0.8.0). The other half, the post-fit return
shapes, is lane `wt-shapes`; sections 10 and 11 say where the two
meet.

The specification is `dev/brmsport-findings.md` section 5 and the
ranked list in `dev/reviews/20260917-brmsport.md`. Seven items were
handed to this lane. Five were live defects and are fixed, one was
already fixed at 0.60.0 and is recorded as dissolved, and one
(`predict(newdata)` without a grouping column) turned out to be two
faults rather than one.

## 0. Where everything is

| script | log | what |
|---|---|---|
| `dev/adefects-install.R` | `install1.txt` to `install4.txt` | roxygenise and install into the lane library `C:/Users/adf44/source/r/adefects-lib`, ONE CHILD PROCESS PER PACKAGE |
| `dev/adefects-repro-base.R` | `repro-base.txt`, `repro-lane.txt` | every defect on the base build and on this one, with brms 2.23.0 on the same call |
| `dev/adefects-evidence.R` | `evidence.txt` | what `fit$data` costs, which refusal an autocorrelation term reaches first, the D5 sweep |
| `dev/adefects-falsealarm.R` | `falsealarm-base.txt`, `falsealarm-lane.txt`, and the two `.tsv` arms | the false-alarm measurement, harvested from brms's own suite and this repository's |
| `dev/adefects-fadiff.R` | `fadiff.txt` | the base-against-lane diff of those two arms, which is what attributes a verdict to this lane |
| `dev/adefects-hypscan.R` | `hypscan.txt` | every `hypothesis()` call in the repository whose string states no relation |
| `dev/adefects-verdicts.R` | `port-record.txt`, `port-ledger-1.txt` | the ported tier's verdict moves, applied idempotently |
| `dev/adefects-run.R`, `dev/adefects-suite.sh`, `dev/adefects-gated.sh` | `suite-*.txt`, `gated.txt` | one R process per test file |
| `dev/adefects-paste-ledger.R` | | pastes the generated ledger summary into `dev/brmsport-findings.md` |
| `dev/adefects-instrument.R` | `instrument.txt` | the lane library IS this worktree's source, function by function |
| `dev/adefects-p1-dollarcost.R` | `p1-dollar.txt` | punch round 1: what a `$.frmtmb_fit` would cost, what the duplicate costs on disk |
| `dev/adefects-p1-d4inverse.R` | `p1-d4inverse.txt` | the four D4 pins, run against the design they must fail on |
| `dev/adefects-all.sh`, `dev/adefects-check.sh` | `runs-complete.txt`, `check-*.txt` | the final pass, and `R CMD check --as-cran` once per package |

**The instrument.** The lane library is built from this worktree:
1045 of 1045 top-level functions in `R/` and 199 of 199 in
`extensions/frmtmb.sample/R/` deparse identically to the installed
namespace (`dev/adefects-instrument.R`, `instrument.txt`). Nothing in
this lane changes a version field, so a version string could not have
seen a stale install.

Data seed 20260917 wherever data are simulated, which is the seed
`dev/brmsport-defects.R` uses, so the numbers below are comparable to
the ones that record carries. The library stack is the lane build, then
the round's reference build `rellib-r3` (read only), then the user
library. The StanHeaders pin is gone; `dev/adefects-run.R` and
`dev/brmsport-run.R` assert the user Makevars `-std=gnu++17` flag and
`tmbstan >= 1.2.1` instead, as `dev/release/run-tests.R` does.

brms's source was re-fetched into `dev/brms-suite/` from the URL in
`dev/brms-suite-audit.md` section 1. Its sha256 is
`b5f5bb5604ec3f87b3ad99f0da4e6a37d8c5488221e8abe459171340992c37a5`,
the recorded value.

## 1. D1: an expression as the time index

`ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` evaluated whatever was
written in the `time` position against the model frame, and the values
became the time points. brms refuses it:
`as_one_variable(deparse(substitute(time)))`, "Cannot coerce 'x + t' to
a single variable name".

Base build, 24 rows, seed 20260917 (`repro-base.txt`):

| term | logLik | control |
|---|---|---|
| `ar(x + t, g, cov = TRUE)` | -30.33421 | -30.88513 for `ar(t, g)` |
| `ar(t - 10 * x, g, cov = TRUE)` | -31.00618 | the same |
| `ma(x + t, g, cov = TRUE)` | -29.84347 | -30.49819 for `ma(t, g)` |
| `arma(x + t, g, cov = TRUE)` | -29.43557 (df 5) | |
| `cosy(x + t, g)` | -30.86409 | |
| `unstr(x + t, g)` | 40.3157, **df 279** | |

The `unstr()` row is the one that shows what was happening: the 24
distinct values of `x + t` became 24 time levels and the term estimated
276 free correlations on 24 observations. brms refuses all five with
the same sentence.

**Fixed** in `autocor_time_name()` (`R/autocor.R`). `time` must be a
symbol; `NA` and an omitted `time` are unchanged. A numeric literal and
a string are refused too, which is also what brms does
(`all_vars("1")` is empty, so `as_one_variable()` rejects it).

The check runs BEFORE the `cov = TRUE` one, which matters for the
ported tier: `brm:106` is written without `cov`, and it now reaches the
grammar refusal rather than the `cov` one (`evidence.txt`).

## 2. D2: the `gr =` grammar

Two faults under one heading.

**`gr = g1/g2` grouped by the quotient.** On numeric group codes the
expression was evaluated, so the series (1, 1) and (2, 2) both had
ratio 1 and merged into one residual block: logLik -30.87554 against
-30.88513, with nothing said. brms refuses the term,
"Illegal grouping term 'g1/g2'. It may contain only variable names
combined by the symbol ':'".

**`gr = g1:g2` was R's sequence operator.** R's `:` gives the
interaction for two factors and a sequence for anything else, so on
numeric codes the whole grouping vector was `1`, of length one. The
base build did not fit it: it reported "time points within groups must
be unique; group '1' has 3 rows at time '1'", with R's warning
"numerical expression has 24 elements: only the first used". Loud, but
for the wrong reason, and the same spelling on factor codes
(`gf1:gf2`) fitted correctly at -30.88513. brms accepts `gr = g1:g2` on
numeric codes and CROSSES the two.

**Fixed** in `autocor_gr_vars()` and `autocor_gr_value()`
(`R/autocor.R`). The grammar is brms's: variable names combined by `:`.
Crossing is done by pasting the components with `_`, which is what
brms's `combine_groups()` does, and each component is checked for
missing values before the paste, because `paste(NA, 1, sep = "_")` is
the string `"NA_1"` and would hide one. A single variable keeps its own
type, so `factor()` still orders numeric codes numerically.

After the fix, on the same data: `gr = g1:g2` and `gr = gf1:gf2` both
give -30.88513, the control's value, and `g1/g2`, `g1 + g2`,
`g1 * g2`, `factor(g)` and `interaction(g1, g2)` are all refused.

**What this breaks.** `gr = factor(g)` and `gr = interaction(g1, g2)`
worked on the base build and gave the right answer (-30.88513 each).
They are refused now, because brms refuses them and because accepting
one call-valued `gr` while refusing another is the grammar that let
`g1/g2` through. The refusal names the replacement.

## 3. D3: a hypothesis with no relation

`hypothesis(fit, "Age")` returned exactly the row of `"Age = 0"`. brms
refuses: "Every hypothesis must be of the form 'left (= OR < OR >)
right'". `brms:::eval_hypothesis()` requires exactly one sign and
exactly two sides.

**Fixed** in `hyp_parse()` (`R/confint.R`). The message contains brms's
own sentence and then names the fix.

**This was a documented frmtmb extension.** The `hypothesis()` page
said "or a bare `"expr"`, which brms does not accept and which is
tested against 0 here". That sentence is gone. No `dev/` document
decided it, and `dev/reviews/20260916-brmsnames.md` line 69 is about a
different "bare spelling" (a parameter NAME written without `b_`).

**What it cost.** `dev/adefects-hypscan.R` parses every `.R` and `.Rmd`
in the repository and walks the syntax tree, so a call split over four
lines counts once. Before the change: 94 `hypothesis()` calls with a
literal hypothesis that states a relation, **48 that do not**, and 15
whose hypothesis is not a literal. All 48 were rewritten, plus three
variables holding a bare string (`icc` in `test-boot.R`, `h` and
`hyps` in `test-multiple-pooling.R`) and three that build one from
`variables()` (`test-review-v29.R`, `test-v19.R`). After: 145 with a
relation and 5 without, of which 4 are the refusals this lane's own
`test-adefects.R` constructs and the fifth is brms's `"b_Age x 0"`
inside the generated tier file. One test expectation changed with it:
the
bootstrap draws are keyed by the hypothesis AS WRITTEN, so
`colnames(attr(h, "draws"))` is now `"exp(x) = 0"` rather than
`"exp(x)"`.

The refusal message that pointed AT the bare spelling was corrected
too: `R/confint.R` told a user to write
`hypothesis(fit, "<name>", method = 'profile')`, and now writes
`"<name> = 0"`.

## 4. D4: `fit$data` was `fit$data2`

A `frmtmb_fit` had no `data` element, and there is no `$.frmtmb_fit`,
so `fit$data` partial-matched `fit$data2`: an empty list on most fits,
and `list(A = <6 x 6 matrix>)` on a fit with `data2 = list(A = A)`.
`fit[["data", exact = TRUE]]` was `NULL`.

**What the element holds, and how it differs from brms's.** It is the
MODEL FRAME, the object `model.frame()` returns: one column per term
the formula names, under the term's own spelling. brms keeps its
VALIDATED RAW DATA there. On `count ~ Trt + Age + offset(Age)` brms's
`$data` is 40 x 7 (`count, Trt, Age, volume, visit, Exp, patient`) and
frmtmb's is 40 x 4 (`count, Trt, Age, offset(Age)`), with a literal
`offset(Age)` column and no column a term does not use. So
`names(fit$data)`, `ncol(fit$data)` and `newdata = fit$data` differ from
brms on any model with a transformed term, and `newdata = fit$data` is
not the brms idiom it looks like. **The first version of this record
said brms keeps a model frame there; it does not.** The model frame is
still the right object to store, because it is the one the package
already had, but the parity claim was wrong and is withdrawn. The
attributes differ too: frmtmb's carries `names`, `row.names`, `class`
and `terms`, brms's carries `drop_unused_levels` and `data_name` as
well, and the missing `data_name` is why `brmsfit-methods:924` is still
a defect (section 10).

**The decision: a real `data` element, not a field served by a `$`
method.** Punch round 1 re-opened this against a measurement
(`dev/adefects-p1-dollarcost.R`, `dev/adefects-log/p1-dollar.txt`),
because the first version of this record claimed the element was free
and it is not.

**The claim that was wrong.** In MEMORY the element is shared: one
`.Internal(inspect())` address for `fit$data` and
`fit$frame[["data_frame"]]`, where an equal copy has a different one.
ON DISK it is not, because R's serializer does not deduplicate an
ordinary shared value. The control, which needs no address reading:

| object | raw bytes |
|---|---|
| `list(mf)` | 320,194 |
| `list(mf, mf)`, the SAME frame twice | 640,314 |
| `list(mf, <an equal copy>)` | 640,314 |

`mf` is 20,000 x 2. Sharing buys nothing at `saveRDS()`. On a real fit,
`y ~ x + (1 | g)`, seed 20260917:

| n | raw with `data` | raw without | delta | gzipped delta |
|---|---|---|---|---|
| 240 | 3,690,464 | 3,684,677 | +5,787 (+0.16%) | +4,825 (+1.49%) |
| 20,000 | 8,029,502 | 7,624,563 | +404,939 (+5.31%) | +313,836 (+14.00%) |

One copy of that 20,000-row frame is 405,431 raw bytes, against a
measured delta of 404,939, so the delta is 99.88 percent of one whole
frame rather than the frame's own figure. The review measured +400,971 and +309,917 on its own
construction; same magnitude, same conclusion.

**What the alternative would cost, measured rather than assumed.** The
first version rejected a `$.frmtmb_fit` on the grounds that `$` is the
commonest read in the package's own code. That is the claim the
coordinator asked to measure, and it is wrong in magnitude. A counting
`$.frmtmb_fit`, registered and verified to be reached from INSIDE
frmtmb's namespace (10 reads during one `predict()`), gives:

| user action | `$` reads on the fit |
|---|---|
| `frm(y ~ x + (1 \| g))` | 9 |
| `predict(fit)` | 10 |
| `predict(fit, newdata, se.fit = TRUE)` | 32 |
| `fitted(fit)` | 12 |
| `simulate(fit, nsim = 5)` | 49 |
| `summary(fit)` | 39 |
| `ranef(fit)` | 4 |
| `vcov(fit)` | 5 |
| `model.frame(fit)` | 1 |
| `emmeans::emmeans(fit, "x")` | 18 |
| `frm_bootstrap(fit, nsim = 20)` | **751** |

One read costs 0.102 us plain and 0.949 us through a method, so
+0.847 us; the CONTROL, a second plain arm built from the same code,
reads 1.0035. At 751 reads, the largest action measured, the method
would add **0.636 ms** to a bootstrap that takes seconds. That is not a
reason to choose either way.

**So the choice is made on behavior, not cost, and brms is the
tiebreaker.** brms's `data` is a REAL element: `names(brmsfit)` contains
it, `fit[["data"]]` returns it, `str()` shows it. A `$` method would
make it a phantom, the shape the family object has and its own page
warns about ("The link fields are not elements of the list, so
`names()`, `[[` and `str()` do not show them"). The residual cost is
0.16% of a 240-row fit and 5.3% of a 20,000-row one, on an object that
is already megabytes because it carries the AD tape, and it is stated
in `?frm`, in `NEWS.md` and here rather than called free.

**The pin, and the inverse case.** Four assertions in
`tests/testthat/test-adefects.R` fail on the alternative design, which
`dev/adefects-p1-d4inverse.R` builds and runs: `"data" %in% names(fa)`,
`!is.null(fa[["data", exact = TRUE]])`, no `$.frmtmb_fit` registered,
and the serialized delta exceeding 0.9 of the frame's columns. On the
alternative all four read FALSE, the delta being exactly 0
(`dev/adefects-log/p1-d4inverse.txt`). The memory sharing is pinned
separately by an ADDRESS comparison with an equal copy as its control,
because the assertion the first version shipped,
`expect_identical(fa$data, fa$frame[["data_frame"]])`, is value equality
and passes on a copy.

The disk assertion is a ratio to the frame's COLUMNS, not to the frame:
a model frame's `terms` attribute carries an environment, and
serializing it from inside a test drags the calling frame in, 3,696,708
bytes against 911 for the columns alone. That is a trap for anyone
re-measuring this.

**What it leaves open.** `data` was the one brms field name that
collided. The other brmsfit element names (`formula`, `family`,
`ranef`, `criteria`, `version`, `algorithm`, `backend`, `basis`,
`stanvars`, `model`) all read `NULL` on a frmtmb fit under `$`'s own
partial matching, which the test asserts with the two names that DO
partial-match as its control.

## 5. D5: `frmtmb.latent::hmm_starts(1)`, dissolved

Not reproducible. It was fixed at 0.60.0 by the conditions lane
(commit `6503fc8`), and `extensions/frmtmb.latent/NEWS.md` says so in
as many words: "`hmm_starts()` and the other `hmm_*()` readers refuse
an object that is not a fit by name. `hmm_starts(1)` used to fail with
`$ operator is invalid for atomic vectors`."

On the reference build `rellib-r3` (frmtmb.latent 0.4.0), all five
exported `hmm_*`/`lca_*` readers, each given `1`, `"a"`, `NULL` and
`list()`, raise a classed
`c("frmtmb_latent_error", "frmtmb_error", "error", "condition")`: 20 of
20 calls, 0 unclassed (`evidence.txt`).

**A false alarm of my own, recorded because it will be made again.**
The first lane run DID see `$ operator is invalid for atomic vectors`.
The cause was the library path: `dev/adefects-repro-base.R` put only
the build under test and the user library on it, so with the lane build
first the reference build dropped out and frmtmb.latent came from the
USER library, at 0.3.0. The script now keeps `rellib-r3` behind
whichever build is under test, and prints the version of every package
it reads.

## 6. D6: `predict(newdata)` without the grouping column

Two faults, and brms decides both.

**With `allow_new_levels = TRUE` brms answers.**
`validate_newdata()`: "grouping factors do not need to be specified by
the user if new levels are allowed". It fills the column with `NA`.
Measured on `brmsfit_example1` with `visit` removed from a 3-row
newdata: `posterior_epred()` returns 25 x 3, and
`brms:::validate_newdata()` fills `visit` with `NA NA NA`
(`repro-lane.txt`). frmtmb stopped at base R's "object 'g' not found",
raised from `eval(comp$bar[[3]], newdata, env)`.

**Without it brms stops too**, and on the same base R error: "object
'visit' not found".

**Fixed** in `fill_new_group_vars()` (`R/predict.R`), called from
`pred_design()`, which is the one place every newdata design is built.
With `allow_new_levels` the missing columns are filled with `NA`;
without it the call is refused by a classed error naming the column,
`allow_new_levels = TRUE` and `re_formula = NA`. brms's refusal is a
bare `simpleError`, but item 2.6e requires frmtmb's refusals to be
classed, and both builds refuse, which is the behavior the user's rule
is about.

Only the ordinary grouping blocks are filled: a factor smooth has its
own named refusal in `smooth_newdata_check()` and an `spde()` block's
levels are mesh row numbers, where a new level means nothing. The
variables filled are `all.vars(comp$bar[[3]])`, which is exactly the
set `eval()` needs; `gr(g, cov = A)` records `bar[[3]]` as `g`, so the
`A` of `data2` is not in it (measured).

On a `y ~ x + (1 | g)` fit with 8 groups, seed 20260917, two newdata
rows:

| call | prediction | se.fit |
|---|---|---|
| `newdata` without `g`, `allow_new_levels = TRUE` | -0.4717254, 0.4963306 | 0.834077, 0.868289 |
| `newdata` with a new level, `allow_new_levels = TRUE` | the same | the same |
| `newdata` without `g`, `re_formula = NA` | the same | 0.3183920, 0.3995182 |

An absent column is an unseen level, an unseen level is the population
value for a maximum-likelihood fit, and the standard error carries the
block's marginal variance where the population one does not. brms's
draws differ between the two because it SAMPLES a new level's effect
per draw; a single parameter vector has nothing to sample.

## 7. D7: `log_lik()` on a fit

`log_lik(fit)` was "could not find function": frmtmb defined no generic
of the name, and `frmtmb.sample` defined its own. `loo()`, `waic()`,
`LOO()`, `WAIC()` and `bayes_R2()` are all in frmtmb with a
`frmtmb_fit` method that refuses and names the route to draws;
`log_lik()` was the one missing.

**Fixed by moving the generic to frmtmb**, which is the rule
`extensions/frmtmb.sample/R/reexports.R` already states: "frmtmb keeps
every generic it still has a method for, and this package registers its
`frmtmb_draws` methods on those rather than defining a second generic
of the same name."

- `R/loo.R` defines `log_lik()` and `log_lik.frmtmb_fit`, registered on
  frmtmb's own generic and on `rstantools::log_lik`.
- `R/generic-owners.R` adds `log_lik = "rstantools"`, so frmtmb's
  binding resolves to rstantools' generic whenever it is loaded.
- `frmtmb.sample` re-exports frmtmb's binding, drops `log_lik` from its
  own owner table, and keeps `log_lik.frmtmb_draws` on the same
  generic. Its doc page for the method is `sample-log_lik`.

Why the lighter route was not taken. Defining the generic in frmtmb
while `frmtmb.sample` kept its own would work only while rstantools is
loaded: without it, `library(frmtmb.sample)` masks frmtmb's binding
with its own fallback generic, whose method table has `frmtmb_draws`
and not `frmtmb_fit`, and `log_lik(fit)` is back to "no applicable
method". rstantools is a Suggests of frmtmb.sample, so that session is
a real one.

`tests/testthat/test-generic-collision.R` gains `log_lik` to its
`brms_shared` list (brms has a `log_lik.brmsfit`), and
`extensions/frmtmb.sample/tests/testthat/test-generic-collision.R`
loses it from `own_generics` and its owner table. Two probes in that
file used `log_lik` as their example of an active binding; they now use
`posterior_epred`, which the package still defines, because the
importer they build imports only `own_generics`.

**What the move costs, found in review.** `S3method(log_lik,
frmtmb_draws)` now resolves through frmtmb's ACTIVE binding at
registration time, so the draws method lands in whichever table that
binding pointed at when frmtmb.sample loaded: rstantools' when
rstantools was loaded, frmtmb's when it was not. So an UNLOAD of
rstantools after that point loses the method:
`library(rstantools); library(frmtmb); library(frmtmb.sample);
unloadNamespace("rstantools")` reaches `log_lik.frmtmb_draws` on the
base build and reports "no applicable method" on this one, and the same
pair starting from `library(brms)` and unloading brms and rstantools
(`dev/adefects-rev-d7.sh`, 26 scenarios per arm; of those, 19 load
frmtmb.sample and the draws method dispatches in 17 of them, and 11
have brms loaded when the probe runs, where brms dispatches in 11 of
11. Counted from the log rather than from the scenario list, which is
where the review's own first figures, 22 of 24 and 12, came from). This
is not
a new CLASS of problem: `dev/adefects-rev-d7-preexisting.R` shows 6 of
14 probed generics, `loo`, `waic`, `bayes_R2`, `as_draws_df`,
`variables` and `ndraws`, already behave that way on the base build,
because they are already re-exported from frmtmb. The move makes
`log_lik` consistent with the majority rather than opening a new hole,
the failure is loud, and it needs an explicit `unloadNamespace()` to
reach. Recorded rather than fixed: fixing it means registering every
re-exported method in both tables, which is a change to
`frm_install_generics()` and to all seven of these names at once.

## 8. False alarms: what the three new refusals cost

Two arms of one harvest, `dev/adefects-falsealarm.R`. Each parses every
`.R`, `.Rmd` and `.Rd` in brms 2.23.0's own tests, vignettes and manual
and in this repository, walks the syntax tree, and asks brms's own
constructor (`brms::ar()` and friends, which validate and do not
evaluate) and frmtmb's parser what each harvested call gets. The base
arm reads the reference build and a `git archive` extract of `f8b45ef`;
the lane arm reads this build and this tree. Diffing the two TSVs is
what attributes a verdict to THIS lane rather than to a divergence the
package already had, such as the documented `cov = FALSE` refusal,
which swamps the raw count (23 of the 57 autocor calls, in BOTH arms).

**Autocor grammar, D1 and D2.** 137 calls harvested, 57 distinct.

| | base arm | lane arm |
|---|---|---|
| both accept | 29 | 29 |
| both refuse | 3 | 3 |
| frmtmb refuses, brms accepts | 23 | 23 |
| frmtmb accepts, brms refuses | 2 | 2 |

**Not one row changed verdict.** The two arms harvest the same 137
calls and the same 57 distinct ones, and every verdict matches. The new
refusals fire on nothing that brms's own suite or this repository
writes: nobody writes an expression as a time index or an operator in
`gr`. The false-alarm rate attributable to this lane is 0 of 57.

`tests/testthat/test-adefects.R` is excluded from the harvest in both
arms. It is built to make each new refusal fire, so counting it would
report this lane's own constructions as its false alarms; with it in,
the lane corpus carried eight autocor calls the base tree does not
have.

The 23 disagreements are all pre-existing and all the same two
recorded divergences: `ma(x)`, `arma(p = 2, q = 2, gr = g)`, `ar(time)`
and 18 more are `cov = FALSE`, which `R/autocor.R:102` and
`dev/feature-gaps.md:228` declare unimplemented; `unstr(week)` needs
both arguments. The 2 misses are `ar(week, subj, p = 2, cov = TRUE)`
and `ar(t, g, p = 2, cov = TRUE)`, where brms limits `cov = TRUE` to
order one and frmtmb deliberately does not.

This measurement says the refusals are cheap, and it does NOT say they
are reachable: no harvested call reaches them. Reachability is proved
by construction in sections 1 and 2 and pinned by
`tests/testthat/test-adefects.R`.

**Hypothesis grammar, D3.** Harvested strings, judged by
`brms::hypothesis()` on `brmsfit_example1`, with brms's grammar refusal
told from its other errors by its text (a refusal of an unknown
PARAMETER counts as grammar-accepted, which is what is being measured).

| | base arm | lane arm |
|---|---|---|
| distinct strings | 98 | 92 |
| both accept | 70 | 88 |
| both refuse | 3 | 4 |
| frmtmb refuses, brms accepts | 0 | 0 |
| frmtmb accepts, brms refuses | **25** | **0** |

The refusal closes 25 distinct strings that brms refuses and frmtmb
answered, and produces no false alarm against brms on 92. The 24 that
exist only in the base corpus are the bare ones this lane rewrote
(`"x"`, `"sigma"`, `"sd_g__Intercept"`, `"exp(x)"`, the ICC
expression, and so on); the twenty-fifth is `"b_Age x 0"`, which
changed verdict in place, and is the ONE row `dev/adefects-fadiff.R`
reports as changed. The 18 that exist only in the lane corpus are the
replacements, every one of them accepted by both.

**A false alarm in the instrument, recorded.** The first run of the
diff reported a second changed row, `ar(t, g, cov = TRUE)` going from
accepted to refused, which would have been a serious regression. It was
`read.delim()`'s default quoting: the harvest writes with
`quote = FALSE` and the corpus contains `ma("log")`, whose double
quotes swallowed the following fields and shifted a verdict into the
wrong row. `dev/adefects-fadiff.R` reads with `quote = ""` and asserts
that neither arm has a duplicate key.

## 9. The tests, seen failing

`tests/testthat/test-adefects.R` and
`extensions/frmtmb.sample/tests/testthat/test-adefects.R`. Both were
run against the reference build before the fix:

```
RESULT frmtmb test-adefects.R pass=16 fail=24 err=4 skip=0 blocks=6
  BLOCK fail=12 err=TRUE  an autocorrelation term refuses an expression as its time index
  BLOCK fail=0  err=TRUE  gr = takes variable names crossed by `:` and nothing else
  BLOCK fail=4  err=FALSE a hypothesis with no relation is refused, as in brms
  BLOCK fail=8  err=FALSE fit$data is the model frame, not a partial match of data2
  BLOCK fail=0  err=TRUE  newdata may omit a grouping column when new levels are allowed
  BLOCK fail=0  err=TRUE  log_lik() on a maximum-likelihood fit says to sample
RESULT frmtmb.sample test-adefects.R pass=0 fail=0 err=2 skip=0 blocks=2
```

(`dev/adefects-log/newtests-base.txt`,
`newtests-base-sample.txt`.) Every block fails. On this build they are
85 and 16 expectations, 0 failures (`suite-core.txt`,
`suite-sample.txt`); punch round 1 added the four D4 pins and the
`allow_new_levels` forwarding block. The review counted 18 passes on
the reference build where this record says 16, on an intermediate
version of the file; the 24 failures, the 4 errors and the six failing
blocks agree exactly.

Where the defect was a wrong NUMBER rather than a missing refusal, the
assertion is on the number: D2 compares the numeric-code fit against
the factor fit and against the control, D6 compares the filled-column
prediction against a genuinely new level and against the population
one, and D4 asserts the identity `fit$data` is `model.frame(fit)` is
`fit$frame[["data_frame"]]`.

## 10. The ported tier and the ledger

`dev/brmsport-ledger.tsv`, regenerated. Totals before and after:

| outcome | round 2 | after this lane |
|---|---|---|
| pass | 192 | 200 |
| defect | 112 | 103 |
| divergence | 34 | 34 |
| pending 2.6d | 12 | 13 |
| cannot transfer | 144 | 144 |

**Rows this lane moved** (`dev/adefects-verdicts.R` applies them
idempotently; `wt-shapes` is regenerating the same ledger in another
worktree, so this is the list to reconcile against):

| id | was | is | why |
|---|---|---|---|
| `brm:106` | defect, silent | pass, own-words | the grammar refusal now fires on brms's own call |
| `brm:108` | defect, silent | pass, own-words | the same |
| `data-helpers:9` | defect, fit-data | pass, own-words | frmtmb refuses the new `visit` level in its own words |
| `brmsfit-methods:179` | defect, fit-data | pass | `fit1$data` is the 40-row frame |
| `brmsfit-methods:417` | defect, misparse | pass | frmtmb's refusal CONTAINS brms's sentence, so brms's own pattern holds |
| `brmsfit-methods:567` | defect, fit-data | pass | `model.frame(fit1) == fit1$data` is now an identity |
| `brmsfit-methods:675` | defect, fit-data | pass | `pp_check(fit1, newdata = fit1$data[1:10, ])` |
| `brmsfit-methods:1035` | defect, fit-data | pass | a real column and a real attribute read, where it used to hold on a `NULL` |
| `brmsfit-methods:345` | defect, fit-data | defect, shape | `fitted()` returns a vector: rule 3, `wt-shapes` |
| `brmsfit-methods:350` | defect, fit-data | defect, shape | the same, ordinal |
| `brmsfit-methods:352` | defect, fit-data | defect, shape | the same |
| `brmsfit-methods:682` | defect, fit-data | defect, internal-error | `pp_check(group = )` is broken, section 11 |
| `brmsfit-methods:764` | defect, fit-data | pending 2.6d, shape | `predict()` returns a vector: item 2.6d |
| `brmsfit-methods:924` | defect, fit-data | defect, output | brms records `attr(data, "data_name")` and frmtmb does not |
| `data-helpers:7` | defect, fit-data | defect, output | R's "contrasts dropped" warning, section 11 |
| `data-helpers:12` | defect, fit-data | defect, different-error | frmtmb refuses `visit`, brms refuses the unused `fac` |

Eight rows became passes; eight kept their verdict and got the class
and the reason the run now reports, because every one of the 13
`fit-data` reasons named a partial match that no longer happens. The
class `fit-data` is gone from the ledger.

The harness itself is unchanged. `dev/brmsport-guards.R` reports 60 of
60 with the case where each guarded thing is absent
(`port-guards.txt`), the own-words specificity check passes all 30
patterns against the other 29 rows' messages with 0 cross-matches, and
the tier runs 15 of 15 files, 631 expectations, 0 failures, errors or
skips (`port-tier.txt`).

**Four harness scripts took a path as an argument.**
`dev/brmsport-run.R`, `dev/brmsport-record.sh`, `dev/brmsport-tier.sh`
and `dev/brmsport-guards.R` named the worktree `frmtmb-wt-brmsport` and
the library `brmsport-lib`, neither of which outlived that lane, and
`dev/brmsport-run.R` asserted the StanHeaders 2.32.10 pin, which was
removed on 2026-09-17. Without this the tier could not be run at all.
`FRMTMB_PORT_LIB` defaults to the round's reference build;
`FRMTMB_PORT_ROOT` has NO default and the two shell scripts refuse to
run without it. A default was the first spelling, and the obvious
default is the main checkout, which `dev/organizer-rules.md` forbids a
lane agent from writing to: an unset variable would have put this
lane's logs in main. **`wt-shapes` is likely to need the same change**,
and its copy will conflict with this one.

## 11. Found and NOT fixed

**Item 7 is the one to act on first.** A partial `re_formula` is
accepted, not honored, and nothing is said; the review found it and it
is the only SILENT one in this list. The rest are loud.

1. **`pp_check(fit, type = "<any>_grouped", group = )` is broken for
   every grouped type.** `pp_check.frmtmb_fit()` passes `...` straight
   to bayesplot's `ppc_*`, so `group = "g"` arrives as a length-one
   CHARACTER STRING where bayesplot wants one value per observation;
   brms resolves the name against the data first. Measured on a plain
   `y ~ x + (1 | g)` fit with 40 rows: `violin_grouped`,
   `stat_grouped` and `dens_overlay_grouped` all die on bayesplot's
   "length(group) must be equal to the number of observations", with
   and without `newdata`, and the same on brms's fixture 1
   (`/tmp` probe, reproduced in `port-record.txt` through
   `brmsfit-methods:682`). This was MASKED before: the row that
   reaches it read `fit1$data`, which used to be the `data2` list, so
   the call failed earlier. It is a LOUD failure of a brms argument,
   not a silent wrong answer, so it is out of this lane's scope; it
   belongs with `brmsfit-methods:694`, `:695`, `:698` and `:699`,
   which are already `argument` and `internal-error` defects.
2. **`predict(newdata = )` leaks R's "contrasts dropped from factor
   Trt" warning**, twice, from `model.matrix()` on a newdata slice;
   brms is silent (`data-helpers:7`). The NUMBERS are not affected:
   the five predictions equal the first five of the full-data ones,
   maximum absolute difference 0 at full precision
   (`evidence.txt`). Noise, not a wrong answer, and not on the silent
   list.
3. **`update(fit, newdata = )` does not record
   `attr(fit$data, "data_name")`**, which brms sets to the name of the
   newdata argument (`brmsfit-methods:924`). Now visible because the
   frame is there; a one-line addition to `update()`, which this lane
   did not make because `update()` is not on its list.
4. **`fitted()` has no `allow_new_levels` argument**, so the D6 fix is
   not reachable through it (`repro-lane.txt`). `fitted()`'s argument
   surface is `wt-shapes`'s. `posterior_epred()` and
   `posterior_predict()` on draws had the same shape of fault and were
   FIXED in punch round 1: they used to swallow the argument in `...`
   and then meet the D6 refusal, whose remedy named it (section 14).
5. **`hypothesis()` output still carries the class name
   `brmshypothesis`** (`R/confint.R:2673`), which contradicts the
   user's rule 2 of 2026-09-17. Filed by the brmsport lane, still
   open, and not on this lane's list.
6. **`frmtmb.latent::hmm(time = )` evaluates an expression the same
   way D1 did** (`extensions/frmtmb.latent/R/hmm.R:795`), and so does
   `car()`/`spde()`'s `gr =` (`R/frame.R:2331`, parsed at
   `R/parse.R:433`). Left alone deliberately: brms has no `hmm()`
   family, so there is no judge for it, `hmm()`'s `time` orders a
   sequence rather than setting a lag (a monotone expression gives the
   same order), and `frmtmb_structure(frame_vars = )` is built to carry
   a compound expression on purpose. Recorded because the next reader
   of D1 will ask.

7. **FIXED in lane wt-reunc** (`dev/reunc-findings.md`): a one-sided
   `re_formula` keeps the terms it names now, as brms's does, and the
   D6 divergence below is gone with it. As filed:
   **A partial `re_formula` is silently ignored**, which is the most
   valuable thing the review found and the one this lane would file
   first. On a `y ~ x + (1 | g) + (1 | h)` fit,
   `predict(fit, newdata, re_formula = ~ (1 | h))` returns a result
   `identical()` to `re_formula = NULL`, the FULL conditional
   prediction, where brms keeps only the `h` block: brms computes the
   required grouping variables from
   `update_re_terms(formula, re_formula)`, so a one-sided formula
   selects blocks there. frmtmb's `re_form_keeps()` (`R/predict.R:63`)
   is a two-way switch, `NULL` against everything else, and
   `?predict.frmtmb_fit` documents only `NULL`, `NA` and `~0`. It is
   the same family as item 2.5e's `fitted(re_formula = NA)` swallow: an
   argument accepted, not honored, and nothing said.
   `dev/adefects-rev-reformula.R` has the construction.

   **It also explains a divergence the D6 fix inherits.** Under
   `re_formula = ~ (1 | h)` the new refusal demands the `g` column,
   because `use_re` is TRUE and every block is kept. That is consistent
   with frmtmb's own semantics and divergent from brms's, where the `g`
   block is dropped and its column is not required. Fixing
   `re_formula` would remove the divergence; fixing D6 alone cannot.
8. **`pp_check(type = "violin")` fails with "object 'ppc_violin' not
   found"** (`dev/adefects-rev-ppcheck.R`). Ungrouped, so it is a
   different defect from item 1: bayesplot has no `ppc_violin`, and
   `pp_check.frmtmb_fit()` builds the function name by pasting
   `"ppc_"` onto `type` without checking that the result exists. Loud,
   like item 1, and the same family as `brmsfit-methods:697`, already a
   defect of class `internal-error`.

## 12. What needs the user

Nothing blocking. Three decisions are recorded here rather than taken:

1. **The bare `hypothesis()` spelling was a documented feature**, and
   removing it is a breaking change to a documented API that 48 call
   sites in this repository used. It was removed because brms is the
   tiebreaker and the brief asked for it. If the user wants it back,
   the one-line inverse is in `hyp_parse()`.
2. **`gr = factor(g)`, `gr = interaction(a, b)` and
   `ar(as.integer(t), g, cov = TRUE)` were correct calls** that are now
   refused, for brms parity. All three fit at -30.88513 on the base
   build, the same value as the bare-name spelling, and brms refuses
   all three. Same shape of decision as 1.
3. **The version floor.** `frmtmb.sample` now needs the frmtmb release
   that exports `log_lik()`, and both NEWS sections are under
   "(development version)" with no number chosen.

## 13. Runs

Every test file in its own R process, on the lane library with the
round's reference build behind it. `dev/adefects-run.R` is a copy of
`dev/release/run-tests.R` with the library swapped; it asserts
`tmbstan >= 1.2.1` and the user Makevars `-std=gnu++17` flag, which is
what replaced the StanHeaders pin, and a whole-suite run in one process
is not used because it has repeatedly hidden state leakage here.

<!-- BEGIN GENERATED: dev/adefects-runs.R -->

| run | files | pass | fail | error | skip |
|---|---|---|---|---|---|
| core suite | 153 | 10302 | 0 | 0 | 133 |
| frmtmb.coupling | 9 | 447 | 0 | 0 | 5 |
| frmtmb.eam | 26 | 1652 | 0 | 0 | 3 |
| frmtmb.latent | 9 | 349 | 0 | 0 | 2 |
| frmtmb.learn | 14 | 408 | 0 | 0 | 13 |
| frmtmb.ode | 10 | 343 | 0 | 0 | 1 |
| frmtmb.sample | 29 | 1727 | 0 | 0 | 3 |
| frmtmb.spline | 14 | 450 | 0 | 0 | 1 |
| gated tier | 23 | 2496 | 0 | 0 | 0 |
| ported brms tier | 15 | 631 | 0 | 0 | 0 |
| **total** | **302** | **18805** | **0** | **0** | **161** |

<!-- END GENERATED -->

**Two ways a run reported nothing and looked clean, both caught by the
summariser rather than by reading the log.**

- **`NOT_CRAN` unset.** The second full core pass reported
  `pass=0 fail=0 err=0 skip=0` for 43 files, which reads exactly like a
  file with nothing in it, because the launcher did not export
  `NOT_CRAN` and every `skip_on_cran()` block was skipped.
  `dev/adefects-suite.sh` now exports it itself, and
  `dev/adefects-run.R` REFUSES to run without it, which is the guard
  failing closed. (The shipped drivers in `dev/release/` all set it
  already; this is a lane-local runner.) `dev/adefects-runs.R` stops on
  any file with zero expectations and zero skips, which is how the 43
  were found.
- **Another R process in the same minutes.** An earlier pass had two
  files produce no RESULT line at all and two produce zero, because the
  false-alarm harvest was running beside the suite.
  `dev/adefects-suite.sh` prints the file name when no RESULT line
  appears, and `dev/adefects-runs.R` stops on it; a runner that sums
  `failed` and not `error` would have printed a clean line for each.

**Do not edit a runner while `Rscript` or `sh` is reading it.** Editing
`dev/adefects-run.R` mid-run to add the `NOT_CRAN` guard also corrupted
the read of `dev/adefects-suite.sh` in the loop that was executing it
(`dev/adefects-suite.sh: line 20: tly: command not found`), which
`dev/lane-rules.md` records for `Rscript` and which is true of `sh`
too. The affected suites were rerun from scratch.

`R CMD check --as-cran`, once per package, on the final source
(`dev/adefects-check.sh`; `--no-manual` is NOT passed, because it skips
the manual sections where an unescaped `%` in Rd surfaces):

| package | status | log |
|---|---|---|
| frmtmb | 1 NOTE | `check-core.txt` |
| frmtmb.sample | OK | `check-sample.txt` |

The NOTE is the pre-existing one on this box: "checking HTML version of
manual ... Skipping checking math rendering: package 'V8' unavailable".
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE` removes the pre-existing
CRAN-incoming WARNING about frmtmb not being on CRAN.

`dev/adefects-log/` holds one log per run. The two partial passes that
preceded them were discarded rather than reported, for the reasons just
above. Punch round 1 reran the core, `frmtmb.sample` and the six other
extension suites, the gated tier, the ported tier, the guards and both
checks on the punched build; the numbers in the table and the two check
statuses are from that pass.

## 14. Punch round 1

`dev/reviews/20260918-adefects.md`, verdict "mergeable after punch", no
blocker. Its scripts are `dev/adefects-rev-*` and its logs
`dev/adefects-rev-log/`. What each item became:

| finding | what was done |
|---|---|
| MAJOR 1, "it costs nothing" is false on disk | FIXED, and the choice re-made against a measurement. The duplicate is KEPT and every claim corrected with the on-disk number, in section 4, `NEWS.md` and the `?frm` roxygen. The alternative was measured rather than assumed: a `$.frmtmb_fit` would add 0.847 us per read and 0.636 ms to the largest action (a 751-read bootstrap), so cost does not decide it; brms's `data` is a real element and a `$` method would make it a phantom, so behavior does. Four assertions now pin the choice and all four fail on the alternative (`dev/adefects-p1-d4inverse.R`) |
| MINOR 1, `fit$data` is the model frame, not brms's raw data | FIXED. Section 4 and `?frm` say what differs, with the `offset(Age)` example; a test pins that the frame carries `offset(z)` and not `z` |
| MINOR 2, the registration target of the draws method | RECORDED in section 7, with the two scenarios and the six sibling generics that already behave that way |
| MINOR 3, two draws methods swallow `allow_new_levels` | FIXED by FORWARDING it, which is what brms does (measured: `posterior_epred()` and `posterior_predict()` answer 25 x 3 on `brmsfit_example1` with `visit` removed and `allow_new_levels = TRUE`, and die on base R's "object 'visit' not found" without it). Every other name in `...` is now refused by name, item 2.5e's rule, which the rest of that file already followed |
| NIT 1, a third now-refused correct call | FIXED. `ar(as.integer(t), g, cov = TRUE)` is named in section 12 beside the other two |
| NIT 2, `gr = a:b` pastes and is not injective | FIXED. One sentence in the `gr` documentation, with the colliding pair and the note that brms merges them too |
| NIT 3, a backticked name is accepted where brms refuses | not changed: frmtmb is the MORE permissive side, so it is a divergence and not a false alarm |
| NIT 4, the two shell scripts default to the main checkout | FIXED. `FRMTMB_PORT_ROOT` has no default and both refuse without it; seen refusing |
| NIT 5, no blank line before `# frmtmb 0.60.0` | FIXED, in both `NEWS.md` files |
| NIT 6, stale text beside corrected numbers | FIXED. `brm:81`'s reason now quotes the named refusal the run records, and the withheld-runs sentence in `ledger-summary.md` is GENERATED from the messages (vacuous against stale) rather than typed, which is why it went stale |
| NIT 7, `rstantools` undeclared | FIXED. Added to frmtmb's Suggests; `NAMESPACE` carries three `S3method(rstantools::...)` directives and `frm_generic_owners` names it as an owner |
| out of scope, a partial `re_formula` | FILED, section 11 item 7, with the reviewer's construction and the note that it explains a divergence the D6 refusal inherits |
| out of scope, `pp_check(type = "violin")` | FILED, section 11 item 8 |

**Two numbers the review re-derived differently, and why.** Its
hypothesis harvest gives 82 distinct strings on the lane tree and 88 on
the base against this record's 92 and 98, because it extracts from a
narrower set of call shapes; the direction, the 24 base-only strings
and the 0 false alarms agree. Its "seen failing" count for
`test-adefects.R` on the reference build is 18 passes where this record
says 16, on an intermediate version of the file; the fail and error
counts, 24 and 4, and the six failing blocks agree exactly.
