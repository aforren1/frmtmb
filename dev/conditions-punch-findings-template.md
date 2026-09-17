
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
@PROBE@
```
<!-- END GENERATED: punch probe -->

"expected subclass" is the extension's subclass for an extension site,
and no subclass for a core site.

### Census plants

`dev/conditions-punch-census.R`:

<!-- BEGIN GENERATED: punch census -->
```
@CENSUS@
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
@SWEEP@
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
@CHECKS@
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
