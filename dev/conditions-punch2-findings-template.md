
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
@PROBE@
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
@OVERREACH@
```
<!-- END GENERATED: punch2 overreach -->

The census test itself now asserts the six legitimate shapes are
clean. All eight real trees pass with no new exemption, and the 26
plants of `dev/conditions-punch-census.R` give the same result as in
round 1:

<!-- BEGIN GENERATED: punch2 census -->
```
@CENSUS@
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
@SUITES@
```
<!-- END GENERATED: punch2 suites -->

No export and no Rd changed, so `R CMD check` was not rerun. The
extensions were roxygenised in a fresh process after core was
installed; their NAMESPACE files did not change.
