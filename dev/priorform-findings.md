# Lane wt-priorform: bf(), the formula grammar and the prior surface

Item 2.6c, findings 4, 6, 7, 12, 13 and 14 of `dev/brms-suite-audit.md`
section 8, section 7 contract 7, and the stray singular-convergence
warning in `tests/testthat/test-prior-compat.R`. Base commit aa9227e,
frmtmb 0.58.0; reference build `C:/Users/adf44/source/r/rellib-r3`;
lane build `C:/Users/adf44/source/r/priorform-lib`. brms 2.23.0 is the
judge throughout. Every count below is pasted from the generated block
of the script named beside it.

brms's suite source is not in the repository (`dev/brms-suite/` is
gitignored) and was absent from both checkouts, so it was fetched again
as `dev/brms-suite-audit.md` section 1 says. The tarball hash matched
the recorded one:
`b5f5bb5604ec3f87b3ad99f0da4e6a37d8c5488221e8abe459171340992c37a5`.

## 1. Before and after, per item

`dev/priorform-probe.R ref|lane` runs one construction per item against
either library. Logs: `dev/priorform-probe-ref-log.txt`,
`dev/priorform-probe-lane-log.txt`.

| item | base (ref) | lane |
|---|---|---|
| 4 `identical(form, bf(form))`, dpar and nl | ERROR "`formula` must be a formula" | TRUE, TRUE |
| 6 `set_prior(..., class = c("b", "sd"))` | ERROR "`class` must be a single non-empty string" | 2 specifications, `$class` = `b sd` |
| 7 `default_prior`, `validate_prior`, `empty_prior`, `as.brmsprior` exported | FALSE for all four | TRUE for all four |
| 7 `print(set_prior("normal(0,1)"))` | `normal(0, 1) class=b` | `b ~ normal(0,1)` |
| 12 `(1 | g) + (x | g)` | accepted | refused by name |
| 13 `y ~ x * cs(g)` | accepted | refused, "The term 'x:cs(g)' is invalid" |
| 14 `bf(y ~ ~ x)` | accepted | refused, names `~~` |

What the base commit DID with the three accepted formulas, rather than
whether it errored, is `dev/priorform-base-behaviour.R` (seed 20260916,
n = 240; logs `dev/priorform-base-behaviour-ref-log.txt` and
`-lane-log.txt`):

- `y ~ ~x` fitted exactly `y ~ x`: fixed effects 0.9497471 and
  0.4817102 and log likelihood -394.7416 in both. Lenient rather than
  wrong, and refused because brms refuses it.
- `y ~ x + (1 | g) + (x | g)` fitted two intercept blocks with standard
  deviations 0.56151 and 0.56093, whose squares sum to the 0.79368^2 of
  `(x | g)` alone, at the same log likelihood -348.814. The split is
  not identified.
- `yo ~ x * cs(c3)` under `sratio()` built a design with the single
  column `x` and NO category-specific term: `cs(c3)` and `x:cs(c3)`
  were dropped without a message. `yo ~ x + cs(c3)` carries `csc3`.
  This one was a silent wrong model.

## 2. What changed

Round 1 of review (`dev/reviews/20260916-priorform.md`) and the user's
two policy decisions changed parts of this lane after its first
report. Section 2P records the punch round; the bullets below are
corrected to the code as it now stands.

### bf() (item 4, item 14, and brms's parameter-formula spellings)

- `bf(form)` on a `frmtmb_formula` returns it untouched when nothing
  else is given. With further arguments it adds dpar formulas or
  constants; a parameter the formula already sets is refused with the
  existing "Duplicated dpar formula" message. brms REPLACES it with a
  message instead. This divergence is deliberate: it matches frmtmb's
  `+ lf()`, which already refuses the same repeat, and replacing in
  silence contradicts the refuse-rather-than-accept policy. Nothing in
  brms's suite asserts the replacement.
- `nl` defaults to `NULL`, brms's default, so a built formula keeps its
  flag. Found in passing: `bf(mvbind(y1, y2) ~ ..., nl = TRUE)` did not
  pass `nl` to the per-response `bf()` calls, so it built linear
  models; it passes it now.
- `refuse_nested_formula()` refuses a top-level `~` on the right-hand
  side in `bf()`, its parameter formulas, `lf()` and `nlf()`. brms checks
  the top level only (`as_formula()`), and so does this: `bf(y ~ x + ~z)`
  is accepted by both.
- brms accepts a one-sided parameter formula named by its parameter,
  `bf(y ~ x, sigma = ~ z)`; frmtmb refused it. `bf()` and `lf()` now
  accept it as the same formula, and refuse an unnamed one with brms's
  wording, "Additional formulas must be named". The invalid-name
  message now says "must not contain dots or underscores", which is
  brms's asserted phrase. These two changes are what turned brms
  assertions F1 to F4 from fail to pass.

### Formula grammar (items 12, 13)

- `refuse_duplicated_re()` in `R/frame.R` runs per linear predictor
  after its bar terms are built, on the pair (grouping factor,
  coefficient name) the design actually has, which is brms's key
  (`group`, `coef`, plus resp/dpar/nlpar, which a per-predictor check
  covers). The first measurement run found 6 misses against brms, all
  `(1 | g:h) + (1 | g/h)` designs, because reformulas expands `g/h` to
  `h:g`. The fix first compared an interaction's factors as a set; the
  reviewer showed that refuses `(1 | g:h) + (1 | h:g)`, which brms
  accepts, so the punch round replaced it (section 2P, F3).
- The first full core suite run (`dev/priorform-core-suite-log1.txt`:
  138 files, 8792 pass, 0 fail, 1 error) found the one design in this
  repository the refusal fires on: `test-review-v29.R`'s animal model,
  `(1 | gr(id, cov = A)) + (1 | id)`. That model IS identified, through
  the relationship matrix, so by this project's own standard it is a
  correct model the check fires on. brms refuses it all the same, and
  accepts the same model with the second term on a copy of the column
  (`dev/priorform-grdup.R`: "Duplicated group-level effects are not
  allowed" against "accepted"). brms is the tiebreaker, so the refusal
  stands; its message now names the copy-the-column spelling when
  either block has a known covariance, the test fits the brms spelling,
  and the help text in `R/confint.R` and `R/methods-fit.R` no longer
  describes the one-column form. The generated set did not contain
  `gr()`, which is why section 4 did not see it. The user decided to
  keep this refusal exactly as brms has it.
- The `cs()` refusal, now `special_term_refusal()` in `R/parse.R`
  and covering brms's whole list (section 2P, F2), judges placement on
  the
  term labels `terms()` expands, as brms's `find_terms()` does. The
  first measurement run refused 12 generated formulas brms accepts,
  all `cs(x) OP cs(x)` self-interactions that `terms()` collapses to
  `cs(x)`; refusing them was a false alarm, and dropping them would
  have been the base defect again, so the `cs()` branch now reads the
  expanded labels and takes each bare `cs()` call. The second run
  found 0 false alarms.

### The prior surface (items 6, 7, contract 7)

- `set_prior()` recycles every argument the way `data.frame()` does,
  one specification per row, and refuses lengths that do not divide the
  longest, which brms also refuses. A `frmtmb_prior` object counts as
  one value. Each row keeps the prior STRING as written, because brms
  prints and returns it verbatim.
- `print.frmtmb_priorlist()` follows brms's `print.brmsprior()`: one
  row prints as brms's `.print_prior()` line, several as a table,
  `show_df` chooses. After the punch round an empty prior prints
  nothing and a row with no density of its own shows its class row's as
  `(vectorized)`, both as brms prints them.
- `as.data.frame()` on a prior object gives brms's `set_prior()`
  columns less `tag`: all character, bounds `NA_character_` when unset,
  `source = "user"`. `$` reads a column of that table and `$<-`
  rebuilds the specifications from the edited table.
- **The object shape was NOT changed** to a data frame. A
  `frmtmb_priorlist` is a list of parsed specifications, and
  `frmtmb.sample` reads it as one in `R/sample.R` and in its tests
  (`unclass(pl)`, `s$dist$kind`). brms's `expect_is(bprior,
  "brmsprior")` is therefore a recorded divergence (P3), and `$`
  carries brms's column reads.
- `default_prior(object, data, family, data2, route)` holds the old
  `get_prior()` body; `get_prior(formula, ...)` is its alias, brms's
  signature. It is a plain function, not an S3 generic: brms's
  `default_prior` is a generic, and registering frmtmb methods on it
  would have needed the owner machinery in `R/generic-owners.R` for a
  gain limited to non-plain-formula inputs. Recorded in section 6.
- The table has brms's columns in brms's order plus `source`, and a
  class `sd` or `cor` row for a block of a distributional parameter's
  predictor now carries that `dpar` (`block_dpar()`). Before, two
  `(1 | g)` blocks in `mu` and `phi` printed as ONE row, because
  `unique()` collapsed identical keys. The sampling-route default for
  those rows is still read class-wide.
- `validate_prior()` runs `resolve_priorlist()` on the assembled model,
  so it refuses exactly what `frm()` refuses with the same message
  (asserted with `identical()` on the two messages), then writes the
  prior into the `default_prior()` table: rows it names get
  `source = "user"`, rows under a class row it set get the class row's
  density and `source = "(vectorized)"`, as brms prints them. The rows
  are written in the order the resolver applies them (section 2P, F1).
- `as.brmsprior()` reads a data frame with a `prior` column into a
  `frmtmb_priorlist` with brms's defaults and column dropping. It
  returns frmtmb's object, not a `brmsprior`, for the reason above.
- `empty_prior()` is a zero-length `frmtmb_priorlist`.
- `as_priorlist()`, the `prior =` boundary, now reads a
  `frmtmb_prior_rows` table as written, skipping flat and vectorized
  rows. Construction in `test-brms-formula-priors.R`: a prior and the
  `validate_prior()` table it produced give `identical()` fixed effects
  and log likelihood.

### The theta row (item 7)

brms's test call `count ~ zBase * Trt + (1 | patient) + (1 + Trt |
visit)`, poisson, gets classes b 4, cor 2, Intercept 1, sd 6 in brms
and b 4, cor 2, Intercept 1, sd 3, theta 5 in frmtmb. **Decision: not a
defect.** Evidence:

1. `theta` is a documented `set_prior()` class naming the internal
   covariance parameters, and every theta row names a parameter the fit
   has: `theta_1` to `theta_4` are the patient log sd and the three
   visit parameters.
2. A theta prior moves the fit. `test-brms-formula-priors.R`, "class
   theta rows name real parameters of the fit", seed 6: a
   `normal(-3, 0.05)` prior on `theta_1` cuts the group sd below half
   its unpenalized value.
3. For a `gp()` length-scale and a `car()` dependence parameter, class
   theta is the ONLY spelling: the refusals of brms's `lscale` and `car`
   classes name it (`brms_prior_class_refusal()`), so dropping the rows
   would leave the table unable to list those parameters.

The row count difference against brms is the theta rows plus the
per-coefficient `sd` rows frmtmb does not list, because its class `sd`
addresses a block and refuses `coef` (an existing, documented
refusal). Both are now written in `?default_prior`, section "Where the
table differs from brms's".

## 2P. Punch round 1

Each finding of `dev/reviews/20260916-priorform.md`, fixed or not, with
its construction. The reviewer's scripts are `dev/priorform-rev-*.R`;
the punch round's are `dev/priorform-punch-*.R`. Every guard added here
was run with its behavior put back (`dev/priorform-punch-revert.R`, log
`dev/priorform-punch-revert-log.txt`), and every one fails:

<!-- priorform-punch-revert:begin -->
```
REVERT bfnl brms-formula-priors pass 169 fail 2 error 0
    test-brms-formula-priors.R : bf() on a built formula adds to it and refuses a repeat 
PUNCH-REVERT none brms-formula-priors pass 171 fail 0 error 0
PUNCH-REVERT order brms-formula-priors pass 167 fail 4 error 0
    a validate_prior() table fits the model its prior fits
    a group-specific cor prior wins over the class-wide one
PUNCH-REVERT slots brms-formula-priors pass 163 fail 4 error 1
    two specifications for one slot are refused (brms)
PUNCH-REVERT twins brms-formula-priors pass 159 fail 0 error 1
    exact twin terms are one block, as brms reads them
PUNCH-REVERT slash brms-formula-priors pass 167 fail 2 error 1
    exact twin terms are one block, as brms reads them
    duplicated group-level effects are refused (brms)
PUNCH-REVERT specials brms-formula-priors pass 163 fail 4 error 1
    a smooth or gp() inside an interaction is refused (brms)
    the whole-term specials are brms's own list
```
<!-- priorform-punch-revert:end -->

### Policy A: the duplicated group-level refusal stays as brms has it

Kept, the animal model included, with the copy-the-column message.

### Policy B: duplicate prior specifications are refused. DONE

`check_prior_slots()` (`R/priors.R`) refuses two specifications with
one slot: class as written, coef (with or without the parentheses of
`(Intercept)`), group, resp, dpar and nlpar. It is called where a user
passes a prior: `frm()`, `validate_prior()`, `frm_simulate()`,
`par_template()`, and `frm_sample()` on the call's own `prior =` before
it is stacked. It is NOT in `resolve_priorlist()`, because
frmtmb.sample stacks the fit's prior, the call's prior and the defaults
there and may repeat a slot on purpose. Exported through the sampling
API for frmtmb.sample.

- Refused, as brms refuses: identical rows, two densities, a density
  followed by bounds only. Accepted: class `b` then `coef = "x"`.
  `test-brms-formula-priors.R`, "two specifications for one slot are
  refused (brms)".
- `test-priors-bounds-grcov.R` writes the tightened box as one
  specification and asserts that the two-specification spelling is
  refused.
- The documentation that stated later-wins or tightening
  (`?set_prior` overlap paragraph, Hard bounds, the broadcast
  paragraph, `?validate_prior`, `fill_prior_table()`, `entry_bounds()`,
  `resolve_priorlist()`) is rewritten.
- The different-slot tests stay green: `test-setprior.R:40`,
  `test-prior-compat.R:885`, `test-lkj.R:292` (section 7).

### F1 MAJOR: specificity wins. FIXED

`prior_specificity_order()` sorts the specifications of each class
stably by specificity (coef, then group, then resp/dpar/nlpar) and
keeps each class in its written positions, so `cor` against the `theta`
hatch still follows written order, which `test-lkj.R:292` asserts. The
resolver applies them in that order, and `fill_prior_table()` writes
the table in the same order, so the table shows what the fit applies.

How fits change (`dev/priorform-punch-order.R`, seed 20260916, base
against lane):

<!-- priorform-punch-order:begin -->
| case | base logLik | lane logLik | fit identical |
|---|---|---|---|
| b_class_then_coef | -409.935633 | -409.935633 | TRUE |
| b_coef_then_class | -261.336784 | -409.935633 | FALSE |
| sd_class_then_group | 215.282086 | 215.282086 | TRUE |
| sd_group_then_class | 211.855093 | 215.282086 | FALSE |
<!-- priorform-punch-order:end -->

Only the specific-first order changed, and it now equals the other
order to the last bit.

The reviewer's 17 round-trip cases (`dev/priorform-rev-roundtrip.R`,
rerun on the punch build, log `dev/priorform-punch-roundtrip-log.txt`):
14 give `identical()` fits and resolved entries. The multivariate case
gives identical resolved entries and a log likelihood 2.27e-13 away
(about 1e-15 relative), because the table lists the rows in a different
order than the prior was written, so the prior terms are summed in a
different order; that is why the test on it compares entries, not bits.
The two remaining cases, "bounds then density" and "density then
bounds" on one coefficient, are now refused by policy B, which is the
intended answer. The round-trip test covers both orders for class `b`
against `coef`, and for class-wide `sd` against a group, and a separate
test covers `cor` against a group.

### F2 MAJOR: smooths and gp() inside an interaction. FIXED

`special_term_refusal()` replaces the `cs()`-only check. Its list,
`brms_whole_term_specials`, is brms 2.23.0's `regex_sp()` for every
type `find_terms()` checks with `complete = TRUE`: `cs`/`cse`,
`s`/`t2`/`te`/`ti`, `gp`, the autocorrelation terms and `mmc`. A test
parses brms's own regular expressions and asserts the two sets are
equal, and fails when the list is cut to `cs` alone. A term that
`terms()` collapses to bare specials (`s(z) * s(z)`) becomes them
before dispatch, so every branch sees a special only as a term of its
own. `x:s(z)`, `x * s(z)`, `I(s(z))`, `x:t2(z, x)` and `x:gp(z)` are
refused; `s(z, by = x)` is accepted. The false-alarm measurement is in
section 4, rows `special`.

### F3 MINOR: twins and group order. FIXED

- `drop_twin_bar_terms()` removes a bar term written twice, compared as
  written before `||` or `/` expand, which is what brms's `terms()`
  does. `(1 | g) + (1 | g)` is one block and fits the single-term log
  likelihood exactly (`identical()` in the test); twins in an nlpar,
  in `sigma`, in one response of an `mvbf()`, twin `gr(cov = )` and
  twin `mm()` are accepted.
- `slash_nested_bars()` flags the bar terms that came from `a/b`, and
  the duplicate key reads those factors in brms's order. So
  `(1 | g:h) + (1 | h:g)` is two blocks, as in brms, and
  `(1 | g:h) + (1 | g/h)` is refused, as in brms.
- Kept: the frmtmb-only covariance pairs on one factor and coefficient
  (`rr + diag`, `us + diag`, `equalto + us`, `ar1 + diag`) are refused.
  brms has no verdict, and the reviewer measured them unidentified on
  base apart from the known-matrix design. Stated in the help of
  `refuse_duplicated_re()` and in NEWS.

### F4 MINOR: guard for the mvbind nl fix. FIXED

`test-brms-formula-priors.R` asserts `$nl` on both responses of
`bf(mvbind(y1, y2) ~ a, a ~ 1, nl = TRUE)`. With the old call put back
(`dev/priorform-rev-revert.R bfnl`) it fails 2 assertions; see the
revert block above.

### F5 MINOR: ledger P15 and P27. FIXED IN THE LEDGER; defects FILED

Both reasons now state the sampling-route gap. The two defects are in
frmtmb.sample, which another lane is editing, and are filed in section
6, not fixed.

### F6 MINOR: the repeated-label guard. FIXED

`test-review-v29.R` rebuilds the guard on `(1 | g)` in `mu` and in
`sigma`: `ranef()` keys both by `g`, which is the repeated-key path the
methods carry, and the test asserts both blocks are reported with their
own values. The comment in `R/methods-fit.R` names that design.

### F7 NITs

- The refusal message keeps the term as written, `gr()` and a covariance
  wrapper included: `(1 | gr(g, dist = "student"))`, `diag(1 + x | g)`.
  FIXED, with a test.
- `print(empty_prior())` prints nothing, as brms. FIXED, with a test.
- A row with no density of its own prints its class row's density with
  source `(vectorized)`, as brms prints `as.brmsprior()` tables. FIXED,
  with a test.
- `cs(0 + time || g)` routing: NOT fixed. It is not a one-line change:
  `expand_double_verts()` splits the `||` term into `diag()` pieces
  before the `cs()` branch sees it, so the fix belongs to how a
  covariance wrapper around `||` is read, which is outside this lane.
  Base does the same.

## 2Q. Punch round 2

The "Recheck, round 1" section of `dev/reviews/20260916-priorform.md`,
scripts `dev/priorform-rev2-*`. The guards added here were put back one
at a time (`dev/priorform-punch-revert.R`, log
`dev/priorform-punch2-revert-log.txt`); every one fails. The controls
are the unmutated runs, `none` on the core file and `samplenone` on the
frmtmb.sample file, and both pass:

<!-- priorform-punch2-revert:begin -->
```
PUNCH-REVERT none brms-formula-priors pass 208 fail 0 error 0
PUNCH-REVERT slotlabel brms-formula-priors pass 205 fail 0 error 1
    two specifications for one slot are refused (brms)
PUNCH-REVERT fillup brms-formula-priors pass 205 fail 3 error 0
    every leaf row of a validate_prior() table is what the fit applies
PUNCH-REVERT mmkey brms-formula-priors pass 207 fail 1 error 0
    multi-membership and smooth twins read as brms reads them
PUNCH-REVERT smtwin brms-formula-priors pass 206 fail 2 error 0
    multi-membership and smooth twins read as brms reads them
PUNCH-REVERT samplenone prior-update pass 7 fail 0 error 0
PUNCH-REVERT update prior-update pass 1 fail 3 error 1
    a call's specification replaces the fit's for the same slot
```
<!-- priorform-punch2-revert:end -->

- **R1: FIXED.** `check_prior_slots()` compares the fields
  (`prior_spec_slot()`: class as written, bare coef, group, resp, dpar,
  nlpar) and keeps brms's printed name for the message only. Pin:
  `set_prior("normal(0, 1)", class = "b", dpar = "sigma") +
  set_prior("normal(0, 2)", class = "b", coef = "sigma")` on
  `bf(y ~ sigma, sigma ~ x)`, and the nlpar analogue. On the round-1
  code the pin errored with "'b_sigma' is given 2 times"
  (`dev/priorform-punch2-r1-before-log.txt`); the `slotlabel` mutant
  above restores the label key and fails again.
- **R2: FIXED.** `fill_prior_table()` fills a row the prior does not
  name from the most specific user row that reaches it, looked for
  anywhere in the table, ranked as the resolver ranks. On the
  reviewer's construction (`dev/priorform-rev2-sdphi.R lane`, log
  `dev/priorform-punch2-sdphi-log.txt`) the row sd / g / phi now reads
  `normal(0, 5)` in both written orders, which is what the fit applies.
  The new test "every leaf row of a validate_prior() table is what the
  fit applies" resolves every row that no other row of its class
  narrows and compares the displayed density with the resolver's entry
  at those parameters, over five prior combinations (sd class, dpar and
  group; b class and dpar; coef, class, lkj and sd on a crossed model),
  and asserts it read more than 15 rows. The `fillup` mutant, which
  looks upward only as before, fails it 3 times.
- **R3: FIXED.** A multi-membership block's duplicate key is brms's
  group name, `mm` plus its member variables, so weights and scale do
  not make a second factor. `(1 | mm(g1, g2)) + (1 | mm(g1, g2,
  weights = cbind(wt1, wt2)))` is refused, naming group `mmg1g2`.
- **R4: FIXED.** A bare whole-term special written twice is dropped to
  one term before dispatch, as brms's `terms()` keeps one copy.
  `y ~ x + s(z) + s(z)` builds one smooth, and its log likelihood is
  `identical()` to `y ~ x + s(z)` (test, seed 14).
- **NIT 6: FIXED.** frmtmb.sample's `drop_superseded()` drops a stored
  specification for any slot the call names, bounds-only included, so
  the call's specification replaces the stored one, as
  `brms:::update.brmsfit` does, and the stacked prior passes
  `check_prior_slots()`. The rule is in `?frm_sample`, rendered and
  read back (`dev/priorform-punch2-rd.R`, log
  `dev/priorform-punch2-rd-log.txt`). New test file
  `extensions/frmtmb.sample/tests/testthat/test-prior-update.R`; the
  `update` mutant, which exempts bounds-only specifications as before,
  fails it.

State of the build when the round was closed. The first session of
this round ended before the fill step: the frmtmb.ode suite log was cut
off after two files, and the fill script stops on a log with no counts.
The second session did not trust the record:

- The lane library matches the worktree
  (`dev/priorform-punch2-libmatch.R`, log
  `dev/priorform-punch2-libmatch-log.txt`): frmtmb 979 functions
  identical, 0 differ; frmtmb.sample 182 identical, 0 differ. The same
  script pointed at the base library finds 14 and 2 that differ, so it
  sees a stale install.
- The revert set above was run again on that library, and
  `dev/priorform-rev2-sdphi.R lane` was run again into
  `dev/priorform-punch2-sdphi-log.txt`. Both give the results recorded
  here.

Filed, not fixed (neither is a one-line change with a test, because
each changes which blocks an existing class `sd` prior reaches and so
changes fits and frmtmb.sample's sampling defaults):

- **FILED: a class `sd` prior with a group and no dpar reaches every
  predictor's block.** `set_prior("normal(0, 5)", class = "sd",
  group = "g")` on `bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta()` puts
  `normal(0, 5)` on the phi block too (`dev/priorform-rev2-sdphi.R`,
  "sd group g only": theta1 and theta2 both `normal(0,5)`). brms
  applies a dpar-less `sd` row to `mu` only, and the phi block keeps its
  default. `block_addressed()` reads an empty `dpar` as every predictor.
- **FILED: in a multivariate model, `sd group = g` with no `resp` is
  accepted and broadcast to every response's block.** brms refuses it
  as matching no parameter (reviewer's recheck, carried item).

## 3. The stray warning (test-prior-compat.R:505)

**A test-data problem, not a defect.** The data generator gave `g` a
random intercept and nothing else, while the model asks for
`(x | g) + (z | h)`, so three of four variance components have true
value zero. `dev/priorform-singular.R ref`, log
`dev/priorform-singular-log.txt`:

- seed 77, the test's data: the warning, with the `g` correlation at
  0.999998 and the `h` correlation at 0.999993; `lme4::isSingular()` on
  the same data and model is TRUE.
- 10 seeds of the old generator: every seed that warned had a
  correlation within 1e-3 of plus or minus 1, and lme4 called all 10
  singular.
- the generator with the components the model asks for, seed 77: no
  warning, correlations -0.0589 and -0.4472, `isSingular()` FALSE, and
  the frmtmb and lme4 log likelihoods agree to 1.15e-08.

<!-- priorform-singular:begin -->
old generator: 8 of 10 seeds warn; new generator: 1 of 10
<!-- priorform-singular:end -->

The one new-generator seed that warns (seed 8) also has a correlation
at -0.999997 and `isSingular()` TRUE, so the warning tracks a boundary
fit rather than the generator. The test now uses the new generator and
asserts `expect_no_warning()` on the fit. Nothing was suppressed.
`dev/priorform-warncount.R`, log `dev/priorform-warncount-log.txt`:
the base file on the base build emits 1 warning, the lane file on the
lane build 0.

## 4. False-alarm rates of the three refusals

`dev/priorform-falsealarm.R`, output `dev/priorform-falsealarm.tsv`
(one row per formula and check), log
`dev/priorform-falsealarm-log.txt`. H is every formula literal in
brms's test suite, purled vignettes and Rd examples, found by walking
the parsed AST. G is generated: all singles, all 153 pairs and 300
seeded triples of an 18-term group-level pool, plus 18 designs putting
each pool term in `mu` and `sigma`; and every two-atom arrangement of 13
atoms (plain variables, `cs()`, `s(z)`, `t2(x, z)`, `s(z, by = x)`,
`gp(z)`, `gp(x, by = g)`, `ar(time)`, `te(x, z)`) over 4 operators,
plus 39 wrapped atoms. T is the set the reviewer showed G lacked: every
pool term written twice, reversed and slash-nested interactions, twin
`gr(cov = )`, and twins in an nlpar, in `sigma` and in one response of
an `mvbf()`. The special check applies the refusal predicate to each
term; the duplicate check assembles the frame.

<!-- priorform-falsealarm:begin -->
seed 20260916; brms 2.23.0; 149 source files (118 Rd); 626 distinct harvested formula literals

| check | set | formulas | brms accepts | false alarms | brms refuses (same reason) | misses | other error, excluded |
|---|---|---|---|---|---|---|---|
| nested | H | 626 | 624 | 0 | 1 | 0 | 1 |
| special | H | 129 | 107 | 0 | 3 | 1 | 19 |
| special | G | 715 | 192 | 0 | 495 | 0 | 28 |
| dup | H | 158 | 86 | 0 | 1 | 0 | 71 |
| dup | G | 489 | 235 | 0 | 254 | 0 | 0 |
| dup | T | 30 | 26 | 0 | 4 | 0 | 0 |
<!-- priorform-falsealarm:end -->

Reading it honestly:

- The instrument sees refusals: the "brms refuses" column is refused
  by frmtmb too in every row, and the first run's 12 false alarms and 6
  misses (section 2) were found by it. The first run's "0 false alarms"
  did NOT cover twins or reversed interactions, where the lane then
  refused 5 designs brms accepts; set T covers them now.
- The one special miss is `y ~ (mmc(w1, w2) * y | mm(g1, g2))`, a
  special inside a BAR term, which the term predicate skips on purpose.
  The model is still refused: frame assembly stops with "mmc() must be
  a term of its own on the left of the bar". So it is a miss of the
  predicate, not of frmtmb.
- The special G "other" rows are the 28 formulas with `te()`, which
  brms refuses as not implemented.
- Of the 158 harvested bar formulas, the "other" column is large. By
  verdict pair (brms/frmtmb): 86 accept/accept, 41 accept/other, 4
  other/accept, 26 other/other, 1 refuse/refuse. The 41 are formulas
  brms accepts that frmtmb stops on for another reason before or during
  frame assembly (terms and families it does not support, on data
  synthesized as numeric columns and six-level factors), so those 41
  never tested the duplicate check. The harvested rate is 0 of 86, not
  0 of 127.
- The data are synthesized, so a factor-valued covariate inside a bar
  is numeric in H. G covers factor coefficients with `f`.

## 5. brms's assertions, ported

`dev/priorform-ledger.R lane|ref` runs every assertion of
tests.brmsformula.R (16) and tests.priors.R (37), and the two
tests.brm.R assertions for items 12 and 13, literally against frmtmb,
and STOPS if a recorded verdict disagrees with the run. Output
`dev/priorform-ledger-lane.tsv` and `-ref.tsv`.

<!-- priorform-ledger-ref:begin -->
55 brms assertions run against frmtmb (ref): 4 pass, 51 fail
ids passing on this build: P21 P22 P23 P24 
<!-- priorform-ledger-ref:end -->

<!-- priorform-ledger-lane:begin -->
55 brms assertions run against frmtmb (lane): 26 pass, 29 fail
verdicts: 26 pass, 8 deliberate divergence, 21 cannot transfer

| id | brms file | block | run | verdict | reason |
|---|---|---|---|---|---|
| F1 | tests.brmsformula.R | validates formulas of non-linear parameters | pass | pass |  |
| F2 | tests.brmsformula.R | validates formulas of non-linear parameters | pass | pass |  |
| F3 | tests.brmsformula.R | validates formulas of non-linear parameters | pass | pass |  |
| F4 | tests.brmsformula.R | validates formulas of auxiliary parameters | pass | pass |  |
| F5 | tests.brmsformula.R | detects use if '~~' | pass | pass |  |
| F6 | tests.brmsformula.R | does not change a 'brmsformula' object | pass | pass |  |
| F7 | tests.brmsformula.R | does not change a 'brmsformula' object | pass | pass |  |
| F8 | tests.brmsformula.R | detects auxiliary parameter equations | fail | cannot transfer | frmtmb has no equating of one dpar to another (bf(y ~ x, sigma1 = "sigma2")); the call is refused as an uninterpretable argument, a feature outside this lane |
| F9 | tests.brmsformula.R | detects auxiliary parameter equations | fail | cannot transfer | frmtmb has no equating of one dpar to another (bf(y ~ x, sigma1 = "sigma2")); the call is refused as an uninterpretable argument, a feature outside this lane |
| F10 | tests.brmsformula.R | detects auxiliary parameter equations | fail | cannot transfer | frmtmb has no equating of one dpar to another (bf(y ~ x, sigma1 = "sigma2")); the call is refused as an uninterpretable argument, a feature outside this lane |
| F11 | tests.brmsformula.R | detects auxiliary parameter equations | fail | cannot transfer | frmtmb has no equating of one dpar to another (bf(y ~ x, sigma1 = "sigma2")); the call is refused as an uninterpretable argument, a feature outside this lane |
| F12 | tests.brmsformula.R | detects auxiliary parameter equations | fail | cannot transfer | frmtmb has no equating of one dpar to another (bf(y ~ x, sigma1 = "sigma2")); the call is refused as an uninterpretable argument, a feature outside this lane |
| F13 | tests.brmsformula.R | update_adterms works correctly | fail | cannot transfer | frmtmb exports no update_adterms(), a formula-editing helper outside this lane |
| F14 | tests.brmsformula.R | update_adterms works correctly | fail | cannot transfer | frmtmb exports no update_adterms(), a formula-editing helper outside this lane |
| F15 | tests.brmsformula.R | update_adterms works correctly | fail | cannot transfer | frmtmb exports no update_adterms(), a formula-editing helper outside this lane |
| F16 | tests.brmsformula.R | update_adterms works correctly | fail | cannot transfer | frmtmb exports no update_adterms(), a formula-editing helper outside this lane |
| P1 | tests.priors.R | default_prior finds all classes for which priors can be specified | fail | divergence | frmtmb lists 5 class theta rows (its internal covariance parameters, a real class) and no per-coefficient sd rows (its class sd addresses a block and refuses coef) |
| P2 | tests.priors.R | default_prior finds all classes for which priors can be specified | fail | cannot transfer | frmtmb's sratio() takes no threshold = "equidistant" and has no cse() alias (family surface, lane wt-famlink) |
| P3 | tests.priors.R | set_prior allows arguments to be vectors | fail | divergence | the object is a frmtmb_priorlist, which holds parsed densities; $ reads brms's columns from it |
| P4 | tests.priors.R | set_prior allows arguments to be vectors | pass | pass |  |
| P5 | tests.priors.R | set_prior allows arguments to be vectors | pass | pass |  |
| P6 | tests.priors.R | print for class brmsprior works correctly | pass | pass |  |
| P7 | tests.priors.R | print for class brmsprior works correctly | pass | pass |  |
| P8 | tests.priors.R | print for class brmsprior works correctly | pass | pass |  |
| P9 | tests.priors.R | print for class brmsprior works correctly | fail | cannot transfer | set_prior(check = FALSE) passes Stan code through, and frmtmb builds no Stan program |
| P10 | tests.priors.R | default_prior returns correct nlpar names for random effects pars | pass | pass |  |
| P11 | tests.priors.R | default_prior returns correct fixed effect names for GAMMs | fail | divergence | a smooth's unpenalized column is named s(x).fx1 where brms names it sx_1 (coefficient naming, lane wt-brmsnames) |
| P12 | tests.priors.R | default_prior returns correct fixed effect names for GAMMs | fail | divergence | as P11, and the nonlinear intercept is listed as (Intercept) where brms writes Intercept; set_prior() takes both |
| P13 | tests.priors.R | default_prior returns correct prior names for auxiliary parameters | fail | divergence | frmtmb lists the class-wide sd row once with no dpar and no per-coefficient sd row, so the phi rows are b, b z, Intercept and sd g |
| P14 | tests.priors.R | default_prior returns correct priors for multivariate models | pass | pass |  |
| P15 | tests.priors.R | default_prior returns correct priors for multivariate models | fail | divergence | default_prior() describes frm() unless route = "sample", and frm() is flat in every slot. With frmtmb.sample attached, route = "sample" ALSO gives (flat) here where brms gives lkj(1): a filed frmtmb.sample defect, not fixed in this lane (dev/priorform-findings.md section 6) |
| P16 | tests.priors.R | default_prior returns correct priors for multivariate models | fail | cannot transfer | frmtmb takes no list of families as `family` (a multivariate model attaches one family per bf()), family surface outside this lane |
| P17 | tests.priors.R | default_prior returns correct priors for multivariate models | fail | cannot transfer | frmtmb takes no list of families as `family` (a multivariate model attaches one family per bf()), family surface outside this lane |
| P18 | tests.priors.R | default_prior returns correct priors for multivariate models | fail | cannot transfer | frmtmb takes no list of families as `family` (a multivariate model attaches one family per bf()), family surface outside this lane |
| P19 | tests.priors.R | default_prior returns correct priors for multivariate models | fail | cannot transfer | frmtmb takes no list of families as `family` (a multivariate model attaches one family per bf()), family surface outside this lane |
| P20 | tests.priors.R | default_prior returns correct priors for categorical models | fail | cannot transfer | categorical() refuses this numeric response ('fewer than two categories'), response-type validation owned by lane wt-famlink |
| P21 | tests.priors.R | set_prior alias functions produce equivalent results | pass | pass |  |
| P22 | tests.priors.R | set_prior alias functions produce equivalent results | pass | pass |  |
| P23 | tests.priors.R | set_prior alias functions produce equivalent results | pass | pass |  |
| P24 | tests.priors.R | set_prior alias functions produce equivalent results | pass | pass |  |
| P25 | tests.priors.R | external interface of validate_prior works correctly | pass | pass |  |
| P26 | tests.priors.R | external interface of validate_prior works correctly | fail | divergence | 10 rows against brms's 9, for the reasons in P1: two class theta rows added, one per-coefficient sd row absent |
| P27 | tests.priors.R | overall intercept priors are adjusted for the intercept | fail | divergence | the fit route is flat. With frmtmb.sample attached, route = "sample" gives student_t(3, 2, 2.5) where brms gives student_t(3, -8, 2.5): the sampling default does not subtract the offset, a filed frmtmb.sample defect not fixed in this lane (dev/priorform-findings.md section 6) |
| P28 | tests.priors.R | as.brmsprior works correctly | pass | pass |  |
| P29 | tests.priors.R | as.brmsprior works correctly | pass | pass |  |
| P30 | tests.priors.R | as.brmsprior works correctly | pass | pass |  |
| P31 | tests.priors.R | as.brmsprior works correctly | pass | pass |  |
| P32 | tests.priors.R | as.brmsprior works correctly | pass | pass |  |
| P33 | tests.priors.R | prior tags are correctly applied | fail | cannot transfer | set_prior() takes no tag, and a brms table carrying one is refused by name, because a tag names a prior inside a Stan program |
| P34 | tests.priors.R | prior tags are correctly applied | fail | cannot transfer | set_prior() takes no tag, and a brms table carrying one is refused by name, because a tag names a prior inside a Stan program |
| P35 | tests.priors.R | prior tags are correctly applied | fail | cannot transfer | set_prior() takes no tag, and a brms table carrying one is refused by name, because a tag names a prior inside a Stan program |
| P36 | tests.priors.R | prior tags are correctly applied | fail | cannot transfer | set_prior() takes no tag, and a brms table carrying one is refused by name, because a tag names a prior inside a Stan program |
| P37 | tests.priors.R | prior tags are correctly applied | fail | cannot transfer | set_prior() takes no tag, and a brms table carrying one is refused by name, because a tag names a prior inside a Stan program |
| B1 | tests.brm.R | brm produces expected errors | pass | pass |  |
| B2 | tests.brm.R | brm produces expected errors | pass | pass |  |
<!-- priorform-ledger-lane:end -->

The 26 passing assertions are ported into
`tests/testthat/test-brms-formula-priors.R`. tests.brmsterms.R was read
and not ported: its seven blocks test brms internals
(`brmsterms()`, `check_re_formula()`, `update_re_terms()`) with no
frmtmb counterpart, except "fixed auxiliary parameters", which is
family validation (lane wt-famlink), and `check_re_formula()` order,
which is section 7 contract 13 and outside this lane.

`test-brms-formula-priors.R` on the base build and on the lane build,
`dev/priorform-run1.R`: <!-- priorform-seefail:begin -->
base build: PRIORFORM brms-formula-priors pass 15 fail 41 error 22 skip 0 in 12 s

lane build: PRIORFORM brms-formula-priors pass 171 fail 0 error 0 skip 0 in 14 s
<!-- priorform-seefail:end -->

## 6. Found and NOT fixed

- **FILED, frmtmb.sample: rescor default on the sampling route.**
  `default_prior(bf(mvbind(y1, y2) ~ x + (x | ID1 | g)) +
  set_rescor(TRUE), dat, family = gaussian(), route = "sample")` gives
  the `rescor` row `(flat)`; brms gives `lkj(1)`. Construction:
  `dev/priorform-rev-route.R` (brms's tests.priors.R data, seed 1),
  log `dev/priorform-rev-route-log.txt`. Ledger P15. Not fixed: another
  lane is editing frmtmb.sample.
- **FILED, frmtmb.sample: offset not subtracted from the default
  intercept location.** `default_prior(y ~ 1 + offset(off), dat,
  route = "sample")` with `y = rep(c(1, 3), each = 5)` and `off = 10`
  gives `student_t(3, 2, 2.5)`; brms gives `student_t(3, -8, 2.5)`,
  because brms subtracts the mean offset from the location. Same
  script and log. Ledger P27. Not fixed, for the same reason.
- **`default_prior` is not an S3 generic.** With brms attached after
  frmtmb, brms's generic masks it and dispatches a frmtmb formula into
  brms. `get_prior()` behaved the same way before this lane.
- **Flat rows read `"(flat)"`** in the table, where brms stores `""`
  and prints `(flat)`. Changing storage breaks the 8 test lines in core
  and frmtmb.sample that compare against `"(flat)"`, for no fit-level
  gain; `frm(prior = )` reads
  both.
- **Coefficient naming** in the table (`s(x).fx1` against `sx_1`,
  `(Intercept)` against `Intercept` for a nonlinear parameter) is the
  naming lane's territory (P11, P12).
- **`categorical()` on a numeric response.** `y2 = c(1, rep(1:3, 3))`
  has three categories and is refused with "the response 'y2' has
  fewer than two categories"; the factor version works. brms accepts
  positive integers. Response-type validation, lane wt-famlink (P20).
- **No list of families** for a multivariate `family =` argument, no
  `threshold = "equidistant"` on `sratio()`, no `cse()` alias, no dpar
  equating (`sigma1 = "sigma2"`), no `update_adterms()`, no prior tags:
  each recorded as cannot transfer in section 5.
- `cs(0 + time || g)` is read as a category-specific term, section 2P
  F7.

## 7. Suites

<!-- priorform-suites:begin -->
frmtmb (core), dev/priorform-core-suite-log.txt:

```
---- GENERATED COUNTS (paste verbatim) ----
package        : frmtmb
files run      : 138
crashed        : 0
pass 8899  fail 1  error 0  skip 133  in 1939 s
files below their 0.57.0 baseline count: 3
                          file pass pass_base skip skip_base
                   test-perf.R    2         3    0         0
           test-prior-compat.R  194       195    0         0
 test-priors-autocor-classes.R   62        63    0         0
baseline files with no run: 0 
---- END GENERATED COUNTS ----
```

frmtmb.coupling, dev/priorform-ext-frmtmb.coupling-log.txt:

```
---- GENERATED COUNTS (paste verbatim) ----
package        : frmtmb.coupling
files run      : 7
crashed        : 0
pass 431  fail 0  error 0  skip 5  in 114 s
files below their 0.57.0 baseline count: 0
baseline files with no run: 0 
---- END GENERATED COUNTS ----
```

frmtmb.eam, dev/priorform-ext-frmtmb.eam-log.txt:

```
---- GENERATED COUNTS (paste verbatim) ----
package        : frmtmb.eam
files run      : 24
crashed        : 0
pass 1633  fail 0  error 0  skip 3  in 1129 s
files below their 0.57.0 baseline count: 0
baseline files with no run: 0 
---- END GENERATED COUNTS ----
```

frmtmb.latent, dev/priorform-ext-frmtmb.latent-log.txt:

```
---- GENERATED COUNTS (paste verbatim) ----
package        : frmtmb.latent
files run      : 7
crashed        : 0
pass 329  fail 0  error 0  skip 2  in 103 s
files below their 0.57.0 baseline count: 0
baseline files with no run: 0 
---- END GENERATED COUNTS ----
```

frmtmb.learn, dev/priorform-ext-frmtmb.learn-log.txt:

```
---- GENERATED COUNTS (paste verbatim) ----
package        : frmtmb.learn
files run      : 12
crashed        : 0
pass 387  fail 0  error 0  skip 13  in 113 s
files below their 0.57.0 baseline count: 0
baseline files with no run: 0 
---- END GENERATED COUNTS ----
```

frmtmb.ode, dev/priorform-ext-frmtmb.ode-log.txt:

```
---- GENERATED COUNTS (paste verbatim) ----
package        : frmtmb.ode
files run      : 8
crashed        : 0
pass 326  fail 0  error 0  skip 1  in 92 s
files below their 0.57.0 baseline count: 0
baseline files with no run: 0 
---- END GENERATED COUNTS ----
```

frmtmb.sample, dev/priorform-ext-frmtmb.sample-log.txt:

```
---- GENERATED COUNTS (paste verbatim) ----
package        : frmtmb.sample
files run      : 19
crashed        : 0
pass 1264  fail 0  error 0  skip 3  in 434 s
files below their 0.57.0 baseline count: 0
baseline files with no run: 0 
---- END GENERATED COUNTS ----
```

frmtmb.spline, dev/priorform-ext-frmtmb.spline-log.txt:

```
---- GENERATED COUNTS (paste verbatim) ----
package        : frmtmb.spline
files run      : 12
crashed        : 0
pass 434  fail 0  error 0  skip 1  in 82 s
files below their 0.57.0 baseline count: 0
baseline files with no run: 0 
---- END GENERATED COUNTS ----
```
<!-- priorform-suites:end -->

Every suite log is later than the last change under `R/` and `tests/`
in core and frmtmb.sample (16:02:56, `test-brms-formula-priors.R`).
The core suite started at 16:05:41, its end time less its duration.
The frmtmb.ode and frmtmb.spline logs were not current: the ode log was
cut off, and the spline log was from 15:41. Both were run again after
the punch round.

The core suite's one failure is `test-perf.R:66`, the bound that a
hundredfold more data costs at most a hundredfold more fit time. It is
a wall-clock bound, and the file passes alone: three reruns gave 3 pass
0 fail each (`dev/priorform-punch2-perf-rerun-log.txt`). That failure
is also why `test-perf.R` is listed below its baseline count.

The two other files below their recorded baseline count dropped by one
assertion each because several `expect_match()` calls on the old print
layout became one `expect_identical()` on the new line, not because an
assertion was lost: `test-prior-compat.R` replaced three with one and
gained one `expect_no_warning()`, `test-priors-autocor-classes.R`
replaced two with one.

`R CMD check --as-cran`, `dev/priorform-check.ps1`, log
`dev/priorform-check-log.txt`, run once in the first pass and once
after the punch round. After the punch round: Status 1 WARNING, 2
NOTEs. The WARNING is "'::' or ':::' import not declared from:
'covr'", from `test-arg-refusal.R` and `test-generic-collision.R`,
both untouched by this lane; the base checkout carries the same
`covr:::` calls. One NOTE is the expected V8 math-rendering one. The
other is an examples-timing NOTE (`VarCorr` 6.42 s, `residuals` 5.33 s
elapsed), taken while the core and extension suites ran in parallel on
the same box; neither example is touched by this lane, and the house
rule is that this NOTE measures load. The first-pass check had no
timing NOTE.

Punch round 2 was not checked again, because it changed no file that
the check reads other than R code and tests. Against the source the
last check built (`dev/priorform-check/frmtmb.Rcheck/00_pkg_src`),
`man/` and `NAMESPACE` are identical, `DESCRIPTION` differs only in the
layout and fields that the build writes, and the other differences are
`R/frame.R`, `R/parse.R`, `R/priors.R` and
`tests/testthat/test-brms-formula-priors.R`. The Rd that round 2 did
change is `?frm_sample`, in frmtmb.sample, and it was rendered and read
back (`dev/priorform-punch2-rd-log.txt`).

The frmtmb.sample suite ran with the pinned StanHeaders library on the
path and with a fresh, empty `FRMTMB_STAN_CACHE`; its sampler tests ran
(`test-sample-direct.R` 136 pass 0 skip), and the cache stayed empty,
because that suite samples through tmbstan and compiles no cached Stan
program.

## 8. Scripts

| script | what it does | output |
|---|---|---|
| `dev/priorform-probe.R` | one construction per item, ref or lane | `dev/priorform-probe-*-log.txt` |
| `dev/priorform-base-behaviour.R` | which model the base fitted for each refused formula | `dev/priorform-base-behaviour-*-log.txt` |
| `dev/priorform-singular.R` | the singular fit behind the stray warning | `dev/priorform-singular-log.txt` |
| `dev/priorform-grdup.R` | brms on the animal model, one column and a copy | stdout |
| `dev/priorform-warncount.R` | warnings one test file emits | `dev/priorform-warncount-log.txt` |
| `dev/priorform-falsealarm.R` | false alarms and misses of the three refusals | `dev/priorform-falsealarm.tsv`, `-log.txt` |
| `dev/priorform-ledger.R` | brms's assertions, run and joined to verdicts | `dev/priorform-ledger-*.tsv` |
| `dev/priorform-fill.py` | pastes every generated block into this file | |
| `dev/priorform-punch-order.R` | fits before and after specificity-wins | `dev/priorform-punch-order-log.txt` |
| `dev/priorform-punch-revert.R` | the punch rounds' guards, reverted, with controls | `dev/priorform-punch-revert-log.txt`, `dev/priorform-punch2-revert-log.txt` |
| `dev/priorform-punch2-libmatch.R` | lane library against the worktree source, core and frmtmb.sample | `dev/priorform-punch2-libmatch-log.txt` |
| `dev/priorform-check.ps1` | R CMD check --as-cran, once | `dev/priorform-check-log.txt` |
| `dev/priorform-install.R` | roxygenise and install into the lane library | |
| `dev/priorform-run1.R`, `dev/priorform-suite.R` | one test file per process | `dev/priorform-suite-*.tsv` |
