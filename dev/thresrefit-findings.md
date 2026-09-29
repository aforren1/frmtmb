# Lane thresrefit: a refit must carry the fitted threshold count

Round of 2026-09-28, punched 2026-09-29. Worktree
`C:\Users\adf44\source\r\frmtmb-wt-thresrefit`, branch `wt-thresrefit`,
private library `C:/Users/adf44/source/r/wt-thresrefit-lib`. "Before" is
the shared reference build of the base commit,
`C:/Users/adf44/source/r/rellib-r3` (frmtmb 0.64.0).

Two defects were assigned. The first is real but reaches ONE path, not
the paths the assignment named. The second is a real latent defect that
no public call reaches, proved by construction in both directions.

**Read "Punch round 1, 2026-09-29" at the end before anything else in
section A.** The reviewer falsified the first version of the fix for
defect A: restoring a factor response's category LABELS is not enough,
because an emptied level also renumbers the codes, and the interior case
was left with the same misaligned row this lane exists to remove. The
response is now RECODED against the fitted categories. The punch round
also made `influence()` count the refits that failed, which fixes a
silent all-NA table that has nothing to do with ordinal thresholds.

## What changed

| file | change |
|---|---|
| `R/thres.R` | `thres_pin_of_fit()` reads a fitted model's threshold layout and category labels; `thres_pin_recode()` puts a factor response's codes back on the fitted model's category scale; `thres_pin_apply()` writes the count into one response's addition-term values as if `thres(x = )` or `thres(x = , gr = )` had been in the formula, and refuses a `thres(gr = )` level that has gained or lost every row |
| `R/frame.R` | `assemble_frame()` takes `thres_pin =`: the recoding runs where the response is extracted, the count where the addition-term values are built, both before `family_finalize()`, which is what reads the count |
| `R/influence.R` | `influence.frmtmb_fit()` passes `thres_pin_of_fit(model)` to every leave-one-out `assemble_frame()` call, and COUNTS the refits that failed: a partly failed table warns with the count and the first reason, an entirely failed one is an error; two doc paragraphs |
| `R/simulate-new.R` | `draw_prior_entry()` refuses an entry whose density covers more than one parameter, instead of writing one draw into all of them; `prior_entry_label()` returns ONE string for such an entry |
| `R/bootstrap.R`, `R/fit.R` | new `?frm_bootstrap` section "Ordinal thresholds"; new paragraphs in the "Ordinal thresholds, thres()" section of `?frm` |
| `tests/testthat/test-thres-refit.R` (new) | 118 expectations in 14 blocks; see Tests |
| `tests/testthat/test-simulate-ergonomics.R` | one comment corrected: it said `frm_simulate()` stops at the natural-scale `newparams` check on an ordinal model, and it does not (see B) |
| `tests/testthat/test-data2.R` | one block asserted the silent all-NA influence table that the failure count replaces; it now asserts the error and its reason |
| `NEWS.md` | new `# frmtmb (development version)` heading with five bullets |
| `dev/thresrefit-probe1.R` .. `probe4.R`, `dev/thresrefit-validate.R`, `dev/thresrefit-priordraw.R`, `dev/thresrefit-agree.R`, `dev/thresrefit-orderedfactor.R`, `dev/thresrefit-of2.R`, `dev/thresrefit-p1-blocking.R`, `dev/thresrefit-p1-controls.R` and their `-log.txt` | evidence |
| `dev/thresrefit-runtest.R`, `dev/thresrefit-suite.sh`, `dev/thresrefit-suite-all.sh`, `dev/thresrefit-block.R`, `dev/thresrefit-install.R`, `dev/thresrefit-news.R`, `dev/thresrefit-rd.R`, `dev/thresrefit-p1-insert.R`, `dev/thresrefit-testlog-all/`, `dev/thresrefit-check/`, `dev/thresrefit-check2/` | harness, test logs and check logs |

### The seam, and why it is where it is

The ordinal threshold count is produced at frame assembly:
`thres_counts()` for a model with `thres()`, and
`ord_tau_init(y, K = max(y))` through `extra_pars` for one without. The
resolved layout is then baked into the FAMILY by `thres_finalizer()`
(`fam[["thres"]]`) and into `frame[["par_template"]]`.

That means a refit that reuses `fit$frame` cannot recount, and a refit
that calls `assemble_frame()` again always does. Following every caller
of `fit_assembled()` and `assemble_frame()`:

| path | rebuilds the frame? | recounts? |
|---|---|---|
| `refit()` (`R/sugar.R`) | no, replaces `frame$y` | no |
| `frm_bootstrap()` | through `refit()` | no |
| `conditional_effects(band = "boot")` | through `frm_bootstrap()` | no |
| `frm_allfit()` | no, reuses `fit$frame` | no |
| `anova(refit = TRUE)` (`anova_refit_ml`) | no | no |
| `confint(method = "profile")` | no | no |
| autoscale pre-fit | no, rescales the frame | no |
| `influence()`, `cooks.distance()`, `dfbeta()`, `dfbetas()` | YES | YES |
| `frm_multiple()` | one `frm()` per imputation | yes, and no fitted model exists to pin to |
| `update(newdata = )` | re-evaluates the stored `frm()` call | yes, by design |
| `frm_simulate()` then `frm()` | two user calls | yes, by design |

So the pin is passed at exactly one call site. `assemble_frame()` grew
the argument rather than `influence()` growing the logic, because the
count is consumed inside frame assembly and because the next caller that
needs it (a future refit that subsets rows) gets it for free.

`thres_pin_apply()` writes `av[["thres"]]`, which is the declared
spelling of a pinned count, so nothing downstream needed a new branch:
the four ordinal families already declare
`accepts_aterms = c("weights", "thres")`, and `family_finalize()`
already resolves a count into the layout. A count the user wrote wins
over the pinned COUNT (`if (!is.null(av[["thres"]])) return(av)`), but
not over the recoding, which is about what a subset's codes mean rather
than about the model. The gate on the family is its own allow-list, so an
ordinal family from an extension that does not declare `thres` is left
alone; `test-thres-refit.R` constructs one.

The pin carries LABELS throughout, never codes, because
`drop.unused.levels = TRUE` and `thres_group_codes()`'s `factor()` both
renumber when a level loses every row. The response's own labels give the
recoding its map; the `thres(gr = )` labels keep level b's thresholds out
of level a's slice, and make "this level has no rows left" a refusal
rather than a silent shift.

## A. Refits that lose a threshold

### The construction

`dev/thresrefit-probe2.R`, log beside it. n = 40, four categories,
logit, slope 0.5, thresholds (-0.6, 0.5, 2.6); data seed 202 gives
`table(y) = 9/16/11/4`. `simulate(fit, nsim = 60, re_formula = NA)` at
seed 11 produces exactly one replicate (number 31) whose response never
reaches category 4: `table = 8/15/17/0` on cumulative and `8/16/16/0`
on sratio.

Grouped: `dev/thresrefit-probe2.R` builds 90 rows over three levels,
seed 303, per-level counts 3/2/3; at simulate seed 13, 17 of 60
replicates have a level that loses its top category.

### What happens today on the reference build

**`frm_bootstrap()` and `refit()` do NOT have the defect.** Measured on
`rellib-r3` (`dev/thresrefit-validate-before-log.txt`, section 3):

| construction | refit logLik | `thres(3)` pinned logLik | relative difference |
|---|---|---|---|
| cumulative, replicate 31 | -41.4680285630979 | -41.4680285623368 | 1.835e-11 |
| sratio, replicate 31 | -41.9541926013013 | -41.9541926022876 | 2.351e-11 |
| grouped cumulative, replicate 6, seed 503 | -99.7389620475353 | -99.738962045463 | 2.078e-11 |

and in every case `length(tau_raw)` is 3 (7 when grouped) on both
sides. `frm_bootstrap(FUN = c(fixef, tau_raw))` returns a 60 x 4 matrix
with column names `x, tau1, tau2, tau3` and NO NA cells, on both arms.
The residual 2e-11 is optimizer tolerance, not a model difference: the
two fits start from different points (a warm start at the fitted
estimates against `frm()`'s own start values).

The entry in `dev/thres-findings.md` ("Refits whose data lose a
category") and the assignment that repeats it are therefore wrong about
`frm_bootstrap()`. Task items 2 and 3 asked for a change there; the
measurement says there is nothing to change, and the property is now
pinned by a test instead. For the record, a fresh `frm()` on replicate
31's data fits 2 thresholds and reaches the SAME log-likelihood
(-41.4680285592851 against -41.4680285623368): the threshold above the
highest observed category contributes probability 0 to every row, so
losing it costs no likelihood and only changes the parameter vector.
That is why the defect is silent.

**`influence()` DOES have the defect**, in two forms.
`dev/thresrefit-validate.R`, data seed 501, n = 50, four categories with
EXACTLY ONE row in category 4 (row 1), logit, slope 0.5.
`table(y) = 20/11/18/1`. Before
(`dev/thresrefit-validate-before-log.txt`, section 1):

| family | NA cells in `inf$fixed` | `cooks.distance()` at row 1 |
|---|---|---|
| cumulative | 1 | NA |
| sratio | 1 | NA |

and the three finite values in that row come from a TWO-threshold model
while the column names still say `tau_raw_1, tau_raw_2, tau_raw_3`.

Grouped is worse, because the loss is positional. Data seed 502, n = 96,
three levels, fitted per-level counts 3/2/2 (`tau_raw` length 7), level
"a" alone reaching category 4 and there on one row (row 10). Before
(section 2), row 10 reads

    0.50921167 -0.91462934 0.32926521 -0.55576307 0.53683315
    -0.66443626 0.79322835 NA

against the `thres(k, gr = g)` reference

    0.50921975 -0.91463523 0.32926979  3.09096292 -0.55579306
     0.53684125 -0.66445177 0.79322886

Every coefficient from position 4 onwards is SHIFTED ONE SLOT EARLIER:
the table reports level b's first threshold in level a's third
threshold's column. The only visible sign is the trailing NA.

### After

`dev/thresrefit-validate-after-log.txt`. Ungrouped: 0 NA cells in the
whole 50 x 4 table, 0 NA in `cooks.distance()`, and
`cooks.distance()[1]` is 10.512458 on cumulative and 83.588491 on
sratio (NA before). Row 1 against the `thres(3)` reference, per
coefficient (`dev/thresrefit-agree.R`, log
`dev/thresrefit-agree-log.txt`):

| family | coefficient | pinned refit | `thres(3)` reference | relative | `|diff| / sd(column)` |
|---|---|---|---|---|---|
| cumulative | x | 0.000627432768754419 | 0.000627417852989218 | 2.4e-05 | 3.1e-07 |
| cumulative | tau_raw_1 | -0.371346624172049 | -0.371346636047913 | 3.2e-08 | 2.5e-07 |
| cumulative | tau_raw_2 | -0.0886355345190902 | -0.0886355354867559 | 1.1e-08 | 2.3e-08 |
| cumulative | tau_raw_3 | 3.11120658611246 | 3.1231163278655 | 3.8e-03 | 4.4e-02 |
| sratio | x | -0.0749764108220967 | -0.0749764109267025 | 1.4e-09 | 2.7e-09 |
| sratio | tau_raw_1 | -0.393802420902748 | -0.393802421005823 | 2.6e-10 | 2.3e-09 |
| sratio | tau_raw_2 | -0.519016999258604 | -0.519016999063165 | 3.8e-10 | 3.3e-09 |
| sratio | tau_raw_3 | 21.6608807211552 | 22.2415759656778 | 2.6e-02 | 2.2e-01 |

The 0.0038 and 0.026 figures printed in `-validate-after-log.txt` are
`max` over ALL coefficients, and the table shows the unidentified
`tau_raw_3` is what produces them. Every IDENTIFIED coefficient agrees
to 2.4e-05 relative at worst (and that worst case is `x`, whose value is
6.3e-04, so the relative figure is inflated by the small denominator);
in units of the spread the deletions themselves produce, the worst
identified coefficient differs by 3.1e-07 of one standard deviation.
The unidentified threshold is the point of the exercise: once the only
row in the top category is gone, that threshold bounds a category
nobody chose, its ML estimate runs off toward infinity, and two
optimizations stop at different places on a flat ridge. The test
therefore compares the identified coefficients against the spread the
deletions themselves produce in each column (a quantity the run
measures) and asserts only that the unidentified one is large in units
of the largest identified threshold, on both sides.

Grouped, after: 0 NA cells, `cooks.distance()` at row 10 is 4.6167013
(NA before), and row 10 is

    0.50922013 -0.91463509 0.32927020 3.06332923 -0.55579326
    0.53684091 -0.66445157 0.79322898

which is the reference position by position. Per coefficient
(`dev/thresrefit-agree-log.txt`), the seven identified coefficients
agree to 1.3e-06 relative at worst and to 1.7e-05 of one column
standard deviation; the unidentified one, level a's third threshold, is
3.06332923400752 against 3.09096291714765 (8.9e-03 relative, 0.14 of
one standard deviation).

### The ordered-factor response, and a regression the first fix caused

Found by construction after the first version of the fix was installed,
not in review: `dev/thresrefit-orderedfactor.R`, logs
`dev/thresrefit-orderedfactor-before-log.txt` and `-after-log.txt`, and
`dev/thresrefit-of2.R` for the diagnosis.

The same data as section 1 with the response as
`factor(y, levels = 1:4, ordered = TRUE)`. `model.frame()` DROPS a factor
level that no row takes, so the subset's response has three levels even
though the column declares four; `frame.R` then refuses a pinned count of
3 by name ("thres(x = 3) asks for 4 categories, and the response is an
ordered factor with 3 levels in the data"). Inside `influence()`'s
per-unit `tryCatch` that meant the WHOLE row came back NA: 4 NA cells
against the base build's 1.

The first fix carried the fitted model's category LABELS with the pin and
restored them only when the subset's remaining labels were an initial
segment of the fitted ones. That closed the top-category case: 0 NA
cells, `cooks.distance()[1] = 10.512458`, and row 1

    0.0006274327688 -0.3713466241720 -0.0886355345191 3.1112065861125

identical to the integer-coded row. **It was the wrong fix, and the
reviewer falsified it.** See "Punch round 1" below: restoring labels is
not enough, because the dropped level also RENUMBERS the codes, and the
initial-segment guard left the interior case with the very misalignment
this lane exists to remove. The final fix recodes the response instead.

### `categorical()` and multinomial-type responses

Checked, and they do NOT have this defect in any in-package refit path.
The reason is structural: a categorical family's category count is
resolved by `resolve_deferred_families()` inside `frm()` and stored on
the SPEC, and every in-package refit reuses `fit$spec`.
`dev/thresrefit-validate.R` section 4, data seed 504, n = 50, three
levels with EXACTLY ONE row in level 3 (row 36): `influence(force =
TRUE)` returns a 50 x 4 table with 0 NA cells on BOTH arms, and
`frm_bootstrap(nsim = 20, seed = 21)` a 20 x 4 matrix with 0 NA cells on
both arms. Deleting row 36 leaves category 3 with no rows, so its
predictor is unidentified and the row reads
`mu3_(Intercept) = -76.75, mu3_x = -23.54`; that is the same
unidentified-parameter phenomenon as the ordinal one, it is identical on
both arms, and no pin can fix it because there is no data. Pinned as a
test (the complement of the ordinal case, checked rather than assumed).

### `frm_multiple()`: measured, and deliberately NOT changed

`dev/thresrefit-probe4.R`, log beside it. Three imputations, the second
with its top category collapsed (tables 17/17/20/6 and 17/17/26/0), and
`frm_multiple(bf(y ~ x + z), family = cumulative())` gives per-fit
`tau_raw` lengths 3, 2, 3. What happens:

- `$pooled` has rows `x` and `z` only. Rubin pooling runs over
  `estimated_coef_names()`, which is `beta` and `betad`; ordinal
  thresholds are not pooled at all, on any imputation set. That is a
  pre-existing gap, not this defect, and it is out of scope.
- `anova(method = "D3")` already REFUSES, by name: "the imputations
  produced different parameter vectors; method = \"D3\" needs one
  common parameterization". The control with three identical
  imputations runs.
- `anova(method = "D1")` and `"D2"` run, and are not wrong for the
  reason above: the absent threshold carries no likelihood.
- `hypothesis()` on the pooled fit runs.

So the one thing that breaks is already named, and there is no fitted
model to take a count from: `frm_multiple()` fits a list of data frames
with no reference fit. Pinning to imputation 1 would refuse a
legitimate set where a later imputation reaches a higher category, and
pinning to the union would need `frm()` to grow a user-facing argument
for an internal purpose. Not changed.

### `update(newdata = )`: measured, and deliberately NOT changed

`dev/thresrefit-probe3.R` section A4: `update(fit, newdata = )` on data
whose top category was collapsed fits 2 thresholds where the original
fit had 3. This is not a refit. `?update.frmtmb_fit` documents it as
"Re-evaluates the stored `frm()` call with the given arguments
replaced", so `update(fit, newdata = d)` must agree with writing
`frm(..., data = d)` by hand; injecting the old fit's count would make
the two differ and would be undiscoverable from the call. Recorded, not
changed. The `?frm` section now says this explicitly.

### A limitation the pin cannot remove

`influence(groups = g)` on a model with `thres(gr = g)` deletes a whole
level. That level then has no rows, and the thresholds the fitted model
gives it have no data at all, so they could only be invented. The unit
is refused by name and its row is NA. What this lane got wrong at first
was the other half: dropping the level instead slid every LATER level's
thresholds one column early, and the reviewer measured it. See "Punch
round 1".

## B. `draw_prior_entry()` on an unordered threshold vector

### Reachability, both directions

`dev/thresrefit-priordraw.R`, logs
`dev/thresrefit-priordraw-before-log.txt` and `-after-log.txt`.

**The defect is real at the function.** On the reference build, with an
entry `list(comp = "tau_raw", idx = 1:3, scale = "internal", dist =
prior_normal(0, 2))` and `set.seed(708)`:

    draw_prior_entry(idx = 1:3) -> returned one number, -2.563726679
    draw_prior_pars(idx = 1:3)  -> tau_raw = -0.7271876417
                                             -0.7271876417
                                             -0.7271876417
                                   | all equal: TRUE

**No public call reaches it.** `draw_prior_pars()` has exactly one
caller. Inside frmtmb, found by deparsing the body of every function in
the namespace (`dev/thresrefit-probe4.R`): `frm_simulate` alone, and
`prior_entry_label()` has two, `check_coverage` and `frm_simulate`.
Outside it, `grep -rn "draw_prior|prior_entry_label"` over every
extension's `R/` and `tests/` matches nothing. And `frm_simulate()`
builds its labels
BEFORE the draw, with
`vapply(entries, function(e) prior_entry_label(frame, slots, e), "")`.
On the reference build `prior_entry_label()` returns one name per index,
so a multi-index entry is length 3 and `vapply` errors first:

    K=3 cratio     -> ERROR: values must be length 1,
                       but FUN(X[[1]]) result is length 3
    K=3 acat       -> same
    K=3 sratio     -> same
    K=3 cumulative -> same

for `frm_simulate(bf(y ~ x), family = ..., data = <40 rows, y = rep(1:4,
10)>, prior = normal(0, 2) on class "Intercept" + normal(0, 1) on class
"b"), seed 702`. A multi-index entry is the only case that can recycle,
so the recycling is unreachable from outside the namespace.

**The census of multi-index entries, corrected.** The first version of
this document said `rescor` was "the other multi-index entry the
resolver can build". That was wrong, and the reviewer found two more,
BOTH reachable through `frm_simulate()`
(`dev/thresrefit-rev-07-multiidx.R`):

| model | entry | base | after |
|---|---|---|---|
| `(1 + x + z \| g)` with `set_prior("lkj(2)", class = "cor")` | `theta` idx 4,5,6 | `simpleError: 'length = 3' in coercion to 'logical(1)'` | named refusal on `theta[4:6]` |
| `ar(time = t, gr = gg, p = 2)` with `set_prior("normal(0, 0.3)", class = "ar")` | `thetaac` idx 1,2 | `simpleError: values must be length 1, but FUN(X[[4]]) result is length 2` | named refusal on `thetaac[1:2]` |

So the rule is about the SHAPE of a prior entry and not about one
family: class `"Intercept"` on an ordinal model with more than one
threshold, class `"cor"` on a block with more than one correlation,
class `"ar"`, `"ma"` or `"cortime"` above order 1, and class `"rescor"`.
`rescor` alone is unreachable here, because `frm_simulate()` refuses
multivariate models by name. The reviewer also measured the fit path on
both of those models as bit-identical between arms (logLik
-175.98786419096427 and -123.63837585947768), so nothing but the
message changed. The NEWS bullet and the comments on
`draw_prior_entry()` and `prior_entry_label()` now describe the rule.

**The complement, so the claim is not "no test reaches it".** An entry
with ONE index draws, correctly, on both arms and bit-identically:

| construction | before | after |
|---|---|---|
| `draw_prior_pars(idx = 1)`, seed 709 | tau_raw = 0.6202043782 | 0.6202043782 |
| `frm_simulate`, two-category cratio / acat / sratio, seeds 703 / 704 | tau_raw_1 = 2.8637285240, b_x = 0.9604788337 | identical |
| `frm_simulate`, `thres(gr = g)` with one threshold per level, seeds 705 / 706, cratio | tau_raw_1 = -1.3510929168, tau_raw_2 = -1.3315166704, b_x = -0.4285386843 | identical |

so the guard added below does not fail closed, and grouped entries that
are one index each still get one draw apiece.

### What was changed, and what was not

The recycling itself is now refused rather than performed, and
`prior_entry_label()` returns one string whatever the entry covers, so
the public call says what it cannot do instead of dying inside
`vapply()`:

    A prior on tau_raw[1:3] is a density on 3 parameters at once, and
    frm_simulate() reports one drawn value per prior. Pin them with
    newparams = tau_raw = instead

`cumulative()` keeps its own, older refusal, which names the ordering
(`scale == "ordthres"`), and that refusal still fires on a
single-threshold cumulative model, unchanged.

Not done: making the draw WORK for the unordered families, which brms
can do (it declares sratio, cratio and acat thresholds as `vector`, so
each element takes the prior independently). That needs one prior entry
per threshold, or a report column per index, and either one reaches the
FIT path: `neg_log_prior_fn()` and `prior_summary()` read the same
entries, and `test-thres.R` pins their objective to
`64 * eps * |objective|`. It is its own change. A capability frmtmb
already declines for `cumulative()` is not worth the risk to a validated
prior objective in this lane.

Also fixed on the way, because it is the same line: on a model with
random effects, `prior_entry_label()`'s sd search evaluates
`e$idx %in% nat_sd_theta_idx(...)` inside `if`, which is an error in R
4.2 and later when `idx` has length above one. The new length-one
branch returns before it.

## Defects found and NOT fixed

- Ordinal thresholds are not pooled by `frm_multiple()` at all (the
  Rubin table covers `beta` and `betad` only). Pre-existing, and
  nothing in this lane's scope reaches it.
- `frm_multiple()` over imputations with different threshold counts
  produces per-fit parameter vectors of different lengths; only
  `anova(method = "D3")` notices, and it refuses by name. See above for
  why it was left.
- `update(newdata = )` recounts. By design; see above.

## Tests

As the file stood at the END of round 0, before the punch round:
`tests/testthat/test-thres-refit.R`, 71 expectations in 8 blocks, seen
to FAIL on the reference build at `pass=43 fail=15 err=2 skip=0`. The
punch round took it to 118 expectations in 14 blocks and the base arm to
`pass=68 fail=23 err=5 skip=0`; those are the current figures, and
`dev/thresrefit-seen-failing-base.txt` holds the current run. The round-0
failures, kept because they are the record of what the first fix pinned:

    a leave-one-out refit keeps the fitted threshold count
      Expected `anyNA(inf$fixed)` to be FALSE. actual: TRUE      (x2)
      Expected `anyNA(cooks.distance(inf))` to be FALSE           (x2)
      Expected `abs(got[["tau_raw_3"]])/scale_id` > 3: NA <= 3    (x2)
    a leave-one-out refit keeps the per-level threshold counts
      Expected `anyNA(inf$fixed)` to be FALSE. actual: TRUE
      Expected `anyNA(cooks.distance(inf))` to be FALSE
      Expected max(|got - want| / spread) < 0.001: NA >= 0.001
      Expected `abs(got[[unident]])/scale_id` > 3: 0.6 <= 3.0
    thres_pin_of_fit() finds nothing to pin off an ordinal model
      Error: object 'thres_pin_of_fit' not found        (the weak form)
    a leave-one-out refit keeps an ordered factor's categories
      Expected `anyNA(inf$fixed)` to be FALSE. actual: TRUE
      Expected `anyNA(cooks.distance(inf))` to be FALSE
      Expected max(|ordered - integer| / spread) < 0.001: NA >= 0.001
    a prior draw never writes one value into a threshold vector
      Expected draw_prior_pars(...) to throw a error
      Expected prior_entry_label(...) to have length 1. Actual: 3
      Error in vapply(...): values must be length 1, but FUN(X[[1]])
        result is length 3

The grouped `0.6 <= 3.0` line is the positional shift caught directly:
on the base build the column that should hold level a's unidentified
third threshold (about 3.06) holds level b's first threshold (-0.56).

No absolute numeric tolerance is written. Coefficient agreement is
`|difference| / sd(column of the influence table) < 1e-3`, both sides
measured by the run; log-likelihood agreement is
`|difference| / |log-likelihood| < sqrt(eps)`; the unidentified
threshold is compared to the largest identified threshold of the same
fit.

## Tests run

The 44 test files the change can reach, one per R process,
`NOT_CRAN=true`, up to 7 processes at a time. `dev/thresrefit-suite.sh`
is the shorter runner used while iterating, and it is the ONE place the
list of files lives: `dev/thresrefit-block.R` reads the list out of that
script and writes the block below into this document, rather than
repeating the list or leaving the counts to be typed. The counts
themselves are the final ones, from the whole-suite logs in
`dev/thresrefit-testlog-all/`. They sum to 5618 passing expectations.

<!-- generated by dev/thresrefit-block.R; do not edit by hand -->
    test-api-spellings.R           pass=36 fail=0 err=0 skip=0
    test-arg-refusal.R             pass=123 fail=0 err=0 skip=0
    test-autoscale.R               pass=46 fail=0 err=0 skip=0
    test-boot.R                    pass=49 fail=0 err=0 skip=0
    test-bracket-access.R          pass=33 fail=0 err=0 skip=0
    test-brms-agreement.R          pass=168 fail=0 err=0 skip=2
    test-brms-families.R           pass=1033 fail=0 err=0 skip=0
    test-brms-shapes-punch.R       pass=55 fail=0 err=0 skip=0
    test-brms-shapes.R             pass=72 fail=0 err=0 skip=0
    test-case-studies.R            pass=30 fail=0 err=0 skip=0
    test-ce-bands.R                pass=167 fail=0 err=0 skip=0
    test-compat.R                  pass=639 fail=0 err=0 skip=0
    test-data2.R                   pass=31 fail=0 err=0 skip=0
    test-diagnostics-ux.R          pass=127 fail=0 err=0 skip=0
    test-edgecases.R               pass=64 fail=0 err=0 skip=0
    test-frame.R                   pass=22 fail=0 err=0 skip=0
    test-importance.R              pass=226 fail=0 err=0 skip=0
    test-influence-plot.R          pass=21 fail=0 err=0 skip=0
    test-input-validation.R        pass=43 fail=0 err=0 skip=0
    test-methods-audit.R           pass=58 fail=0 err=0 skip=0
    test-multiple-pooling.R        pass=33 fail=0 err=0 skip=0
    test-mv-gaps.R                 pass=46 fail=0 err=0 skip=0
    test-numerical-robustness.R    pass=709 fail=0 err=0 skip=0
    test-ordinal-fitted.R          pass=129 fail=0 err=0 skip=0
    test-ordinal.R                 pass=109 fail=0 err=0 skip=0
    test-osa-inference.R           pass=34 fail=0 err=0 skip=0
    test-predfix.R                 pass=92 fail=0 err=0 skip=0
    test-review-fixes.R            pass=18 fail=0 err=0 skip=0
    test-review-v29.R              pass=137 fail=0 err=0 skip=0
    test-setprior.R                pass=27 fail=0 err=0 skip=0
    test-simulate-density.R        pass=459 fail=0 err=0 skip=0
    test-simulate-ergonomics.R     pass=50 fail=0 err=0 skip=0
    test-simulate-newdata.R        pass=44 fail=0 err=0 skip=0
    test-spectral.R                pass=113 fail=0 err=0 skip=0
    test-sratio-thresholds.R       pass=29 fail=0 err=0 skip=0
    test-structure.R               pass=135 fail=0 err=0 skip=0
    test-tabular-inputs.R          pass=9 fail=0 err=0 skip=0
    test-thres-refit.R             pass=118 fail=0 err=0 skip=0
    test-thres.R                   pass=67 fail=0 err=0 skip=0
    test-unpinned-seams.R          pass=30 fail=0 err=0 skip=0
    test-v07.R                     pass=25 fail=0 err=0 skip=0
    test-v14.R                     pass=81 fail=0 err=0 skip=0
    test-v15.R                     pass=41 fail=0 err=0 skip=0
    test-verbose.R                 pass=40 fail=0 err=0 skip=0

## Whole core suite

Run once, at the end: `dev/thresrefit-suite-all.sh`, log
`dev/thresrefit-suite-all-log.txt`, per-file logs in
`dev/thresrefit-testlog-all/`. Every one of the 180 `test-*.R` files in
one R process of its own, 7 at a time, ungated.

Rerun in full after the punch round, on the final code.

<!-- generated by dev/thresrefit-suite-all.sh; paste verbatim -->
    launched 180 files
    files=180 pass=11997 fail=0 err=0 skip=152
    files with a nonzero fail or err: none
    aborted: 0

The count rose from 11949 to 11997 with the punch round's 47 new
expectations in `test-thres-refit.R` and its one net new expectation in
`test-data2.R`.

The 152 skips are the gated tiers, which this run does not enable:
test-brms-likelihood.R 37, test-brms-methods.R 46, test-drmtmb-agreement.R
13, test-brms-priors.R 12, the thirteen test-bcm-*.R files 39 between
them, test-brms-agreement.R 2, test-rl-example.R 2, test-fuzz.R 1.

The extensions' own suites were NOT run. The change is core, and no
extension reaches either seam: `grep -rn
"assemble_frame|draw_prior|prior_entry_label|thres_pin"` over
`extensions/*/R/` and `extensions/*/tests/` matches three COMMENTS in
frmtmb.eam and no code. `thres_pin_of_fit()` also returns `NULL` for any
response whose family type is not `"ordinal"`, so an extension family
cannot be pinned.

## R CMD check

**The full run**, once: `R CMD build` then
`R CMD check --as-cran frmtmb_0.64.0.tar.gz`, log
`dev/thresrefit-check/check.log`. `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`
and `_R_CHECK_FORCE_SUGGESTS_=false`.

    Status: 1 ERROR, 1 WARNING, 2 NOTEs

Every content section is OK, including `checking Rd files`,
`checking R code for possible problems` (28s), `checking examples` (31s),
`checking examples with --run-donttest` (37s), `checking tests` (366s)
and `checking re-building of vignette outputs` (224s). All four non-OK
items are this machine, not the package:

- `checking PDF version of manual ... WARNING` and
  `checking PDF version of manual without index ... ERROR`:
  `pdflatex is not available`. There is NO TeX on this box.
  `where /R C:\Users\adf44 pdflatex.exe` returns nothing and
  `tinytex::tinytex_root()` is `""`, so the TinyTeX path
  `lane-rules.md` gives no longer exists. The release harness's own
  `dev/release/check.log` has `checking PDF version of manual ... OK`,
  so this is a change in the machine.
- `checking for non-standard things in the check directory ... NOTE`:
  `frmtmb-manual.tex`, which the failed `texi2pdf()` left behind. A
  consequence of the item above.
- `checking HTML version of manual ... NOTE`: "Skipping checking math
  rendering: package 'V8' unavailable". This is the one expected NOTE
  `lane-rules.md` names.

`_R_CHECK_FORCE_SUGGESTS_=false` is needed because `drmTMB` is not
installed on this machine at all (a `Suggests` at `>= 0.7.0`); without
it the check stops at `checking package dependencies ... ERROR` after
30 seconds and runs nothing. `test-drmtmb-agreement.R` skips its 13
expectations for the same reason.

**A partial second run.** The full run above was started before the
ordered-factor regression was found, and the punch round changed more
code after it, so it checked a tree that differs from the delivered one
in `R/thres.R`, `R/influence.R`, two places in `R/frame.R`, the doc
sections of `R/bootstrap.R` and `R/fit.R`, all three touched Rd files,
`test-thres-refit.R` and `test-data2.R`. Rather than spend a second
20-minute run, the sections that delta can affect were rechecked on the
FINAL tree with the long sections skipped
(`dev/thresrefit-check2/check.log`, rerun after the punch round):

    R CMD check --as-cran --no-tests --no-examples --no-vignettes \
                --no-build-vignettes --no-manual
    Status: 2 WARNINGs, 1 NOTE

`checking Rd files`, `Rd metadata`, `Rd line widths`,
`Rd cross-references`, `Rd \usage sections`, `Rd contents`,
`code/documentation mismatches`, `S3 generic/method consistency`,
`use of S3 registration`, `dependencies in R code` and `R code for
possible problems` (26s) are all OK. Both WARNINGs and the NOTE are
`--no-build-vignettes` itself ("Files in the 'vignettes' directory but no
files in 'inst/doc'", "Directory 'inst/doc' does not exist", "no prebuilt
vignette index"). The delta's own evidence is elsewhere:
`test-thres-refit.R` passes 118 of 118 in its own process, the whole
180-file suite passes on the final code, and all three changed Rd
sections were verified by RENDERING them with `Rd2txt` and scanning the
rendered text (`dev/thresrefit-rd.R`, log
`dev/thresrefit-rd-log.txt`: "style problems: none"). The widest
rendered line of `man/frm.Rd` is 73 columns.

**An instrument of my own, wrong twice and recorded rather than
deleted.** That scanner reported "1 rendered line over 80" in the punch
round, on the section heading "The Laplace approximation, and how to
check it:", and I excused it by checking that `git show HEAD:man/frm.Rd`
reported the same one line. Both numbers were artifacts of the scanner.
`Rd2txt` underlines a heading the teletype way, with an underscore, a
BACKSPACE and the character (`_\bT_\bh_\be`), and my stripper removed the
underscores and left the backspaces, so a 46-column heading measured 85.
Two symptoms should have caught it sooner: the "wide" line was always a
heading and never prose, and a `grepl()` for the heading's own text
failed against the supposedly stripped string. The stripper now removes
the `_\b` pairs, no line of any of the three files exceeds 80, and the
exemption is gone. The scanner is checked in the other direction too: a
95-character line is still flagged.


## Punch round 1, 2026-09-29

Review: `dev/reviews/2026-09-29-thresrefit.md`, verdict NOT MERGEABLE.
Every finding was reproduced before it was fixed. Scripts:
`dev/thresrefit-p1-blocking.R` and `dev/thresrefit-p1-controls.R`, each
with a `-before-log.txt` (arm `base`, 0.64.0) and an `-after-log.txt`
(arm `lane`).

### Blocking 1. The first fix was the wrong fix

The reviewer falsified the claim that an interior loss is "left
unpinned". Leaving it unpinned means leaving it MISALIGNED, which is the
signature this lane exists to remove.

The diagnosis the first fix missed: `assemble_frame()` builds its model
frame with `drop.unused.levels = TRUE`, so an emptied factor level does
not merely disappear from the labels, it RENUMBERS every category above
it. With four categories and an empty second, a row truly in the third is
coded 2. Restoring the labels cannot repair that; the codes have to be
put back on the fitted scale.

So the fix is now `thres_pin_recode()`, which runs on the response
immediately after `extract_y()` and before anything reads it, and maps
each subset code through the fitted model's own category labels. An
interior category then simply has no observations, exactly as it has none
when the same data are coded as integers.

Option (a) of the punch instruction, chosen. Option (b) is used only
where a refit genuinely cannot represent the fitted model, and never
produces a shifted row.

**Ungrouped, the reviewer's construction.** Seed 901, n = 60,
`table(y) = 36/1/13/10`, one row (row 1) in the interior category 2.
Deleting it empties category 2. The `thres(3)` reference on that subset
is

    0.356690563  0.534975255  -21.292348854  0.162248682

| arm | coding | NA in the deleted row | the row |
|---|---|---|---|
| base | integer | 0 | 0.356690563  0.534975255  -21.449994726  0.162248682 |
| base | ordered factor | 1 | 0.356690854  0.534977288  0.162248984  NA |
| base | character | 0 | 0.356690563  0.534975255  -21.449994726  0.162248682 |
| lane | integer | 0 | 0.356690563  0.534975255  -21.449994726  0.162248682 |
| lane | ordered factor | 0 | 0.356690563  0.534975255  -21.449994726  0.162248682 |
| lane | character | 0 | 0.356690563  0.534975255  -21.449994726  0.162248682 |

On `base` the ordered factor's `tau_raw_2` column holds `tau_raw_3`'s
value (0.162248984 against the reference's 0.162248682) and
`identical(ordered, integer)` over the WHOLE table is FALSE. On `lane`
all three codings are bitwise identical, `identical()` TRUE on the whole
60 x 4 table, the emptied category's unidentified threshold sits in its
own column (-21.45 against the reference's -21.29, a flat ridge), and
`cooks.distance()` has no NA.

**Grouped.** My first grouped construction (seed 902 here) did not
reproduce the reviewer's, and the reason is worth recording: a category
absent from ONE `thres(gr = )` level is still present in the response, so
the factor keeps its level and even `base` is aligned (seed 902, base and
lane both 0 NA and bitwise identical). The level only disappears when the
category is absent GLOBALLY. Constructed at seed 904, n = 96, category 2
present on exactly one row, in level a, `nthres = 3/2/2`, `tau_raw`
length 7:

| arm | coding | NA | the deleted row |
|---|---|---|---|
| base | integer | 0 | 0.371900709 0.138229085 -20.760649158 0.633160001 2.009694234 -20.096186900 0.223434960 -29.178326137 |
| base | ordered factor | 3 | 0.371900550 0.138229648 0.633157812 2.009681458 0.223434410 NA NA NA |
| lane | every coding | 0 | 0.371900709 0.138229085 -20.760649158 0.633160001 2.009694234 -20.096186900 0.223434960 -29.178326137 |

Five coefficients one slot early with three trailing NA on `base`,
exactly the reviewer's signature; bitwise identical across codings on
`lane`.

**What cannot be pinned is refused, never shifted.** A `thres(gr = )`
level with no rows left has thresholds with no data at all.
`thres_pin_apply()` now refuses that unit by name ("cannot drop a
level") rather than dropping the level, which is what slid the later
levels' thresholds. `influence(groups = )` on a model whose grouping
factor IS the threshold grouping factor therefore errors, where `base`
returned a 3 x 8 table with 7 of 24 cells NA and the rest misaligned.

### Blocking 2. influence() now counts its failures

`influence.frmtmb_fit()` counts the deletion refits that returned an
error and keeps the first message and unit. A partly failed table warns;
a table whose every refit failed is an error carrying that message.

| construction | base | lane |
|---|---|---|
| `data = ` with two rows in a category the fit never saw (seed 1101, n = 50, fit `n_tau = 2`) | `subscriptOutOfBoundsError`, "subscript out of bounds" | error: "influence(): all 50 deletion refits failed ... The first was observation 1: Number of thresholds is smaller than required by the response ... reaches category 4" |
| `data = ` with a fourth `thres(gr = )` level "d" (seed 1103, k = 96) | `subscriptOutOfBoundsError` | error: "all 96 deletion refits failed ... level(s) 'd' ... are not in the fitted model ... cannot add a level" |
| `groups = ` deleting a whole `thres(gr = )` level | 3 x 8 table, 7 of 24 cells NA, silent, misaligned | error: "all 3 deletion refits failed ... The first was 'idg' level a: ... have no rows left, so the 3 threshold(s) the fitted model gives them cannot be estimated" |
| ONE row in an unseen category (seed 1202, n = 40), so one unit refits | `subscriptOutOfBoundsError` | warning: "39 of 40 deletion refits failed and their rows are NA. The first was observation 1: ...", and the table comes back with row 7 filled |

The findings text that said "that unit is skipped" was wrong by a factor
of the sample size, and the second construction was not recorded at all.
Both are above.

**A pre-existing silent all-NA table, also fixed, and not ordinal.**
`tests/testthat/test-data2.R` block "data2 outlives its calling
environment across saveRDS" ASSERTED the old behavior:
`expect_true(all(is.na(i_env$fixed)))` on a fit whose structural matrix
was lost with its environment. Every refit fails there for a reason that
has nothing to do with thresholds, and 0.64.0 returned a 6 x 3 matrix of
NA without a word. That block now asserts the error and its reason
("all 6 deletion refits failed", "cannot find 'A'"), and both new
expectations are SEEN TO FAIL on `base`
(`pass=29 fail=2 err=0 skip=0`).

### Finding 3. A hand-written thres(K) on a factor response

`thres_pin_apply()` returned at `!is.null(av[["thres"]])` before the
label work, so the one coding with a user-written count came back empty.
The recoding is now a separate step that runs first, which is right for a
second reason: a count the user wrote is a statement about the MODEL, and
the recoding is about what a subset's codes MEAN. Seed 501, the deleted
top row of `bf(y | thres(3) ~ x)`:

| arm | integer | ordered factor | character |
|---|---|---|---|
| base | 0 NA, row filled | 4 of 4 cells NA | identical to integer |
| lane | 0 NA | identical to integer, bitwise | identical to integer |

### Finding 4. The defect-B census

Corrected in place, under "B. ... The census of multi-index entries,
corrected". The NEWS bullet now describes the rule (any prior entry
covering more than one parameter, naming `Intercept` on an ordinal
model, `cor`, `ar`, `ma`, `cortime` and `rescor`) instead of the ordinal
case alone, and so do the comments on `draw_prior_entry()` and
`prior_entry_label()`.

### Nits

- `?frm_bootstrap` no longer pins a number. It said "about 2e-11
  relative"; the reviewer measured 8.98e-12 to 6.09e-11 over four
  families, grouped and ungrouped, so the sentence now says "to
  optimizer tolerance" and the range lives here and in the review.
- `prior_entry_label()` kept `length(e$idx) != 1L`, deliberately, so a
  zero-length index cannot reach the vector comparison below it either;
  it now returns `comp[]` for that case rather than `comp[NA:NA]`. The
  comment says which classes reach the branch.
- The `accepts_aterms` half of the `thres_pin_of_fit()` gate now HAS a
  construction that makes it false, in `test-thres-refit.R`: an ordinal
  `frmtmb_family()` from outside the package that declares only
  `weights` gives `NULL`, and the same family declaring `thres` gives
  `nthres = 3` with the fitted labels. Keeping the half is right, because
  pinning writes `thres` into the addition-term values and a family that
  does not accept it would then be refused by name.

### Recorded as an observation, not changed

`simulate()` followed by `frm()`, written out by hand, fits K-1
thresholds: 2 instead of 3 ungrouped and 6 instead of 7 grouped, at the
same maximized log-likelihood to 7.4e-11 to 1.3e-10 relative
(the reviewer's `dev/thresrefit-rev-01-boot.R`). That is documented
behavior of `frm()`, which counts the thresholds of whatever data it is
given, and `thres(K)` pins it. `?frm_bootstrap` and the `?frm` section
now say so, because a reader of either will meet it.

### What the controls say, so the new refusals do not fail closed

`dev/thresrefit-p1-controls.R`, log `-after-log.txt`.

- `influence(groups = "id")` on `bf(y | thres(gr = g) ~ x + (1 | id))`,
  seed 1201, n = 120, where `id` is NOT the threshold grouping factor:
  8 x 8 table, 0 NA cells, 0 warnings, `cooks.distance()` with no NA, on
  both arms. So the "cannot drop a level" refusal does not touch the
  ordinary case; the smallest number of rows any `g` level keeps across
  the eight `id` deletions is 35.
- The pin stays inert: `influence(force = TRUE)` on cumulative, sratio,
  cratio and acat with nothing to lose, and on gaussian, bernoulli and
  poisson, gives 0 NA cells and 0 warnings, and
  `thres_pin_of_fit()` is `NULL` on all three non-ordinal fits.
- Nothing wrong warns about nothing: the control at seed 1101 without a
  doctored `data` reports 0 warnings and 0 NA cells.

### Test counts for this round

`test-thres-refit.R` grew from 71 to 118 expectations in 14 blocks.

    test-thres-refit.R (lane)            pass=118 fail=0  err=0 skip=0
    test-thres-refit.R (FRMTMB_LIB=base) pass=68  fail=23 err=5 skip=0
    test-data2.R (lane)                  pass=31  fail=0  err=0 skip=0
    test-data2.R (FRMTMB_LIB=base)       pass=29  fail=2  err=0 skip=0

The new pins for blocking 1 and blocking 2, all SEEN TO FAIL on `base`
and all behavioral (`dev/thresrefit-seen-failing-base.txt`):

    an emptied INTERIOR category does not slide the row over
      Expected `anyNA(inf$fixed)` to be FALSE
      Expected `anyNA(cooks.distance(inf))` to be FALSE
      Expected `got$ordered` to be identical to `got$integer`
    an emptied category does not slide a grouped row over
      Expected `anyNA(inf$fixed)` to be FALSE
      Expected `got$ordered` to be identical to `got$integer`
    a hand-written thres(K) is refit on every response coding
      Expected `anyNA(inf$fixed)` to be FALSE
      Expected `got$ordered` to be identical to `got$integer`
      Expected `anyNA(got$ordered[itop, ])` to be FALSE
    influence() counts the refits that failed and says so
      <subscriptOutOfBoundsError> in `[<-`(*tmp*, i, names(fe_i), ...)
    a refit cannot add or drop a thres(gr = ) level
      <subscriptOutOfBoundsError> in `[<-`(*tmp*, i, names(fe_i), ...)
    data2 outlives its calling environment across saveRDS (test-data2.R)
      Expected `suppressWarnings(influence(r_env, groups = "g"))` to
        throw a error                                            (x2)
## NEWS entry

See `NEWS.md`, `# frmtmb (development version)`.
