# conditional_effects fixes: lane ledger (`wt-ce`)

Worktree `C:\Users\adf44\source\r\frmtmb-wt-ce`, branch `wt-ce`, from
`2210aa1` (frmtmb 0.51.0). Private library
`...\scratchpad\ce-lib`, private Stan cache `dev/stan-cache` (copied
from the main checkout, 53 programs plus `makevars-cxx17.mk`).

Written incrementally: the scratchpad is wiped without warning.

## The contract, and how the brief numbers it

The brief's ten defects are numbered by the REVIEW's 22-row table
(`dev/reviews/2026-09-05-brms-methods.md` section 2), while the
`Class:` verdicts live in `dev/brms-methods-tests.md` under a different
numbering. The map, once, so nothing below is ambiguous:

| brief item | review row | doc finding | what |
| --- | --- | --- | --- |
| 1 | 1 | 1b | ce plots `exp(eta)`, not the expected response |
| 2 | 2 | 1c | mixture `theta` response scale is the predictor |
| 3 | 3 | 1d | ce refuses a mixture with `theta ~ x` |
| 4 | 4 | 2 | `signif(mean +/- sd, 3)` as the held VALUE |
| 5 | 5 | 3 | `int_conditions =` accepted and ignored |
| 6 | 6 | 4 | frame lacks `y`, `z`, `cond__`, `effect1__` |
| 7 | 8 | 6d | `re_formula = NULL` takes the first observed level |
| 8 | 9 | 6b | ordinal layout, ignored `categorical =`, effect key |
| 9 | 10 | 17 | `mo()` gets a 100-point grid |
| 10 | 18, 21 | 11b, 15 | `ranef()`/`coef()` keys; `se()`'s unused sigma |

**One contradiction in the brief, recorded rather than resolved
silently.** Item 3 is spelled out as a defect to fix ("brms returns the
x panel"), and the same review row 3 appears in the brief's list of
design choices that are NOT the contract. Row 3 is class **C** in both
the review and the doc. I fixed it: it is the enabling half of item 2
(without it a mixture with `theta ~ x` has no default panel at all to
compare), and the fix is a fallback that changes nothing for any model
whose selected dpar already has something to plot.

## Exclusion table, before

`brms_exclusions()` at `tests/testthat/helper-brms-methods.R:551`,
32 rows, 11 of class D:

| list | key | finding | class |
| --- | --- | --- | --- |
| brms_ce_shapes | r2 | 17 | D |
| brms_ce_shapes | r12a..r12e | 6b | D (5 rows) |
| brms_ce_shapes | r13 | 6b | D |
| brms_ce_shapes | r16, rC16 | 1b | D (2 rows) |
| brms_dpars_of | r17:theta1 | 1c | D |
| brms_dpars_of | r14c:sigma | 15 | D |
| brms_ce_shapes | r5 | 18 | P |
| brms_ce_shapes | r17 | 1d | C |
| brms_linpred_shapes | r12e, r13, r17 | 13 | C (3 rows) |
| brms_meanlink_shapes | 10 keys | 14 | P (10 rows) |
| brms_resid_shapes | 6 keys | omissions | P (6 rows) |

## Measurements

(filled in as they are taken)

### Baseline, from the installed core alone (no brms, no Stan)

`scratchpad/ce-before.R`, worktree core at 2210aa1 in `ce-lib`:

| defect | measured before |
| --- | --- |
| 1 zero-inflated ce | `ce[1] = 0.68557372533` vs `predict(response)[1] = 0.54967822329`, ratio 1.2472 at grid point 1 and 1.9538 at 100; `ce == predict(type = "conditional")` to 0 |
| 1 hurdle ce | ratio 0.86193188603 at point 1, 3.37134980095 at point 100 (the doc's numbers, reproduced) |
| 4 signif | held at `-1.0500, -0.0107, 1.0300` vs exact `-1.0545225911165, -0.0106839128603, 1.0331547653958`, max abs 4.52e-3 |
| 5 int_conditions | `ignoring unknown arguments to predict(): int_conditions` |
| 6 columns | `x, estimate__, se__, lower__, upper__` |
| 7 re_formula NULL | Subject 308 to 5.68e-14, 86.5480358829 away from the population curve, no `Subject` column |
| 8 ordinal | keyed `"x"`, 300 rows, `categorical =` a no-op |
| 9 mo() | 100 rows, grid `0, 0.0303, 0.0606, ...` |
| 2 mixture theta | `predict(response, theta1)` `identical()` to the link scale, range [-1.2197, 2.7993], 37.5% of rows outside [0, 1], `max|response - plogis(eta)| = 1.857` |
| 3 mixture ce | `No plottable predictors found for dpar 'mu1'` |
| 10 ranef/coef | `names(ranef) = "Days | Subject"`, `names(coef) = "Subject"`, `ranef(fit)$Subject` NULL |
| 10 se() sigma | `predict(response, dpar = "sigma") = 1`, link scale 0 (and `sigma(fit)` already 0) |

### After, same script

Every row above now reads: ZIP ratio 1 and 1 (`ce == predict(type = "response")` exactly), hurdle
1 and 1, held values exact to 0, no warning, columns
`x, y, cond__, effect1__, estimate__, se__, lower__, upper__`,
`re_formula = NULL` equal to the population curve with a `Subject`
column of NA, ordinal keyed `x:cats__`, `mo()` 4 rows at 0 1 2 3,
theta1 in [0.228, 0.9426] and `max|response - plogis(eta)| = 2.22e-16`,
the mixture default drawing the `x` panel.

### The Wald band for the expected response, measured against a bootstrap

`bf(y ~ x, zi ~ x) + zero_inflated_poisson()`, 300 rows, 11 grid points,
`band = "boot"` with 200 refits (`seed = 4`):

- estimates identical between the two bands (max abs 0).
- 95% width, Wald `0.35914, 0.33004, 2.50673` at points 1, 6, 11;
  bootstrap `0.36938, 0.33034, 2.78366`. Max relative difference in
  width over the grid **0.1278** - which is ONE DRAW, not a bound. The
  review swept it: median 0.125 over eight seeds, range 0.063 to 0.198,
  falling to 0.041 at 800 refits. The honest statement is that the gap
  is the percentile bootstrap's own Monte Carlo error at 200 refits and
  converges away, except at the ENDS of a steep shape, where the Wald
  band stays about 28% wider at 3000 refits because it is symmetric on
  its own scale and the percentile band is not. The help page's band
  section now says this.
- every Wald lower bound is strictly positive (the band is symmetric on
  the mu link's scale, `log` here, so it cannot cross zero).
- on a family whose mean IS `linkinv(eta_mu)` the new route is not
  taken at all, and the gaussian band is bit-identical to the old one
  (max abs difference 0).

## What landed

| # | fix | site |
| --- | --- | --- |
| 1 | expected response + joint delta band on the default display | `R/conditional-effects.R` `mean_display`/`pred_dpar`, `ce_band_scale()`. The joint delta method itself, `predict_mean_se()`, is PRE-EXISTING at `2210aa1` and untouched; what is new is the route to it from the default display |
| 2 | mixture mixing weights report the softmax | `R/families.R` `post$dpar_response` in `mixture()` and `mixture_mvn()`; `R/predict.R` `dpar_report_hook()`, `dpars_natural()` |
| 3 | the default effect search falls back across dpars | `ce_plot_vars_any()`, used in `ce_grids_build()` |
| 4 | the moderator is held at the exact `mean +/- sd` | `ce_second_values()`, label in `ce_effect2()` |
| 5 | `int_conditions =` implemented | new formal, `ce_int_cond()`, `ce_grid_values()` |
| 6 | brms's columns | `ce_frame()` |
| 7 | `re_formula = NULL` is a new level | `ce_group_vars()`, `na_vars` in `ce_build_nd()` |
| 8 | `categorical =` honored, `"x:cats__"` key, dots checked | `ce_display_kind()`, `ce_finalize(cats_key =)`, `ce_dots()`, `ord_prob_se(weights =)` |
| 9 | `mo()` steps by one | `ce_step_vars()`, `ce_grid_values(stepwise =)` |
| 10 | `ranef()` keyed by the factor; `se()`'s sigma reports 0 | `R/methods-fit.R`; `R/predict.R` `se_unused_sigma()` |
| C | brms's two-way row order, brms's `method =` names | `ce_build_nd()`, `ce_profile_eta_ci()`, `ce_method()` |
| C | `trials()` held at 1 on the grid, with a message | `ce_trial_vars()` in `ce_grids_build()` |

## Exclusion table, after

**22 rows: 0 defect, 18 paradigm, 4 design choice** (was 32: 11/17/4).

| list | key | finding | class |
| --- | --- | --- | --- |
| brms_ce_shapes | r5 | 18 | P |
| brms_ce_shapes | r12e | 13 | C |
| brms_ce_shapes | r13 | 9 | P |
| brms_linpred_shapes | r12e, r13, r17 | 13 | C |
| brms_meanlink_shapes | 10 keys | 14 | P |
| brms_resid_shapes | 6 keys | omissions | P |

Dropped: `r2` (17), `r12a`-`r12d` (6b), `r16`, `rC16` (1b), `r17` (1d),
`r17:theta1` (1c), `r14c:sigma` (15) - eleven defect rows, and `r13`
came back as a paradigm row for a different reason (see below).

Two rows are NOT defects and are new:

- `r12e` under finding 13: a `cs()` term stores its evaluated column
  and the label `"csz"`, not the variable, so no display can enumerate
  a `cs()` predictor and brms's effect list carries a `"z:cats__"`
  panel frmtmb's does not.
- `r13` under finding 9: a NOMINAL family's category probabilities have
  no thresholds, so the ordinal delta method does not apply and the
  wald band is refused by name. Its own block covers the shape under
  `band = "boot"`, where the estimate is the fit's own and agrees with
  brms elementwise.

The one-way loop and the guard's probe now pass `categorical = TRUE` to
BOTH packages for a polytomous family (`brms_ce_args()`): brms refuses
its own default for a nominal family and warns against it for an
ordinal one, so default-against-default compared a refusal with a
curve.

## The three unpinned seams

- `R/conditional-effects.R`, `method = "predict"` on a family that
  draws the response whole: PINNED, and repaired first. The two family
  checks were raised INSIDE the grid loop, after the point estimate had
  been computed, so on `mixture_mvn()` the call died on
  `replacement has 200 rows, data has 100` instead of reaching its own
  message. They are hoisted out of the loop. `tests/testthat/test-unpinned-seams.R`
  asserts `simulate()` accepts the fit, `method = "predict"` refuses by
  name, the expected-response display refuses for the other reason (the
  mean is an n by D matrix), and `dpar = "theta1"` still draws.
- `R/par-template.R:155` and `R/priors.R:719`
  (`carry_finalized_responses()`): NOT pinnable from outside, and the
  measurement is in the test file rather than a claim. I installed a
  copy of the package with BOTH lines removed
  (`scratchpad/ce-src2` -> `ce-lib2`) and ran `get_prior()` and
  `par_template()` against a family whose `family_finalize()` replaces
  a link, and again against one that swaps which dpar is primary:
  **identical output, both libraries, every call**. The reason is
  structural: the response spec's own `primary_dpars` is not refreshed
  by finalize (the frame refreshes `resp$dpars[[i]]$link` only), and
  neither function reads a dpar's link. What the test pins instead is
  the helper's own contract - that it carries the finalized links and
  family, and passes a family without `family_finalize` through
  untouched - plus that both call sites still answer on such a family.

## Verification

| what | result |
| --- | --- |
| `test-brms-methods.R`, warm, gated | 47 tests, **950** assertions, 0 fail, 0 skip, 154.6 s (was 46 / 872) |
| `test-brms-likelihood.R`, warm, gated | 32 tests, **352** assertions, 0 fail (351 before, plus the one `attr(ranef(fit)[[1]], "term")` assertion this round added) |

## The design choices left open, with a cost each

Adopted this round because each was a few lines and the divergence was
a porting stumble with nothing to recommend it: review row 3 (doc 1d,
the mixture's covariate), row 7 (doc 5, the two-way row order), row 12
(doc 6c, brms's `method =` names), and the `trials()` grid rule brms
documents and messages about.

Left open, for the user to decide:

| finding | choice | cost of adopting brms's side |
| --- | --- | --- |
| review 14 / doc 7 | `hypothesis()` returns the flat `frmtmb_hypothesis` table, not a `brmshypothesis` list | a `$hypothesis` alias is about five lines, but it COLLIDES with the table's own `hypothesis` column, so `h$hypothesis` would have to mean two things; the brms idiom fails loudly today |
| review 19 / doc 12 | `coef()`'s return type is a function of the model (vector, list, or list of data frames) | always returning a list breaks the `glm`-style `coef(fit)` vector that the `stats::coef()` fallback exists to give; the honest fix is one paragraph on the help page, not a change |
| review 22 / doc 13 | `cs()` columns are unreachable through `dpar =` | needs the cs() term to store its VARIABLE (frame.R, a sibling's file) plus a dpar spelling; it would also let `conditional_effects()` enumerate a `cs()` predictor and put `r12e` back in the agreement loop, which is the only exclusion row that a repair could still remove |
| doc 6 | `conditional_effects(dpar =)` enumerates that dpar's own predictors, brms enumerates the model's | matching brms means drawing flat panels for predictors the dpar does not contain; frmtmb's is the better display and the cost is that a loop over the returned list runs a different number of times |
| doc 6e | `predict()` spells it `re.form`, `conditional_effects()` `re_formula` | accepting `re_formula` as an alias on `predict()` is a few lines; today it warns by name, which is right for a direct call and is no longer misattributed from `conditional_effects()` |

## Not done, and why

- `frmtmb.sample`'s `conditional_effects.frmtmb_draws()`
  (`extensions/frmtmb.sample/R/conditional-effects-draws.R`) still
  builds the OLD column set: it shares `ce_grids_build()`,
  `ce_boot_one()` and `ce_finalize()` (whose signatures only gained
  arguments with defaults, so it is unaffected) but assembles its own
  per-effect frame with `g$nd[g$ev]`. Giving it brms's columns means
  calling `ce_frame()`, which would have to join the exported
  extension API in `R/sampling-api.R` - a sibling lane's file this
  round. Its grid row order follows `ce_build_nd()` and so did change
  with the rest.
- `predict_mean_se()`'s refusal for a family whose mean is not one
  number per observation (`mixture_mvn()`) surfaces from `predict()`
  rather than naming `conditional_effects()`. The message names the
  family and the reason, and the test pins it; a display-level refusal
  would be a second message saying the same thing.
| full core suite, one file per process, by name | **110 of 110 files ran, 1097 tests, 6168 assertions, 0 failed, 0 errors**, 88 skips (the three brms-gated tiers, run separately with the gate: 45 + 29 + 11, plus 2 in `test-brms-agreement.R` and 1 in `test-fuzz.R`) |
| `R CMD check --as-cran` | **Status: OK** - no ERROR, no WARNING, no NOTE. Tests OK in 547 s, vignettes rebuilt. `_R_CHECK_CRAN_INCOMING_=false`, pandoc from RStudio's quarto tools on PATH |
| lint | `R/conditional-effects.R` 0 lints, the same as its baseline; the pre-existing counts in `R/predict.R` (9), `R/families.R` (6) and `R/methods-fit.R` (2) are all outside the edited hunks |
| roxygen | `roxygenise()` rewrote `man/conditional_effects.Rd` and `man/ranef.Rd`; a second run wrote nothing (idempotent) |

The core suite ran as one sequential stream plus three parallel ones
(the box carries five lanes); 28 files ran twice as a result, which the
audit deduplicates by name. No file is missing and none failed.

# Punch round, 2026-09-05

Against `dev/review-ce.md` (GO-WITH-FIXES, verdict at `:1003`). P1 and
P2 are the blockers, P3-P5 the cheap fixes, P6-P9 the follow-ups the
review listed. Every number below was measured here after the fix,
against the same shapes the review used where it named one.

## P1 (blocking): `band = "boot"` dropped the group variance

**What it was.** `re_formula = NULL` sets the grid's grouping column to
`NA` (`R/conditional-effects.R:1550`), so every bootstrap refit
predicted a NEW level's ZERO modes and the percentile band carried
coefficient uncertainty only. The interval was the population interval
under another name, and nothing on the returned object said so.

**The fix.** The wald band adds an unseen level's marginal variance
(`lp_extra_var()`); a bootstrap has replicates instead of a variance,
so it now DRAWS that level's effects once per replicate. The mechanism
is a placeholder level: the boot grid carries an observed level so the
design maps it (`ce_boot_grids()`), and each replicate overwrites that
level's coefficients (`ce_draw_new_levels()`), which is exactly
`z_i' u` through the ordinary `Z`. The covariance drawn from is
`lp_extra_var()`'s own - `covstruct_registry[[...]]$vcov()` with the
Student-t factor - and the blocks it skips (`gr_cov`, `gr_prec`,
`car`, `spde`) are ZEROED here for the same reason, so the two bands
make the same assumption about the same models.

New helpers in `R/conditional-effects.R`: `ce_new_level_spec()`,
`ce_boot_grids()`, `ce_draw_new_levels()`; wired into the
`band == "boot"` call site and into `ce_boot_draws()`'s `FUN`.

sleepstudy, `Reaction ~ Days + (Days | Subject)`, resolution 6,
`boot = 100`, `seed = 11`, band WIDTH:

| band | grid points 1..6 |
| --- | --- |
| wald, population | 26.00 26.69 31.20 38.21 46.59 55.73 |
| wald, new group | 96.78 107.89 131.89 163.18 198.34 235.65 |
| boot, population | 24.61 26.81 31.32 38.33 46.55 55.89 |
| boot, new group, BEFORE | the population band, ratio 1.00 throughout |
| boot, new group, AFTER | 92.69 95.39 122.07 159.63 194.38 220.26 |

`boot new / boot population` **3.77 3.56 3.90 4.16 4.18 3.94** (was
1.00 throughout); `boot new / wald new` 0.96 0.88 0.93 0.98 0.98 0.93.
The point estimate is unchanged (max abs 0 against the population
curve) and `Subject` is still all `NA`.

**The ordinal case, which is why this blocked.** `y ~ x + (1 | g)`,
`cumulative()`, 25 groups, n = 400: the per-category wald band is
refused BY NAME with `re_formula` and the message sends the user to
`band = "boot"`, so boot is the only band available there.

| | before | after |
| --- | --- | --- |
| estimate, new vs population | max abs 0 | max abs 0 (unchanged, correctly) |
| WIDTH, new vs population | **max abs 0**, bit-identical | ratio 1.80 to 4.57, wider at every one of the 15 rows |

Pinned by "band = 'boot' carries a new group's variance, ordinal
included" in `tests/testthat/test-ce-bands.R`, which asserts the
gaussian and the ordinal case, the identical estimate, the `NA`
grouping column, and that the band is not the population band.

## P2 (blocking): the `ranef()` re-key broke `frmtmb.sample`

Reproduced first: **5 failures**, `test-sample-direct.R:67` and
`test-sampling-ported.R:560,561,562,566`, all the block key.

Fixed on the SAMPLE side, because core's re-key is the deliberate
brms/lme4 parity the review confirms: the five assertions now read the
factor key and the `"term"` attribute. One code fix went with them, and
it is a defect the re-key created rather than a test edit:
`ranef.frmtmb_draws()` in
`extensions/frmtmb.sample/R/methods-draws.R` assembled its arrays by
NAME, so for a model with two blocks on one factor - which now share a
name - `out[[tn]] <- st` wrote the same name twice and the SECOND block
was DROPPED: the list came back with one entry instead of two. (An
earlier draft of this ledger said the first block's draws were returned
twice; the review reproduced the assembly against both cores and it is
the drop, which is worse.) It indexes by position now, as core's own
`print()` and `as.data.frame()` do, and carries the `"term"` attribute
through.

`extensions/frmtmb.sample/NEWS.md` gains a BEHAVIOR CHANGE bullet: its
own `ranef()` output is re-keyed, `[["1 | g"]]` returns `NULL`, and the
two-blocks-on-one-factor assembly changed.

| suite | before | after |
| --- | --- | --- |
| `frmtmb.sample`, 10 files, one process each | 136 tests, 881 assertions, **5 failures** | 136 tests, **888 assertions, 0 failures**, 2 skips |

888 = 881 + the 5 that now pass + 2 new `"term"` assertions.

## P3: `effect2__` collapsed distinct moderator values

`ce_effect2()` took its levels from `round(v, 2)`, so two values that
rounded together became ONE level - `plot()` groups on that column, so
two curves were drawn as one series - and the user's names were dropped
with them because the name-count guard stopped matching. The levels now
come from the distinct VALUES and the rounding is the label only;
labels widen past two decimals only as far as they must to stay
distinct, because assigning a duplicated label to `levels<-` merges the
levels again.

| call | before | after |
| --- | --- | --- |
| `int_conditions = list(z = c(lo = 0.001, hi = 0.002))` | 1 level, `"0"`, names dropped | **2 levels, `"hi" "lo"`** |
| same, unnamed | 1 level | **2 levels, `"0.002" "0.001"`** |
| default `mean +/- sd` | `1.05 0.04 -0.97` | **unchanged** |

Pinned by "effect2__ keeps distinct moderator values distinct" in
`tests/testthat/test-ce-bands.R`.

## P4: the NEWS family list

Enumerated through `mean_is_mu()` here rather than copied: **14**
user-facing families move (`asym_laplace`, `beta_binomial`, `binomial`,
`cox`, `hurdle_gamma`, `hurdle_lognormal`, `hurdle_poisson`,
`lognormal`, `shifted_lognormal`, `zero_inflated_asym_laplace`,
`zero_inflated_beta`, `zero_inflated_binomial`,
`zero_inflated_negbinomial`, `zero_inflated_poisson`), plus any
`trunc()` response and any `mixture()`; 23 do not. The old bullet named
only the zero-inflated, hurdle, `trials()`, truncated and mixture
cases and then said "every other family is unchanged, bit for bit",
which told `lognormal`, `shifted_lognormal`, `asym_laplace` and `cox`
users their plots had not moved. They had.

**The lognormal case, measured here.** `y ~ x`, `lognormal()`, response
range 0.077 to 7.809:

```
old curve (the mu dpar on its identity link) : -1.4886 -0.7529 -0.0171  0.7186  1.4543
new curve (the mean, exp(mu + sigma^2/2))    :  0.2859  0.5966  1.2452  2.5987  5.4234
new == predict(type = "response")            : max abs 0
```

The old curve was NEGATIVE over most of the grid for a strictly
positive response. NEWS now names all fourteen and says so.

## P5: the `ranef()` NEWS bullet now says what BREAKS

It listed what is kept (the `"term"` attribute, the `grp` column,
`VarCorr()`'s key - all three verified by the review). It now leads
with the migration hazard: `ranef(fit)[["Days | Subject"]]` returns
`NULL`, and a `NULL` flows into arithmetic as `numeric(0)` rather than
erroring. It also points at `frmtmb.sample`, whose draws surface is
re-keyed with it.

## P6: the mixing weight's SE is one-predictor, documented and pinned

Measured here on a three-component fit (`theta1 ~ x, theta2 ~ x`,
n = 600), against the analytic joint delta method over BOTH theta
predictors' coefficients using `vcov(fit)`:

```
reported vs the one-predictor rule : max rel 6.66e-16   (it IS that rule)
reported vs the JOINT delta, theta1: 6.6% to 25.9% WIDER
reported vs the JOINT delta, theta2: 5.5% to 26.1% WIDER
wider at every row                 : TRUE
K = 2, reported vs exact           : max rel < 1e-10
```

The review measured 8.9-17.3% on its own fit; the range is data
dependent and the DIRECTION is not. Documented in `?mixture` and
`?predict.frmtmb_fit` and pinned by "a mixing weight's response-scale
SE is the one-predictor rule" in `tests/testthat/test-unpinned-seams.R`,
which asserts it IS the one-predictor rule, that it is never narrower
than the joint one, and that K = 2 is exact.

## P7: the band help no longer implies the three bands agree

The band section now says what they do: wald and bootstrap agree in the
middle of a grid and the residual there is the bootstrap's own Monte
Carlo error (13% at 200 refits, 4% at 800 on a zero-inflated fit),
while at the ENDS of a steep shape they need not converge at all (the
wald band stays about 28% wider at 3000 refits, being symmetric on its
own scale where the percentile band is not).

## P8: the `frmtmb.sample` draws method, recorded not fixed

`extensions/frmtmb.sample/R/conditional-effects-draws.R:41,45,81,87,91,98,102`.
It shares `ce_grids_build()`, so it SILENTLY INHERITED four grid
changes (brms's row order, exact `mean +/- sd`, `mo()` stepwise grids,
`trials()` held at 1 with its new `message()`) while keeping the old
frame shape, the old `cond__` rule, the old ordinal key, and the old
`re_formula` / `categorical` semantics. Migrating it needs `ce_frame()`
in the exported extension API, which lives in `R/sampling-api.R`, a
sibling lane's file this round; the coordinator's reach for the
extension was the `ranef()` fix. Both halves are now stated in that
package's NEWS so no user meets them silently.

## P9: three claim-wording corrections

- `predict_mean_se()` is PRE-EXISTING at `2210aa1` and untouched; only
  the route to it from the default display is new. The "What landed"
  table says so.
- "within 0.128 of a 200-refit bootstrap" was one draw, not a bound.
  The measurement section now gives the review's sweep (median 0.125
  over eight seeds, 0.063 to 0.198, 0.041 at 800 refits) and says the
  gap is the bootstrap's Monte Carlo error except at the ends of a
  steep shape.
- `tests/testthat/test-unpinned-seams.R`: a `family_finalize()` cannot
  swap which dpar is PRIMARY (`primary_dpars` is never refreshed); the
  available variant is reordering `family$dpars`, which is what the
  review ran.

## Not adopted from the review

- Section 1b's tail gap and section 12's `R CMD check` residuals
  (LaTeX and V8 missing on the reviewer's box) are recorded, not
  "fixed": the first is a documented property of the two band methods,
  the second is an environment gap. My own as-cran run on this box was
  `Status: OK`.
- The review's correction that main is at `9b9011a`, not `2210aa1`, is
  right and changes nothing: the three commits on top touch
  `_pkgdown.yml`, `docs/**` and a spline test only.

## Punch-round verification

| what | result |
| --- | --- |
| `test-brms-methods.R`, gated, warm | 47 tests, **950** assertions, 0 fail, 0 skip, 121.2 s |
| `test-brms-likelihood.R`, gated, warm | 32 tests, **352** assertions, 0 fail |
| `test-ce-bands.R` | 21 tests, **159** assertions, 0 fail (was 19 / 139: +2 tests, P1 and P3) |
| `test-ce-facets.R` | 6 / 92 / 0 |
| `test-effects.R` | 9 / 54 / 0 |
| `test-ordinal-fitted.R` | 13 / 87 / 0 |
| `test-mvn-mixture.R` | 9 / 277 / 0 |
| `test-unpinned-seams.R` | 4 tests, **28** assertions, 0 fail (was 3 / 23: +1 test, P6) |
| `test-message-uniqueness.R` | 1 / 6 / 0 |
| `test-boot.R`, `test-simulate-ergonomics.R`, `test-v15.R`, `test-methods.R`, `test-bracket-access.R` | 43, 44, 41, 27, 8 assertions, 0 fail |
| **`frmtmb.sample`, all 10 files** | 136 tests, **888** assertions, **0 fail**, 0 error, 2 skips |
| roxygen | 3 pages rewritten, second run silent (idempotent) |
| full core suite after the punch round, one file per process, by name | **110 of 110 files, 1100 tests, 6193 assertions, 0 failed, 0 errors**, 88 skips (the gated tiers). Up from 1097 / 6168 by the three blocks this round added. `test-perf.R` failed once under three-way contention on its wall-clock envelope (`t_large < 100 * max(t_small, 0.01)`, 2.6 vs 1.0) and passes alone; it touches nothing this lane changed |
| `R CMD check --as-cran`, re-run on the final code | **Status: OK** - 0 ERROR, 0 WARNING, 0 NOTE. `checking tests [548s] OK`, examples OK, `--run-donttest` OK, vignettes rebuilt [394s] |
| lint | `R/conditional-effects.R` 0, unchanged from its baseline |

## Re-check residuals, same punch round

Against `dev/review-ce.md` from `:1075` (the re-check; updated verdict
at `:1360`). R1 was blocking-class, R2 to R4 are one-liners.

### R1 (blocking): `boot =` reuse could not tell a population bootstrap
from a new-group one

**Why the key could not see it.** `ce_boot_key()` hashed the grids'
`nd` columns, and the two grids are BYTE IDENTICAL: `ce_ref_value()`
holds an unvaried factor at `levels(col)[1L]` and the new-level
placeholder is `bk[["levels"]][1L]`, which on sleepstudy are both
`"308"`. Two independent choices happening to coincide, so a population
call and a `re_formula = NULL` call keyed the same. The reuse path is
one the help page recommends, and the failure reopened P1 exactly:

```
population band width : 22.65 24.29 25.77 36.95 48.12
new-level  band width : 83.40 101.60 121.43 156.04 221.99
```

**The fix.** `ce_new_level_key()` (new, in `R/conditional-effects.R`)
reduces the new-level spec to what identifies the DRAWS (which blocks
are redrawn, where they are written, whether they are drawn or zeroed)
and `ce_boot_key()` takes it as a sixth argument. `ce_boot_draws()`
keeps it beside the key on the returned object (`bs$ce_new`) so a
mismatch can say WHICH of the two it is, and the reason is checked
BEFORE the grid reason, since the grids can be identical.

Both directions are now refused by name, with the reason:

| reuse | before | after |
| --- | --- | --- |
| population draws in a `re_formula = NULL` call | ACCEPTED, returned the POPULATION band (22.78 to 59.89 where the correct band is 92.51 to 269.33) | REFUSED: "predictions for a different group: this call conditions on a NEW group (re_formula = NULL), and those draws do not carry that group's effects, so their percentiles are the population band" |
| new-level draws in a population call | ACCEPTED, band four times too wide | REFUSED: "... those draws carry a NEW group's effects ... so their percentiles are too wide for it" |
| population draws in a population call | worked | works, reproduces its own band to **max abs 0** |
| new-level draws in a new-level call | worked | works, reproduces its own band to **max abs 0** |
| a genuinely different grid | refused | still refused, and still with the GRID reason |

Pinned by "a bootstrap cannot be reused across the new-group boundary"
in `tests/testthat/test-ce-bands.R`, which asserts both refusals by
message, both matching reuses reproducing their band, and that the
grid-mismatch reason is still reachable and distinct.

### R2: a Student-t block is drawn gaussian, deliberately

One line added to the band section of `?conditional_effects`: a
Student-t random-effect block enters EITHER band as a gaussian with the
t variance (`nu / (nu - 2)` times the scale matrix), because the delta
method has no other shape to offer and the bootstrap draws the same way
so the two bands stay comparable. It is the right variance around a
heavier-tailed truth, not the right quantile. `R/predict.R:1538-1541`
already carried the honest comment on the wald side; the choice is now
recorded where a user reads it.

### R3: the old by-name assembly DROPPED the second block

Corrected in both places. `out[[tn]] <- st` writes one name twice, so
for two blocks on one factor the second was dropped entirely and the
list came back with ONE entry, not with the first block's draws twice.
The review reproduced it against both cores (main by-name: 2 entries,
correct; lane by-name: 1 entry; lane by-position: 2 entries, correct).
Worse than this ledger and the sample NEWS first described; the
by-position fix handles it either way.

### R4: the P5 NEWS bullet leads with the breakage

`NEWS.md` now opens the bullet with **"what breaks is
`ranef(fit)[["Days | Subject"]]`, which now returns `NULL`"** and the
`numeric(0)` hazard, then gives the rationale and what is kept.

While there: the prose spaced hyphens I had introduced in this
section, and in the roxygen this round added, are gone; the house
style keeps `-` for arithmetic.

### Runs after the residuals

| what | result |
| --- | --- |
| `test-ce-bands.R` | **22 tests, 167 assertions, 0 fail** (was 21 / 159: +1 test for R1) |
| `test-brms-methods.R`, gated, warm | **47 tests, 950 assertions, 0 fail, 0 skip**, 145.6 s |
| `test-message-uniqueness.R` | 1 / 6 / 0 |
| `test-boot.R` | 7 / 43 / 0 |
| `test-unpinned-seams.R` | 4 / 28 / 0 |
| roxygen | 3 pages rewritten, second run silent |
| lint | `R/conditional-effects.R` 0 |
