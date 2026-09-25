# Lane wt-simnewdata: findings

Lane of 2026-09-23/24, off main `51bfaa4e` (frmtmb 0.62.0, frmtmb.sample
0.10.0). Base build: `rellib-r3`, read only. Lane build:
`C:/Users/adf44/source/r/simnewdata-lib`. Every number below names the
script that produced it and the log it is in, under
`dev/simnewdata-log/`.

## Answer first

All four items are done. Punch round 1 (see its section) made the
newdata draw about 6 times cheaper per replicate with identical draws,
made a missing covariate draw NA as `predict()` does, kept
`frm_bootstrap()` bitwise at 0.62.0 by the user's decision, restored
the released NEWS section, and corrected this record on the residual
correlation: frmtmb and brms differ at seen times too.

1. `simulate()` takes `newdata` and `allow_new_levels`, and
   `pp_check(newdata = )` on a fit answers instead of refusing. The two
   ported rows `brmsfit-methods:675` and `:682` pass. Bin 1 went from
   250 to 252 of 494, and the ledger stopped on exactly those two stale
   verdicts and on nothing else.
2. `simulate()` reads `re_formula` through the same resolution as
   `predict()` (`re_keep_plan()`, split out of `re_resolve()`). `~1` and
   `~0` now draw what `NA` draws. A partial formula keeps the named
   terms and redraws the rest. `~(1 | nosuch)` and a partial formula
   beside a factor smooth are refused in `predict()`'s words.
3. The multinomial `error_binned` refusal says "a set of counts over
   categories", on the fit and on draws.
4. Neither frmtmb.sample test file leaves `Rplots.pdf` now: measured by
   deleting the file and running each file in its own process.

**Two things the brief did not anticipate, both found by measuring:**

- **brms fixture 1 is a structured draw.** Its `arma(visit, patient,
  cov = TRUE)` term makes `sim_is_structured()` TRUE
  (`dev/simnewdata-probe1.R`, `probe1-base.txt`). The feasibility probe
  was run on a plain gaussian, so it never met this. If `simulate()`
  refused every structured draw at newdata, as `predict()` and the
  draws method do, `:675` and `:682` would stay defects. So the lane
  rebuilds a residual correlation term on the newdata rows
  (`autocor_for_newdata()`). A structured FAMILY is still refused.
- **`simulate(re_formula = NA)` redrew population smooths.** Every
  release through 0.62.0 redrew the penalized coefficients of `s()`,
  `gp()` and `hsgp()` from the smoothing prior under `NA`. On
  `y ~ s(x)` (200 rows, seed 21) the draws spread with a per-row sd of
  2.1305 against a residual sd of 0.2829, and the row means of the
  draws sat 4.8183 RMS from the fitted curve (`probe2-base.txt`).
  `pp_check()` on a fit defaults to `NA`, so `pp_check()` on every
  smooth model compared the data with random curves. `predict(NA)`
  keeps a population smooth, so this is in item 2's scope, and it is
  fixed. See item 2.

**Files outside the listed ones that the change forced:** frmtmb.learn's
duplicated-payoff refusal said "simulate() takes no newdata", and its
`test-counterfactual.R` asserted that `simulate()` has no `newdata`
formal. The release suite caught it (69 of 70 with 1 failure). Both now
say what is true; see "Decisions, and notes for the merge". frmtmb.sample's draws
method carried the `error_binned` sentence of item 3 and changed with
it.

## Item 1: `simulate(newdata = )` and `pp_check(newdata = )`

### Construction

`simulate.frmtmb_fit()` (R/predict.R) gains `newdata` and
`allow_new_levels`, after `censored`, so no positional call moves. At
newdata, one replicate:

1. draws `b` by the `re_formula` plan (item 2), exactly as in sample;
2. writes that `b` into a copy of the fit and evaluates every dpar at
   newdata with `frm_linpred(type = "link")` and the link inverse
   (`sim_newdata_dpars()`, R/simulate-newdata.R). This is the
   `re_formula = NA` composition the feasibility report named: the
   redraw is composed with the newdata design by evaluating the design
   at the drawn `b`. An unseen level under `allow_new_levels = TRUE`
   gets `predict()`'s own draw (`predict_new_level_spec()` and
   `predict_new_level_draw()`, called, not copied). An ordinal `cs()`
   term's offsets are rebuilt on the new rows (`ord_cs_values()`).
3. evaluates the addition terms on newdata (`aterms_for_newdata()`),
   so `trials()`, `se()` and `trunc()` follow the rows;
4. rebuilds a residual correlation on the newdata rows, if the fit has
   one;
5. draws with `sim_draw(sim_context(...))`, the one implementation.

`pp_check.frmtmb_fit()` passes `newdata` to `simulate()`, takes `y` from
the newdata's response (`pp_check_newdata_y()`, as the draws method's
`draws_response_values()` does), and reads `group` and `x` from the
newdata.

### The identity it must meet

At the same seed, `simulate(newdata = <the fitted rows>)` must draw what
`simulate()` draws, because only the source of the dpars differs. Run on
13 designs under `NULL` and `NA` (`dev/simnewdata-invariant.R`,
`invariant.txt`):

    gaussian        re_formula = NULL identical
    gaussian        re_formula = NA   identical
    poisson         re_formula = NULL identical
    poisson         re_formula = NA   identical
    binomial        re_formula = NULL identical
    binomial        re_formula = NA   identical
    zip             re_formula = NULL identical
    zip             re_formula = NA   identical
    cumulative      re_formula = NULL identical
    cumulative      re_formula = NA   identical
    cumulative_cs   re_formula = NULL identical
    cumulative_cs   re_formula = NA   identical
    categorical     re_formula = NULL identical
    categorical     re_formula = NA   identical
    trunc           re_formula = NULL identical
    trunc           re_formula = NA   in-sample ERROR: trunc(): rejection sa
    smooth          re_formula = NULL max rel diff 2.09e-13, 3497 of 3500 ce
    smooth          re_formula = NA   max rel diff 2.09e-13, 3497 of 3500 ce
    distributional  re_formula = NULL identical
    distributional  re_formula = NA   identical
    mixture         re_formula = NULL identical
    mixture         re_formula = NA   identical
    student_ar      re_formula = NULL identical
    student_ar      re_formula = NA   identical
    offset          re_formula = NULL identical
    offset          re_formula = NA   identical

The `smooth` rows differ at 2.09e-13 relative: `mgcv::PredictMat()`
rebuilds the basis on newdata and the product differs in the last bits.
The `trunc` row under `NA` is a pre-existing rejection-sampling refusal
that the base build gives identically (`invariant-base.txt`), because
the redrawn effects put one row's mean where the bound excludes nearly
all of its mass. The `binomial`, `cumulative_cs`, `mixture` (the
reporting softmax of `theta` must not reach the simulator) and
`student_ar` (the rebuilt residual correlation) rows are the ones a
wrong construction would break.

### The `re_formula = NA` composition, at a count that can see it

`test-simulate-newdata.R`, "re_formula = NA at newdata shares one group
draw across a level": `y ~ x + (1 | g)`, `sd_g = 1`, `sigma = 0.5`, two
newdata rows in each of two levels, 2000 replicates. The implied
within-level correlation is `sd_g^2 / (sd_g^2 + sigma^2)`, near 0.8
here, and the test asserts the observed mean within-level correlation
within 5 standard errors of it, with the standard error
`(1 - rho^2) / sqrt(2000)` computed in the test; the across-level
correlation within `5 / sqrt(2000)` of 0; and, as the control, the
conditional draws (`NULL`) within `5 / sqrt(2000)` of 0 for the
within-level pairs. A design with rho near 0.8 was chosen over the
probe's 0.116 so the test is many standard errors from both 0 and the
target at an affordable count.

### brms evidence

`dev/simnewdata-brms.R`, `brms.txt`, brms 2.23.0 on its own
`brmsfit_example1`:

- `pp_check(fit1, newdata = fit1$data[1:10, ])`: a ggplot, 10 observed
  values.
- `violin_grouped` with `group = "visit"` and that newdata: a ggplot
  whose groups are the newdata's (`1`).
- newdata without the response: refused, "Response variables must be
  specified in 'newdata'", for `ppc` and for `ppd` alike on this fit.
  frmtmb refuses for `ppc` and answers for `ppd`. brms's `ppd` refusal
  here comes from its `arma()` term, which is fitted with `cov = FALSE`
  and so needs the lagged response; frmtmb's fixture fits `cov = TRUE`,
  whose draw needs no response. This is recorded, not tested against
  brms on a design without `arma()`.
- an unseen grouping level: refused unless `allow_new_levels = TRUE`,
  as frmtmb refuses it.

`dev/simnewdata-brms-arma.R`, `brms-arma.txt`: brms 2.23.0 and frmtmb
on `y ~ 1 + ar(time, gr = g, cov = TRUE)`, 40 groups of 6, AR(1) 0.7,
newdata rows at times 2 and 3 of one new group and time 2 of another:

    rows 1-2 (one group, lag 1): cor 0.592; rows 1-3 (two groups): cor 0.027; draws 750
    time 9, never seen: answers: cor of the two rows 0.620 
    frmtmb: rows 1-2 cor 0.609; rows 1-3 cor 0.002; fitted R[2,3] 0.608
    frmtmb, time 9: ERROR: simulate(newdata = ): the residual correlation term ar(time, gr = g, cov = TRUE) counts lags in the fitted time levels (6 of them), and newdata puts 1 row(s) at a time outside that set (time = 9) 

### Which families refuse newdata, and why

Every family with a structured simulator (`sim_ctx`), because its draw
walks the frame it was fitted on: `mixture(groups = )` draws one class
per FITTED group (`mixture_sim_groups()` reads `ctx$block`), the hidden
Markov family of frmtmb.latent walks the fitted state sequence, and the
learning families of frmtmb.learn walk the fitted trial sequence
(`ln_recurse(ctx$block, ...)`). `mixture_mvn()` goes through the same
slot and is refused with them, although its draw is row by row; it is
also what `predict(newdata = )` refuses. Admitting it is a one-line
change to `sim_newdata_refuse_structure()` and would want its own test.

### Tests seen failing on the base build

`test-simulate-newdata.R` (new, 21 blocks, 40 assertions) and the new
blocks of `test-pp-check-types.R`. On rellib-r3
(`seefail-base-test-simulate-newdata.R.txt`,
`seefail-base-test-pp-check-types.R.txt`):

    RESULT test-simulate-newdata.R pass=5 fail=8 err=14 skip=0 blocks=21
    RESULT test-pp-check-types.R pass=182 fail=1 err=6 skip=0 blocks=21
    RESULT test-draws-methods.R pass=147 fail=0 err=1 skip=0 blocks=22
    RESULT test-counterfactual.R pass=66 fail=1 err=1 skip=0 blocks=15

The five assertions of `test-simulate-newdata.R` that pass on the base
build are controls, not behavior changes: the two row-index identities
of the partial-formula block, its across-term correlation (zero under
both readings), the "a formula naming every term conditions on them"
block, and `expect_false(identical(NULL, NA))`. No assertion sits after
an erroring one in a way that hides a behavior: each block that errors
on base errors at its first `simulate(newdata = )` call, and every
assertion after it in that block needs the same argument.

## Item 2: `simulate()`'s `re_formula`

### The decision it implements, and the one reading it needed

The user decided (2026-09-23) that `~1` means no group effects and a
partial formula is honored. In `simulate()` "no group effect" cannot
mean zero, because `NA` has always REDRAWN the effects there (lme4's
unconditional simulation), and the user's decision named `NA` as the
reading `~1` joins. So a term that is not kept is redrawn, and a kept
term enters at its estimate. That answers the question wt-reunc filed
("redraw the dropped terms, as `NA` does for all of them?"): yes.

brms has no `simulate()`. Its `posterior_predict(re_formula = NA)`
drops the group effects rather than redrawing them. That difference is
0.62.0's and is not moved here; `simulate()`'s `NA` and `predict()`'s
`NA` answered different questions before this lane and still do, and
`?simulate.frmtmb_fit` now says so in a section of its own. What brms
settles is the reading of the formula: `~1`, `~0` and `NA` give
identical draws (`brms.txt`: "~1 identical to NA: TRUE", "~0 identical
to NA: TRUE", "NULL identical to NA: FALSE"), and `~(1 | visit)` on a
`(1 + Trt | visit)` fit differs from both `NULL` and `NA`. brms drops
`~(1 | nosuch)` silently (equal to `NA`); frmtmb refuses it, by the
user's decision of 2026-09-23.

### Construction

`re_resolve()` is split: `re_keep_plan()` returns which columns of
which named terms a formula keeps, with every refusal `re_resolve()`
had, and `re_resolve()` builds its view from that plan. Nothing else in
`re_resolve()` changed, and the bitwise battery below covers the
partial-formula paths of `predict()`, `fitted()` and `frm_linpred()`.

`sim_re_plan()` turns the plan into what one replicate redraws:

- `NULL`, or a formula keeping every column of every term: nothing.
- `NA`, `~0`, `~1`: every block `predict(re_formula = NA)` drops, which
  is every block except a population smooth, `gp()` or `hsgp()`. A
  factor smooth is group-level and is redrawn.
- a partial formula: a block whose columns are all dropped is redrawn
  whole; a block with some columns kept has its dropped columns drawn
  given the kept ones at their estimates,
  `b_D | b_K ~ N(V_DK V_KK^-1 b_K, V_DD - V_DK V_KK^-1 V_KD)`. This is
  refused for `rr`, `gr_cov`, `gr_prec`, `car`, `spde` and Student-t
  blocks, where this code has no conditional law.

`draw_b()` runs once per replicate whenever anything is redrawn whole,
over every block, as 0.62.0's `NA` did, and the positions not redrawn
are overwritten with the estimates. So on a fit without a population
smooth, `NA` consumes the caller's stream exactly as before and returns
the same vector; the bitwise battery shows it.

A `re_formula` that is not `NULL`, `NA` or a formula (`"g"`) used to
condition on everything with nothing said; `check_re_form()` now
refuses it, as it does in `predict()`.

### What the two callers do now

`dev/simnewdata-callers.R` (`callers-base.txt`, `callers-lane.txt`,
`callers-compare.txt`), a `(1 | g)` gaussian and a `s(w)` gaussian,
seeded:

    base: smooth pp_check, sd of replicated values 4.050, of y 1.478
    base: smooth dharma(re_formula = NA), KS p against uniform 1.65e-25
    base: mixed dharma(re_formula = ~1), KS p 0.813
    lane: smooth pp_check, sd of replicated values 1.469, of y 1.478
    lane: smooth dharma(re_formula = NA), KS p against uniform 0.581
    lane: mixed dharma(re_formula = ~1), KS p 0.0101
    dharma mixed default               identical
    dharma mixed NA                    identical
    dharma mixed ~1                    differs
    dharma smooth NA                   differs
    pp_check mixed default             identical
    pp_check mixed NULL                identical
    pp_check mixed ~1                  differs
    pp_check smooth default            differs

- `dharma_residuals()` defaults to `NULL`: identical. With `NA`:
  identical on the mixed model. With `~1`: it now draws what `NA`
  draws, so it moved. Checked from the saved values: on the base build
  `~1` is identical to the default `NULL`, on the lane build identical
  to `NA`. Its KS p-value moved from 0.813 to 0.0101 for that reason:
  a marginal DHARMa residual on a fit with real group effects is not
  uniform, which is `NA`'s long-standing behavior, not a new defect.
- `pp_check()` on a fit defaults to `NA`: identical on the mixed model,
  moved on the smooth model, where the replicated values' sd went from
  4.050 to 1.469 against the data's 1.478, and DHARMa's KS p-value
  under `NA` from 1.65e-25 to 0.581.
- `frm_bootstrap()` keeps 0.62.0's whole-model bootstrap (user
  decision, 2026-09-24; see "Punch round 1").

### Tests seen failing on base

In `test-simulate-newdata.R`: "re_formula = ~1 and ~0 redraw the group
effects, as NA does" (2 failures on base), "a partial re_formula keeps
the named term and redraws the other" (2), "a dropped slope is drawn
given the kept intercept" (2, behavioral since punch round 1: mean
max |z| 33.2 and covariance max |z| 37.4 against the exact conditional
law, `seefail-base-test-simulate-newdata.R.txt`), the three refusals
(3 failures), "re_formula = NA keeps a population smooth curve" (1). In
`test-pp-check-types.R`: "pp_check(re_formula = ~1) redraws the group
effects, as NA does" (1).

### The NEWS exception

The development section says what changed and what breaks. wt-reunc's
paragraph "`simulate()` is the exception and does not change" stays in
the released 0.62.0 section, where it was true; the first pass of this
lane removed it and punch round 1 restored it byte for byte (the
0.62.0 and older sections of `NEWS.md` are identical to `51bfaa4e`).

## Item 3: the multinomial `error_binned` wording

`pp_check_polytomous_noun()` gives "a set of counts over categories"
for multinomial and "a category" otherwise, which is how the
`residuals()` refusal beside it says it. The draws method in
frmtmb.sample carried the same sentence ("two methods should not
describe one case two ways"), so it changed too; that file is outside
this lane's listed files, and the change is the one string. brms
refuses `pp_check()` outright on a multinomial fit ("'pp_check' is not
implemented for this family", `brms.txt`), so there is no brms wording
to match; the refusal and its class were already right.

Assertions, each seen failing on base: `test-pp-check-types.R`,
"error_binned on a multinomial fit calls it counts, not a category"
(core); and one assertion added to the multinomial block of
frmtmb.sample's `test-draws-methods.R`
(`seefail-base-test-draws-methods.R.txt`: 1 error, the old sentence).
There was no multinomial `error_binned` assertion before; the existing
one covers ordinal and categorical and still passes.

## Item 4: `Rplots.pdf`

`dev/simnewdata-rplots.sh` deletes
`extensions/frmtmb.sample/tests/testthat/Rplots.pdf`, runs one file in
its own process with the gated environment on, and checks whether the
file came back.

    SIMNEWDATA_LIB=base
    DRAWS  test-draws-methods.R  RESULT test-draws-methods.R pass=147 fail=0 err=0 skip=0 blocks=22
    DRAWS  test-brms-suite-methods.R  RESULT test-brms-suite-methods.R pass=61 fail=0 err=0 skip=0 blocks=18
    SIMNEWDATA_LIB=lane
    clean  test-draws-methods.R  RESULT test-draws-methods.R pass=148 fail=0 err=0 skip=0 blocks=22
    clean  test-brms-suite-methods.R  RESULT test-brms-suite-methods.R pass=61 fail=0 err=0 skip=0 blocks=18

- `test-draws-methods.R`: the `mcmc_plot()`/`pairs()` block opens
  `pdf(NULL)` and closes it on exit, as the block at line 504 does.
- `test-brms-suite-methods.R` is generated, so the GENERATOR changed:
  `dev/brmsport-gen.R` emits `grDevices::pdf(NULL)` and
  `withr::defer(grDevices::dev.off(), teardown_env())` at the top of
  both packages' methods files. withr is in both packages' Suggests.

## The ledger, before and after

`FRMTMB_PORT_ROOT=<worktree> FRMTMB_PORT_LIB=<lane lib>`,
`dev/brmsport-record.sh`, then `dev/brmsport-ledger.R`. The first
ledger run stopped (`ledger-run1.txt`):

    HOLDS but has manual verdict 'defect': brmsfit-methods:675 in frmtmb
    HOLDS but has manual verdict 'defect': brmsfit-methods:682 in frmtmb
    Error: 2 problem(s); the ledger was not written

Those two rows were removed from `dev/brmsport-verdicts-manual-fit.tsv`
and nothing else was edited by hand.

Before (`before-ledger-summary.md`, the committed 0.62.0 summary) and
after (`dev/brmsport-log/ledger-summary.md`, `ledger-final.txt`):

    | `tests.brmsfit-methods.R` | 223 | 116 | 47 | 30 | 30 |
    | **total** | **494** | **250** | **64** | **35** | **145** |
    | defect | refuses-accepted | 24 |
    Bin 1 passes: 250 of 494 (50.6%).
    | `tests.brmsfit-methods.R (sample)` | 61 | 13 | 6 | 0 | 42 |
    | **total** | **165** | **63** | **7** | **1** | **94** |
    ->
    | `tests.brmsfit-methods.R` | 223 | 118 | 45 | 30 | 30 |
    | **total** | **494** | **252** | **62** | **35** | **145** |
    | defect | refuses-accepted | 22 |
    Bin 1 passes: 252 of 494 (51.0%).
    | **total** | **165** | **63** | **7** | **1** | **94** |

## Nothing else moved

`dev/simnewdata-bitwise.R`: eight designs (gaussian, poisson, a
distributional model, `(1 + x | g)`, `(1 | g) + (1 | h)`, a mixed
poisson, a distributional mixed model, and `s(w) + (1 | g)`) under ML,
`REML = TRUE` and `frmtmb_control(profile = TRUE)`: `vcov()`, `fixef()`,
`logLik()`, `summary()`, `predict()` and `fitted()` at seed 7 with the
caller's `.Random.seed` after each, `predict(re_formula = NA)`,
`fitted(re_formula = ~1)`, `simulate()` under `NULL`, with and without
`seed =`, and under `NA`; on the two-term design the partial
`predict()`, `fitted()`, `frm_linpred(newdata =, se.fit = TRUE)` and the
`~(1 | nosuch)` refusal; `pp_check()`'s plotted data under its default
and under `~1`. Compared with `identical()`:

            identical moved_as_intended           DIFFERS 
                  453                63                 0 

Every quantity outside the items is identical. The moved quantities are
the ones this lane targets: `simulate(re_formula = ~1)`, the partial
`simulate()`, `pp_check(re_formula = ~1)`, `simulate(re_formula = NA)`
on the smooth design, and the caller's stream after each of the first
three, because a draw that now redraws group effects takes more of the
stream. `vcov()` errors on poisson under `REML = TRUE` and `profile`,
identically in both builds; that is the filed defect "vcov() on a REML
poisson fit with no random effects".

## The release tiers

Every log below postdates the last edit to the files it covers: the
last core or frmtmb.sample source edit is 00:33 (R/predict.R,
R/simulate-newdata.R), the lane library was installed at 00:42, the
frmtmb.learn edits are 04:18 to 04:43 and its install 04:43, and the
generated brms-suite files were last written at 04:53.

- **Ungated suite**, all eight packages, one file per process
  (`dev/simnewdata-run-suite.ps1`, a copy of the release driver with
  the lane library first; `suite.log`, 2026-09-24 04:15, summed by
  `dev/simnewdata-suite-compare.R` into `suite-compare.txt`):

      files 280 pass 16458 fail 1 err 0 skip 176 unparsed 0 
                   pkg                  file pass fail err skip
      252 frmtmb.learn test-counterfactual.R   69    1   0    0
      
      files whose count moved against the 0.62.0 baseline:
                 pkg                    file pass pass.base skip skip.base
              frmtmb   test-pp-check-types.R  193       183    0         0
              frmtmb test-simulate-newdata.R   40        NA    0        NA
        frmtmb.learn   test-counterfactual.R   69        70    0         0
       frmtmb.sample    test-draws-methods.R  148       147    0         0

  The one failure is frmtmb.learn's `test-counterfactual.R`, the
  assertion that `simulate()` has no `newdata` (see "What needs the
  user"). Its package was fixed and its 14 files re-run after the edit
  (`dev/simnewdata-run-learn.sh`, `suite-learn.log`,
  2026-09-24 04:56): 14 files, pass 407, fail 0, err 0, skip 13.
  Against the 0.62.0 baseline no file fell except that one: the three
  that rose are this lane's (+10 `test-pp-check-types.R`, +1
  `test-draws-methods.R`, and the new `test-simulate-newdata.R`
  with 40). `test-counterfactual.R` is 69 against the baseline's 70
  because the lane removed the `formals()` assertion and changed one
  other in place.
- **Gated tier**, all environment variables of `run-gated.ps1`
  (`dev/simnewdata-run-gated.ps1`; `gated.log`, 2026-09-24 05:19,
  re-run after the final regeneration): GATED ran 39 of 39;
  39 files, pass 3283, fail 2, err 0. The 2 failures are the stale
  skew-normal pin in `test-drmtmb-agreement.R`, identical on the base
  build (`drmtmb-base.txt`); the release gated log skipped that file
  (13 skips, drmTMB tier off), which is why 0.62.0 did not see it. The
  brms suite files all pass with the new verdicts, including
  `test-brms-suite-methods.R` at 162.
- **R CMD check --as-cran**, built with vignettes
  (`dev/simnewdata-run-check.ps1`): frmtmb Status: 1 NOTE, the
  environmental V8 math-rendering NOTE, in-check tests
  [ FAIL 0 | WARN 16 | SKIP 169 | PASS 10924 ]
  (`check-frmtmb-testthat.Rout`). frmtmb.learn Status: 1 NOTE,
  examples timing on
  `bandit2arm_delta` at 5.9 s elapsed against 2.24 s user, which this
  lane did not touch and which measures load here; in-check tests
  [ FAIL 0 | WARN 0 | SKIP 15 | PASS 399 ].
  frmtmb.sample: the first check ERRORED after its tests reported
  [ FAIL 0 | WARN 15 | SKIP 11 | PASS 1811 ], at "cannot open file
  '.../RLIBS_3f0c32fe5a92/cli/R/sysdata.rdb'", a file in the check's own
  temporary library that disappeared mid-run, with C: at 16 GB free
  (99 percent used); no library went hollow (user 410 of 410, rellib-r3
  9 of 9, lane 3 of 3). Re-run alone (`dev/simnewdata-run-check-sample.ps1`):
  Status: OK, in-check tests
  [ FAIL 0 | WARN 15 | SKIP 11 | PASS 1811 ].

## Punch round 1

The reviewer's verdict was MERGEABLE with no blockers. Its probes are in
`dev/simnewdata-review/`. The user decided two things on 2026-09-24
that changed the round: `frm_bootstrap()` keeps 0.62.0's behavior
(item 7), and `re_formula = NA` must keep every smooth (item 8,
recorded only).

### 1. The residual correlation against brms

Rewritten under "Decisions, and notes for the merge", with the
reviewer's numbers. The
difference is not limited to unseen times. The user decided to keep
frmtmb's reading; `?simulate.frmtmb_fit` and NEWS now say how a lag
is counted and that brms differs.

### 2. The conditional-draw test is behavioral now

"a dropped slope is drawn given the kept intercept" no longer calls
`sim_re_plan()` or `sim_draw_b()`. It simulates the FITTED rows of one
level under `re_formula = ~ (1 | g)` on the reviewer's
`(1 + x1 + x2 | g)` design, 4000 replicates, and compares the three
rows' sample mean and covariance with the exact conditional Gaussian.
The tolerance is 5 standard errors of a sample mean and of a sample
covariance under normality, computed in the test from the law being
tested. Two further assertions check that the instrument can tell the
conditional law from the marginal one (both more than 5 SE away).
Using the fitted rows means the base build answers the same call
instead of failing on a missing argument, and it fails for a real
reason: base conditions on every column, and gives mean max |z| 33.2
and covariance max |z| 37.4 (`seefail-base-test-simulate-newdata.R.txt`).
The reviewer's own end-to-end probe gave 1.79 against the conditional
and 150.9 against the marginal (`dev/simnewdata-review/log/conditional.txt`).

### 4. The newdata design is built once per call

`sim_newdata_design()` builds each linear predictor's newdata pieces
(`lp_eta_design()`) and the `cs()` offsets once; a replicate only
recombines them with its `b`, in `lp_eta_design()`'s own order. A
nonlinear predictor still goes through `frm_linpred()` per replicate.
The aterms and the residual-correlation structure were already built
once.

- Draws unchanged: 37 cases (9 designs; `NULL`, `NA`,
  `allow_new_levels`, `NA` with an unseen level, a partial formula)
  saved from the build before the change and after it, compared with
  `identical()` (`dev/simnewdata-cache-bitwise.R`, `cache-compare.txt`):

      identical 37 of 37 

  The fitted-rows identity (`invariant.txt`) therefore holds as before.
- Cost, `dev/simnewdata-timing.R`: the reviewer's design (500 rows,
  `sigma ~ x`), 2000 replicates per block, three arms interleaved in
  one process, minimum of 5 rounds, with `simulate()` on the fitted
  rows as the control:

      before: per-replicate ms, min of 5: fitted rows 0.635, newdata 6.205, newdata NA 6.625; newdata / fitted 9.77
      after: per-replicate ms, min of 5: fitted rows 0.890, newdata 0.970, newdata NA 1.000; newdata / fitted 1.09

### 5. Unusual newdata

`dev/simnewdata-edge-predict.R`, measured on the base build first
(`edge-predict-base.txt`) and then on the lane (`edge-predict-lane.txt`):

- A missing covariate: `predict()` leaves the row's cells NA (and warns
  in its overflow wording), and `frm_linpred()` gives NA. `simulate()`
  gave NaN with base R's "NAs produced" warning from the family's
  generator (`dev/simnewdata-review/log/edge.txt`). Now the rows whose
  parameters are not finite at the estimates are found once, draw NA,
  and one classed warning names them; the family draws the other rows.
  A residual correlation is built on the finite rows. A row that a
  redrawn effect makes non-finite is NA in that replicate, with one
  warning counting the cells.
- A grouping column absent from newdata: `predict()` and
  `frm_linpred()` DO read a same-named object from the formula
  environment, deliberately (`fill_new_group_vars()` follows
  `model.frame()`). With a global `g` of the right length both answer
  with it; with the wrong length both die on "non-conformable arrays".
  So `simulate()` keeps the environment lookup, as the instruction
  allowed when `predict()` deliberately does otherwise, and refuses the
  wrong length by name instead of dying.
- Tests: "a newdata row with a missing covariate draws NA, as
  predict()" and "a grouping factor read from the environment must fit
  newdata". On the base build they fail only because `newdata` does
  not exist there. The pre-round lane behavior they pin is the
  reviewer's measurement above; that build no longer exists to rerun.

### 6. Rplots, rerun on the final files

`rplots-lane.txt`, written after the regenerated
`test-brms-suite-methods.R` of 04:53:

    SIMNEWDATA_LIB=lane
    clean  test-draws-methods.R  RESULT test-draws-methods.R pass=148 fail=0 err=0 skip=0 blocks=22
    clean  test-brms-suite-methods.R  RESULT test-brms-suite-methods.R pass=61 fail=0 err=0 skip=0 blocks=18

### 7. `frm_bootstrap()`, as the user decided it

`frm_bootstrap()` keeps 0.62.0's whole-model bootstrap: under its
default `re_formula = NA` every replicate redraws the group effects AND
the penalized coefficients of every smooth, `gp()` and `hsgp()` term.
`simulate()`'s public `NA` still holds the curve. The bootstrap reaches
the same body through an internal setting
(`sim_fit_draws(redraw_smooths = TRUE)`, `sim_re_plan(smooths = TRUE)`)
that `simulate()` does not expose.

- Bitwise against 0.62.0 (`dev/simnewdata-boot.R`, 40 replicates,
  seed 1; `boot-compare.txt`): the default on a smooth fit and on a
  mixed fit, and `re_formula = NULL` on each, compared with
  `identical()` against rellib-r3:

      smooth default               identical
      mixed default                identical
      smooth NULL                  identical
      mixed NULL                   identical

  Before this fix the lane's default on the smooth fit DIFFERED
  (`boot-compare-before.txt`), because it had inherited the new `NA`.
- The two settings on `y ~ s(x)` (200 rows, seed 21; `boot-lane.txt`),
  bootstrap sd of the intercept / of the smooth's standard deviation
  `exp(theta)`:

      lane: smooth default, estimate -0.2833 / 5.2011; bootstrap sd 0.0385 / 1.4376; converged 40 of 40
      lane: smooth NULL, estimate -0.2833 / 5.2011; bootstrap sd 0.0432 / 0.3693; converged 40 of 40
      lane: intercept Wald se 0.0460

  The whole-model default spreads the smooth's standard deviation 3.9
  times as far as the conditional bootstrap (1.4376 against 0.3693);
  the intercept's spread is about the same (0.0385 against 0.0432,
  Wald 0.0460).
- **No new argument was added**, although the decision asked for one
  "that conditions on the fitted random effects and smooths". That
  argument already exists: `re_formula = NULL`, the package's name for
  it everywhere, conditions on both, in 0.62.0 and now (bitwise,
  above). A second argument that says the same thing would need a rule
  for when the two disagree. `?frm_bootstrap` now documents both
  settings, and NEWS says so. The user accepted this on 2026-09-24.
- Test: "frm_bootstrap() redraws a smooth by default, and NULL holds
  it". It compares the refitted curve at x = 0.25 and asserts that the
  spread under the default is more than 5 times the spread under
  `NULL` (measured 1.80 against 0.054, ratio 33.7). Seen failing on the
  lane build before this fix, installed from a copy of the tree with
  the old bootstrap call (`boot-test-prefix.txt`, ratio 1.0); passes on
  the lane (`boot-test-lane.txt`) and on rellib-r3, whose behavior it
  pins.
- Not BREAKING: nothing about `frm_bootstrap()` changed against
  0.62.0. The NEWS line says it is unchanged and names `NULL`.

### 8. Recorded: `NA` must keep every smooth

Filed in `dev/test-backlog.md` as decided and waiting for a lane, with
the reviewer's evidence (`rv-smooth.R`, `log/smooth.txt`).

### 3. NEWS

The released 0.62.0 section is restored byte for byte; see "The NEWS
exception" under item 2.

### The tiers, rerun for the round

The machine crashed at about 09:24 during the first rerun; that
`gated.log` was truncated and is void. After the crash the lane library
was checked against the source (`dev/simnewdata-libcheck.R`,
`libcheck-after-crash.txt`): every function the changed files define
matches the installed namespace (frmtmb 163, frmtmb.sample 93,
frmtmb.learn 13, none differ), no library has a hollow package, and a
fit produces a number. Every log below postdates the last source edit
(R/bootstrap.R 08:10 and R/simulate-newdata.R 08:09, installed at 08:14; test files
by 08:30).

- **Gated tier**, all environment variables (`gated.log`, 10:47):
   GATED ran 39 of 39 (39 files: pass 3283, fail 2, err 0, skip 0; summed by `dev/simnewdata-suite-compare.R`).
  The 2 failures are the skew-normal pin in `test-drmtmb-agreement.R`,
  the same on the base build. It is flipped on main in `cad68e21`, after
  this branch's base, so it is expected here until merge.
- **R CMD check --as-cran**, frmtmb and frmtmb.sample, built with
  vignettes (`dev/simnewdata-run-check-r2.ps1`):
  frmtmb: Status: 2 NOTEs; in-check tests [ FAIL 0 | WARN 16 | SKIP 169 | PASS 10928 ] (`check-r2-frmtmb-testthat.Rout`).
  frmtmb.sample: Status: OK; in-check tests [ FAIL 0 | WARN 15 | SKIP 11 | PASS 1811 ] (`check-r2-frmtmb.sample-testthat.Rout`).
  frmtmb's second NOTE is examples timing (`dharma_residuals` 12.25 s
  elapsed, `profile.frmtmb_fit` 5.61 s, `VarCorr` 5.52 s). None of the
  three is this lane's code except through `simulate()`, and the
  example's `simulate()` call is not slower: 8 calls of 250 replicates
  took 0.40 and 0.44 s on the base build against 0.54 and 0.46 s on
  the lane, with a fixed arithmetic control of 0.38 to 0.43 s in the
  same processes (`dev/simnewdata-dharma-timing.R`,
  `dharma-timing.txt`). The NOTE measures load, as it has before on
  this machine.

### The closing Rd change

The sentence on how a lag is counted at newdata (user decision of
2026-09-24) changed `R/predict.R` (12:32) and
`man/simulate.frmtmb_fit.Rd` (12:33, rendered in `rd.txt`). R CMD
check --as-cran on frmtmb afterwards (`dev/simnewdata-run-check-r3.ps1`):

    Status: 2 NOTEs; in-check tests [ FAIL 0 | WARN 16 | SKIP 169 | PASS 10928 ]
    (`check-r3-frmtmb-testthat.Rout`, 2026-09-24 13:22)

The NOTEs are V8 and examples timing, this time on
`profile.frmtmb_fit` alone (9.31 s elapsed), which this lane does
not touch; the previous check named three other examples. It measures
load. frmtmb.sample is unchanged since its check of 11:59.

## What I did not do, and why

- **`predict(newdata = )` still refuses a residual correlation term**,
  and so does `posterior_predict(newdata = )` on draws, while
  `simulate(newdata = )` now answers. `predict()` belongs to lane
  wt-predfix and the draws method to frmtmb.sample. The construction is
  `autocor_for_newdata()` and would serve both.
- **`predict()` at newdata with a `cs()` term.** `predict_dpar_values()`
  calls `with_cs_offsets()`, which builds the offsets from the FITTED
  rows' values and writes them into a per-dpar list under the response
  name. `simulate()` does not use it (it rebuilds the offsets on the new
  rows). Not measured on `predict()`, because `predict()` is another
  lane's; recorded for it.
- **`pp_check(type = "error_binned")` on a fit still bins simulated
  responses**, where brms bins `posterior_epred()` draws (wt-correct's
  filed defect 6). Not in this lane's items.
- **An unseen level of a factor smooth under `NA` at newdata** takes
  the population curve under `allow_new_levels = TRUE`, as `predict()`
  does, rather than a fresh curve. The newdata design has no route to
  draw one. Documented in `?simulate.frmtmb_fit`, and `NA` does not
  imply `allow_new_levels` on such a fit, so it is never silent.
- **`dev/brmsport-gen.R` and `dev/brmsport-ledger.R` write CRLF on
  Windows.** Every generated and ledger file was normalized to LF after
  it was written, so the diff shows content only.
- `mixture_mvn()` at newdata is refused though its draw is row by row;
  see item 1.

## Decisions, and notes for the merge

- **`simulate(re_formula = NA)` and brms's `NA` answer different
  questions** (redraw against drop). That is 0.62.0's state and the user
  decision of 2026-09-23 keeps `NA` as it was; it is recorded so that
  nobody reads the new `~1 == NA` identity as brms parity of the draws
  themselves.
- **DECIDED by the user on 2026-09-24: `frm_bootstrap(re_formula =
  NULL)` is the way to condition on the fitted random effects and
  smooths, and no separate argument is added** (see item 7 of punch
  round 1).
- **A residual correlation at newdata reads time differently from
  brms, at seen times as well as unseen ones** (corrected in punch
  round 1; the first pass said the difference was limited to times the
  fit never saw, which was wrong). frmtmb places each newdata row at
  its fitted time LEVEL and takes the fitted correlation between those
  levels, as the likelihood does at fit time. brms sorts each group's
  newdata rows by time and counts lags by POSITION, ignoring the time
  values. The reviewer measured both on `ar(time, gr = g, cov = TRUE)`
  with brms's posterior mean `ar[1]` 0.608 (phi^2 0.369), correlation
  of two rows of one group
  (`dev/simnewdata-review/rv-brms-time.R`,
  `dev/simnewdata-review/log/brms-time.txt`):

  | newdata times | brms | frmtmb |
  |---|---|---|
  | 2, 3 | 0.619 | 0.622 |
  | 2, 4 (seen, a gap) | 0.656 | 0.389 |
  | 4, 2 (the same rows reversed) | 0.615 | 0.389 |
  | 1, 6 (seen, far apart) | 0.634 | 0.103 |
  | 2, 9 (9 never seen) | 0.613 | refused |
  | 2, 9, 3: times 2 and 9 | 0.390 | refused |

  So frmtmb answers differently from brms at SEEN times with a gap, and
  says nothing there, because under its reading nothing is wrong. brms's
  reading is also not consistent under marginalization: the correlation
  of the rows at times 2 and 9 is 0.613 alone and 0.390 once a time-3
  row is added to the same newdata. DECIDED by the user on
  2026-09-24: frmtmb keeps its time-level reading, as a documented
  departure from brms, because brms's position-based reading is not
  consistent under marginalization. `?simulate.frmtmb_fit` (Structured
  draws) now says how a lag is counted at newdata and that brms
  differs, and the NEWS entry for `simulate(newdata = )` gives the
  2, 4 numbers; no code changed.
- **This lane edited frmtmb.learn, which is outside its listed files.**
  The release suite caught it: `test-counterfactual.R` asserted
  `"newdata" %in% names(formals(simulate.frmtmb_fit))` is FALSE, and
  the duplicated-payoff refusal in `R/family.R` (and
  `?bandit2arm_delta`) told the user "simulate() takes no newdata",
  which item 1 makes false. The message now says `simulate()` refuses
  newdata for a learning family, the test asserts that sentence and the
  structured refusal, and frmtmb.learn has a NEWS entry. Seen failing on
  base (`seefail-base-test-counterfactual.R.txt`). If lane wt-phase3a or
  wt-phase3b touches frmtmb.learn, this is a merge point.
- **The frmtmb.sample draws method's `error_binned` sentence** changed
  with the fit method's; also outside the listed files, one string.
- **`test-drmtmb-agreement.R` fails 2 assertions on the base build and
  the lane build alike** (`drmtmb-base.txt`, `drmtmb-lane.txt`): its
  "frmtmb's default start stalls" block pins the skew-normal stall as a
  KNOWN DEFECT to be flipped when the start is fixed, and 0.62.0 fixed
  it. The pin is flipped on main in commit `cad68e21`, after this
  branch's base, so on this branch the 2 gated failures are expected
  until merge. This lane did not edit that file.
