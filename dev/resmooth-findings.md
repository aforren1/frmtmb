# Lane wt-resmooth: `re_formula = NA` keeps every smooth

Worktree `C:\Users\adf44\source\r\frmtmb-wt-resmooth`, branch
`wt-resmooth`, private library
`C:/Users/adf44/source/r/wt-resmooth-lib`. The "before" arm is the
0.64.0 reference build in `C:/Users/adf44/source/r/rellib-r3`; the
scripts that have a before/after arm take `RESMOOTH_LIB=base` to select
it, and section 10 says which those are.

## 1. The defect, reproduced

`dev/resmooth-probe.R` (seed 5, n = 300, the data of
`dev/simnewdata-review/rv-smooth.R`), logs `dev/resmooth-before.txt` and
`dev/resmooth-after.txt`. `max|NA-NULL|` is the largest row difference
between `fitted(re_formula = NA)` and `fitted(re_formula = NULL)`;
`rowsd/sigma` is the median over rows of the per-row standard deviation
of 1500 `simulate(re_formula = NA)` draws over the fitted sigma.

                          blocks NA keeps  max|NA-NULL|    rowsd/sigma
    construction           before  after  before   after  before after
    ---------------------- ------ ------  ------ -------  ------ -----
    s(x)                        1      1       0       0   1.000 1.000
    s(x) + s(g, bs = "re")      1    1,2   0.906  4.4e-16  1.260 1.000
    s(x, g, bs = "fs", k=5)  none  1,2,3    2.08  4.4e-16  1.741 1.000
    t2(x, g, c("cr","re"))   none    1,2    1.84  4.4e-16  1.919 1.000
    s(x) + (1 | g)              2      2   0.906    0.906  1.258 1.258
    sigma ~ s(x)              1,2    1,2       0        0  1.000 1.000
    t2(x, z)                1,2,3  1,2,3 4.4e-16  4.4e-16  1.000 1.000
    gp(x)                       1      1       0        0  1.000 1.000

The in-sample and the `newdata` route agree in both arms (the probe
reports the kept blocks separately for each). `4.44e-16` is not zero and
is not claimed to be: the `NA` path sums one product per smooth block
while the `NULL` path takes one product over the whole `Z`, so the two
differ in summation order. `t2(x, z)` already showed that figure before
the change.

`conditional_effects(re_formula = NA)` on the same fits, from the same
probe: on `y ~ s(x, g, bs = "fs", k = 5)` the drawn curve ran from
0.4901 to 0.4901, a FLAT line, because the term was the whole model and
was removed; it now runs from -0.1814 to 1.7905.

`predict.frmtmb_fit()` SIMULATES the predictive distribution, so its
`Estimate` cannot be used to test a design: on `y ~ x` it sits 0.089
from `fitted()` (`dev/resmooth-ndcheck.R`, reference build). The probe
uses `frm_linpred()` for that comparison instead. This cost one probe
run and is recorded so the next reader does not repeat it.

## 2. What brms does (item 3 and item 4 of the task)

`dev/resmooth-brms.R`, log `dev/resmooth-brms.txt`, four brms 2.23.0
fits on the same data (`chains = 2, iter = 1000, seed = 7`), cached in
`dev/resmooth-brms-fits.rds`.

**brms keeps every smooth under `re_formula = NA`, bitwise.**

    fit                               identical(NA, NULL)  max|NA-NULL|
    --------------------------------  -------------------  ------------
    y ~ s(x) + s(g, bs = "re")        TRUE                 0
    y ~ s(x, g, bs = "fs", k = 5)     TRUE                 0
    y ~ t2(x, g, bs = c("cr","re"))   TRUE                 0
    y ~ s(x) + (1 | g)                FALSE                1.966

`identical()` is on the whole draws matrix and `max|NA-NULL|` with it;
the last row's 1.966 is over draws and rows, and the same difference on
the posterior mean per row is 0.898 (section 2's parity table).

**A new level of the grouping factor.** brms refuses it for all three
smooth constructions, with "New factor levels are not allowed", at
`re_formula = NA` and at `NULL`, and with `allow_new_levels = TRUE` as
well as without it. The grouping factor of a smooth is an ordinary
predictor to brms, not a group-level grouping factor, so it never
reaches the `allow_new_levels` machinery. A missing column is
"The following variables can neither be found in 'data' nor in
'data2': 'g'". On `y ~ s(x) + (1 | g)` brms behaves as frmtmb does: `NA`
answers at a new level, `NULL` needs `allow_new_levels`.

mgcv itself, at an unseen level (same log):

- `s(g, bs = "re")`: `PredictMat()` returns a 5 x 11 matrix where the
  fit has 10 columns, silently, so a caller that does not check the
  level mis-indexes the coefficient vector. frmtmb refuses first.
- `s(x, g, bs = "fs", k = 5)`: a 5 x 50 matrix of ZERO rows, which is
  the population curve.
- `t2(x, g, bs = c("cr","re"))`: `PredictMat()` errors,
  "non-conformable arguments".

**Decision, recorded:** frmtmb MATCHES brms by refusing, by name, both a
missing grouping column and an unseen level at every `re_formula`, and
keeps one opt-in brms does not have: `allow_new_levels = TRUE` gives an
unseen level of an `fs` term the population curve, which is mgcv's own
zero row. An `re` basis or margin is refused under
`allow_new_levels = TRUE` too, because it has no zero row. That opt-in
predates this lane and is more permissive than brms; it is kept rather
than removed because it answers a real question (the population curve)
that `re_formula = NA` no longer answers.

**The parity check.** brms's coefficients cannot be pushed into a frmtmb
fit: the two rotate the basis differently (brms 19, 50, 50 and 9 design
columns against frmtmb's 2 + 18, 1 + 50, 1 + 50 and 2 + 8). The task's
stated fallback is used instead, and sharpened into a numerical
statement: least squares of brms's posterior-mean `posterior_epred(NA)`
on frmtmb's `re_formula = NA` design, residual relative to the centred
target (`dev/resmooth-parity.R`, log `dev/resmooth-parity.txt`).

Columns: the residual of brms's epred(NA) on frmtmb's NA design now, on
the design 0.64.0's rule produced, and on brms's own design (the floor);
then the reverse direction, frmtmb's eta(NA) on brms's design.

    fit                    now       0.64.0    brms's own  reverse
    ---------------------  --------  --------  ----------  --------
    s(g, bs = "re")        7.81e-14  5.923e-1  2.53e-15    8.27e-14
    s(x, g, bs = "fs")     4.52e-15  1.000e+0  1.78e-15    4.44e-15
    t2(x, g, cr+re)        1.52e-15  1.000e+0  1.42e-15    1.63e-15
    s(x) + (1 | g)         9.36e-14  9.360e-14 3.10e-15    9.98e-14

The two designs now span each other's answer at the level of the linear
algebra. Under 0.64.0's rule the `fs` and `t2` designs explained NOTHING
of the centred variation (residual 1.000: what was left was the
intercept), and the `s(g, bs = "re")` design left 59 percent of it. The
`s(x) + (1 | g)` row is unchanged by the fix, which is the control: `NA`
drops the intercept effect and keeps the smooth, before and after.

## 3. What changed in the package

- `R/predict.R`, `lp_eta_design()`: the `setdiff(sm_ids,
  smooth_group_block_ids(lp))` under `!use_re` is gone. Every smooth
  block is kept at population level.
- `R/predict.R`, `pred_design()`: `drop_grp` is gone, so the `newdata`
  route rebuilds every smooth's basis at every `re_formula`.
- `R/predict.R`, `smooth_newdata_check()`: the `use_re` argument is
  gone, because the checks now apply at every `re_formula`. The
  missing-column and unseen-level messages no longer offer
  `re_formula = NA` as the way out, and say why. The `re`-basis refusal
  names "an `re` factor in its basis" rather than `bs = "re"`, so a
  `t2()` with an `re` margin is not called a `bs = "re"` term.
- `R/simulate-newdata.R`, `sim_group_block_ids()`: every smooth, `gp()`
  and `hsgp()` block is excluded, so `simulate(re_formula = NA)` redraws
  none of them.
- `R/simulate-newdata.R`, `sim_re_plan()`: `every` (may an unseen level
  at `newdata` be invented?) is now FALSE whenever the FIT has a
  group-indexed smooth, rather than whenever a REDRAWN block is a
  smooth, which after the change would never be true. `redraw` is
  `integer(0)` rather than `unlist()`'s `NULL` when nothing is redrawn.
- `R/interop.R`, `emm_grid_vars()`: a group-indexed smooth's grouping
  factor joins the emmeans reference grid, whatever `use_re` says. It was
  in neither source the grid reads: `ce_plot_vars()` leaves it out on
  purpose (it is not a curve to draw) and `emm_group_vars()` reads bar
  terms, `car()` and `spde()` only. Without this the change turned
  `emmeans()` on such a fit from an ANSWER into a refusal for want of the
  column. Measured on `y ~ f + s(x, g, bs = "fs", k = 5)`, `g` crossed
  with `f` (`dev/resmooth-emm.R`, logs `dev/resmooth-emm-before.txt` and
  `-after.txt`): `emmeans(fit, ~f)` gave 0.1389 and 0.8801 on 0.64.0 with
  the term dropped, ERRORED on this lane's first pass, and gives 0.1424
  and 0.8836 now, averaged over `g`'s levels. The CONTRAST across `f` is
  0.7412 in both arms, as it has to be: the smooth names no `f`, so the
  difference is the fitted coefficient whichever way the smooth is
  treated and only the level of the means moves. On the `t2()` fit the
  means go from 0.0766 and 0.8108 to 0.1494 and 0.8836, again with the
  contrast 0.7342 unchanged. The same change makes
  `emmeans(re_formula = NULL)` work on such a fit, where 0.64.0 refused
  it for exactly the same missing column.

  `g` has to be CROSSED with `f` to measure this, and the first
  construction was not. `g` cycling 1..10 against an alternating `f` puts
  both on period 2, so `g`'s parity fixed `f`, emmeans read the design as
  nested ("g %in% f") and contrasted something else: the contrast came
  out 0.9 percent from the coefficient instead of at round-off. The test
  and the probe both use `rep(1:10, each = 30)` now, and say why.
- `R/predict-brms.R`: new `smooth_b_idx()`; see section 5.
- `R/predict.R`, `fitted_point_se()` and `R/autocor.R`,
  `autocor_cond_fd_se()`: the finite-difference standard error
  differences every smooth's coefficients at every `re_formula`. See
  section 5.
- `R/re-formula.R`: the partial-formula refusal's REASON is corrected.
  Behavior unchanged; see section 4.
- Docs: `?frm_linpred` ("What `re_formula = NA` drops"), the
  `allow_new_levels` parameter of the same block, `?simulate.frmtmb_fit`
  (group-level terms, new data), `?pp_check`, `?frm_bootstrap`,
  `?frm_lp_basis`, and `vignette("case-studies")`.
- `NEWS.md`: a new `# frmtmb (development version)` heading above
  `# frmtmb 0.64.0`, with a "Breaking changes" section that says which
  default `pp_check()` results change.

`re_governed_b()` was NOT changed, though its name says "the blocks
`re_formula` governs" and `re_formula` no longer governs a group-indexed
smooth. Removing the group smooth from it would have regressed a
measured fix of lane wt-reunc (`Est.Error` 68 percent wrong on a
cumulative fit, `dev/reunc-log/fdsmooth-prefix.txt`). Its doc block now
says why the smooth is there.

`vignette("brms-migration")` has no row on the `re_formula` RULE, only
two rows on the argument's spelling, so it needed no change. The row
that did need one was in `vignette("case-studies")`, and it was not
prose alone: the `fosr-fig` chunk called
`frm_linpred(ffs, newdata = xg, re_formula = NA)` with NO `subject`
column, which after this change is a named refusal. The chunk now
predicts at an unseen subject with `allow_new_levels = TRUE`, which is
the population coefficient function, and the closing bullet says so.

## 4. The partial-formula refusal: a proposal, not a change

Its stated reason is now false, so this section says so and proposes;
it does not decide.

`re_keep_plan()` refuses a PARTIAL `re_formula` on a fit that also has
group-level content a formula cannot name. Its stated reason was
"`re_formula = NA` drops it with every other group-level term and NULL
keeps it, but a partial formula cannot say which". For a factor smooth
that reason is now FALSE: `NA` keeps it and `NULL` keeps it, so there is
nothing for a partial formula to resolve. For a `car()` or `spde()`
field the reason still holds exactly.

The 0.63.0 handoff settled that a partial `re_formula` stays refused on
a fit with a factor-smooth term, and this lane was told to say so and
propose rather than decide. So:

- The refusal is UNCHANGED. `frm_linpred(fit, re_formula = ~(1 | h))` on
  `y ~ s(x, g, bs = "fs", k = 5) + (1 | h) + (1 | f)` is still an error
  of class `frmtmb_error` matching "cannot name"
  (`dev/resmooth-edge-before.txt`, `dev/resmooth-edge-after.txt`, same
  in both arms).
- Its MESSAGE is corrected, because it stated something the code no
  longer does. It now says a partial formula cannot say whether that
  content stays, and points at `NULL` and at `NA`, naming what each
  does. No test pins the old sentence; two pin the fragment "cannot
  name" and both still pass.
- **Proposal, for the user to decide:** drop the factor-smooth entry
  from `re_fit_components()$unnamed`, so that a partial `re_formula`
  beside a factor smooth is ACCEPTED and keeps the smooth, as `NULL`
  and `NA` both do, while `car()` and `spde()` keep the refusal. That
  would make `re_formula` mean one thing for smooths at every value.
  The cost is that a user who wrote a partial formula expecting the
  smooth to go would now silently get it, which is the situation
  `s(x)` is already in.

**The whole accept-and-refuse surface, the facts the decision needs.**
Measured by the review on
`gaussian(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | f) + (1 | h)))`,
identical on both arms (`dev/reviews/2026-09-29-resmooth.md` section 5,
script `dev/resmooth-rev-probe2.R`):

    re_formula          result             smooth kept?
    ------------------  -----------------  -----------------------
    NULL                OK                 yes
    NA                  OK                 yes, max|-NULL| 0.529
    ~0                  OK, equals NA      yes
    ~1                  OK, equals NA      yes
    ~(1|f) + (1|h)      OK, equals NULL    yes
    ~(1|f)              REFUSED            unanswerable by fiat
    ~(1|nosuch)         REFUSED by name    n/a

and on a fit with ONE bar term, `~(1 | f)` is accepted and equals `NULL`
exactly. So the smooth is kept at EVERY `re_formula` the package accepts,
and the refusal fires on exactly one shape: a formula that drops some
columns of a bar term on a fit that also carries a factor smooth. Nothing
about the smooth is being guessed at that shape, which is why the
message no longer says it is.

The message took two passes to get right, and both faults were the same
kind: it described a condition the code does not test.

As shipped in 0.64.0 it said a partial formula "cannot say whether that
content stays" and then, in the same sentence, that `NULL` and `NA` both
keep it, which is a contradiction; and it called a smooth "group-level
content", which this lane's own `?frm_linpred` says it is not.

The first rewrite said the formula "drops some columns of a group-level
term" and called that "the one shape where a formula reaches inside a
term". That is false of the case that actually fires most: `~(1 | f)` on
`y ~ s(x, g, bs = "fs") + (1 | f) + (1 | h)` drops the WHOLE `(1 | h)`
term and reaches inside nothing. Its advice, "name every column of the
terms you want", is also what that user had already done. The condition
`re_keep_plan()` tests is `!all(vapply(keep, all, NA))`: at least one
column of at least one addressable component is not kept, which a
wholly unnamed term satisfies as much as a partly named one.

The message now states that condition as the code has it, "names fewer
group-level terms, or fewer columns of one, than this fit has", and its
advice is "name every group-level term and every column of each, or use
re_formula = NULL to keep them all, or NA to drop them". The refusal
itself is unchanged through all of this, and the two tests that pin the
fragment "cannot name" still pass.

`re_formula = ~(1 | nosuch)` stays refused, unchanged and untouched, as
the 0.63.0 handoff settled.

## 5. The finite-difference standard error left out every smooth

Found because the fix put a term into the ESTIMATE at
`re_formula = NA` whose uncertainty the finite-difference route then
left out. It turned out to be older and larger than that.

`fitted()` on an ordinal, categorical or multinomial fit reports
`Est.Error` per category through `fit_fd_se()`, and
`autocor_cond_fd_se()` does the same for a `cov = FALSE` autocor
quantity. Both took their `b` positions from `re_governed_b()`, which
excludes a POPULATION smooth always and returns nothing at all at
`re_formula = NA`. So the differenced function was "move the thresholds
with the fitted curve FROZEN", which is not the standard error of
anything. The analytic route (`lp_delta_A()`) has always carried every
smooth's columns, so the two routes disagreed.

Measured against Monte Carlo over the same joint covariance the delta
method uses: draw the outer parameters AND the smooth's coefficients
jointly, recompute the category probabilities, take the standard
deviation. `dev/resmooth-mcse.R` (4000 draws, seed 4), log
`dev/resmooth-mcse.txt`. The Monte Carlo standard error of an sd at 4000
draws is 1.1 percent, so a few percent is agreement and 50 to 75 percent
is not.

Median over the nine reported cells of the relative gap from the Monte
Carlo standard deviation, with the smooth's coefficients differenced and
without them:

    cumulative fit                        with     without
    ------------------------------------  -------  -------
    s(x, k = 8), 6 coefficients           0.0209   0.7447
    s(x, g, bs = "fs", k = 5), 30 coefs   0.0495   0.4990

**The residual, stated as a range and not as one number.** "Within 2
percent" is true of the POPULATION-smooth fit alone. The table above
gives 2.09 percent there and 4.95 percent on the factor-smooth fit, and
the review's own nine fits put the worst cell at 7.85 percent, on
`s(x) + (1 | g)` at `re_formula = NULL`
(`dev/reviews/2026-09-29-resmooth.md` section 1). So the honest statement
is **within 2 percent on a population smooth and within 8 percent beside
a group effect**. That residual is delta-method curvature, not a missing
term: it is present on BOTH arms and is the pre-existing accuracy of a
first-order method applied to a curved function of a random intercept.
The claim is corrected in these words in `NEWS.md` too.

The error is not one-signed, which is why it was worth measuring rather
than reasoning about: the variance carries the cross term
`2 Jd' Vdb Jb`, so dropping the `b` columns moved the answer UP on the
population-smooth fit (0.37825 against a Monte Carlo 0.06509 on the
first row, a factor of 5.8) and DOWN on the factor-smooth fit (0.02607
against 0.06270). `dev/resmooth-popse-before.txt` and
`dev/resmooth-fdse.txt` hold the two.

Fixed by `smooth_b_idx()`, used in `fitted_point_se()` and
`autocor_cond_fd_se()`: every smooth, `gp()` and `hsgp()` block's `b`
positions are differenced at every `re_formula`, and `re_governed_b()`'s
group-level positions join them when `re_formula` keeps the group
effects. After the fix the shipped value equals the reference that
differences those coefficients to 27.9 ulp on the population-smooth fit
and to 1.6 ulp on the factor-smooth one, which is the same arithmetic in
a different summation order.

Two tests pin it, both SEEN TO FAIL on the reference build (see section
7). Both read the smooth's `b` positions off the frame rather than
calling `smooth_b_idx()`, so the failure on a build without the fix is
behavioural and not "that symbol does not exist".

## 6. What was deliberately NOT done

- **`frm_bootstrap()` is untouched.** Its default is
  `re_formula = NA` with `redraw_smooths = TRUE`, which takes the
  `seq_along(blocks)` branch of `sim_re_plan()` and is therefore
  unaffected by the change to `sim_group_block_ids()`. Measured:
  `sim_re_plan(fit, NA, smooths = TRUE)` returns blocks 1,2,3 on
  `y ~ s(x, g, bs = "fs", k = 5)` before and after, while
  `sim_re_plan(fit, NA)` returns none after and 1,2,3 before
  (`dev/resmooth-edge-before.txt`, `dev/resmooth-edge-after.txt`). The
  0.63.0 decision holds: the whole-model parametric bootstrap redraws
  the smooths, `simulate()` does not.
- **The extensions needed no change.** `extensions/frmtmb.sample/R` has
  no use of `smooth_group_block_ids`, `sim_group_block_ids`,
  `sim_re_plan`, `lp_eta_design` or `re_form_keeps`: its draws methods
  pass `re_formula` to core's `frm_linpred()`, so they inherit the fix.
  `ce_group_vars()` in core already excludes smooth, `gp()` and `hsgp()`
  blocks, and `frmtmb.sample`'s `conditional_effects` draws method calls
  it, so neither blanks a factor smooth's grouping column. Their suites
  were run anyway; section 7 has the counts.
- **`re_governed_b()` was not narrowed**; see section 3.
- **The partial-formula refusal was not removed**; see section 4.
- **`allow_new_levels = TRUE` on an `fs` term was not made to refuse**,
  which would have matched brms exactly. It is the only remaining way
  to ask for the population curve of such a model, and removing it
  would leave the question unanswerable. Recorded as a deliberate
  difference from brms in `?frm_linpred` and in NEWS.

## 7. Tests

`dev/resmooth-run-one.R <file> [base]` runs one file in one R process
and prints `RESULT <file> lib=<lane|base> pass= fail= error= skip=
warn=`, summed over the file's blocks, with each failing block's
message. `base` selects the 0.64.0 reference build, which is how each
new assertion was SEEN TO FAIL.

### Seen to fail on the reference build

    file                            base arm            lane arm
    ------------------------------  ------------------  ---------
    test-smooth-population.R        26 pass 11 fail 1 error  40 pass
    test-simulate-newdata.R         49 pass  6 fail 0 error  55 pass
    test-pp-check-types.R          196 pass  2 fail 0 error 198 pass
    test-predict-re-uncertainty.R   68 pass  2 fail 0 error  78 pass
    test-emmeans.R                  41 pass  1 fail 1 error  53 pass
    (frmtmb.spline) test-curve.R    ran green, see below     65 pass

Every lane-arm figure above is 0 fail, 0 error, 0 skip.

**Two figures in this table were wrong when first written, and the
review caught both.** `test-smooth-population.R` is 26 pass 11 fail, not
25/12: the 25/12 came from a run made BEFORE the three cross-package
tolerances in that file were loosened from 1e-9 to 1e-5, so one assertion
that then failed now passes. `dev/resmooth-seenfail-test-smooth-
population.R.txt` is that superseded run and is kept as the record of it;
`dev/resmooth-seenfail2-test-smooth-population.R.txt` is the base arm of
the file as shipped. `test-predict-re-uncertainty.R` at 68 pass 2 fail is
right for the file as shipped, and the 67 pass 1 fail in
`dev/resmooth-seenfail-fdse.txt` is a snapshot from when only one of its
two standard-error blocks existed;
`dev/resmooth-seenfail2-test-predict-re-uncertainty.R.txt` re-measures
it. Both errors are the same kind: a count copied from a log that a later
edit superseded. The generated-counts discipline was applied to the suite
totals and not to this hand-built table, which is where it should have
been applied too.

`extensions/frmtmb.spline/tests/testthat/test-curve.R` is the one file
whose base arm is not a useful measurement: its old block asserted that
`frm_curve(re_formula = NA)` on an `fs`-only fit gives a CONSTANT, which
was true of 0.64.0 and is the thing the fix removes. The block was
rewritten, so the base arm would be measuring a test that no longer
exists. What was seen instead is the failure the change caused: that file
reported 57 pass with 1 error against the lane build before the block was
rewritten, on the named refusal for the missing grouping column
(`dev/resmooth-ext/frmtmb.spline-test-curve.R.txt` at the time). It now
asserts what the new rule says, including the refusal and that two
subjects give different curves.

The failing blocks on the base arm, by name:

- `re_formula = NA keeps the fs smooth and drops (1 | g)`
- `an fs term needs its grouping column at every re_formula`
- `an unseen fs level errors at every re_formula`
- `s(g, bs = 're') is a smooth, so re_formula = NA keeps it`
- `a t2() with an re margin is kept by re_formula = NA`
- `conditional_effects draws the reference level's fs curve`
- `re_formula = NA keeps a GROUP-INDEXED smooth curve too`
- `pp_check reproduces an fs fit's spread at the default re_formula`
- `Est.Error at re_formula = NA carries a kept smooth's uncertainty`
- `Est.Error of a category probability carries a POPULATION smooth`
- `a group-indexed smooth puts its factor in the reference grid`

The emmeans block is the one whose base-arm failure is on its CAUSE and
not its symptom: on the base build `emmeans()` answers on such a fit, so
the assertion that fails there is the one on the reference grid's
variables, which is what the answer was wrong about.

### The pp_check assertion, and why it is a ratio

`dev/resmooth-ppcheck.R` (seed 5, n = 300, `nsim = 400`, seed 3), logs
`dev/resmooth-ppcheck-before.txt` and `-after.txt`. `pp_check()` calls
`simulate(re_formula = NA)`, so the statistic is taken from that.

                           rowsd/sigma  sd(yrep)/sd(y)  identical(NA,
    construction           before after  before  after  NULL) bef/aft
    ---------------------- ------ -----  ------ ------  --------------
    s(x) + s(g, bs = "re") 1.1363 1.0023 1.0013 0.9930  FALSE  TRUE
    s(x, g, bs = "fs", k5) 2.9348 1.0023 1.3359 0.9845  FALSE  TRUE
    t2(x, g, c("cr","re")) 3.2380 1.0023 1.3530 0.9848  FALSE  TRUE
    s(x) + (1 | g)         1.1395 1.1395 1.0017 1.0017  FALSE  FALSE
    s(x)                   1.0023 1.0023 0.9951 0.9951  TRUE   TRUE

The test asserts both ratios against `6 / sqrt(2 * (nsim - 1))`, which
is six times the sampling standard deviation of ONE row's sd estimate at
`nsim` draws. At `nsim = 400` that is 0.2124: the before values 1.9348
and 0.3359 from 1 are outside it and the after values 0.0023 and 0.0155
are well inside. The bound is derived from the run's own `nsim` and is
deliberately generous for the pooled statistic, whose sampling error is
smaller. The complement is asserted too: on `y ~ s(x) + (1 | g)` the
ratio must be `sqrt(1 + tau^2/sigma^2)` from the fit's own estimates,
1.1395 measured.

`P(sd(yrep) >= sd(y))` was measured as well (0.9750 before against
0.2450 after on the `fs` fit) and is NOT asserted: on
`s(x) + s(g, bs = "re")` it was 0.5000 before, so it does not
discriminate, and a tail probability near 0.975 is not outside a
conventional band. The dissolved statistic is recorded rather than
deleted, because it is a result about the statistic.

### The guard with its positive condition absent

`sim_re_plan()`'s `every` flag decides whether `simulate(newdata = )`
may invent an unseen level. Written the obvious way after the fix (over
the REDRAWN blocks, which now hold no smooth) it would read TRUE, and
then a `newdata` row at an unseen `fs` level would be given the
population curve while the draw claimed a fresh level. Both arms are
built in `dev/resmooth-edge.R` and recorded in
`dev/resmooth-edge-after.txt`:

    design with the guard absent (nl_ok = TRUE)  OK
    design with the guard present (nl_ok = FALSE) ERROR: New levels in
      the factor-smooth term s(x,g): 99 ...

### The whole core suite and the two extension suites

Emitted by `dev/resmooth-summarise.R` into
`dev/resmooth-suite-counts.txt` and pasted verbatim, because a typed
count is a sentence no verifier checks. The summariser REFUSES to print
totals until each driver log says DONE, and that refusal was observed:
a redirect creates a result file before the R process writes anything,
so the file count reached 179 while four gated brms files were still
running, and an earlier run of the summariser reported 12304 passes over
175 usable files and named the four.

    == GENERATED BLOCK, paste verbatim ==
    core frmtmb suite: 179 files
      pass=14040 fail=0 error=0 skip=14 warn=17
      no file failed, errored or failed to report
      files with a skip: 2
        test-drmtmb-agreement.R.txt skip=13
        test-fuzz.R.txt skip=1
    frmtmb.sample + frmtmb.spline suites: 50 files
      pass=2623 fail=0 error=0 skip=2 warn=13
      no file failed, errored or failed to report
      files with a skip: 2
        frmtmb.sample-test-scale.R.txt skip=1
        frmtmb.spline-test-scale.R.txt skip=1
    == end ==

Both files with a skip were then run with their own gate set, and the
two `test-scale.R` skips are performance tiers that gate themselves.

The core suite above ran against the library as it stood before this
lane's last three edits (the `emm_grid_vars()` fix and two comment
edits). Those are covered instead by `R CMD check --as-cran` on a
tarball built from the FINAL source, whose `checking tests` step passed
in 319 s, and by a rerun of the reachable files one per process against
a second private library and then against the canonical one; the
per-file counts are in `dev/resmooth-t2-*.txt` and
`dev/resmooth-final-*.txt`.

### R CMD check --as-cran, once, on the final source

    Status: 1 NOTE
    * checking examples ... [29s] OK
    * checking examples with --run-donttest ... [32s] OK
    * checking tests ... [319s] OK
    * checking re-building of vignette outputs ... [163s] OK
    * checking PDF version of manual ... OK
    * checking HTML version of manual ... [20s] NOTE
      Skipping checking math rendering: package 'V8' unavailable

That NOTE is the one `dev/lane-rules.md` records as expected on this
box. The vignette step matters here: `vignette("case-studies")` has a
chunk this change would otherwise have broken, and it re-built.

## 8. Defects found and NOT fixed

- **A `t2()` with an `re` margin cannot be predicted at an unseen level
  at all.** `mgcv::PredictMat()` errors with "non-conformable
  arguments" (`dev/resmooth-brms.txt`), and frmtmb refuses earlier with
  its own named message. Not a frmtmb defect; recorded so the refusal is
  not mistaken for a frmtmb limitation that could be lifted.
- **`s(g, bs = "re")` at an unseen level makes `mgcv::PredictMat()`
  return a matrix of the WRONG WIDTH** (5 x 11 where the fit has 10
  columns) rather than erroring. frmtmb's check fires first, so nothing
  reaches it here, but any other caller of `PredictMat()` on a stored
  `random.effect` basis is exposed. Recorded, upstream.
- **The machine cannot compile a fresh Stan program with the shared
  library as it stands.** `StanHeaders` 2.39.1 against `rstan` 2.32.7
  dies in `compileCode()` at `make: *** Error 1`, which is exactly the
  failure `dev/lane-rules.md` documents, and
  `C:/Users/adf44/source/r/pinlib` does not exist on this machine. The
  brief for this lane said Stan compiles fresh and works; it does not.
  Worked around by installing `StanHeaders` 2.32.10 from the CRAN
  archive into THIS LANE's private library only
  (`C:/Users/adf44/source/r/wt-resmooth-lib`, tarball and log beside
  it), after which all four brms fits compiled and sampled in about 40
  seconds each. The shared `pinlib` should be rebuilt; this lane did not
  touch anything outside its own library.

## 9. Version bump

This changes the ANSWER of `predict()`, `fitted()`, `frm_linpred()`,
`simulate()`, `pp_check()`, `conditional_effects()` and
`posterior_epred()` at `re_formula = NA` on any fit with a
group-indexed smooth, changes the `Est.Error` of a category probability
on any fit with any smooth, and turns two silent answers into named
refusals. It is a breaking change of behavior, not a patch: a MINOR
bump. The number is not chosen here.

## 10. Where the numbers come from

Every script is in `dev/`. The five that have a before/after arm take
`RESMOOTH_LIB=base` to run against the 0.64.0 reference build instead of
this lane's library: `resmooth-probe.R`, `resmooth-edge.R`,
`resmooth-ppcheck.R`, `resmooth-fdse.R` and `resmooth-popse.R`.
`resmooth-run-one.R` takes `base` as its second argument for the same
purpose. The rest are lane-only: `resmooth-ndcheck.R` measures the
reference build directly, `resmooth-brms.R` uses no frmtmb at all, and
`resmooth-parity.R` and `resmooth-mcse.R` need the fix to be present.

    dev/resmooth-probe.R        section 1     before/after logs
    dev/resmooth-ndcheck.R      section 1     predict() simulates
    dev/resmooth-brms.R         section 2     brms's rule, new levels
    dev/resmooth-parity.R       section 2     design-space parity
    dev/resmooth-edge.R         sections 4,6  the refusals, the guard
    dev/resmooth-ppcheck.R      section 7     the pp_check statistic
    dev/resmooth-fdse.R         section 5     the fs Est.Error gap
    dev/resmooth-popse.R        section 5     the s(x) Est.Error gap
    dev/resmooth-mcse.R         section 5     the Monte Carlo reference
    dev/resmooth-run-one.R      section 7     one test file, one process
    dev/resmooth-summarise.R    section 7     the generated counts block
    dev/resmooth-emm.R          section 3     emmeans on such a fit
    dev/resmooth-msg.R          section 4     the refusal on both shapes
    dev/resmooth-vignette-chunk.R  section 3  the edited vignette chunk

`dev/resmooth-brms-fits.rds` is 10 MB of cached brms fits and
`dev/resmooth-brms-design.rds` 99 KB of design pieces taken from them.
Neither is needed to read this document; the consolidating session may
want to drop the first or add it to `.gitignore`.

## 11. How the suite was run, and what that cost

`dev/resmooth-run-one.R` sets `NOT_CRAN=true` AND
`FRMTMB_BRMS_FIT_TESTS=true`, so the run is the GATED tier, not the
ordinary one. That is deliberate (the gated tier is where the brms
comparisons live, and this change is about matching brms) and it is why
two files dominate the wall clock: `test-brms-likelihood.R` and
`test-brms-methods.R` fit brms models, brms has no compile cache, and
every fit compiles a fresh Stan program at about a minute of `cc1plus`.
In main's own ungated log those two report 17 pass with 37 skip and 16
pass with 46 skip; here they report 428 pass 0 skip and 979 pass 0 skip.
`test-brms-likelihood.R` alone took about four hours of wall clock, most
of it in `cc1plus`, against a machine also running another lane's
suite.

Two more files skip under the flags the runner sets, and both were run
separately with their own flag rather than left as a skip:

    test-drmtmb-agreement.R  FRMTMB_DRMTMB_FIT_TESTS=true  131 pass
    test-fuzz.R              FRMTMB_FUZZ=true FRMTMB_FUZZ_N=60  2 pass

The fuzz tier was run at size 60 rather than its default 300, which is a
smoke run of the greedy cover and not the full plan. It is reported as
what it is.

## 12. The cost of differencing every smooth (nits round)

The review measured a wall-clock regression the first pass did not
report, up to 28.9x on `fitted()` for an ordinal `gp()` fit, and
attributed it to `re_b_batches()` returning `NULL` and abandoning
batching for the whole fit. **The number reproduces; the attributed cause
does not.** Counted, not reasoned about (`dev/resmooth-batchdiag.R`, log
`dev/resmooth-batchdiag.txt`): `re_b_batches()` returns `NULL` on none of
the three fits, and the bar blocks keep their batches throughout. On
`s(x, k = 8) + (1 | g)` at 40 levels it returns 7 batches, one covering
all 40 levels of the `us` block and six for the smooth's six columns.

The real cause is one line up. Batches are cut per column POSITION,
`for (k in seq_len(D))` with `D = bk[["dim"]]`, which is right for a
`(1 + x | g)` term whose row loads its level's intercept and its slope.
But for a smooth, `gp()` or `hsgp()` block `dim` IS the coefficient
count, so the split produces one batch per coefficient and saves nothing.

The instrument is the review's clock rule (10 ms tick, blocks grown past
1.2 s, minimum of three rounds, a no-smooth CONTROL carried) plus one it
did not use: the EXACT count of model evaluations the route makes, which
does not move with machine load. The clock is reported beside it because
the control moved by up to 1.36x between arms, so any wall-clock ratio
inside about 1.4 here is noise.

    cell                                  evals before/after   batches
    -----------------------------------   ------------------   -------
    A gp(x) 160 coefs, newdata 3 NEW pos     329 -> 329        160 -> 160
    B gp(x) 160 coefs, in sample             329 ->  11        160 ->   1
    C gp(x) 160 coefs, 3 observed rows        15 ->  11          3 ->   1
    D s(x,k=8)+(1|g) 40 lv, NULL               25 ->  25          7 ->   7
    E s(x,k=8)+(1|g) 40 lv, NA                 23 ->  23          6 ->   6
    F hsgp(x) k=12, in sample                  33 ->  33         12 ->  12
    CONTROL x+(1|g) 40 lv, NULL                11 ->  11          1 ->   1

Wall clock on three arms, one process each
(`dev/resmooth-cost3.R`, logs `dev/resmooth-cost3-{base,before,after}.txt`):

    cell                             0.64.0   before    after
    ------------------------------   ------   -------   -------
    A gp newdata 3 NEW positions     0.0497   1.9500    1.5600
    B gp in sample                   0.0195   1.2500    0.0202
    D s(x,k=8)+(1|g) NULL            0.0806   0.1725    0.1369
    E s(x,k=8)+(1|g) NA              0.0756   0.1444    0.1119
    CONTROL x+(1|g) NULL             0.0731   0.0644    0.0537

**The fix.** `re_b_batches()` now tries the whole kept set of a block as
ONE batch before falling back to the per-position split, and takes it
when no row loads two of the block's columns. An exact `gp()` in sample
has an indicator design, so that test passes and 160 batches become 1.
A `(1 | g)` block reaches the same single batch by a shorter route; a
`(1 + x | g)` block FAILS the test and falls through to the split
unchanged; a multi-membership block still returns `NULL`. Cell B's
regression is gone: 329 evaluations to 11, which is what 0.64.0 itself
spent, and 1.2500 s to 0.0202 s against 0.64.0's 0.0195 s.

**It is bit-neutral.** `dev/resmooth-batchident.R` saves the `Est.Error`
array of 44 cells on each arm, over the review's own constructions A, B,
C, D, F and I plus an `hsgp()`, a `gp()` at observed rows, a
`(1 + x | g)` fit, a multi-membership fit and a no-smooth control, at
`NA` and `NULL`, in sample and on newdata. All 44 are `identical()`, max
absolute difference 0.000e+00 (`dev/resmooth-batchident.txt`). A test
pins the batch count and the agreement with the unbatched route, and
fails on the pre-nits build with "Expected `bt` to have length 1. Actual
length: 160".

**What is NOT fixed, with its numbers.** Cell A, the review's headline,
is unchanged at 329 evaluations and about 31x over 0.64.0. Its three
newdata rows are OFF the fitted `gp()` positions, so each kriging row
loads all 160 columns and no batch of two columns can be attributed by
row. Cells D, E and F are unchanged for the same reason at 1.5 to 1.7x:
every row loads every basis column of a smooth or an `hsgp()`. Under
finite differences one evaluation pair per such column is irreducible.
Both are filed in `dev/test-backlog.md` with the route that would fix
them (the chain rule: a category probability depends on `b` only through
`eta`, which is linear in `b`, so `dp/db = (dp/deta) Z` needs one pair
per ROW rather than per coefficient, and `Z` is already built), and with
the one cheaper partial measure the round did not take: an `fs` block's
columns partition by LEVEL, so its 50 columns could be 5 batches rather
than 50, which the whole-block test cannot see.

### The suites after the nits round

Both were rerun in full against the build the nits round produced, one
file per process, with the same gated tier as the first pass. Generated
by `dev/resmooth-summarise.R` into `dev/resmooth-suite2-counts.txt`
(`RESMOOTH_SUITE_DIR=dev/resmooth-suite2`,
`RESMOOTH_EXT_DIR=dev/resmooth-ext2`) and pasted verbatim:

    == GENERATED BLOCK, paste verbatim ==
    core frmtmb suite: 179 files
      pass=14063 fail=0 error=0 skip=14 warn=17
      no file failed, errored or failed to report
      files with a skip: 2
        test-drmtmb-agreement.R.txt skip=13
        test-fuzz.R.txt skip=1
    frmtmb.sample + frmtmb.spline suites: 50 files
      pass=2623 fail=0 error=0 skip=2 warn=13
      no file failed, errored or failed to report
      files with a skip: 2
        frmtmb.sample-test-scale.R.txt skip=1
        frmtmb.spline-test-scale.R.txt skip=1
    == end ==

14063 against the first pass's 14040: the 23 are this round's new
assertions. The two gated brms files report the same 428 and 979 with no
skip. The two files that skip under the runner's flags were rerun with
their own gate against this build too: `test-drmtmb-agreement.R` 131
pass, `test-fuzz.R` 2 pass at size 60.

`R CMD check --as-cran` was run once, on the first pass, and is NOT rerun
here: the round's code change is `re_b_batches()`, which the check's own
`checking tests` step covers through the suite, and the Rd changes were
verified by rendering instead (`dev/resmooth-rdcheck.R` renders the six
touched Rd files and greps the RENDERED text, the only unescaped `%` in
any of them being roxygen's own header comment).
