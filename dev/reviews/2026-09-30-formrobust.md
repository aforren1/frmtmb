# Review of lane wt-formrobust, 2026-09-30

Adversarial review. Worktree `C:\Users\adf44\source\r\frmtmb-wt-formrobust`,
base `57c25589` (frmtmb 0.66.0), uncommitted. Lane build
`C:/Users/adf44/source/r/wt-formrobust-lib`; base build `rellib-r4`;
brms 2.23.0. Every script is `dev/formrobust-rev-*` and every log is in
`dev/formrobust-rev-log/`. No package file was edited. The Stan cache
was copied to the session scratchpad before any sampling.

## Verdict

**NOT MERGEABLE.** One blocker: the BREAKING `emmeans()` offset change
goes the wrong way. brms 2.23.0 run with emmeans includes a linear
predictor's offset at the grid's value, and the lane removes it, with
a test that pins the wrong number. The rest of the lane holds up well.
The addition-term expressions match brms's `standata()` bitwise on
every shape I tried. The ARMA fill equals brms's own `.predictor_arma()`
to 1.5 ulp at the same seed. 749 of 780 base-against-lane items are
identical, and the 10 that differ are all explained.

## The instrument

`formrobust-rev-verify-lib.R` sourced every R file of the worktree and
compared each function by deparse against the installed namespace:
frmtmb 1428 functions, 0 differ, 0 missing, NAMESPACE identical;
frmtmb.sample 235, 0, 0, identical (`verify-lib.txt`).

## BLOCKER 1. `emmeans()` drops the offset, and brms keeps it

The lane read `emm_basis.brmsfit()` (`posterior_linpred(offset =
FALSE)`) and concluded that brms leaves the offset out. It did not run
brms. brms's grid carries emmeans's own `.offset.` column, and emmeans
adds it back.

`formrobust-rev-brms-emm.R` and `-emm2.R` (data seed 21, Stan seed 7,
1 chain, 2000 iterations; logs `brms-emm.txt`, `brms-emm2.txt`),
`yc ~ x + f + offset(log(time))`, poisson:

| quantity | value |
|---|---|
| brms emmean, level a | 1.022449 |
| median of b0 + bx mean(x) | 0.3196746 |
| difference | 0.7027742 |
| log(mean(time)) | 0.7027742 |
| `ref_grid(fb)@grid$.offset.` | 0.7027742 |
| brms emmean at `time = 1` | 0.320 (the offset is log 1 = 0) |
| `posterior_linpred(offset = FALSE)` minus b0 | 0 (the basis alone does drop it) |

The same holds for `offset(time)` (emmean 0.859 = b0 + bx mean(x) +
mean(time)) and for `sigma ~ 1 + offset(log(time))` with
`dpar = "sigma"` (0.9003 = 0.1975 + 0.7028). An offset inside a
nonlinear parameter's formula is left out (nl emmean 1.021405 equals
the median of a + b mean(x) exactly), because emmeans cannot see it
there. `epred = TRUE` includes the offset once, at log(mean(time))
(2.779994 both ways), and the lane matches that.

frmtmb lane, same data (`formrobust-rev-offset.R`, `offset-after.txt`):
the emmean minus (b0 + bx mean(x)) is 0 for the poisson model, and the
emmean minus b_sigma is 0 for the sigma model. Both are 0.7027742 away
from brms. The nonlinear case and `epred = TRUE` agree with brms.

Why this blocks: this is a BREAKING bullet in NEWS that claims brms
parity. brms contradicts it, and it reverses emmeans's own convention
for glm, which brms inherits. `test-offset-grid.R` ("emmeans() leaves
the offset out of the linear predictor") asserts the wrong value to
64 ulp, and findings section 2 and the `?frmtmb-emmeans` text repeat
the claim.

The fix is one line, and I measured it.
`formrobust-rev-offset-fix.R` (`offset-fix.txt`, same data) keeps the
`offset` attribute in `emm_terms()` as a runtime mutant, and
`emm_drop_offsets()` stays as it is. emmeans then adds `.offset.` as
it does for brms, and all five cases agree with brms:

| case | result |
|---|---|
| poisson emmean minus (lin + log(mean(time))) | 0 0 |
| at `time = 1`, emmean minus lin | 0 0 |
| `epred = TRUE` over exp(lin + log(mean(time))) | 1 1 |
| sigma emmean minus (b_sigma + log(mean(time))) | 0 |
| nonlinear emmean minus (a + b mean(x)) | 0 (offset still left out, as brms) |

The same mutant (M05b) turns only two assertions of
`test-offset-grid.R` red: the emmeans value, and the
`e1 - (e0 - 0.1)` identity, which assumes the offset is left out.
The NEWS bullet, `?frmtmb-emmeans` and findings section 2 need the
same correction.

## Minors

1. **`refit()` of a recoded bernoulli fit with the original values
   fits garbage.** `refit(fn, d$yn)` on the `-1/-2` fit gives logLik
   3.048205e+304, with only two convergence warnings
   (`formrobust-rev-bern.R`, seed 51, `bern-after.txt`). The lane
   records this as not recoded, but the result is not a refusal. Code
   `newresp` with the stored `bin_levels`, or refuse values other than
   0 and 1.
2. **The fraction guard refuses a plausible correct model.** Effect
   coding `-0.5/0.5` is refused ("frmtmb refuses fractions"), while
   brms codes it `0/1` (`standata` Y equals the 0/1 column). With many
   distinct fractions brms refuses too ("only two different values").
   So the guard differs from brms only for exactly two fractional
   values, and `-0.5/0.5` is the common case there. **Needs the user:**
   keep the divergence, narrow it, or drop it. The all-`TRUE`
   divergence is right: brms codes an all-`TRUE` response as all
   failures.
3. **The `fitted()` and `posterior_epred()` divergence under an NA
   response is larger than the write-up shows, and NEWS and the compat
   row do not call it a divergence.** On frmtmb.sample draws
   (`formrobust-rev-sample-arma.R`, sampler seed 3, 2000 draws, the
   lane's construction), brms's `.predictor_arma()` run on the SAME
   draws gives an epred sd of 0.582 / 0.670 / 0.687 on rows 6 to 8.
   frmtmb.sample gives 0.097 / 0.103 / 0.122. The Laplace `fitted()`
   side is defensible, since it has no draws of the missing past. On
   draws, brms's reading, where the missing response is a latent that
   the epred draws integrate over, is coherent and easy to implement
   (reuse the draw closure). **Needs the user.** At minimum, NEWS and
   the compat `fitted` row should say that this diverges from brms.
4. **Test gaps (mutation testing).** 25 mutants were run against the
   lane's tests (`formrobust-rev-mut/`, `mut.log`; controls all
   green). 22 were caught. These 3 survive:
   - M05: `emm_drop_offsets()` made the identity. No test checks a
     grid-route offset value, only finiteness.
   - M07: the fill residual omits the MA term
     (`err = y - mu` in place of `y - mu - sma`). The arma test checks
     only `!anyNA`.
   - S01: frmtmb.sample fills with the expected value in place of a
     draw. The draws test checks only `!anyNA` and the rows before the
     first NA. My measurement separates the two cleanly: row 6
     predictive sd is 1.192 with the draw fill and 1.049 under the
     mutant, against 1.200 from brms on the same draws.

   The code itself is right on both ARMA points (see Notes). The tests
   would not notice if it regressed.
5. **An NA produced by an expression is not refused by name.**
   `trunc(lb = ifelse(x > 1, NA, -5))` reaches the optimizer and fails
   with "NA/NaN gradient". `trials(ifelse(x > 1, NA, c))` fails with
   "response must be integer counts". brms refuses all of these at
   `standata()`. Base dropped such rows through `na.omit`, so this is
   a behavior change (`formrobust-rev-aterm2.R` section D, seed 101).
6. **Formula-environment constants.** `weights(wt * k)` and
   `trials(k + 0)` with a scalar `k` are accepted, but `trials(k)`
   (bare) is refused. brms refuses all three ("can neither be found in
   'data' nor in 'data2'"). Findings section 1 describes brms as
   reading such constants. The extension is harmless. The bare-name
   inconsistency is not (`formrobust-rev-refusals.R`, seed 91).
7. **Wrong brms statements in findings section 1.** brms does not
   recycle a single value for `cens()`'s interval bound ("Argument 'y2'
   needs to have length equal to the number of data rows"), and
   neither does frmtmb. frmtmb's message is misleading, though: "upper
   bounds must not be NA on interval-censored rows"
   (`aterm2-after.txt` section C).
8. `weights(t * 2)` with no column `t` fails with "non-numeric argument
   to binary operator": `t` resolves to `base::t`. The variable is not
   named.
9. A factor inside an offset, `offset(log(as.numeric(ef)))`, lets
   "variable 'ef' is not a factor" escape from
   `conditional_effects()` and `emmeans(epred = TRUE)`
   (`offset-after.txt`).
10. A pooled `update()` stores the evaluated bf() object in `fit$call`.
    It deparses to about 3 KB of family closures
    (`update-after.txt`).
11. House style: `vignettes/brms-migration.Rmd` has an added line of
    89 columns. The frmtmb.sample NEWS names "0.66.0", and the round
    rule says no version numbers. In `test-update-pool.R`,
    `expect_message(suppressMessages(...), NA)` is vacuous. In
    `test-formrobust-draws.R`, `frm_sample()` is wrapped in
    `suppressWarnings()` and not in `allow_warnings()`.
12. `update()` on a multivariate fit with a complete formula still
    replaces everything and drops `sigma ~ z`, where brms refuses
    ("Updating formulas of multivariate models is not yet possible").
    This is unchanged from base, but it is now inconsistent with the
    univariate BREAKING change.

## Notes: verified, no action

- **Addition-term expressions** (`formrobust-rev-aterm.R`,
  `-aterm2.R`, seed 101). brms `standata()` is bitwise equal on 12
  shapes: `weights(wt*2)`, `weights(wt, scale = TRUE)`, `se(s/2)`,
  `rate(t*2)`, `trials(c+1)`, `trials(max(c)+1)`,
  `trunc(lb = min(y)-1)`, `trunc(ub = max(y)*2)`,
  `weights(log(wt)+1)`, `weights(wt^2)`, `se(sqrt(s))` and
  `weights(1 + 0*wt)`. Variables named `t` and `c` work.
  Multivariate with `subset()` on one response and expressions on
  either response is bitwise equal to brms (`lb_y`, `weights_y`,
  `weights_y2`, `se_y2`, and scale = TRUE under a subset), with N_y 30
  and N_y2 60 as brms gives when NA falls outside the subset. `cens()`
  with an interval expression, with and without NA rows, and
  `index(id*3)` are identical to the precomputed columns. Readers of
  the slimmer frame are identical to the precomputed-column fit:
  residuals, influence (numerics), refit, bootstrap, predict(newdata),
  conditional_effects, update(newdata), and frmtmb.sample draws,
  log_lik, loo, posterior_predict and posterior_epred
  (`formrobust-rev-sample-aterm.R`). The one exception is
  `trunc(min(y))` on newdata, which brms also evaluates on newdata's
  rows.
- **ARMA fill against brms** (`formrobust-rev-arma.R`, seeds 77 to
  80). At the frmtmb MLE, with brms's own `.predictor_arma()` and a
  minimal brmsprep, one group at the same seed: 20 cases (ar1, ar2,
  arma11, ma1, student; four missing patterns each: rows 5 to 8, the
  first row and row 3, all rows, and alternate rows). 16 are identical
  and 4 are within 1.5 ulp. Three groups with 40000 draws each: KS
  p >= 0.129 on every row with spread, sd ratio 0.9865 to 1.0060, and
  fitted() Estimate within |z| 2.40 of brms's mean over fills. On
  frmtmb.sample draws against brms on the same draws: predictive sd
  ratio 0.984 to 1.017, KS p >= 0.150. The lane's worry about 1.30
  against 1.22 was Monte Carlo noise in brms's 1000 draws.
- **update() pooling** (`formrobust-rev-update.R`, seed 41).
  `bf(y ~ x2, sigma ~ w)` and `sigma ~ 1` replace, with brms's message.
  A plain formula keeps `sigma ~ z`, and its logLik equals the direct
  fit bitwise. Family order is as brms's. A dpar constant is pooled.
  The nonlinear body replacement keeps `a` and `b`.
- **Exports** (`formrobust-rev-exports.R`). Attaching frmtmb masks
  `stats::ar` with the usual message. No R code in core or any
  extension calls `ar()`. `autocor()` works with brms unloaded, after
  `loadNamespace("brms")`, and with brms attached. `bf(autocor =)` and
  `acformula()` fits equal the in-formula fits bitwise, also for
  `mvbind`. `drop_unused_levels = FALSE` gives identical logLik, and
  NA with a warning on the unused level where brms predicts from the
  prior (a documented divergence).
- **No regression** (`formrobust-rev-regress.R`, `-regress-cmp.R`,
  seed 71). 39 models (gaussian to mixture, mv, rescor, subset, mi,
  ARMA both forms, nl, smooth, mo, cens, trunc, se, rate, trials,
  weights, bernoulli 0/1) and 20 items each: fn, gr at the optimum and
  at a shifted point, parameters, logLik, fixef, fitted (in sample and
  newdata), predict (seed), residuals, simulate (seed),
  conditional_effects, emmeans three ways, and default_prior. 749 are
  identical and 21 fail identically in both builds. The 10 that differ
  are all on the two offset models: 8 are base errors that now answer,
  and 2 are `names(fit$data)` gaining `time`.
- **Base failures of the new tests** are behavioral for 6 of 7 files
  (lane logs `before-dbg-*`). The weak form is limited to
  `autocor()`/`acformula()` not existing in `test-brms-api-formrobust.R`.
  The tests have no absolute tolerances: the arma bound is on a
  dimensionless sd ratio.

## Counts (generated by `formrobust-rev-counts.R` from the per-file logs)

Ungated, NOT_CRAN=true, one file per process, package attached
(`suite-counts.md`). Every one of the 334 files loaded frmtmb from
`wt-formrobust-lib`:

| package | files | with RESULT | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| frmtmb | 199 | 199 | 13646 | 0 | 0 | 163 | 0 |
| frmtmb.coupling | 11 | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 29 | 29 | 1743 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 10 | 10 | 359 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 15 | 429 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 11 | 547 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 44 | 44 | 2223 | 0 | 0 | 4 | 0 |
| frmtmb.spline | 15 | 15 | 553 | 0 | 0 | 1 | 0 |
| total | 334 | 334 | 20042 | 0 | 0 | 192 | 0 |

Gated tier (`gated-counts.md` plus `gated2.log` for the two files that
also need `FRMTMB_DRMTMB_FIT_TESTS` and `FRMTMB_FUZZ`): frmtmb 31 files,
2802 pass, 8 fail, 0 err, 0 skip; frmtmb.sample 5 files, 201 pass,
0 fail. The 8 failures are exactly the stale "now HOLDS" rows of
findings section 7: brm:112, brmsfit-methods:112, :113 and :747, and
standata:75, :112, :970 and :974. These match the lane's counts.

The lane's R CMD check logs agree with its claims: frmtmb 1 NOTE (HTML
manual), frmtmb.sample OK. I did not rerun check.

## For the user

- Blocker 1 needs a fix before merge. The direction is above.
- Decide minors 2 and 3: the fraction guard, and whether
  frmtmb.sample's epred should integrate the fill as brms does.
- Likely consolidation overlaps: `R/frame.R`, `R/parse.R`,
  `R/conditional-effects.R`, `R/interop.R`, `R/families.R`, `NEWS.md`.

## Re-check after punch round 1

I read findings section 12 and the diff. I had no snapshot of the
reviewed tree, so I found the change set by modification time since
this review. It touches 10 R files in core (`autocor.R`, `bf.R`,
`compat.R`, `confint.R`, `families.R`, `frame.R`, `interop.R`,
`predict.R`, `sampling-api.R`, `sugar.R`), frmtmb.sample's
`methods-draws.R`, 7 test files, NEWS and the vignette. The lane
library again matches the source: 1431 core functions and 235
frmtmb.sample functions, 0 differ (`rev2-verify-lib.txt`).

### Verdict: NOT MERGEABLE

The original blocker is fixed, and I measured the fix against brms on
every shape I could build. Punch round 1 introduced one new blocker:
a pooled `update()` silently resets a family option.

### NEW BLOCKER. A pooled `update()` drops a family's non-link options

`bform_call()` writes the pooled formula into the call, and
`family_call_of()` turns the family into a constructor call. The
constructor call is written whenever it returns the same family name
and links. That test cannot see an option that is not a link.

`huber(k = 3)` inside `bf()` becomes `frmtmb::huber(link = "identity",
link_sigma = "log")` in the call. The update evaluates that call, so
the refit uses the default `k = 1.345`, with no word
(`formrobust-rev2-famcall.R`, seed 95, `rev2-famcall.txt`):

| fit | logLik |
|---|---|
| `update(fit, y ~ x + z)` | -231.49610424 |
| direct fit, default k | -231.49610424 (identical) |
| direct fit, k = 3 | -237.749592515 |

A family given as `frm(family =)` is not affected, because the call
keeps the original expression. Other constructors with options besides
links are exposed the same way: `whittle(tapers =)` and
`cox(df =, degree =, intercept =)`. `categorical(levels =)` falls back
to the stored object and is safe (`rev2-famargs.txt`).

The fix: write the constructor call only when every formal of the
constructor is a link argument, and store the object otherwise. Then
add a `huber(k = 3)` case to `test-update-pool.R`. No current test
catches this: the mutant R08, which always stores the object, is
caught only by the assertion on the call's text.

### 1. emmeans() offsets: fixed and brms-exact

Both brms scripts reran with identical output (`rev2-brms-emm*.txt`).
The lane's own five cases are exact (`rev2-offset-lane.txt`). For the
new shapes, brms ran at Fixed_param draws equal to frmtmb's estimates
(the test helper `brms_fixed_fit()`, my Stan cache copy), in
`formrobust-rev2-offset.R`, `-offset-mv.R` and `-offset-nl.R` (data
seed 21). All of these agree with brms to at most 0.9 ulp:

- **Offset plus factors:** `~ f` with `by = "h"`, `~ f | h`,
  `at = list(time = c(1, 2))` averaged, `~ f | time` at `c(1, 2)`,
  `epred` at `c(1, 2)`, and `~ x` at `x = c(-1, 1)`. The hand formula
  lin + mean(log(c(1, 2))) also gives 0 0.
- **Two offsets** (`offset(log(time)) + offset(z)`), link and epred;
  the hand formula gives 0 0.
- **An offset in mu and one in sigma:** mu, mu epred, `dpar = "sigma"`,
  and sigma at `time = c(1, 2)`.
- **Nonlinear:** the body, `nlpar = "a"` (whose offset is included, as
  brms includes it), `nlpar = "b"`, epred, and `nlpar = "a"` at
  `z = 0`.
- **Multivariate, one response selected:** `resp = "yc"`,
  `resp = "y2"`, `resp = "yc"` at `time = 1`, `resp = "yc2"`, and
  epred.

One multivariate case diverges, and frmtmb is right. With one family
and the offset on one response only, emmeans over both responses in
frmtmb gives each response its own offset. brms has one `.offset.`
column for the whole grid, so it adds yc's offset to yc2 as well:

| case | brms minus frmtmb |
|---|---|
| default grid | log(mean(time)) / 2 = 0.3514 (1.545 against 1.194) |
| `~ f \| time` at `time = 2` | log(2) / 2 |
| `~ f \| time` at `time = 1` | 0 |

frmtmb's value is the correct one; brms's is a defect. **Minor:** NEWS
and `?frmtmb-emmeans` do not record this as a divergence from brms.
The mutant R02, which removes the grid-route switch, is caught by
`test-offset-grid.R`. With mixed families and no `resp`, both packages
refuse.

### 2. ARMA epred fill on draws: brms-exact by value

`formrobust-rev2-sample-epred.R` compares frmtmb.sample with brms's
`.predictor_arma()` on the same 2000 posterior draws (sampler seed 3,
RNG seed 2). brms runs one draw at a time, which is frmtmb.sample's
RNG order.

- **By value:** `posterior_epred(newdata)` and
  `posterior_predict(newdata)` are both `identical()` to brms, with a
  maximum difference of 0.
- **By law** (brms on all draws at once, a different RNG order): the
  epred sd on rows 6 to 8 is 0.568 / 0.678 / 0.685 against brms's
  0.579 / 0.668 / 0.708, with KS p of 0.436 / 0.075 / 0.460.
- The rows before the first NA are identical to the epred with the
  full response.

The mutants S01 and R01 (expected fill, fill off) are now caught.

### 3. New refusals and their guard-absent cases

`formrobust-rev2-misc.R`, seed 93, `rev2-misc.txt`.

**NA from an expression (BREAKING).**
- **Refused:** `weights(ifelse(x > 1.5, NA, w))` ("NA on 6 of 60
  rows").
- **Guard absent:** an NA in the response, in a plain predictor, or in
  the expression's own variable still drops the rows through
  `na.omit`. nobs is 58 in each case, and the logLik is identical to
  the precomputed-column fit.
- **Guard absent, other shapes:** an NA only on rows outside a
  `subset()`, where brms also accepts it; an `mi()` model with an NA in
  the expression's variable; an `mi()` response with `weights(w * 2)`;
  and an interval bound that is NA on non-interval rows. All fit.

**`refit()`.**
- **Refused:** the value 2, and `"yes"/"no"` on a character-response
  fit, each by name.
- **Guard absent:** 0/1, logical, 0/1 on the character-response fit, a
  gaussian refit with any values, and binomial counts all fit.
- A `newresp` with NA passes the guard and ends in an optimizer
  failure. That is a note only; I did not check whether base behaves
  the same.

**Single-bound `cens()`.**
- **Refused:** `cens(cc, 10)` and `cens(cc, max(y) + 10)`, with the
  length stated, as brms requires.
- **Guard absent:** `cens(cc, ub)` as a column, and `cens(cc)` with no
  interval rows.

**`weights(t * 2)` with no column `t`.**
- **Refused, by name.** `offset(log(t))` is refused the same way.
- **Guard absent:** a column `t`, and `t` as a numeric vector of the
  formula environment.
- The guard also fires on `weights(vapply(w, round, 1))`, a function
  passed as an argument. brms refuses that too ("can neither be
  found"), so this is a note, not a false alarm against brms.

**Multivariate `update()`.**
- **Refused:** a complete formula, and a univariate fit updated with
  an `mvbf()`, both in brms's words.
- **Guard absent:** `update(newdata =)` and `update(control =)` on a
  multivariate fit both work. The delta refusal predates this lane.

Each new refusal is caught by its mutant: R03 (cens), R04 (function
variable), R05 (multivariate update), and P03 and P04.

### 4. Formula-environment constants

The worker's claim holds (`rev2-misc.txt`, `rev2-envinfl.txt`).
`I(x * ck)`, `offset(z * ck)` and `sigma ~ I(z * ck)` all read `ck`
from the environment. `trials(k)` equals `trials(12)` bitwise.

After `k` or `ck` changes, `trials(k)`, `I()`, `offset()` and the dpar
formula all behave the same way:
- `fitted(newdata =)` moves with the new value.
- `update()` refits with the new value.
- `influence()` refits re-read it. This holds for `offset(x * ck)` as
  for `trials(k)`.
- In-sample `fitted()`, `residuals()` and `simulate()` keep the
  stored values.

That is R's convention, and it applies alike to addition terms,
predictors and offsets. The divergence from brms is documented in
`?bf` and NEWS.

A pre-existing defect, not this lane's: `influence()` of
`y ~ I(x * ck) + (1 | g)` fails on base and lane alike ("The model
uses `x`, which is not a column of `data`"). The refits read the model
frame, which holds `I(x * ck)` and not `x` (`rev2-envinfl*.txt`).

### 5. `eval(fit$call)` after a pooled update

`eval(u$call)` gives an identical logLik and identical coefficient
names in all 6 cases: a dpar formula, `student(link_sigma =
"softplus")`, `Gamma(link = "identity")`, a dpar constant, a
nonlinear model, and the `huber` case. The stored call reads as, for
example, `frmtmb::bf(y ~ x + z, sigma ~ z, family =
stats::gaussian(link = "identity"))`.

This agreement holds by construction, since the update itself
evaluates that call. That is exactly why it cannot reveal the blocker
above: the call and the fit are wrong together.

### 6. Mutants and counts

The run is `formrobust-rev-run-mut2.sh`, logged in `mut2.log`, with
`FRMTMB_BRMS_FIT_TESTS=true`. All 8 control files are green.

- **The review's 25 mutants:** 23 still apply, and all 23 are caught.
  S01 is now caught. M07 is caught by the new hand-recursion test. M05
  is caught, but by a signature error. Two no longer apply, and the
  helper refuses them by name: M05b (the attribute strip is gone; P01
  is its replacement) and M19 (the fraction guard became a warning;
  P02 covers it).
- **The worker's P01 to P06:** all caught, with failure counts equal
  to the worker's `p1-mut.log`.
- **New mutants R01 to R08 for the round-1 code**, each run on all 8
  lane test files: all caught. R01 (epred fill off) by
  `test-formrobust-draws.R`. R02 (multivariate grid switch) by
  `test-offset-grid.R`. R03 (cens length), R04 (function variable) and
  R07 (bare constant sent to the frame) by `test-aterm-expr.R`. R05
  (multivariate update allowed), R06 (object stored in the call) and
  R08 (family object stored) by `test-update-pool.R`.

The ungated suite, which I reran with one file per process
(`suite2-counts.md`), equals the worker's `counts-p1.md` exactly.
Every file loaded frmtmb from `wt-formrobust-lib`.

| package | files | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 199 | 13682 | 0 | 0 | 164 | 0 |
| frmtmb.coupling | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 29 | 1743 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 10 | 359 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 429 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 547 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 44 | 2224 | 0 | 0 | 4 | 0 |
| frmtmb.spline | 15 | 553 | 0 | 0 | 1 | 0 |

The one new core skip is the gated brms comparison in
`test-offset-grid.R`; in the mutation controls, with the gate on, that
file passes 31 with 0 skipped. I spot-checked the worker's gated
logs: 37 RESULT lines. `test-offset-grid.R` passes 31 with 0 skipped,
and `test-brms-suite-standata.R` has the 4 stale rows as before.

The 39-model regression battery reran against base with the same
result as round 0: 749 identical, 21 fail in both builds, and 10
differ, all explained, all on the offset models
(`rev2-regress-cmp.txt`).

### For the user

- The new blocker needs the family-call fix above and a test.
- Recording the multivariate emmeans offset divergence from brms in
  NEWS and `?frmtmb-emmeans` is a minor.
