# Lane formrobust: failures on correct brms code in formula handling

Base: 57c25589 (frmtmb 0.66.0). Worktree `frmtmb-wt-formrobust`, branch
`wt-formrobust`, library `C:/Users/adf44/source/r/wt-formrobust-lib`.
"Before" is the read-only base build `rellib-r4`; every repro script
takes `FORMROBUST_LIB=""` for that arm. brms is 2.23.0. The brms
functions this lane read are printed into `dev/formrobust-brms-src/` by
`dev/formrobust-brms-src.R`.

## 1. Addition terms given an expression

**brms.** `get_ad_values()` evaluates the term's argument with
`eval2(expr, data)` on the data rows after `subset_data()`. The model
frame holds only `all_vars()` of the term, and `validate_data()`
refuses a variable that is not in `data` or `data2` ("can neither be
found in 'data' nor in 'data2'"), a constant of the environment
included. `data_response()` recycles a single value for `trials()`,
the censoring code of `cens()`, `trunc()`, `vint()` and `vreal()`, but
not for the interval bound `y2` of `cens()` ("Argument 'y2' needs to
have length equal to the number of data rows"). It scales
`weights(x, scale = TRUE)` as `weights / sum(weights) * length(weights)`.
(Corrected in punch round 1: this paragraph said `cens()` recycles and
that brms reads environment constants; see section 12.)

**Defect.** `assemble_frame()` put the term's expression into the
model-frame formula as a term. There `wt * 2` is an interaction and
`s / 2` a nesting, so `weights(wt * 2)`, `rate(time * 2)`, `se(s / 2)`
and any multivariate model holding one of them stopped with "invalid
model formula in ExtractVars". `trunc(lb = min(y) - 1)` stopped with
"variable lengths differ (found for 'min(y)')", and `weights(wt * k)`
with `k` a constant of the formula environment stopped the same way.
`trials(n + 1)`, `cens(c == 1)`, `trunc(lb = lb - 1)` and
`subset(x > 0)` worked by accident: their expressions happen to be
valid formula terms. `weights(w, scale = TRUE)` was refused ("takes
exactly one argument"). Before and after:
`dev/formrobust-log/repro1-before.txt`, `repro1b-before.txt`,
`repro1b-after.txt` (`dev/formrobust-repro1.R` seed 1,
`dev/formrobust-repro1b.R` seed 11).

**Change.** A call-valued term puts only its variables in the frame
(`expr_frame_vars()` in `R/frame.R`): a column of `data`, or an object
of the formula environment with one value per row. A constant is read
when the term is evaluated. The term is evaluated on the response's
rows, a single value is recycled, and a value of any other length is
refused by name. A literal keeps its old path. Since punch round 1 an
expression that gives `NA` is refused by name, and a single-value `y2`
of `cens()` is refused as brms refuses it (section 12). Since round 3,
by the user's decision, every variable an addition term reads must be
a column of the data, as in brms, and one found only in the formula
environment (a scalar `k`, a full-length vector, or `t` found only as
the function `t()`) is refused by name (section 14). This supersedes
the round-1 reading of environment constants.
`weights(w, scale = TRUE)` is parsed into the expression
`w / base::sum(w) * base::length(w)`, so the frame, a refit on a subset
and every reader of the value take brms's scaling on the response's
rows. The newdata paths (`aterms_for_newdata()`, `ce_aterms()`)
already evaluated the expression on newdata; they did not change.

**Evidence.** In `dev/formrobust-repro1b.R` (seed 11), each expression
term and its precomputed column give identical `logLik()`, `fitted()`
and `predict()` on newdata, and identical `simulate(newdata = )`. This
holds for `weights`, `rate`, `trials`, `se`, `se(sigma = TRUE)`,
`trunc`, a multivariate `weights` + `se` model, a multivariate
`rate` + `subset` model, `weights` + `mi()`, a constant `k` from the
environment, and `index(id * 3)` with `mi(x, idx = )`. brms's
`make_standata()` gives the same `weights`, `se`, `denom` and `lb`
vectors bitwise (`test-aterm-expr.R`, "the values are brms's own").

`mi(x, idx = )` in a predictor takes one variable name in brms
("'mi' only accepts single untransformed variables"), and frmtmb keeps
that refusal. `index()` evaluates its expression, as brms does.

## 2. `offset()` in `conditional_effects()` and `emmeans()`

**brms.** `prepare_conditions()` holds every variable of
`bterms$allvars`, `time` included, at its mean. `posterior_epred()` at
that grid then includes the offset at `log(mean(time))`.
`get_all_effects()` does not list an offset's variable as a display
(`dev/formrobust-log/brms-effects.txt`). `emm_basis.brmsfit()` takes
`posterior_linpred(offset = FALSE)` on the link scale, but the terms
its `recover_data()` hands emmeans carry the offset, so emmeans puts it
in the grid's `.offset.` column and adds it back at the grid's value.
brms's emmean therefore INCLUDES a predictor's offset, at
`log(mean(time))` or at the `at =` value, and leaves out only an offset
inside a nonlinear parameter's formula, which its grid cannot see. With
`epred = TRUE` it includes every offset once. (The round-0 text of this
section read the basis alone and said brms leaves the offset out. That
was wrong, and the review measured brms: section 12.)

**Defect.** The model frame held the column `offset(log(time))` and not
`time`. Both grids are built from the frame's columns, so neither had
`time`. `log(time)` then found `stats::time()`, and
`conditional_effects()` stopped with "non-numeric argument to
mathematical function" (`offset(time)` gave "invalid type (closure)").
`emmeans()` stopped with "undefined columns selected", and on the nl
model with "Variable 'z' missing from newdata"
(`dev/formrobust-log/repro2-before.txt`, `dev/formrobust-repro2.R`
seed 21).

**Change.** The frame also holds the variables of every `offset()`
term of every linear predictor (`offset_vars()`). An offset's variable
is not a default display (`ce_lp_vars(rsv = FALSE)`). After punch round
1, `emmeans()` keeps the `offset` attribute in the design route's terms,
so emmeans adds `.offset.` as it does for brms. On the grid route the
selected predictor's prediction carries its own offset, and
`emm_drop_offsets()` removes only the offsets of the other predictors,
that is, of the parameters a nonlinear body reaches. Both responses of
a multivariate fit on one grid take the grid route, so that each gets
its own offset rather than one shared `.offset.` column. Against
rellib-r4 this is a bug fix, not a breaking change: every one of these
calls errored before.

**Evidence** (`dev/formrobust-repro2b.R`, seed 21):
`conditional_effects()` estimate / `exp(eta + log(mean(time)))` = 1.
The default displays are `x`, `f` (`dev/formrobust-repro2c.R`). All of
`offset(time)`, `offset(log(time) + 0.1)`, `sigma ~ offset(log(time))`
and a nonlinear parameter's `a ~ f + offset(z)` answer both calls
(`dev/formrobust-log/repro2-after.txt`). The emmeans values after punch
round 1 are in section 12.

## 3. `predict(newdata = )` with an NA response under `cov = FALSE`

**brms.** `validate_newdata(check_response = FALSE)` fills a missing
response column with `NA`. `.predictor_arma()` walks the rows in order.
For each row it adds the MA and AR terms of the earlier residuals, and
when `Y[n]` is `NA` it sets `y` to `posterior_predict_<family>(n)` at
the current shifted mean. It then stores `err = y - eta_before_ar`.
`posterior_predict()` then draws every row again at the shifted mean.
So a fill is a draw per posterior draw, and `posterior_epred()` carries
the fill's randomness too.

**Change.** `autocor_cond_newdata_y()` returns `NA` for a missing value
and all `NA` for an absent column. It refuses a non-numeric response
only. `autocor_cond_mu_fill()` is the recursion on plain numerics with
a fill. `predict()` fills with a draw from the family at the shifted
mean, using the replicate's own random stream
(`arma_cond_fill_dpars()`, exported for frmtmb.sample, whose
`posterior_predict()` does the same per draw). `fitted()` and
`frm_linpred()` fill with the expected value, the shifted mean itself.
For gaussian and student, the only families `cov = FALSE` takes, the
shift is linear in the residuals, so this is the expectation of brms's
filled recursion at fixed parameters. With no response missing,
`fitted()` keeps the old `autocor_cond_mu()` path, to the bit.

frmtmb.sample refused `posterior_predict(newdata = )` of every
`cov = FALSE` model as a structured draw. Core's `predict()` already
exempted it. The exemption now reaches the draws too, for newdata and
for `re_formula`.

**Evidence.** `dev/formrobust-repro3.R` (seed 31, 30 groups of 8,
AR(1) at 0.6, `ar` estimated at 0.543, `sigma` at 1.024, 4000 draws,
`propagate_error = FALSE`). Row 5 of group 1 reads observed residuals:
Est.Error / sigma = 1.0115. Row 6 reads row 5's fill:
Est.Error / (sigma sqrt(1 + ar^2)) = 1.00315, and Est.Error / sigma
= 1.1415, where an expected-value fill would give 1. The `fitted()` of
row 6 / the hand recursion = 1. On draws (`dev/formrobust-repro3c.R`,
2000 draws) the sd of row 6 over row 5 is 1.0937.

brms itself (`dev/formrobust-brms-arma.R`, one chain of 2000,
`dev/formrobust-log/brms-arma.txt`), rows 5 to 8 of group 1:

| quantity | brms | frmtmb |
|---|---|---|
| predict Est.Error | 1.0367 1.1896 1.1825 1.3003 | 1.0222 1.1751 1.2223 1.2206 |
| fitted Estimate | 1.99687 0.98368 1.43549 1.72239 | 1.99637 0.94744 1.40844 1.69977 |
| fitted Est.Error | 0.082312 0.589030 0.645964 0.705964 | 0.081575 0.093792 0.101448 0.121045 |

The predictive spreads agree to Monte Carlo error. The fitted means
agree within about two Monte Carlo standard errors of brms's randomly
filled mean (0.59 / sqrt(1000) = 0.019). **Deliberate divergence:**
brms's `fitted()` Est.Error on a row after a missing response includes
the spread of the random fill (0.59 against 0.094 on row 6).
frmtmb's `fitted()` is the conditional expectation given the observed
rows, and its Est.Error is parameter uncertainty alone. brms mixes
the predictive noise of an unobserved response into an expected value.
The reason for the divergence, as the coordinator's decision of punch
round 1 words it: a maximum likelihood fit has no draws to propagate the
fill through. It is documented in NEWS, `?fitted.frmtmb_fit` and the
compat `fitted` row. On draws there is no such reason, and since punch
round 1 frmtmb.sample's `posterior_epred()` fills each draw with a
sampled value, as brms does (`arma_cond_fill_epred()`; section 12).

## 4. `update()` keeps the stored parameter formulas

**brms** (`update.brmsfit()`, `update.brmsformula()`,
`dev/formrobust-brms-update.R`, output
`dev/formrobust-log/brms-update.txt`). A complete formula or `bf()`
replaces the location formula. `flist = c(object$pforms, object$pfix,
formula.$pforms, formula.$pfix)` pools the parameter formulas, and a
later one replaces an earlier one with "Replacing initial definitions
of parameters". On a nonlinear fit the mode is "replace", with the
message "Argument 'formula.' will completely replace the original
formula in non-linear models", and `nl` is taken from the fit when the
formula is plain. `bf(count ~ a + b, nl = TRUE)` is a valid brms object.
"No non-linear parameters specified" is raised at `brm()`.

**Defect.** frmtmb's `bf()` refused `nl = TRUE` without a parameter
formula at the call (ledger row `brmsfit-methods:955`). A complete
formula in `update()` replaced the whole `bf()`, so
`update(fit, y ~ x2)` on `bf(y ~ x, sigma ~ z)` lost `sigma ~ z`
(`dev/formrobust-log/repro45-before.txt`: coefficients `Intercept x z`).

**Change.** `bf()` no longer refuses. `parse_spec()`'s existing
refusal ("A nonlinear formula needs at least one nonlinear-parameter
formula") fires in `frm()`, and `bf(...) + lf(a ~ 1)` now works.
`update_complete_bform()` pools a complete formula with the stored
parameter formulas, constants, equations and `nlf()` formulas, with
brms's two messages and brms's family order (the formula's family, then
the `family` argument, then the stored one). A model with nothing to
pool keeps the old path, so its call keeps the formula as written.
**BREAKING** for a complete formula on a model with dpar formulas.

## 5. `bernoulli()` on two values that are not 0 and 1

**brms** (`data_response.brmsframe()`,
`dev/formrobust-log/brms-binary.txt` from
`dev/formrobust-brms-binary.R`). `bin_levels` is
`frame$basis$resp_levels`, stored with the fit, else
`levels(as.factor(Y))`. A numeric response with one value becomes
`c(0, value)`, or `c(0, 1)` when that value is 0. Then
`Y = as.integer(factor(Y, levels)) - 1`, and a third value is refused
("requires responses to contain only two different values"). Measured:
`-1, -2` codes 1, 0; `1, 2` codes 0, 1; all 5 codes 1; `"yes"/"no"` and
`TRUE/FALSE` code by level; all-`TRUE` codes 0.

**Change.** `extract_y()` codes a bernoulli response with
`bernoulli_levels()` and `bernoulli_code()`. The coding goes onto the
response's family (`bin_levels`), and `carry_finalized_responses()`
carries it into the fit's spec. So `influence()` refits, which assemble
from `model$spec`, code against the fit's levels. `frm_bootstrap()`
refits take the simulated 0/1 codes as they are. A response read from
newdata is coded with the same levels (`response_codes_newdata()`) in
`residuals(newdata = )`, `pp_check(newdata = )` and frmtmb.sample's
`predictive_error(newdata = )`. `simulate()` and the draws give 0/1,
as brms's `posterior_predict()` does. The points of
`conditional_effects()` are the raw response, as brms's
`make_point_frame()` plots them.

**Evidence.** `test-bernoulli-coding.R` (seed 5): `-1/-2`, a factor, a
logical, `1/2` and a character response give a `logLik()` and
`fixef()` identical to the 0/1 fit. The same holds for `fitted()`,
`residuals()`, `predict()`, `simulate()`, `AIC()`, `influence()`,
`frm_bootstrap()` and `conditional_effects()`. brms's `standata()` `Y`
equals frmtmb's codes for six types. A newdata holding only `-2`, which
coded afresh would read as the 1 of `c(0, -2)`, gives the 0/1 fit's
residuals. frmtmb.sample (`test-formrobust-draws.R`, sampler seed 3):
the draws, `log_lik()`, `posterior_predict()` and
`predictive_error(newdata = )` on the `-1/-2` fit are identical to the
0/1 fit's.

**Deliberate divergences.** (a) An all-`TRUE` logical codes 1, where
brms codes it 0 because `as.factor(TRUE)` has one level; a logical is
read as its number. (b) Withdrawn in punch round 1: two fractions
were refused in round 0. They are now coded as brms codes them, with a
warning when every value lies strictly inside (0, 1), and effect coding
such as -0.5 and 0.5 does not warn. `test-open-issues.R` keeps the
lme4#682 guard as an asserted warning.

## 6. Small API items

- **`0 + intercept`** (rows `standata:970`, `:974`). brms's
  `data_rsv_intercept()` warns "Reserved variable name 'intercept' is
  deprecated. Please use 'Intercept' instead.", fills
  `data$intercept` with ones and keeps the column name `intercept`.
  frmtmb now does the same: the warning, a column of ones in the data
  (`rsv_lower_fill()`) and in newdata (`pred_design()`), and brms's
  refusal of a non-ones column. The design equals brms's `standata()`
  `X` bitwise. The column is not a default display of
  `conditional_effects()`.
- **`frm(drop_unused_levels = )`** (row `standata:1133`) reaches
  `model.frame()` and is kept for `influence()` refits. With `FALSE`,
  an unused level of a predictor gives a zero column, which frmtmb
  drops as rank deficient. brms keeps it for its prior, so the row
  stays a divergence. A grouping factor loses the unused level either
  way (`dev/formrobust-dul-probe.R`: `ranef()` lists a, b, c for both,
  and the logLik relative difference is 0).
- **Autocorrelation terms as objects** (row `brm:112`). brms's
  `arma()`, `ar()`, `ma()`, `cosy()` and `unstr()` return an `ac_term`.
  `bf() + <ac_term>` is refused: "Autocorrelation terms can only be
  specified on the right-hand side of a formula, not added to a
  'brmsformula' object." The row asserts exactly that refusal. The
  task brief said brms "allows" adding the object; brms's source says
  otherwise. What brms allows outside the formula is `bf(autocor = )`
  and `+ acformula()`, and frmtmb now reads both. It adds the terms to
  the location formula, which is brms's reading, and refuses an
  `acformula()` without `resp` on a multivariate formula, as brms does.
  Exporting `ar` masks `stats::ar` when frmtmb is attached, as brms's
  does.
- **`autocor()`** (rows `brmsfit-methods:112`, `:113`) returns `NULL`
  with brms's "Method 'autocor' is deprecated and will be removed in
  the future.", which is brms's answer for every fit since 2.11. It
  validates `resp` with the shared `brms_validate_resp()`. The generic
  is registered on brms's (`frm_generic_owners`, `S3method(brms::
  autocor, frmtmb_fit)`), and `test-generic-collision.R` lists it.

## 7. The ported brms suite

The gated files (`FRMTMB_BRMS_FIT_TESTS=true`) report these rows as
stale, which is expected: `dev/brmsport-verdicts.tsv` is not edited
here.

| row | recorded | now |
|---|---|---|
| `brm:112` | cannot transfer | holds |
| `brmsfit-methods:112` | cannot transfer | holds |
| `brmsfit-methods:113` | cannot transfer | holds |
| `brmsfit-methods:747` | defect | holds |
| `standata:75` | defect | holds |
| `standata:112` | pass, own words | holds as brms wrote it (brms's message) |
| `standata:970` | cannot transfer | holds |
| `standata:974` | cannot transfer | holds |

Rows whose reason changes and whose outcome does not:

- `brmsfit-methods:955`: the refusal is gone and `update()` refits. The
  row asserts `is(up, "brmsfit")`, so it is now a class-name divergence
  (rule 2 of 2026-09-17), as `:953`.
- `standata:1133`: `drop_unused_levels = FALSE` is accepted. frmtmb
  drops the all-zero `xc` column that brms keeps, so it is a divergence
  and no longer "cannot transfer".

## 8. Test runs and checks

Every test file ran in its own R process with `library(<pkg>)`
attached (`dev/formrobust-run-tests.R`, driven by
`dev/formrobust-run-par.sh`). The seven unchanged extensions loaded
from `rellib-r4` behind the lane's core, as the library lines below
show. The block below is generated by `dev/formrobust-counts.R` and
`dev/formrobust-counts-new.sh` from the per-file logs; the per-file
RESULT lines are in `dev/formrobust-log/suite2.log` and
`dev/formrobust-log/gated1.log`.

The eight gated failures are the stale "now HOLDS" rows of section
7, and nothing else (`dev/formrobust-log/dbg-gated-*.txt`). The gated
tier ran before the last edit, which changed only the reason text of
two `R/compat.R` rows and a vignette paragraph; the ungated suite and
the checks ran on the final code.

### Final ungated suite (NOT_CRAN true) (`dev/formrobust-log/suite2-files`)

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

Libraries the files loaded (package path):

- 199 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb
- 11 x C:/Users/adf44/source/r/rellib-r4/frmtmb.coupling
- 29 x C:/Users/adf44/source/r/rellib-r4/frmtmb.eam
- 10 x C:/Users/adf44/source/r/rellib-r4/frmtmb.latent
- 15 x C:/Users/adf44/source/r/rellib-r4/frmtmb.learn
- 11 x C:/Users/adf44/source/r/rellib-r4/frmtmb.ode
- 44 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb.sample
- 15 x C:/Users/adf44/source/r/rellib-r4/frmtmb.spline

### Gated tier (FRMTMB_BRMS_FIT_TESTS true) (`dev/formrobust-log/gated1-files`)

| package | files | with RESULT | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| frmtmb | 31 | 31 | 2802 | 8 | 0 | 0 | 0 |
| frmtmb.sample | 5 | 5 | 201 | 0 | 0 | 0 | 0 |

Files with a failure, an error, a warning or no RESULT line:

- frmtmb test-brms-suite-brm.R: pass=22 fail=1 err=0 skip=0 warn=0
- frmtmb test-brms-suite-methods.R: pass=159 fail=3 err=0 skip=0 warn=0
- frmtmb test-brms-suite-standata.R: pass=83 fail=4 err=0 skip=0 warn=0

Libraries the files loaded (package path):

- 31 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb
- 5 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb.sample

### The new test files, before (rellib-r4) and after (lane library)

| file | before | after |
|---|---|---|
| frmtmb--test-arma-na-newdata | pass=0 fail=0 err=4 skip=0 warn=0 | pass=13 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-aterm-expr | pass=0 fail=0 err=8 skip=0 warn=0 | pass=49 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-bernoulli-coding | pass=0 fail=0 err=5 skip=0 warn=0 | pass=37 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-brms-api-formrobust | pass=1 fail=0 err=6 skip=0 warn=0 | pass=26 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-offset-grid | pass=1 fail=1 err=3 skip=0 warn=0 | pass=15 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-update-pool | pass=1 fail=1 err=2 skip=0 warn=0 | pass=11 fail=0 err=0 skip=0 warn=0 |
| frmtmb.sample--test-formrobust-draws | pass=0 fail=0 err=2 skip=0 warn=0 | pass=10 fail=0 err=0 skip=0 warn=0 |

`R CMD check --as-cran` (built and checked in
`dev/formrobust-check/<pkg>/`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`):

- frmtmb: Status: 1 NOTE. The NOTE is the V8 math-rendering note on
  the HTML manual, which lane-rules.md expects.
- frmtmb.sample: Status: OK.
- No other extension changed, so none was checked.

After the core check, one roxygen paragraph of `R/interop.R` was
rewrapped to 80 columns (`man/frmtmb-emmeans.Rd` changes only its line
breaks) and the package was reinstalled; nothing else changed.

### Punch round 1

The final runs of the round, on the final code. The blocks below are
generated by `dev/formrobust-counts.R` from the per-file logs
(`suite3.log`, `gated2.log`). The gated tier now includes
`test-offset-grid.R`, whose brms comparison needs Stan. Its eight
failures are the same stale "now HOLDS" rows as round 0 (section 7).

### Punch round 1: final ungated suite (NOT_CRAN true) (`dev/formrobust-log/suite3-files`)

| package | files | with RESULT | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| frmtmb | 199 | 199 | 13682 | 0 | 0 | 164 | 0 |
| frmtmb.coupling | 11 | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 29 | 29 | 1743 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 10 | 10 | 359 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 15 | 429 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 11 | 547 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 44 | 44 | 2224 | 0 | 0 | 4 | 0 |
| frmtmb.spline | 15 | 15 | 553 | 0 | 0 | 1 | 0 |

Libraries the files loaded (package path):

- 199 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb
- 11 x C:/Users/adf44/source/r/rellib-r4/frmtmb.coupling
- 29 x C:/Users/adf44/source/r/rellib-r4/frmtmb.eam
- 10 x C:/Users/adf44/source/r/rellib-r4/frmtmb.latent
- 15 x C:/Users/adf44/source/r/rellib-r4/frmtmb.learn
- 11 x C:/Users/adf44/source/r/rellib-r4/frmtmb.ode
- 44 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb.sample
- 15 x C:/Users/adf44/source/r/rellib-r4/frmtmb.spline

### Punch round 1: gated tier (FRMTMB_BRMS_FIT_TESTS true) (`dev/formrobust-log/gated2-files`)

| package | files | with RESULT | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| frmtmb | 32 | 32 | 2833 | 8 | 0 | 0 | 0 |
| frmtmb.sample | 5 | 5 | 201 | 0 | 0 | 0 | 0 |

Files with a failure, an error, a warning or no RESULT line:

- frmtmb test-brms-suite-brm.R: pass=22 fail=1 err=0 skip=0 warn=0
- frmtmb test-brms-suite-methods.R: pass=159 fail=3 err=0 skip=0 warn=0
- frmtmb test-brms-suite-standata.R: pass=83 fail=4 err=0 skip=0 warn=0

Libraries the files loaded (package path):

- 32 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb
- 5 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb.sample

Mutation runs (`dev/formrobust-log/p1-mut.log`, verbatim):

```
CONTROL--test-arma-na-newdata  RESULT test-arma-na-newdata.R pass=14 fail=0 err=0 skip=0 warn=0 tests=5
CONTROL--test-aterm-expr  RESULT test-aterm-expr.R pass=62 fail=0 err=0 skip=0 warn=0 tests=11
CONTROL--test-bernoulli-coding  RESULT test-bernoulli-coding.R pass=46 fail=0 err=0 skip=0 warn=0 tests=5
CONTROL--test-formrobust-draws  RESULT test-formrobust-draws.R pass=11 fail=0 err=0 skip=0 warn=0 tests=2
CONTROL--test-offset-grid  RESULT test-offset-grid.R pass=31 fail=0 err=0 skip=0 warn=0 tests=7
M05_emm_keep_offsets--test-offset-grid MUTATED frmtmb::emm_drop_offsets  RESULT test-offset-grid.R pass=5 fail=0 err=5 skip=0 warn=0 tests=7
M06_fill_expected_in_predict--test-arma-na-newdata MUTATED frmtmb::arma_cond_fill_dpars  RESULT test-arma-na-newdata.R pass=13 fail=1 err=0 skip=0 warn=0 tests=5
M07_err_omits_ma--test-arma-na-newdata MUTATED frmtmb::autocor_cond_mu_fill  RESULT test-arma-na-newdata.R pass=13 fail=1 err=0 skip=0 warn=0 tests=5
M07b_fill_at_unshifted_mean--test-arma-na-newdata MUTATED frmtmb::autocor_cond_mu_fill  RESULT test-arma-na-newdata.R pass=12 fail=2 err=0 skip=0 warn=0 tests=5
P01_emm_terms_drop_offset--test-offset-grid MUTATED frmtmb::emm_terms  RESULT test-offset-grid.R pass=25 fail=6 err=0 skip=0 warn=0 tests=7
P02_bern_fraction_nowarn--test-bernoulli-coding MUTATED frmtmb::extract_y  RESULT test-bernoulli-coding.R pass=45 fail=1 err=0 skip=0 warn=0 tests=5
P03_refit_no_guard--test-bernoulli-coding MUTATED frmtmb::refit.frmtmb_fit  RESULT test-bernoulli-coding.R pass=45 fail=1 err=0 skip=0 warn=2 tests=5
P04_na_expr_allowed--test-aterm-expr MUTATED frmtmb::assemble_frame Error: The optimizer failed on this model (gaussian, ML, nlminb): NA/NaN gradient evaluation. The likelihood was undefined or unbounded somewhere the optimizer stepped. Refit with verbose = TRUE to see which stage broke, try another optimizer (frmtmb_control(optimizer = "optim")), or start the fit nearer the optimum with the `start` argument of frm()  RESULT test-aterm-expr.R pass=58 fail=0 err=1 skip=0 warn=1 tests=11
P05_grid_drop_all_offsets--test-offset-grid MUTATED frmtmb::emm_drop_offsets  RESULT test-offset-grid.R pass=29 fail=2 err=0 skip=0 warn=0 tests=7
S01_draws_fill_expected--test-formrobust-draws MUTATED frmtmb::arma_cond_fill_dpars  RESULT test-formrobust-draws.R pass=10 fail=1 err=0 skip=0 warn=0 tests=2
P06_emm_keep_all_offsets--test-offset-grid MUTATED frmtmb::emm_drop_offsets RESULT test-offset-grid.R pass=30 fail=1 err=0 skip=0 warn=0 tests=7
```

`R CMD check --as-cran` of the round:

- frmtmb: Status: 1 NOTE. The NOTE is the V8 math-rendering
  note on the HTML manual, which lane-rules.md expects.
- frmtmb.sample: Status: OK.

## 9. Decided not to do, and why

- `conditional_effects()` refuses `trunc()` and `se()` terms whose
  variables are not pinned in `conditions` (`ce_aterms()`, pre-existing
  and deliberate). brms holds such a variable at its mean, and a
  `min(y) - 1` bound at `mean(y) - 1`. Not changed.
- frmtmb.sample's `frm_sample()` takes no `drop_unused_levels`. Its
  `...` goes to tmbstan. Fit with `frm(drop_unused_levels = FALSE)`
  and sample the fit.
- `refit(newresp = )` takes a bernoulli response as 0/1 codes, which
  is what `simulate()` returns. Since punch round 1 any other value is
  refused by name. It is not recoded, because with levels `1` and `2`
  a vector of 1s would be ambiguous between a code and a value.
- The variables of an addition term's expression or of an offset now
  enter `na.action`. So `weights(ifelse(is.na(w), 1, w))` drops the
  rows where `w` is `NA`. brms's `allvars` does the same.

## 10. Defects found and not fixed

- `vignettes/brms-migration.Rmd` still says that `ar()`, `ma()`,
  `arma()`, `cosy()` and `unstr()` "are supported in their covariance
  form only", where `R/autocor.R` fits brms's default `cov = FALSE`
  too. Not this lane's text, so it is reported and not edited.

## 11. Versions

Core: minor (new exported functions `ar()`, `ma()`, `arma()`, `cosy()`,
`unstr()`, `acformula()`, `autocor()`, a new `frm()` argument, and
breaking changes to `update()`, `bernoulli()` and the NA refusal of
addition-term expressions).
frmtmb.sample: patch, with its floor raised to the next frmtmb, for
`arma_cond_fill_dpars()`, `arma_cond_fill_epred()` and
`response_codes_newdata()`.

## 12. Punch round 1 (`dev/reviews/2026-09-30-formrobust.md`)

The review found one blocker and 12 minors. Every item was worked, and
the coordinator decided minors 2 and 3.

**Blocker: `emmeans()` offsets.** brms keeps a predictor's offset
(section 2 corrected). `emm_terms()` keeps the offset attribute, and
`emm_drop_offsets()` drops only the offsets of the predictors that the
selected one does not own, as described in section 2.
`dev/formrobust-p1-offset-check.R` runs the review's five cases on the
lane build (`dev/formrobust-log/p1-offset-check.txt`): the poisson
emmean minus `lin + log(mean(time))` is 0 0; at `time = 1`, the emmean
minus `lin` is 0 0; `epred / exp(lin + lmt)` is 1 1; the sigma emmean
minus `b_sigma + lmt` is 0; and the nonlinear emmean minus
`a + b mean(x)` is 0. `test-offset-grid.R` now asserts these values. It
also compares with brms's `emmeans()` at fixed parameters (Fixed_param
draws equal to frmtmb's estimates, `brms_fixed_fit()`, gated) for
poisson at the grid mean, at `time = 1`, `epred = TRUE`, and the sigma
dpar, all to `expect_exact_num()`. It covers the grid route
(`re_formula = NULL`), which includes the offset once, and a
multivariate grid, where each response gets its own offset.

**m1.** `refit()` refuses a bernoulli value other than 0 and 1 by name.

**m2 (decided: match brms).** Two values are coded 0/1 by level order.
When every distinct value lies strictly inside (0, 1), frmtmb warns
that this looks like a proportion. `-0.5/0.5` and `0/0.5` do not warn.
brms's `standata()` agrees for `-0.5/0.5` and `0.1/0.9` too.

**m3 (decided: fill with draws on draws).** frmtmb.sample's
`posterior_epred()` fills a missing `cov = FALSE` response with a
draw per draw (`arma_cond_fill_epred()`, exported). Core's `fitted()`
keeps the expected-value fill, which is now documented as a divergence
in NEWS, `?fitted.frmtmb_fit`, the compat `fitted` row and the
migration vignette.

**m4: mutants.** `dev/formrobust-p1-run-mut.sh`, through the
reviewer's runner, with the gated tier on. The log is
`dev/formrobust-log/p1-mut.log`, and the controls are green. The
review's M05 (`emm_drop_offsets()` made the identity) now fails 5
blocks, but for a signature reason, so P06 is the same mutant at the
new signature: it fails the nonlinear-offset assertion. M07 (the fill
residual omits the MA term) fails the new hand-recursion test of
`test-arma-na-newdata.R`. S01 (frmtmb.sample fills with the expected
value) fails the new epred-spread assertion of
`test-formrobust-draws.R`. Five more mutants also fail:
- P01, the round-0 offset code: 6 failures, including both gated brms
  comparisons.
- P02, no fraction warning.
- P03, no `refit()` guard.
- P04, no NA refusal.
- P05, the grid route dropping the selected predictor's offset.

The review's M06 and M07b still fail.

**m5.** An addition-term expression that gives `NA` is refused by name,
as brms refuses it. At 0.66.0 `na.omit` dropped those rows, because
the expression was a frame column. NEWS lists this under breaking
changes.

**m6 (REVERSED in round 3 by the user's decision; see section 14. The
text below is the round-1 reasoning, kept as the record.)** brms
refuses `trials(k)`,
`trials(k + 0)` and `weights(wt * k)` for a scalar `k` of the
environment, because it reads variables from `data` and `data2` alone.
frmtmb reads the formula environment everywhere, which is a
model.frame convention. A literal `trials(10)` is accepted by both, so
`k <- 10; trials(k)` meaning the same is the consistent reading. The
bare name now takes the same rule as the expression. This is a
recorded divergence, and the test asserts all three against the
literal.

**m7.**
- A single-value `y2` of `cens()` is refused with the length, as brms
  requires, where it read as "must not be NA on interval-censored
  rows".
- `weights(t * 2)` without a column `t` is refused by name ("R finds
  only the function t()").
- A factor read inside an offset no longer lets model.frame()'s
  "is not a factor" warning escape: `xlev_for()` keeps the whole-term
  variables only.
- A pooled `update()` stores `frmtmb::bf(y ~ x + z, sigma ~ z, family
  = stats::gaussian(link = "identity"))` in the call. The family call
  is written only when evaluating it gives back the same family and
  links; otherwise the family object is stored.
- The 89-column vignette line is rewrapped.
- The frmtmb.sample NEWS names no version.
- The vacuous `expect_message(..., NA)` is replaced by an assertion on
  the "Replacing initial definitions" message.
- The sample test uses `allow_warnings()`.
- `update()` with a complete formula on a multivariate fit is refused
  with brms's message.

## 13. Punch round 2 (the re-check of the review)

**Blocker: a pooled `update()` dropped a family's non-link options.**
`family_call_of()` wrote `huber(k = 3)` as `frmtmb::huber(link = ,
link_sigma = )`, and the update refitted at the default k. The fix is
in two parts:
- The constructor call is written only when every formal of the
  constructor is `link` or `link_<dpar>`.
- The evaluated call must also be `identical()` to the stored family
  on everything but environments.

Otherwise the family object is stored. The review's
`dev/formrobust-rev2-famcall.R` (seed 95) now gives an update logLik of
-237.749592515, the direct fit at k = 3, where it gave -231.49610424,
the default k (`dev/formrobust-log/p2-famcall-after-full.txt`).
`test-update-pool.R` adds an end-to-end `huber(k = 3)` update and
`family_call_of()` checks on `huber(k = 3)`, `whittle(tapers = 4)` and
`cox(df = 6)`. The guard-absent case, `student(link_sigma =
"identity")`, is still written as its call. On the round-1 build this
block failed 3 assertions: the huber logLik and two of the three
families. The round-1 logic, reproduced in
`dev/formrobust-p2-oldfamcall.R`, wrote `huber` and `cox(df = 6)` as
calls, dropping `k` and `df`; `whittle(tapers = 4)` already fell back to
the object (`dev/formrobust-log/p2-oldfamcall.txt`).

**Minor: the multivariate emmeans divergence** is recorded in NEWS and
`?frmtmb-emmeans`. Over both responses with no `resp =`, frmtmb gives
each response its own offset. brms's single `.offset.` column adds one
response's offset to every response.

**Optional: `refit()` with NA.** At 0.66.0 an NA in `newresp` gave
"NA/NaN function evaluation" (`dev/formrobust-p2-refitna.R`,
`dev/formrobust-log/p2-refitna-before.txt`). It is now refused by
name, in `test-bernoulli-coding.R`.

Runs: the whole suite of every package once, one file per process,
package attached (`dev/formrobust-log/suite4.log`). The counts are
generated by `dev/formrobust-counts.R` into
`dev/formrobust-log/counts-p2.md`:

### Punch round 2: whole suite (NOT_CRAN true) (`dev/formrobust-log/suite4-files`)

| package | files | with RESULT | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| frmtmb | 199 | 199 | 13689 | 0 | 0 | 164 | 0 |
| frmtmb.coupling | 11 | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 29 | 29 | 1743 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 10 | 10 | 359 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 15 | 429 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 11 | 547 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 44 | 44 | 2224 | 0 | 0 | 4 | 0 |
| frmtmb.spline | 15 | 15 | 553 | 0 | 0 | 1 | 0 |

Libraries the files loaded (package path):

- 199 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb
- 11 x C:/Users/adf44/source/r/rellib-r4/frmtmb.coupling
- 29 x C:/Users/adf44/source/r/rellib-r4/frmtmb.eam
- 10 x C:/Users/adf44/source/r/rellib-r4/frmtmb.latent
- 15 x C:/Users/adf44/source/r/rellib-r4/frmtmb.learn
- 11 x C:/Users/adf44/source/r/rellib-r4/frmtmb.ode
- 44 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb.sample
- 15 x C:/Users/adf44/source/r/rellib-r4/frmtmb.spline


`R CMD check --as-cran` of frmtmb after this round: Status: 1 NOTE (the
V8 math-rendering note on the HTML manual). frmtmb.sample did not change
this round; its round-1 check was Status: OK.

## 14. Round 3: addition terms read the data alone (user decision)

This round reverses m6 of section 12. By the user's decision, which
matches brms, every variable that an addition term reads must be a
column of `data`. A variable found only in the formula environment is
refused by name, with a message that tells the user to put it in the
data. The reason: predictions on newdata and refits read an addition
term again, and an outside value can be missing or changed by then,
for example after `saveRDS()`/`readRDS()` in another session. On the
response side, as `trials(k)`, that change would be silent. The
divergence of section 12 is now parity with brms.

**Change.**
- `check_aterm_data_vars()` in `R/frame.R` checks every free variable
  of every non-literal addition term, `cens()`'s `y2` included. It runs
  after the existing `se()` refusal, so that refusal still speaks first.
- The frame takes only the term's variables that are in the data.
- Function calls such as `log()` stay allowed, and a literal stays a
  constant.
- `weights(t * 2)` without a column `t` falls under the same rule. The
  message adds "(R finds only the function t() from the formula)".
- On newdata, `aterms_for_newdata()` no longer lets a variable that
  newdata lacks resolve in the formula environment.
- Predictor formulas and offsets are not changed; they read the
  environment.

**Before and after** (`dev/formrobust-p3-envvec.R`, seed 3; logs
`p3-envvec-before.txt` and `p3-envvec-after.txt`):

| term | 0.66.0 | now |
|---|---|---|
| `trials(k)` | refused (model.frame's scalar check) | refused, names `k` |
| `trials(k + 0L)` | refused (model.frame's scalar check) | refused, names `k` |
| `weights(x^2 + k)` | refused (model.frame's scalar check) | refused, names `k` |
| `weights(w_out)`, a full-length vector outside the data | fitted | refused (BREAKING, in NEWS) |
| `trials(n)` with `n` a column | fitted | fitted |

**Tests** (`test-aterm-expr.R`):
- Refused: `weights(wt * k)`, `rate(time * k)`, `trials(k)`,
  `trials(k + 0L)` and `weights(w_out)`.
- Guard-absent cases:
  - the same constants as data columns fit, identical to the literal
    `trials(12)` and to the precomputed weights;
  - `weights(log(wt) + 1)` equals the precomputed column;
  - `I(x * kx)` in a predictor, with `kx` in the environment, equals
    the precomputed column.
- On the round-2 build the four `k` refusals failed:
  `dev/formrobust-log/p3-aterm-expr-before.txt`, 4 of 68 assertions.

**Docs.** The addition-term paragraph of `?bf` and NEWS say this.

**Runs.**
- The whole suite of every package ran once, one file per process,
  package attached (`dev/formrobust-log/suite5.log`). The counts are
  generated by `dev/formrobust-counts.R` into `counts-p3.md`:
  - frmtmb: 199 files, 13695 pass, 0 fail, 1 err, 164 skip.
  - All seven extensions: 0 fail, 0 err, 0 warn.
- The one error was the absent case of `test-brms-parity-defects.R`
  (lane defects' file). It expected check_frame_variables()'s wording
  for a missing `se()` column, and the addition-term rule now names
  that column first. Its expectation was updated to the new message,
  and the file reran at 106 pass, 0 fail.
- The gated brms-suite, agreement and port files (`gated3.log`, 12
  files) show only the same 8 stale "now HOLDS" rows.
- `R CMD check --as-cran` of frmtmb, run after that test edit:
  Status: 1 NOTE (the V8 math-rendering note).
