# Lane formula2: brms formula grammar

Worktree `C:\Users\adf44\source\r\frmtmb-wt-formula2`, branch
`wt-formula2`, base `1f40800d` (frmtmb 0.65.0, frmtmb.sample 0.13.0).
Private library `C:/Users/adf44/source/r/wt-formula2-lib`; the "before"
arm reads `C:/Users/adf44/source/r/rellib-r3`. Nothing is committed.

brms 2.23.0's source was read, not recalled: `dev/formula2-dump-brms.R`
and `dev/formula2-dump-brms2.R` deparse `brmsformula()`, `lf()`,
`validate_formula.*()`, `brmsterms.brmsformula()`,
`expand_dot_formula()`, `validate_terms()`, `get_model_matrix()`,
`terms_fe()`, `terms_re()`, `check_fdpars()`, `update.brmsfit()` and
`update.brmsformula()` into `dev/formula2-brms-src/`. brms's behaviour
on each case is in `dev/formula2-probe-equate.log` (bf() rules, the
generated Stan code) and `dev/formula2-brms-equate-fit.log` (a brms fit:
`variables()`, `summary()`).

## 1. What changed

### 1.1 Equating a dpar to another: `bf(y ~ x, sigma1 = "sigma2")`

brms reads a named character string in `bf()` as an equation and keeps
it in `pfix` beside the constants. Its Stan program declares `sigma2`
alone and sets `sigma1 = sigma2` in `transformed parameters`
(`dev/formula2-probe-equate.log`). frmtmb now does the same:

- `bf()` (R/bf.R, `bf_dots()`) takes a named string and applies brms's
  five rules with brms's messages (`check_dpar_equations()`): not
  itself, not class `mu`, same class (`sub("[[:digit:]]*$", "", dp)`,
  brms's `dpar_class()`), right side not fixed, right side not
  predicted. The rules run on every route: `bf()`, `bf()` on a built
  formula, `+ lf()`, `+ nlf()`.
- `theta1 = "theta2"` is refused by name. brms stops on it with R's
  `invalid 'type' (character) of argument` from the sum it takes to
  normalize fixed mixing proportions (`dev/formula2-probe-equate.log`).
- `lf()` takes named constants and equations, as brms's `lf()` does
  (its entries go through `bf()`); `nlf()` refuses them by name.
- The parser (R/parse.R) refuses what brms refuses once the family is
  known: an unknown target with brms's "Parameter 'sigma3' cannot be
  found.", a target with a formula (the family's `default_forms`
  included). The equated dpar is parsed with the TARGET's link, so its
  value is the target's exactly.
- The frame (R/frame.R, a second pass per response) gives the equated
  dpar the target's linear predictor under its own name: the same
  design, the same `idx` into `betad`. It owns no parameter, so every
  path that computes a dpar from `X %*% betad[idx]` (objective,
  fitted, predict, simulate, frm_linpred, the draws) reads the target's
  value with no code of its own.
- The paths that attach something PER PARAMETER skip it: default
  prior rows and prior resolution (R/priors.R), autoscale planning
  (R/autoscale.R), start values and stationary-escape checks
  (R/fit.R), brms's coefficient table (R/brms-names.R), and
  frmtmb.sample's default priors. A prior on `sigma1` is refused:
  "sigma1 is equated to sigma2 in this model, so it holds no parameter
  of its own".
- Names follow brms. `fixef()`, `get_prior()`, `par_template()` hold
  `sigma2` alone. `variables()` and `summary()` list `sigma1` beside
  `sigma2` with its value, as brms lists its transformed parameter:
  brms's `brms_coef_table()` carries an `equated` attribute and
  `brms_natural_rows()` places the name in predictor order. The draws
  of frmtmb.sample carry a `sigma1` column copied from `sigma2`, added
  after the sampled columns (as `thetaK` is) so no sampled column moves.

### 1.2 `cmc`, cell-mean coding

brms's `cmc` is NOT about category-specific or nonlinear designs, as
the task statement glossed it. Its Rd: "Indicates whether automatic
cell-mean coding should be enabled when removing the intercept by
adding 0 to the right-hand of model formulas." `validate_terms()` sets
the terms' intercept to 1 and `int = FALSE` when a formula has no
intercept and `cmc = FALSE`, and `get_model_matrix()` then drops the
`(Intercept)` column: `y ~ 0 + g` keeps `g`'s treatment contrasts. It
reaches the group-level terms through `terms_re()`. `bf(cmc =)` sets
the location formula only; `lf(cmc =)` sets each formula it holds.

frmtmb: `bf(cmc =)` and `lf(cmc =)` are flags (`check_flag()`); the
frame builds the fixed design with the intercept and drops its column,
keeping the intercept in the stored terms so a prediction builds the
same columns and selects `param_colnames`. A group-level term is given
an intercept for `mkReTrms()` (`bar_with_intercept()`), and its
intercept rows are dropped from the level-major `Zt`; the component
carries `cmc_intercept` so `pred_design()` rebuilds it the same way.

### 1.3 `.` in a formula

brms's `expand_dot_formula()` runs `stats::terms(formula, data = data)`
on the formula and on each parameter formula, keeps the original
attributes, and leaves a formula it cannot expand as written.
`expand_dot_bform()` (R/bf.R) is the same function, called in the four
entry points after `as_bform()`: `frm()`, `get_prior()` (and so
`default_prior()`), `par_template()`, `frm_simulate()`. The frame's
refusal stays for a `data` with no columns (an environment) and its
message now says that; it said "frmtmb does not expand `.`".

`update(fit, bf(~ ., family = acat()))`: brms's `update.brmsfit()`
passes a `bf()` delta to `update.brmsformula()`, which updates the
location formula with `stats::update()` (keeps it when the model is
nonlinear and the delta is only `.`), pools `pforms` and `pfix` with a
later one replacing an earlier one ("Replacing initial definitions of
parameters", a message), and takes the delta's family over the
`family` argument. `update_delta_bform()` (R/confint.R) does the same.
A plain formula delta is unchanged.

### 1.4 A list of families

brms's `validate_formula.mvbrmsformula()` recycles a single family and
takes a list positionally, one per response, refusing a length
mismatch ("If 'family' is a list, it has to be of the same length as
the number of response variables."). Each entry goes through
`validate_formula.brmsformula()`, which sets it only when the `bf()`
has no family. `as_bform()` (R/bf.R) now does the same; a list on a
univariate formula is refused. This agrees with frmtmb's rule for a
single family argument and does not touch the `+ family` rule the user
confirmed on 2026-09-29.

## 2. Validation

### 2.1 brms's log density at frmtmb's estimate

Four new gated rows of `tests/testthat/test-brms-likelihood.R`, each
`brms_lp_check()`: brms's compiled program (flat priors) evaluated at
frmtmb's maximum likelihood estimate, against `logLik()` (`joint =
TRUE` for row 24, which has a group-level term).
`dev/formula2-run-rows.R` runs the four by name with
`frmtmb.brms_lp_report = TRUE`; its lines, pasted from
`dev/formula2-rows-summary.txt`:

```
LPCHECK const 1.136868377e-13 grad 0.000323 ours -610.260451
RESULT row 17c: sigma1 = "sigma2" equates the components' sigma pass=4 fail=0 err=0 skip=0
LPCHECK const 5.684341886e-14 grad 4.66e-15 ours -306.5260617
RESULT row 24: cmc = FALSE on a population- and a group-level term pass=3 fail=0 err=0 skip=0
LPCHECK const -2.842170943e-14 grad 3.7e-05 ours -102.9563961
RESULT row 25: y ~ . expands against the data as in brms pass=2 fail=0 err=0 skip=0
LPCHECK const -2.842170943e-14 grad 0.000118 ours -185.9433169
RESULT row 26: a list of families, one per response pass=2 fail=0 err=0 skip=0
```

`const` is brms's log density minus frmtmb's at the same parameters,
`grad` the largest absolute entry of brms's gradient there; the rows
assert `|const| < 1e-6 * max(1, |logLik|)` and `grad < 1e-3`. Data
seeds: row 17c 11, row 24 5, row 25 6, row 26 8. Row 17c also asserts
that brms's `parameters` block declares `sigma2` and not `sigma1`.
Three new programs went into `dev/stan-cache` (121 files before, 124
after). The whole file, gated, after the change: pass=451 fail=0 err=0
skip=0 (`dev/formula2-g-test-brms-likelihood.log`).

### 2.2 The equated mixture is the hand-written model

`tests/testthat/test-dpar-equate.R` and
`dev/formula2-probe-equate-sample.R` (data seed 11, n = 300): an RTMB
objective written by hand with one log sigma for both components and
theta2 as the reference. Its log:

```
HAND nll at frm estimates 610.505291582506, -logLik(fit) 610.505291582506
HAND rel diff at same point 0.000e+00
HAND optimum -610.505291585714, frm optimum -610.505291582506, rel -5.254e-12
HAND max |est diff| / se: 6.350249e-05
```

At the same parameters the two objectives return the same double: the
difference is exactly 0, computed, not a rounded print. The
hand-written model optimized from its own start reaches the same
optimum to 5.3e-12 relative, and the estimates agree to 6.4e-5
standard errors, which is optimizer stopping. `logLik()` has df 6, one
fewer than the unequated mixture's 7 (`dev/formula2-probe-mix.R
before`). `bf(y ~ x) + lf(sigma1 = "sigma2")` gives the identical
`logLik()`.

### 2.3 Names against brms

brms 2.23.0's fit of the same model (`dev/formula2-brms-equate-fit.R`,
data seed 11, sampler seed 1, log beside it): `variables()` is
`b_mu1_Intercept b_mu2_Intercept b_mu1_x b_mu2_x sigma1 sigma2 theta1
theta2 Intercept_mu1 Intercept_mu2 lprior lp__`, and `summary()` lists
`sigma1` and `sigma2` under Further Distributional Parameters with the
same values. frmtmb's `variables()` is the same set less
`Intercept_mu1`, `Intercept_mu2`, `lprior` and `lp__`, which
`?variables` documents as absent from an ML fit. The order differs as
it already did for every mixture (frmtmb puts `b_mu1_x` before
`b_mu2_Intercept`). `summary()$spec_pars` lists `sigma1` then
`sigma2`, as brms does. frmtmb.sample's draws
(`dev/formula2-probe-equate-sample.R sample`, seed 3) carry `sigma1`,
and its draws are `identical()` to `sigma2`'s. `get_prior()` has a
`sigma2` row and no `sigma1` row, as brms's `default_prior()`
(`dev/formula2-probe-equate.log`).

### 2.4 cmc, dot and family list against brms's Stan data

`tests/testthat/test-formula-cmc.R`, `test-formula-dot.R` and
`test-family-list.R` compare the frame with `brms::standata()`,
`brms:::validate_formula()` and `brms::default_prior()`, which compile
nothing. Covered: brms's own `standata:742` case (`X` is the single
column `gb`); a three-level factor with `lf(sigma ~ 0 + g + (0 + g |
h), cmc = FALSE)` (X_sigma, and Z_1_sigma_1 and Z_1_sigma_2 row for
row, no third coefficient); five dot formulas (`y ~ .`, `y ~ . - x2`,
`y ~ 0 + .`, `y ~ . + (1 | g)`, `y ~ x1 * .`: the expanded formula and
the design columns both equal brms's); `sigma ~ .`; and
`y | weights(w) ~ .`. `y ~ 0 + g` with `cmc = FALSE` gives the
IDENTICAL `logLik()` of the same columns written by hand, and a
prediction on new data rebuilds the group-level columns.

### 2.5 Seen to fail before

Each new file against `rellib-r3` (`dev/formula2-run-test.R frmtmb
<file> before`, logs `dev/formula2-tb-*.log`):

```
RESULT frmtmb test-conditions.R pass=149 fail=0 err=1 skip=0 warn=0
RESULT frmtmb test-dpar-equate.R pass=0 fail=0 err=5 skip=0 warn=0
RESULT frmtmb test-family-list.R pass=0 fail=0 err=3 skip=0 warn=0
RESULT frmtmb test-formula-cmc.R pass=0 fail=0 err=4 skip=0 warn=0
RESULT frmtmb test-formula-dot.R pass=0 fail=0 err=3 skip=0 warn=0
```

Each error is the old behaviour, not a missing symbol: "Cannot
interpret bf() argument 'sigma1'" and "'cmc'", "lf() takes two-sided
formulas", "A formula with `.` is not supported", "The model formula
needs a response (left-hand side)" (the update delta), and "Cannot
interpret `family` of class list". The cmc file gained three
misspelling assertions after this run (section 3).

### 2.6 The ported brms suite

Gated, both arms, one file per process (`dev/formula2-g-*.log` after,
`dev/formula2-gb-*.log` before):

```
after   test-brms-suite-brmsformula.R pass=11 fail=5
        test-brms-suite-standata.R    pass=86 fail=1
        test-brms-suite-priors.R      pass=36 fail=0
        test-brms-suite-methods.R     pass=162 fail=0
before  brmsformula 16/0, standata 87/0, priors 36/0, methods 162/0
```

frmtmb.sample's copy of `test-brms-suite-brmsformula.R` fails the same
five rows (2.7). Every failure is the harness's "now HOLDS but is
recorded as 'cannot transfer'", so these are the verdicts to change;
lane `defects` owns the ledger and I did not touch `dev/brmsport-*`.

| row | now | why |
|---|---|---|
| brmsformula:30, :32, :34, :36, :38 | HOLDS, both tiers | brms's messages, verbatim |
| standata:744 | HOLDS | `bf(cmc = FALSE)`'s X is brms's |
| standata:738, :739, :740 | still no | frmtmb's `cumulative()` has no `disc`, so `lf(disc ~ ...)` is refused as a dpar the family lacks. `lf(cmc =)` itself works: the same case on `sigma` equals brms (2.4) |
| standata:745 | still no | reads `Z_1_1`, which the harness's standata view does not carry; frmtmb's Z equals brms's (2.4) |
| standata:928 | still no, for a new reason | `y ~ .` now expands to `y ~ x1 + x2`, as brms, but brms's data has `x1 == x2` and frmtmb drops the aliased `x2` as lm() does ("rank deficient; dropping column(s): x2", `dev/formula2-probe-928.R`). The verdict should become a divergence on rank deficiency, not a defect on the dot |
| priors:91 to :94 | still no, for a new reason | the family list works, but brms's `y2` (1, 1, 2, 3, ...) is not in (0, 1), and frmtmb's `default_prior()` assembles the frame and refuses the Beta response ("beta: response must lie strictly in (0, 1)"), where brms's does not look at it. With `y2 / 4` frmtmb's table has `sigma` and `ar` for y1, `phi` for y2 and no `ar` for y2: all four assertions (`dev/formula2-probe-priors91.R`) |
| brmsfit-methods:959 | still no, and cannot | `update(fit3, bf(~ ., family = acat()))` now fits (family acat, the old formula; `dev/formula2-probe-959.R`), but the row asserts `is(up, "brmsfit")`, which rule 2 of 2026-09-17 forbids. The verdict should become `divergence`, as row 957's is |

brmsfit-methods:953 (`bf(. ~ ., a + b ~ 1, nl = TRUE)`) is not in this
lane's list. The delta now keeps the nonlinear body and replaces `a`
and `b`, but that row asserts `is(up, "brmsfit")` too.

### 2.7 Whole suites

Both packages, every test file gated (`FRMTMB_BRMS_FIT_TESTS=true`,
`NOT_CRAN=true`), one file per process, 16 at a time
(`dev/formula2-suite.sh`), counted by `dev/formula2-summarize.R` into
`dev/formula2-suite-summary.txt`:

```
SUITE frmtmb files=186 pass=14462 fail=8 err=0 skip=14 warn=0
SUITE frmtmb.sample files=36 pass=2154 fail=5 err=0 skip=1 warn=0
NOT CLEAN:
                                                   log           pkg
        frmtmb-gated-test-brms-suite-brmsformula.R.log        frmtmb
           frmtmb-gated-test-brms-suite-standata.R.log        frmtmb
            frmtmb-gated-test-message-uniqueness.R.log        frmtmb
 frmtmb.sample-gated-test-brms-suite-brmsformula.R.log frmtmb.sample
                          file pass fail err skip warn
 test-brms-suite-brmsformula.R   11    5   0    0    0
    test-brms-suite-standata.R   86    1   0    0    0
     test-message-uniqueness.R    5    2   0    0    0
 test-brms-suite-brmsformula.R   11    5   0    0    0
WITH SKIPS:
           pkg                    file skip
        frmtmb test-drmtmb-agreement.R   13
        frmtmb             test-fuzz.R    1
 frmtmb.sample            test-scale.R    1
```

The brms-suite failures are the flips of 2.6. The two
`test-message-uniqueness.R` failures were two templates this lane had
written twice ("lf() sets ''..." and the multivariate delta refusal).
Both were rewritten; after reinstalling, that file reads pass=6
fail=0, and the 12 test files that call `+ lf(` or `update(` rerun
clean, 814 pass and 0 fail (`dev/formula2-rerun/`). The skips are the
opt-in tiers (drmTMB, fuzz, scale). The grammar fuzz tier was run on
its own with `FRMTMB_FUZZ=true`: pass=2 fail=0
(`dev/formula2-t-fuzz.log`). No warning escaped (warn=0 everywhere).

## 3. The argspell hole

`bf()` gains the formal `cmc` and reads a named string as an equation;
`lf()` gains `cmc` and reads named constants and strings. A misspelled
argument is still refused, asserted in `test-formula-cmc.R`:
`bf(y ~ 0 + g, cmcc = FALSE)` ("Cannot interpret bf() argument
'cmcc'"), `lf(sigma ~ 0 + g, cmcc = FALSE)` ("lf() takes parameter
formulas ..."), and `bf(y ~ x, famly = "gaussian")`. Since punch
round 1 the last one is named as a misspelling first: "bf() has no
argument `famly`. Did you mean `family`? Read as an equation, it fails
brms's rule: Can only equate parameters of the same class ..."
(`nearest_formal()`, the argspell helper, over `bf()`'s own arguments;
`lf()` the same over its own). A string between two names of one class
that the family lacks is refused at the parse, "dpar(s) not available".

## 4. What I did not do, and defects found

- **priors:91 to :94 stay blocked** on `default_prior()` validating the
  response (2.6). Whether `default_prior()` should skip the response
  checks, as brms's does, is for the owner of the prior table.
- **frmtmb.sample puts no default prior on a mixture's `sigma1` or
  `sigma2`**: `default_priors_for()` matches `dpar == "sigma"` only, so
  a mixture samples its sigmas flat where brms puts
  `student_t(3, 0, 2.6)` on them (the sampling log of 2.3 lists only
  the two intercept priors). Pre-existing; not changed.
- **`update(fit, bf(y ~ x2))` with a complete bf() replaces the model**,
  as frmtmb documents, while brms's `update.brmsformula()` pools the
  old dpar formulas into it. Only the delta spelling now reads as
  brms's.
- `summary()`'s Formula line does not print `sigma1 = sigma2`, which
  brms prints under the formula (brms prints every parameter formula
  there; frmtmb prints the location formula alone, before this lane
  too). Cosmetic; not changed.
- An equation between two `hmm()` transition cells: see section 8, M1.
  frmtmb.latent reads `dp$constant` in two places; an equated cell is
  neither constant nor free there, and its suite passes (8.5).
- brms REPLACES an equation with a later formula on the same parameter
  (`bf(y ~ x, sigma1 = "sigma2", sigma1 ~ z)`, with its message
  "Replacing initial definitions of parameters"); frmtmb refuses it as
  a duplicate on all four routes (the reviewer's
  `dev/formula2-rev-equate.log`, E4). Recorded, not changed: the
  refusal is loud, and the `update()` delta already pools with brms's
  replacement rule.
- The task statement glossed `cmc` as "combine main effect columns in
  category-specific or nonlinear designs"; brms's own Rd and code say
  cell-mean coding when the intercept is removed (1.2), and that is
  what was built. brms's separate `gp(..., cmc =)` is not touched.

## 5. Where the edits are

Formula parsing is shared with lanes `aterms2` and `defects`. There is
no `R/formula.R`; `bf()` and `lf()` live in `R/bf.R`, which carries
most of the change. `R/parse.R`, five hunks in `parse_one_response()`
and `print.frmtmb_spec()`: the cmc attribute in `lin_dpar()` and after
the `center` block, the equation checks after `lin_dpar()`, the
equated branch of the dpar loop, the spec print. `R/frame.R`: the dot
message in `check_frame_variables()`; `bar_with_intercept()` before
`grp_structural_ops`; in `assemble_frame()` the equated placeholder
and second pass, the cmc fixed design, the cmc group-level rows, and
the `cmc_intercept` field of a block's components. `R/predict.R`: one
hunk in `pred_design()`. Punch round 1 added: `mm_member_designs()`
and the mm branch of `assemble_frame()`, and the cmc refusal on the
level-read structures (R/frame.R); `re_term_cnms()`, `re_lp_cmc()` and
`re_keep_plan()` (R/re-formula.R); `coef()` (R/methods-fit.R); the
call `frm()` keeps (R/fit.R). Also `R/confint.R` (update, hyp_env_vals),
`R/brms-names.R`, `R/methods-fit.R`, `R/priors.R`, `R/autoscale.R`,
`R/fit.R`, `R/par-template.R`, `R/simulate-new.R`, and frmtmb.sample's
`R/draws-brms.R` and `R/sample.R`.

## 6. Version

A minor bump for frmtmb (new grammar: `cmc`, equations, dot expansion,
family lists, the `update()` delta). A patch bump for frmtmb.sample;
its change reads an attribute only the new frmtmb sets, so with an
older frmtmb nothing changes and no floor is forced. Raising its
frmtmb floor to the release that carries this is what makes `sigma1`
appear in the draws.

## 7. R CMD check --as-cran

Built and checked in `dev/formula2-check/<pkg>/`
(`dev/formula2-check.ps1`), `R_LIBS` the lane library, then rellib-r3,
then the user library, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`:

- frmtmb: `Status: 1 NOTE`, the expected V8 math-rendering note on the
  HTML manual.
- frmtmb.sample: `Status: OK`.

Both were checked again after punch round 1, on the fixed sources, with
the same result: frmtmb `Status: 1 NOTE` (V8), frmtmb.sample
`Status: OK`.

## 8. Punch round 1 (review `dev/reviews/2026-09-29-formula2.md`)

### 8.1 B1: `update(newdata =)` expanded `.` again

`frm()` now expands `.` once and, when that changed the formula,
stores the expanded formula (before any family is attached) as the
call's `formula` (R/fit.R). `update()` evaluates that call, so new
data fits the same model, as brms's `update.brmsfit()` reuses
`object$formula`. A fit without a `.` keeps its call as written. Test
`test-formula-dot.R`, "update(newdata =) keeps the formula `.`
expanded to" (seed 6): a new column is not a predictor, `logLik()` is
identical, `update(data =)` the same, new data without `x2` is refused
by name, and a `bf(y ~ ., sigma ~ x1)` fit keeps its dpar formula. On
the round-0 build it failed 5 of its assertions: fixef
`Intercept x1 x2 extra`, logLik -85.68, and no error on the missing
column (`dev/formula2-p1pre-test-formula-dot.log`). The reviewer's
`dev/formula2-rev-dot.R` rerun on the fixed build
(`dev/formula2-p1-rev-dot.log`): `update(newdata = d2 with an extra
column) fixef: Intercept x1 x2 g2 g3 g4 g5 g6`.

### 8.2 B2: `cmc = FALSE` on a multi-membership term

The mm spec carries `cmc = FALSE` from its predictor, and
`mm_member_designs()` builds each member's design with the intercept
and drops its column, at fit and at prediction (the component keeps
the spec). The other bar shapes: `gr(by =)`, `|ID|` and `||` reach
`mkReTrms()` as plain bars and were already covered (reviewer's log:
`gr(by)` cnms `fb fc fb fc`); `cs()` is not available inside a group
term in frmtmb. The structures that read one coefficient per level,
`ar1()`, `hetar1()`, `cs()`, `homcs()`, `toep()`, `homtoep()`, `ou()`,
`exp()`, `gau()`, `mat()`, are now refused under `cmc = FALSE` by name:
treatment contrasts would leave a level out of an AR or distance
structure, and brms has none of them. Tests (`test-formula-cmc.R`,
seeds 21 and 22): the mm cnms are `fb fc`, the design equals brms's
`W_1_1 * Z_1_k_1 + W_1_2 * Z_1_k_2` for k = 1, 2 with no third
coefficient, the prediction at newdata equals the in-sample one, and
`ar1(0 + t | g)` is refused. Round-0 build: cnms `fa fb fc`, both
design comparisons failed, and `ar1()` fitted without a word.

### 8.3 M1: an equation onto a link-scale dpar

`brms_coef_table()` records an equated dpar only when its target's row
is natural. A link-scale target (an `hmm()` transition cell) is a
coefficient, and the cell equated to it is then no parameter at all,
so it is not listed. The reviewer's `dev/formula2-rev-hmm.R` (seed
4001) on the fixed build (`dev/formula2-p1-rev-hmm.log`): logLik
-843.8734 df 5 against the mixture's -843.8734 df 5 (relative 8.4e-13);
`variables()` `b_mu1_Intercept b_mu2_Intercept b_tr22_Intercept
sigma1 sigma2`; `summary()` returns. `hypothesis(h1, "tr22_Intercept =
0")` works, and `"tr12_Intercept = 0"` is refused as a parameter that
cannot be found (`dev/formula2-p1-hmm-hyp.R` and log). The core test
("an equation onto a link-scale dpar lists no natural value") sets the
family's `link_scale_dpars` by hand on an equated mixture, because core
does not suggest frmtmb.latent; on the round-0 build it stopped with
"attempt to apply non-function".

### 8.4 M2 and minors

- M2: `re_term_cnms()` takes `cmc`, and `re_keep_plan()` reads each
  predictor's setting (`re_lp_cmc()`), so `re_formula = ~ (0 + f | h)`
  names the `fb fc` term of a `cmc = FALSE` fit. Test in
  `test-formula-cmc.R`; round-0 build: "matches no group-level term".
- `coef()` has no `$sigma1` block (brms's `coef()` covers group-level
  terms only and has no such block either). Test in
  `test-dpar-equate.R`; round-0 build listed it.
- `print()` of a `bf()` or `lf()` shows `(cmc = FALSE)` and
  `(center = FALSE)`; it already showed `sigma1 = sigma2`.
- The misspelling message, section 3.

### 8.5 Counts

New and changed test files, after the fix (`dev/formula2-p1-test-*.log`):

```
RESULT frmtmb test-conditions.R pass=150 fail=0 err=0 skip=0 warn=0
RESULT frmtmb test-dpar-equate.R pass=36 fail=0 err=0 skip=0 warn=0
RESULT frmtmb test-family-list.R pass=7 fail=0 err=0 skip=0 warn=0
RESULT frmtmb test-formula-cmc.R pass=31 fail=0 err=0 skip=0 warn=0
RESULT frmtmb test-formula-dot.R pass=25 fail=0 err=0 skip=0 warn=0
RESULT frmtmb test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
```

The same three files on the round-0 build, with the new tests
(`dev/formula2-p1pre-*.log`):

```
RESULT frmtmb test-dpar-equate.R pass=32 fail=1 err=1 skip=0 warn=0
RESULT frmtmb test-formula-cmc.R pass=22 fail=6 err=2 skip=0 warn=0
RESULT frmtmb test-formula-dot.R pass=20 fail=5 err=0 skip=0 warn=0
```

Both whole suites again, gated, one file per process
(`dev/formula2-suite/`, round 0's logs moved to
`dev/formula2-suite-r0/`), `dev/formula2-suite-summary-p1.txt`:

```
SUITE frmtmb files=186 pass=14486 fail=6 err=0 skip=14 warn=0
SUITE frmtmb.sample files=36 pass=2154 fail=5 err=0 skip=1 warn=0
```

The 11 failures are the ported-row flips of 2.6 and nothing else
(core brmsformula 5, standata 1; frmtmb.sample brmsformula 5). The
skips are the opt-in tiers. frmtmb.latent's suite with this core and
rellib-r3's frmtmb.latent: 10 files, pass 359, fail 0, err 0, skip 2
(scale tier), warn 0 (`dev/formula2-suite/frmtmb.latent-*.log`).

## 9. Punch round 2

### 9.1 B3: the cmc refusal refused correct models

The refusal on the ten level-read structures fired on every bar
without an intercept. `cmc_changes_columns()` (R/frame.R) now builds
the bar's left-hand side on the model frame both ways; when the
columns are the same (numeric slopes, `f:x` without `f`), the bar is
neither rewritten nor refused. The refusal, and its "would leave one
level out", now fire only where a factor's first level really goes.
If the columns cannot be built there, the old path (rewrite, refuse)
is kept.

Test `test-formula-cmc.R`, "cmc = FALSE is accepted where it changes
no column" (seed 41): `cs`, `ar1`, `hetar1`, `homcs`, `toep` and
`homtoep` over `0 + x1 + x2`, and `cs(0 + f:x1 | g)`, give the same
cnms under `cmc = FALSE` as under the default. On the round-1 build it
stopped at the first, "cmc = FALSE does not apply to cs(0 + x1 + x2 |
g)" (`dev/formula2-p2pre-test-formula-cmc.log`). The refused case,
`ar1(0 + t | g)` over a factor, stays refused.

### 9.2 Minors

- `fit$call` of a plain `y ~ .` formula now carries the expanded
  formula itself, so it prints `frm(formula = y ~ x1 + x2, ...)`
  (asserted in `test-formula-dot.R`; round-1 build: the deparsed
  list). A `bf()` with a `.` still stores the expanded `bf()` object.
- Recorded, not changed: a `y ~ . - w` fit refit on data without `w`
  is refused ("The model uses `w`"), because the expanded formula is
  `y ~ (x1 + w) - w`. brms stores the same formula and its `standata()`
  without `w` works. The refusal is loud, and frmtmb refuses
  `y ~ x1 + w - w` without `w` the same way when there is no dot.

### 9.3 Counts

`test-formula-cmc.R` pass=44 and `test-formula-dot.R` pass=26, fail,
err, skip and warn 0. The 27 core files that reach the guard or the
changed call (covariance structures, `mm()`, `cmc`, `update()`, `y ~
.`), one per process (`dev/formula2-rerun2/`): pass 2810, fail 0, err
0, skip 0, warn 0.

Every extension's ungated suite against this build (core and
frmtmb.sample from the lane library, the other extensions from
rellib-r3), one file per process (`dev/formula2-suite/<pkg>-plain-*`):

```
frmtmb.coupling logs 11 results 11 pass 542 fail 0 err 0 skip 5 warn 0
frmtmb.eam logs 29 results 29 pass 1733 fail 0 err 0 skip 3 warn 0
frmtmb.latent logs 10 results 10 pass 359 fail 0 err 0 skip 2 warn 0
frmtmb.learn logs 15 results 15 pass 429 fail 0 err 0 skip 13 warn 0
frmtmb.ode logs 11 results 11 pass 547 fail 0 err 0 skip 1 warn 0
frmtmb.sample logs 36 results 36 pass 1946 fail 0 err 0 skip 4 warn 0
frmtmb.spline logs 15 results 15 pass 553 fail 0 err 0 skip 1 warn 0
```

The six extensions other than frmtmb.sample total 4163 passes and 25
skips, and frmtmb.sample 1946 and 4, the reviewer's counts for both
arms.

`R CMD check --as-cran` on frmtmb after punch round 2
(`dev/formula2-check/frmtmb/`): `Status: 1 NOTE`, the V8
math-rendering note.
