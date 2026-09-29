# Lane wt-csfactor: `cs()` on a discrete predictor, and `y ~ x + cs(x)`

Two defects, both pre-existing through 0.64.0. A was a silent wrong
answer; B was an unidentified model that was fitted and reported.

Base commit: 0.64.0. Reference build for every "before" number:
`C:/Users/adf44/source/r/rellib-r3`. Lane library:
`C:/Users/adf44/source/r/wt-csfactor-lib`. brms 2.23.0, rstan 2.32.7,
StanHeaders 2.39.1 from the user library.

## What was wrong

**A.** `assemble_frame()` stored one `cs()` term as
`as.numeric(eval(cexpr, mf, env))`, a single numeric column with one
`bcs<j>` vector of `K - 1` coefficients. On a factor that is the INTEGER
CODES, so the model was linear in 1, 2, 3 instead of carrying one
coefficient per level, and `ord_cs_values()` re-evaluated the same
expression against `newdata` alone, so a `newdata` factor was recoded
from ITS OWN levels. On a character column `as.numeric()` produced `NA`s
and the fit died.

**B.** Nothing checked whether a `cs()` column was already spanned by the
population-level design. `y ~ x + cs(x)` adds `-x * b` to every category
and `-x * bcs[k]` to category `k`, so `b -> b + t`, `bcs[k] -> bcs[k] - t`
leaves the density unchanged: one flat direction, reported as an estimate
with a standard error.

## What changed

| file | change |
|---|---|
| `R/frame.R` | new `cs_term_design()`, `cs_term_newdata()`, `check_cs_identified()`, and from punch round 1 `cs_alias_cols()`, `cs_alias_phrase()`, `cs_alias_advice()`, `mo_indicator_basis()`; the `cs_info` loop builds one entry per design COLUMN and stores the term's model matrix as `lp[["cs_mm"]]` |
| `R/influence.R` | punch round 1: a deletion refit that succeeds with FEWER coefficients warns, naming the unit and the `cs()` level, instead of leaving a silent `NA` in `cooks.distance()`; new `influence_coef_labels()` |
| `R/predict.R` | `ord_cs_values()` takes the newdata columns from the new `cs_newdata_columns()`, which rebuilds each term's model matrix through the FIT's terms, `xlevels` and `contrasts` |
| `R/compat.R` | the `cs_pred()` x `group:ordinal_cs` registry note says what the columns are and that the both-sides form is refused |
| `R/fit.R` | new `?frm` section "Category-specific effects, cs()" |
| `R/confint.R` | `?variables` documents `bcs_<column>[k]` |
| `vignettes/brms-migration.Rmd` | the departure from brms for `y ~ x + cs(x)`, in the "What changes" list |
| `NEWS.md` | a `# frmtmb (development version)` heading, three bullets |
| `tests/testthat/test-cs-factor.R` | NEW, 84 assertions in 15 blocks |
| `tests/testthat/test-brms-likelihood.R` | six new rows in the row-12 block |
| `tests/testthat/test-v17.R`, `test-review-v29.R` | both fitted `y ~ x + cs(x)`; each now asserts the refusal and fits `cs(x)` alone, which spans the same set of distributions |
| `dev/test-backlog.md`, `dev/sratio-findings.md` | the two entries marked done, with the corrected number for B |

Nothing in `extensions/` needed a change; see "frmtmb.sample" below.

Every consumer reads `lp[["cs"]]` as a flat list of entries and only
needs `par` and `vals`, so one entry per COLUMN instead of per TERM
carried the fix through `objective.R`, `brms-shapes.R`, `confint.R`,
`priors.R`, `importance.R`, `predict.R` and `simulate-newdata.R`
unchanged. `lp[["cs_mm"]]` is the per-term model matrix spec (`terms`
with `predvars`, `xlevels`, `contrasts`, `colnames`, `label`); each
column entry carries `tid` and `col` into it. A frame with no `cs_mm`
falls back to the old expression path, so an object assembled by an
older build still predicts.

## A: before and after, seed 405, n = 500

Script `dev/csfactor-repro.R` (the wt-predfix reviewer's construction
from `dev/predfix-review2/r2-cs2.R`). Logs `dev/csfactor-log/before.txt`
and `dev/csfactor-log/after.txt`. Model `frm(bf(yo ~ cs(fc)),
family = sratio())`, `fc` a three-level factor.

| quantity | before (0.64.0) | after | the model written with hand dummies |
|---|---|---|---|
| `logLik` | -446.9233 | -446.6361 | -446.6361 |
| coefficients | 2 | 4 | 4 |
| `fixef()` rows | `fc[1]`, `fc[2]` | `fcb[1]`, `fcb[2]`, `fcc[1]`, `fcc[2]` | `fb[1]`, `fb[2]`, `fcc[1]`, `fcc[2]` |
| `variables()` | `bcs_fc[k]` | `bcs_fcb[k]`, `bcs_fcc[k]` | |
| `default_prior()` class `b` coef rows | `fc` | `fcb`, `fcc` | |
| `fitted()` at `factor("c")` | 0.1876, 0.6725, 0.1399 | 0.7518, 0.0922, 0.1560 | 0.7518, 0.0922, 0.1560 |
| `fitted()` at level a | 0.1876, 0.6725, 0.1399 | 0.1944, 0.6722, 0.1333 | 0.1944, 0.6722, 0.1333 |
| `cs()` on a character column | optimizer died, "NA/NaN gradient evaluation", after a "NAs introduced by coercion" warning | -446.6361, 4 coefficients | |
| `cs()` on `factor(1:3)` | -446.9233, 2 coefficients | -446.6361, 4 coefficients | |

The 0.64.0 answer at `factor("c")` is level a's row, which is the
recoding half of A: a one-row `newdata` factor has one level, `as.numeric`
made it 1, and level 1 is the reference. The corrected row equals the
empirical shares of level c in the data (0.7518, 0.0922, 0.1560), which
is what a saturated three-level model must reproduce; so does every other
row. The in-sample `fitted()` matrix equals the hand-dummy fit's
BITWISE (max absolute difference 0, `dev/csfactor-log/consumers.txt`).

`cs(factor(1:3))` is the case where the defect was invisible in the
numbers: the codes and the levels coincide, so only the coefficient COUNT
and the names gave it away.

## brms's parameterization, read and not guessed

`dev/csfactor-brms.R`, log `dev/csfactor-log/brms.txt`. `standata()` and
`stancode()` only; nothing sampled.

| formula | `X` columns | `Xcs` columns | `get_prior()` class `b` coefs |
|---|---|---|---|
| `yo ~ cs(fc)` | none | `fcb`, `fcc` | `fcb`, `fcc` |
| `yo ~ cs(fch)` (character) | none | `fchb`, `fchc` | `fchb`, `fchc` |
| `yo ~ cs(fnum)` (`factor(1:3)`) | none | `fnum2`, `fnum3` | `fnum2`, `fnum3` |
| `yo ~ x + cs(x)` | `x` | `x` | `x` |
| `yo ~ fc + cs(fc)` | `fcb`, `fcc` | `fcb`, `fcc` | `fcb`, `fcc` |
| `yo ~ cs(x)` | none | `x` | `x` |

So brms builds treatment dummies for a factor AND for a character
column, names them exactly as `model.matrix()` does, and declares
`matrix[Kcs, nthres] bcs`: one coefficient row per dummy per threshold.
frmtmb now produces the same columns in the same order and names them
`bcs_<column>[k]`, brms's own spelling.

For `y ~ x + cs(x)` brms declares BOTH `vector[Kc] b` and
`matrix[Kcs, nthres] bcs` over the same column, so it fits the
unidentified model.

The `bcs_` NAME was read out of brms's namespace rather than inferred
(`dev/csfactor-brms-names.R`, log `dev/csfactor-log/brms-names.txt`).
`brms:::rename_cs()` builds `csenames <- paste0("bcs", p, "_",
csenames)` from `bframe$frame$cs$vars`, which is the `Xcs` column set,
and `brms:::prepare_predictions_cs()` reads them back with the regex
`paste0("^bcs", p, "_", escape_all(csef), "\\[")`. So the name is
`bcs_<Xcs column>[k]`, and for a factor the Xcs columns are the dummies.
`brm(empty = TRUE)` was tried first and is no use here: an empty
`brmsfit` has no draws, so `variables()` returns nothing.

## Log-density identity against brms's compiled program

`dev/csfactor-lp.R`, seed 405, n = 400, log `dev/csfactor-log/lp.txt`.
It reuses `tests/testthat/helper-brms.R`'s `brms_lp_check()`:
`rstan::log_prob(sf, upars, adjust_transform = FALSE)` at frmtmb's
estimates against `logLik(fit)`, plus the gradient there.

| row | measured constant | max abs gradient | `logLik` |
|---|---|---|---|
| `sratio`, `y ~ x + cs(fc)` | 0.000e+00 | 7.132e-04 | -336.0989836 |
| `cratio`, `y ~ x + cs(fc)` | 0.000e+00 | 7.132e-04 | -336.0989836 |
| `acat`, `y ~ x + cs(fc)` | 0.000e+00 | 1.653e-05 | -335.6000751 |
| `sratio`, `y ~ x + cs(fch)` | 0.000e+00 | 7.132e-04 | -336.0989836 |
| `sratio`, `y ~ x + cs(fnum)` | 0.000e+00 | 7.132e-04 | -336.0989836 |
| `sratio`, `y ~ cs(fc) + cs(x)` | 0.000e+00 | 7.994e-05 | -331.9562306 |

The constant is zero to the printed precision AND the comparison inside
`brms_lp_check()` is `abs(lp - ours) < 1e-6 * max(1, abs(ours))`, so this
is an identity check and not a tolerance. `sratio` and `cratio` land on
the same log-likelihood here because with three categories and a
category-specific term on every boundary the two families span the same
set of distributions; both agree with brms, which is the claim.

These six rows are now in `test-brms-likelihood.R`'s row-12 block. They
also pin the column ORDER: brms builds its cs design from `~ f + z`, so
the factor's dummies come before `z`, and a translator that reversed them
would fail on the two-term row.

The gated file was then run whole, with `FRMTMB_BRMS_FIT_TESTS=true` and
`NOT_CRAN=true`, compiling 64 Stan programs from scratch because this
machine has no `dev/stan-cache` (`dev/csfactor-log/gated-likelihood.txt`,
cache written to `dev/csfactor-stan-cache`, gitignored):

    RESULT test-brms-likelihood.R pass=440 fail=0 error=0 skip=0 warn=0

Nothing skipped, which is what the gated tier exists for.

## B: what was measured, and the decision

`yo ~ x + cs(x)`, seed 405, n = 500, `dev/csfactor-log/before.txt`:

    Intercept[1] -0.248400  se 9.0249e-02
    Intercept[2]  0.804003  se 1.3090e-01
    x             0.110717  se 2.709529e+05
    x[1]         -0.129316  se 2.709529e+05
    x[2]          0.240033  se 2.709529e+05
    logLik -513.249191

`dev/sratio-findings.md` recorded NaN standard errors for this model on
the 0.64.0 build. On THIS construction they are finite and 2.7e5, so the
NaN belongs to that lane's data and not to the model; both are the same
flat direction reaching `vcov()` differently. The record has been
corrected in place rather than contradicted elsewhere.

The project tiebreaker is to match brms unless it is clearly obvious brms
should do it the other way. brms fits the model. The refusal is the
departure the task asked for, and the reason it is the right one here is
that frmtmb's `frm()` is maximum likelihood: there is no prior to hold
the ridge, the optimizer stops at an arbitrary point on it, and the
standard error frmtmb prints is a number with no meaning. `cs(x)` alone
spans exactly the same set of distributions, so nothing is lost.

The check is a RANK test, not a name test, so `y ~ I(x) + cs(x)` and
`y ~ x1 + x2 + cs(x1 + x2)` are caught too:

    rank(cbind(1, X, Z[, 1:j])) < rank(cbind(1, X)) + j

for each `cs()` column `j` in order, with every column scaled to unit
norm before `qr()`. It compares against the RANK of `cbind(1, X)` and not
its column count, because the stored `X` carries all-zero placeholder
columns for `mo()`, `mi()` and `me()` terms that the objective fills
later; comparing counts would have fired on every such model. The
constant column `1` is in the test because an ordinal family has no
intercept in `X` and the thresholds are one intercept per boundary, so a
constant `cs()` column is unidentified too.

False alarms, measured on the shapes the suite and the vignettes use
(`tests/testthat/test-cs-factor.R`, "the identifiability check does not
fire on a separable model", plus the 28 test files below):
`y ~ z + cs(x)`, `y ~ fc + cs(x)`, `y ~ x + cs(fc)`, `y ~ cs(fc) + cs(x)`,
`y ~ cs(x)`, `y ~ cs(x) + cs(x:z)`, `rating ~ period + carry + cs(treat)`,
`rating ~ x1 + cs(x2) + (1 + x2 || subject)` and
`y ~ x + cs(z)` all fit. Nothing in the 28 files newly refused. The
random-slope case is worth naming: `(1 + x2 || subject)` beside `cs(x2)`
is NOT refused, and should not be, because a random slope is centered at
zero with a variance rather than carrying a free mean.

Two files DID fit `y ~ x + cs(x)` and had to change: `test-v17.R`'s
"cs() category-specific effects match direct ML" and
`test-review-v29.R`'s "cs() terms enter the ordinal predictions". Both
now assert the refusal and fit `cs(x)` alone. In `test-v17.R` the
generating model itself writes `x` twice (`-0.3 * x` plus
category-specific `g_t`), so only the sums `0.3 + g_t[k]` ever entered
the density; the hand-written reference `optim()` was rewritten with six
parameters instead of seven and the identity still holds, which is direct
evidence that the two models span the same distributions.

## Every consumer, followed

`dev/csfactor-consumers.R`, logs `dev/csfactor-log/consumers.txt` (after)
and `consumers-before.txt` (before).

| consumer | result |
|---|---|
| fitting (`objective.R`) | one `bcs<j>` per column, log density equals brms's |
| `fitted()`, `frm_linpred(type = "response")` | correct per level; in sample identical to the hand-dummy fit (max abs difference 0) |
| `predict()` (20000 draws, `propagate_error = FALSE`) | 0.74555, 0.09660, 0.15785 at `factor("c")` against the fitted 0.7518, 0.0922, 0.1560, inside the draw count's own Monte Carlo standard error |
| `simulate()` | same, at 20000 draws |
| `conditional_effects(effects = "fc")` | one row per level per category, equal to the hand-dummy fit's predictions to 1e-10 relative |
| `summary()`, `fixef()`, `variables()`, `hypothesis()` | `fcb`/`fcc` rows, `bcs_fcb[k]` names |
| `default_prior()` | one class `b` row per dummy, matching brms's `get_prior()` |
| `set_prior(class = "b", coef = "fcc")` | reaches that dummy's pair only: (-2.5295, 2.1438) to (-0.1050, 0.0322) while `fcb`'s pair does not move. Before, the same call was "Prior target not found (class=b, coef=fcc)" |
| `set_prior(class = "b")` | covers all four |
| a `newdata` level the fit never saw | refused by name. Before, `factor("zz")` was silently recoded and returned level a's row |
| `emmeans()` | UNCHANGED and still refuses; see below |
| `residuals()` | `type = "response"` is refused for `sratio`, as before and by design |
| `frmtmb.sample::frm_sample()` | no change needed; see below |

### frmtmb.sample

`dev/csfactor-sample.R` and `dev/csfactor-sample-before.R`, n = 250,
2 chains, 600 iterations, seed 7. Logs `dev/csfactor-log/sample.txt`
and `sample-before.txt`.

The extension reads `cs()` only through core's exported
`cs_offsets_add()` and `brms_coef_table()`, and its fixef regex
(`^b((|s|cs|sp|mo|me|mi|m))_`) is name-agnostic, so NOTHING in
`extensions/frmtmb.sample/` was changed. Verified by installing it into
the lane library and sampling:

| quantity | before | after |
|---|---|---|
| `fixef(ds)` rows | `Intercept[1]`, `Intercept[2]`, `fc[1]`, `fc[2]` | ... `fcb[1]`, `fcb[2]`, `fcc[1]`, `fcc[2]` |
| `posterior_epred()` at `factor("c")` | 0.1583, 0.6567, 0.1850 (level a's) | 0.8386, 0.0494, 0.1120 |
| the same row of the three-level grid | 0.7991, 0.0803, 0.1207 | 0.8386, 0.0494, 0.1120 |
| ML fit at the same rows | 0.7960, 0.0808, 0.1232 | 0.8395, 0.0494, 0.1111 |

So the wrong answer reached the posterior surface too, and the fix
reaches it without an extension change.

## What I decided NOT to do

- **No `droplevels()` call of my own on a `cs()` factor.** WITHDRAWN AS
  FIRST WRITTEN, and the correction is on this page rather than only in
  conversation. The bullet said brms leaves an unused level's all-zero
  dummy in place and that frmtmb's constant-column refusal then catches
  it. Both halves are false, and the punch-round reviewer measured both.
  `assemble_frame()` builds its model frame with
  `drop.unused.levels = TRUE` (`R/frame.R`, two call sites), UPSTREAM of
  `cs_term_design()`, so frmtmb never sees a zero column: on a
  three-level factor with `table(f) = 114 126 0` the fit succeeds with
  ONE cs column, `fixef()` rows `fb[1]`, `fb[2]`, and no refusal
  (`dev/csfactor-rev-unused.R`, seed 405, n = 240; reproduced here as
  logLik -257.5667, `dev/csfactor-log/p1.txt`). brms drops it too:
  `standata()` gives `Kcs = 1`, `Xcs = fub`
  (`dev/csfactor-rev-brms.R`). The two packages AGREE, for a reason that
  has nothing to do with the new check, and `test-cs-factor.R`'s "an
  unused cs() factor level drops its column, as in brms" now pins it,
  `nlevels(f) == 3` asserted first so the level really is declared.
- **No escape hatch for `y ~ x + cs(x)`.** The refusal is at frame
  assembly, so it also stops `frmtmb.sample::frm_sample()`, where a
  proper prior WOULD identify the model as it does in brms. That is a
  real narrowing and it is recorded here rather than worked around: an
  argument to disable an identifiability check is a thing that gets
  passed by habit. If the Bayesian case is wanted later, the right shape
  is a `frm_sample()`-only relaxation gated on a proper prior covering
  the block, not a `frm()` argument.
- **No change to `emmeans()`.** It refuses `~ fc` on a `cs(fc)`-only
  model with "No variable named fc in the reference grid", on the
  reference build and on this one, because a `cs()` term is not a column
  of `X` and the reference grid is built from `X`. Pre-existing, out of
  scope, filed below.
- **No change to the `bcs<j>` component numbering.** A factor now spends
  several slots (`bcs2`, `bcs3` for a three-level factor whose family
  already took slot 1 for its thresholds). The numbering was already
  documented as "read the terms in order rather than parsing the number"
  for `mo()`'s `zeta<j>`, and the brms-facing names are the `bcs_<column>`
  ones.
- **No `hypothesis()` change.** `hypothesis(fit, "bcs_fcc[1] = 0")` needs
  `class = NULL`, as every non-`b_` name does; that convention is
  pre-existing and documented.

## Defects found and NOT fixed

- **`variables()` on a `frm_sample()` draws object of an ORDINAL fit
  lists the RAW internal names.** Filed first as a `cs()` defect; the
  punch-round reviewer showed it is WIDER than that, so the filing is
  corrected here and in `dev/test-backlog.md`. It gives `bcs2_1`,
  `bcs2_2` where `variables(fit)` gives `bcs_fcb[1]`, `bcs_fcb[2]`, AND
  `tau_raw_1`, `tau_raw_2` where `variables(fit)` gives
  `b_Intercept[1]`, `b_Intercept[2]` - and it does the `tau_raw` half for
  an ordinal fit with NO `cs()` term at all (`sratio, yo ~ x` gives
  `b_x tau_raw_1 tau_raw_2 lp__`, `dev/csfactor-rev-docs2.R`). A gaussian
  fit's draws are correct. Measured on both builds
  (`dev/csfactor-log/sample.txt`, `sample-before.txt`), so pre-existing
  and unchanged. `fixef(ds)` DOES carry brms's rows, through
  `draws_fixef_ordinal()`, so the gap is in `variables()` alone.
  `?variables` said the draw columns "follow the same convention"; that
  sentence now records the exception instead. Not fixed: it is a
  frmtmb.sample naming question, not a wrong answer.
- **`fitted(mv, newdata = )` errors on any multivariate ordinal fit**,
  with R's own `vapply()` message "values must be length 1, but
  FUN(X[[1]]) result is length 0". Found by the punch-round reviewer with
  `cs()` and without it, same message either way, and in-sample
  `fitted(mv)` works. Not about `cs()`, pre-existing, filed in
  `dev/test-backlog.md` with the construction.
- **`mi(x) + cs(x)` is NOT settled.** The rank test cannot see an `mi()`
  or `me()` column, both being zero placeholders at assembly, so both
  pairs are accepted. `me(x, sdx) + cs(x)` looks identified: 0 of 5
  standard errors NaN on both builds (`dev/csfactor-log/p1.txt` and the
  reviewer's own run). `mi(x) + cs(x)` was not measured, because the
  observed half of an `mi()` column IS the data column, so a partial
  aliasing is plausible and the question needs its own construction with
  a controlled missingness fraction. `?frm` and `NEWS.md` say the gap
  exists rather than implying the check is exhaustive.
- **`emmeans()` cannot reach a `cs()` predictor.** Above. A reference
  grid over a variable that appears only inside `cs()` is not built, so
  `emmeans(fit, ~ fc)` errors with emmeans's own message rather than a
  designed refusal. Pre-existing on 0.64.0.
- **A fresh Stan compile on this machine needs
  `R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win` set by the
  caller**, because `path.expand("~")` here is `C:\Users\adf44` and R
  reads only `<HOME>/.R/Makevars.win`. Without it `rstan::stan_model()`
  fails with "invalid connection", which names neither the compiler nor
  the flag. `dev/release/run-tests.R` already sets it; any lane script
  that compiles Stan must too. Recorded because the failure message is
  misleading, not as a package defect.

## Tests

`dev/csfactor-run-tests.R` runs one file per R process and sums
`passed`, `failed`, `error` and `skipped` from the result frame.
`dev/csfactor-batch.ps1` runs up to seven of them at a time.

`tests/testthat/test-cs-factor.R`, SEEN TO FAIL on the reference build
before it passed here (`dev/csfactor-log/test-before.txt`,
`test-after.txt`):

| build | pass | fail | error | skip |
|---|---|---|---|---|
| rellib-r3 (0.64.0) | 7 | 27 | 3 | 0 |
| lane | 46 | 0 | 0 | 0 |

The 27 failures cover both halves of A (the log-likelihood gap, the
coefficient count, the names, the `newdata` recoding, the character
column, the prior coverage, the stored frame) and all seven refusals of
B. Two of its tests were rewritten after the first "before" run showed
them PASSING on the unfixed build: `predict()`/`simulate()` and
`conditional_effects()` compared the factor fit against its own
`frm_linpred()`, and both halves of A moved together, so the comparison
agreed while both sides were wrong. They now compare against the
hand-dummy fit and fail on the reference build.

No absolute numeric tolerance appears in the file. Each one is either a
ratio to the log-likelihood's own magnitude, a ratio to the largest
probability in the matrix being compared, the Monte Carlo standard error
`sqrt(p(1-p)/nsim)` of the draw count actually used, a ratio to the
TREATMENT-coded fit's own distance from the same target (the
ordered-factor test, where how close a saturated fit comes to the
empirical shares is the optimizer's business and not the contrast
coding's), or `expect_identical()` where the quantity is exact (a 0/1
dummy column). The ordered-factor bound was written as
`1e-6 * max(E)` first and FAILED at 3.0e-6 against 8e-7 on the fixed
build, which is what sent it to the measured reference.

### The whole core suite

`dev/csfactor-log/suite-pre.txt`, one file per R process, seven at a
time, `NOT_CRAN=true`, against the lane library:

    RESULT lines: 180 of 180
    pass=11923 fail=0 error=0 skip=152

Thirteen files report `pass=0`. Every one is a TOP-LEVEL gate, not a
file that asserted nothing: the eleven `test-brms-suite-*.R` files and
`test-brms-priors.R` call `skip_unless_brms_suite()` or
`skip_unless_brms_fit()` above the first `test_that()`, which needs
`FRMTMB_BRMS_FIT_TESTS=true`, and `test-drmtmb-agreement.R` (13 skips)
needs drmTMB, which is not in the lane library. All thirteen are
pre-existing for an ungated lane run.

### Ordered factors

`dev/csfactor-ordered.R`, log `dev/csfactor-log/ordered.txt`.
`model.matrix()` gives an ordered factor `contr.poly`, and so does brms:
both build `fo.L` and `fo.Q`, and frmtmb's stored columns are brms's
`Xcs` values (-0.707107, 0.408248 on the first row in both). The names
are `bcs_fo.L[k]` and `bcs_fo.Q[k]`. No extra `brms_lp_check()` row was
added for this: the polynomial and the treatment codings go through the
same `model.matrix(tt, mfc, contrasts.arg = )` call, and a row would cost
another Stan compile to test the contrast machinery of `stats` rather
than anything in frmtmb.

### After the final install

Every file this change touches, rerun against the reinstalled lane
library (`dev/csfactor-log/batch-final.txt`):

    RESULT test-cs-factor.R       pass=46  fail=0 error=0 skip=0
    RESULT test-compat.R          pass=639 fail=0 error=0 skip=0
    RESULT test-bracket-access.R  pass=33  fail=0 error=0 skip=0
    RESULT test-v17.R             pass=36  fail=0 error=0 skip=0
    RESULT test-review-v29.R      pass=138 fail=0 error=0 skip=0
    RESULT test-brms-agreement.R  pass=168 fail=0 error=0 skip=2
    RESULT test-brms-shapes.R     pass=72  fail=0 error=0 skip=0
    RESULT test-ordinal.R         pass=109 fail=0 error=0 skip=0

`test-brms-agreement.R`'s two skips are the gated distributional-gaussian
and `mo()` rows, neither about `cs()`.

### R CMD check --as-cran, once

`dev/csfactor-log/rcheck.txt`, tarball
`dev/csfactor-check/frmtmb_0.64.0.tar.gz`, pandoc and TinyTeX on PATH,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`:

    Status: 1 NOTE
    * checking HTML version of manual ... NOTE
      Skipping checking math rendering: package 'V8' unavailable
    * checking tests ... Running 'testthat.R' [391s] OK
      [ FAIL 0 | WARN 17 | SKIP 307 | PASS 9124 ]
    * checking re-building of vignette outputs ... [205s] OK
    * checking examples with --run-donttest ... [43s] OK

The V8 NOTE is the one lane-rules.md says to expect. The 17 warnings are
the same 17 the 180-file run reports, in `test-autocor-cond.R` (1),
`test-brms-families.R` (12), `test-brms-names.R` (3) and
`test-nl-rtmb-scope.R` (1); none of those files touches `cs()`.

The first attempt aborted at `* checking package dependencies ... ERROR
Package suggested but not available: 'drmTMB'`, which stops the check
before it does any work. drmTMB is on this machine only in `rellib-r3`,
so the check was rerun with `rellib-r3` LAST on `R_LIBS` (read-only, and
`R CMD check` installs into its own `.Rcheck` library, never into an
`R_LIBS` entry) rather than with `_R_CHECK_FORCE_SUGGESTS_=false`, so
`test-drmtmb-agreement.R` and the sibling-extension tests actually ran.
That is why the check reports 9124 passes where the ungated 180-file run
reports 11923: the check run is `NOT_CRAN=false`, so every
`skip_on_cran()` tier skips (307 skips against 152).

The first `R CMD build` of this lane wrote `frmtmb_0.64.0.tar.gz` into
`C:/Users/adf44/source/r/`, the SHARED parent of every lane's worktree,
because `R CMD build` writes to the working directory. Every lane would
write that same name there. The build was redone inside
`dev/csfactor-check/`, and the stray tarball in the shared directory was
left alone rather than deleted, since another lane may be reading a file
of that name right now (rtmb-pitfalls.md item 20).

### Gated brms parity files

Every gated brms file, run with `FRMTMB_BRMS_FIT_TESTS=true` after the
likelihood file had warmed the Stan cache, four processes at a time
(`dev/csfactor-log/batch-gated.txt`). These are the files whose UNGATED
run reports `pass=0` or a large skip count, so an ungated suite says
nothing about them:

    RESULT test-brms-suite-standata.R        pass=87  fail=0 error=0 skip=0
    RESULT test-brms-suite-brm.R             pass=23  fail=0 error=0 skip=0
    RESULT test-brms-suite-priors.R          pass=36  fail=0 error=0 skip=0
    RESULT test-brms-suite-brmsterms.R       pass=5   fail=0 error=0 skip=0
    RESULT test-brms-suite-brmsformula.R     pass=16  fail=0 error=0 skip=0
    RESULT test-brms-suite-families.R        pass=84  fail=0 error=0 skip=0
    RESULT test-brms-suite-data-helpers.R    pass=6   fail=0 error=0 skip=0
    RESULT test-brms-suite-brmsfit-helpers.R pass=3   fail=0 error=0 skip=0
    RESULT test-brms-suite-emmeans.R         pass=11  fail=0 error=0 skip=0
    RESULT test-brms-suite-methods.R         pass=162 fail=0 error=0 skip=0
    RESULT test-brms-priors.R                pass=103 fail=0 error=0 skip=0
    RESULT test-brms-methods.R               pass=979 fail=0 error=0 skip=0
    RESULT test-brms-agreement.R             pass=187 fail=0 error=0 skip=0

    BATCH produced a RESULT line for 13 of 13 files

No file skipped anything. `test-brms-methods.R` is the one that matters
most here: ungated it is 16 pass and 46 skips, and gated it fits the brms
counterpart of every shape in its matrix, the `y ~ x + cs(z)` row among
them, and compares `fixef()`, `ranef()`, `posterior_epred()`,
`conditional_effects()`, `hypothesis()` and `residuals()` against brms's
own. All 979 pass.

# Punch round 1, 2026-09-29

Review `dev/reviews/2026-09-29-csfactor.md`, verdict NOT MERGEABLE on
three RECORD defects, with the code holding on every construction it
tried. Measurements below are `dev/csfactor-p1.R`, log
`dev/csfactor-log/p1.txt`, on the reviewer's seeds.

## 1. `?frm` said `cs()` was not available in a multivariate fit. False.

Removed. `bf(yo ~ x + cs(f)) + bf(yo2 ~ z + cs(f))` fits, and no such
refusal exists anywhere in `R/frame.R`; the sentence was written from
nothing. Measured here on my own construction (seed 77, n = 300, both
responses pure noise, so the log likelihood is not the reviewer's):

    npar 14
    fixef():     yo_Intercept[1..2] yo2_Intercept[1..2] yo_x yo2_z
                 yo_fcb[1..2] yo_fcc[1..2] yo2_fcb[1..2] yo2_fcc[1..2]
    variables(): bcs_yo_fcb[1..2] ... bcs_yo2_fcc[1..2]
    logLik(mv)                -650.9878
    logLik(u1) + logLik(u2)   -650.9878
    difference                 5.273137e-09

With no `rescor` the multivariate likelihood IS the product of the two, so
the last line is an IDENTITY and the 5.3e-09 is the optimizer's noise; the
reviewer's independent run gave -1.3e-09 on its own construction. One bad
predictor out of two is still refused, naming `'yo.mu'`.

`?frm` now says the opposite of what it said, with the coefficient
spelling, and `test-cs-factor.R`'s "cs() works in a multivariate fit, in
every response" pins the names, the sum identity and the one-bad-response
refusal, so the page cannot drift back.

## 2. The rank test was blind to `mo()`. Fixed, not only documented.

`yo ~ mo(m) + cs(m)`, `m` an ordered factor with 4 levels, fitted with
3 extra degrees of freedom buying 6.4e-09 of log likelihood, two negative
covariance eigenvalues and 9 of 9 NaN standard errors (the reviewer's
`dev/csfactor-rev-gap.R`, seed 1907). Pre-existing, and exactly the
symptom defect B exists to remove.

The reason is structural: `X`'s `mo()` column is an all-zero placeholder
the objective fills later, so it contributes nothing to
`rank(cbind(1, X))` and every `cs()` column then raises the rank.

The fix is a second test rather than a name-based special case, because
the span the `mo()` column can occupy IS known at frame time. The
objective puts `b * D * cumsum(zeta)[codes + 1]` there. Over the simplex
those vectors are the functions of `codes` that vanish at code 0, which is
the treatment-dummy basis of the codes, times the interaction multiplier
for `mo(x):z`. So `mo_indicator_basis()` builds that basis and
`check_cs_identified()` refuses when

    rank(cbind(1, X, Z, Dmo)) == rank(cbind(1, X, Z))   and
    rank(cbind(1, X, Dmo))    >  rank(cbind(1, X))

The first condition says no simplex escapes the cs span. The second keeps
the refusal ABOUT `cs()`: a `mo()` term the population-level design alone
already absorbs is a different, pre-existing problem and is left to report
itself. Treating the whole indicator basis AS the mo column would have
been the false-alarming version, and the no-refusal rows below are what
separate the two.

| shape | now | why |
|---|---|---|
| `yo ~ mo(m) + cs(m)` | REFUSED | `cs(m)` spans every function of `m` |
| `yo ~ mo(m)` alone | FITTED, logLik -329.0959 | |
| `yo ~ cs(m)` alone | FITTED, logLik -324.8122 | |
| `yo ~ mo(m) + cs(x)` | FITTED, -327.0273 | `x` is continuous |
| `yo ~ mo(m) + cs(f)` | FITTED, -328.3364 | `f` is unrelated to `m` |
| `yo ~ mo(m) + cs(mc)`, `mc` a 2-level COARSENING of `m` | FITTED, -328.6165 | `cs(mc)` spans 2 of the 4 values, so the monotone shape stays free |
| `yo ~ mo(m):z + cs(m)` | FITTED, -324.699418765, df 11, 0 of 9 non-finite standard errors, smallest vcov eigenvalue 0.003151 | `z` times a function of `m` is not a function of `m`. CORRECTED after punch round 2: the evidence that this one is identified is the finite standard errors and the positive-definite covariance, NOT the likelihood gain. An `mo()` term costs `D - 1` simplex coordinates plus its coefficient, so against `cs(m)` alone (-324.812229967, df 8) the interaction buys 0.11 units over 3 df, an LR of 0.23, which says nothing either way. The first write-up called it 1 df and read the 0.11 as a real effect |

`mi()` and `me()` are NOT covered and cannot be by this route: their
column is filled with observed-or-latent VALUES rather than a function of
a data column, so there is no basis to test against.
`me(x, sdx) + cs(x)` is accepted and looks identified (logLik -776.7281,
0 of 5 standard errors NaN); `mi(x) + cs(x)` is accepted and UNSETTLED.
`?frm` and `NEWS.md` both say so, in place of implying the check is
exhaustive.

## 3. The `droplevels` bullet was wrong in both directions

Rewritten in "What I decided NOT to do" above, with the withdrawal on the
page rather than only in conversation. The behavioral consequence the
reviewer measured is now LOUD.

`influence()` fills its table by name, so a deletion refit that succeeds
with one column fewer leaves those cells `NA` and `cooks.distance()` `NA`
for that unit, silently. Measured before and after on seed 90291, n = 120,
`yo ~ x + cs(f)` with level `c` holding exactly one row (row 120):

    before: 1 of 120 cooks.distance() NA, no warning and no message
    after:  1 of 120 NA, and

    Deleting unit '120' left a design without 'bcs3_1 (cs fc, boundary 1)',
    'bcs3_2 (cs fc, boundary 2)', so those cells are NA and
    cooks.distance() is NA for that unit. A factor level holding one row
    of its own disappears from the deletion subset, and the refit then
    estimates one coefficient fewer.

The table's own column names are the internal ones, and `bcs3_1` does not
say which LEVEL vanished, so `influence_coef_labels()` maps the `bcs<j>`
part back through the frame's cs entries and the warning names both.

The check is on the COEFFICIENT NAMES a refit came back with, so it is
not about `cs()` at all: any factor level a deletion empties triggers it.
The re-check reviewer measured `yo ~ x + f`, a plain factor with no
`cs()` term, warning once where 0.64.0 was silent with one `NA`. That is
the intended reach and `NEWS.md` says so in a clause; `cs()` only makes
the case easier to meet, because a level there carries `K - 1`
coefficients instead of one.

The warning is deliberately a WARNING and nothing more. Lane wt-thresrefit
is adding a failed-refit counter to `influence()` and a `thres_pin`
argument to `assemble_frame()`; neither reaches this case, because the
refit does not FAIL and `thres_pin` pins thresholds rather than `xlevels`.
Building the same machinery here would collide at consolidation. The
FULLER fix is to carry the fit's `xlevels` into a deletion refit the way
that lane carries the threshold count, and it has a catch worth
recording: a pinned `xlevels` puts the absent level's all-zero dummy back
into the design, which the new constant-column refusal would then reject,
so a pinned refit would have to be exempted from
`check_cs_identified()`. That is a decision for whoever merges the two
lanes, not a change to make blind from here.

The influence hunk is three lines in the loop, one warning block and one
helper, all inside `R/influence.R`, so it should not conflict with a
failed-refit counter in the same function.

## Nits

- **The refusal advice is now conditional on what the formula holds.**
  Before, every aliased-column refusal said "Write 'x + cs(x)' as 'cs(x)'
  alone", which on `poly(x, 2) + cs(x)` names a term that is not in the
  formula and, if followed, deletes the quadratic. Now the message names
  the aliased population-level COLUMNS, from the least-squares weights of
  the cs column on `cbind(1, X)` (exact, because the column is in that
  span), and branches three ways:

  | formula | columns named | advice |
  |---|---|---|
  | `x + cs(x)` | `'x'` | write it as `cs(x)` alone |
  | `poly(x, 2) + cs(x)` | `'poly(x, 2)1'` | no term spells `x` alone; `I(x^2) + cs(x)` reaches the same maximum with one parameter fewer |
  | `s(x) + cs(x)` | `'s(x).fx1'` | that column is a SMOOTH's unpenalized linear part, no spelling keeps both, drop one |
  | `x + z + cs(I(x + z))` | `'x', 'z'` | as the poly branch |

  The smooth branch keys on the `.fx<j>` suffix `R/frame.R` gives a
  smooth's unpenalized columns. `?frm` and `NEWS.md` record that
  `s(x) + cs(x)` has no rewrite keeping both, which is a genuine narrowing
  and was undocumented.
- **The rank tolerance is documented.** `?frm` now says the rank is
  `qr()`'s at its default `tol = 1e-7` over unit-norm columns, so a merely
  ill-conditioned design is refused past a condition number of about 1e7.
  The reviewer measured the boundary: `eps = 1e-6` (condition number
  2.2e6) fits and `eps = 1e-7` (2.3e7) is refused.
- **`?variables` no longer claims the draws object's columns follow the
  same convention.** It names the exception, `tau_raw_1` and `bcs2_1`, and
  says `fixef()` on the draws object does report brms's rows. The filed
  defect is widened above and in `dev/test-backlog.md`.
- **The one absolute constant in `test-cs-factor.R` is gone.**
  `expect_gt(max(abs(P1[1, ] - P3[1, ])), 0.3)` is now
  `> 0.5 * max(abs(P3[3, ] - P3[1, ]))`, a ratio to the separation between
  level a and level c that the run itself measures. Every other bound in
  the file is a ratio to a measured quantity or a SIGMA COUNT
  (`max(abs(pr - p) / mcse) < 5`, dimensionless), so the claim above is
  now literally true. The reviewer was right that the original was
  defensible as a discrimination floor and that the sentence was not.
- **`dev/test-backlog.md` line 970 reflowed**, and the two pre-existing
  defects the review found are filed there: the multivariate
  `fitted(newdata = )` error with its construction, and the widened
  `variables()` filing.

## Punch round 1 counts

`test-cs-factor.R` grew from 46 assertions to 74, in five new blocks
(multivariate, the `mo()` refusal with its six no-refusal controls, the
message branches, the unused level, the `influence()` warning). SEEN TO
FAIL again on the reference build with the final file:

| build | pass | fail | error | skip |
|---|---|---|---|---|
| rellib-r3 (0.64.0) | 17 | 37 | 5 | 0 |
| lane | 74 | 0 | 0 | 0 |

14 of the 15 blocks fail on the reference build. The one that passes on
both is "the identifiability check does not fire on a separable model",
the false-alarm control, which must pass on both.

The whole core suite was rerun, because `R/influence.R` is reached by
every path that refits and the `mo()` branch is reached by every `cs()`
frame. One file per R process, eight at a time, `NOT_CRAN=true`
(`dev/csfactor-log/suite-p1.txt`):

    RESULT lines: 180 of 180
    pass=11954 fail=0 error=0 skip=152 warn=17

Diffed line by line against the pre-punch run, ONE file's count changed
and it is the new one:

    < RESULT test-cs-factor.R pass=43 fail=0 error=0 skip=0 warn=0
    > RESULT test-cs-factor.R pass=74 fail=0 error=0 skip=0 warn=0

Every other file is identical, skips and warnings included, so neither the
`influence()` warning nor the `mo()` refusal moved anything else. (The
pre-punch run of `test-cs-factor.R` reads 43 rather than 46 because it
predates the ordered-factor block.)

The gated log-density file was rerun too, on the now-warm Stan cache, to
confirm the new refusal branches do not fire on any of its 60 shapes
(`dev/csfactor-log/gated-likelihood-p1.txt`):

    RESULT test-brms-likelihood.R pass=440 fail=0 error=0 skip=0 warn=0

the same 440 as before the punch round, with nothing skipped.

### R CMD check --as-cran, after the punch round

Rebuilt and rechecked (`dev/csfactor-log/rcheck-p1.txt`), same environment
as before, `rellib-r3` last on `R_LIBS` so drmTMB is present:

    Status: 1 NOTE
    * checking HTML version of manual ... NOTE
      Skipping checking math rendering: package 'V8' unavailable
    * checking tests ... Running 'testthat.R' OK
      [ FAIL 0 | WARN 17 | SKIP 307 | PASS 9124 ]
    * checking re-building of vignette outputs ... [174s] OK
    * checking examples with --run-donttest ... [25s] OK

Identical to the pre-punch check in every line, the V8 NOTE included.
The pass count does not move because the check runs with
`NOT_CRAN=false` and `test-cs-factor.R` opens with `skip_on_cran()`, so
the whole file skips there; the 180-file run above is what exercises it.

# Punch round 2, nits only, 2026-09-29

Verdict on re-check was MERGEABLE. Five nits, all closed. Measurements are
`dev/csfactor-p2.R`, log `dev/csfactor-log/p2.txt`, seed 1907, n = 300.

1. **The fallback advice branch no longer advertises a `poly()` that is
   not there.** On `x + z + cs(I(x + z))` it named `('x', 'z')` correctly
   and then pasted the polynomial worked example. `cs_alias_advice()` now
   has four branches, each keeping its own example: a plain column, a
   `.fx<j>` smooth column, a `poly(` column, and a generic fallback that
   names nothing but the columns already listed ("No population-level term
   spells 'I(x + z)' on its own, so there is nothing to delete by that
   name: drop either the cs() term or the population-level term(s)
   carrying the column(s) named above."). The message-branch test now pins
   the fallback text and asserts that `poly` and `I(x^2)` are ABSENT from
   it, and that the poly branch does not carry the fallback wording.

2. **The `mo()` refusal names the spelling that keeps an interaction.**
   `mo(m) * z + cs(m)` is refused, because `*` expands to the main effect
   plus the interaction and the main effect is the unidentified part;
   `z + mo(m):z + cs(m)` fits, df 12, logLik -324.125678297, 0 of 10
   non-finite standard errors, smallest vcov eigenvalue 0.002812. The
   message now says to write the second rather than the first, and the
   test asserts both the refusal of `mo(m) * z + cs(m)` and the fit of
   `z + mo(m):z + cs(m)`.

3. **The `mo(m):z` degrees of freedom were wrong in the round-1 table**
   and are corrected in place above: an `mo()` term costs `D - 1` simplex
   coordinates plus its coefficient, so the interaction buys 0.11 units
   over 3 df, an LR of 0.23, which is no evidence at all. What actually
   shows that row is identified is the finite standard errors and the
   positive-definite covariance (smallest eigenvalue 0.003151), and the
   table now says so.

4. **`y ~ mo(m) + m`, with no `cs()` anywhere, is filed.** It fits in
   silence and it is not identified: df 8 and logLik -328.534572784
   against `m` alone at df 5 and -328.534572790, so three extra parameters
   buy 6e-09; smallest covariance eigenvalue -6441, 6 of 6 standard errors
   non-finite, and the only sign of it is R's own
   `In sqrt(diag(V)) : NaNs produced`. With `m` written as an UNORDERED
   factor the eigenvalue is -6724 and the standard errors happen to stay
   finite, so nothing at all is visible. Measured identically on both
   builds. The brief said 4 non-finite standard errors; on this
   construction it is 6 of 6 with the ordered factor and 0 of 6 with the
   unordered one, and the invariant across both is the negative
   eigenvalue, so that is what the backlog entry leads with.
   `check_cs_identified()`'s comment said such a case "is left to report
   itself"; it now says nothing reports it today and points at the
   backlog. Not refused here on purpose: it is a statement about `mo()`,
   the `mo()` column is a zero placeholder at assembly, and the fix
   belongs with a rank check over the FILLED design.

5. **The `influence()` warning is general, and said so.** It keys on the
   coefficient names a refit returned, so any factor level a deletion
   empties triggers it; the reviewer measured `yo ~ x + f` with no `cs()`
   warning once where 0.64.0 was silent with one `NA`. Recorded above and
   in `NEWS.md` in one clause.

Files touched this round: `R/frame.R` (`cs_alias_advice()` branches, the
`mo()` message, the `check_cs_identified()` comment), `man/frm.Rd` (from
roxygenising), `NEWS.md`, `dev/test-backlog.md`,
`tests/testthat/test-cs-factor.R`, this file, plus `dev/csfactor-p2.R`
and `dev/csfactor-log/p2.txt`.

## Punch round 2 count

`test-cs-factor.R` grew from 74 assertions to 84, in the two blocks the
nits touched (the message branches and the `mo()` controls):

    RESULT test-cs-factor.R pass=84 fail=0 error=0 skip=0 warn=0
