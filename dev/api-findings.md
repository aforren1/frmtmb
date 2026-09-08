# wt-api: four nonlinear usability defects

Worktree `C:/Users/adf44/source/r/frmtmb-wt-api`, branch `wt-api`, from
0.53.0. Private library
`.../scratchpad/ap-lib`; every script prints `packageVersion("frmtmb")`
before it does anything.

Reproduction data set (task brief): 40 subjects, `w` 1..45 by 0.25,
`logf = b_i - chi_i log w + 1.1 exp(-0.5 ((w - mu_i) / 1.6)^2)`,
`I ~ rexp(rate = 1 / exp(logf))`; 7080 rows. The reference fit

```r
frm(bf(I ~ apo - chi * logw + exp(lamp) * exp(-0.5 * ((w - (10 + pk)) / exp(lsig))^2),
       apo ~ 1 + (1 | id), chi ~ 1 + group + (1 | id), pk ~ 1 + (1 | id),
       lamp ~ 1, lsig ~ 1, nl = TRUE),
    family = exponential(link = "log"), data = d)
```

fits in **2.32 s** and recovers the truth (chi group effect 0.4389
against 0.35, apo 1.062 against 1.2, pk near 0 against a centre of 10).

---

## Defect 1: reserved nonlinear parameter names

### Reproduced (0.53.0)

| spelling | result |
| --- | --- |
| `bf(I ~ mu - chi * logw, mu ~ 1 + (1\|id), chi ~ 1 + group, nl = TRUE)` | `Error: object 'mu' not found` from `eval(predvars, data, env)` |
| same, through `par_template()` | same error, so it is a SPEC-time fault, not a fit-time one |
| `bf(I ~ b - chi * logw, b ~ 1 + (1\|id), ...)`, no `start` | **fits normally** |
| same plus `start = list(b = 1)` | `Error: start$b must have length 40` |

`theta`, `beta`, `betad`, `sigma`, `shape`, `nu`, `phi`, `miss`, `z`
all fit normally as nonlinear parameter names under
`exponential(link = "log")`, which has only `mu` as a dpar. That is the
clue: the fatal set is the family's own dpars, nothing else.

### Cause, measured

`R/parse.R:1503`

```r
nlpars <- if (length(nl_dpars)) {
  setdiff(c(names(pforms), nl_dpars), fam[["dpars"]])
} else character(0)
```

A name that is one of the family's dpars is subtracted out, so
`mu ~ 1 + (1 | id)` is read as a formula for the dpar `mu`, not as a
nonlinear parameter. `nl = TRUE` then also gives `mu` a BODY (the
response right-hand side, `parse.R:1494`). In `nl_dpar()`
(`R/parse.R:1296`) the body's free names split three ways, and
`nl_dpar_refs` deliberately excludes the parameter's own name:

```r
nl_dpar_refs = setdiff(intersect(vars, dparnames), c(name, pars)),
datavars     = setdiff(vars, pars),
```

so `mu` inside `mu`'s own body lands in `datavars` and is asked of the
model frame. Hence "object 'mu' not found", naming nothing.

The `b` case is different in kind. `b` IS a nonlinear parameter there
and its coefficients are `beta["b_(Intercept)"]`. What collides is only
the `start`/`newparams` namespace, whose keys are the par-template
COMPONENT names.

### The unsafe set, enumerated

* **Fatal, spec time** (refused after this lane): a family dpar name
  used inside its own nonlinear body. On `nl = TRUE` that is always the
  location dpar, so `mu` for every family that `nl = TRUE` accepts
  (`nl = TRUE` already requires a single-`mu` location family,
  `parse.R:1489`). Through `nlf()` it is any dpar given a body that
  names itself.
* **Silently reinterpreted, not a defect**: a family dpar name used in
  ANOTHER dpar's body is an `nl_dpar_ref` and reads that parameter's
  per-row value. `bf(y ~ sigma * x + a, a ~ 1, nl = TRUE)` under
  `gaussian()` fits, with `sigma` the dpar's value, not a new
  parameter. This is the documented variance-function extension
  (`parse.R:1284`) and must keep working.
* **Ambiguous at `start =` / `newparams =`**: the par-template
  component names, measured across model shapes to be exactly
  `beta`, `betad`, `b`, `theta`, `thetaac`, `thetar`, `miss`
  (`R/frame.R:2382-2420`). A nonlinear parameter may carry any of these
  names and fit; only the start list misroutes.

### Namespacing versus refusal, measured

Namespacing an nl parameter called `mu` means renaming it internally
and mapping back on every user-facing surface. Measured cost: **57
lines across 13 files key on a dpar name** (`grep '\[\["dpar"\]\]\|\$dpar'
R/*.R`) and **39 lines across 9 files read `nlpars`**. The parsed spec's
`dpars` list is itself keyed by parameter name, so `mu` the bodied dpar
and `mu` the nonlinear parameter cannot both be in it: the rename would
have to survive priors, `par_template`, `start`, `fixef`, `coef`,
`ranef`, `predict`, `conditional_effects` and the two sampling
extensions, and any surface missed returns a silent wrong name.

The reference implementation settles it. brms 2.23.0, same collision:

* `bf(y ~ mu - 2 * x, mu ~ 1, nl = TRUE)` **silently discards the
  nonlinear body** and generates `mu = rep_vector(0.0, N); mu +=
  Intercept;` -- an intercept-only model, no warning.
* `bf(y ~ sigma - 2 * x, sigma ~ 1, nl = TRUE)` under `gaussian()`
  errors with `No non-linear parameters specified.`

So brms has no namespacing to be compatible with, and on the `mu` case
it is strictly worse than an error. **Decision: named refusal at parse
time.** Namespacing is refused as disproportionate and as a divergence
from brms with no upside.

---

## Defect 2: two spellings of one coefficient name

### Spelling census, measured on the reference fit

| surface | spelling |
| --- | --- |
| `names(fixef(fit))` | `apo`, `chi`, `pk`, `lamp`, `lsig` (dpar KEYS, not coefficient names) |
| `names(unlist(fixef(fit)))` | `apo.(Intercept)`, `chi.group`, ... |
| `rownames(vcov(fit))`, `vcov(full = TRUE)` | `apo_(Intercept)`, `chi_group`, ... |
| `rownames(confint(fit))` | underscore |
| `names(par_template(fit)$beta)` | underscore |
| `names(fit$estimates$beta)` | underscore |
| `outer_par_names()`, `estimated_coef_names()` | underscore |
| `hypothesis(fit, "chi_group = 0")` | underscore (the dot spelling is refused) |
| `set_prior(nlpar = , coef = )` | dpar and bare column, separately |
| `summary(fit)$coefficients` | a LIST keyed by dpar, bare column rownames |
| `coef(fit)` | list keyed by dpar, then by grouping factor, bare columns |
| `insight::find_parameters()` | dpar keys |
| `marginaleffects::get_coef()` | underscore, from `estimated_coef_names()` (`R/interop.R:131`), and `get_vcov()` is `vcov()`, so both halves of a marginal effect agree |
| `emmeans` | never sees a nonlinear parameter: `emm_mu_linpred()` refuses a fit whose mu has a body (`R/interop.R:114`), and the compat registry records that refusal (`R/compat.R:1373`) |
| `frmtmb.sample::fixef.frmtmb_draws()` | underscore, parenthesis-free (`par_name_bare(estimated_coef_names(fit))`), matching the draws matrix, `variables()` and `hypothesis()` |
| `frmtmb.sample::ranef.frmtmb_draws()` | inherits core's keys, indexes by position |

**Only one surface produces the dot spelling, and frmtmb does not
choose it**: base R's `unlist()` composes a list key with an element
name using `.`. `fixef()` returns a LIST keyed by dpar whose elements
are named by DESIGN COLUMN.

### Why a separator change cannot reconcile them

The canonical name (`R/frame.R:2108-2116`) drops the location dpar:

```r
cn <- colnames(X)
if (!identical(dp[["name"]], "mu")) cn <- paste(dp[["name"]], cn, sep = "_")
if (length(spec$responses) > 1) cn <- paste(resp$resp_name, cn, sep = "_")
```

Measured on a plain `y ~ x + (1 | g)` gaussian fit:

```
unlist(fixef(f)):  mu.(Intercept) | mu.x | sigma.(Intercept)
rownames(vcov(f)): (Intercept)    | x    | sigma_(Intercept)
```

They differ in the separator AND in whether `mu` is named at all. So
this is not one name with two separators; it is a base-R composite key
against a frmtmb coefficient name.

### Blast radius of moving `unlist(fixef())`, measured

`unlist()` IS an internal generic and does dispatch (verified: a
`unlist.<class>` method on a classed list is called), so the change is
technically possible. Its cost is not.

* 240 literal `[["a.b"]]` index sites in the monorepo.
* **19 files carry both `unlist(fixef(...))` and literal dot-name
  indexing**: core `tests/testthat/test-mm.R` (11), `test-sparsex.R`
  (5), `test-vignette-wiener.R` (6), `vignettes/case-studies.Rmd` (7),
  and five extensions -- frmtmb.eam (`test-defects.R` 3,
  `test-family.R` 8, `test-gddm-recovery.R` 9, `test-sampling.R` 1,
  `test-surface.R` 10, `test-variability.R` 6, `ddm.Rmd` 9,
  `gddm.Rmd` 12), frmtmb.latent (`test-hmm.R` 10, `latent.Rmd` 7),
  frmtmb.learn (`learning.Rmd` 5), frmtmb.ode (`test-ode-events.R` 12,
  `test-ode-tv.R` 3, `test-ode.R` 6), frmtmb.sample
  (`test-simulators.R` 4).

Roughly 134 index sites across five extensions. The lane's scope is
core only, with the smallest possible extension change. **Refused: the
`unlist()` spelling does not move.**

### What does change

The defect that is real and fixable is that the package POINTS the user
at the mismatch and emits it:

1. `?fixef` says `unlist(fixef(fit))` is "the vector `confint()` and
   `hypothesis()` name their rows by". That is false on every model.
2. `frm_bootstrap()`'s default statistic is
   `function(f) unlist(fixef(f))`, so a bootstrap standard error cannot
   be lined up with `sqrt(diag(vcov(fit)))` by name -- which is the
   whole point of the comparison. Measured blast radius of moving THIS
   one: a single core assertion (`tests/testthat/test-boot.R:27`) and
   nothing in any extension or vignette (every vignette passes its own
   `FUN`).
3. There is no supported route from `fixef()` to the canonical vector.
   Added `fixef(object, flatten = TRUE)`.

---

## Defect 3: ranef() blocks are not addressable

### Reproduced

Three `(1 | id)` terms on three nonlinear parameters:

```
names(ranef(fit)):   id | id | id
  [[1]] dim 40x1  term "apo: 1 | id"
  [[2]] dim 40x1  term "chi: 1 | id"
  [[3]] dim 40x1  term "pk: 1 | id"
ranef(fit)[["id"]]   -> the apo block, silently
```

### History that constrains the fix

`dev/reviews/2026-09-06-sample-ce.md` was read and contains **no mention
of ranef** (zero hits): that review is about the conditional-effects
draws band, `estimate__` on non-identity families, and `hypothesis()`'s
affine probe. The constraint therefore comes entirely from the NEWS
entry.

`NEWS.md` 0.52.0 carries a BEHAVIOR CHANGE entry: `ranef()` is keyed by
the GROUPING FACTOR, as brms and lme4 key it, deliberately away from
the block label, because `ranef(fit)$Subject` used to be `NULL` where
`coef(fit)$Subject` was a data frame. The block label rides in the
`"term"` attribute, is `VarCorr()`'s key and the `grp` column of
`as.data.frame()`. `R/methods-fit.R:640` appends then names, precisely
so a duplicate label cannot drop a block.

`extensions/frmtmb.sample/R/methods-draws.R:154` indexes core's list BY
POSITION and copies `names(per[[1]])` wholesale, with a comment saying
why. So it follows whatever core names the list, and needs no change as
long as the list stays a positionally indexed list with names.

---

## Defect 4: a flat nonlinear term reports the wrong diagnosis

### Reproduced, with a correction to the brief

Dropping the `10 +` from the peak centre. The fit does NOT stay at the
start point: `exp(-0.5 * ((w - 0) / exp(0))^2)` at `w = 1` is 0.607, so
there IS a gradient at the start. The optimizer drives `lsig` down and
`pk` away, the bump underflows, and the model becomes flat DURING the
optimization. Measured at the reported optimum:

```
gradient  apo_(Intercept)  -5.18e-04
          chi_(Intercept)   1.14e-03
          chi_group         7.09e-04
          pk_(Intercept)    0.00e+00   <- flat
          lamp_(Intercept)  0.00e+00   <- flat
          lsig_(Intercept)  0.00e+00   <- flat
          theta_1          -9.90e-04
          theta_2          -1.28e-04
          theta_3           0.00e+00   <- flat (the pk group sd)
```

`summary()` reports NaN standard errors on ALL SIX coefficients.

### What diagnose() reports today

```
Max |gradient|: 0.001138 at chi_(Intercept)
Hessian positive definite: FALSE
Non-finite standard errors: apo_(Intercept), chi_(Intercept), chi_group,
  pk_(Intercept), lamp_(Intercept), lsig_(Intercept), theta_1, theta_2, theta_3
```

Every parameter is named, so none is. The four flat directions are not
distinguished from the five identified ones, and `worst_grad` points at
`chi_(Intercept)`, which is fine.

The fit-time warning is `R/utils.R:300` `warn_nonfinite_cov()`: "the
covariance could not be recovered from the Hessian and the model is
probably overparameterized". The model is not overparameterized; four
directions are flat.

`fit$obj$he()` is unavailable here ("Hessian not yet implemented for
models with random effects"), so an empty-Hessian-row test has to come
from the gradient.


---

# What landed

## Defect 1

* `R/parse.R:1524-1563` `parse_one_response()`: a name carrying BOTH a
  parameter formula and a nonlinear body is refused, naming the family
  and what it reserves. No data is needed, so `par_template()` refuses
  it too. CORRECTION, from the review: the second branch of that check,
  for a name that is NOT a family dpar, is UNREACHABLE and was
  described here as delivered behavior. Re-measured on all four routes
  that could reach it -- `bf(nl) + lf(a) + nlf(a)`,
  `bf(nl) + nlf(a) + lf(a)`, `bf(nl, a ~ 1) + nlf(a)`, and a duplicate
  inside one `bf()` -- and every one is refused earlier, the first
  three by `R/bf.R:338`, `:351`, `:358` ("nlf() sets 'a', which the
  bf() it is added to already sets") and the fourth by "Duplicated dpar
  formula: 'a'". So the only member `both` can hold is the location
  dpar, and the first branch always takes it. The branch is KEPT as a
  guard, with a comment at `R/parse.R:1546` saying it is unreachable
  while those refusals stand: without it, a non-dpar name reaching this
  point would be told it is a distributional parameter, which would be
  false. No test covers it, because none can.
* `R/frame.R:799-843` `check_nl_self_reference()`, called from
  `assemble_frame()` at `R/frame.R:1057` before
  `resolve_nl_dpar_refs()`: a nonlinear body naming its OWN parameter
  with no column of that name in `data` is refused, with the
  family-dpar sentence added when it applies. The check is made where
  the data is known so that a real column of that name still wins, as
  it does for a dpar reference.
* `R/fit.R:1879-1917` `nl_start_collision_msg()`, used at
  `R/fit.R:1961-1979` (`make_start()`) and `R/simulate-new.R:594-610`
  (`frm_simulate(newparams = )`): when a `start` component name is also
  a nonlinear parameter of the model, the component reading still
  applies (it is the only one the template can express) but the message
  names both readings and gives the spelling that reaches the
  parameter. A length error carries it; a length that happens to match
  warns instead of setting the wrong object silently.

## Defect 2

* `R/methods-fit.R:531-604` `fixef(object, flatten = TRUE)`: the
  canonical named vector, read straight off `fit$estimates`, which
  already carries the template's names, so the dpar-prefix rule is not
  restated anywhere it could drift.
* `R/methods-fit.R:531-549`: the `?fixef` claim that
  `unlist(fixef(fit))` is "the vector confint() and hypothesis() name
  their rows by" is replaced by what the two spellings are and why.
* `R/bootstrap.R:35` `frm_bootstrap()`'s default `FUN` is
  `fixef(f, flatten = TRUE)`. BEHAVIOR CHANGE, blast radius measured:
  one core assertion (`tests/testthat/test-boot.R:27`), nothing in any
  extension, nothing in any vignette (all pass their own `FUN`).
* REFUSED: moving `unlist(fixef())`. Measured above.

## Defect 3

* `R/methods-fit.R:658-716` `ranef_pick()` plus `[[.ranef_frmtmb` and
  `$.ranef_frmtmb`. The list keeps its 0.52.0 grouping-factor keys and
  its order; the KEY gains a second accepted form, the block label. A
  factor name several blocks share is refused and names the labels.
  Positional indexing takes a fast path, so `print()`,
  `as.data.frame()` and frmtmb.sample's position-indexed
  `ranef.frmtmb_draws()` never touch the new code.
* frmtmb.sample UNCHANGED, and verified: `methods-draws.R:154` indexes
  core's list by position and copies `names(per[[1]])` wholesale, and
  all three of its `ranef` assertions
  (`test-sample-direct.R:69-70`, `test-sampling-ported.R:570`,
  `test-reparam.R:545`) are single-block models.

## Defect 4

* `R/confint.R:806-895` `diagnose_flat()`, `flat_pars()` (memoized on
  the fit cache) and `flat_par_note()`.
* `R/confint.R:1004-1012`: `diagnose()` gains `$flat`, computed only
  when the Hessian is not positive definite or a standard error is not
  finite, so a healthy fit spends no gradient evaluations.
* `R/confint.R:1078-1094`: the printed report gains a "Flat directions"
  block.
* `R/utils.R:300-317` `warn_nonfinite_cov()` and `R/utils.R:356-365`
  `solve_joint_precision()` take the fit and name the flat directions
  in place of "probably overparameterized"; `R/fit.R:2050-2064`
  `check_convergence()` does the same at fit time. Callers threaded at
  `R/methods-fit.R:391`, `R/methods-fit.R:403`, `R/confint.R:1787`,
  `R/predict.R:27`.

## Documentation

* `vignettes/frmtmb.Rmd`, nonlinear section: starting values for
  bump-shaped terms, and the reserved-name rule.
* `vignettes/diagnostics.Rmd`, new "NaN standard errors" subsection:
  the two causes and which one `diagnose()` reports.
* `NEWS.md`: a development-version section naming both behavior
  changes.

## Tests

* `tests/testthat/test-nl.R` (+6 tests): one per reserved name, the two
  refusals, the two cases that must keep fitting (a real data column,
  a dpar reference), and every template component name through
  `nl_start_collision_msg()`.
* `tests/testthat/test-api-spellings.R` (+5): the flattened names equal
  `vcov()`'s rows on two model shapes (the review counted correctly;
  "three" here was wrong, and the reviewer confirmed set equality on
  twelve shapes including the multivariate one where only the ORDER
  differs), the list shape is unchanged,
  the bootstrap statistic matches the covariance, and the `unlist()`
  composite is PINNED so the refusal cannot drift into looking
  canonical.
* `tests/testthat/test-methods.R` (+4): label addressing, the refusal
  of an ambiguous factor name, the single-block model unchanged
  (including base `$` partial matching), and position-indexed
  `print()` / `as.data.frame()`.
* `tests/testthat/test-diagnostics-ux.R` (+3): the flat set is named
  and is a strict subset of the NaN standard errors, the same model
  with the centre on the data has none, and a healthy fit never
  populates the cache.
* `tests/testthat/test-boot.R:27`: the moved bootstrap names.

No new test FILE: the suite stays at 131.


---

# Verification

## Core suite, one file per process

131 files on disk, 131 RESULT lines, 131 unique names, no duplicates,
no file without a line (audited with `comm -23` against `ls`).

Totals: **0 failed, 5 errored, 133 skipped, 7393 passed.**

All 5 errors are the SAME harness limitation and were proved
pre-existing by installing the unmodified main checkout (`54d4d92`)
into a second library and running the same files there:

| file | errors | baseline | cause |
| --- | --- | --- | --- |
| `test-compat-register.R` | 1 | 1, identical | `local_mocked_bindings()` needs pkgload; a `library()`-loaded namespace has none |
| `test-influence-plot.R` | 2 | 2, identical | same |
| `test-structure.R` | 2 | 2, identical | same |

`R CMD check` runs the same suite through `tests/testthat.R` and reports
`checking tests ... [347s] OK`, so these three are artifacts of the
per-file runner only.

## Gated tier

`FRMTMB_BRMS_FIT_TESTS=true`, `NOT_CRAN=true`,
`FRMTMB_STAN_CACHE` on a private copy of the union cache
(`ap-stan-cache`, 112 files; the original was copied, never written),
`R_MAKEVARS_USER` on a copy of `dev/stan-cache/makevars-cxx17.mk`. Run
because `test-brms-methods.R` and `test-brms-likelihood.R` read
`ranef()` and `diagnose()`.

| file | failed | error | passed | seconds |
| --- | --- | --- | --- | --- |
| `test-brms-methods.R` | 0 | 0 | 950 | 147 |
| `test-brms-likelihood.R` | 0 | 0 | 394 | 68 |
| `test-brms-agreement.R` | 0 | 0 | 187 | 256 |
| `test-brms-priors.R` | 0 | 0 | 103 | 14 |

**1634 passed, 0 failed, 0 errored.**

## Extension suites

No extension was modified. frmtmb.sample was nonetheless run against
the changed core, because `ranef()` changed shape and its
`ranef.frmtmb_draws()` reads core's list:

| file | failed | error | skipped | passed |
| --- | --- | --- | --- | --- |
| `test-sample-direct.R` | 0 | 0 | 0 | 136 |
| `test-sampling-ported.R` | 0 | 0 | 1 | 208 |
| `test-reparam.R` | 0 | 0 | 0 | 260 |
| `test-draws-methods.R` | 0 | 0 | 0 | 97 |

## Roxygen

Idempotent: two consecutive `roxygenise()` runs give byte-identical
`NAMESPACE`, `man/fixef.Rd`, `man/ranef.Rd`, `man/frm_bootstrap.Rd`
(md5 compared). Only those four files differ from HEAD.

## R CMD check --as-cran, with the manual

`_R_CHECK_CRAN_INCOMING_=false`, `_R_CHECK_FORCE_SUGGESTS_=false`,
pandoc from the RStudio quarto tools directory, TinyTeX prepended to
PATH, `R_LIBS` with `ap-lib` first, and **no `--library=`** and **no
`--no-manual`**. 58 checks.

```
Status: 1 NOTE
* checking HTML version of manual ... [37s] NOTE
Skipping checking math rendering: package 'V8' unavailable
```

Measured cause: `requireNamespace("V8")` is `FALSE` on this machine.
This is the one expected NOTE. `frmtmb-manual.pdf` was built
(542,679 bytes), so the manual half of the check ran.

No ERROR, no WARNING. Every other line is OK, including:

```
* checking S3 generic/method consistency ... OK   <- [[, $, fixef(flatten=)
* checking examples ... [49s] OK
* checking examples with --run-donttest ... [52s] OK
* checking tests ... [348s] OK
* checking re-building of vignette outputs ... [400s] OK
* checking PDF version of manual ... OK
```

The check was run twice: once mid-lane, and again on the FINAL tree
after `diagnose()`'s roxygen was corrected from "four checks" to five.
Both runs give the same one NOTE.

## Diff

`git -C C:/Users/adf44/source/r/frmtmb-wt-api diff --name-only`, 22
files:

```
NAMESPACE  NEWS.md
R/bootstrap.R  R/confint.R  R/fit.R  R/frame.R  R/methods-fit.R
R/parse.R  R/predict.R  R/simulate-new.R  R/utils.R
man/diagnose.Rd  man/fixef.Rd  man/frm_bootstrap.Rd  man/ranef.Rd
tests/testthat/test-api-spellings.R  tests/testthat/test-boot.R
tests/testthat/test-diagnostics-ux.R  tests/testthat/test-methods.R
tests/testthat/test-nl.R
vignettes/diagnostics.Rmd  vignettes/frmtmb.Rmd
```

Untracked: `dev/api-findings.md` (this file) and nothing else. No
`dev/*.log`, no scratch script in the worktree; every script of this
lane lives in the scratchpad under the `ap-` prefix. DESCRIPTION is
untouched. Nothing committed.


---

# Punch round, against dev/reviews/2026-09-06-api.md

Verdict PUNCH, five items plus three optional. Every item was
re-measured here before it was acted on, not taken on the review's word.

## Item 1 (BLOCKER): the flat-direction message asserted a nonlinear cause

The reviewer is right and the reproduction is exact. On the fixture
`frm(y ~ mo(mo) + x)` with `set.seed(6)`, 300 rows, all four `mo` levels
populated (70, 75, 67, 88):

* `nlpars` is `character(0)` -- no nonlinear anything.
* `zeta1_2` is genuinely flat: the objective is
  `220.205072942557450` at the optimum, `220.205072942557933` at
  `+1` and `220.205072942557166` at `-5`. So the DETECTOR is right, and
  it generalizes past the nonlinear case it was built for.
* The message did not generalize with it. BEFORE:

  > The likelihood is FLAT in 1 direction: zeta1_2 ... **A nonlinear
  > term that has left its own support does this; give it a starting
  > value that puts it back (see par_template())**

  A false lead on a model with no nonlinear term, where main shipped
  only a vague one. That is Defect 4's own failure mode relocated.

Fixed by gating the CAUSE, not the diagnosis. `R/confint.R:884`
`fit_has_nlpars()`; used at `R/confint.R:925` (`flat_par_note()`) and
`R/confint.R:1151` (`diagnose()`'s print block). AFTER, same fit:

> ... not identified AT THIS POINT rather than over-parameterized:
> moving it does not change the likelihood at all

and the printed block ends "Moving it does not change the likelihood at
all, so no amount of optimization will pin it down: what the fit reports
for that parameter is wherever the optimizer stopped. par_template()
names it." No `nonlinear`, no `bump`. The nonlinear fit keeps its own
explanation, verified on the peak model.

Regression test at `tests/testthat/test-diagnostics-ux.R`, three cases:
the `mo()` fixture asserts `zeta1_2` is flagged, that the objective does
not move when it moves, and that neither the printed report nor the
`vcov()` warning contains "nonlinear" or "bump"; a second asserts the
nonlinear fit still gets the nonlinear sentence; a third pins the
`theta` label.

## Item 2 (MEDIUM): `flatten = TRUE` is not `hypothesis()`'s spelling

Re-measured, and the reviewer is right on every branch:

```
hypothesis("(Intercept) = 0")        OK   (an accident of R parsing)
hypothesis("x = 0")                  OK
hypothesis("sigma_(Intercept) = 0")  REFUSED: could not find function "sigma_"
hypothesis("`sigma_(Intercept)` = 0") OK
hypothesis("sigma_Intercept = 0")     OK   (hypothesis's own spelling)
```

`hypothesis()` dropped from the claim in three places and replaced by
what is true, with a pointer to the three-vocabulary comment the package
already carried at `R/confint.R:61-66`: `R/methods-fit.R:534-545`
(description), `R/methods-fit.R:578` (the example comment, which still
said "confint() and hypothesis() name their rows by" above a
`flatten = TRUE` call), and `NEWS.md:37`. `man/fixef.Rd` follows from
roxygen.

## Item 3 (MEDIUM): which models the ranef() change breaks

Re-measured. All three shapes give two blocks on one factor and all
three now refuse `ranef(fit)$g`:

| model | block labels |
| --- | --- |
| `y ~ x + (1 + x \|\| g)` | `1 \| g`, `0 + x \| g` |
| `y ~ x + (1 \| g) + (0 + x \| g)` | `1 \| g`, `0 + x \| g` (identical: `\|\|` desugars to this) |
| `bf(y ~ x + (1 \| g), sigma ~ (1 \| g))` | `1 \| g`, `sigma: 1 \| g` |
| `y ~ x + (1 + x \| g)` | ONE block, `$g` unchanged |

Named in `NEWS.md` and in `man/ranef.Rd`'s paragraph
(`R/methods-fit.R:620-630`). Three fixtures added to
`tests/testthat/test-methods.R`: the `\|\|` model (with the two-bar
spelling asserted to give identical labels), the distributional random
effect, and the correlated slope as the control that must NOT change.

## Item 4 (LOW): the dead branch

Verified unreachable myself on all four routes; see the corrected
paragraph above. Kept as a guard with a comment at `R/parse.R:1546`
saying so, and `dev/api-findings.md` corrected.

## Item 5: verified as a non-issue, not skipped

The coordinator re-measured: the merge base of main and `wt-api` is
`028a862`, and
`git diff --name-only 028a862 wt-api -- dev/frequency-domain-todo.md`
is empty, so the branch never touches that file and a three-way merge
keeps main's copy. The reviewer's warning holds for a checkout or a
reset, not for the merge that will actually run. No merge was made, no
commit, and the branch is left where it is on instruction.

## Optional items, all three taken

* **1.2**: `nl_start_collision_msg()` gains an `arg` argument
  (`R/fit.R:1889`, default `"start"`), and `R/simulate-new.R:600,606`
  passes `"newparams"`. A `newparams` error now reads
  `newparams$b must have length 40. \`newparams$b\` sets the
  random-effect vector ...` where it used to answer with `start$b`.
  Both spellings pinned in `tests/testthat/test-nl.R`.
* **2.2**: `@param flatten` (`R/methods-fit.R:556-559`) now says that
  dividing by a `vcov()` diagonal is `NA` for a constant dpar and gives
  the selection that avoids it,
  `intersect(names(cf), rownames(vcov(fit)))`.
* **4.2**: `flat_par_display()` (`R/confint.R:897`) labels a bare
  `theta_k` from `log_sd_theta_index()`, so the peak fit now reports
  `theta_3 (sd of pk: 1 | id (Intercept))` in both the warning and the
  printed block.

## Items recorded and not acted on

* **4.1** (worst case `n + 1` gradient evaluations, 0.28 s on a
  33-parameter fit) and **4.3** (the 1e-8 threshold is absolute). Both
  are INFORMATIONAL in the review and both are already argued in the
  `diagnose_flat()` roxygen: the check is gated on a covariance that has
  already failed, and a parameter with a Hessian column norm near 1e-8
  carries a standard error near 1e4 and is not usefully identified. A
  candidate cap would trade a measured bound for an unmeasured
  heuristic, so it is left alone.
* **2.3** (`R/allfit.R:122` also calls `unlist(fixef())`): confirmed
  name-blind -- it feeds a `vapply()` reduced to a scalar spread, so no
  name reaches the user. Nothing to change; recorded so the next lane
  does not rediscover it.

## Punch-round counts

Core suite, one file per process: 131 files on disk, 131 RESULT lines,
131 unique names, none missing, none twice. **0 failed, 5 errored, 133
skipped, 7426 passed.** The 5 errors are the same three pre-existing
pkgload files (`test-compat-register.R` 1, `test-influence-plot.R` 2,
`test-structure.R` 2), unchanged from the pre-punch run and from main.

Passing count 7395 before the punch round (the 7393 reported earlier was
computed before the `test-nl.R` fixture correction landed, which is the
two-assertion gap the review noted; the reviewer's 7395 was right) and
7426 after: +31, from `test-diagnostics-ux.R` 104 to 116,
`test-methods.R` 44 to 59, and `test-nl.R` 40 to 44.

### Gated tier, re-run after the punch changes

`diagnose()`'s printed block and `ranef()` are both read by this tier,
so it was run again rather than assumed unaffected. Same environment as
before (private `ap-stan-cache`, `R_MAKEVARS_USER` on a copy of the
makevars).

| file | failed | error | passed |
| --- | --- | --- | --- |
| `test-brms-methods.R` | 0 | 0 | 950 |
| `test-brms-likelihood.R` | 0 | 0 | 394 |
| `test-brms-agreement.R` | 0 | 0 | 187 |
| `test-brms-priors.R` | 0 | 0 | 103 |

**1634 passed, 0 failed, 0 errored**, unchanged.

### Roxygen, re-checked

Idempotent after the punch edits: a further `roxygenise()` reproduces
`NAMESPACE`, `man/diagnose.Rd`, `man/fixef.Rd`, `man/frm_bootstrap.Rd`
and `man/ranef.Rd` byte for byte (md5 compared before and after).

### R CMD check --as-cran, re-run on the punched tree

Same invocation: no `--library=`, no `--no-manual`, TinyTeX and pandoc
on PATH. **58 checks. Status: 1 NOTE. No ERROR, no WARNING.**

```
* checking HTML version of manual ... [17s] NOTE
Skipping checking math rendering: package 'V8' unavailable
```

the one expected environment NOTE, cause measured
(`requireNamespace("V8")` is `FALSE`). `frmtmb-manual.pdf` built at
543,474 bytes. The phases this round touches are all OK:

```
* checking S3 generic/method consistency ... OK
* checking Rd files ... OK
* checking for code/documentation mismatches ... OK
* checking examples ... [20s] OK
* checking examples with --run-donttest ... [21s] OK
* checking tests ... [112s] OK
* checking re-building of vignette outputs ... [145s] OK
* checking PDF version of manual ... OK
```

### Tree after the punch round

22 modified files, unchanged in membership from the pre-punch list.
Untracked: `dev/api-findings.md` and `dev/reviews/2026-09-06-api.md`
(the review itself). No `dev/*.log`, no scratch file in the worktree,
`DESCRIPTION` untouched, the main checkout untouched, nothing committed,
and `wt-api` still at `028a862` as instructed.
