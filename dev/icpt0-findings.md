# Lane icpt0: brms's reserved `Intercept` and `bf(center = FALSE)`

Lane `icpt0`, 2026-09-25, worktree `/home/user/wt/icpt0`, branch
`lane/icpt0`. Base: frmtmb 0.63.0 (`91cbf355`). brms 2.23.0.

## What brms does

brms 2.23.0 source read from the installed namespace
(`has_rsv_intercept()`, `terms_fe()`, `validate_terms()`,
`get_model_matrix()`, `data_rsv_intercept()`, `stan_center_X()`,
`frame_fe()`, `prior_fe()`, `reorder_pars()`), and every case below
reproduced with `make_standata()` and `get_prior()`. Script
`dev/icpt0-brms-probe.R`, output `dev/icpt0-brms-probe.out`.

1. `Intercept` (or the deprecated `intercept`) is reserved only in a
   fixed-effect formula WITHOUT an intercept. brms then sets the terms'
   `intercept` attribute back to 1, builds `model.matrix()`, drops
   `(Intercept)`, and fills a data column `Intercept` with ones. So
   factors get TREATMENT contrasts (`0 + Intercept + f` gives
   `Intercept, fb, fc, fd`; `0 + f` alone gives `fa ... fd`), and the
   `Intercept` column sits where its term is written
   (`0 + x + Intercept` gives `x, Intercept`).
2. Such a predictor is not centered (`stan_center_X()` is FALSE), so
   the intercept is an element of the `b` vector: `get_prior()` lists
   class `b`, coef `Intercept`, and no class `Intercept`. A class
   `Intercept` prior is refused ("do not correspond to any model
   parameter").
3. `1 + Intercept + x` and `Intercept + x` are not reserved: brms looks
   for a data column and stops ("can neither be found in 'data'").
4. A data column `Intercept` (or `intercept`) that is not all ones is
   refused in any model that uses the reserved name ("Variable name
   'Intercept' is reserved in models without a population-level
   intercept"). In a model with an intercept, a covariate column named
   `Intercept` collides with the renamed `(Intercept)` and brms stops
   ("Internal renaming led to duplicated names").
5. Distributional parameters (`sigma ~ 0 + Intercept + z`, `hu ~ ...`),
   nonlinear parameters (`a ~ 0 + Intercept + x`, identical to
   `a ~ 1 + x`, because brms never centers a nonlinear parameter),
   mixtures and multivariate responses all follow the same rule, one
   formula at a time.
6. Ordinal families refuse it: "Cannot remove the intercept in an
   ordinal model." `disc ~ 0 + Intercept + z` is accepted.
7. `mo()`, `s()`, `cs()`: nothing special. They sit outside `Xc`, and
   with no centering there is nothing to interact with. `cs()` needs an
   ordinal family, which refuses `0 + Intercept`.
8. `(0 + Intercept | g)` works in brms only when a population-level
   formula of the model reserves the name (the data column then
   exists); otherwise it is "not found".
9. `bf(y ~ x, center = FALSE)` is the same mechanism: the intercept
   column is kept, not centered, and is class `b`. It applies to the
   location formula only; `lf(sigma ~ z, center = FALSE)` sets it for a
   parameter formula. On an ordinal family it leaves the thresholds as
   class `Intercept` but uncentered.
10. `reorder_pars()` reorders parameter BLOCKS and puts a block named
    `..._Intercept` first; a `b_Intercept` that is an element of `b`
    keeps its column's place. So `fixef()` of `0 + x + Intercept` is
    `x, Intercept`, and of `bf(y ~ x, sigma ~ 0 + Intercept + z)` is
    `Intercept, x, sigma_Intercept, sigma_z`.

## What changed

The design choice: `0 + Intercept` is REWRITTEN at parse time into the
formula with an intercept, and the linear predictor is marked
`center = FALSE`. `bf(center = FALSE)` and `lf(center = FALSE)` set the
same mark. Every numeric path (design, start values, autoscale, rank
check, prediction, emmeans, the tape) then sees exactly the model
`1 + x`, and only the prior layer and brms's coefficient order read the
mark. No `Intercept` column is ever added to `data` or `newdata`.

- `R/parse.R`: `strip_rsv_intercept()` and `rsv_intercept_fixed()`
  detect brms's rule 1 on the fixed part and rewrite it; the rewrite is
  verified by `terms()` (intercept 1, the remaining term labels
  identical and in order) and refuses anything else by name.
  `parse_linpred()` returns `center = FALSE` and `rsv_intercept` (the
  number of terms before `Intercept`). `parse_one_response()` refuses
  the ordinal case and applies `bform$center` and an `lf()` formula's
  `center` attribute.
- `R/bf.R`: `bf(center = )` and `lf(center = )`, validated with
  `check_flag()`, carried through `bf_update()` and the `mvbind()`
  path. Roxygen: new `@param`s and a section "brms's reserved
  `Intercept`".
- `R/frame.R`: `rsv_intercept_order()` moves the intercept column to
  brms's position (rule 1) using `assign`; `check_rsv_intercept_data()`
  refuses a non-ones `Intercept`/`intercept` column (rule 4); the
  not-a-column refusal for `Intercept` now says where the name is
  reserved (rule 3, 8); every linear predictor record carries
  `center`.
- `R/priors.R`: with `center = FALSE`, class `Intercept` picks nothing
  (and the "not found" refusal says to use class `b`, coef
  `Intercept`), class `b` picks the intercept too, the prior table lists
  `b`/`Intercept` instead of class `Intercept`, `dpar_has_predictor()`
  is TRUE (so `sigma ~ 0 + Intercept` is a predictor, not the scalar a
  class `sigma` prior names), and `ordinal_center_offset()` returns
  NULL (rule 9). `?set_prior` updated.
- `R/brms-shapes.R`: `brms_fixef_rows()` does not lead with an
  uncentered intercept (rule 10).
- `R/compat.R`: grammar feature `0 + Intercept`, with `* works`,
  `group:ordinal refused` and `prior works` rules.
- `extensions/frmtmb.sample/R/sample.R`: `default_priors_for()` gives
  an uncentered intercept no default, as brms leaves class `b` flat.
- `extensions/frmtmb.sample/R/evidence-ratio.R`: the flat-prior advice
  names class `b`, coef `Intercept` for it; and see the defect below.
- `vignettes/brms-migration.Rmd`: one row in "What ports directly".
- Port rows: `dev/brmsport-verdicts.tsv`, `dev/brmsport-ledger.tsv`
  and `tests/testthat/test-brms-suite-standata.R` (see "Ports").

## Validation

Script `dev/icpt0-validate.R`, output `dev/icpt0-validate.out`, seed 1,
n = 80. Every number here is copied from that output.

A. Design against `brms::make_standata()`: 12 formulas (`0 + Intercept
+ x`, `0 + x + Intercept`, `-1 + Intercept + x`, `x + Intercept - 1`,
`0 + Intercept + f`, `0 + Intercept + x * f`, `0 + x * f + Intercept`,
`0 + Intercept + x:f`, `0 + Intercept`, `X_sigma` of
`sigma ~ 0 + Intercept + z`, `X_a` of `a ~ 0 + Intercept + x`, and
`bf(y ~ x * f, center = FALSE)`). Column names identical in all 12,
`max|dX| = 0` in all 12.

B. `b`/`Intercept` rows of the prior table against `brms::get_prior()`:
identical for `0 + Intercept + x` (3 rows), `0 + Intercept + f` (5),
`sigma ~ 0 + Intercept + z` (6), `center = FALSE` (3),
`center = FALSE` + `lf(sigma ~ z, center = FALSE)` (6), and a `mvbf`
(identical).

C. ML. `0 + Intercept + x` and `1 + x` have one likelihood; this is an
IDENTITY, not a measurement, because the rewritten formula builds the
same design and the same tape. Numerical check at full precision:
`dlogLik 0` and bitwise-identical coefficients for gaussian `x`,
gaussian `x + f`, `sigma ~ 0 + Intercept + z`, poisson `x * f`, the
nonlinear `a ~ 0 + Intercept + x`, and `center = FALSE`. Written last
(`0 + x + Intercept`) the column moves, the optimizer path changes, and
`max|dcoef| = 4.44e-16`. Against the classical references:
`max|coef - lm()| = 1.69e-07` (coefficient scale 1.33); against
`glm()`, `max|coef| 3.24e-06`, `dlogLik -1.43e-11`.

D. MAP against an independent log posterior written from brms's Stan
code (section F: `lprior += normal_lpdf(b | 0, s)` over the whole `b`,
including `b_Intercept`):

Each line: model and prior; frmtmb (intercept, slope); reference;
max|d|/se.

- `0 + Intercept + x`, class b normal(0, 0.1): (0.302231, 0.870160);
  (0.302231, 0.870160); 1.4e-06.
- `center = FALSE`, class b normal(0, 0.1): (0.302231, 0.870160);
  (0.302231, 0.870160); 1.4e-06.
- `1 + x`, class b normal(0, 0.1), slope only: (2.918390, 0.221095);
  (2.918390, 0.221095); 1e-06.
- `0 + Intercept + x`, class b coef Intercept: (0.062234, 1.106836);
  (0.062234, 1.106836); 3.1e-07.
- `1 + x`, class Intercept, centered: (-1.961279, 0.700155);
  (-1.961290, 0.700159); 8.3e-06.

The first and third lines are the whole point: one prior, two models,
intercepts 0.30 and 2.92. A class `Intercept` prior on
`0 + Intercept + x` is refused.

E. `cumulative()` under `bf(center = FALSE)`: the class `Intercept`
density sits at the thresholds, not at `tau - mean(x) * b`. Reference:
the cumulative logit written out, plus the log-Jacobian of the
`(tau_1, log increments)` storage that `?set_prior` documents for the
thresholds. center = TRUE: `max|d| 3e-07`; center = FALSE:
`max|d| 3.6e-06`; the two fits differ (first threshold 2.27 against
-0.61), so the offset is really dropped. A first reference without
the Jacobian disagreed by 0.094 in both arms alike; adding it closed
both, which is how the documented Jacobian was confirmed rather than
assumed.

F. brms's generated Stan code (in the output): for `0 + Intercept + x`
and for `center = FALSE`, `vector[K] b` and `lprior +=
normal_lpdf(b | 0, 1)`, with no `Intercept` parameter and no
`means_X`.

## Downstream methods

Checked in `/tmp/lanes/icpt0/icpt0-play2.R`, `icpt0-play3.R` and pinned
in `tests/testthat/test-rsv-intercept.R`:

- `predict()`/`fitted(newdata = )`: no `Intercept` column needed;
  `fitted(newdata)` identical to the `1 + x` fit.
- `simulate()`, `summary()`, `vcov()`, `confint()`, `hypothesis(
  "Intercept = 0")`, `variables()`, `coef()`, `model.matrix()`: as for
  `1 + x` (the internal column is `(Intercept)`, so names do not
  change); `fixef()`/`summary()`/`vcov()` take brms's order (rule 10).
- `conditional_effects()`: effects `x`, `f`; `Intercept` is not a
  variable. `emmeans()`: identical to the `1 + x` fit.
- `prior_summary()`/`set_prior()`/`get_prior()`/`default_prior()`:
  brms's rows (B); `route = "sample"` defaults match
  `brms::default_prior()` row for row (sample test).
- `sigma()`: `sigma ~ 0 + Intercept` returns the constant sigma.
- `REML`, `sparse_x`, `autoscale`, mixtures, `mvbf()`/`mvbind()`,
  `mo()`, `s()`, offsets, `poly()`: fitted; offset and `poly()` fits
  bitwise equal to the `1 + x` fits.
- `frm_sample()`: not run (tmbstan chains are slow); its default
  priors and its flat-prior advice are pinned directly.

Speed: nothing was added to the objective. The tape is the one `1 + x`
builds, which C shows as a bitwise-equal log-likelihood; no timing was
needed.

## Refusals added

Each names the replacement or the reason:

- `0 + Intercept` on an ordinal family (brms refuses it too).
- `Intercept` inside another term (`Intercept:x`, `I(2 * Intercept)`).
  brms reads these as products with a column of ones; frmtmb refuses
  rather than guess ("Intercept:x is x").
- `0 + intercept` (brms's deprecated lower-case spelling, which brms
  still accepts with a warning).
- A non-ones `Intercept`/`intercept` data column in a model that
  reserves the name (brms refuses it too).
- Class `Intercept` on an uncentered intercept: the existing "Prior
  target not found" now says to use class `b`, coef `Intercept`.
- `1 + Intercept`, `(0 + Intercept | g)` with no such data column:
  the existing not-a-column refusal now says where the name is
  reserved and that a group-level intercept is `(1 | g)`.

## Pins seen failing on the base build

- `test-rsv-intercept.R`, block "a data column of ones named Intercept
  gives brms's model": on base (`FRMTMB_LIB=base`) with a data column
  `Intercept = 1`, `0 + Intercept + f` fitted R's no-intercept coding,
  "rank deficient; dropping column(s): fd", coefficients `Intercept,
  fa, fb, fc` (the `d` cell and contrasts against it, not brms's
  model), and `fitted(newdata)` stopped with "Variable 'Intercept'
  missing from newdata". Base: `fail=1 err=9` for the file.
- `frmtmb.sample` `test-evidence-ratio.R`, "the flat-prior advice for
  a centered intercept has no stray field": fails on base (defect
  below).
- `frmtmb.sample` `test-default-priors-brms.R`, the new block: with the
  lane core and the BASE `frmtmb.sample`, `default_priors_for()` put
  brms's student-t on the uncentered intercept (2 failures), which is
  the counterfactual for the `sample.R` change.
- `test-brms-suite-standata.R` on base: `fail=1` (row 976 does not
  hold); on the lane, 87 pass.

## Defect found and fixed

`er_slot_spelling()` (frmtmb.sample) wrote
`set_prior("normal(0, 1)", class = "Intercept",  = "")` for the
intercept of a univariate model: `paste0()` over zero fields still
pastes its constants. Pinned above.

## Decided not to do, and divergences left

- `(0 + Intercept | g)` in a model that reserves the name: brms accepts
  it (the data column exists). frmtmb refuses it and names `(1 | g)`,
  the same group-level term. Supporting it means adding a column to
  `data` and to every `newdata` path, for a spelling of `(1 | g)`.
- Lower-case `intercept`: refused, not warned. brms names that
  coefficient `b_intercept`; accepting it under frmtmb's `Intercept`
  name would change a name silently.
- A `newdata` column `Intercept` that is not all ones is ignored; brms
  refuses it. The prediction is right either way (the intercept is
  ones by definition).
- Pre-existing, not changed: a nonlinear parameter's prior-table row is
  coef `(Intercept)` where brms lists `Intercept`, and
  `fixef()` of a nonlinear model leads with every intercept where brms
  keeps each parameter's block together. An ordinal `yo ~ 0 + x` is
  silently fitted as `yo ~ x`; brms refuses it ("Cannot remove the
  intercept in an ordinal model").
- `dev/test-backlog.md` still lists "`0 + Intercept`" under "Open -
  medium"; not edited here because neighbouring entries belong to
  other lanes this round. The consolidator can close it.

## Ports

`dev/brmsport-gen.R` cannot run here: it needs
`dev/brms-suite/brms_2.23.0.tar.gz` (sha256-checked), which is not in
the tree, and neither CRAN nor a v2.23.0 git tag is reachable. So the
generated file was edited by hand to what the generator emits (the
reason strings rendered with the generator's own `emit_reason()`), and
the diff touches only these rows:

- `standata:976` cannot transfer -> pass.
- `standata:970`, `standata:974` stay cannot transfer, with the new
  reason (the lower-case spelling is refused by name).

`dev/brmsport-ledger.tsv` was updated by hand the same way (outcome,
reason, held, and the messages from a recording run). The count tables
in `dev/brmsport-findings.md` were NOT regenerated: for
`tests.standata.R` they move from pass 42 / cannot transfer 40 to
43 / 39, and the "cannot transfer | absent" total from 106 to 105.

## Tests run

One file per process, lane library first. RESULT lines are in the
final report.

## NEWS entry

* `0 + Intercept` is brms's reserved intercept: `y ~ 0 + Intercept + x`
  fits the model `y ~ 1 + x` with the intercept as an ordinary class
  `"b"` coefficient, not centered, so a class `"b"` prior reaches it
  and class `"Intercept"` is refused, as in brms. It works in the
  location, distributional and nonlinear-parameter formulas; factors
  take treatment contrasts; `newdata` needs no `Intercept` column.
  `bf(center = FALSE)` and `lf(center = FALSE)` are the same
  mechanism. Ordinal families refuse `0 + Intercept`, as brms does.
  The maximum likelihood fit is unchanged; a MAP fit matches brms's
  prior placement (`dev/icpt0-findings.md`).
* `frm_sample()`'s flat-prior advice no longer reads
  `class = "Intercept",  = ""` for a univariate intercept.
