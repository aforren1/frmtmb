# Lane sampfront: frmtmb.sample readers and front end (round of 2026-11-11)

Base e8abed0a, frmtmb 0.81.0, frmtmb.sample 0.28.0. Brief:
`C:/Users/adf44/source/r/round-briefs/r082-sampfront.md`. Status: fix
round of review `dev/reviews/2026-11-11-sampfront.md` (FIX) done,
uncommitted; see "Fix round" at the end, which supersedes the earlier
sections where they differ. This file is the record; every number below
is pasted from the log named beside it (`dev/sampfront-out/`, local).

Summary: all twelve backlog items of the brief are closed. The fix
round added one core function, `integrated_design_zero()`, for the
structural 1530 rule, so frmtmb.sample needs the frmtmb release that
carries it (the floor moves from `frmtmb (>= 0.81.0)` to that release;
DESCRIPTION is left at 0.81.0 for the consolidation to set). Bumps:
frmtmb minor (a new export), frmtmb.sample minor (new exported functions
`pp_average()`, `model_weights()`, `prior_samples()`, a new argument of
`frm_sample()`, and a BREAKING acceptance of `newdata`).

## Arms and scripts

- Before: `rellib-r20` (frmtmb 0.81.0, frmtmb.sample 0.28.0). After:
  the lane library `wt-sampfront-lib` (frmtmb.sample only; core is
  unchanged) in front of `rellib-r20`.
- `dev/sampfront-repro.R [base]`: one reproducer per item; data seed 11
  (items 1 to 5, 7), seed 77 (item 6), sampler seed 1, `lap_pair()`
  seed 1. Logs `repro-base.txt`, `repro-lane1.txt`, `repro-lane2.txt`.
- `dev/sampfront-brms-probe.R`: brms 2.23.0 on the parity items, data
  seed 11, brms seed 1, 1 chain, 400 iterations. Log `brms-probe.txt`.
- `dev/sampfront-frame-probe.R`: what a fit's frame holds for the
  subset blockers (`frame-probe.txt`). `dev/sampfront-dul-probe.R`:
  what `drop_unused_levels = FALSE` leaves in a fit.
- `dev/sampfront-probe-attack.R [base]`: the sampfix reviewer's attack
  on the laplace probe (`dev/sampfix-rev-01-probe.R`) run on each arm.
- `dev/sampfront-run-test.R`, `dev/sampfront-runjobs.sh`: one test file
  per R process, at most 3 at once, lane or base arm.

## Item 1, backlog 4717: a subset of the fitted rows as `newdata`

brms 2.23.0 (`brms-probe.txt`, y ~ x + (1 | g), 48 rows):
`log_lik(newdata = d[1:10, ])` against the fitted columns 1:10, max
difference 0; rows `c(5, 3, 3, 20)`, 0; `loo(newdata = d[1:10, ])`
elpd -15.71843 against the sum of the fitted pointwise rows 1:10,
-15.71843, pointwise max difference 0; `waic` pointwise 0. With x moved
in row 2 brms recomputes that row (0.9335582) and leaves the others (0):
that is new data, which stays refused here by the policy of 0.81.0.
On `ar(time = t, gr = g)`: rows with `t >= 3`, max difference 0.4603141
from the fitted columns; rows 1:10 (whole leading runs), 0. So for an
autocorrelation term brms's subset answer is not the fitted column.

Change (`R/methods-draws.R`, `R/loo.R`, `R/loo-predict.R`):
`draws_newdata_rows()` places a subset by its row names (a repeated row,
`"3.1"` from `make.unique()`, by its stem), then requires every variable
the model reads to hold its fitted value at the placed row
(`draws_newdata_matches()`, the fitted-rows check of 0.81.0 cut to the
rows). A wrong placement cannot pass the comparison, and on the models
it serves a row equal in every variable has the same column.
`draws_subset_blocker()` refuses, naming the term, on an autocorrelation
term, a `me()` term, an `mi()` term, `subset()` rows, or a family whose
likelihood has no column per row. `log_lik()` computes the matrix and
cuts the columns; `loo()`, `waic()`, `psis()` and `loo_model_weights()`
(and so `model_weights()` and `posterior_average()`) hand the subset to
`log_lik()`. Every matrix still goes through `loo_check_finite()`, and
the reader is still `row_lpdf(exact = TRUE)` (`draws_row_loglik()`,
unchanged). `loo()`'s `yhash` is the subset's response.

Evidence: `repro-base.txt` refuses all five calls; `repro-lane1.txt`:
`log_lik(newdata = d[1:10, ])` dim 100 x 10, identical to the columns;
`d[c(5, 3, 20), ]` identical; `loo()` pointwise max difference 0 over
10 rows; `waic` 10 rows; `psis` 100 x 10. New test file
`test-loo-newdata-subset.R`: base pass=13 fail=11 err=5, lane pass=43
fail=0 (`t-subset-base.txt`, `t-subset-lane.txt`). The guard's absent
case is in it: a model without the term takes the same subset.

`pp_check()`'s `loo_*` types hand `newdata` to `log_lik()` and draw
their predictions on `newdata`, so they take a subset too, with no code
of their own. The plain tier showed it: `test-ppcheck-loo.R` asserted a
refusal of `dd[-1, ]`, which is a subset, and failed=1 on the lane. The
test now asserts a refusal of moved rows, and that on a subset the
weights are those rows' PSIS weights and the response is the subset's
(`run2/`: lane pass=35; base err=1). Tests that pinned the subset
refusal were rewritten to use moved rows: `test-loo.R`,
`test-loo-newdata-fitted.R`, `test-loo-predict.R`,
`test-loglik-autocor.R` (its no-term control).

Not done: the `loo_*()` functions still refuse a subset (filed). The
brief scoped `loo()` and `log_lik()`.

## Item 2, backlog 4584: `log_lik(combine = FALSE)` without rescor

brms 2.23.0 (`brms-probe.txt`): array 200 x 48 x 2, dimnames
`NULL, NULL, c("y", "z")`, sum over the layers equals the combined
matrix (0), layer y equals `resp = "y"` (0); one response gives the
matrix; `loo(combine = FALSE)` errors in brms. Change: `log_lik()` builds
each response's matrix (`draws_loglik_matrix()`) and stacks them with
`draws_resp_stack()`, brms's names (`ya` for `y_a`). `subset()` models
keep their refusal. Test `test-loglik-combine.R`: base err=2, lane
pass=13 (`run1/`).

## Item 3, backlogs 3215 and 4223: GEV from a formula

`gev_xi_chart_start()` (core) reads `mu` and `sigma` through
`frm_linpred()`, which refuses the formula route's object by its class
`frmtmb_unfitted`. On that route `fit$estimates` are the starting values,
which is what the backlog's fix names. Change (sample only):
`mo_sample_pars()` hands core a copy without the marker, as
`draws_base_fit()` does. Base: "frm_linpred() needs a fitted model";
lane: samples (`repro-lane1.txt`: b_Intercept 0.181, sigma 1.164, xi
0.086 on 100 rows simulated at 0, 1, 0.1, seed 4). Test in
`test-gev-tmp-xi.R`: base err=1, lane pass=16.

## Item 4, backlog 3594: `posterior_linpred()` and an unknown argument

`frm_check_dots(..., .allow = draws_new_level_args)`, as
`posterior_epred()`. Base returned the 100 x 48 matrix for `junk = 1`;
lane refuses. Test in `test-draws-spellings.R`: base fail=1, lane pass.

## Item 5, backlog 4370: `prior_samples()`

brms 2.23.0: warning "'prior_samples' is deprecated. Please use
'prior_draws' instead.", result identical to `prior_draws()`. Defined
as a plain function, as brms's is, that warns and calls `prior_draws()`
(the owner's generic while brms is loaded). Test in
`test-prior-draws.R`: base err=1, lane pass.

## Item 6, backlogs 1456 and 1530: the laplace probe and watch

1456: the watch took the first draw's pattern of missing cells as the
reference. Repro (`repro-base.txt`, the sampfix construction, c1 = -1
and k = 0 on draw 1): laplace `posterior_epred()` 480 of 480 cells
non-finite, full draws 80 of 480. 1530: `(0 + x | g)` at `x = 0` was
refused, because `NA * 0` is `NA`.

Change: a quantity reads the integrated values unless every cell that
is finite at fill `NA` is the same at fill 0, and every cell that is
`NA` at fill `NA` is the same at fills 0 and 0.5772156649015329
(`draws_fill_cells()`). The probe uses it on the first draw; the watch
now evaluates each draw itself (`watch(at_k, part)`), and a draw with a
missing cell is evaluated at the two finite fills, on the same random
numbers, and compared with itself. A cell that is the same at both
finite fills takes the value at 0. A draw with no missing cell costs no
extra evaluation (by construction); a draw with one costs two. Every
caller (`posterior_epred()`, `posterior_linpred()`, `posterior_predict()`
and its joint rescor path, `hypothesis()`, `pp_mixture()`,
`conditional_effects()`) passes its evaluator.

Evidence: `repro-lane1.txt`: 1456 refused; 1530 identical to the full
draws; the control at `x = 1` refused. The sampfix attack
(`attack-base.txt`, `attack-lane.txt`): base PASS 24 FAIL 3, lane PASS
25 FAIL 2; the one change is "epred newdata x=0 default re" FAIL to
PASS. The two remaining FAILs are the same on both arms (a new grouping
level on draws is refused by the draws-wide new-level refusal, before
the probe). The script stops at its cumulative section on both arms
(`lap_pair()`'s label count no longer matches on 0.81.0); the sections
before it ran. Tests in `test-laplace-draws.R`: base fail=3 err=1, lane
pass=107.

## Item 7

- 1547: each accessor (`as.array()`, `as.matrix()`, `as.data.frame()`)
  checks its own dots against `draws_accessor_dots`. Base: "draws_
  accessor_args() has no argument `chains`"; lane: "as.array() has no
  argument `chains`. Did you mean `chain`?". Test in
  `test-draws-methods.R`: base err=1, lane pass.
- 1655: `frm_sample(drop_unused_levels =)` on the formula route; on a
  fit it is refused by name with the other assembly arguments, and
  `update()` refuses it too. Base: "passing unknown arguments". Note:
  `frm(drop_unused_levels = FALSE)` keeps the level in the frame, but
  the grouping factor's block holds the 6 observed levels of 7
  (`dev/sampfront-dul-probe.R`, core behavior, filed, brms not
  measured). The test asserts the frame. Test in
  `test-stan-control.R`: base err=2, lane pass=51 after the assertion
  was corrected (the first lane run asserted 7 block levels and
  failed=2; that assertion was wrong, not the change).
- 1559: `posterior_average()` refusals quote names (`pa_quote()`, brms's
  `collapse_comma()`); brms 2.23.0 gave "Parameters 'b_x', 'nope' cannot
  be found in all of the models." `model_weights()` and `pp_average()`
  defined, exported, and added to the owner table (owner brms) and to
  `test-generic-collision.R`'s tables. `pp_average()` takes each model's
  draws as brms's `validate_draw_ids()` does (a random sample in drawn
  order), refuses `draw_ids` and a model of another response in brms's
  words, and summarizes with `posterior_summary()`. Tests in
  `test-postfit-draws.R`: base err=3, lane pass=89;
  `test-generic-collision.R`: base fail=43, lane pass=94.
- 1533: a message when `laplace = TRUE` samples a REML model with
  priors. Test in `test-laplace-draws.R` (with the ML control, which
  gets no message).

## Decisions and deviations

- `model_weights(weights = "loo"/"waic")` does not suppress the loo
  warnings as brms's `SW()` does; `posterior_average()` did not before
  either.
- brms's "it not supported" typo in `pp_average(method =)` is written
  "is not supported".

- Not run under OpenBLAS: no numerical code changed. The new
  exactness assertions compare a matrix with columns of a matrix that
  the same code computed in the same process (`log_lik(newdata = sub)`
  against `log_lik()[, rows]`, PSIS per column), or laplace draws with
  full draws, as the existing laplace tests do.
- `R CMD check` ran once, on the final source, as the lane-env asks.

## Test runs

Every log's `lib:` line shows frmtmb.sample from `wt-sampfront-lib` and
frmtmb 0.81.0 and the other extensions from `rellib-r20` (core is not
changed). Totals summed from the RESULT lines of the per-file logs.

- Plain tier (`plain/`, NOT_CRAN=true, one file per process, the build
  before the last refactor): frmtmb.sample 89 files pass=4150 fail=4
  err=1 skip=18 warn=0; and every file of core and the extensions that
  names frmtmb.sample: frmtmb 29 files pass=2934 fail=0 err=0 skip=13;
  frmtmb.eam 3 files pass=130 skip=2; frmtmb.latent 1 file pass=13;
  frmtmb.ode 1 file pass=15; frmtmb.spatial 2 files pass=238 skip=4.
  The four failures and the error were in three frmtmb.sample files,
  all from this change: `test-ppcheck-loo.R` and `test-loglik-autocor.R`
  pinned the refusal of a subset (rewritten, see item 1), and
  `test-message-uniqueness.R` found two refusal templates that
  `pa_split()` and `pp_average()` repeated (moved into
  `pa_refuse_unnamed()` and `pa_check_ndraws()`). After the fix
  (`run2/`): those three pass=110, 6, 35 with fail=0 err=0, and the six
  other touched files pass as well.
- Gated tier (`gated/`, FRMTMB_BRMS_FIT_TESTS=true, final build): every
  frmtmb.sample file, 89, pass=4521 fail=0 err=0 skip=1 warn=0 (the one
  skip is the scale tier's own gate, FRMTMB_SCALE_TESTS); and the seven
  files of core and the extensions that skipped on the gate in the
  plain tier: frmtmb 5 files pass=214, frmtmb.eam 1 file pass=13,
  frmtmb.spatial 1 file pass=126, all fail=0 err=0 skip=0.
- Sessions (`dev/release/sampfront-sess.log`, final build):
  `SESSION frmtmb.sample forward files=84/89 pass=4157 fail=0 err=0
  skip=18 warn=0`, the same line for reverse and preload. Every LEAK
  line carries `[with load of ...]`; none is untagged.
- Before and after per touched file (`run1/`, `run2/`, base arm
  rellib-r20): test-loo-newdata-subset base pass=13 fail=11 err=5, lane
  pass=43; test-loglik-combine base err=2, lane pass=13;
  test-laplace-draws base fail=3 err=1, lane pass=107;
  test-draws-spellings base fail=1, lane pass=100; test-draws-methods
  base err=1, lane pass=177; test-prior-draws base err=1, lane pass=93;
  test-postfit-draws base err=3, lane pass=89; test-gev-tmp-xi base
  err=1, lane pass=16; test-stan-control base err=2, lane pass=51;
  test-loo-newdata-fitted base err=1, lane pass=77; test-generic-
  collision base fail=43, lane pass=94; test-ppcheck-loo base err=1,
  lane pass=35; test-loo and test-loo-predict pass on both arms (their
  edits replace a subset by moved rows, which both arms refuse).
- `R CMD check --as-cran` of frmtmb.sample (`dev/sampfront-check.sh`,
  tree deleted after): "Status: 1 NOTE", the expected "Skipping checking
  math rendering: package 'V8' unavailable"; tests `[ FAIL 0 | WARN 0 |
  SKIP 124 | PASS 1079 ]`.

## Files touched

frmtmb.sample: `R/conditional-effects-draws.R`, `R/draws-brms.R`,
`R/generic-owners.R`, `R/loo-predict.R`, `R/loo.R`, `R/methods-draws.R`,
`R/postfit-draws.R`, `R/prior-draws.R`, `R/sample.R`, `NAMESPACE`,
`NEWS.md`, `man/` (frm_sample, prior_draws, sample-log_lik, sample-loo,
sample-loo_predict, new model_weights and pp_average). Tests: new
`test-loo-newdata-subset.R`, `test-loglik-combine.R`; edited
`test-draws-methods.R`, `test-draws-spellings.R`,
`test-generic-collision.R`, `test-gev-tmp-xi.R`, `test-laplace-draws.R`,
`test-loglik-autocor.R`, `test-loo-newdata-fitted.R`,
`test-loo-predict.R`, `test-loo.R`, `test-postfit-draws.R`,
`test-ppcheck-loo.R`, `test-prior-draws.R`, `test-stan-control.R`.
`dev/test-backlog.md` (12 CLOSED lines, 2 filed items); `dev/sampfront-*`.
Shared areas: `sample.R` (priortab and names lanes also edit it) in
four small hunks: `mo_sample_pars()`, `sample_assemble()`, the
`frm_sample()` formals and inert list and REML message, and roxygen;
`methods-draws.R`, `loo.R` (names, ppdraw lanes).

Files that `sed -i` edited were written back with LF line ends where the
checkout has CRLF; git normalizes them (`core.autocrlf=true`), and the
diff carries only the content changes.

## Fix round (review `dev/reviews/2026-11-11-sampfront.md`, FIX)

The coordinator's six points, in order. Scripts of the review were rerun
on the fix build; their logs are in `dev/sampfront-out/fix/` (the
review's own logs in `dev/review-sampfront-out/` are left as they were).

### 1. BLOCKER: 1530 decided by the design, not by a second fill

The second scalar fill (0.5772) is gone. A cell that is `NA` at fill
`NA` may skip the integrated values only when every design entry that
multiplies one at its row is exactly 0, read per column of the design.
That needs the design at new rows, which only core builds
(`pred_design()`, internal), so core exports one function,
`integrated_design_zero(fit, newdata, resp, allow_new_levels)` (new file
`R/integrated-design.R`): `TRUE` for a row where every group-level
design entry (`re_parts$mm`) and every smooth and Gaussian-process
basis entry (`sm_parts$Xr`, and no kriged row) is exactly 0, in every
linear predictor of the response; `FALSE` for every row on a nonlinear
body, a `me()` or `mi()` term, a latent-residual block, an
autocorrelation term (whose conditional mean reads the fitted rows'
residuals), a design that fails to build or whose rows are not
newdata's. A missing entry is not a 0. Nothing in core calls it; grep
of `R/` and `extensions/*/R` finds only frmtmb.sample.

frmtmb.sample (`R/methods-draws.R`): `draws_fill_cells(a, b, skip)`
takes a cell that is `NA` at `NA` as reading nothing only when it is
`NA` at 0 too (its own `NA`) or its row's `skip` is `TRUE`; a finite
cell must be the same at 0. The probe and the watch both use it; the
watch evaluates a draw with a missing cell once more, at 0, and returns
that value. `draws_design_skip()` reads the skip at the call's newdata
for `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`
(and the joint rescor draw, the AND over responses) and
`conditional_effects()` (per grid, per cell, with or without a level
plan; a grid whose design fails skips nothing). Without newdata
(fitted rows), `hypothesis()` and `pp_mixture()` there is no skip, which
is 0.81.0's rule. `conditional_effects(method = "posterior_predict")`
now watches the expected values, as its probe does, not the simulated
ones. The watch holds back each evaluation's warnings and raises only
the returned one's (through `frm_warning()`, the condition kept): a
simulator at the `NA` fill warned "NAs produced" six times on an
answered call in the first run of the new test.

1456 still holds: the all-NaN first draw is `NA` at both fills (its
own), and draws 2 to 6 are `NA` at `NA` and finite at 0 with no skip
(a nonlinear body), so the call is refused (test passes).

Evidence (fix build):
- `dev/review-sampfront-probe.R` (`fix/probe-lane.txt`): every case of
  the review's table is refused: A epred, linpred, predict; A control;
  B sum-coded; C `||`; D and D2 `step()`; E `conditional_effects()`.
- `dev/review-sampfront-probe-real.R` (`fix/probe-real-lane.txt`): the
  laplace call at `x = -1` in group 1 is refused (review lane build:
  returned mean -0.028 against the full draws' 0.8774).
- `dev/review-sampfront-same.R` (`fix/same-compare.txt`, 28
  quantities, base rds of the review against the fix build): every
  full-draws quantity and every laplace quantity both arms answer is
  identical (max |diff| 0); the refusals match; the laplace 1530 epred
  and seeded predict at `x = 0` are answered by the fix build (base
  refuses) and are identical to the full draws' (`identical()` TRUE on
  both, from `fix/same-lane-fix.rds`).
- `dev/sampfront-probe-attack.R` (`fix/attack-lane.txt`): PASS 25 FAIL
  2, as before; the two FAILs are the draws-wide new-level refusal.
- Tests: `test-laplace-draws.R` gains "a group-level design that
  cancels or a step() still refuses" (cases A, B, C, D2, E of the
  table and `conditional_effects(method = "posterior_predict")` at
  `x = -1`, which must refuse, and the `(0 + x | g)` at `x = 0` epred
  and predict-mode curve, which must equal the full draws). On the
  review's lane build the refusals passed silently (the review's
  table); on the fix build 118 pass, 0 fail (`fix/t-laplace-lane.txt`,
  `run4/`). Core test `tests/testthat/test-integrated-design.R`: 10
  pass (`run3/`), with the cancelling, sigma-term, missing-entry,
  autocorrelation, nonlinear, smooth and unbuildable cases.

A core defect found on the way and filed: a prediction rebuilds a
group-level term's design under the session's contrasts, so case B
errors ("does not match the fitted model") if the contrasts option is
reset after the fit; the test keeps the option set for the prediction.

### 2. MAJOR: rows without the fitted row names are placed by value

`draws_newdata_rows()` places a row by its row name when that names a
fitted row holding the same values; otherwise each row is placed on a
fitted row equal in every variable the model reads (an exact key per
row: labels as text, numbers by `sprintf("%a")`, then the full
comparison of `draws_newdata_matches()` on the placed rows). A row that
equals no fitted row is refused with its position ("its row 2 equals no
fitted row in every variable the model reads"); a lacking variable is
named. The blockers (autocorrelation, `me()`, `mi()`, `subset()` rows,
no column per row) are unchanged.

Evidence: `dev/review-sampfront-subset.R` and `-subset2.R` (51 cases,
`fix/subset-lane.txt`, `fix/subset2-lane.txt`): 30 answered, 21
refused; every answered case with a meant row has max difference 0
from the fitted columns of that row. Five outcomes changed from the
review's lane run, all from refused to answered: `tibble[11:20, ]`
(difference 0) and the four "names reset" calls (`log_lik`, `loo`,
`psis`, `loo_model_weights` of `d[11:20, ]` renamed 1:10), which the
script expected refused under the old rule and has no meant row for;
the test "rows without the fitted row names are placed by their
values" asserts `log_lik()` of that construction is identical to
columns 11:20, and likewise for `which(g == "2")` and a repeat. Every
refusal now gives a true reason. `test-loo-newdata-subset.R`: 53 pass
(`run3/`).

### 3. minor: the `me()` term as written

The refusal reads the term's `written` label: "the model has the term
me(x, sdx)" (`fix/subset-lane.txt`); the test asserts `me(x, sx)`.

### 4. minor: drop_unused_levels

The filed backlog item now carries brms 2.23.0's numbers from the
review (`dev/review-sampfront-dul.R`): N_1 = 7 levels and K = 3 columns
in brms, 6 effects and 2 coefficients in `frm()`. `?frm_sample` and
NEWS now say what the argument does today (the frame only) and that
brms keeps the level in both.

### 5. minor: why sessions count 84 of 89 files

The five `test-brms-suite-*.R` files (brmsformula, brmsterms, families,
helper-copy, methods) call `skip_unless_brms_suite()` at top level,
which skips the whole file without `FRMTMB_BRMS_FIT_TESTS=true`, and
testthat records no result for it, so an ungated session counts 84.
They pass in the per-file gated tier, and a gated forward session was
run (below).

### 6. minor: the subset's cost

Filed: `log_lik()` on a subset computes every fitted column, then cuts.

### Runs on the fix build

Every log's `lib:` line names frmtmb and frmtmb.sample from
`wt-sampfront-lib` and the other extensions from `rellib-r20`. Totals
are summed from the RESULT lines of the per-file logs.

- Plain tier (`fixplain/`, one file per process, 3 at once): frmtmb.sample
  89 files pass=4175 fail=1 err=0 skip=18 warn=0; frmtmb 37 files (the
  29 that name frmtmb.sample, the new `test-integrated-design.R`, and
  the files that read the namespace or the exports: core-boundary,
  arg-refusal, generic-collision, bracket-access, registry-helper,
  api-spellings, no-superassign, adefects, dpar-eta-api, guardrail1,
  nl-rtmb-scope, negbinomial2) pass=3415 fail=0 skip=13; frmtmb.eam 3
  pass=130 skip=2; frmtmb.latent 1 pass=13; frmtmb.ode 1 pass=15;
  frmtmb.spatial 2 pass=238 skip=4. The one failure was
  `test-conditions-census.R`: the watch's re-raise was a bare
  `warning()`; it is `frm_warning()` now, and the file passes (11 pass,
  `run4/`, final build), with `test-laplace-draws.R` (118) and
  `test-message-uniqueness.R` (6).
- Gated tier (`fixgated/`, FRMTMB_BRMS_FIT_TESTS=true, final build):
  frmtmb.sample 89 files pass=4540 fail=0 err=0 skip=1 (the scale
  tier's own gate) warn=0; frmtmb 5 files pass=214, frmtmb.eam 1 pass=13,
  frmtmb.spatial 1 pass=126, all skip=0.
- Sessions (`dev/release/sampfront-fixsess.log`, final build): forward,
  reverse and preload each `files=84/89 pass=4176 fail=0 err=0 skip=18
  warn=0` (the five top-level-skipped brms-suite files, point 5). Gated
  forward session (`fix-sess-gated-forward.txt`): `SESSION
  frmtmb.sample forward files=89/89 pass=4540 fail=0 err=0 skip=1
  warn=0`. No LEAK line is untagged in any of the four.
- `R CMD check --as-cran` (`dev/sampfront-check.sh <package>`, trees
  deleted after): frmtmb "Status: 1 NOTE", tests `[ FAIL 0 | WARN 0 |
  SKIP 506 | PASS 19152 ]`; frmtmb.sample "Status: 1 NOTE", tests
  `[ FAIL 0 | WARN 0 | SKIP 124 | PASS 1100 ]`. Each NOTE is the
  expected "Skipping checking math rendering: package 'V8' unavailable".
- Not run: the other extension suites in full. Core's change is one
  function no core path and no extension other than frmtmb.sample calls
  (grep of `R/` and `extensions/*/R`), so it reaches no other suite; the
  files of the other extensions that name frmtmb.sample ran above, and
  core's check ran its whole suite. OpenBLAS not run: no numerical code
  changed.

### Files touched in the fix round

Core: `R/integrated-design.R` (new), `NAMESPACE`,
`man/integrated_design_zero.Rd` (new), `NEWS.md`,
`tests/testthat/test-integrated-design.R` (new). frmtmb.sample:
`R/methods-draws.R`, `R/conditional-effects-draws.R`, `R/sample.R`
(roxygen), `R/loo.R` (roxygen), `NEWS.md`, `man/frm_sample.Rd`,
`man/sample-log_lik.Rd`, `tests/testthat/test-laplace-draws.R`,
`test-loo-newdata-subset.R`, `test-loo-newdata-fitted.R`.
`dev/test-backlog.md` (the drop_unused_levels item rewritten; the
subset-cost and contrasts items filed), `dev/sampfront-runjobs.sh` (R
gets no stdin: it read the job list, which truncated two summaries;
the per-file logs were complete), `dev/sampfront-check.sh` (takes the
package), `dev/sampfront-warn-probe.R`.
