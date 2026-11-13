# Review of lane sampfront (round of 2026-11-11)

Reviewer: adversarial review, not the lane's author. Tree:
`frmtmb-wt-sampfront` (branch `wt-r82-sampfront`, uncommitted, base
e8abed0a). Lane library `wt-sampfront-lib` (frmtmb.sample built
2026-11-11 13:40 UTC, after the last source edit at 13:31 UTC); base
`rellib-r20` (0.81.0). Every script below is in `dev/` of the lane tree
and every log in `dev/review-sampfront-out/` (local, gitignored). The
probe and subset scripts take `base` as their first argument for the
before arm; each script's header gives its call.

## Verdict: FIX

One BLOCKER: the new two-fill laplace probe and watch (backlogs 1456 and
1530) let a quantity that reads the integrated random effects through,
computed at random effects of 0, without a word. 0.81.0 refused every
one of these cases. The other items hold up under attack. The subset
route of item 4717 gave no wrong column in 51 constructed cases (25
answered, each with difference 0 from the meant columns; 26 refused),
but it refuses tibble and dplyr subsets that it could answer (MAJOR,
finding 2).

## Findings

### 1. BLOCKER: the second finite fill is one scalar, so a read that cancels or is flat passes as "reads nothing"

`draws_fill_cells()` (R/methods-draws.R) accepts a cell that is `NA` at
fill `NA` when it is the same at fill `0` and fill `0.5772...`. The fill
is ONE value written into EVERY integrated component
(`draws_fit_at()`: `est[[cp]][] <- fill`). Any cell whose value depends
on the integrated values only through a combination whose weights sum
to zero, or through a function flat on [0, 0.5772], is the same at both
fills, and the watch then returns its value at fill 0: the prediction
with every random effect set to 0. That is the 1530 case generalized
far past "the fill times 0":

- `(1 + x | g)` at `x = -1`: the group term is `c * 1 + c * (-1) = 0`
  at every constant fill `c`.
- `(1 + f | g)` with sum-coded `f`, at the level coded `-1`. This is
  the common psych design (sum contrasts, by-subject slopes).
- `(1 + x || g)` at `x = -1`.
- A breakpoint model, `a + bb * (x - om) * step(x - om)` with
  `om ~ 1 + (1 | g)`, at an `x` below the population breakpoint and
  above the group's: `step()` is 0 at both fills.

Reproducer `dev/review-sampfront-probe.R` (draws from the lane's own
`lap_pair()` construction, data seed 77, draw seed 1; logs
`probe-lane.txt`, `probe-base.txt`):

| case | base (0.81.0) | lane |
|---|---|---|
| A. `posterior_epred`, `(1 + x \| g)`, x = -1 | refused | returned, max abs diff from full draws 1.01166 |
| A. `posterior_linpred`, same | refused | returned, 1.01166 |
| A. `posterior_predict`, same, seeded | refused | returned, 1.01166 |
| A control, x = -0.5 | refused | refused |
| B. `(1 + f \| g)` sum-coded, f = b | refused | returned, 1.73599 |
| C. `(1 + x \|\| g)`, x = -1 | refused | returned, 1.01099 |
| D2. `step()` breakpoint, group 6 | refused | returned, 0.304048 |
| E. `conditional_effects(re_formula = NULL)`, x = -1 | refused | returned, 0.932268 |

On real sampler output (`dev/review-sampfront-probe-real.R`,
`probe-real-lane.txt` / `probe-real-base.txt`; data seed 77, sampler
seed 1, 1 chain, 600 iterations, `y ~ x + (1 + x | g)`,
`posterior_epred(newdata = data.frame(x = -1, g = 1))`):

```
true group-1 value at x = -1: 1.022638
full draws, group 1, x = -1: mean 0.8774
full draws, re_formula = NA, x = -1: mean 0.0184
laplace draws, group 1, x = -1: mean -0.028 (returned without a word)   # lane
laplace draws: REFUSED: posterior_epred() needs the random effects, ... # base
```

The lane's laplace answer is the population-level prediction, labeled as
the group's. The same hole is in the probe (`draws_probe_reads()`) and
the watch (`draws_laplace_watch()`), so all six callers carry it.

A related weakness: `conditional_effects(method = "posterior_predict")`
watches the simulated values (`part = identity`), while
`posterior_predict()` watches the distributional parameters
(`part = function(o) o$dp`). Discrete simulated values coincide at two
fills by chance, so that caller's check is weaker even after the fix
below. The probe of that caller already reads the expected values; the
watch should too.

Proposed fix. The NA-cell test has to show that the cell does not
depend on the integrated values, and two scalar fills cannot show that.
Either:

1. (preferred) evaluate at fill 0 and at two or more VECTOR fills whose
   entries differ per component, with both signs and a spread of
   magnitudes (for example a fixed quasi-random vector such as
   `4 * sin(seq_along(v) * 2.399963)` and its negation, generated
   without touching `.Random.seed`). A linear cancellation then needs
   the design to be orthogonal to an irrational vector, and a step or
   threshold is crossed by one of the signs. Accept a cell only when it
   is the same at every finite fill; or
2. decide 1530 structurally: a cell reads nothing integrated when every
   design entry that multiplies an integrated component at that row is
   exactly 0, and keep 0.81.0's refusal for every other `NA` cell.

Add cases A, B and D2 above as tests that must refuse, with the x = 0
case of 1530 as the guard's absent case, and make the
`conditional_effects()` watch check the expected values.

### 2. MAJOR: a subset without the fitted row names is refused, though the lane's own invariant makes it answerable

Rows are placed by row names. A tibble has none, and dplyr verbs reset
them (`dplyr::filter(d, g == "2")` has row names 1 to 8), so
`tibble[11:20, ]`, `dplyr::filter(d, g == "2")` and
`dplyr::slice(d, 9:16)` are refused (`subset-lane.txt`,
`readers.txt`), while `subset(d, g == "2")` and `d[d$g == "2", ]` of
the same rows are answered. The reason given is "its rows do not hold
the fitted values of the rows their names give", which is wrong for a
tibble: it gives no names. brms's `log_lik(newdata =)` evaluates the
newdata rows themselves and reads no row names (its subset answer equals
its fitted columns, `dev/sampfront-out/brms-probe.txt`); I did not run
brms on a tibble. lane-env-r082.md: "a refusal or warning that fires on
a correct model is worse than none".

The lane's own comment (`draws_newdata_rows()`) says "on the models it
serves a row equal in every variable has the same column". I confirmed
it: a refit with rows 41:43 copied from rows 1:3 and newdata named 1:3
gives columns with max difference 0 from rows 41:43's
(`subset-lane.txt`). So placement by VALUE is safe on exactly the models
that pass `draws_subset_blocker()`: any fitted row equal in every
variable the model reads has the same column, and the ambiguity the
0.81.0 note worried about only arose with `me()`, `mi()` and
autocorrelation, which stay refused.

Proposed fix: when name placement fails, place each newdata row at any
fitted row equal in every variable the model reads
(`draws_newdata_matches()` per row, or a keyed match on the compared
variables), and refuse only a row that equals no fitted row. Say in the
reason that a row was not found, not that a name was wrong. Add the
tibble and `dplyr::filter()` cases as tests that answer.

### 3. minor: the `me()` refusal names the term by its mangled label

`log_lik(newdata = d[1:10, ])` on `y ~ me(x, sdx)` says "the model has
the term mexsdx, whose latent value per row brms does not carry to
newdata" (`subset-lane.txt`). The brief asks for a reason that names
the term. Proposed fix: name it as the user wrote it, `me(x, sdx)`
(the frame keeps the call; core's me term carries its expression).

### 4. minor: `frm_sample(drop_unused_levels = FALSE)` is accepted but changes nothing in the sampled model; brms keeps the level

The lane's open question, answered with brms 2.23.0
(`dev/review-sampfront-dul.R`, `dul.txt`; `make_standata()`, no
compile; `g` with 7 levels, 6 used; `f` with 3 levels, 2 used):

```
brms drop_unused_levels = TRUE : N_1 (levels of g) = 6  max J_1 = 6
  y ~ f: K = 2  colnames X: Intercept fb
brms drop_unused_levels = FALSE : N_1 (levels of g) = 7  max J_1 = 6
  y ~ f: K = 3  colnames X: Intercept fb fc
frm drop_unused_levels = TRUE : nlevels(fit$data$g) = 6  length(b) = 6  ranef rows = 6
  y ~ f: coefficients: Intercept fb
frm drop_unused_levels = FALSE : nlevels(fit$data$g) = 7  length(b) = 6  ranef rows = 6
  y ~ f: coefficients: Intercept fb
```

brms keeps the unused level in the group-level block (7 levels, so the
draws hold `r_g[7,Intercept]`, a draw from the prior) and as an all-zero
fixed-effect column. `frm()` keeps it in the frame only. So the new
`frm_sample()` argument passes a value that has no effect on the model,
and the lane's test asserts only the frame. The defect is core's, not
this lane's, but the filed backlog item says "brms's block was not
measured"; it is measured now. Proposed fix: update the filed item with
these numbers (both the block and the fixed design differ), and do not
describe `drop_unused_levels` in `?frm_sample` as working "As in
frmtmb::frm()" without saying what it changes today.

### 5. minor: sessions ran 84 of 89 files because five gated files skip at file level

The five files with a FILE line and no RESULT line in each of the three
session logs (`dev/release/sampfront-sess-sessions/*.txt`) are
`test-brms-suite-brmsformula.R`, `-brmsterms.R`, `-families.R`,
`-helper-copy.R` and `-methods.R`. Each calls
`skip_unless_brms_suite()` at top level, which skips the whole file
without `FRMTMB_BRMS_FIT_TESTS=true`, and testthat records no result for
a file skipped there. The sessions ran ungated, so those five never ran
in one session. They hide no failure in the per-file gated tier (pass
16, 4, 84, 36, 61; fail 0 err 0; `dev/sampfront-out/gated/`). I ran the
gated forward session (`dev/release/run-session.R frmtmb.sample
forward`, NOT_CRAN=true, FRMTMB_BRMS_FIT_TESTS=true,
FRMTMB_RELLIB=wt-sampfront-lib; `sess-gated-forward.txt`): `SESSION
frmtmb.sample forward files=89/89 pass=4521 fail=0 err=0 skip=1
warn=0`, every LEAK line tagged `[with load of ...]`. So nothing is
hidden. Proposed fix: the lane's record should say why the count is
84, and a lane whose change reaches the brms-suite files should run one
gated session as well.

### 6. minor: a subset costs the whole matrix

`log_lik(newdata = subset)` computes every fitted column and then cuts
(`draws_loglik_matrix(... cols)`): `out[, cols]` after the per-draw
`draws_row_loglik()` over all rows. `loo(newdata = d[1:10, ])` on a
100 000-row fit costs the full model per draw. brms evaluates the
subset only. Not a correctness issue; worth a backlog line, since the
cut could move inside the per-draw call.

## Checks that passed

- **4717 placement and refusals** (`dev/review-sampfront-subset.R`,
  `subset-lane.txt`, `subset-base.txt`; `dev/review-sampfront-subset2.R`,
  `subset2-lane.txt`). On the lane, every accepted subset gave max
  difference 0 from the fitted columns of the rows meant: `d[1:10, ]`,
  reordered `d[c(20, 5, 3), ]`, duplicated `d[c(3, 3, 3), ]`, the frame's
  own rows, rows equal in value to other rows (renamed to them),
  `droplevels()`, reordered factor levels, `g` as character, an extra
  column, custom row names, rows after an `na.omit` drop, data that has a
  row named "31" from `rbind`, `weights()`, `offset()`, `cens(cc, y2)`,
  `trunc()`, `poly()`, `scale()`, `gp()`, `s()`, `sigma ~ z`, `mvbf`
  with and without `rescor`. Refused, as they should be: names reset
  (`d[11:20, ]` named 1:10), x moved by 1e-12, the response moved, a
  weight, offset, `y2`, `trials()` or a `sigma`- or nl-formula variable
  moved, a dropped NA row, `ar`, `ar(cov = TRUE)`, `arma`, `cosy`,
  `unstr`, `fcor`, frmtmb.spatial's `sar(lag)` and `sar(error)`,
  `me()`, and `subset()` multivariate models. `g` as integer is refused
  (a label against a number): conservative, not wrong. The base arm
  refuses all 49 calls. `loo()`, `psis()` and `loo_model_weights()`
  refuse the renamed rows with the same reason.
- **Readers** (`dev/review-sampfront-readers.R`, `readers.txt`, 40
  draws, newdata `d[1:10, ]`): `log_lik`, `logLik`, `loo`, `waic`,
  `psis`, `loo(x, x0, newdata =)`, `loo_model_weights`,
  `model_weights` (stacking and loo), `posterior_average(weights =
  "loo", control = list(newdata =))` and `pp_check(type =
  "loo_pit_overlay")` all take the subset. `loo_R2`, `loo_predict`,
  `loo_epred`, `loo_linpred` and `loo_predictive_interval` refuse it,
  as the lane records and files; `kfold()` refuses everything.
  `loo_compare()` takes no `newdata` (brms's neither). One consequence
  to know, not a defect: `add_criterion(x, "loo", newdata = s)` stores
  the 10-row criterion, a later `loo(xa)` returns it, and
  `loo_compare(xa, x0a)` stops in the loo package on the differing
  observation counts. brms 2.23.0 does the same (`use_stored` in
  `R/loo.R` of the brms tree under `dev/brms-suite/`).
- **Six watch callers on ordinary draws**
  (`dev/review-sampfront-same.R`, `same-compare.txt`): 28 quantities
  (epred, linpred, predict seeded, hypothesis, conditional_effects both
  methods, the 1530 newdata, pp_mixture, multivariate predict with and
  without rescor, and laplace draws at `re_formula = NA`) are identical
  between base and lane on full draws and on the laplace calls that
  both arms answer; the refusals match. The lane's 1530 answers on
  laplace draws equal the full draws exactly (epred and seeded predict).
- **4584 layout**: the lane's brms 2.23.0 probe (`dev/sampfront-out/
  brms-probe.txt`) shows `200 x 48 x 2`, dimnames `NULL, NULL,
  c("y", "z")`; the lane's array (`repro-lane2.txt`) is `100 x 48 x 2`
  with the same dimnames, and the layers sum to the combined matrix
  (difference 0). Same dim order and names.
- **3215 GEV** (`dev/review-sampfront-gev.R`, `gev.txt`; data seed 4,
  n = 150, sampler seed 1, 2 chains of 1000): the formula route and the
  fit route agree on `y ~ 1`, `y ~ x` and `y ~ x + (1 | g)`; the largest
  rough z of a mean difference is 0.88 (b_Intercept, `y ~ x`).
- **New exports and owners** (`dev/review-sampfront-owner.R`,
  `owner.txt`): with brms attached before or after frmtmb.sample,
  `model_weights` and `pp_average` resolve to brms's generics and a
  brmsfit reaches brms's method; `prior_samples()` warns once in brms's
  words and reaches brms's `prior_draws()`. The 43 base failures of
  `test-generic-collision.R` are the lane's new table entries checked
  against the base build, which has neither generic
  (`run*/frmtmb.sample--test-generic-collision-base.txt`: "`actual` is
  length 21, `expected` is length 23", and the child's
  "object 'model_weights.frmtmb_draws' not found").
- **Line endings**: `git ls-files --eol` shows `i/lf` for every touched
  file (working copies are a mix of `w/lf` and `w/crlf`);
  `git diff --stat` and `git diff --stat --ignore-cr-at-eol` both give
  "30 files changed, 1256 insertions(+), 239 deletions(-)". The diff is
  content only.
- **Plain tier on the final build** (`dev/review-sampfront-out/plain/`,
  one file per process, 2 at once, NOT_CRAN=true, lane library, the
  lane's own runner `dev/sampfront-runjobs.sh` with my job list
  `dev/review-sampfront-plain-jobs.txt`): 89 files, pass=4157 fail=0
  err=0 skip=18 warn=0; every `lib:` line names `wt-sampfront-lib`.
  The lane's plain tier was on the build before its last refactor
  (fail=4 err=1 there); this closes that gap. The 18 skips are the
  brms-fit gate.
- **Gated tier** totals regenerated from the lane's logs: 89 files,
  pass=4521 fail=0 err=0 skip=1, every `lib:` line the lane library.

## Not checked

- OpenBLAS: not run; no numerical code changed.
- `R CMD check`: not rerun; the lane's "Status: 1 NOTE" stands unverified.

## Short check of the fix round

Reviewer's scripts are `dev/reviews/sampfront-short-*.R`. Logs are in
`dev/reviews/sampfront-short-out/`. Lane library at the check: frmtmb
built 16:07 UTC and frmtmb.sample built 16:32 UTC, both after their
last source edits (10:54 and 11:32 local).

### Verdict: SHIP

Two conditions, neither in the lane's code: the consolidation must set
frmtmb.sample's floor to the core release that carries
`integrated_design_zero()`, and the item in S2 below must be filed.

### S1. BLOCKER 1 is fixed

- The review's table, rerun by the lane (`dev/sampfront-out/fix/probe-lane.txt`,
  `fix/probe-real-lane.txt`): every case refuses.
- New attack (`sampfront-short-design.R`; `design-lane.txt`,
  `design-base.txt`): `y ~ x + (0 + x | g) + T` at newdata `x = 0`.
  Each `T` is another integrated term, so the term must be seen.
  - Refused, correctly, for every `T`: `s(z)`, exact `gp(z)`,
    `gp(z, k = 6)`, `s(z, by = w)`, `(1 | mm(g1, g2))`, frmtmb.spatial
    `car(A, gr = loc)`, a `sigma ~ (1 | h)` term (the reader is
    conservative: a gaussian's epred does not read sigma),
    `zi ~ (1 | h)`, a shared `(1 | i | g)` id across mu and sigma, and
    `mvbf` with both responses.
  - Answered, each identical to the full draws: the control
    (`(0 + x | g)` alone), and `mvbf` with `resp = "y"`, where the other
    response's correlated term is not read.
  - `(1 | gr(h, by = f))`: my data were invalid for it (the fit
    stopped).
  - New levels: the draws-wide new-level refusal stops the call before
    the probe.
- I found no case where the reader says "zero" and a laplace answer
  differs from the full draws.
- By construction, the integrated parameters are `b` and `miss`
  (`R/fit.R`). `miss` and the latent and autocorrelation blocks in `b`
  return `FALSE`. The remaining blocks of `b` reach predictions through
  the `re_parts` and `sm_parts` that the reader inspects.

### S2. minor, not a regression: `conditional_effects(method = "posterior_predict")` gives an all-NA band where a predictive draw reads a dpar's group term

Watching the expected values is not sound for every family. A gaussian
or student epred does not read `sigma`, and a negbinomial epred does not
read `shape`. With `sigma ~ (1 | h)` or `shape ~ (1 | h)` and
`re_formula = NULL`, the probe and the watch both pass, and the draws are
simulated at the NA fill.

Reproducer `sampfront-short-ce.R` (`ce-lane.txt`, `ce-base.txt`): every
estimate and bound is NA, with 6 "NAs produced" warnings. The base arm
gives the identical result, so this is not a regression and not a
plausible wrong number. `posterior_predict()` on the same model refuses,
because it watches `dp` (`design-lane.txt`).

Zero-inflated, hurdle and mixture models are fine. Their epred reads
`zi`, `hu` and `theta`, so an integrated value in any of them is caught.
Censoring does not enter epred or the draw's parameters.

Proposed fix: watch the distributional parameters in that path, as
`posterior_predict()` does. If not fixed now, file it.

### S3. MAJOR 2 is fixed

Reproducer `sampfront-short-place.R` (`place-lane.txt`). Each answer
below has difference 0 from the meant fitted columns:

- duplicate fitted rows (41:43 copies of 1:3), value-placed either way:
  `match()` takes the first, and the columns are identical, so it does
  not matter;
- an integer `k` given as double, with and without names;
- `g` with another level set;
- a `Date` and a `POSIXct` read through `as.numeric()`;
- a row named 12 holding row 30's values: name placement fails the
  value check, so it falls back to value placement on row 30.

Refused with a true reason:

- an NA in `x` ("its row 2 equals no fitted row ...");
- the same clock time in another time zone.

The key's separators could in principle collide on labels holding
`\r` or `\t`. The full comparison by `draws_newdata_matches()` on the
placed rows catches that.

### S4. Minor items are done

- `me(x, sdx)` is named as written (`fix/subset-lane.txt`).
- The drop_unused_levels item carries brms's N_1 = 7 and K = 3.
- 84 of 89 is explained.
- A gated session ran 89 of 89 files: pass=4540 fail=0, no untagged
  LEAK line (`fix-sess-gated-forward.txt`).
- The subset-cost item and the contrasts item are filed.

### S5. The core export

`integrated_design_zero()` has `@keywords internal` and its NEWS entry
says "(internal, exported for frmtmb.sample)". That fits core's other
exports for extensions. The name says what it returns, and the contract
is conservative: `FALSE` wherever the design does not say. Its only
caller is frmtmb.sample.

frmtmb.sample still says `frmtmb (>= 0.81.0)` and calls the function
through `import(frmtmb)`. Against 0.81.0 it fails at run time on laplace
draws with newdata. The floor must move to the release that carries it
(core minor bump). The lane records this as left to the consolidation;
it is required, not optional.

### S6. Runs, regenerated from the logs

- Gated frmtmb.sample: 89 files, pass=4540 fail=0 err=0 skip=1. All 96
  `lib: frmtmb` lines name the lane library.
- Sessions, three orders: 84/89, pass=4176 fail=0.
- Gated forward session: 89/89, pass=4540.
- No untagged LEAK line.
- The new tests fail on base: `test-loo-newdata-subset.R` fail=13
  err=6, `test-laplace-draws.R` fail=3 err=2 (`run3/`). The lane build
  passes: pass=53 and pass=118.

Not needed:

- Other extension suites in full. The core change is one new export
  that no core path or other extension calls. The files of other
  extensions that name frmtmb.sample ran, and core's check ran its
  whole suite.
- OpenBLAS. No numerical code changed. The reader tests design entries
  for exact 0. A BLAS difference in a smooth basis product could only
  turn an exact 0 into a nonzero, which means refusing, never answering
  wrongly.
