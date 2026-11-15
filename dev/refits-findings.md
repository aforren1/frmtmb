# Lane refits: response coding, refit plumbing, core refusals

Round of 2026-11-13, base ec8178d3 (frmtmb 0.82.0), worktree
`frmtmb-wt-refits`, library `wt-refits-lib`. Brief:
`round-briefs/r083-refits.md`. Status: DONE, uncommitted.

Scripts (all in `dev/`): `refits-lib.R` (shared setup; `ARM=before`
reads rellib-r21 alone), `refits-repro.R` (one section per item, logs
`refits-repro-before.log` and `refits-repro-after.log`), `refits-brms.R`
and `refits-brms2.R` (brms 2.23.0, no compile; logs `refits-brms.log`,
`refits-brms2.log`), `refits-brms-me.R` (brms fit plus marginaleffects),
`refits-probe-*.R` (single questions). Logs are local (gitignored);
every number below is pasted from them.

## Summary

Closed: 1259, 1285, 1567, 1687, 2939, 3188, 3427, 3430, 3434, 3436,
4375, and 4018 in the part frmtmb can fix. Not done: 3225 (group-level
`cs()`, much larger than its size; re-filed with the sizing) and 4350
(upstream only; recorded as marginaleffects-1). Bump: minor (two
BREAKING behavior changes of `update()`; a refusal lifted). Floor:
frmtmb.latent's changed `test-hmm-starts.R` needs this frmtmb (its
`update()` takes `na.action` from the fit); no extension's code changed.

## Files touched

- `R/autocor.R`: `autocor_outer_sigma()` (4375).
- `R/bf.R`: `expand_dot_bform()`, `dot_drop_removed()` (1567).
- `R/compat.R`: the notes of two rules, `r("gr_prec", "*", ...)` and
  `r("gr_cov", "*", ...)`, gain one sentence on the `||` term (2939).
  No status changes.
- `R/conditional-effects.R`: `ce_base_frame()` (3436).
- `R/confint.R`: `confint()` check order (4375); `update()` takes
  `refit_args`, records `update_shown_formula()` (3427, 3430); the
  `?update` paragraph.
- `R/fit.R`: `frm()` stores `refit_args` (3430); `fit_assembled()` sets
  `marginaleffects_modeldata` (`me_modeldata()`, 3188).
- `R/frame.R`: `extract_y()` (1259); `gr_cov_check_within(split = )`,
  `gr_known_split()`, the split and the `col_pick` field (2939);
  `assemble_frame(smooth_map = )`, `smooth_map_of()`,
  `smooth_frozen_of()`, `smooth_random_frozen()`, `spec_k` (3434).
- `R/influence.R`: passes `smooth_map` (3434); doc.
- `R/interop.R`: `recover_data.frmtmb_fit()` (4018).
- `R/parse.R`: `expand_double_verts()` for `gr(cov/prec)` (2939).
- `R/predict.R`: `fitted_mv()` (1285), `col_pick` in the newdata
  design (2939). Not the contrasts or CE hunks.
- `R/sandwich.R`: `resolve_cluster()` (3430).
- `R/sugar.R`: `refit()`, `refit_bernoulli_codes()` (1687).
- `man/refit.Rd`, `man/update.frmtmb_fit.Rd`,
  `man/influence.frmtmb_fit.Rd` (roxygen).
- Tests: new `tests/testthat/test-refit-plumbing.R`; changed
  `test-bernoulli-coding.R`, `test-double-bar.R`, `test-formula-dot.R`
  and `extensions/frmtmb.latent/tests/testthat/test-hmm-starts.R`.
- `NEWS.md`, `dev/test-backlog.md` (CLOSED lines, new section),
  `dev/upstream-bugs.md` (emmeans-1 workaround text, marginaleffects-1
  and -2), this file, and the `dev/refits-*` scripts.

## Items

### 1285: `fitted(mv, newdata = )` on a multivariate ordinal fit

Misdiagnosed in part: it fails only for a ONE-row `newdata`, not on
every fit. `fitted_mv()` stacked `p[, , k]`, which drops a 1 x 4 layer
to a vector, and `nrow()` of a vector is `NULL`. Fixed with
`drop = FALSE` (`R/predict.R`). `refits-repro.R` section B, seed 77:

```
before: B cs 1 row       ERROR values must be length 1, but FUN(X[[1]])
                         result is length 0
before: B no cs 2 rows   2 4 6
after:  B cs 1 row       1 4 6
after:  B no cs 1 row    1 4 6
after:  B 1 row equals row 1 of 2 rows   1 0   (same names; max diff 0)
```

### 1259: character ordinal response with non-numeric labels

brms 2.23.0 refuses every character ordinal response, the digits
`"1".."4"` too, with "Family 'cumulative' requires either positive
integers or ordered factors as responses." (`refits-brms.log`, A).
frmtmb now refuses every character ordinal response in those words
plus the response's name and the fix (`extract_y()`, `R/frame.R`).
First round: digit text kept fitting; the fix round (coordinator, the
standing rule to match brms) refuses it too, BREAKING, and drops the
"character vector" coding from `?influence.frmtmb_fit` and from
`test-thres-refit.R`'s codings. Seed 2501, n = 60:

```
before: ERROR [simpleError]: missing value where TRUE/FALSE needed
        (after the warning "NAs introduced by coercion")
after:  ERROR [frmtmb_error]: Family 'cumulative' requires either
        positive integers or ordered factors as responses. 'y' is a
        character vector with values that are not numbers ('moderate',
        'severe', 'mild', 'none'), ...
both:   ordered factor over the same labels, and the integer codes:
        Intercept[1..3] -2.0394597 0.1030346 2.3889422, x 2.5289973
fix round, reviewer's probe2 (dev/refits-f2-probe2-after.log):
  'yd' "1".."4": first round fitted (Intercept[1] -1.171250 ...);
  now ERROR Family 'cumulative' requires either positive integers or
  ordered factors as responses. 'yd' is a character vector ('1', ...
  sratio() the same; "2.5" and "n/a" refused in the same sentence
```

### 1687: `refit(newresp = )` of a recoded bernoulli fit

`refit()` now takes the response's own two values as well as the 0/1
codes (`refit_bernoulli_codes()`, `R/sugar.R`).

The first round read a vector of only 0 and 1 as codes. The review
(F1) found that silent-wrong: on a fit coded 1/2 the number 1 was a
success in `refit(fit, rep(1, n))` and a failure in `refit(fit, c(2,
1, 1, ...))` and in `residuals(newdata = )`. The fix round applies the
reviewer's rule: values first; else 0/1 is codes; else refuse; a
vector both readings accept with different results (always constant)
is refused, naming how to pass it (a factor states values, a logical
states codes); `boot_refits()` passes its draws as logicals, so they
never reach the guess. Evidence (seed 9, n = 40, the reviewer's
probe5 data; `dev/refits-f1-readers.R`, log `refits-f1-readers.log`;
reviewer's scripts rerun, logs `dev/refits-f1-bern-after.log`,
`dev/refits-f1-probe5-after.log`):

```
refit(f12, c(2, 1, 1, ...)) = refit(f01, c(1, 0, 0, ...))   TRUE
refit(f12, rep(1, n))           refused: "reads two ways on this fit"
residuals(f12, newdata y = 1) = residuals(f01, newdata y = 0)  TRUE
pp_check(f12, newdata y = 1): the y plotted is 0
bayes_R2(draws, newdata) f12 (1s) = f01 (0s)                   TRUE
brms 2.23.0 standata(empty fit, newdata y = 1)$Y: 0
log_lik(draws, newdata = ): refused by frmtmb.sample (no newdata
  argument), so not a reader here
reviewer's bern script, against its run on the first-round build: the
  only changed lines are the three constant ambiguous vectors (1/2,
  -1/0 and character 1/2), now refused; {0,1}, logical, factor,
  {-1,1}, {0,2} identical; bootstrap 0/1 vs 1/2: identical, max|diff| 0
```

First round, seed 9, n = 100:

```
before: ERROR refit(): a bernoulli response is refitted on its 0/1
        codes ... newresp holds '-2', '-1'
after:  refit(fb, d7$yb): max |coef - fit coef| = 0
after:  factor "no"/"yes" response refitted on itself: 0
both:   refit on simulate() codes: Intercept 0.2562342, x 0.9650831
```

### 1567: `y ~ . - w` refit on data without `w`

`terms()` expands the dot to `y ~ (x + z + w) - w`, which names `w`.
brms 2.23.0 stores that same formula but reads only the terms that
stay: `standata(y ~ (x + z + w) - w)` on data without `w` gives columns
Intercept, x, z (`refits-brms2.log`). `expand_dot_bform()` now leaves a
removed column out of the dot (`dot_drop_removed()`, `R/bf.R`):
`y ~ . - w` becomes `y ~ (x + z)`. A removed column read elsewhere in
the formula (`. - w + w:x`) stays. Seed 3:

```
before: y ~ (x + z + w) - w; update(newdata = without w):
        ERROR The model uses `w`, which is not a column of `data`
after:  y ~ (x + z); update(newdata = without w): coefficients
        identical to update(newdata = with w); Intercept, x, z
```

Not done: `y ~ x + w - w` (no dot) is still refused on data without
`w`; brms accepts it (columns Intercept, x). Filed as a new LOW item.

### 3427: `getCall()` after a delta update

`update()` records the delta as a call (`update_shown_formula()`,
`R/confint.R`): the plain formula when the resolved `bf()` holds
nothing else and the family is the default or the call's `family`
argument, else `frmtmb::bf(...)` from `bform_call()`. Seed 7:

```
before: frm(formula = list(formula = z ~ x + (1 | g), pforms = ...
        3125 characters
after:  frm(formula = z ~ x + (1 | g), data = dk)   41 characters
after:  frm(formula = frmtmb::bf(z ~ x, sigma ~ x), data = dk)
both:   update() of the updated fit gives the same fixef (identical)
```

### 4375: `confint(parm = "x")` on a saved `ar(cov = TRUE)` fit

`confint()` refuses a fit saved by 0.76.0 only when the asked
parameters include `sigma`'s coefficients or the method is `"boot"`
(`autocor_outer_sigma()`, `R/autocor.R`). Neither the estimate nor the
inverse-Hessian block of another parameter changes when `sigma` alone
is reparameterized (the Jacobian rows of the others are unit rows).
Seed 4, the saved fit made by removing `sigma_is`:

```
before: confint(old, parm = "x"): ERROR ... built by frmtmb 0.76.0
after:  confint(old, parm = "x"): x 0.1900523 0.6936 0.4418262
        confint(fit, parm = "x"): x 0.1900523 0.6936 0.4418262
both:   confint(old): ERROR ... built by frmtmb 0.76.0
```

### 3430: `vcov_cluster()` outside the model; `update()` at the caller

- `vcov_cluster(cluster = ~ cl)` with `cl` not in the model reads the
  call's data. Those data are now checked to still hold the fit's rows
  and values (`refit_data_matches()`), and the cluster is aligned by
  row name, not by position (`resolve_cluster()`, `R/sandwich.R`).
  When the data are gone, the refusal says so. A changed value of `cl`
  alone cannot be detected: the fit never read `cl`.
- Fix round (review F2, silent-wrong, pre-existing): `~ cl` was
  evaluated in the fit's frame with the formula's environment as
  enclosure, so a workspace `cl` of length n replaced the data's
  column. Now the variables are read from the fit's own columns when
  it holds them all, then from the call's data (checked as above), and
  from the formula's environment only when the data lack them, as
  `sandwich::vcovCL()` reads them. Reviewer's probe2 V10 (seed 1, 120
  rows, a global `cl <- factor(rep(1:2, length.out = 120))`; logs
  `reviews/refits-review-probe2-after.log`, `refits-f2-probe2-after.log`;
  every other V row unchanged):

```
first round: 0.0162685 0.0645226 0.0215026
fix round:   0.3053286 0.1815976 0.1108744  (= dd$cl's)
```
- `update()` took `na.action`, `se` and `verbose` from the stored call
  and evaluated them where it was called. `frm()` now records the
  values on the fit (`refit_args`), and `update()` uses them
  (`R/fit.R`, `R/confint.R`). Seed 1:

```
before: update() of a fit made with na.action = my_na inside a function:
        ERROR object 'my_na' not found
before: update() of a fit made with se = s inside a function, a global
        s of length 100 present: ERROR `se` must be TRUE or FALSE, not
        an integer of length 100
after:  both refit; residual lengths 480 480 (na.exclude padding kept)
```

### 3188: marginaleffects `predictions()` on a matrix response

marginaleffects 1.0.0 puts the model data in a data.table, which
cannot hold a matrix column (upstream, `dev/upstream-bugs.md`
marginaleffects-2). A fit whose data hold a matrix column now carries
the data with that column split (`Y.a`, `Y.b`, `Y.c`) in the
`marginaleffects_modeldata` attribute, which marginaleffects reads
before insight (`me_modeldata()`, `R/fit.R`). insight's `get_data()`
and `get_response()` keep the matrix. Seed 5, n = 200:

```
before: predictions(fm): ERROR Column 3 ['Y'] is length 600 but column
        1 is length 200; malformed data.table.
after:  predictions(fm): 600 x 15; first row a 3.167451 se 0.1390591
both:   avg_predictions(): a 2.860000, b 4.105000, c 3.034999 (= the
        column means of fitted())
```

brms could not be checked: brms 2.23.0's `predictions()` needs the
`collapse` package, which the user library lacks
(`refits-brms-me.log`). The review checked it against
`nnet::multinom()` instead (max relative difference 4.3e-6).

Fix round (review LOW): `refit(newresp = )` of a matrix-response fit
now writes the new response into the frame's `data_frame`, so the
`marginaleffects_modeldata` of the refit holds it (`R/sugar.R`); it
held the old response. Test in `test-refit-plumbing.R`, seen failing on
the first-round build (`dev/refits-tests/plumbing-prefix.log`).

### 4018: `recover_data()` errors lose their class

Partly fixed. emmeans's `ref_grid()` calls the method inside `try()`
when `data` is missing, and stops with a returned text. No returned
value can carry the class: `tryCatch()` turns a condition whose type is
character into a `simpleError` (R 4.6.1 `tryCatchOne()`; tried and
reverted, `refits-probe-emm2.R`). When `data =` is given, `ref_grid()`
calls the method outside `try()`, and the method now raises the
condition as it is (`R/interop.R`):

```
before: ref_grid(fit, dpar = "nope", data = d): simpleError
after:  ref_grid(fit, dpar = "nope", data = d): frmtmb_error
both:   ref_grid(fit, dpar = "nope"): simpleError, the message itself
```

### 3434: `s(x)` in an `influence()` deletion refit

A deletion refit now keeps the fit's smooth basis: `assemble_frame()`
takes `smooth_map` and builds each smooth from the fitted one at the
subset's rows (`smooth_basis_at()`), as `predvar_map` keeps `poly()`
(`R/frame.R`, `R/influence.R`). A fit made before records which term a
smooth came from (`spec_k`) refits as before. Exact reference: the
deletion of row 1 written as the fit's basis at the 59 rows, the
wiggly columns as one `homdiag()` block and the null-space column as a
fixed one (seed 8, n = 60, `y ~ z + s(x, k = 5)`):

```
            table_fx    written_fx  table_z    written_z
before:   -0.1869716  -0.1858258  0.2915599  0.2915267
after:    -0.1858256  -0.1858258  0.2915267  0.2915267
max |dfbetas| (Intercept): before 0.357, after 0.263
```

### 3436: `ce_base_frame()` repeats `fit_refit_data()`

Done: it calls `fit_refit_data()` (`R/conditional-effects.R`).

### 4350: the marginaleffects argument warning (upstream)

Nothing frmtmb can change: `sanity_dots()` keys its white list on
`class(model)[1]` inside the function. Recorded as
marginaleffects-1 in `dev/upstream-bugs.md`. Not posted.

### 2939: a factor in a `||` term over `gr(cov = )` / `gr(prec = )`

brms 2.23.0 fits `(0 + f || gr(g, cov = A))` as one sd per column over
the known matrix, `r_1_k = sd_1[k] * (Lcov_1 * z_1[k])`
(`refits-brms.log`, E), and codes `(1 + f || gr(g, cov = A))` with the
formula's contrasts: Z columns Intercept, fb, fc; sd rows Intercept, fb,
fc (`refits-brms3.log`). frmtmb now fits both:

- `expand_double_verts()` (`R/parse.R`) writes a `gr(cov = )` or
  `gr(prec = )` term under `||` as `diag()` over the whole term, brms's
  coding, not as lme4's one-term pieces. `mm(cov = )` keeps lme4's
  pieces (filed: a factor there is still refused).
- The frame splits such a term into one single-column `gr_cov` or
  `gr_prec` component per column (`gr_known_split()`, `R/frame.R`). Each
  keeps the term's bar, with `col_pick` naming its column, which the
  newdata design reads (`R/predict.R`). No new covariance structure:
  every reader of a one-column `gr_cov` block reads these.

Found on the way, and fixed by the coding change: with lme4's pieces and
the factor refusal lifted, `(1 + f || gr(g, cov = A))` had fitted
Intercept, fa, fb and fc, one sd more than brms, logLik -115.108 where
brms's model gives -128.694 (`refits-2939b.R`, seed 11). Never shipped:
on 0.82.0 that term was refused.

Exact references, `refits-2939.R` and `refits-2939b.R`, seed 11, 12
levels, A[i, j] = 0.5^|i - j|, n = 120 (log `refits-2939-after.log`):

```
(0 + f || gr(g, cov = A)), (0 + fa + fb + fc || ...), cor = FALSE:
    -115.107965226781 (all three)
exact gaussian marginal ML (optim):   -115.107965226781
frmtmb - exact:                        -8.527e-14
gr(prec = solve(A)), factor and numeric: -115.107965226781 (both)
(1 + f || ...) and (1 + fb + fc || ...): -128.694123077274 (both)
    exact ML -128.694123077258; sds 0.694863 0.789940 1.277404
fitted(newdata), factor against numeric: max |diff| 0
before (0.82.0): ERROR ... is not available yet
```

### 3225: group-level `cs()`: not done, larger than M

brms 2.23.0 fits `(cs(1) | g)` on `sratio()` with a correlated block of
one effect per threshold added to `mucs[n, k]`, `cumulative()` with its
"experimental" warning, and `(1 + cs(x) | g)` as six correlated effects
(`refits-brms.log`, F). frmtmb has no random-effect design per
threshold: the objective, `fitted()`, `predict()` with new levels,
frmtmb.sample's `lp_Z()`, brms's names and the priors would each need
it. The refusal stays (it names the term and the alternative). Filed
with this sizing under "Filed in the round of 2026-11-13".

## Tests

New file `tests/testthat/test-refit-plumbing.R`, 12 tests, one per item
fixed (`dev/refits-test1.R`, logs `dev/refits-tests/plumbing-*.log`):

```
RESULT test-refit-plumbing.R (rellib-r21) pass=4 fail=6 error=10 skip=0
       warn=2 tests=12   (every test has a failure or an error)
RESULT test-refit-plumbing.R (lane)       pass=45 fail=0 error=0 skip=0
       warn=0 tests=12
```

Changed pinned tests (each pinned the behavior this lane changed):
`test-bernoulli-coding.R` (refit with the values was refused; now it
equals the refit with the codes, and a value that is neither is
refused), `test-double-bar.R` (the factor refusal over `gr(cov = )`
becomes the brms coding; the "named as written" refusal is now checked
on `mm(cov = )`), `test-formula-dot.R` (`y ~ . - x2` deparses to
`~(x1 + g)`, where brms writes `~(x1 + x2 + g) - x2`; the design
columns are still brms's).

Core, one file per process, 3 at a time (`dev/refits-par.sh core`,
log `dev/refits-core.log`, every file's `lib:` the lane library):

```
CORE ran 331 of 331; pass=24434 fail=1 err=0 skip=240 warn=0
the one failure: test-formula-dot.R, the pinned formula text above,
fixed in the test after the run
```

Extensions, every one installed from the worktree into the lane
library (`dev/refits-install-ext.sh`), one file per process
(`dev/refits-par.sh ext`, log `dev/refits-ext.log`, all 199 `lib:`
lines the lane library):

```
EXT ran 199 of 199
frmtmb.coupling pass=551  fail=0 err=0 skip=5  warn=0
frmtmb.eam      pass=1851 fail=0 err=0 skip=5  warn=0
frmtmb.latent   pass=696  fail=1 err=0 skip=2  warn=0
frmtmb.learn    pass=474  fail=0 err=0 skip=13 warn=0
frmtmb.ode      pass=549  fail=0 err=0 skip=1  warn=0
frmtmb.sample   pass=4238 fail=0 err=0 skip=18 warn=0
frmtmb.spatial  pass=529  fail=0 err=0 skip=6  warn=0
frmtmb.survival pass=575  fail=0 err=0 skip=0  warn=0
```

The one failure pinned the old `update()`: frmtmb.latent's
`test-hmm-starts.R` expected `hmm_starts()` to refuse a fit made with a
local `na.action`, which `update()` now takes from the fit. The test
now builds the refusal with an empty local `data2`, which the fit does
not store, and asserts that the local `na.action` refits.

Rerun of every changed test file (`dev/refits-rerun.log`):

```
test-formula-dot.R      pass=27  fail=0 err=0 skip=0 warn=0
test-refit-plumbing.R   pass=45  fail=0 err=0 skip=0 warn=0
test-bernoulli-coding.R pass=48  fail=0 err=0 skip=0 warn=0
test-double-bar.R       pass=155 fail=0 err=0 skip=0 warn=0
test-hmm-starts.R       pass=95  fail=0 err=0 skip=0 warn=0  (latent)
```

Gated tier: the 38 core files that skipped above, with
`FRMTMB_BRMS_FIT_TESTS=true` (`dev/refits-gated.log`): GATED ran 38 of
38, pass=3457 fail=0 err=0 skip=15 warn=0. The 15 skips are other
gates: drmTMB (13), the scale tier (1) and the fuzz tier (1). Run on
their own: `test-drmtmb-agreement.R` with `FRMTMB_DRMTMB_FIT_TESTS`
pass=131 fail=0 (it holds the `gr(cov = A)` animal model), and
`test-fuzz.R` with `FRMTMB_FUZZ` pass=2 fail=0. The scale tier was not
run (a cost row of the spline curves; no path of this lane).

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
4 threads; `dev/refits-tests/openblas-*.log`; `dev/ciharden-out/`
deleted after):

```
test-refit-plumbing.R   pass=45  fail=0 error=0 skip=0 warn=0
test-double-bar.R       pass=155 fail=0 error=0 skip=0 warn=0
test-bernoulli-coding.R pass=48  fail=0 error=0 skip=0 warn=0
test-formula-dot.R      pass=27  fail=0 error=0 skip=0 warn=0
```

The new numeric assertions over 20 perturbations y * (1 + 2e-15 k),
k = 1..20 (`dev/refits-ulp.R`; logs `refits-ulp.log`, reference BLAS,
and `refits-ulp-openblas.log`, the same maxima):

```
max relative |logLik diff|: factor vs numeric 0, against the exact
  density 7.4e-16, (1 + f ||) vs numeric 0, prec vs cov 4.9e-16
  (test tolerances 1e-10, 1e-9, 1e-10, 1e-8)
influence ratio |table - written| / |deletion move|: max 1.87e-06
  (test bound 1e-3)
```

Every suite in one session, three ways, every extension installed from
the worktree (`dev/refits-sessions.sh sess 2`, log `dev/refits-sess.log`):
SESS ran 27 of 27, fail 0 and err 0 in all 27, and no LEAK line without
`[with load of ...]`. Core forward/reverse/preload files=318/331,
pass=24436/24601/24704, skip=240/240/242; frmtmb.sample files=86/91,
pass=4238 each (the same 13 and 5 files without a result as the 0.82.0
record, `dev/round-20261111.md` line 393); the other extensions as in
the per-file tier, latent pass=698 with its changed test.

`R CMD check --as-cran`, built with vignettes inside
`dev/refits-check/` (deleted after; summaries
`dev/refits-check-*-summary.log`):

```
frmtmb:        Status: 2 NOTEs  (V8 math rendering; example timing
               frm_hazard_reads 5.04 s elapsed, an example this lane
               did not touch: load, per dev/lane-rules.md). Tests 38 min OK.
frmtmb.latent: Status: 1 NOTE   (V8 math rendering)
```

No brms-suite ledger row was edited; none of the hand-written verdict
rows names a construct this lane changed. The ledger is regenerated at
consolidation.

## Questions for the user

Questions 1 and 2 of the first round are settled by the fix round
(digit text refused as brms refuses it; the reviewer's rule for
`refit()`). Question 3 (`collapse` for brms's marginaleffects route)
is moot: the review checked 3188 against `nnet::multinom()`.

## Fix round (review FIX, 2026-11-13)

Done: F1 (1687), F2 (3430), F3 (NEWS BREAKING marks: getCall() after a
delta, influence() tables of s()/t2(), formula(fit) of `y ~ . - w`,
vcov_cluster(), plus 1259 digit text and the refit() reading), 1259
(digit text refused), the 3188 LOW (refit() refreshes the data), and
the two backlog filings (gp() deletion refits, MEDIUM; the 1567 gaps,
added to the LOW item). Files touched beyond the first round:
`R/bootstrap.R` (`boot_refits()`), `R/sugar.R`, `R/sandwich.R`,
`R/frame.R` (`extract_y()`), `R/influence.R` (doc), `NEWS.md`,
`dev/test-backlog.md`, `man/refit.Rd`, `man/influence.frmtmb_fit.Rd`,
`tests/testthat/test-refit-plumbing.R`, `tests/testthat/test-thres-refit.R`
(the character coding dropped from `refit_codings`).

New and changed tests, seen failing first:

```
test-refit-plumbing.R on the first-round build (plumbing-prefix.log):
  pass=52 fail=5 (digit text x2, the 1/2 reading, the workspace
  cluster, the refit data)
test-refit-plumbing.R on rellib-r21 (plumbing-before2.log):
  pass=5 fail=7 error=11 skip=0 warn=2 tests=14
test-refit-plumbing.R fixed (plumbing-after2.log): pass=57 fail=0, 14
```

Core, one file per process (`dev/refits-core2.log`): CORE2 ran 331 of
331, pass=24440 fail=0 err=0 skip=240 warn=0; in it
test-refit-plumbing 57/0, test-thres-refit 110/0,
test-bernoulli-coding 48/0, test-double-bar 155/0, test-formula-dot
27/0. Extension files that read the changed paths (`dev/refits-ext2.log`):
latent test-hmm-starts 95/0; sample test-add-criterion 41/0,
test-incl-thres-draws 43/0, test-loo 138/0 skip 2, test-prior-draws
93/0, test-readers-resp-me 106/0; spatial test-car-postfit 14/0.

OpenBLAS 0.3.32 with LAPACK, `OPENBLAS_CORETYPE=Haswell`, 4 threads
(`dev/refits-tests/ob32-*.log`; `dev/ciharden-out/` deleted after):
test-refit-plumbing 57/0, test-double-bar 155/0, test-bernoulli-coding
48/0, test-formula-dot 27/0, test-thres-refit 110/0, latent
test-hmm-starts 95/0; `refits-ulp.R` the same maxima as above
(`refits-ulp-ob32.log`).

## Merge note

After the parity2 merge, run `influence()` on an `s(exp(x))` fit:
parity2 makes the frame evaluate a smooth's term, and the frozen basis
(`smooth_random_frozen()`) calls `PredictMat()` on the frame, so the
two must agree on the column the smooth reads (review, merge dry runs).
The trees are CRLF; convert to LF before `git merge-file`.

## Fix round: R CMD check

`R CMD check --as-cran` of core after the fix round, built with
vignettes in `dev/refits-check/` (deleted after; summary
`dev/refits-check-core2-summary.log`): Status: 1 NOTE (V8 math
rendering of the HTML manual). Examples 90 s OK, tests 34 min OK.

## Fix round 2 (short check FIX, 2026-11-13)

`dev/reviews/2026-11-13-refits-short.md`. The short check's extra probe
rerun (`dev/reviews/refits-short-extra.R`, log
`dev/refits-n12-extra-after.log`, against
`dev/reviews/refits-short-extra-after.log`):

- **N1 (false refusal from the 3188 refresh).** `refit_data_matches()`
  takes `ignore_response`; `resolve_cluster()` sets it, so the call's
  data are compared on every column but the response, which a refit
  replaced. `update()` keeps the full comparison (a response changed
  in the call's data after the fit is still another model there).
  M1, seed 5, n = 200, `multinomial()`:
  first fix round: ERROR "... no longer hold the data ...";
  now `vcov_cluster(refit, ~ cl)` = 0.0523188 0.0514621 0.0493927
  0.0617483 = `vcov_cluster(refit, dd$cl)`.
- **N2 (silent, pre-existing).** `resolve_cluster()` falls back to the
  formula's environment only when the call's data were read (and lack
  the variable) or the call names no data. X1 (fit made in a function
  on local data, a global `cl` of length n): before 0.0162685 0.0645226
  0.0215026 (the global `cl`), now refused "could not find `cl` ...
  (d_local, which are gone)". The call's data are read with
  `fit$call[["data"]]`, not `$data`, which partial-matches `data2`.
- **`influence()` of a univariate refit (silent, pre-existing).** Fixed
  rather than filed: a deletion frame takes the fit's own response on
  the kept rows (`influence_set_response()`, `R/influence.R`) when
  `data` is not given and the fit has one response. U1/U2, seed 8,
  n = 30: max |diff| of the table against a fresh fit on the new
  response 1.850286 -> 3.33e-07; against the original fit's table
  1.06e-05 -> 1.850289.
- Filed LOW: the raw `Ya`, `Yb`, `Yc` columns of a `cbind()` response
  in the attribute stay old after `refit()`.

Unchanged against the short check's logs: `refits-short-f2.R` (S1 to
S11, `dev/refits-n12-f2-after.log`, identical) and
`refits-review-probe2.R` (V0 to V12, identical); `refits-short-3188.R`
differs only in the refit's `vcov_cluster(~ cl)`, which now computes
(0.052319 0.051462 0.049393 0.061748).

Tests: three new tests in `test-refit-plumbing.R` (gone data refused;
refit of a matrix response; influence of a refit). On the first fix
round's build: pass=57 fail=2 error=1
(`dev/refits-tests/plumbing-prefix2.log`, each new test fails); now
pass=60 fail=0, 16 tests (`plumbing-after3.log`).

Core, one file per process. The run started 3 at a time; after the
machine-wide notice (memory commit limit) I stopped my own driver and
its two running files after 152 of 331 had a result, and ran the 179
others 2 at a time (`dev/refits-core3b.log`; the 152 kept in
`dev/refits-core3-keep/`, all results in `dev/refits-core3-all.txt`).
No file failed for memory or disk. 331 of 331: pass=24443 fail=0 err=0
skip=240 warn=0; test-refit-plumbing 60/0, test-sandwich 58/0,
test-influence-plot 21/0, test-refit-data 122/0, test-thres-refit
110/0, test-bernoulli-coding 48/0.

OpenBLAS 0.3.32 with LAPACK, `OPENBLAS_CORETYPE=Haswell`, 4 threads,
the same six files (`dev/refits-tests/ob32b-*.log`; `dev/ciharden-out/`
deleted after): 60/0, 58/0, 21/0, 122/0, 110/0, 48/0.

`R CMD check --as-cran` of core after fix round 2 (built with vignettes
in `dev/refits-check/`, deleted after; summary
`dev/refits-check-core3-summary.log`): Status: 1 NOTE (V8 math rendering
of the HTML manual). Examples 46 s OK, tests 25 min OK.
