# Review: lane brmsnames, plan item 2.6c

Adversarial reviewer pass on `frmtmb-wt-brmsnames`, base `aa9227e`
(frmtmb 0.58.0, frmtmb.sample 0.6.0). The brief was to FALSIFY.

Libraries. BASE arm: the shared reference build
`C:/Users/adf44/source/r/rellib-r3`, read only. LANE arm: the worker's
`C:/Users/adf44/source/r/brmsnames-lib`, read only. Nothing was
installed anywhere. Before any measurement,
`dev/brmsnames-rev-verify.R` compared every top-level function of the
worktree sources against the installed namespaces by deparsed body:
frmtmb 968 functions, 1 differs (`covstruct_registry`, a list the
script reads only in its first assignment, so a false positive);
frmtmb.sample 193, 0 differ; frmtmb.learn 44, 0 differ. So the lane
build is the worktree. Real brms fits used brms 2.23.0 with `pinlib`
ahead of the user library (StanHeaders 2.32.10 asserted), compiled
fresh; no `FRMTMB_STAN_CACHE` hit is involved in any brms number here.
No git write operation. No package file edited.

Scripts added by this review, all under `dev/`, prefix
`brmsnames-rev-`: `verify.R`, `brmsgrep.R`, `data.R`, `brms.R`,
`frm.R`, `shapes.R`, `silent.R`, `varcorr-se.R`, `mutants.R`,
`guards.R`, `collide.R`, `collide-brms.R`, `colon.R`, `rcontent.R`,
`excluded.R`, `excluded-base.R`, `interop.R`, `r2mv.R`. Logs are in
`dev/brmsnames-rev-log/`. Fitted objects are in `dev/stan-cache/`
(`brmsnames-rev-*.rds`, gitignored).

Verdict: **not mergeable**. One BLOCKER (a silent wrong answer on
brms's own spelling, new in this lane for the `b_` form) and four
MAJORs, two of which falsify the lane's headline claims (names are
brms's; collisions are unreachable). The rest of the work holds up
well under independent checks: every number is the same number under
the renaming, the `VarCorr()` standard errors are right to 5e-9, and
the correlation order is right at K = 4.

---

## Findings

### BLOCKER 1. `hypothesis()` reads an interaction name as R's `:` operator

Construction: `dev/brmsnames-rev-colon.R base|lane`, data seed 8,
`y ~ x * f`, n = 400; and the C1 draws against the C1 brmsfit.

On a fit, `fixef()` gives `x` 0.2415 and `x:fe` 1.5794:

| call | base | lane |
| --- | --- | --- |
| `hypothesis(fit, "x:fe > 0")` | 0.2415 | 0.2415 |
| `` hypothesis(fit, "`x:fe` > 0") `` | 1.5794 | 1.5794 |
| `hypothesis(fit, "b_x:fe > 0")` | error "object 'b_x' not found" | **0.2415** |

On draws, same model family (C1, sampler seed 3):
`hypothesis(ds, "x:fe > 0")` gives Estimate 0.2787, which is the mean
of `b_x`; the mean of `b_x:fe` is 0.0455. On the brms fit of the same
model and data, `hypothesis(b, "x:fe > 0")` gives 0.0719, the mean of
`b_x:fe`.

Why it is silent: `x:fe` parses as `b_x:b_fe`, a sequence from `b_x` to
`b_fe`, which has length 1 whenever the two differ by less than 1, so
the length check in `hypothesis.frmtmb_fit()` passes.

brms rewrites `:` to `___` before evaluating
(`brms:::eval_hypothesis()`, `rename(h, c(":", "[", "]", ","), ...)`).
`draws_hypothesis_coef()` in this lane already does that for
`scope = "ranef"`; the standard scope on fits and draws does not.
`variables()` now lists `b_x:fe`, and the lane tells a user to write
brms's strings, so the `b_` spelling is a new route to the wrong
answer. The bare spelling was already wrong at base.

### MAJOR 1. The names are not brms's once a name has a character brms renames

The 107 of 107 was measured on four models whose names contain only
letters and underscores. `par_name_bare()` drops parentheses and
nothing else. brms applies `brms:::rename()` (space, parentheses,
brackets, comma and quotes dropped; `+ - * / ^ =` become
`P M MU D E EQ`), replaces whitespace in group levels with `.`, pastes
the levels of an interaction group with `_`, and drops `_` and `.` from
response names.

Construction: `dev/brmsnames-rev-data.R` (data seed 2026), real brms
fits from `dev/brmsnames-rev-brms.R` (chains 2, iter 60, seed 1),
frmtmb.sample draws from `dev/brmsnames-rev-frm.R` (seed 3), compared
by `dev/brmsnames-rev-shapes.R`.

| term | brms | frmtmb (lane) |
| --- | --- | --- |
| `I(x^2)` | `b_IxE2` | `b_Ix^2` |
| factor level `c-d` | `b_fcMd`, `b_x:fcMd` | `b_fc-d`, `b_x:fc-d` |
| group level `lvl 1` | `r_h[lvl.1,Intercept]` | `r_h[lvl 1,Intercept]` |
| group `g:h2` | `sd_g:h2__Intercept`, `r_g:h2[1_p,Intercept]` | `sd_gh2__Intercept` (fit), `r_g:h2[1:p,Intercept]` |
| mv response `y_a` | `b_ya_Intercept`, `sd_g__ya_Intercept`, `r_g__ya[1,Intercept]`, `bayes_R2` row `R2ya` | `b_y_a_Intercept`, `sd_g__y_a_Intercept`, `r_g__y_a[1,Intercept]`, `R2y_a` |

The divergence reaches `fixef()`, `coef()`, `ranef()` dimnames,
`VarCorr()` dimnames, `posterior_summary()`, `as_draws_df()` and
`bayes_R2()`. `I(x^2)` and `poly(x, 2)` are ordinary formulas, and a
name such as `b_Ix^2` cannot be written unquoted in a hypothesis.

### MAJOR 2. Name collisions are still reachable, and silent on the fit

The lane removed the shadow-note machinery "on the argument that `b_`
prefixes make collisions unreachable", and the `hyp_env_vals()` comment,
`?hypothesis` and NEWS say no collision is reachable. Finding 10 of
`dev/brmsnames-findings.md` says "brms has the same collision". Both
are false.

Construction: `dev/brmsnames-rev-collide.R base|lane` (data seed 12)
and `dev/brmsnames-rev-collide-brms.R` (brms's own parse, no compile).

| model | `hypothesis()` string | returns | intended | brms |
| --- | --- | --- | --- | --- |
| K1 `y ~ sigma_Intercept + (1 \| g)` | `sigma_Intercept = 0` | 1.0733 (log sigma intercept) | -1.4435 (the covariate) | no collision: brms's residual SD is `sigma`, the covariate is `b_sigma_Intercept` |
| K4 `mvbf(y ~ x_z, y_x ~ z)` | `y_x_z = 0` | 0.4634 (`y_x`'s `z`) | 0.7899 (`y`'s `x_z`) | no collision: resp `y_x` is `yx`, so `b_y_x_z` and `b_yx_z` |
| K5 `y ~ Intercept + x` | `Intercept = 0` | 0.1518 (the covariate) | 1.2357 (the intercept) | refused: "Internal renaming led to duplicated names" |
| K3 `(1 \| g:h2) + (1 \| gh2)` | `sd_gh2__Intercept` (class NULL) | 2.2602 (the `g:h2` SD) | 0.0003 (the `gh2` SD) | no collision: `sd_g:h2__Intercept` |
| K2 `bf(y ~ sigma_z + (1 \| g), sigma ~ z)` | `sigma_z = 0` | 0.0668 (sigma's `z`) | 2.1032 (mu's `sigma_z`) | brms fits it and keeps both, as `b_sigma_z` and `b_sigma_z__1` (C6 brmsfit) |

Every one of these numbers is identical at base, so on the fit this is
not a regression. What the lane adds: the draws labels now duplicate
(`brms_par_labels()` returns `b_sigma_Intercept`, `b_sigma_z`,
`b_y_x_z` and `b_Intercept` twice in K1, K2, K4 and K5), and on draws
every accessor then fails with posterior's "Duplicate variable names
are not allowed" (C6: `fixef`, `coef`, `posterior_summary`, `summary`,
`as_draws_df`, `as.array`, `as.matrix`, `hypothesis(scope = "coef")`),
where brms's fit of the same model answers all of them. The rule to
port is brms's: build every name through one renamer and refuse a
duplicate (`rename(check_dup = TRUE)`), or suffix it as brms's stanfit
repair does.

### MAJOR 3. frmtmb.sample returns a log-sigma row that brms does not have

Construction: `dev/brmsnames-rev-shapes.R`, C1, C2, C5, C7, C8 (scalar
`sigma`; C8 is the plain `y ~ x + (1 | g)`) and C4 (mv, no rescor).

In brms, a residual SD that is not modeled is `sigma`, which is not in
`fixef()` or `coef()`. frmtmb.sample names the sampled log sigma
`b_sigma_Intercept`, so `fixef(ds)` has 3 rows where brms has 2 (C2),
`coef(ds)$g` has a `sigma_Intercept` coefficient slice brms has not,
and `posterior_summary(ds)` and `as_draws_df(ds)` carry the extra
column. On C4 the extra rows are `sigma_y_a_Intercept` and
`sigma_y2_Intercept`, where brms has `sigma_ya` and `sigma_y2`. This is
the reason K1 in MAJOR 2 collides.

The pinning test cannot see it. `bo_shim()` in
`extensions/frmtmb.sample/tests/testthat/test-brms-output.R` builds the
brms side from `brms::bf(y ~ x + (1 + x | g), sigma ~ 1)` unless
`simple_sigma = TRUE`, which makes brms name its log sigma
`b_sigma_Intercept` too. NEWS says `coef()` "broadcasts every
population-level coefficient, as brms's does"; brms's `coef()` on the
model the user wrote has no sigma column.

### MAJOR 4. `ranef(ds)` and `coef(ds)` silently return NA where 0.6.0 returned values

Construction: `dev/brmsnames-rev-excluded.R` (lane) and
`dev/brmsnames-rev-excluded-base.R` (base), data seed 21, sampler seed
4, chains 1, iter 100.

| model | lane | base |
| --- | --- | --- |
| `y ~ x1 + rr(x1 + x2 \| g, d = 1)` | `ranef(ds)$g` 36 of 36 estimates NA; `coef(ds)$g` 144 of 144 cells NA | values, for example level 1 `x1` Estimate -2.218 |
| `(1 \| gr(id, cov = A)) + (1 \| id)` | second block `ranef(ds)[["1 \| id"]]` 12 of 12 NA; `coef(ds)` 48 NA | values |

No warning or refusal. `draws_ranef_layout()` fills NA for any block
without `r_` labels, following brms's NA fill for missing draws, but
brms has neither of these models (brms has no `rr()`, and it refuses
the duplicated animal-model term). The ML `ranef(fit)` answers both.
Either compute them per draw as base did, or refuse by name.

### MINOR 1. The guards do not catch a coefficient mislabel or a correlation reorder

Construction: `dev/brmsnames-rev-mutants.R none|rlev|rcoef|corord`, in
memory through `assignInNamespace()`, running
`test-brms-output.R`, `test-draws-methods.R` and `test-brms-names.R`,
plus an independent check (data seed 5, sampler seed 8,
`y ~ x + (1 + x + z + w | g)`) that compares each `r_` column and each
`VarCorr(ds, summary = FALSE)` correlation with the ML `ranef()` and
`cov2cor(varcorr_matrices())` of a fit set to the same draw.

| mutation | independent check | brms-output | draws-methods | brms-names |
| --- | --- | --- | --- | --- |
| none | 0 of 120 r_ wrong, 0 of 36 cor wrong | 22 / 0 | 145 / 0 | 56 / 0 |
| levels reversed | 120 of 120 r_ wrong | 18 / **4** | 145 / 0 | 56 / 0 |
| coefficients reversed | 120 of 120 r_ wrong | 22 / 0 | 145 / 0 | 56 / 0 |
| correlations column-major | 12 of 36 cor wrong | 22 / 0 | 145 / 0 | 56 / 0 |

(PASS / FAIL.) A coefficient mislabel passes because the shim takes its
names from `brms_par_labels()` and `bo_brms_order()` then reorders by
the label, so brms's positional reshape reads the mislabeled columns in
the order that gives the right coefficient. That is claim 1's
laundering, demonstrated. The current code is correct on both counts;
the tests would not notice if it stopped being.

Two more guards with no pin, from `dev/brmsnames-rev-guards.R`:
zeroing the derivative of `residual__` (Est.Error prints 0 against a
delta-method 0.05244) passes `test-brms-names.R` 56 of 56; and no test
calls `log_lik(pointwise = TRUE)` or `add_point_estimate = TRUE`,
which the lane added as refusals.

### MINOR 2. Old positional calls on a fit change answer without an error

Construction: `dev/brmsnames-rev-silent.R`, six models. On every one,
`fixef(fit, TRUE)` was the flattened numeric vector at base and is the
per-dpar list now; `ranef(fit, TRUE)` carried `condSD` at base and does
not now. This is brms's slot meaning and NEWS says `flatten` and
`condVar` must now be named, so it is recorded, not ranked higher.
`test-brms-names.R` line 110 asserts the new meaning.

### MINOR 3. Absolute tolerances in the new test file

`tests/testthat/test-brms-names.R` lines 74, 76, 81, 85, 87, 160, 163,
165 and 171 use fixed `tolerance =` values (1e-12, 1e-4, 1e-6, 1e-8,
1e-10), which `dev/lane-rules.md` forbids. The line 81 comparison
(`VarCorr()` Est.Error against `hypothesis()` at 1e-4) is also not an
independent check: both use the same finite-difference step rule.

### MINOR 4. Multivariate `VarCorr()` has no `residual__`

C4: brms's `VarCorr()` has `residual__` with rows `ya` and `y2`;
frmtmb's has only `g`, on the fit (by construction of
`varcorr_layout()`) and on draws. A shape difference.

Seen, pre-existing and unchanged, not ranked: `bayes_R2(ds)` on the
two-response C4 model returns one row, `R2y_a`, where brms returns
`R2ya` and `R2y2` (`dev/brmsnames-rev-r2mv.R`: identical at base and
lane, 0.1088).

### NIT

- `hyp_eval()` still says "a natural-scale name a coefficient has taken
  over carries a leading dot" (`R/confint.R` about line 2284); the dot
  aliases are gone.
- `ranef_pick()`'s refusal calls the term label "VarCorr()'s key"
  (`R/methods-fit.R` about line 920); `VarCorr()` is keyed by group now.
- `?fixef` says `hypothesis()`'s own spelling is `sigma_Intercept` and
  reads `sigma_(Intercept)` backquoted (`R/methods-fit.R` lines
  686 to 689); it is `b_sigma_Intercept` and `` b_sigma_(Intercept) ``.
- A smooth's null-space coefficient is `b_sx.fx1`; brms's is
  `bs_sx_1`.

---

## Claims I failed to break

1. **Same numbers under the renaming (claim 6).**
   `dev/brmsnames-rev-silent.R base|lane|compare`, data seed 31,
   sampler seed 5, six models (distributional with a sigma group term,
   mv with `|p|`, nonlinear, cumulative, `s(x) + (1 | g)`, `gp(x)`):
   `fixef()`, `fixef(flatten = TRUE)`, `ranef()`, `coef()`, the
   per-block matrices (`VarCorr()` at base, `varcorr_matrices()` in the
   lane), `confint()` and the whole draws matrix are `identical()`
   between arms on all six. Every lane hypothesis-vocabulary value
   equals a base value; the only lane name without one is
   `sd_g__sigma_Intercept`, the SD the base could not reach. The old
   positional `hypothesis(fit, h, 0.1)` now errors on `class`, loudly.
2. **Each `r_` column holds what its name says.**
   `dev/brmsnames-rev-rcontent.R`: draws 1, 50 and 100 of C1, C3
   (`|ID|` across mu and sigma), C4 (mv), C5 (nonlinear) and C7 (two
   blocks on `g` plus `g:h2`), against the ML `ranef()` at that draw:
   597 cells, 0 wrong, 0 names missing. On the silent-change models, 180
   more, 0 wrong. The mutants above show the check has power.
3. **`VarCorr()` Est.Error on a fit (claim 2).**
   `dev/brmsnames-rev-varcorr-se.R`, data seed 77. Independent of the
   lane in three ways: numDeriv's Richardson Jacobian, the covariance
   from `vcov(fit, full = TRUE)` by row name, and lane entries matched
   to the reference by ESTIMATE value, not by name. Worst relative
   difference: `us` K = 4 plus a second group, ML 1.1e-10, REML 4.3e-11;
   correlation near +1 with a near-zero SD 4.5e-9; `diag` 2.9e-11;
   `cs` 7.2e-11; `ar1` 4.0e-10; mv `|p|` 5.3e-11; sigma `|q|` 3.4e-11.
   0 unmatched entries.
4. **brms's correlation order at K = 4 (claim 5).** 0 of 36 cells wrong
   (MINOR 1, "none" row); 12 of 36 wrong under the column-major
   mutation, so the check can say no. The real brms fits of K = 3 (C1)
   and K = 4 (C2) give `VarCorr()` dimnames identical in shape and names
   to frmtmb.sample's.
5. **Accessor shapes against a real brmsfit (claim 1), where the names
   agree.** `dev/brmsnames-rev-shapes.R` on C1 to C8: `ranef()` and
   `ranef(summary = FALSE)` dims and dimnames, `VarCorr()` and
   `VarCorr(summary = FALSE)` element names and dimnames, the
   `hypothesis()` object's elements and columns, including
   `scope = "ranef"` and `"coef"`, `summary()` columns, and
   `conditional_effects()` columns match brms on every model, except
   where MAJOR 1, MAJOR 3 or MINOR 4 applies. Variable ORDER differs,
   as recorded in the lane's "not fixed" list.
6. **Excluded blocks and the duplicate fallback (claim 4).** `rr()` and
   `s()` blocks get `b[i]` and no `r_`; the animal model gets 12 `r_`
   and 12 `b[i]`, no duplicate label, and `VarCorr()` keys the second
   block by term label (`dev/brmsnames-rev-excluded.R`). brms refuses
   the animal model, so there is no brms name to miss. The consequence
   for `ranef(ds)` is MAJOR 4.
7. **Blast radius outside tests (claim 7).** No extension `R/` code
   calls `VarCorr()`, `variables()` or `hypothesis()`, or reads a draws
   column by an old name. Core `R/` reads `varcorr_matrices()` for
   `print()`, `summary()` and `confint_varcorr()`.
   `dev/brmsnames-rev-interop.R base|lane|compare`, data seed 3:
   `insight::find_parameters`, `get_parameters`, `get_variance`,
   `get_statistic`, `find_random`, `emmeans::emmeans`, `print(fit)` and
   `print(summary(fit))` are `identical()` between arms, and
   `marginaleffects::avg_slopes()` differs only in the attached model
   object. parameters, performance and broom.mixed are not installed.
   The lane's own `R CMD check` logs show vignettes re-built and
   examples with `--run-donttest` OK for frmtmb and frmtmb.sample.
8. **Guards that do fail closed (claim 8).** `fit_refuse_draws_args()`
   made a no-op fails 6 assertions in 2 blocks of `test-brms-names.R`.
   `posterior_summary()` on a fit sent to the default method gives
   "is.atomic(x) is not TRUE" in a standalone process, which the
   test's "needs posterior draws" pattern rejects. (Registering the
   mutated method did not take effect inside `test_file()` in my run, so
   this one is shown by construction, not by a failing count.)

## What to fix before merge

1. BLOCKER 1: rename `:` (and brms's other characters) in the hypothesis
   string and the evaluation environment together, as
   `brms:::eval_hypothesis()` does, on fits, `frm_multiple()` and draws.
   Pin it with a model where `abs(b_x - b_x:fe) < 1`, so the old
   behavior returns a number rather than an error.
2. MAJOR 1: one renamer used by `brms_coef_names()`, `brms_re_rnames()`,
   `brms_r_labels()` and the group and response names, ported from
   `brms:::rename()`, level whitespace to `.`, interaction levels joined
   by `_`, and `_` and `.` dropped from response names. Extend
   `dev/brmsnames-naming.R` with `I(x^2)`, `poly(x, 2)`, a level with a
   space and a hyphen, `(1 | g:h)` and a response with an underscore.
3. MAJOR 2: refuse a duplicate brms name at `frm()` with brms's message,
   or suffix it as brms does; correct the comment in `hyp_env_vals()`,
   `?hypothesis`, NEWS, and finding 10 of the findings file.
4. MAJOR 3: decide the name of an unmodeled log sigma in draws so that
   `fixef()` and `coef()` on draws match brms's for the model the user
   wrote, and make the `test-brms-output.R` shim use that model.
5. MAJOR 4: compute or refuse; do not return NA silently.
6. MINOR 1 and 3 in the same pass: a K = 4 and a coefficient-mislabel
   pin that does not route through `brms_par_labels()`, a
   `residual__` Est.Error pin, a `log_lik(pointwise = TRUE)` refusal
   test, and ratio tolerances.
