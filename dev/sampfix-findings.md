# Lane sampfix: findings

Base 1f40800d (frmtmb 0.65.0, frmtmb.sample 0.13.0). The reference
build is `C:/Users/adf44/source/r/rellib-r3`; the lane build is
`C:/Users/adf44/source/r/wt-sampfix-lib` ahead of it. Every number
below is in `dev/sampfix-log/10-summary.txt`, which
`dev/sampfix-10-summary.R` generates from the logs; each item names the
script and its seed.

## What changed

### 1. Laplace draws were read in the wrong layout (silent wrong answer)

**Why they were NaN.** `frm_sample(laplace = TRUE)` draws hold the
outer parameters alone, `b_Intercept b_x sigma theta_1 lp__` for
`y ~ x + (1 | g)`. Every draws method mapped a draw back onto the
parameter template with `draws_par_index()`, which assumes the full
layout: `beta`, then `b` (six columns), then `betad`, then `theta`. So
`b` was read from `sigma`, `theta_1`, `lp__` and three positions past
the end (`NA`), and `sigma` and `theta` were read past the end. Where a
prediction reached an `NA` it was NaN; where it reached a real column it
was finite and wrong.

Measured (`dev/sampfix-01-repro.R`, seed 1212, sampler seed 3, N = 30 in
6 groups):

- reference: `posterior_epred()`, `posterior_linpred()` and
  `posterior_predict()` gave 4500 non-finite cells of 4500 on
  `y ~ x + ar(t, g) + (1 | g)` and 3000 of 4500 on `y ~ x + (1 | g)`.
- the finite case is worse. At `newdata` with rows in groups 1 and 2
  (`dev/sampfix-06-probes.R`), the reference returned finite values, and
  each row's "group effect" was exactly another column: group 1's was
  `theta_1` (to 2.2e-16) and group 2's was `lp__` (to 7.1e-15).
- `fitted()`, `predict()`, `residuals()`, `predictive_error()` and
  `bayes_R2()` returned the same kind of numbers (the last with only a
  base `STATS is longer` warning); `pp_check()`,
  `predictive_interval()`, `hypothesis()` on `sd_g__Intercept` and
  `pp_mixture()` stopped in `quantile()` or bayesplot errors
  (`dev/sampfix-03-sweep.R`).
- an `mi()`-only model (no group-level block): the reference's
  `posterior_epred(resp = "y")` gave 300 non-finite cells of 4500.

**Decision: read laplace draws in their own layout, refuse by name what
needs the integrated values, compute the rest exactly.** brms has no
Laplace route (its sampler always visits `r_`), so brms's semantics give
no answer to copy; the choice is about what is right.

Two ways to compute the refused quantity were measured
(`dev/sampfix-04-laplace-validate.R`: data seed 2024, 12 groups of 10,
4 chains of 1000 post-warmup draws, sampler seed 11 on both routes,
Laplace-conditional draws seed 5), against the full draws' in-sample
`posterior_epred()`:

| construction | max mean diff / sd | spread ratio |
|---|---|---|
| each draw's conditional modes `b_hat(theta)` | 0.0631 | 0.2305 to 0.7368 |
| `b` drawn from `N(b_hat, H^-1)` | 0.0650 | 0.9476 to 1.0518 |

The conditional modes are wrong: they drop the uncertainty of `b`
given `theta` and understate the spread by up to a factor of four. The
Laplace-conditional draw is right to within Monte Carlo error on this
gaussian model, where the Laplace approximation of `b` given `theta` is
exact (for another family it is itself an approximation). It is NOT
shipped. **The user decided this on 2026-09-29:** `laplace = TRUE`
samples an approximate posterior whose use is checking the
approximation, and anyone who wants predictions samples without it. Two
further reasons: `posterior_epred()` would depend on the random seed,
and `log_lik()`, `ranef()`, `coef()` and `loo()` refuse laplace draws
by design, so filling `b` in for the predictive methods alone would make
the methods disagree about what a draw contains.

What the lane does instead:

- `draws_par_index(fit, laplace =)` and `draws_index(x)` give the
  layout the draws actually have; every reader goes through it.
- `draws_fit_at()` fills what the draws do not hold with `NA`.
- `draws_laplace_probe()` evaluates the call at the first draw with
  those values at `NA` and at `0`; a result that changes read them, and
  the call is refused by name. The probe keeps the caller's random
  stream (tested). Wired into `posterior_epred()`,
  `posterior_linpred()`, `posterior_predict()` (so `fitted()`,
  `predict()`, `residuals()`, `predictive_*()`, `bayes_R2()` and
  `pp_check()`), `hypothesis()` where it needs the fit, `pp_mixture()`
  and `conditional_effects()`.
- `draws_require_b()` (`log_lik()`) now asks whether anything was
  integrated, not whether the model has a group-level block.

Why a probe and not a rule on the layout: `re_formula = NA` drops
`(1 | g)` and keeps a smooth, whose coefficients are also in `b`, and a
distributional parameter with no group-level term reads no `b` at all.
The probe answers exactly for each call. The guard-absent cases are
tested and hold bit for bit against full draws with the same outer
values: `posterior_epred()`, `posterior_linpred()` and (same seed)
`posterior_predict()` at `re_formula = NA`, with and without `newdata`;
`posterior_epred(dpar = "sigma")`; `hypothesis()` on
`sd_g__Intercept`; `VarCorr()`; `conditional_effects()` at its default
and at `re_formula = NULL` with a seed; `posterior_epred(resp = "x")` on
the `mi()` model. On real draws against full draws of the same model
(`dev/sampfix-04-laplace-validate.R`), the population-level
`posterior_epred(re_formula = NA)` agrees to 0.0203 of a posterior sd
with spread ratio 0.9512 to 1.0068, and `conditional_effects()` to
0.0465 of its se (full-draws bulk ESS on `b_Intercept` 1466).

`conditional_effects()` used to refuse every laplace draws object with
a group-level block; it now draws its default curves and refuses a
smooth.

### 2. `laplace = TRUE` with nothing to integrate

tmbstan evaluates `-obj$env$random`, which is `NULL` there. A Laplace
approximation over nothing is the model itself, so the call now says so
in a message and samples on the full route; the draws are identical to
those of the call without `laplace` (tested, `identical()`).

Found while there: on a REML fit sampled as it stands (`prior =
"flat"`), the objective integrates `beta` too, the draws hold
`betad theta lp__`, and the reference labeled them
`b_Intercept b_x sigma` (`dev/sampfix-06-probes.R`). That is refused
now. With the default priors the objective is rebuilt with `beta`
sampled, and it still samples (tested).

### 3. `pp_check()`'s `loo_*` types

Built as `brms:::pp_check.brmsfit()` builds them: `log_lik()` on the
same draws as the predictions, `r_eff` from `loo::relative_eff()` on the
chain structure when every draw is used and on one chain otherwise (with
brms's warning for an `NA`), `loo::psis(-ll, r_eff)`, then `lw` as
`weights(log = TRUE)` or the psis object itself. Where `log_lik()`
refuses (laplace draws, `newdata`, `re_formula`, a likelihood that does
not factor by row) the type refuses before any prediction, naming
`pp_check()` and quoting `log_lik()`.

Validated against brms 2.23.0 (`dev/sampfix-05-ppcheck-brms.R`, data
seed 31 as `dev/arcovsample-rev-07-ppcheck.R`, brms sampler seed 31, 2
chains of 500): a frmtmb draws object built from brms's own draws and
brms's own stanfit, so `log_lik()` agrees exactly (max abs diff 0).
bayesplot's `ppc_loo_*` functions were traced to capture their inputs.
For `loo_pit_overlay`, `loo_pit_qq`, `loo_intervals` and `loo_ribbon`,
over all 1000 draws and over `draw_ids = seq(1, 1000, by = 7)`: the same
argument is passed, `y` is identical, and the log weights and the
Pareto k values differ by exactly 0.

`loo_pit` is not in `bayesplot::available_ppc()` (bayesplot deprecated
`ppc_loo_pit()`), and brms refuses it as an invalid type. The draws
method did not check the type; it now refuses as brms does, in brms's
words, as core's fit method already did.

### 4. `@` behind `%||%`

The sweep found three reads of `x$stanfit@sim$chains %||% 1L`:
`draws_raw_array()`, `nchains.frmtmb_draws()` (which
`draws_derived_matrix()` calls) and `draws_chain_id()` (already guarded
at 0.65.0). All three now go through `draws_nchains()`. The only other
`@` read in the package is `check_stan_draws()`, on a stanfit it was
just handed. On the reference, 26 calls of the 75 in the sweep died on
`@` for a `stanfit = NULL` object, and 2 more (`nuts_params()`,
`log_posterior()`) on bayesplot's "no applicable method"; on the lane,
0 died on `@`. (The first version of this file said 28, counting those
2 as `@`; corrected in the nits round, `dev/sampfix-10-summary.R`.)
`nuts_params()` and
`log_posterior()`, which read the stanfit itself, refuse by name
instead of bayesplot's "no applicable method".

### 5. Ordinal and `cs()` draw names

The draws stored `tau_raw_k` (cumulative: the first threshold and log
increments) and `bcs<j>_k`. They now store brms's names and brms's
values: `b_Intercept[k]` holding the thresholds (`b_Intercept[g,k]`
under `thres(gr = )`, `b_<resp>_Intercept[k]` in a multivariate model)
and `bcs_<column>[k]`.

- Core: each ordinal family declares the inverse of its threshold map,
  `post$ord_thresholds_raw` (`ord_raw_from_tau()`), beside the forward
  map; `thres(gr = )` declares its per-group inverse; `brms_fixef_rows()`
  carries it as `inv` on each block (`identity` for `cs()`).
- frmtmb.sample: `draws_ordinal_cols()` reads those blocks;
  `draws_to_natural()` maps and renames them when the draws are stored
  and maps them back for the model; `draws_outer_cols()` and
  `draws_fixef_ordinal()` read the new names. A family with no declared
  inverse keeps the internal names, so no column carries a name its
  values do not have. frmtmb.sample's floor rises to the frmtmb that
  carries this at consolidation (other lanes need new core exports), so
  the names do not depend on which core is installed.

Checked on four fits (`dev/sampfix-07-ordinal.R`, data seed 405, sampler
seed 3): cumulative, `sratio` with `cs(fc)`, cumulative with
`thres(gr = grp)`, and a two-response cumulative model.
`variables(ds)` is `variables(fit)` plus `lp__` on all four. Against the
reference build at the same seed, the stored draws mapped the lane's
way and `fixef()` differ by exactly 0, and `posterior_epred()` by 0 on
these four. The reviewer's 11 fits (`dev/sampfix-rev-04-ordinal.R`)
found `posterior_epred()` 1 ulp (2.2e-16) apart on 4 of them, all with
cumulative thresholds, which is consistent with the round trip through
`c(tau1, log(diff(tau)))` not returning the sampled log increments bit
for bit (not measured separately); so the claim is
"within 1 ulp" for `posterior_epred()`, and exactly 0 for the stored
draws and `fixef()`. The reference's own draws object (old names) read
by the lane build gives the same `posterior_epred()` and `fixef()`,
difference 0.

Readers by name that were checked: `posterior_summary()`,
`as_draws_*()`, `hypothesis()` (including `"Intercept[1] <
Intercept[2]"`), `summary()`, `mcmc_plot()`, `fixef()`, `ranef()` (no
ordinal columns), `VarCorr()`. What stores the names: a saved draws
object. It keeps working (tested); code that selected `tau_raw_1` from
the draws must now select `b_Intercept[1]` (NEWS says so).

`?variables` loses its exception sentence and now names the one real
gap left (see "Found and not fixed").

### 6. One-parameter model

tmbstan's model is `vector[N] y`, and its init sanitizer passes a
numeric init as `y`. rstan reads a length-one numeric as a scalar, so
every chain died "no more scalars to read" (reference: 2 of 3 init
routes fail, `init = "random"` worked; `dev/sampfix-02-onepar.R`, seed
1212). `stan_init_arrays()` passes each numeric init as a
one-dimensional array, which the sanitizer keeps; all three routes
sample on the lane.

## The sweep of every exported draws method

`dev/sampfix-03-sweep.R` calls 75 methods on laplace draws of
`y ~ x + (1 | g)` and on full draws with `stanfit = NULL` (seed 1212,
sampler seed 3). 52 of 150 results changed (45 before the nits round,
which reworded the refusals). Logs:
`dev/sampfix-log/03-ref.txt`, `dev/sampfix-log/03-lane.txt`. On the lane
no call returns a non-finite number except two that are right:
`hypothesis(ds, "x > 0")`'s `Evid.Ratio` is `Inf` because 300 of 300
draws are above 0 (brms's convention), and `loo_compare(ds, ds)`'s
`p_worse` is `NA` for a model compared with itself
(`dev/sampfix-09-nonfinite.R`). No call ends in an internal error.

## Tests

New in frmtmb.sample, each seen to fail on the reference build (the new
test file run against `rellib-r3`, `dev/sampfix-runtest.R ref`):

Whole suites on the lane build, one file per R process, with
`NOT_CRAN=true` and `FRMTMB_BRMS_FIT_TESTS=true`
(`dev/sampfix-suite.sh`, logs under `dev/sampfix-log/suite-*`): core,
182 files, 14382 pass, 0 fail, 0 error, 0 warn, 14 skip (13 in
`test-drmtmb-agreement.R`, gated on `FRMTMB_DRMTMB_FIT_TESTS`, 1 in
`test-fuzz.R`, gated on `FRMTMB_FUZZ`); frmtmb.sample after the nits
round, 40 files, 333 tests, 2324 pass, 0 fail, 0 error, 0 warn, 1 skip
(`test-scale.R`, gated). The last change of the round (the probe's
re-raise going through `frm_stop()`, which `test-conditions-census.R`
requires) was followed by a rerun of the 8 files it can reach, not of
all 40;
frmtmb.eam `test-sampling.R`, 97 pass.

`R CMD check --as-cran` (frmtmb before the nits round, which did not
touch core; frmtmb.sample rerun after it), built with vignettes inside
`dev/sampfix-check/` (`dev/sampfix-check.ps1`): frmtmb `Status: 1 NOTE`
(the environmental V8 math-rendering NOTE on the HTML manual);
frmtmb.sample `Status: OK`.

| file | reference | lane |
|---|---|---|
| `test-laplace-draws.R` | 14 pass, 18 fail, 10 error | 90 pass |
| `test-ppcheck-loo.R` | 4 pass, 4 fail, 4 error | 27 pass |
| `test-draws-no-stanfit.R` | 0 pass, 2 error | 21 pass |
| `test-ordinal-draws-names.R` | 0 pass, 3 fail, 4 error | 27 pass |

Changed: `test-draws-methods.R` (laplace-shaped `conditional_effects()`
now equals the full draws' curve instead of refusing),
`test-sample-direct.R` and `test-sratio-draws.R` (select the thresholds
as `b_Intercept[k]`).

The frmtmb.sample brms suite (`test-brms-suite-*.R`) passes with the
change, so no ledger row there flips. One ledger row's REASON is stale:
`brmsfit-methods:685` (core tier, `pp_check(type = "loo_pit_qq")` on an
ML fit) cites "dev/brmsnames-findings.md not fixed 8", the draws
method's missing `lw`, which this lane fixes. The row stays non-pass,
because the fit method refuses by design (an ML fit has no posterior
draws); its class should move from defect to a no-draws divergence.
`dev/brmsport-*` was not edited.

## Found and not fixed

1. **`conditional_effects()` on draws ignores an observed level named in
   `conditions` under `re_formula = NULL`.** Pre-existing, full draws as
   well as laplace ones. With `conditions = data.frame(g = factor(1))`
   or `g = 2`, shifting `r_g[1,Intercept]` or `r_g[2,Intercept]` by 10
   in every draw leaves the curve unchanged to five digits, and the two
   levels give the same curve: the method draws a NEW level's effect
   per draw even where the condition names an observed one. The fit
   method uses the level (its `re_formula = NULL` curve minus the `NA`
   one is 0.0432, `ranef(fit)` for level 1), and so does brms.
   `dev/sampfix-08-ce-level.R`, seed 9. The cause is that
   `ce_new_level_spec()` builds a new-level spec for every grouping
   variable in `na_vars`, including one `conditions` names. Silent
   wrong answer; it belongs to a conditional-effects lane.
2. **The covariance parameters in the draws keep internal names and the
   unconstrained scale** (`theta_1` where brms has `sd_g__Intercept`).
   `VarCorr()` and `hypothesis()` compute the brms quantities, so
   nothing is wrong, but `variables(ds)` and `variables(fit)` still
   differ there. Renaming needs a per-covariance-structure map (an
   `us()` block's theta is jointly the sds and the correlations); not
   done. `?variables` now states this gap instead of the ordinal one.
3. Under REML with the default priors, `laplace = TRUE` samples an
   objective rebuilt WITHOUT the REML integral over `beta`, so the
   draws are those of the ML Laplace posterior. It is not wrong, but it
   is not what "a REML fit" suggests. Not changed.
4. `?frm_sample`'s `prior` argument says list names are "parameter names
   as in the draws". Named-list priors resolve against the internal
   outer names (`tau_raw`, `beta`, ...), which already differed from the
   draws for `sigma` and now also for the thresholds. Not changed.
5. **M2 (reviewer): a harmless false alarm.** `(0 + x | g)` at `newdata`
   with `x = 0` is refused on laplace draws, because the probe's `NA`
   times 0 is `NA`, although the prediction there does not depend on
   `b`. The user is told to sample without `laplace = TRUE`, which
   works. `dev/sampfix-rev-01-probe.R`, rerun on the final build in
   `dev/sampfix-log/11-probe-lane.txt`.
6. **M7 (reviewer), pre-existing:** `conditional_effects(re_formula =
   NULL)` on draws draws a NEW level from `N(0, Sigma(theta))` (core
   `ce_draw_new_levels()`), where brms 2.23.0 calls its predictions with
   `allow_new_levels = TRUE` and the default `sample_new_levels =
   "uncertainty"`, which draws from the existing levels' effects. Full
   draws as well as laplace ones (`dev/sampfix-rev-08-brmsce.R`). Item 1
   above is the same code path seen from `conditions`.
7. **R2 (reviewer), a NaN the watch can still pass.** In
   `bf(yn ~ log(c1) * x + exp(a)^k, c1 ~ 1, a ~ 1 + (1 | g), k ~ 1,
   nl = TRUE)`, with `c1 = -1` AND `k = 0` on draw 1 only, draw 1 is NaN
   on every row for its own reason, so the probe sees NaN at both fills
   and the watch takes the all-`NA` pattern as its reference; draws 2 to
   6, which read `a`'s group effects, match it. `posterior_epred()`
   returns 480 of 480 cells non-finite where the full draws give 80.
   It needs two exact parameter values at one draw and the result is
   NaN, never a finite wrong number. Not fixed; comparing each draw with
   itself at fill 0 would close it at the cost of a second evaluation
   per draw. `dev/sampfix-rev-r2-01-watch.R`, section 1a.

## Last round (re-check R1 to R3)

- **R1, fixed.** `draws_laplace_watch()` spent nearly all its time in
  `unlist(v)` building one name per cell. It now uses
  `unlist(v, use.names = FALSE)`. Measured with the reviewer's
  `dev/sampfix-rev-r2-03-cost.R` (`y ~ x + (1 | g)`, 20000 rows in 200
  groups, 1000 draws, data seed 5, draws seed 1, `re_formula = NA`,
  3 repetitions, medians), laplace over full draws:

  | run | `posterior_predict()` | `posterior_epred()` | watch, 1000 calls |
  |---|---|---|---|
  | before (`13-cost-before.txt`) | 1.964 | 2.000 | 3.670 s |
  | after (`13-cost-after.txt`) | 1.037 | 1.050 | 0.160 s |
  | after, again (`13-cost-after2.txt`) | 1.067 | 1.024 | 0.130 s |

  The machine was loaded during the before run (full-draws
  `posterior_predict()` 5.56 s against 3.28 s after), so read the ratios
  and the watch-alone time, not the absolute seconds; the reviewer's
  own before run gave ratios 1.946 and 1.794 and 2.160 s. Results are
  `identical()` between laplace and full draws in every run.
- **R3, fixed.** The watch compared `!is.finite()` patterns, so a later
  draw that overflowed to `Inf` on its own was refused as a read. It
  compares `is.na()` patterns now. Checked first on every watched path
  (`dev/sampfix-12-fillpaths.R`, `dev/sampfix-log/12-fillpaths.txt`,
  data seed 77, draws seed 1): with the probe off and the watch replaced
  by a recorder, every read of the `NA` fill came out `NA` or NaN and
  none came out `Inf`, on 38 combinations. The paths were
  `posterior_epred`, `posterior_linpred` and `posterior_predict`'s
  dpars, in sample and at `newdata`, on gaussian, poisson, bernoulli,
  lognormal with `(1 | g)` in mu and sigma, cumulative, a smooth, and a
  nonlinear `exp(a)^k` body; plus `conditional_effects()` at its default
  and at `re_formula = NULL` on a smooth plus `(1 | g)`, the `mi()` model,
  `pp_mixture()` and `posterior_epred()` on a grouped mixture. The
  `hypothesis()` path had no read to observe (an `sd_` needs none). A
  test (`test-laplace-draws.R`: lognormal, `b_Intercept = 720` on draw 3)
  was seen to fail on the build before the change
  (`dev/sampfix-log/r3-laplace-before.txt`) and passes now.

## Version floor

frmtmb.sample needs no new core EXPORT. It reads a new family field,
`post$ord_thresholds_raw`, through the existing `brms_fixef_rows()`.
Settled at consolidation (M5): frmtmb.sample's floor rises to the next
frmtmb anyway, because other lanes need new core exports, so the draws'
names never depend on which core is installed; NEWS says so. The
fallback to the internal names stays in the code for a family that
declares no inverse. Bumps: minor for frmtmb.sample (a visible rename
of draw columns), patch for frmtmb.

## Nits round (after review, dev/reviews/2026-09-29-sampfix.md)

- **S1.** A laplace refusal suggests `re_formula = NA` only when the
  function the user called takes it and the same call at
  `re_formula = NA` passes the probe on this model
  (`draws_laplace_refuse()`); otherwise it says only "Sample without
  laplace = TRUE". Rerun of the reviewer's message script on the final
  build (`dev/sampfix-log/11-messages-lane.txt`): `bayes_R2()`, the
  smooth-only model and `pp_mixture()` carry no hint now.
- **S2.** Every refusal names the function the user called:
  `draws_as_caller()` records the outermost caller, and `fitted()`,
  `predict()`, `residuals()`, `predictive_error()`,
  `predictive_interval()`, `bayes_R2()`, `pp_check()`, `loo()`,
  `waic()`, `psis()`, `loo_compare()`, `hypothesis(scope =)` and
  `coef()` set it. In the rerun every refusal names the called function;
  the reviewer's check flags the second `pp_check()` line only because
  that message begins `pp_check(type = 'loo_pit_overlay')`, which names
  it with its type.
- **S3.** `?frm_sample` leads with the user's reason and no longer cites
  a dev/ script.
- **M1.** `draws_laplace_watch()`: on laplace draws, a draw whose result
  has non-finite cells where the first draw's has none is refused with
  the same message. Compared with "any non-finite output", this keeps a
  row that is `NA` at every draw (a missing covariate) from being taken
  for a read. The reviewer's construction (`c0 + exp(a)^k`, `k = 0` on
  draw 1) is a test in `test-laplace-draws.R`, seen to fail on the build
  before this round (`dev/sampfix-log/nits-laplace-prenits.txt`, the
  failure at line 233) and passing now; the probe rerun no longer lists
  it among its FAILs.
- **r_eff.** `test-ppcheck-loo.R` pins brms's rule: the chain structure
  when every draw is used, one chain for a subset, with the two rules
  shown to differ on the same draws.
- **Found by R CMD check in this round:** the rewrite of the probe block
  dropped `is_na_re_form()`, which `conditional_effects()` calls at
  `re_formula = NULL` with a group-level term. The probe then failed at
  both fills, took that for an error unrelated to the fill, and let the
  call through, and a smooth's curves came back `NA` without a word.
  Restored, and a test (`test-laplace-draws.R`, a curve at
  `re_formula = NULL` that reads a smooth) was seen to fail on the build
  without it (`dev/sampfix-log/nits-laplace-noisna.txt`). The probe
  itself was the hole: an error at both fills let the call through. It
  now re-raises that error, which the call's first draw would raise
  anyway, so a probe that fails on its own account cannot wave a call
  through again.
- **S4.** The corrections above (26, not 28; within 1 ulp);
  `dev/sampfix-log/05-lane.txt` regenerated on the final build (all
  four types, all draws and the subset: max |lw diff| 0, max |k diff| 0).
