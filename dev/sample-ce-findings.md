# SAMPLE-CE lane findings

Worktree `C:/Users/adf44/source/r/frmtmb-wt-sample-ce`, branch
`wt-sample-ce`, base 68d6782 (core 0.52.0, `frmtmb.sample` 0.2.0).
Written incrementally; never commit.

## Item 1: the conditional_effects() draws frame

### The gap, as the review records it

`dev/reviews/2026-09-05-conditional-effects.md` sections D-E and P8 name
seven line numbers in
`extensions/frmtmb.sample/R/conditional-effects-draws.R`. Confirmed by
reading the file at 68d6782:

| line | what it does | what core's fit method does |
| --- | --- | --- |
| `:41` | `ce_cats_display(rspec, dpar)` | `ce_display_kind(rspec, dpar, categorical)`; `categorical =` is not even a formal of the draws method |
| `:45` | `ce_grids_build(..., conditions, data)` | `..., int_conditions, na_vars`; so `int_conditions` is unimplemented and `re_formula = NULL` takes the first observed level |
| `:81`, `:91` | `d <- g$nd[g$ev]` | `ce_frame(nd, ev, g$v2, cond, cats)` |
| `:87` | `cats__` appended after the band columns | `ce_frame()` puts it before `effect1__` |
| `:98` | `cond__` only when `length(cond_sets) > 1L`, as character | always, as a factor over every condition level |
| `:102` | `ce_finalize(...)` with `cats_key` omitted | `cats_key = categorical`, so the ordinal key is `"x:cats__"` |

So the returned columns are the varied predictor(s), `estimate__`,
`se__`, `lower__`, `upper__` and sometimes `cats__`/`cond__`. Core's
0.52.0 frame carries, in order: varied predictor(s), every other model
variable at its held value, `cond__`, `cats__` (ordinal only),
`effect1__`, `effect2__` (two-way or ordinal), then the band columns.

### What has to be exported from core

`R/sampling-api.R` already exports `ce_grids_build`, `ce_boot_one`,
`ce_finalize`, `ce_cats_display`, `ce_structure_check`, `ce_re_formula`.
Missing for parity, all in `R/conditional-effects.R`:

- `ce_frame` (`:309`) - the column set itself
- `ce_group_vars` (`:437`) - the `na_vars` a `re_formula = NULL` call blanks
- `ce_display_kind` (`:457`) - `categorical =` and the `cats_mean` display
- `ce_new_level_spec` (`:672`), `ce_boot_grids` (`:702`),
  `ce_draw_new_levels` (`:725`) - a NEW group's effects, drawn per draw
  from THAT draw's theta, which is what makes `re_formula = NULL` mean
  a new level on the draws side too

`ce_draw_new_levels(f, nspec)` reads `f$estimates$b` and
`f$estimates$theta` only, so `draws_fit_at()`'s per-draw fit feeds it
unchanged. On the draws side the new group is drawn per POSTERIOR DRAW
rather than per bootstrap replicate, which is brms's own construction.

### What landed

`R/sampling-api.R` gains six exports in the conditional-effects seam:
`ce_frame`, `ce_display_kind`, `ce_group_vars`, `ce_new_level_spec`,
`ce_boot_grids`, `ce_draw_new_levels`.
`extensions/frmtmb.sample/R/conditional-effects-draws.R` is rewritten
against them: the frame comes from `ce_frame()`, `cond__` is always a
factor over every condition level, the ordinal key is `cats_key`, and
the grid takes `int_conditions` and `na_vars` like the fit method.
`categorical =` and `seed =` are new formals; unknown dots now warn
instead of vanishing.

`re_formula = NULL` draws the new group per POSTERIOR DRAW from that
draw's own theta (`ce_draw_new_levels()`), where the fit method draws
it per bootstrap replicate. That is brms's construction and it is why
the band widens instead of reproducing the population band.

### The column and grid comparison (sm-ce-parity.R)

Every case below: `identical(names(draws_frame), names(fit_frame))` is
TRUE, and every non-band column (the whole grid: varied predictors,
held values, `cond__`, `effect1__`, `effect2__`, `cats__`) is
`all.equal()`-identical between the two.

| case | columns, in order | key | grid identical |
| --- | --- | --- | --- |
| `effects = "x"`, `re_formula = NA` | x, y, z, g, cond__, effect1__, estimate__, se__, lower__, upper__ | `x` | TRUE |
| default effect list | same, per effect (`x`, then `z`) | `x`, `z` | TRUE |
| `conditions = data.frame(z = c(-1, 1))` | same, cond__ with two levels | `x` | TRUE |
| `int_conditions = list(z = c(lo = -1, hi = 1))` on `x:z` | x, z, y, g, cond__, effect1__, effect2__, + band | `x:z` | TRUE; `effect2__` levels `hi,lo` on both |
| `re_formula = NULL` | same, `g` column all `NA` on both | `x` | TRUE |
| ordinal, default | x, y, cond__, cats__, effect1__, effect2__, + band | `x:cats__` on both | TRUE |
| ordinal, `categorical = FALSE` | x, y, cond__, effect1__, + band | `x` on both | TRUE |

Before this change the draws frame was `x, estimate__, se__, lower__,
upper__` in case 1 (no held values, no `cond__`, no `effect1__`), the
ordinal key was `x`, `int_conditions` was not a formal, and
`categorical = FALSE` could not be asked for.

Two numbers that are not column parity and are the point of the
`re_formula` fix, same model and seed, `resolution = 6`:

- mean band width, population curve: **1.917**
- mean band width, new level: **3.929**

The old code returned the population band under the new-level name.
And the ordinal expected-category-number display agrees with the fit
method's delta-method estimate to three decimals over the grid (fit
1.220 to 2.817, draws 1.220 to 2.818).

---

## Item 2: the Savage-Dickey evidence ratio

### What landed

`extensions/frmtmb.sample/R/evidence-ratio.R` is new;
`R/methods-draws.R:190` (`hypothesis.frmtmb_draws`) gains two columns,
`evid_ratio` and `post_prob`, and an `"evid_ratio_mcse"` attribute.
Core was not touched for this: the prior-density accessor reads the
resolved specification through the already-exported
`resolve_prior_input()`, and evaluates the density in this package.

The construction, per row:

- DIRECTIONAL (`<`, `>`): the posterior odds of the claim, counted from
  the draws. No prior enters. Always available.
- POINT (`=`): the Savage-Dickey ratio, kernel density of the drawn
  quantity at the point over the PRIOR density there.
- Written with no comparison at all (this package's summary spelling,
  which `hyp_parse_all()` also reports as `two.sided`): no ratio, no
  complaint. Nothing there to weigh.

The numerator is brms's own construction, rebuilt from
`brms:::density_ratio()`: widen the range so it contains the point,
`stats::density()` on 4096 grid points, spline read at the point,
floored at zero. Bandwidth is `density()`'s default `nrd0`, Silverman's
`0.9 min(sd, IQR / 1.34) n^(-1/5)`. Matching brms's smoother is what
makes the two numbers comparable at all.

The denominator is NOT brms's. brms estimates the prior density from
`sample_prior = "yes"` draws; this evaluates it from the specification
the model was sampled under. That is the half of the ratio where no
Monte Carlo error is necessary.

### Validation against brms (sm-brms-validate.R)

brms 2.23.0, rstan 2.32.7, StanHeaders 2.39.1. rstan needed
`R_MAKEVARS_USER` with `CXX17FLAGS = ... -std=gnu++17`: rstan emits
`Rcpp::plugins(cpp14)`, whose `-std=c++1y` lands after Makeconf's
CXX17STD and wins, and StanHeaders 2.39.1 needs C++17. Same finding as
the sibling lanes' `bc-`/`bm-makevars-cxx17.mk`. tmbstan is NOT broken
on this machine right now: `tmbstan_build_broken()` returns FALSE and a
short chain recovers the slope, so the frmtmb side is real sampling.

Both shapes: 4 chains, 3000 iterations, 1000 warmup, 8000 kept draws
each side, seed 31, `prior(normal(0, 1), class = "b")` on both sides.
The rest of the stack is each package's own defaults. Those were
printed side by side per model rather than taken on trust, since the
whole comparison depends on their being equal:

| slot | brms `prior_summary()` | frmtmb `get_prior(route = "sample")` |
| --- | --- | --- |
| b | `normal(0, 1)` (passed on both sides) | passed at `frm_sample(prior =)` |
| Intercept, S1 | `student_t(3, 0.4, 2.5)` | `student_t(3, 0.4, 2.5)` |
| Intercept, S2 | `student_t(3, 0.2, 2.5)` | `student_t(3, 0.2, 2.5)` |
| sigma | `student_t(3, 0, 2.5)` | `student_t(3, 0, 2.5)` |
| sd (S2) | `student_t(3, 0, 2.5)` | `student_t(3, 0, 2.5)` |

(`get_prior(route = "sample")` reports the DEFAULTS, so class `"b"`
prints `(flat)` there; the `normal(0, 1)` reached the chain through
`frm_sample(prior =)`.) The likelihood identity behind these two
shapes is `dev/brms-likelihood-tests.md`'s row 1 family and its
random-effect row 7 construction.

**S1, `y ~ x`, gaussian, n = 150.** Posteriors first: brms mean 0.20446
sd 0.07731, frmtmb mean 0.20542 sd 0.07682.

| hypothesis | brms | frmtmb | frmtmb MCSE | abs diff |
| --- | --- | --- | --- | --- |
| `x = 0` | 0.425087 | 0.468621 | 0.0550 | 0.0435 |
| `x > 0` | 295.296 | 249.000 | 54.50 | 46.3 |
| `x < 0` | 0.0033864 | 0.0040161 | 0.00103 | 0.00063 |

**S2, `y ~ x + (1 | g)`, gaussian, n = 240, 12 groups.** brms mean
0.06579 sd 0.06860, frmtmb mean 0.06514 sd 0.06813.

| hypothesis | brms | frmtmb | frmtmb MCSE | abs diff |
| --- | --- | --- | --- | --- |
| `x = 0` | 9.18031 | 9.58529 | 0.1957 | 0.4050 |
| `x > 0` | 4.99251 | 4.93032 | 0.0788 | 0.0622 |
| `x < 0` | 0.200300 | 0.202827 | 0.00331 | 0.00253 |

S1's difference is inside the combined Monte Carlo error of two
independent 8000-draw runs (0.0435 against a reported 0.0550, 0.79
sigma). **S2's is not, read that way**: 0.4050 against 0.1957 is 2.07
reported MCSE. It is covered only once the gap is DECOMPOSED, which is
what the next section does and what this sentence originally elided:
1.42% of 9.18 is 0.130 of systematic denominator error, leaving 0.275
of sampler noise against a combined MCSE of
`sqrt(0.1957^2 + 0.2373^2)` = 0.3075, i.e. 0.89 sigma. The number is
right; the claim about it was loose.

### The sharp comparison: the denominator, isolated

Independent chains confound the method with the sampler. So the same
ratio was computed on BRMS's OWN posterior draws with this package's
numerator, which makes the numerator input byte-identical and leaves
the denominator as the only difference:

| | S1 | S2 |
| --- | --- | --- |
| exact prior density `N(0, 1)` at 0 | 0.39894228 | 0.39894228 |
| brms prior-draw KDE at 0 (n = 8000) | 0.38605907 | 0.40460660 |
| brms `Evid.Ratio` | 0.42508665 | 9.18031049 |
| ours on brms's own draws | 0.41135915 | 9.31065567 |
| relative gap, denominator only | **3.229%** | **1.420%** |

So the two constructions agree exactly on the numerator, and differ by
the error in brms's kernel estimate of a density it could have
evaluated: 3.2% low on S1, 1.4% high on S2. That error does not shrink
with the posterior draw count, only with the prior draw count, and it
carries straight into the reported Bayes factor.

The directional ratios are the control, and they are EXACT on shared
draws:

| | brms | ours |
| --- | --- | --- |
| S1 greater | 295.29629630 | 295.29629630 |
| S1 less | 0.00338643 | 0.00338643 |
| S2 greater | 4.99250936 | 4.99250936 |
| S2 less | 0.20030008 | 0.20030008 |

### Monte Carlo error at the draw count used

The estimator is the between-CHAIN spread of the ratio,
`sd(per-chain) / sqrt(chains)`. A bootstrap over rows would understate
it: the draws are autocorrelated. It is slightly conservative, since a
per-chain estimate carries a wider `nrd0` bandwidth than the pooled one.

| shape | draws | chains | ratio | MCSE | chain range |
| --- | --- | --- | --- | --- | --- |
| S1, on brms draws | 8000 | 4 | 0.4114 | 0.0376 | 0.1621 |
| S1, on frmtmb draws | 8000 | 4 | 0.4686 | 0.0550 | 0.2568 |
| S2, on brms draws | 8000 | 4 | 9.3107 | 0.2373 | 0.9606 |
| S2, on frmtmb draws | 8000 | 4 | 9.5853 | 0.1957 | 0.9294 |

**Read that before believing a third digit.** At 8000 draws the point
ratio is good to roughly 5% (S1) and 2% (S2), and the S1 chain range
of 0.16 on a ratio of 0.41 is 39% of the value. This is why the number
is exposed as `attr(h, "evid_ratio_mcse")` rather than left implicit.

### What refuses, and why

Each refusal names the parameter and its own reason, leaves that row's
`evid_ratio` at `NA`, and does not disturb the other rows.

| case | reason reported |
| --- | --- |
| flat class `"b"` (the DEFAULT, as in brms) | `x has no proper prior (frm_sample() leaves class "b" flat, as brms does); write one with set_prior(class = "b") and resample` |
| `sd_g__Intercept = 0` | not a population-level coefficient; its prior sits on the log standard deviation, and a point null at zero is on the boundary |
| `x^2 = 0` | not an affine function of the coefficients it names |
| `Intercept = 0` | the prior is about the intercept at the predictor MEANS (the entry carries an `offset` over the other coefficients), so its density is not the density of the tested quantity |

What is accepted beyond a bare coefficient: an affine map of ONE
coefficient (`2 * x = 0` divides the density by `|b|`), and an affine
combination of SEVERAL when every one carries a normal prior, since the
convolution is then normal (`x - w = 0` under two independent
`normal(0, 1)` priors is `normal(0, sqrt(2))`; asserted in
`test-evidence-ratio.R`).

### Left out, deliberately

- No prior-DRAW route. A general nonlinear hypothesis could get a
  denominator by sampling the prior and smoothing it, which is what
  brms does. It would put back exactly the Monte Carlo error the exact
  denominator removes, and it needs an RNG per density kind that core
  does not have. Refusing by name is the better answer.
- `evid_ratio` is not added to the FIT method's `hypothesis()`. A Wald
  or profile interval has no posterior density to put over a prior one.
- `extensions/frmtmb.sample/R/loo.R:665` says `hypothesis()` "gives the
  posterior probability of a directional claim". That is now
  incomplete: it also gives a point Bayes factor. The file belongs to
  the wt-protocol lane this round and was not touched. One sentence for
  whoever holds it next.

---

## Verification

One process per test file, counts audited by name. The runner puts the
package NAMESPACE in the test environment's parent so an unqualified
internal resolves the way `devtools::test()` would.

### frmtmb.sample suite: 949 pass, 0 fail, 2 skips (12 files)

| file | pass | fail | error | skip |
| --- | --- | --- | --- | --- |
| test-conditional-effects-draws.R | 37 | 0 | 0 | 0 |
| test-draws-methods.R | 97 | 0 | 0 | 0 |
| test-draws-spellings.R | 48 | 0 | 0 | 0 |
| test-evidence-ratio.R | 24 | 0 | 0 | 0 |
| test-loo.R | 70 | 0 | 0 | 1 |
| test-message-uniqueness.R | 6 | 0 | 0 | 0 |
| test-parallel-chains.R | 4 | 0 | 0 | 0 |
| test-prior-route.R | 9 | 0 | 0 | 0 |
| test-reparam.R | 260 | 0 | 0 | 0 |
| test-sample-direct.R | 136 | 0 | 0 | 0 |
| test-sampling-ported.R | 208 | 0 | 0 | 1 |
| test-simulators.R | 50 | 0 | 0 | 0 |

888 at 0.2.0 with 2 skips, so **+61**, which is exactly the two new
files (37 + 24). No pre-existing file's count moved.

Two entries needed a harness correction, both verified:

- `test-parallel-chains.R:54` calls `testthat::local_mocked_bindings()`
  with no `.package =`, which needs a pkgload-loaded namespace. Under a
  plain `test_file()` it errors; run with `pkgload::load_all()` the file
  is 4 pass, 0 fail, 0 error. Nothing to do with this lane.
- `test-simulators.R` crashed once with a `frmtmb_register_compat()`
  `.onLoad` error, a race against the `R CMD check` running in the same
  private library at that moment. Run alone: 50 pass, 0 fail.

### frmtmb core suite: 6535 pass, 0 fail, 92 skips (117 files)

117 tracked `tests/testthat/test-*.R` at 68d6782 and none added by this
lane. (The brief said 116; `git ls-files` says 117.) The named files
the brief asks for:

| file | pass | fail | error | skip |
| --- | --- | --- | --- | --- |
| test-ce-bands.R | 167 | 0 | 0 | 0 |
| test-ce-facets.R | 92 | 0 | 0 | 0 |
| test-bracket-access.R | 8 | 0 | 0 | 0 |
| test-message-uniqueness.R | 6 | 0 | 0 | 0 |

(core has no `test-conditional-effects*.R`; `test-ce-bands.R` and
`test-ce-facets.R` are that suite.)

Two errors and one warning in the totals, neither from this lane:

- `test-influence-plot.R` 13 pass / 2 error under the plain harness, on
  the same `local_mocked_bindings()`-without-`.package` limitation.
  Under `pkgload::load_all()`: **21 pass, 0 fail, 0 error**.
- `test-prior-compat.R` carries one warning in the test "coef and group
  narrow the classes that read them". Pre-existing; this lane's only
  core change is roxygen text and an export list.

### R CMD check --as-cran

`_R_CHECK_CRAN_INCOMING_=false`, pandoc on PATH from
`/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools`.

- `frmtmb.sample`: **Status: OK** over 52 checks. No NOTE, no WARNING.
  Vignettes rebuilt, tests run (45 s). Re-run on the final source
  after the last comment edit: OK again.
- `frmtmb` core: the first run stopped at `checking package
  dependencies ... ERROR`, "Packages suggested but not available:
  'brokenstick', 'frmtmb.spline'". Both are environment gaps, not code:
  `frmtmb.spline` was installed into the private library from
  `extensions/`, and `brokenstick` is a CRAN Suggests that is not on
  this machine at all, so the rerun sets
  `_R_CHECK_FORCE_SUGGESTS_=false`. With those two settled:
  **Status: OK** over 56 checks, no NOTE and no WARNING, tests
  run and vignettes rebuilt (408 s).

### Toolchain note

`R CMD check` and every brms run needed
`R_MAKEVARS_USER` pointing at a file containing

    CXX17FLAGS = -O2 -Wall -mfpmath=sse -msse2 -mstackrealign -std=gnu++17

rstan 2.32.7 emits `// [[Rcpp::plugins(cpp14)]]`, whose `-std=c++1y`
lands after Makeconf's `CXX17STD`, and the last `-std` wins, so a Stan
unit compiles as C++14 while StanHeaders 2.39.1 needs C++17. Without
it `brm()` dies with a bare `Error in sink(type = "output")`. Same
finding as the sibling lanes' `bc-`/`bm-makevars-cxx17.mk`; it is
recorded under the toolchain heading of `dev/brms-likelihood-tests.md`.

## Touched files

Core:

- `R/sampling-api.R` - the conditional-effects seam section, six
  `@aliases`, six names in the `@rawNamespace export()` list
- `NAMESPACE`, `man/frmtmb-sampling-api.Rd` - roxygen output
- `NEWS.md` - a "(development version)" heading and one bullet
- `vignettes/brms-migration.Rmd` - one list item under "When you still
  want brms"

`extensions/frmtmb.sample`:

- `R/conditional-effects-draws.R` - rewritten against core's engine
- `R/evidence-ratio.R` - NEW
- `R/methods-draws.R` - `hypothesis.frmtmb_draws()` only
- `tests/testthat/test-conditional-effects-draws.R` - NEW
- `tests/testthat/test-evidence-ratio.R` - NEW
- `vignettes/sampling.Rmd` - two new sections
- `NEWS.md` - a "(development version)" heading and five bullets
- `man/frmtmb.sample-package.Rd` - roxygen picked up a stale URL that
  `DESCRIPTION` already carried; not a change of this lane's making,
  and leaving it would make `roxygenise()` non-idempotent

Not touched, on purpose: `R/priors.R` and `frm_sample()`'s prior
defaults (wt-priors2), `extensions/frmtmb.sample/R/loo.R`
(wt-protocol).


---

# Punch round, 2026-09-06

Against `dev/reviews/2026-09-06-sample-ce.md`. Nine items; all nine
addressed, none disputed. Numbers below were re-measured after the
machine restart, on the same private library.

## 1. BLOCKER. The draws method plotted mu, not the expected response

The reviewer is right, and the reproduction stands. `ce_boot_one()`
branches on `dpar = NULL` to predict the expected response, and the
draws method passed its resolved `dpar` (`"mu"`) unconditionally where
core's fit method passes `pred_dpar`.

**The fix keeps ONE rule rather than two copies.** `ce_pred_dpar()` is
new at `R/conditional-effects.R:530`; core's own fit method now calls
it at `:1734` and derives `mean_display <- is.null(pred_dpar)` at
`:1738`, which is what that variable always was, since `dpar` is
resolved by then and can never be `NULL`. The draws method calls it at
`extensions/frmtmb.sample/R/conditional-effects-draws.R:96` and passes
the result to `ce_boot_one()` at `:145` and `:153`. `mean_is_mu()`
stays private, as the review asked.

Measured, five-point grid, 2 chains x 800 iterations, `draws / fit`:

| family | before (reviewer) | after |
| --- | --- | --- |
| zero-inflated poisson | 1.3527 to 1.3692 | 0.9985 to 1.0181 |
| hurdle poisson | sign-changing | 0.9998 to 1.0070 |
| `trunc(lb = 0)` gaussian | **-0.9919** at point 1 | 0.9985 to 1.0097, none negative |
| mixture(gaussian, gaussian) | one component's mu | 0.9959 to 1.0047 |

Largest relative deviation from the fit curve after the fix: 1.4%
(zero-inflated), 1.0% (truncated), 0.5% (mixture). Before: 21.7% on the
zero-inflated fit measured the same way. `dpar = "theta1"` on the
mixture, which was already correct, still agrees to 0.0009.

## 2. Tests that can see it

`extensions/frmtmb.sample/tests/testthat/test-conditional-effects-draws.R`
gains five blocks (`:235`, `:253`, `:269`, `:287`, `:310`) and a
three-family fixture at `:191`. Each family carries a gated
`expect_lt(rel_gap(...), 0.12)` against the fit method AND at least one
assertion that needs no chain agreement at all:

| family | ungated assertion | holds under old code |
| --- | --- | --- |
| zero-inflated | `all(estimate__ < dpar = "mu" display)` | **FALSE** (they were the same curve) |
| `trunc(lb = 0)` | `all(estimate__ > 0)`, `all(lower__ > 0)` | **FALSE** (-0.6179 at point 1) |
| mixture | `min(abs(estimate__ - mu1 display)) > 0.5` | **FALSE** (identical) |

Demonstrated rather than asserted: the old default display IS the
`dpar = "mu"` display, so asking for that explicitly reproduces the
defect exactly (`scratchpad/sm-oldcode.R`). Measured on the truncated
fixture, new default `0.6498 1.0815 1.8644 2.9478 4.1233` against old
`-0.6179 0.5673 1.7525 2.9377 4.1229`. The gated tolerance separates
cleanly too: relative gap to the fit curve is 0.0136 now and 0.2167
under the old path, against a 0.12 gate.

`ce_pred_dpar()` itself is pinned at `:235`, including that a named
`dpar` and a category display are left alone.

## 3. The affine probe now runs both signs

`extensions/frmtmb.sample/R/evidence-ratio.R:195-202`. The fitting
points were `0` and the unit vectors, all non-negative, so `abs()` was
indistinguishable from the identity on every one of them and on the
single non-negative check point. `abs(x) = 0` reported a Bayes factor
exactly 2x too large, because the prior of `|X|` for `X ~ N(0, 1)` is a
half-normal with twice the density at zero. The loop now probes `chk`
and `-chk`; cost is one extra `hyp_eval()` per hypothesis.
`test-evidence-ratio.R:146` asserts the refusal, and asserts that four
legitimate affine hypotheses (`x = 0`, `2 * x = 0`, `x - w = 0`,
`-x = 0`) still pass without warning.

## 4. The MCSE has at least four blocks

`extensions/frmtmb.sample/R/evidence-ratio.R:314`, `nc <- max(nc, 4L)`. Two chains gave `sd()` of two numbers over `sqrt(2)`,
which the reviewer measured understating a real 0.46 gap as 0.18.
Splitting each chain in half keeps blocks consecutive, so the
within-block autocorrelation the estimator exists to respect is
preserved; it is the split-chain construction the sampler diagnostics
already use. More than four chains keep their own boundaries.
`test-evidence-ratio.R:161` asserts four blocks on the 2-chain fixture,
and `:80` now pins the whole estimator by recomputing it from the same
blocks rather than asserting `mcse < evid_ratio`, which any MCSE under
400% of the value satisfied.

## 5. The migration vignette no longer overpromises

`vignettes/brms-migration.Rmd:468-479`. "for any parameter that carries
a proper prior" is replaced by the truthful scope, and the list of what
refuses is named inline so the page a brms user reads to decide
matches what `sampling.Rmd` already said.

## 6. The brms agreement is a test

`test-evidence-ratio.R:176`, behind `skip_unless_brms()`, no Stan.
`identical(er_kde_at(v), density_ratio(v, point = 0))` on three fixed
vectors: one inside the drawn range, one where the range-widening
branch fires, and one read at a point inside the range. Reached through
`utils::getFromNamespace()` rather than `:::` so as-cran stays quiet.
All three `identical()`, not "close".

## 7. allow_new_levels is accepted

`conditional-effects-draws.R:64` now calls core's `ce_dots()` in place
of a hand-rolled leftover check, and `:109` folds its result into
`anl`. `ce_dots()` is the eighth seam export
(`R/sampling-api.R:151-157, 263, 279`); exporting it rather than adding
two names to a `setdiff()` is the same argument as `ce_pred_dpar()` -
the set of accepted arguments is a rule, and a rule written twice
drifts. `test-conditional-effects-draws.R:310` asserts both spellings
warn-free and that a genuinely unknown argument still warns.

## 8. The sigma refusal says what sigma is

`evidence-ratio.R:214-232`. Every name reaching that branch resolved in
the hypothesis environment, so a name that is not a coefficient is one
of the natural-scale summaries that environment also carries. One
message now covers all of them and gives the file's own roxygen reason
(the prior sits on the internal parameter with a change of variables in
between) instead of "not a population-level coefficient", which was
true of the name and false of the model.

## 9. The trailing slash

`extensions/frmtmb.sample/DESCRIPTION:19` gains it, and
`man/frmtmb.sample-package.Rd:14` follows from re-roxygenizing. This
clears the as-cran URL WARNING the reviewer saw with
`_R_CHECK_CRAN_INCOMING_` on; the first round ran with it off and so
reported OK, which the review is right to qualify.

## Nothing disputed

All nine are real. Two review observations outside the punch list were
also acted on, since they were cheap and correct:

- `dev/sample-ce-findings.md`'s "Every difference is inside the combined
  Monte Carlo error" was looser than its own data on S2 (2.07 reported
  MCSE). The sentence now decomposes the gap, as the review says it
  should. The number is unchanged.
- `test-conditional-effects-draws.R`'s ordinal assertion
  `all(estimate__ > 1 & < 3)` held for any probability vector. It is
  replaced at `:156` by an equality against the per-category display of
  the SAME draws, which is exact because the map is linear.

## Verification, punch round

Same runner as before: one process per test file, test environment
parented on the package namespace.

### frmtmb.sample: 978 pass, 0 fail, 2 skips (12 files)

| file | pass | fail | error | skip |
| --- | --- | --- | --- | --- |
| test-conditional-effects-draws.R | **57** | 0 | 0 | 0 |
| test-draws-methods.R | 97 | 0 | 0 | 0 |
| test-draws-spellings.R | 48 | 0 | 0 | 0 |
| test-evidence-ratio.R | **33** | 0 | 0 | 0 |
| test-loo.R | 70 | 0 | 0 | 1 |
| test-message-uniqueness.R | 6 | 0 | 0 | 0 |
| test-parallel-chains.R | 4 | 0 | 0 | 0 |
| test-prior-route.R | 9 | 0 | 0 | 0 |
| test-reparam.R | 260 | 0 | 0 | 0 |
| test-sample-direct.R | 136 | 0 | 0 | 0 |
| test-sampling-ported.R | 208 | 0 | 0 | 1 |
| test-simulators.R | 50 | 0 | 0 | 0 |

949 before this round, so **+29**: +20 in the conditional-effects file
(37 to 57) and +9 in the evidence-ratio file (24 to 33). No other
file's count moved. `test-parallel-chains.R` again reports 3 pass /
1 error under a plain `test_file()` and **4 pass / 0 error** under
`pkgload::load_all()`, the `local_mocked_bindings(.package =)`
limitation section (f) of the review documents; it is not a regression
and not this lane's.

The new brms assertion is NOT skipped here: brms is installed, so
`er_kde_at()` was compared against `brms:::density_ratio()` on three
fixed vectors and `identical()` on all three.

### frmtmb core: 6537 pass, 0 fail, 91 skips (117 files)

`R/conditional-effects.R` changed this round, so the whole core suite
ran rather than the touched files alone. `test-ce-bands.R` 167,
`test-ce-facets.R` 92, `test-bracket-access.R` 8,
`test-message-uniqueness.R` 6, all clean.

Two deltas from the first round's 6535/92, both explained and neither
this lane's:

- `test-influence-plot.R` 13 pass / 2 error under the plain harness,
  21 pass / 0 error under `pkgload::load_all()`. Same limitation as
  above.
- `test-ps.R` moved from 49 pass / 1 skip to 51 pass / 0 skip. Its only
  conditional skip is `skip_if_not_installed("brokenstick")` at `:249`,
  and `brokenstick` is now present on this machine where it was absent
  before the restart. Nothing in this lane touches it.

`test-prior-compat.R` still carries its one pre-existing warning.

### Roxygen

`roxygenise()` on both packages, twice: the second pass writes nothing
for either. Idempotent.

### R CMD check --as-cran, frmtmb.sample

`_R_CHECK_CRAN_INCOMING_=false`, pandoc on PATH: **Status: OK** over 52
checks, no NOTE and no WARNING, tests run and vignettes rebuilt.

Item 9 needed the other setting to verify, so a second run with
`_R_CHECK_CRAN_INCOMING_=TRUE` and network available (no tests, no
vignettes; those are covered above). The URL block the review reported

```
Found the following (possibly) invalid URLs:
  URL: https://aforren1.github.io/frmtmb/frmtmb.sample
       (moved to https://aforren1.github.io/frmtmb/frmtmb.sample/)
    From: DESCRIPTION
          man/frmtmb.sample-package.Rd
    Status: 301
```

is **gone**: `grep -ci "invalid url|moved to https"` over the check log
returns 0. What remains under incoming feasibility is one WARNING that
no edit can clear, and that is not this lane's:

```
New submission
Strong dependencies not in the CRAN or BioC software repositories:
  frmtmb
Suggests or Enhances not in mainstream repositories:
  frmtmb.latent
```

An extension in a monorepo whose core is not on CRAN reports that
whatever it does.
