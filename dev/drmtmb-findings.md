# drmTMB lane: verified surface, agreement, and what to take from it

**STATUS 2026-09-22: reviewed MERGEABLE (punch round 1 recheck), held
uncommitted in the worktree by the user's decision until lanes
`wt-reunc`, `wt-correct` and `wt-skewinit` land, so all four merge
together. Two nits fixed after the recheck: the over-long README line
reflowed, and the upstream issue's last line wrapped in `as.numeric()`
so its printed output matches its comment.**

Lane `wt-drmtmb`, 2026-09-22. frmtmb 0.61.0 (base 2d70b33, no `R/`
change), drmTMB 0.7.0 (CRAN Windows binary), RTMB 2.0, TMB 1.9.25,
R 4.6.1. drmTMB is in the private library
`C:/Users/adf44/source/r/drmtmb-lib`; frmtmb is the release build in
`rellib-r3`. Every number below comes from a script in `dev/drmtmb-*.R`
and its log in `dev/drmtmb-log/`, named beside it.

## Answer

1. **The two packages fit the same likelihood wherever both fit a
   model.** 27 comparisons agree, from gaussian location-scale to
   Tweedie, animal-model and `mi()` fits. The log-likelihoods
   agree to within 2.8e-8. The two objectives, evaluated at one shared
   point, agree to within 3.9e-11 of |logLik|. Over those 27 rows,
   estimates agree to within 1.5e-4 standard errors, and standard
   errors to within 9.9e-5 relative, the worst on the near-boundary
   bivariate row. The test file's own 24 comparisons are a subset on
   the same seeds and are tighter on standard errors, 1.8e-5; each
   figure names the set it covers, and `dev/drmtmb-test-margins.R`
   prints the test file's. Every parameter map is written out (table
   below).
2. **One frmtmb defect.** `skew_normal()` can stop at alpha = 0, a
   stationary point, with convergence code 0 and no warning. In the
   reproduction it is 24.7 log-likelihood units below the maximum. On
   a designed case it happens in 28 of 40 seeds. Root cause: the alpha
   start takes its sign from the skewness of the raw response, not of
   the residuals. Starting from the residual skewness is NOT sufficient
   to remove it; a two-start refit was, on 80 of 80 fits. Not fixed
   here, because this lane changes no `R/`.
3. **One drmTMB inconsistency, drafted as an upstream issue and not
   filed.** Under `REML = TRUE`, adding a random effect to `sigma`
   also integrates `beta_sigma` out. The REML criterion therefore
   changes with the model, and a model reports a restricted
   log-likelihood 4.7 to 5.5 units below the model it contains (10 of
   10 seeds). Attributed exactly (to 10 digits) by rebuilding drmTMB's
   own template.
4. **Every other difference is a parameterization**, and each is now a
   stated map: drmTMB caps every correlation at 0.999999, rescales a
   phylogeny to unit height, uses `sigma = 1 / sqrt(precision)` for
   beta, beta-binomial, NB2 and Gamma, `sigma^2 = phi` for Tweedie, and
   `nu = 2 + exp(eta)` for Student-t.
5. **Part of drmTMB's website surface is not in the CRAN release.**
   `heritability()`, `icc()`, `repeatability()`, `chibar_pvalue()`,
   `lrt_boundary()` and `objective_at()` are absent from drmTMB 0.7.0,
   and its `anova()` refuses every call. They are exported on GitHub
   `main` (0.7.1, commit 54df129, 2026-09-21). I read that source; I
   did not run it.
6. The test file passes with its gate set (13 tests, 128 expectations)
   and skips all 13 in each of the other three states: gate unset,
   drmTMB absent, `NOT_CRAN` unset. It is a fit tier behind
   `FRMTMB_DRMTMB_FIT_TESTS`, like the brms fit tier, because
   otherwise every check job on every platform would install drmTMB and
   run about 50 fits. Four mutated parameter maps each make it fail.
   README, CONTRIBUTING and the submission draft now include drmTMB,
   with only claims that a run in this lane supports.

## How the agreement is measured

`dev/drmtmb-agree-lib.R` holds the engine. For each model it takes
both fits and a map from drmTMB's outer parameter vector to frmtmb's,
one row per parameter: the function, its inverse and its derivative.
It then reports:

- the two log-likelihoods at their own optima;
- frmtmb's objective minus drmTMB's, at drmTMB's optimum and at
  frmtmb's optimum (both `obj$fn`, so this is the Laplace or REML
  objective of each package at one point);
- the same gap one standard error away from drmTMB's optimum in every
  coordinate, with a fixed alternating sign pattern rather than a draw,
  so that the engine leaves the RNG stream alone;
- estimate differences divided by frmtmb's standard error, and the
  ratio of standard errors after the delta-method map.

Gaps at the optima alone are weak evidence for the map. The gradient
vanishes there, so a map that is slightly wrong still passes. The
off-optimum point is the check that caught the correlation cap. With
the uncapped map `tanh` for `rho12`, the gap there was 8.8e-6 to
1.5e-5 in the four bivariate models. With the capped map it is below
4e-11.

Standard errors: frmtmb's come from `vcov(fit, full = TRUE)`, which is
identical to `RTMB::sdreport(...)$cov.fixed` on the outer vector
(`dev/drmtmb-vcovcheck.R`: maximum |log ratio| 0 on two models).
drmTMB's come from `fit$sdr$cov.fixed`.

## Agreement, per model

Scripts: `dev/drmtmb-agree.R` (seeds `sim_grouped(101)`,
`sim_pedigree(202)`, `sim_meta(303)`, tree seed 404) and
`dev/drmtmb-agree2.R` (`sim_batch2(505)`, missingness seed 506).
Table generated by `dev/drmtmb-agree-print.R` into
`dev/drmtmb-log/agree-table.md` and pasted unchanged. Gaps are
frmtmb minus drmTMB, in log-likelihood units. Group sizes: 40 groups of
10, 140 animals (20 founders, 3 generations of 40) with 3 records
each, 50 species with 4 records each, 60 studies.

| Model | logLik frmtmb | logLik drmTMB | diff | gap at drm opt | gap at frm opt | gap off opt | max diff/SE | max SE ratio - 1 |
|---|---|---|---|---|---|---|---|---|
| a gaussian mu RE, ML | -549.419427 | -549.419427 | 5.6e-09 | -1.8e-12 | -3.3e-12 | -2.7e-12 | 9.5e-05 | 1.1e-05 |
| a gaussian mu RE, REML | -553.082193 | -553.082193 | -2.6e-12 | -2.0e-12 | -2.6e-12 | -2.8e-12 | 1.1e-06 | 4.2e-08 |
| b gaussian mu+sigma RE, ML | -543.081472 | -543.081472 | 2.8e-10 | 1.0e-09 | -6.7e-12 | -2.0e-12 | 1.6e-05 | 2.4e-06 |
| b gaussian mu+sigma RE correlated labelled, ML | -542.208327 | -542.208327 | 2.8e-08 | 2.1e-08 | 2.0e-09 | -1.2e-11 | 8.8e-05 | 8.4e-06 |
| c beta mu RE, ML | 282.847258 | 282.847258 | 6.5e-09 | 6.2e-11 | 1.7e-10 | 4.0e-11 | 8.5e-05 | 8.9e-06 |
| d nbinom2 mu RE, ML | -866.176500 | -866.176500 | 7.1e-09 | 8.1e-12 | -1.3e-11 | 3.4e-13 | 8.6e-05 | 8.4e-06 |
| d nbinom2 sigma RE, ML | -882.168472 | -882.168472 | 1.5e-09 | 2.7e-10 | -3.1e-10 | -1.3e-12 | 4.8e-05 | 7.6e-06 |
| e bivariate gaussian rho12, ML | -1102.488311 | -1102.488311 | -1.3e-09 | -1.1e-11 | -1.2e-11 | -1.7e-11 | 3.9e-05 | 1.1e-06 |
| e bivariate gaussian rho12 + correlated RE labelled, ML | -1048.911121 | -1048.911121 | 1.5e-08 | -2.5e-11 | -1.7e-11 | -6.8e-12 | 1.5e-04 | 5.7e-06 |
| e bivariate gaussian rho12 + correlated RE labelled, REML | -1056.365117 | -1056.365117 | -2.8e-09 | 2.8e-11 | -8.0e-12 | 1.6e-11 | 6.9e-05 | 2.8e-06 |
| e boundary: bivariate RE correlation at 1, ML | -1056.209243 | -1056.209243 | 1.4e-11 | 1.1e-12 | 2.5e-11 | -6.4e-12 | 5.0e-06 | 9.9e-05 |
| f animal model A, ML | -581.986715 | -581.986715 | 1.4e-08 | 4.5e-13 | -1.3e-11 | -6.8e-13 | 1.2e-04 | 2.4e-06 |
| f animal model A, REML | -584.903247 | -584.903247 | 6.9e-12 | 6.9e-12 | 5.2e-12 | 6.3e-12 | 2.0e-07 | 6.3e-09 |
| f2 phylo vs gr(cov = vcv(tree)), ML | -215.266820 | -215.266820 | 1.5e-12 | -1.1e-12 | -1.3e-12 | 1.7e-13 | 1.7e-06 | 1.9e-07 |
| f2 phylo vs gr(cov = vcv(tree, corr = TRUE)), ML | -215.266820 | -215.266820 | -1.5e-12 | -1.3e-12 | -1.5e-12 | -5.7e-14 | 1.9e-06 | 3.8e-07 |
| g student mu RE, ML | -581.220334 | -581.220334 | 8.9e-12 | 8.6e-12 | 1.5e-11 | -1.4e-11 | 5.8e-07 | 1.4e-07 |
| h cumulative logit mu RE, ML | -511.844440 | -511.844440 | 9.7e-10 | 4.5e-13 | -2.5e-12 | 6.4e-11 | 3.2e-05 | 1.8e-05 |
| i skew_normal, ML | -449.564542 | -449.560314 | -4.2e-03 | -2.0e-12 | -2.6e-12 | 4.0e-13 | 2.8e-02 | 7.7e+00 |
| j meta-analysis meta_V vs se(sigma = TRUE), ML | -34.816109 | -34.816109 | -9.9e-14 | 7.1e-14 | 1.1e-13 | 6.4e-14 | 5.6e-07 | 4.7e-07 |
| k beta_binomial mu RE, ML | -871.528684 | -871.528684 | 5.2e-09 | 3.8e-11 | -3.1e-11 | -2.3e-11 | 7.4e-05 | 7.4e-06 |
| l tweedie mu RE, ML | -666.884883 | -666.884883 | 1.1e-08 | 0.0e+00 | -3.4e-13 | -1.8e-12 | 1.2e-04 | 8.3e-06 |
| m lognormal mu RE, ML | -387.463244 | -387.463244 | 3.1e-09 | -8.6e-12 | -7.7e-12 | -7.1e-12 | 4.9e-05 | 3.6e-06 |
| m lognormal sigma RE, ML | -453.985219 | -453.985219 | -2.2e-10 | -7.3e-11 | -4.5e-12 | -1.8e-11 | 1.3e-05 | 3.4e-06 |
| n Gamma(log) mu RE, ML | -414.189105 | -414.189105 | -7.0e-10 | 2.7e-10 | 4.2e-12 | -1.0e-10 | 3.5e-05 | 1.8e-06 |
| n Gamma(log) sigma RE, ML | -500.569756 | -500.569756 | -2.3e-10 | -2.1e-10 | 1.7e-12 | -2.2e-12 | 3.0e-06 | 8.5e-07 |
| o zero-inflated nbinom2, ML | -789.969363 | -789.969363 | -1.0e-11 | 9.3e-12 | -1.6e-11 | -2.2e-11 | 3.8e-06 | 5.2e-07 |
| p sd(g) ~ w vs nl exp(lsd) * zz, ML | -511.319520 | -511.319520 | -9.9e-12 | -4.5e-12 | 5.9e-12 | -2.5e-12 | 3.1e-06 | 4.7e-08 |
| q mi() gaussian missing predictor, ML | -991.448708 | -991.448708 | 3.5e-10 | 9.1e-12 | 1.9e-11 | 3.0e-12 | 2.1e-05 | 1.1e-07 |

Rows a to h, j to q: 27 comparisons, all agreeing. Row i is the
frmtmb skew-normal defect: the objectives agree to 2.6e-12 at a
common point, so it is the same likelihood, and the optima differ.
The two f2 rows reuse one drmTMB fit; the boundary row (e) sits at
group correlation 0.9937. The meta-analysis row has a third opinion:
`metafor::rma(method = "ML")` gives logLik -34.8161092509 and
tau 0.19110194, against sigma 0.19110031 (drmTMB) and 0.19110027
(frmtmb).

Three further models agree outside the engine (`dev/drmtmb-surface2.R`):
`phylo_interaction(1 | plant:poll, tree1, tree2)` equals frmtmb's
`(1 | gr(pair, cov = kronecker(C2, C1)))` with correlation-form trees,
to 4.0e-13 in logLik. The labelled all-four block over `mu1`, `mu2`,
`sigma1`, `sigma2` gives -1033.43072358 (drmTMB) and -1033.43072871
(frmtmb); both report singular convergence, because the data have no
such structure. drmTMB's own `animal(pedigree = )` builder and the
additive matrix built in `dev/drmtmb-agree-lib.R` give logLik
difference 0 exactly.

### Parameter maps (drmTMB value d, frmtmb value f)

| Quantity | drmTMB | frmtmb | Map | How established |
|---|---|---|---|---|
| gaussian, lognormal, skew-normal location and log sigma | `beta_mu`, `beta_sigma` | `beta`, `betad` | identity | all gaps < 1e-10 |
| random-effect log SD | `log_sd_*` | `theta` | identity | same |
| beta precision | `beta_sigma`, sigma = phi^(-1/2) | `betad` (log phi) | f = -2 d | help page; off-optimum gap |
| beta-binomial precision | same as beta | `betad` (log phi) | f = -2 d | same |
| NB2 dispersion | sigma = shape^(-1/2) | `betad` (log shape) | f = -2 d | help page; glm.nb theta 1.188793 vs 1.188789 |
| NB2 or Gamma sigma random intercept SD | `log_sd_sigma` | `theta` of log shape | f = d + log 2 | off-optimum gap |
| Gamma CV | sigma = shape^(-1/2) | `betad` (log shape) | f = -2 d | off-optimum gap |
| Tweedie dispersion | sigma = phi^(1/2) | `betad` (log phi) | f = 2 d | help page; gap |
| Tweedie power | 1 + plogis(eta) | `power12`, same | identity | gap |
| Student-t nu | 2 + exp(d) | 1 + exp(f), as brms | f = log(1 + exp(d)) | gap; drmTMB cannot reach nu in (1, 2] |
| skew-normal shape | `nu`, mean parameterization | `alpha`, as brms | identity | gap 2.6e-12 |
| cumulative logit cutpoints | first raw, then log increments | `tau_raw`, same | identity (reordered) | gap |
| any correlation (`rho12`, random effects) | 0.999999 tanh(d) | f / sqrt(1 + f^2) | see `rho_map()` | boundary, below |
| phylogeny | tree scaled to unit height | raw `cov` | f = d - log(h) / 2 | f2 rows |
| random intercept SD in `sd(g) ~ w` | `beta_sd_mu` | nl: SD of `zz` and `lsd` slope | intercept to log SD, slope to slope | row p |
| meta-analysis residual SD | sigma with `meta_V()` | sigma with `se(, sigma = TRUE)` | identity | row j |

## Disagreements and their root causes

### 1. frmtmb: `skew_normal()` stalls at alpha = 0 (frmtmb defect)

**Observed.** Row i: frmtmb's logLik is 4.228e-3 below drmTMB's, while
the objectives agree to 2.6e-12 at a common point. frmtmb returned
alpha = -0.0034 with standard error 12.9; drmTMB returned 0.358.

**Root cause.** alpha = 0 is a stationary point of the skew-normal
likelihood in the mean parameterization. `fam_skew_normal()` already
knows this and starts away from zero, but it takes the sign of the
start from the skewness of the RAW response:
`2 * sign(m3) + 0.5 * m3`, with `m3` computed on `y`
(`R/families.R`, `init_dpars$alpha`). When covariates or group effects
make the raw response skewed the other way from the residuals, the
start is on the wrong side of zero. nlminb then walks to the stationary
point and reports relative convergence. Seed 101: raw skewness
-0.016, residual skewness +0.020, start -2.008. frmtmb's own profile
over alpha, from its own objective, rises monotonically from -450.30
at alpha = -1 to -449.5645 at 0 and -449.5603 at 0.36, so the stall
is not a local maximum. Started at +2, frmtmb reaches drmTMB's
optimum: both print logLik -449.560313567
(`dev/drmtmb-log/skewnormal.log`).

**Rate.** Design: `y = xs + 1.5 * (|z| - sqrt(2/pi)) + N(0, 0.3^2)`,
`xs = -3|N(0,1)|`, n = 200, seeds 1 to 40, one `set.seed(s)` per seed
and no other draw. frmtmb ends below drmTMB by more than 1e-6 in
**28 of 40 seeds, by 6.406930 to 35.137329** log-likelihood units. All
28 have the raw and residual skewness of opposite sign; all 40 seeds
do, so opposite sign is necessary in this design and not sufficient.
This is a designed case, not a field rate.

The first version of this record said 33 of 40, 6.107917 to 32.692143.
That came from `dev/drmtmb-skewnormal.R` drawing an unused `rnorm(n)`
before the covariate, which the prose did not mention, so the stated
design did not reproduce the stated rate. The dead draw is gone from
that script and the rate above is from the design as written.
`dev/drmtmb-skewstart.R` runs BOTH streams, because the second one is
where the sufficiency result below was measured:

| Stream | default start | residual-skewness start | two-start +2 and -2 |
|---|---|---|---|
| as stated | 28 of 40 stall | 0 of 40 | 0 of 40 |
| with the extra draw | 33 of 40 stall | 1 of 40 (seed 34) | 0 of 40 |

**Reproduction, frmtmb only.** `dev/drmtmb-repro-skewnormal.R`
(seed 1 of that design): default start gives logLik -284.0261427,
alpha -0.000837 (SE 5.46), convergence 0, message "relative convergence
(4)", no warning and no message; `start = list(betad = c(0, 2))` gives
-259.3389354, alpha 7.506. Gap 24.69 units. The maximum gradient at the
stalled point is 1.4e-5, so a gradient check does not catch it.

**Pinned** in the test file as a known defect, with a comment to flip
it to `expect_same_model()` when the start is fixed.

#### 1.1 The residual-skewness start is NOT sufficient (sizing note)

The obvious fix, taking frmtmb's own rule `2 * sign(m3) + 0.5 * m3` on
the residuals of a `mu`-only fit instead of on the raw response, does
not close the defect. Measured with `frm(start = )`, moving alpha alone
and leaving frmtmb's own `log(sd(y))` sigma start in place
(`dev/drmtmb-skewstart.R`):

- as stated: 0 of 40 stall. With the extra draw: **1 of 40 stalls,
  seed 34**, at logLik -276.9013 against -251.5943, a gap of 25.31,
  from a correctly signed start of 2.5524.
- At that seed, alpha starts of 0.5, 1, 2, 3, 4, 5, 8, 10 and 16 all
  reach -251.5943, and only 2.5524 falls into alpha = 0. So alpha = 0
  attracts a correctly signed start too, and the sign of the start is
  not what decides the outcome.
- **Two starts, +2 and -2, keeping the better optimum: 0 of 40 on both
  streams, 80 of 80 fits.**

One measurement detail is worth keeping, because it cost a wrong
conclusion here first: moving the sigma start at the same time changes
the answer. With `sigma` started at 0 rather than at `log(sd(y))`,
seed 34's residual start reaches the optimum, and the residual start
then looks sufficient at 0 of 40 on both streams. An in-family fix
would move alpha alone, so that is the measurement that counts.

**Route and size.** Two starts and keep the better optimum costs one
extra fit per skew-normal model, about 0.5 lane-day in
`R/families.R` plus tests. The residual-skewness start is cheaper but
is measured insufficient. **The user decides the route**; nothing is
implemented here. Acceptance test: both streams of
`dev/drmtmb-skewstart.R` at 0 of 40, plus the pinned test flipped to
`expect_same_model()`.

### 2. drmTMB: REML criterion changes with a sigma random effect

**Observed.** With `sigma ~ z + (1 | g)` under `REML = TRUE`, drmTMB's
outer vector holds only the two log SDs, while frmtmb also keeps
`betad`. The criteria differ, so there is nothing to compare.

**Root cause, exact.** drmTMB marks `beta_sigma` as random, so
integrated, when `sigma` has a random effect, and not otherwise
(`fit$obj$env$random`: A = `beta_mu u_mu`, B = `beta_mu beta_sigma
u_mu u_sigma`). `dev/drmtmb-repro-reml-sigma.R` rebuilds drmTMB's own
ML template for model A with `beta_sigma` random and evaluates it at
B's fitted `log_sd_mu`. Result: -351.2997895, equal to B's reported
REML logLik -351.2997895, with B's sigma SD at 1.6e-5.

**Consequence.** `sigma ~ z` is the SD -> 0 limit of
`sigma ~ z + (1 | g)`, but in drmTMB its REML logLik is higher by 4.67
to 5.51 units, in 10 of 10 seeds (`dev/drmtmb-reml-nesting.R`, seeds
1 to 10, 30 groups of 8, no sigma group variation in the data). ML
nests in drmTMB (0 of 10). frmtmb REML nests (0 of 10), and agrees
with drmTMB on model A to at most 2.9e-9. A REML likelihood-ratio
statistic for the sigma variance component is therefore about -10 in
drmTMB. The development version's `lrt_boundary()` refuses REML pairs
whose `mu` coefficient names differ; from its source, it compares only
those names, so it would accept this pair and clamp the negative
statistic to p = 0.5. I did not run it.

Neither criterion is "the" REML for a nonlinear scale model. The
defect is the switch: the criterion depends on whether a variance
component is present.

**Draft upstream issue, not filed:**

> Title: REML criterion changes when `sigma` gains a random effect,
> so nested REML fits do not nest
>
> Thank you for drmTMB. With `REML = TRUE`, adding `(1 | g)` to the
> `sigma` formula also integrates `beta_sigma` out, alongside
> `beta_mu`. Without a `sigma` random effect, only `beta_mu` is
> integrated. So `sigma ~ z` and `sigma ~ z + (1 | g)` are scored on
> different criteria, and the larger model can report the lower
> restricted log-likelihood even when its extra standard deviation
> collapses to zero.
>
> ```r
> library(drmTMB)
> set.seed(1)
> ng <- 30; n <- ng * 8
> g <- factor(rep(seq_len(ng), each = 8))
> x <- rnorm(n); z <- rnorm(n)
> y <- 1 + 0.5 * x + rnorm(ng, 0, 0.6)[g] + rnorm(n, 0, exp(-0.2 + 0.3 * z))
> d <- data.frame(y, x, z, g)
>
> fA <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z), data = d, REML = TRUE)
> fB <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), data = d,
>              REML = TRUE)
>
> as.numeric(logLik(fA))                     #> -346.0610
> as.numeric(logLik(fB))                     #> -351.2998   (5.24 lower)
> exp(fB$opt$par[["log_sd_sigma"]])          #> 1.578e-05   (the added SD)
>
> unique(names(fA$obj$env$par)[fA$obj$env$random])
> #> "beta_mu" "u_mu"
> unique(names(fB$obj$env$par)[fB$obj$env$random])
> #> "beta_mu" "beta_sigma" "u_mu" "u_sigma"
>
> # ML nests as expected
> mA <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z), data = d)
> mB <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), data = d)
> as.numeric(logLik(mB)) - as.numeric(logLik(mA))   #> -2.9e-08
> ```
>
> Rebuilding model A's own ML template with `beta_sigma` marked random
> reproduces model B's REML value exactly, which is how we traced it:
>
> ```r
> e <- mA$obj$env
> objJ <- TMB::MakeADFun(data = e$data, parameters = e$parList(),
>                        map = e$map,
>                        random = c("beta_mu", "beta_sigma", "u_mu"),
>                        DLL = e$DLL, silent = TRUE)
> as.numeric(-objJ$fn(fB$opt$par[["log_sd_mu"]]))   #> -351.2998 = logLik(fB)
> ```
>
> We saw the same sign in 10 of 10 seeds of this design, with the drop
> between 4.67 and 5.51. A consistent choice either way would make REML
> log-likelihoods comparable across `sigma` random effects. Tested on
> drmTMB 0.7.0, R 4.6.1, Windows.

A second, minor note for the same issue or a separate one: `bf(y ~ x
+ ar(p = 1))` fails with "argument "x" is missing, with no default"
from `stats::ar()`, where other unsupported terms get a named refusal
(`dev/drmtmb-log/surface.log`).

### 3. Correlations at the boundary: drmTMB caps at 0.999999

**Observed.** In the lane's first run, with the second response
sharing the first one's group effect exactly (`sim_v0()` in
`dev/drmtmb-boundary.R`), frmtmb's logLik was higher by 9.46e-6 and
the objectives differed by 2.80e-5 at drmTMB's optimum.

**Root cause.** drmTMB's correlation link is `0.999999 * tanh(eta)`.
Evaluated at drmTMB's optimum under candidate maps c * tanh, the gap is
2.80e-5 for c = 1, 2.55e-5 for c = 1 - 1e-7, 2.2e-9 for c = 1 - 1e-6
and -2.5e-4 for c = 1 - 1e-5. drmTMB's objective along its own
correlation coordinate plateaus at -1068.4817701448 from eta = 14,
and frmtmb reached 1 - rho = 6.7e-7, beyond the cap. 27 drmTMB R
functions contain the literal 0.999999. The same cap applies to
`rho12` (off-optimum gap 8.8e-6 to 1.5e-5 with c = 1, below 4e-11 with
the cap) and to the mu-sigma correlation (the cap-removed mutant of the
test file fails there at 6.1e-10 of |logLik|, against 2.2e-14 with the
cap). Not a defect in either package; a boundary difference of order
1e-5 log-likelihood units.

### 4. Things that looked like disagreements and were not

- drmTMB's `phylo()` rescales the tree to unit height. Against the raw
  `ape::vcv(tree)` (height 0.7615) the log SD differs by log(h) / 2 =
  -0.136 and the logLik is identical; with the correlation form it is
  the identity.
- `drmTMB()` and `drm_formula()` find a matrix or tree (`A`, `tree`)
  through the calling frame, so a plain wrapper function hides them.
  The test file evaluates its calls in the test's own frame.

## drmTMB's surface, verified by running it

Logs: `dev/drmtmb-log/surface.log`, `surface2.log`, `ns.log`,
`frm-gaps.log`, `grby.log`. "Refused" means `drmTMB()` stopped with an
error; the message is in the log.

| Capability | Result in drmTMB 0.7.0 |
|---|---|
| gaussian: RE in mu; RE in sigma; random slope in sigma; mu-sigma correlated via `|p|` | all fitted |
| gaussian: correlated `(1 + x | g)` in mu | fitted |
| Student-t: RE in mu | fitted |
| Student-t: RE in sigma; RE in nu | both refused |
| beta: RE in mu | fitted |
| beta: RE in sigma; correlated slope block; `|p|` across mu and sigma | all refused |
| NB2: RE in mu; RE in sigma | both fitted |
| NB2: RE in mu and sigma together; `|p|` across mu and sigma | both refused |
| NB2: RE in `zi` | refused |
| lognormal and Gamma: RE in mu; RE in sigma | fitted (row m, n) |
| cumulative logit: RE in mu | fitted |
| Poisson: correlated `(1 + x | g)` | fitted |
| `sd(g) ~ w`, gaussian, w constant within g | fitted; log sd_g = c0 + c1 w_g (row p) |
| `sd(g) ~ x` with x varying within g | refused, with a clear message |
| `sd(g) ~ w` for Poisson or NB2 | refused |
| `sd(g) ~ w` with REML | refused |
| REML: gaussian with RE in mu, `sigma ~ z`, RE in sigma, correlated slopes, fixed only; bivariate gaussian; binomial `(1 | g)` | fitted |
| REML: Student-t, beta, NB2, Poisson, cumulative logit | refused: "implemented for univariate/bivariate Gaussian and binomial models" |
| three responses (`mvbind` or three families) | refused: one and two responses only |
| mixed-family pair (gaussian with Poisson, gaussian with NB2) | refused |
| nonlinear formula (`nl = TRUE`) | refused; `y ~ exp(x)` fits as a fixed transform |
| custom family (`quasi()`, a user family object) | refused |
| `s()` smooth | refused, pointing to mgcv and gamlss |
| `mo()`, `trials()`, `gr()` | not found or refused |
| `ar(p = 1)` | error from `stats::ar()`, not a refusal |
| `rho12 ~ x` | fitted |
| `phylo(tree =)`, `animal(pedigree =)`, `animal(A =)`, `phylo_interaction(tree1, tree2)` | fitted; agree with frmtmb given the matrix |
| `meta_V(V = vi)` (diagonal) and `meta_V(V = M)` (full) | fitted; the full case equals a hand-coded N(Xb, M + s^2 I), logLik -13.82506657 |
| zero-one beta, hurdle NB2, zero-truncated NB2 | fitted (`frm-gaps.log`) |
| binary missing predictor, `impute_model(family = binomial())` | fitted |
| missing response, `miss_control(response = "include")` | fitted |
| `worm_plot()`, `centile_chart()`, `check_drm()` | return ggplots or a diagnostic table |
| `estimator = "mspl"` for binomial with RE | fitted; `logLik()` refused by design |
| `anova()` between two fits | refused for every call |
| `heritability`, `icc`, `repeatability`, `chibar_pvalue`, `lrt_boundary`, `objective_at` | absent from the 0.7.0 namespace (`ns.log`) |
| MCMC or sampling | nothing in the namespace; `simulate()` is the only draw method |

Parameterizations and links, from fits in `surface.log`: beta, NB2,
beta-binomial and Gamma `sigma` are inverse square-root precisions or
shapes. NB2 theta from exp(-2 sigma) is 1.188789, against
`MASS::glm.nb` 1.188793. Student-t nu is 2 + exp(eta). Skew-normal is
the mean parameterization, and its `nu` is brms's `alpha`. Cumulative
logit uses P(y <= k) = logit^-1(theta_k - eta), with the location
intercept dropped and the cutpoints stored as first value then log
increments. `rho12` is 0.999999 tanh(eta).

## What frmtmb does that drmTMB does not, tested in both

frmtmb side in `dev/drmtmb-log/surface2.log`, drmTMB side in
`surface.log`.

| Feature | frmtmb | drmTMB 0.7.0 |
|---|---|---|
| nonlinear formula | fitted, logLik -588.67, convergence 0 | refused |
| three responses with `rescor` | fitted | refused |
| mixed-family multivariate, gaussian + NB2 with shared `|p|` | fitted, convergence 0 | refused |
| `|p|` across mu and shape (NB2) | fitted | refused |
| `|p|` across mu and phi (beta); correlated beta slopes | fitted (singular convergence: no such structure in these data) | refused |
| Student-t sigma RE; `zi` RE; NB2 mu and shape RE together | fitted | refused |
| REML for beta and NB2 | fitted | refused |
| `s()` smooth | fitted | refused |
| custom family as an R log-density | fitted (logistic), convergence 0 | refused |
| brms spellings `trials()`, `gr(cov =)` | accepted (rows k and f use them) | not found or refused |
| sampling | companion package `frmtmb.sample` (not run in this lane) | none in the namespace |

## What drmTMB does that frmtmb does not, tested in both

| Feature | drmTMB 0.7.0 | frmtmb 0.61.0 | In brms 2.23.0? |
|---|---|---|---|
| random-effect SD formula `sd(g) ~ w` | direct grammar | no grammar; the same gaussian model through `nl` (row p) | no grammar; `gr(g, by = f)` for a factor |
| `gr(g, by = f)`, SD per level of a group-level factor | as `sd(g) ~ f`, logLik -584.117229864 | refused: "gr() supports (x \| gr(g, cov = A)) or (1 \| gr(g, prec = Q))" | yes |
| tree or pedigree as direct input | `phylo(tree =)`, `animal(pedigree =)` | refused a tree in `gr(cov =)`: "cov must be a square matrix" | no |
| full known sampling covariance, Cov = V + sigma^2 I | `meta_V(V = M)` | no route; `fcor()` refused | `fcor(M)` exists but means sigma^2 M (Stan code in `brms-parity.log`) |
| zero-one-inflated beta | `zero_one_beta()` | refused | yes |
| hurdle NB2 | `truncated_nbinom2()` + `hu` | refused | yes |
| zero-truncated NB2 | `truncated_nbinom2()` | `trunc()` refused: negbinomial has no CDF | yes, via a CDF |
| binary missing predictor | `impute_model(binomial())` | refused: "mi() responses need a gaussian or student model" | no, same limit |
| worm plot, centile chart | yes | no export of that kind | no |
| `rho12 ~ x` | fitted | not tried | no |
| capability ledger grading each route's evidence | in its documentation | no equivalent | no |

## Feature candidates (recorded, not built)

Sizes are lane-days for one worker including tests and a findings
file. brms parity is the project's first rule, so the brms column
decides order.

1. **Three brms families: `hurdle_negbinomial`,
   `zero_one_inflated_beta`, and a negbinomial CDF for `trunc()` and
   `cens()`.** brms has all three. Files: `R/families.R` (the
   `hurdle_poisson` and `zero_inflated_beta` definitions are the
   templates), the family registry, `R/links-brms.R` if names need
   mapping, and a NB lcdf through the regularized incomplete beta.
   Size: 2 lane-days for all three, 1 of them for the CDF. Validate
   against drmTMB (`truncated_nbinom2`, `zero_one_beta`, now installable
   and agreeing elsewhere to 3e-8), glmmTMB (`truncated_nbinom2`, and
   hurdle through `ziformula`), and `brms::make_standata()` for the
   design.
2. **`gr(g, by = f)`.** brms has it; frmtmb refuses it. It is the brms
   spelling for covariate-dependent random-effect SDs. Files:
   `R/parse.R` (the `gr()` parser), `R/covstruct.R`,
   `R/par-template.R`, `VarCorr()` and `variables()` naming
   (`sd_g__Intercept:f`), and a frame check that `f` is constant within
   `g`, which brms and drmTMB both enforce. Size: 2 to 3 lane-days.
   Validate against drmTMB `sd(g) ~ f` for an intercept (logLik already
   known, -584.117229864 on `sim_grouped(101)`), against lme4 with one
   indicator-coded intercept per level of `f`, and against
   `brms::make_standata()` (`Jby`, `Nby`).
3. **Continuous `sd(g) ~ x`.** Not in brms. frmtmb already fits the
   gaussian random-intercept case through `nl` and agrees with drmTMB
   (row p). Recommendation: document that recipe with a test (0.5
   lane-day) rather than add grammar brms lacks.
4. **Boundary-corrected variance-component LRT.** Not in brms (no LRT
   at all). frmtmb's `anova()` already documents that its chi-square
   p-value is conservative at the boundary (`R/confint.R`, the
   `anova.frmtmb_fit` examples). drmTMB's development version has
   `chibar_pvalue()` and `lrt_boundary()` for q = 1 and 2 independent
   components (Self and Liang 1987), from its source on GitHub `main`;
   I did not run it. Files: `R/confint.R` (an argument to `anova()` and
   the detection of which components were dropped), tests beside
   `test-confint-anova.R`. Size: 1 lane-day for q = 1 and independent
   q = 2; 2 more for general q with information-based weights.
   Validate against `RLRsim::exactRLRT()` for one gaussian component,
   against a parametric bootstrap null built with `frm_simulate()` and
   `frm_bootstrap()`, and against drmTMB's `chibar_pvalue()` once it
   reaches CRAN.
5. **Heritability, ICC, repeatability.** Not in brms (brms users take
   `performance::icc()` or `VarCorr()` arithmetic). A grep of `R/`
   finds no `get_variance` method, so `performance::icc()` probably
   does not reach frmtmb; not run. drmTMB's development version computes
   sigma^2_focal over a total, with a delta-method interval on the
   log-SD scale. Size: 1 lane-day for the `insight` methods, which is
   the brms-parity route, and 1 more for an own accessor with a profile
   interval. Validate against `performance::icc()` on the same lme4
   fit, and against drmTMB's development `icc()` on the animal model
   in row f.
6. **Full known covariance for meta-analysis.** brms has `fcor(M)`
   (Cov = sigma^2 M). drmTMB's `meta_V(V = M)` is Cov = M + sigma^2 I,
   which in brms is `fcor(M)` with sigma fixed at 1 plus an
   observation-level random intercept. Implementing `fcor()` gives
   brms parity and, with a fixed sigma, drmTMB's model. Files:
   `R/parse.R` (`fcor` is listed there and refused),
   the residual-correlation code beside `unstr()`. Size: 1 to 2
   lane-days. Validate against the hand-coded likelihood in
   `dev/drmtmb-surface2.R` (-13.82506657), drmTMB, and
   `metafor::rma.mv(V = M)`.
7. **Skew-normal start (the defect above).** 0.5 lane-day for the
   two-start route, which is the one that measured 0 of 40 on both
   streams; the residual-skewness route is cheaper and insufficient
   (section 1.1). Not a feature, listed here so the size is recorded.
8. **Lower priority, not brms features:** a tree accepted by
   `gr(cov =)` (convert with `ape::vcv(corr = TRUE)`, 0.5 lane-day);
   worm plots from the existing quantile residuals (0.5 lane-day);
   non-gaussian `mi()` predictors (brms has the same limit as frmtmb).

## The gated test file

`tests/testthat/test-drmtmb-agreement.R`: 13 tests, 128 expectations.

**Three gates, and why the third one exists.** `skip_on_cran()`,
`skip_if_not_installed("drmTMB", "0.7.0")`, and
`FRMTMB_DRMTMB_FIT_TESTS = "true"`. The variable follows
`helper-brms.R`'s `skip_unless_brms_fit()` and its
`FRMTMB_BRMS_FIT_TESTS`, because this file is a fit tier, not a
structural one: every test fits models, about 50 in all, in two
packages. `drmTMB (>= 0.7.0)` stays in DESCRIPTION Suggests, because
the file calls `drmTMB::` and the declaration is what keeps
`R CMD check` quiet about it.

**The cost of not gating**, stated plainly: `R-CMD-check.yaml` names no
`dependencies:`, so `setup-r-dependencies@v2` installs all Suggests,
and `check-r-package@v2` sets `NOT_CRAN=true`. Without the variable,
every check job would install drmTMB and run about 50 fits: four jobs
over three platforms (windows, macos, ubuntu release, ubuntu devel).
The r-devel job would most likely build drmTMB from source, because
the public RSPM binaries are built for the release version, and drmTMB
compiles a TMB template. Measured locally, the fits alone are 47 s.

**Where it runs.** Locally, through the release gated tier:
`dev/release/run-gated.ps1` now exports `FRMTMB_DRMTMB_FIT_TESTS` and
lists `test-drmtmb-agreement.R` beside the brms tiers, one R process
per file, which is where a local gated run belongs.

**In CI: nothing added, deliberately.** The only tier workflow that
could host it, `brms-likelihood.yaml` (job name `brms-tiers`), is
scoped to the two tiers that compile Stan, says so in its header, and
its name is referenced by branch protection; a drmTMB job there would
need its own dependency install anyway, so it would not be cheaper
than a separate workflow. A separate workflow would cost one
ubuntu-latest job per push or pull request touching `R/`, `tests/` or
`DESCRIPTION`: checkout, setup-r and a Suggests install (drmTMB is a
binary from public RSPM on release), then about 1 minute of fits, so
roughly 3 to 5 minutes per run and a new required check to maintain.
It is not added here. If the user wants the tier in CI, that estimate
is the input to the decision.

Runs, one file per process with `dev/drmtmb-run-test.R`. All four
states were verified:

| State | Result | Log |
|---|---|---|
| gate set, drmTMB installed, `NOT_CRAN=true` | 13 tests, 128 expectations, 0 failed, 0 error, 0 skipped, 83 s (57.8 s on a quieter run; this box times as load) | `test-with.log` |
| gate unset, drmTMB installed, `NOT_CRAN=true` | 13 skipped | `test-nogate.log` |
| gate set, drmTMB absent, `NOT_CRAN=true` | 13 skipped | `test-without.log` |
| gate set, drmTMB installed, `NOT_CRAN` unset | 13 skipped | `test-cran.log` |

Mutation runs, each with ONE map changed in a copy of the file:

| Mutation | Failures | Log |
|---|---|---|
| correlation cap removed (cc = 1) | 4: the off-optimum gap in every one of the four comparisons that carry a correlation map, at 6.1e-10 to 2.2e-8 of \|logLik\| against a bound of 1e-11 | `test-mut1.log` |
| phylo map without log(h) / 2 | 3 | `test-mut2.log` |
| NB2 sigma-RE map without log 2 | 6 | `test-mut3.log` |
| Student-t nu map replaced by the identity | 4 | `test-mut4.log` |

Tolerances are ratios, as `dev/lane-rules.md` requires.
`expect_same_model()` makes five assertions, and all five margins are
measured by `dev/drmtmb-test-margins.R` over its 24
`expect_same_model()` calls:

| Assertion | Bound | Worst measured | Margin |
|---|---|---|---|
| gap between the two optima, over \|logLik\| | 1e-9 | 5.22e-11 | 19.1x |
| gap at drmTMB's optimum, over \|logLik\| | 1e-9 | 3.79e-11 | 26.4x |
| gap one SE off the optimum, over \|logLik\| | 1e-11 | 1.39e-13 | 72.1x |
| estimate gap, over the standard error | 1e-2 | 1.49e-04 | 67.2x |
| \|log standard-error ratio\| | 1e-3 | 1.83e-05 | 54.8x |

The off-optimum point moved from a seeded random sign pattern to a
fixed alternating one, so that the function leaves the RNG stream
alone. Its margin went from 23x to 72x, and the cap mutant now fails in
all four correlation comparisons rather than three: with the seeded
point and the earlier 1e-9 bound, the mu-sigma map did not fail.

What the file pins as a difference rather than agreement: the
skew-normal default-start stall (flip when fixed), and, for REML with a
sigma random effect, frmtmb's nesting plus drmTMB's integrated
`beta_sigma` (it fails if drmTMB changes its criterion, which is the
point).

drmTMB fits are read through `$obj`, `$opt` and `$sdr`, which are not
drmTMB API. A drmTMB release that changes them breaks the file loudly.

## R CMD check

`dev/drmtmb-check.sh`, run once on the final tree: `R CMD build` with
vignettes, then `R CMD check --as-cran` on the tarball, with
`drmtmb-lib` first on `R_LIBS`, `pinlib` before the user library,
pandoc and TinyTeX on `PATH`, and
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`. Logs:
`dev/drmtmb-log/check-build.log`, `check-console.log`,
`check-00check.log`.

**Status: 2 NOTEs.** The release check on main
(`dev/release/check.log`) records 1 NOTE for core, the V8 one, so the
second NOTE is this run's and is not called pre-existing here:

- `checking HTML version of manual ... NOTE`, "Skipping checking math
  rendering: package 'V8' unavailable". This is the NOTE the release
  check also carries.
- `checking examples ... [76s] NOTE`, listing `residuals.frmtmb_fit`
  at 5.95 s elapsed and `pp_check` at 5.23 s. `dev/lane-rules.md`
  documents an examples-timing NOTE on this box as a measurement of
  LOAD, with a fixed-arithmetic control that swings a factor of 4.7 on
  identical work, and it records the same example at 1.05 s quiet
  against 5.28 s busy. Two other lanes were running during this check.
  This lane changed no `R/` file and no example, so it cannot have
  moved these times; what it cannot do is prove the NOTE would be
  absent on a quiet machine without rerunning the check there, which
  the once-per-lane budget does not allow.

The lines the Suggests addition had to clear: `checking package
dependencies ... OK`, `checking for unstated dependencies in 'tests'
... OK`, `checking DESCRIPTION meta-information ... OK`, `checking
tests ... [553s] OK`. The vignettes rebuilt (488 s). Under
`--as-cran`, `NOT_CRAN` is unset, so the new file's 13 tests skip
inside the check; they are run separately, above.

## What I did not do, and why

- **drmTMB's development version was not installed or run.** The brief
  pins the CRAN binary. The claims about `heritability()`,
  `lrt_boundary()` and `objective_at()` come from reading
  `R/heritability.R`, `R/lrt-boundary.R` and `R/objective-at.R` at
  commit 54df129. `objective_at()` was not needed: evaluating each
  package's TMB objective at the other's estimates is what the engine
  does.
- **No `spatial()` comparison.** frmtmb's spatial routes (SPDE, CAR, GP)
  and drmTMB's `spatial(coords =)` would need their covariance
  functions matched first; not attempted.
- **No brms fits.** brms appears only through its namespace and
  `make_stancode()` (`dev/drmtmb-brms-parity.R`), for the parity
  column.
- **frmtmb's full test suite was not run.** No `R/` code changed. The
  new file and `helper-reference.R` (a comment edit) are the only test
  changes, and the new file ran as above.
- **No CI workflow for the tier.** The reasoning and the cost estimate
  are under "The gated test file"; the local route is
  `dev/release/run-gated.ps1`.
- **Coverage, not only agreement, was out of scope.** Agreement says
  the two packages compute the same numbers, not that the intervals
  cover.
- **The skew-normal rate is for one designed case.** I did not estimate
  how often field data put raw and residual skewness on opposite sides.
- **Neither skew-normal fix is implemented.** Section 1.1 measures
  both; the route is the user's decision.

## Files

Changed: `README.md`, `CONTRIBUTING.md`, `DESCRIPTION` (Suggests),
`dev/submission-draft.md`, `dev/release/run-gated.ps1` (the gate and
the file), `tests/testthat/helper-reference.R` (srr comment).
New: `tests/testthat/test-drmtmb-agreement.R`, this file, 21 R
scripts `dev/drmtmb-*.R`, `dev/drmtmb-check.sh`, and `dev/drmtmb-log/`.
