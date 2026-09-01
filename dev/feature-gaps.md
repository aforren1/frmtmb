# Feature gaps vs lme4 / glmmTMB / brms vignettes

Source: sweep of the installed vignettes of all three packages
(2026-08-31, post v0.13). Items already rejected or delivered are at
the bottom.

## Tier 1: cheap, real value — DONE in v0.14

All items below shipped in v0.14 (plus `frm_simulate()` from Tier 2
and the `anova()` duplicate-label fix). Notes: `profile = TRUE` is the
glmmTMB-style approximation (betas ride inside the Laplace step), and
`condVar` follows the TMB/glmmTMB convention, not lme4's
purely-conditional `postVar`.

- `ranef(fit, condVar = TRUE)`: conditional SDs of the modes from
  `sdreport`'s `diag.cov.random` (already computed for smooth edf).
  Both lme4 (`postVar` attribute, caterpillar plots) and glmmTMB
  (`condsd` column) provide this; we provide nothing.
- `as.data.frame()` methods for `ranef` (grp/term/condval/condsd) and
  `VarCorr` (grp/var1/var2/vcov/sdcor): the tidy shapes broom.mixed
  and plotting workflows expect.
- `se()` addition term: known per-observation sampling SD
  (meta-analysis). Pairs with `gr(g, cov = A)` for phylogenetic
  meta-analysis; brms treats this as core. Small: fixed additive
  component in the sigma of gaussian-type families.
- Binomial proportion response + `weights = size` (the glmer idiom).
  We currently require integer counts + `trials()`.
- Family quick wins: `bernoulli` (alias), `geometric` (negbinomial
  with shape fixed), `exponential` (Gamma shape 1), `weibull`,
  `shifted_lognormal` (ndt dpar), `hurdle_gamma`, `hurdle_lognormal`,
  `zero_inflated_binomial`, `zero_inflated_beta`, `asym_laplace`
  (quantile regression: fixed-quantile dpar).
- `frm_allfit()`: refit under every registered optimizer and tabulate
  logLik/fixef/theta agreement (lme4 `allFit`). Trivial given
  pluggable optimizers + `update()`.
- `MakeADFun(profile = "beta")`: profile fixed effects out of the
  inner problem (glmmTMB `glmmTMBControl(profile = TRUE)`, glmer
  `nAGQ = 0` relative). Speedup for many-coefficient models; expose
  as a `frmtmb_control()` flag after benchmarking.

## Tier 2: moderate effort - DONE in v0.15 except where noted

- ~~De novo simulation~~ DONE in v0.14 (`frm_simulate()`).
- ~~mo() monotonic effects~~ DONE (standalone additive terms;
  interactions and group-level mo() remain out).
- ~~sratio / cratio / acat~~ DONE. ~~cs()~~ DONE in v0.17
  (sratio/cratio/acat; refused for cumulative).
- ~~influence() + cooks.distance~~ DONE.
- ~~hetar1 / homcs / homtoep / exp / gau / mat~~ DONE. ~~equalto~~
  DONE in v0.17 (zero-parameter registry entry).
- ~~rr(..., d = k)~~ DONE in v0.16 (two-space b: factors in the
  parameter vector, expand_b() maps to the coefficient space the Z
  matrices span; exact vs glmmTMB; se.fit deferred).
- ~~conditional_effects(method = "predict") + data-frame
  conditions~~ DONE.
- ~~vint() / vreal()~~ DONE (Wiener diffusion showcase in
  test-v15.R).

## Missing data (assessed 2026-08-31)

- `frm_multiple()` (Rubin pooling over imputations, mice interop):
  DONE in v0.15; pooled model comparison (D1/D2/D3) added on top of
  it, see the Tier 3 entry below.
- ~~In-model `mi()`~~ DONE in v0.16 for continuous predictors
  (gaussian/student imputation models; validated against the
  closed-form linear-gaussian marginal). Discrete predictors stay
  impossible (no discrete latents in the Laplace class - same wall
  Stan has). `me()` measurement error is a small follow-on (mi with
  a known-error observation model). mi() in dpar formulas, smooths
  of mi() variables, and mi() interactions remain out.

## Mixture models (assessed 2026-08-31)

- ~~Observation-level `mixture(fam1, fam2, ...)`~~ DONE in v0.17
  (logsumexp lpdf via RTMB::logspace_add; suffixed component dpars;
  mixing weights are multinomial-logit dpars with full linear
  predictors, so mixture-of-experts is free; quantile-spread mean
  inits; exact vs direct ML). Multimodality is inherent to mixture
  ML - use frm_allfit / bounds.
- ~~Group-level mixtures~~ DONE: v0.18 latent-class regression,
  v0.19 WITH continuous random effects (growth mixtures). The
  earlier "needs per-group per-class Laplace / crossed impossible"
  analysis was wrong: swapping sum and integral, the marginal is
  exactly the integral of a per-group class-mixture conditional on
  b, so ONE Laplace over b handles it and crossed designs are
  structurally fine. What remains approximate is the Laplace of the
  (non-Gaussian) mixture integrand: ~0.1 logLik bias in validation,
  parameters to 0.01; quadrature exact for univariate per-group
  integrands, approximate when class-specific intercepts couple.

## Multivariate mixture components (noted 2026-08-31, vs clustTMB)

clustTMB (Havron) does model-based clustering of MULTIVARIATE
gaussian observations (mclust-style per-cluster mean vectors and
covariance matrices, VVV etc.), with covariates and spatial GMRFs in
the gating and expert parts. DONE in v0.20 as `mixture_mvn(K, D)`:
per-class D-dim means as full linear predictors (covariates, REs) +
unstructured per-class covariance; gating on covariates via the
theta dpars. Validated vs hand-rolled ML to 1.3e-10. Remaining
sub-gaps, logged in ?mixture_mvn: covariances are `us`-only and
covariate-free (no EII..VEV taxonomy, no covariance regression); no
cens/trunc, no mvbf components (mixture over mvbf with rescor is the
general form; unscheduled), and simulate() needs a simulator
interface that can see family-level extras.

## Tier 3: positioning decisions

- ~~`frm_multiple`~~ DONE in v0.15; v0.20 adds varcorr pooling on
  transformed scales and Rubin-pooled hypothesis(). Fixed effects
  match mice::pool to 2.1e-6. Pooled model comparison is now DONE
  too: `anova.frmtmb_multiple(method = c("D3", "D1", "D2"))`, with
  `use = c("likelihood", "wald")` for D2 and a `constraint =` form
  for the Wald rules. D3 is the default and is the reason our
  architecture helps: the Meng-Rubin pooled-estimate leg needs each
  imputation's likelihood at one common parameter vector, which is
  `obj$fn(pbar)` on the stored objective - one Laplace solve per
  imputation, no refit. Parameters are pooled on the optimizer's own
  (internal) scale, so variance components pool as log/Cholesky,
  consistent with `$pooled`.
  Validation (tests/testthat/test-pooled-anova.R): a poisson GLM has
  no dispersion parameter, so our ML coefficients, inverse-information
  vcov and `df.residual` all equal `glm()`'s, and D1 / D2(wald) match
  `mice::D1` / `mice::D2` end to end to 1e-5 relative (8.5e-8 on the
  statistic; the residue is optimizer noise, BFGS at reltol 1e-14),
  D2(likelihood) to 1e-8 (3.8e-12 on the statistic). Two convention
  differences worth remembering. (1) `df.residual()` counts a gaussian
  fit's dispersion parameter and `lm()` does not, so our default
  `dfcom` is one smaller than mice's; `dfcom =` overrides it. (2)
  `mice::D3()` does not plug the pooled estimates in - `fix.coef()`
  makes them an offset and refits `. ~ 1`, which re-estimates the
  dispersion AND frees an intercept. On a gaussian pair that put mice
  10% away on the statistic and 77% on the p-value, while our value
  is within 5.8e-4 of `mitml::testModels(method = "D3")` (the plug-in
  form; the rest is mitml pooling sigma^2 on the variance scale where
  we pool log sigma) and within 1.1e-13 of the Meng-Rubin formula
  written out by hand.
- ~~Sparse X option~~ DONE in v0.20 (`frmtmb_control(sparse_x =
  TRUE)`, glmmTMB `sparseX` analog). Identical estimates (gradient
  3e-14), 13.8% tape memory saved; dense fallback for NA-factor
  newdata rows.
- ~~`autoscale` internal predictor scaling~~ DONE in v0.20
  (lme4 >= 1.1.37 analog): scaled pre-fit + exact back-transform +
  warm start, optimizer/convergence/Hessian in scaled units. 12
  scaling warnings to 0, logLik identical to manual standardization.
- Sandwich/robust SEs (`vcovHC`, `bread`/`estfun`): still skipped;
  glmmTMB does cluster-level scores. Revisit only on demand.
- `car(M, gr, type = "icar"/"bym2"/"escar")` (brms spelling) and an
  SPDE-Matern covstruct (noted 2026-08-31, fMRI/spatial discussion).
  The backend is ready: sparse GMRFs are TMB's home turf (dgmrf +
  sparse Laplace Hessian; sdmTMB/VAST scale to 1e5+ latent nodes),
  and `gr(g, prec = Q)` already tapes an advector-scaled sparse Q.
  Missing pieces are grammar and hyperparameters, not capability:
  (1) Q(theta) assembled on the tape as an AD-weighted linear
  combination of fixed sparse matrices - ICAR from adjacency
  (needs sum-to-zero constraint + the TMB normalize trick for the
  parameter-dependent constant), BYM2 mixing, SPDE
  tau^2(kappa^4 C + 2 kappa^2 G1 + G2); (2) `gr(prec=)` beyond
  intercept-only. Mesh/adjacency construction stays out of scope
  (fmesher/spdep are preprocessing, same posture as HRF
  convolution). vs references: brms has car() under full MCMC
  pricing; glmmTMB has no CAR; sdmTMB owns SPDE but is
  fixed-likelihood. Deferred since v0.1; schedule on spatial
  demand.

## Method-surface residue (v0.21 audit)

- ~~Deviance residuals (`residuals(type = "deviance")`,
  lme4/glmmTMB)~~ DONE in v0.24 as a per-family `post$dev_fn` unit
  deviance. Covers gaussian, poisson, binomial, bernoulli, Gamma,
  exponential, inverse.gaussian, negbinomial/nbinom2, nbinom1,
  geometric, beta, tweedie; every other family is refused by name
  (ordinal/mixture/multinomial/zi/hurdle have no standard unit
  deviance), as are `trunc()` and `cens()` responses. Exact vs
  `stats::glm` (0 ulp at a shared optimum) for the five glm families,
  exact vs glmmTMB's `dev.resids` for both nbinom parameterizations
  (nbinom1 follows their convention: the size stays at the fitted
  `mu / phi`), and vs the saturated log-likelihood for beta and
  tweedie, which glmmTMB returns NA for. `deviance()` is unchanged:
  still `-2 logLik` (lme4). Dunn-Smyth stays covered via
  `dharma_residuals()` and OSA.
- ~~`se.fit` for the expected response of non-identity-mean
  families~~ DONE in v0.24: the delta method runs jointly over every
  dpar linear predictor against the joint coefficient covariance, so
  zi/hurdle, lognormal, trials-binomial and `trunc()`ed responses all
  answer. Gradients `dm/deta_k` are central differences (one step per
  predictor, relative), validated against the analytic zi-poisson and
  lognormal forms to 1e-8 and against glmmTMB's own response-scale
  delta method to ~1e-5 relative. The identity-mean path is
  bit-identical to before.
- `getME()`: small accessor vocabulary (X, Z, theta, beta, b) on the
  lme4 generic; low value until a concrete downstream consumer
  appears (design access exists via `fit$frame`).

## Already ahead (no action)

- Custom families are plain R functions with AD verification vs
  brms's stanvars/Stan-code injection and glmmTMB's C++ enum-edit +
  recompile chapter (their "hacking" vignette is our `custom_family()`
  one-liner).
- Priors: brms-style `set_prior` with natural-sd-scale Jacobian;
  glmmTMB's priors vignette lists LKJ/elementwise as "planned".
- Post-hoc MCMC: `frm_sample`/`check_laplace` vs glmmTMB's DIY
  Metropolis vignette.
- OSA residuals, DHARMa, pp_check, quadrature matching GLMMadaptive,
  hypothesis wald/profile/boot, warm-start bootstrap.
- `||`, nesting/crossing, offsets, `numFactor`, dispformula-equivalent
  (dpar formulas), zi/hurdle with full formulas: all present.

## Explicitly out (locked or moot)

- ~~`me()` measurement error~~ DONE in v0.17 as `x | mi(sdx)` (the
  brms spelling since me() was folded into mi). Discrete-predictor
  mi(): impossible.
- ~~`gp()`~~ DONE in v0.18: exact (dense SE kernel + nugget) and
  Hilbert-space (`k =`) forms. Extended post-v0.19: up to 3
  dimensions (per-dimension lengthscales by default, `iso = TRUE`
  for a shared one) and kriging prediction for the exact form at
  unseen positions.
- ~~mo()/mi() interactions~~ DONE in v0.18 (two-way `:`/`*` with
  numeric terms; shared simplex per mo variable).
- OpenMP objective parallelism (RTMB limitation; benchmarks fine
  without it). brms threading, glmmTMB parallel vignettes are moot.
- Compilation management (precompile, cmdstanr backends): moot,
  nothing ever compiles.
- Group-level mixtures for crossed designs (see the mixture section:
  out of the Laplace class).
