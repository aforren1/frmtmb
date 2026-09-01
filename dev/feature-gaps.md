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
theta dpars. Validated vs hand-rolled ML to 1.3e-10.

The covariance taxonomy landed next as
`mixture_mvn(K, D, model = )`, in mclust's own model-name vocabulary:
EII, VII (spherical), EEI, VEI, EVI, VVI (diagonal), EEE (one shared
full covariance) and VVV (the default, free per class). The extras
shrink with the model - a single log-SD for EII, one `us` block for
EEE, `sigmavol*`/`sigmashape*` for the volume-shape pair - so
confint() and frm_sample() label exactly the parameters a model has.
Validated against mclust::Mclust() with intercept-only means (where
the two models are the same likelihood): starting our ML at mclust's
EM solution, logLik agrees to 1e-12, class means and covariances and
posterior z to 1e-15, on faithful (K=2, D=2) and a simulated K=3,
D=3 set; from our own cold start we reach mclust's optimum to 1e-8.
The eight likelihoods are bit-identical at a shared spherical-equal
parameter point, and the fitted maxima respect every nesting chain.

Remaining sub-gaps, logged in ?mixture_mvn: the six models with a
class-varying eigenvector basis (EEV, VEV, EVE, VEE, VVE, EVV) need
constrained-orientation machinery and are absent; covariances take no
linear predictor (no covariance regression); no cens/trunc, no mvbf
components (mixture over mvbf with rescor is the general form;
unscheduled), and simulate() needs a simulator interface that can see
family-level extras.

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
- ~~`car(M, gr, type =)` (brms spelling) and an SPDE-Matern
  covstruct~~ DONE in v0.28, together with `gr(prec=)` beyond
  intercept-only. All three pieces shipped.

  (1) `car(M, gr = g, type =, con_sd =)` is a predictor special (like
  `gp()`), not a bar term, and carries brms's whole argument list and
  all four types. The block looks like `(1 | g)` - one intercept per
  location, `dim = 1`, a synthetic bar - so `ranef()`, `predict()`
  with `newdata`, `VarCorr()`, `simulate()` and `draw_b()` ride the
  existing machinery; only `confint_varcorr()` needed new rows
  (`sd(car)` plus brms's `car` / `rhocar` on a new `"prop"`
  back-transform, logit rather than log or Fisher-z).

  No normalize trick is needed anywhere, because every
  log-determinant is analytic in the parameters and the fixed part is
  computed once at frame time. `escar`: log|Q| = n log tau +
  sum log d_i + sum log(1 - rho e_i), with e the eigenvalues of
  D^-1/2 W D^-1/2 (fixed). Intrinsic: the graph Laplacian is rank
  n - c for c connected components, and the constrained precision
  tau (L + sum_j kappa_j s_j s_j') is full rank with log|Q| = n log
  tau + log|K| - so the density costs one sparse matrix-vector
  product and NO on-tape factorization at all. `bym2` is the one
  exception: the Riebler mixture has no sparse precision, so it uses
  the dense marginal covariance sd^2[(1 - rho) I + (rho/scale) K^-1]
  with brms's `.car_scale` (reproduced exactly, tested against the
  formula).

  CONSTRAINT. The intrinsic types take brms's soft sum-to-zero
  constraint with its precision riding on tau (as in brms's
  non-centered zcar), applied per connected component. That keeps the
  density proper - which is what makes ranef/predict/simulate defined
  on the block - and keeps log|Q| exact. `con_sd` (default brms's
  1e-3) is the constraint sd relative to the field sd; the fit
  converges quadratically onto the hard-constrained (esicar)
  likelihood as it shrinks. Measured on a 4 x 4 lattice against a
  hand-rolled hard sum-to-zero ML: 1e-3 off by 4.7e-4 in logLik
  (3.6e-5 relative in sdcar), 1e-4 by 4.7e-6 (3.6e-7), 1e-5 by 4.5e-8
  (3.7e-9), 1e-6 by 9.2e-10, 1e-7 lost to roundoff (1.1e-4). The
  default stays at brms's value: the bias is four orders below the
  parameter's own SE, and tighter settings cost optimizer robustness
  (over 25 lattice refits nlminb reported false convergence 0 times
  at 1e-3, once at 1e-4, 6 times at 1e-5). `esicar` selects the same
  density as `icar` - the brms difference between them is a Stan
  parameterization detail that ML does not see.

  VALIDATION (tests/testthat/test-car-spde.R), all against a
  hand-rolled marginal-gaussian direct ML with the same block
  covariance, which the Laplace approximation reproduces exactly:
  icar 4.8e-9 on the logLik, escar 1.1e-11, bym2 2.9e-12, spde
  2.6e-11, plus a disconnected two-component graph (the rank
  correction and the per-component constraint) and simulate-recover
  over 15 lattice replicates for sdcar and rhocar. No opportunistic
  cross-check was possible: neither CARBayes, spaMM, spdep, fmesher
  nor INLA is installed here.

  (2) `gr(g, prec = Q)` now takes correlated slopes. Precision-side
  Kronecker: inv(A (x) Sigma) = A^-1 (x) Sigma^-1, so the block
  precision is Q (x) Sigma^-1, assembled on the tape as an AD-weighted
  sum of the d(d+1)/2 fixed sparse matrices Q (x) E_ab and staying as
  sparse as Q. Exact against the dense `gr(cov = solve(Q))`
  equivalent: 0 ulp on the logLik, 3e-15 on theta, 5e-16 on b.
  `RTMB::solve`, not base's - the S4 advector method is not imported.

  (3) `spde(fem, gr = node)` takes the mesh's finite-element triple
  (fmesher `c0`/`g1`/`g2` or INLA `M0`/`M1`/`M2`) as fixed data and
  tapes Q = tau^2 (kappa^4 M0 + 2 kappa^2 M1 + M2) as an AD-weighted
  sparse sum through `RTMB::dgmrf`; theta = (log tau, log kappa), with
  sd and range reported through the planar alpha = 2 identities.
  Validated on a 1-D chain whose FEM matrices are known tridiagonals
  (fmesher is not installed).

  RESIDUE. Mesh and adjacency construction stay out of scope
  (fmesher/spdep are preprocessing, same posture as HRF convolution).
  `spde()` maps observations to mesh nodes through a grouping factor,
  not through a general barycentric projector matrix `A`, which is
  what a real fmesher workflow wants next. Any sum-to-zero constraint
  puts a dense rank-c update in the block's Laplace Hessian, and
  `bym2` is dense outright, so the practical field size is in the low
  thousands rather than the 1e5 nodes sdmTMB/VAST reach. `car()`
  refuses unseen locations at prediction time (as brms does).

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
- ~~`getME()`: small accessor vocabulary (X, Z, theta, beta, b) on the
  lme4 generic~~ DONE in v0.28: `getME.frmtmb_fit` registered on
  `lme4::getME` (Suggests, so the registration is delayed) with the
  vocabulary X, Z, Zt, beta/fixef, b, theta, lower, sigma, flist,
  n_rtrms, n_rfacs; anything else errors listing it. X and Z agree
  with `lmer`'s on sleepstudy; `b` is coefficient space (rr blocks
  expanded), `theta` is the internal unconstrained scale and `lower`
  is therefore all `-Inf`, both documented as such. Design extractors
  need `resp=` on a multivariate fit; the scalar names answer without
  one.

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

## ODE models via RTMBode (user request 2026-09-01)

RTMBode (kaskr, r-universe only, NOT CRAN) wraps deSolve's ode() with
adjoint AD and taped dynamics: ode(y, times, func, parms, method,
...) differentiates through BOTH the dynamics parameters and the
initial states (its Lotka-Volterra example estimates pars, yini, and
an observation SD through MakeADFun + nlminb + sdreport). brms has
NO ODE grammar (stanvars + hand-written Stan is its answer), so this
is union-coverage territory (pharmacokinetics, epidemic models,
nlmixr-adjacent).

Path: (1) feasibility probe - nl = TRUE bodies are arbitrary R
evaluated on the tape, so RTMBode::ode() inside an nl body may
already nearly work (population PK: nlpars with REs feeding the
dynamics parameters, mu read off the solved trajectory); establish
what breaks (times/data plumbing, matrix indexing on the tape,
sdreport). (2) If viable, sugar later: an ode-aware term or an
odefun = argument on bf(); keep deSolve method/atol/rtol
passthrough.

Packaging (the non-CRAN dependency): the established CRAN mechanism
is Suggests + Additional_repositories, used ON CRAN today by brms,
bayesplot and bridgesampling for cmdstanr (Additional_repositories:
https://stan-dev.r-universe.dev/) and by the INLA ecosystem
(inlabru, sdmTMB). So: Suggests: RTMBode (and RTMBp if the parallel
backend ever ships), Additional_repositories:
https://kaskr.r-universe.dev, every ODE feature behind
requireNamespace("RTMBode") with an informative install message,
tests skip_if_not_installed, compat registry rows conditional.
CRAN's check farm does NOT install additional-repo packages, so all
examples/vignettes must skip cleanly without it - the brms/cmdstanr
posture exactly. The RTMB/RTMBp split (variant under a different
package NAME) is the other pattern, for when the variant replaces
rather than extends; not needed here.
