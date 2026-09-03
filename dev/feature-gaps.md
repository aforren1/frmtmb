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
unscheduled).

The simulator sub-gap is CLOSED in v0.36. The family-level extras a
`mixture_mvn()` draw needs are reachable through the structured
simulator contract (`fam$sim_ctx(ctx)`, R/families.R): one
implementation per family, called identically by `simulate()`,
`posterior_predict()` and `frm_simulate()`. The same contract carries
`hmm()`'s per-sequence chain walk and `mixture(groups =)`'s per-group
class draw, and the frame-level residual-correlation draw that
`ar()`/`ma()`/`cosy()` need.

The MEASUREMENT-MODEL case of that general form is DONE as `lca(K)`
in v0.35: a mixture over a vector of conditionally independent
CATEGORICAL responses, which is poLCA's latent class analysis. It
takes the same contained shape mixture_mvn() took (one matrix
response, gating on covariates, no random effects) and for the same
reason: the item profiles are family extras rather than K component
family objects, so none of it generalizes to arbitrary component
families with rescor. The general mixture-over-mvbf stays
unscheduled. Validated against poLCA on carcinoma (K = 3, logLik to
4.2e-8, profiles to 2.8e-8) and on the election latent class
regression (logLik to 1.1e-7, gating coefficients to 9.9e-7).

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
- ~~brms R-side autocorrelation (`ar()`, `ma()`, `arma()`, `cosy()`,
  `unstr()`)~~ DONE in v0.34. The residual of a group becomes one
  multivariate draw, `y_g ~ N(mu_g, D R D)` with `D` the diagonal of
  that group's `sigma` values, so a `sigma ~ x` distributional model
  enters through the diagonal and a scalar `sigma` is the homogeneous
  special case. `R` is UNIT-DIAGONAL, so `sigma` stays the marginal
  residual SD everywhere in the package; brms instead parameterizes
  `ar`/`ma`/`arma` by the innovation SD (`cholesky_cor_ar1()` divides
  by `1 - ar^2`) while `cosy`/`unstr` use the marginal one, an
  inconsistency we deliberately do not reproduce.

  Where the pieces live: a `thetaac` component of the parameter
  template (outer, so REML integrates only the mu fixed effects, as
  with `theta`); one entry per response in `frame$autocor`; the
  density replaces `fam$lpdf` in `build_objective()`. The block is
  evaluated one PATTERN at a time (a distinct set of present time
  levels), so ragged and unsorted groups are handled by construction:
  a balanced design costs one on-tape Cholesky per gradient, ragged
  data one per distinct pattern, and nothing ever builds an `n x n`
  matrix.

  ALL FIVE STRUCTURES SHIPPED, including `ma()` and `arma()`. The
  ARMA autocorrelation function is exact rather than a truncated
  MA(inf) expansion - psi weights up to lag q, then the
  `max(p, q) + 1` linear system of Brockwell & Davis 3.3.1 solved with
  `RTMB::solve` - and agrees with `stats::ARMAacf` to 1e-16.
  Stationarity and invertibility come from the Monahan/Jones partial
  autocorrelation transform (`nlme::corARMA`'s), so no parameter value
  leaves the parameter space. Higher orders work; brms limits
  `cov = TRUE` to order one.

  Validated against `nlme::gls`/`lme` under BOTH ML and REML. Worst
  log-likelihood gap over the suite: AR(1) 3e-12, AR(2) 4e-10, MA(1)
  4e-13, ARMA(1,1) 8e-10, corCompSymm 2e-11, corSymm 3e-9, ragged
  AR(1) 1e-10, AR(1) + varIdent 3e-10, `lme(random = ~ 1 | subj,
  correlation = corAR1())` 9e-10. Correlation parameters agree to
  1e-7 or better except `unstr` (2.6e-5) and the ARMA/AR(2) ML fits
  (3e-6), where the surface is flat enough that both optimizers stop
  inside a region of that width - the log-likelihood agreement is the
  binding evidence there.

  Refused, each because the likelihood no longer factorizes over
  rows: `weights()`, `cens()`, `trunc()`, `se()`, `mi()`,
  `rescor = TRUE`, mixtures, `quadrature = TRUE`, non-linear (`nl`)
  formulas, `frm_simulate()` and `residuals(type = "osa")`. brms
  refuses the same core set. Also refused: brms's `cov = FALSE`
  default, which is the residual-REGRESSION formulation and a
  different likelihood; and every family but gaussian and student,
  where brms silently switches to a latent AR process on the linear
  predictor - which is a random effect and is spelled
  `ar1(factor(week) + 0 | subj)` here.

  Divergence from brms worth remembering: the lag is the distance
  between the rows' positions in the GLOBAL set of time levels, so a
  group missing a time point gets the wider lag. brms indexes
  `ar`/`ma`/`arma`/`cosy` by position WITHIN the group (its Stan code
  slices `chol_cor[1:nobs[i], 1:nobs[i]]` and carries no time index;
  only `unstr` has `Jtime_tg`), so brms treats a missing row as no
  gap. Ours is nlme's reading and is what the ragged gls agreement
  test pins down.

  Remaining, in rough order of value: (1) `se(x, sigma = TRUE)`, which
  brms supports here by adding `diag(se^2)` to the group covariance -
  a small addition to `autocor_loglik()` plus a relaxation of the
  aterm guard; (2) `set_prior()` targeting the correlation parameters
  (bounds on `thetaac_*` already work); (3) `frm_simulate()` support,
  which needs the de novo path to draw per group; (4) `rescor` +
  a time structure, whose joint covariance is the Kronecker product of
  the two; (5) a `cov = FALSE` residual-regression form, if anyone
  actually wants brms's default rather than the marginal likelihood.

- ~~Sandwich/robust SEs (`vcovHC`, `bread`/`estfun`)~~ DONE as
  `vcov_cluster()` / `cluster_scores()`; see `dev/sandwich.md` for the
  design record. The per-OBSERVATION half of this note stands: a
  marginalized objective has no per-observation contribution, so there
  is still no `estfun()` method. Per-CLUSTER scores are exact, and are
  what shipped.
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
- ~~`data2 = list()`: the five sites that captured structural matrices
  from the calling environment~~ DONE in v0.29. All five go through
  one helper, `lookup_structural(expr, data2, data, env, what)`:
  `gr(g, prec = Q)`, `gr(g, cov = A)`, `equalto(..., V)`,
  `car(M, ...)` and `spde(fem, ...)`. Resolution order is data2 (bare
  name), then the expression evaluated with data2 in front of the data
  mask, then the historical data-then-formula-env path. The middle
  step is a deliberate divergence from brms, whose rule is name-only:
  it makes `gr(g, cov = solve(Q))` work with `Q` in data2. The objects
  are stored as `fit$data2` and threaded through every re-assembly
  path (`refit()`, `influence()`, `update()`, `drop1()`,
  `anova(refit = TRUE)`, `frm_allfit()`, and `frm_multiple()` per
  imputation), so a `saveRDS()`ed fit refits in a session where the
  environment that built the matrix is gone;
  `tests/testthat/test-data2.R` proves that against the same fit
  built without data2, whose deletion refits all fail.
  `frm_bootstrap()` and the autoscale pre-fit reuse the assembled
  frame and needed nothing. data2 stays optional and the fallback path
  is unchanged, so pre-v0.29 code keeps working.

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

PROBE DONE 2026-09-01, branch wt-odeprobe. Full write-up in
dev/ode-feasibility.md; scripts in dev/ode/. Headline: it already
works, with NO package change. A population PK model (12 subjects,
one-compartment oral absorption, BSV on log ka and log ke) fits
through the plain nl spelling

  bf(conc ~ pk_ode(exp(lka), exp(lke), exp(lV), time, id, dose),
     lka ~ 1 + (1|id), lke ~ 1 + (1|id), lV ~ 1, nl = TRUE)

where pk_ode() is a user function looping RTMBode::ode() once per
subject. The nl body sees row-wise nlpar advectors plus raw data
columns and resolves helpers through the formula environment, which
is all a per-group solve needs (a nlpar constant within group is read
off the group's first row). Verified four ways: numeric vs closed
form 2e-8; identical objective to a hand-rolled MakeADFun; IDENTICAL
logLik to the same model with the analytic mu at n_id 6/12/25/50 (so
the Laplace approximation is right THROUGH the adjoint node); and
agreement with nlmixr2 FOCEi to 3 decimals on all six quantities
(frmtmb finds a marginally better optimum, and needs no C compiler).

Downstream all works: sdreport SEs, wald and profile confint,
predict() in-sample AND on newdata (dense grids, new levels,
re.form = NA - the nl branch re-evaluates the body against newdata's
own columns, so the helper re-derives its groups), simulate() and
refit, REML, ranef/VarCorr/coef. Only predict(se.fit = TRUE) fails,
on the pre-existing "not supported for the nonlinear predictor" nl
gap. Edge cases pass: t = 0 observations, duplicate times, ragged and
row-shuffled designs, NA rows.

Two hard findings that shape the design:

1. NEVER stack subjects into one big system. The Laplace inner
   Hessian through a single ode() node degrades with the state count:
   lsoda goes silently NaN above ~8 states, lsode at 8, and
   ode45/rk4/euler CRASH the R process (exit 127) at 32. adams
   survives to 48 then hangs. Solve one small system per group -
   which is the natural pharmacometric regime anyway. Minimal
   frmtmb-free reproduction in dev/ode/probeE8-laplace-limit.R;
   worth an upstream RTMBode issue. Related: a NaN ODE tape poisons
   later MakeADFun objects in the same R session, so size sweeps and
   tests need one process per case.
2. Fixed-step integrators (rk4, euler) give a WRONG likelihood
   (-63.93 vs -60.46). Every adaptive one agrees to 7 digits; lsoda
   is fastest.

Cost: linear in subjects, ~0.4 s per subject for a 2-state system
(12 subjects 3.5 s, 50 subjects 20 s, all-in including sdreport),
32x-181x the same model with a closed-form mu.

Remaining work is ergonomics only. Recommended: ONE exported helper,
frm_ode(dynamics, init, times, parms, group, method, ...), called
from the nl body, owning the ADoverload("[<-")/("c") locals, the t0
prepend, the within-group sort and scatter, an assertion that every
parms/init column is constant within group (a within-group covariate
is currently unidentified - coefficient pinned at its start, non-PD
Hessian, NaN SEs), a fixed-step refusal, a state-count warning, and
a tryCatch penalty for failed solves. Do NOT build an odefun = or
ode() bf() grammar: it needs new machinery in parse/frame/objective/
predict/simulate for capability the nl body already has, and would
be a worse nlmixr2. Boundary: dosing EVENT tables (multiple doses,
infusions, evid/amt) are out - that needs deSolve events = plumbed
through RTMBode and was not probed. Sizing: 3-4 days total (helper
1-1.5d, tests 0.5d, packaging 0.5d, PK vignette 0.5-1d).

A custom_family() whose lpdf calls ode(), with times/ids/doses in
vreal()/vint(), also works and reaches the identical objective, but
is slower (optimize 5.4 s vs 3.3 s) and its predict() defaults to the
primary dpar's link scale rather than the concentration. Keep it as a
documented alternative for models where the ODE does not enter the
mean.

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

## Hidden Markov models (user request 2026-09-02)

Union-coverage territory: brms has NO HMM grammar. The reference the
user pointed at (cas-bioinf/covid19retrospective, manuscript/hmm.Rmd)
injects a hand-written forward algorithm into brms through `stanvars`
with a custom family - raw Stan, not grammar.

Statistical position: HMM latent states are discrete, so the Laplace
approximation does not apply to them and does not need to. The forward
algorithm marginalizes the state sequence EXACTLY, as a per-group
recursion of small matrix products, which tapes in RTMB. This is
hmmTMB's architecture. Composition with random effects is the v0.19
insight again: sum the discrete states inside (forward algorithm, exact
given b), integrate the Gaussian b outside (Laplace).

PROBE DONE 2026-09-02, branch wt-hmmprobe. Full write-up in
dev/hmm-feasibility.md; scripts in dev/hmm/. Headline: the thesis holds
and rung 1 already works, with NO package change.

Rung 1 - `custom_family()` + `vint(g, t)` payloads, the lpdf running
the forward recursion per sequence and scattering each sequence's log
likelihood onto its first row with a constant sparse indicator - fits a
gaussian HMM with covariate-dependent transitions and per-state random
effects. `check_custom_family()` PASSES unmodified. Validated four
ways: identical to a hand-rolled MakeADFun (4.2e-9 on logLik, 8.4e-6 on
the coefficients); depmixS4 to 2.1e-8 (single sequence) and 1.3e-6 (30
sequences with transition covariates, every coefficient to 4 decimals);
hmmTMB to 8.7e-10 on logLik and 2.8e-7 on the parameters (fixed
effects, stationary initial distribution); and, with `(1|g)` on a state
mean, adaptive Gauss-Hermite quadrature over the per-group scalar b,
which puts the Laplace bias at 0.126 in logLik (8.9e-5 relative) and
4.4e-4 in the parameters - an order of magnitude better than v0.19's
group mixtures. Categorical emissions (the covid model class) work
identically: exact against the numeric forward (1.1e-13) and an
independent BFGS (1.5e-8 on the parameters).

Cost is linear in ROWS and free in the number of sequences: 5000 rows
cost the same cut as 1x5000 or 1000x5 (tape 543 vs 653 ms, fn 0.70 vs
0.88 ms). Tape build is the part that grows slightly faster than linear
and is what becomes annoying above T ~ 20000 (1.9 s per build against
5 ms to evaluate). Log-space via `logspace_add` and Zucchini scaling
both tape and agree to 2.8e-13; scaling is 1.8x faster on the gradient,
log-space is the robust default. The `[<-` gotcha costs only 1.5x on
the tape build here, because the assigned vector is length K, not n.

Three sharp edges, all measured. (1) The post-fit surface SILENTLY
LIES: `fitted()`, `residuals()` and `predict(type = "response")` return
state 1's mean at every row, because `response_mean()` falls back to
`dpars$mu`; no per-row `mean_fn` could be right, since E[y_t] needs the
state-occupancy probability from forward-backward. (2) The cold start
converges to a LOCAL optimum 8.1 logLik below the global one on the
random-effect model, with convergence code 0, max|grad| 3.5e-4, a PD
Hessian and `diagnose()` reporting nothing. (3) REML runs and should
not: it integrates only `mu`'s fixed effects, leaving the transition
and SD dpars outer - a partial restricted likelihood matching no
standard definition. Also: a start with both state means at the median
is a fixed point of the label symmetry and stalls at the one-state
solution; sequences of length 1 reproduce `mixture()` to 3.4e-13 but
report 2 unidentified transition parameters in `df`; NA responses must
be masked with a `vint` flag rather than dropped by `na.omit`, which
silently shortens the chain; and the stationary initial distribution
tapes fine through `RTMB::solve` (4.6e-13 vs numeric) but only when the
transitions are constant.

Found while validating, worth an upstream note: hmmTMB silently treats
a data column named `state` as KNOWN STATES and maximizes the
complete-data likelihood (-1247.25 against -1216.40 on the same data),
with no message. Drop the column before comparing.

Rung 2 (`hmm(K, family, time =, group =, init =)`) is recommended and
sized at 7-9 days in dev/hmm-feasibility.md. Most of it is `mixture()`
machinery verbatim - suffixed per-state dpars with the full formula
grammar, the multinomial-logit weight block (K copies of it, one per
source state), logsumexp, quantile-spread inits, the `mix_g` group
structure, and the sum-inside-integral objective branch. The genuinely
new work is the time/group contract in the frame, the taped forward
recursion, and forward-backward/Viterbi post-processing - the last
being both the largest piece and the reason to build rung 2 at all,
since it is what makes `fitted()` mean something. Refuse initially:
weights/cens/trunc/mi, rescor/mvbf, quadrature, OSA, REML, and
`predict(se.fit = TRUE)` on the response scale.

## t2() newdata prediction defect (found and FIXED 2026-09-02)

`predict(newdata = )` on any `t2()` smooth failed from the day the
prediction path was written until v0.35.1. Fixed; kept here because the
root cause is a genuinely obscure mgcv contract.

Two separate things were wrong.

1. `smooth2random(type = 2)` returns no `trans.U` for a t2 smooth, so
   the s() path's `PredictMat %*% U` multiplied by NULL. The t2 split is
   reported as `pen.ind` instead (penalty index per original basis
   column, 0 = the unpenalized null space): scale `sm$X` by `trans.D`,
   then stable-sort the columns as penalty 1, penalty 2, ..., penalty 0
   and you get `cbind(rand[[1]], ..., Xf)` exactly.
2. That mapping alone was still not enough, because
   `mgcv::PredictMat(sm, data)` did not reproduce `sm$X` on the TRAINING
   rows (max diff 6.2), while for s() it did. ROOT CAUSE: a t2 smooth
   carries TWO identifiability constraints. `smooth.construct.t2.smooth.spec`
   sets `object$C` (the fit constraint, which keeps the penalty blocks
   separable so `pen.ind` means anything) and, separately,
   `object$Cp <- matrix(colSums(X), 1, ncol(X))` (the conventional
   sum-to-zero PREDICTION constraint). `smoothCon(absorb.cons = TRUE)`
   absorbs `C` into `X`/`S` and `Cp` into `Xp`/`Sp`, and `PredictMat`
   honors the `Cp` parameterization. s() has no `Cp`, hence no
   discrepancy there.

   gamm4 hits the same wall and solves it the other way round: it fits in
   the `C` parameterization and then maps the fitted coefficients into
   the `Cp` parameterization before predicting -- `gamm4.r`,
   `object$coefficients <- G$P %*% object$coefficients`, commented "If
   prediction parameterization differs from fit parameterization,
   transform now... (important for t2 smooths, where fit constraint is
   not good for component wise prediction s.e.s)". frmtmb has no
   component-wise smooth s.e. to protect, so it drops `Cp` instead:
   `smoothCon(..., modCon = 3)`, documented as "set fit and predict
   constraint to fit constraint". `modCon >= 3` only does
   `sm$Cp <- NULL`, so `X` and `S` come out bit-identical (diff 0) and
   the fit, including logLik, does not move.

FIX (v0.35.1): `R/frame.R` builds smooths with `modCon = 3` and stores
the pen.ind permutation on `sm_info` via `smooth_pen_order()`, which
verifies the scale-and-permute identity on the training rows rather than
asserting it. `R/predict.R` takes the t2 branch
`sweep(PredictMat, 2, D, "*")[, ord]`. The refusal that remains is
narrow: only if `smooth_pen_order()` returns NULL, i.e. some future
smooth class reaches the no-`trans.U` path with a `pen.ind` that does
not mean what t2's means.

## Robustness items (user discussion 2026-09-02)

Four distinct items under the word "robust". Items 1 and 3 are
delivered (below); items 2 and 4 remain queued.

1. NUMERICAL robustness audit - DELIVERED, see the next section.
2. ~~t-distributed random effects~~ DELIVERED as
   gr(g, dist = "student", dist_nu = 5); probes in dev/tre/ and
   dev/tre-feasibility.md (Laplace bias = per-group constant, argmax
   unmoved; quadrature = TRUE exact escape hatch; nu not ML-estimable,
   fixed 5 by minimax; brms's diag t-block is still MVT).
3. ~~huber() family~~ DELIVERED, see the numerical-robustness section
   below (k fixed at 1.345, working-likelihood caveats documented).
4. ~~Cluster-robust (sandwich) vcov~~ DELIVERED, `dev/sandwich.md`.
   The claim held: per-CLUSTER scores are exact and cheap, from one
   tape carrying a per-cluster weight parameter (`obj$gr` at a cluster
   indicator IS that cluster's score, because a data-free cluster
   integrates its Gaussian prior to exactly 1). CR0/CR1/CR1p/CR1S ship
   in the clubSandwich spelling; CR2/CR3 are refused, because
   Bell-McCaffrey is defined through a hat matrix a Laplace-marginal
   likelihood does not have.
robustlmm-style bounded-influence estimating equations (DAStau) are
NOT a likelihood and stay out of scope; item 2 is our answer to the
same concern.

## Numerical robustness audit (item 1, delivered)

The concern was right and the damage was wider than the two RTMB
robust densities cover. RTMB 1.9 exports exactly `dbinom_robust(x,
size, logit_p)` and `dnbinom_robust(x, log_mu, log_var_minus_mu)`;
RTMBdist ships no `_robust` form at all. Everything else had to be
rewritten in log space by hand.

What the probe found at the BASE commit (per-observation nll taped as a
function of one dpar's linear predictor; value finiteness plus gradient
against a central difference):

- every logit-scale dpar - `mu` for bernoulli / binomial /
  beta_binomial / beta / zero_inflated_binomial, and every `zi`, `hu`
  and `quantile` gate - was already wrong in the second decimal at
  eta = +30 and returned -Inf (NaN gradient) at +40. Asymmetric,
  because only `1 - plogis(eta)` cancels.
- `cloglog` was far worse: `1 - exp(-exp(eta))` is exactly 1 from
  eta = 4, so a cloglog binomial had NO gradient past a single-digit
  linear predictor.
- negbinomial and geometric gave NaN at eta = -40: `dnbinom2` forms
  `p = mu / var`, which is 0 / 0 once exp(eta) underflows.
- cumulative and sratio died at eta = -40, cratio at +40 (mirrored
  floating-point failure of what is, for the logit link, the same
  density); cumulative with near-coincident thresholds died at -30.
- cumulative(probit) was the worst survivor: a difference of two
  saturated `pnorm` values is 0 from |eta| = 8.
- acat / categorical / multinomial summed `exp(eta)` and overflowed at
  |eta| = 709 / (K - 1).

The eta-scale design. The dpar contract hands the lpdf the
INVERSE-LINKED value, which is what every numeric post-fit path wants,
and recomputing eta with `qlogis()` inside the lpdf would just
reintroduce the round trip. So `build_objective()` now stores the
linear predictor beside each dpar under a reserved `.eta_<dpar>` name
(the `.cs` precedent), and the accessors `robust_logit()`,
`robust_logmu()`, `log_dpar()`, `gate_logs()` and `mu_pair()` in
R/families.R read it. Off the tape the entries are simply absent, every
accessor returns NULL, and the family falls back to the plain form -
correct there, because nothing off the tape is differentiated. No
family contract changed and no post-fit path was touched. Whether a
robust form EXISTS is a property of the link, so R/links.R grew
optional `logit_eta` (logit, cloglog) and `log_eta` (log) fields; a
link without one leaves the family on the plain path automatically,
which is why `identity`-link binomials still work.

Why not hand-roll a log-space Poisson while we were there: `RTMB::dpois`
dispatches on `osa` and `simref` objects, so replacing it with plain
arithmetic would silently break `oneStepPredict()` residuals and
tape-side `simulate()`. The RTMB `_robust` densities keep that dispatch,
which is the other reason to prefer them. dpois is finite and correct
through |eta| = 700 anyway.

Known remaining limits, all documented and pinned by tests:
- `cumulative(probit)` keeps the plain CDF difference. There is no
  exact log-space form: `RTMB::pnorm` carries no `log.p`, and the
  mirrored `pnorm(-b) - pnorm(-a)` trick needs a branch on the sign of
  a parameter-dependent quantity. Use the logit link.
- `nbinom1` at |eta| > 709: its negative-binomial size is `mu / phi`,
  so the size itself underflows with the mean. Intrinsic to the
  parameterization, unchanged by this work.
- The ordinal families' `osa` re-tape branches keep the plain forms.
  `oneStepPredict()` runs at the fitted values, where they are fine.

## huber() (item 3, delivered)

`k` is a fixed argument of the constructor, `huber(k = 1.345)`, not a
dpar. It states where the analyst draws the line between a residual and
an outlier, which is a modelling choice; `MASS::rlm()` treats it the
same way, and estimating it would let the likelihood buy fit by
widening the gaussian core. The normalizer
`sqrt(2 pi) (2 Phi(k) - 1) + (2 / k) exp(-k^2 / 2)` agrees with a
numeric integral to 1.8e-14 and the density integrates to one to the
same order. `rho` is written branch-free as
`(u^2 - max(|u| - k, 0)^2) / 2` with `max(x, 0) = (x + |x|) / 2`,
because the regime switch has to happen on the tape.

`rlm()` fixes the scale at a MAD-type estimate and iterates the
location; `huber()` estimates `sigma` by ML jointly with `mu`. That is
the whole difference: coefficients agree to 5e-3 .. 9e-3 on ordinary
data (7e-2 under 5% gross contamination, where the ML scale doubles),
and to 6e-6 .. 1.3e-5 when `sigma` is held at `rlm`'s own scale with
`bf(y ~ x, sigma = rl$s)`. `rho`'s kink stops any optimizer around a
maximum absolute gradient of 1e-4, so the sharp check is Huber's own
estimating equations, `X' psi(u) = 0` and `sum(u psi(u)) = n`, which
the fit satisfies to the same order.
## `nlf()`: recorded, not implemented (v0.34.x, wt-polish)
## `nlf()`: DONE (v0.36, wt-nlf)

The v0.34 note was right about the shape and wrong about the cost. The
work was not "a nonlinearity flag per dpar plus a nested body in the
objective" bolted beside the existing branch: it was deleting the
branch. `nl = TRUE` is now sugar - "mu carries a body, and the body is
the response formula's right-hand side" - and `parse_one_response()`
has ONE path that builds every dpar, with `nl_body` / `nl_pars` on the
dpar entries that have one and a dependency sort over the whole set.
Everything downstream was already keyed on `lp$nl_body` rather than on
"the model is nonlinear", so the generalization touched five call sites
(evaluation scope in `objective.R`, `eval_dpars()` and `predict()`,
`ce_plot_vars()`'s recursion, `check_ode_constancy()`'s body walk,
`drop_nl_lexical_datavars()`'s loop) and no algorithm.

Scope delivered:

- `nlf(dpar ~ body)` on ANY dpar of any family, including the ones
  `nl = TRUE` could never reach: a nonlinear `sigma` (or `shape`, or a
  mixture's `theta`) with a linear `mu`.
- Bodies chain to any depth, computed in dependency order (a stable
  topological sort that keeps declaration order among independent
  dpars, which is why the existing nl coefficient ordering is
  unchanged). Cycles are refused by name, as brms refuses them.
- The composition identity `bf(y ~ a) + nlf(a ~ exp(b * x)) +
  lf(b ~ 1)` == `bf(y ~ exp(b * x), b ~ 1, nl = TRUE)` holds exactly
  (logLik, coefficients, vcov and dimnames, predictions all at 0).
- `frm_ode()` composes inside an `nlf()` body with no change: the
  within-group constancy check now walks every body, not just mu's.
- Two documented supersets of brms. (1) brms insists on `nl = TRUE`
  when the response formula names a nonlinear parameter; here `nlf()`
  has already declared the name, so the flag is optional. (2) A body
  may name another dpar of the same response and read its per-row
  VALUE; in brms such a name is a data column and the model is refused
  when no column has it. A column still wins here, so no brms body
  changes meaning.

The dpar reference closes the one hole in the `varFunc` translation
table: `varPower(form = ~ fitted(.))` is
`nlf(sigma ~ ls + th * log(abs(mu))) + lf(ls ~ 1, th ~ 1)`. The
likelihoods are identical (frmtmb's objective at `gls`'s own estimates
reproduces `logLik(gls)` to 9e-10); the ESTIMATORS are not. `gls`
alternates between the mean fit and the variance function against
iteratively updated fitted values and stops without a gradient check,
so on the reference data it lands 1.1e-4 lower in logLik with a
gradient of 0.09, where the joint maximization here reaches 6.6e-4.
Expect agreement to about 1e-2 on the coefficients, not 1e-10 - unlike
`varPower(~v)`, `varExp`, `varIdent` and `varFixed`, whose variance
functions do not depend on the fitted values and therefore agree
exactly. The migration vignette's variance-function section (added on
main after this branch) should replace its "no frmtmb spelling"
sentence for `varPower(form = ~ fitted(.))` with this one, and carry
the estimator caveat.

Not in scope, and deliberately: `nlf(resp =)` and `nlf(dpar =)` (both
deprecated in brms; placement after the right `bf()` says the same
thing), `flist =`, and `loop =` (accepted and ignored - frmtmb always
evaluates a body once over whole vectors, which is brms's
`loop = FALSE`, and an elementwise body has the same value either way).
An `nlf()` left-hand side naming several parameters is refused rather
than split: one shared body makes them the same function of the data,
which is an aliased model. brms refuses it too.

## Multi-membership `mm()` / `mmc()` (delivered v0.35)

Was absent entirely, which the v0.34 audit flagged: a ported brms
formula `(1 | mm(g1, g2))` died with "could not find function mm",
several stages away from the grammar it belongs to.

The design claim is narrow and it held. `mm()` changes the Z MATRIX
and nothing else: the block is an ordinary `us`/`diag` block over the
pooled level set, and each observation's design row puts weight `w_j`
on member `j`'s columns instead of a 1 on one column. So the registry
nll, the Laplace step, `ranef()`, `VarCorr()`, `ngrps()`, `simulate()`
and `residuals()` needed no multi-membership branch at all; the work
was the parse (`parse_mm_call()`, `split_mmc_lhs()`), the Z builder
(`mm_pooled_levels()`, `mm_index_weights()`, `mm_member_designs()`,
`mm_local_Z()`), and one newdata hook in `predict.R` that emits one
`re_parts` entry per member with its design pre-scaled by that
member's weight, so `re_eta()` and `re_design_matrix()` accumulate the
weighted sum with no branch of their own.

brms semantics reproduced exactly (verified against
`make_standata()`'s `J_*`, `W_*`, `Z_*` arrays, 0 difference on five
formulas): levels are the members' own level sets concatenated in
written order and deduplicated (`frame_re()`); the default weights are
`1/J` and are NOT rescaled; `scale = TRUE` divides a SUPPLIED weight
matrix by its row sums (`data_gr_local()`); `mmc(x1, x2)` is ONE
coefficient whose covariate is member specific.

Refused, and why: any covariance structure but `us`/`diag` (the pooled
levels have no order, coordinates, or relationship matrix for one to
be defined on), `gr(mm(...), cov =)`, an `|ID|` key over an `mm` term
(a merged block indexes one level set per row), and brms's `by =`,
`cor =`, `id =`, `pw =`, `cov =`, `dist =` (each refusal names the
spelling that replaces it: `cor = FALSE` is `diag(x | mm(...))`).

One thing the lane had to reconcile rather than add. The two standard
error paths grouped new-level variance differently, and multi-membership
is the first term that makes one block contribute SEVERAL entries to one
row, so the disagreement became reachable: the linear-predictor path
added one quadratic form per entry (`w1^2 S + w2^2 S`), the
expected-response path keyed on the BLOCK and summed the entries first
(`(w1 + w2)^2 S`). Neither is right on its own - which one applies
depends on whether the row's two members name the same unseen level
(one draw, perfectly correlated) or two different ones (two independent
draws). Both now go through `extra_var_blocks()`, which keys on
(block, level) and sums design rows inside a key before taking one
quadratic form. Single-membership terms key on a constant, so they keep
the one-level-per-row grouping; the linear-predictor path additionally
gains the `|ID|` cross-covariance terms it was dropping, which the
expected-response path already had.

Left open. `mm() x quadrature` is declared untested rather than
claimed: an `(1 | mm(g1, g2))` block passes the scalar-intercept guard,
but nothing checks that the Gauss-Kronrod rule marginalizes a design
whose rows load several levels at once. `getME("flist")` skips mm
blocks, because lme4's flist has no representation for a row that
belongs to several levels. `by =` (level-specific variances) and
`pw =` are the two brms arguments with no analog anywhere in the
package yet, and `by =` would be the next one worth having.
## Three family gaps closed (v0.35, 2026-09-02)

`categorical()`, `von_mises()` and `cox()` were the brms families with
no frmtmb spelling. All three shipped. What each cost, and what stayed
out.

### categorical()

The likelihood was already here as `multinomial(K)` on a count-matrix
response; what was missing was the brms spelling on a bare factor.
Shipped as a `type = "categorical"` family over a `1..K` category
index: `K - 1` linear predictors, category 1 the reference, dpars named
`mu<Level>` as brms names them (confirmed from
`make_stancode(bf(y ~ x), family = categorical())`, which emits `b_mub`
/ `b_muc` for levels `a < b < c` and calls
`categorical_logit_glm_lpmf`). The main formula applies to every
category unless a dpar formula overrides one, again as brms does.
Validated to 4e-10 in the log-likelihood against `nnet::multinom` and
to 1e-8 against `multinomial(K = 3)` on the one-hot response.

The one place brms and frmtmb genuinely differ is WHEN the categories
are known. brms builds Stan code after seeing data; `parse_spec()` here
is deliberately data-free, so a dpar named `muc` would be rejected as
unknown before any response is read. Resolved with a `defer` hook on
the family and one line in `frm()` (mirrored in `get_prior()` and
`frm_simulate()`): a bare `categorical()` is a placeholder that
`resolve_deferred_families()` swaps for the concrete family once the
response is in hand. `categorical(levels =)` and `categorical(K =)`
build it directly for the paths that have no data.

Not done: `conditional_effects()` and `dharma_residuals()` branch on
`type == "ordinal"` and take the wrong path for a nominal response.
Both live outside this lane; the compat registry leaves them at the
untested default rather than claiming they work. `residuals()` is
refused outright - a nominal response has no scale for one.

### von_mises()

Nearly free, and the Bessel blocker in the brief turned out not to
exist. `RTMB` exports an S4 `besselI` method for advectors:
`MakeTape(function(x) besselI(x, 0))` tapes and differentiates, and
`RTMBdist::dvm()` already builds the von Mises density on it in the
exponentially scaled form (`log I0(k) = log besselI(k, 0,
expon.scaled = TRUE) + k`, which is what keeps a large concentration
from overflowing). So no series approximation, no accuracy caveat, and
no refusal: the family is three lines of lpdf plus the `tan_half` link
(`linkfun = tan(mu/2)`, `linkinv = 2 atan(eta)`), brms's
parameterization exactly.

The one thing that had to be written was the SIMULATOR.
`RTMBdist::rvm()` delegates to `circular::rvonmises()`, which takes
scalar parameters only, so a distributional `kappa ~ x` could not
simulate. Replaced with a vectorized Best-Fisher (1979) rejection
sampler over per-row `mu` and `kappa`, with the `kappa -> 0` uniform
case split out (the rejection constants divide by it).

Validated against a hand-rolled RTMB likelihood to 1e-6 (both for
constant and for distributional kappa) and against
`circular::mle.vonmises()` to 1e-5 on mu and 1e-2 on kappa (its kappa
comes from a different root finder). `residuals()` runs but the
differences are NOT wrapped; `residuals(type = "osa")` is refused
upstream by `dvm()` itself.

### cox()

brms's exact construction, read off `make_stancode()` for
`bf(t | cens(c) ~ x + (1|g))` with `family = cox()`:

    bhaz  = Zbhaz  * sbhaz      // M-spline basis, simplex weights
    cbhaz = Zcbhaz * sbhaz      // I-spline basis, the SAME weights
    cox_log_lpdf  = log(bhaz) + eta - cbhaz * exp(eta)
    cox_log_lccdf =                  - cbhaz * exp(eta)

with `simplex[Kbhaz] sbhaz`, a `dirichlet(1)` prior, and the basis from
`brms:::bhaz_basis_matrix()`: `splines2::mSpline` / `iSpline`, default
`bhaz(df = 5, intercept = TRUE)` (cubic), internal knots on response
quantiles and boundary knots at
`c(max(min(y) - diff(range(y))/50, 0), max(y) + diff(range(y))/50)`.

Reproduced without taking `splines2` as a dependency: the M-spline and
I-spline bases are built here from a hand-rolled Cox-de Boor recursion
(`bspline_basis()`), with `I_j` the reverse cumulative sum of the
order-(degree + 2) B-spline basis on the once-more-repeated boundary
knots. They agree with `splines2` to 1e-16 in both the intercept and
no-intercept cases, and the tests also check `I_j` against a numeric
integral of `M_j` and that each M-spline integrates to one. The simplex
rides in the parameter vector as `sbhaz_raw`, its `Kbhaz - 1` softmax
coordinates with the first pinned at zero.

Two small contract changes carried it:

- `family$aterm_data(y, aterms)`, an optional hook in `assemble_frame()`
  returning family-level DATA that no addition term supplies. The
  spline bases are a function of the validated response, so they are
  built once there and ride with the addition-term values the objective
  already bakes into the tape. Nothing else uses the hook yet.
- `fam_lcdf()`, an arity dispatcher, because the Cox survivor function
  needs the family-level extra parameters and the three-argument
  `lcdf(q, dpars, aterms)` contract had no room for them. Every other
  family keeps the old signature.

`cens()` then needed no new machinery at all: the objective already
replaces a censored row's density with a windowed CDF difference, which
for this family is exactly the survivor function. Right, left and
interval censoring all run. Frailty models come free through Laplace,
which is the point of the family: `time | cens(c) ~ x + (1 | g)`
recovers a frailty SD of 0.8 within 0.3 and agrees with
`coxph(... + frailty(g, distribution = "gaussian"))` to 0.1 on the
coefficient.

Validated exactly (1e-6 in the log-likelihood, 1e-4 in the
coefficients and the baseline simplex) against a hand-rolled M-spline
PH likelihood built independently in the test, and approximately
(2e-2) against `survival::coxph`, which is the right claim: coxph
leaves the baseline fully nonparametric while this spends `df` spline
weights on it.

KNOWN AND DOCUMENTED: ML routinely drives one or more baseline weights
to the simplex boundary, so their softmax coordinates run to minus
infinity along a flat ridge and the optimizer reports singular
convergence. The gradient is zero there and the regression
coefficients are at their optimum - a test asserts both - but the
warning is real and brms does not meet it, because its Dirichlet prior
keeps the weights interior. Lowering `df` is the remedy. A penalized
or bounded baseline would remove the warning; it would also stop being
maximum likelihood, so it was not done.

Refused, deliberately: `fitted()` and `predict(type = "response")` (a
survival time has no mean the censored rows identify - brms refuses the
same question, `posterior_epred` has no cox method) and `simulate()`
(no quantile function for the cumulative baseline). `trunc(lb = )` runs
as delayed entry through the same log-CDF but is declared conditional,
not verified: there is no external left-truncated reference in the
suite yet.
DONE-rung-2, branch wt-hmm. `hmm(K, family, time =, group =, init =,
trans =)` is a first-class family in R/hmm.R, wired through additive
hooks in parse/frame/objective/predict. Everything the probe predicted
transferred: suffixed per-state dpars with the full formula grammar,
the multinomial-logit block (K copies), the logspace_add fold,
quantile-spread inits, and a third objective branch beside `rescor`
and `mix_g`.

What was built differently from the sketch. The recursion is sliced by
TIME rather than by sequence: sorting the sequences by decreasing
length makes the set still running at step s a prefix of that order, so
one step is K^2 vectorized `logspace_add` calls over all live sequences
at once. The tape then holds O(Tmax K^2) nodes rather than O(n K^2),
and the sparse scatter the rung-1 recipe needed disappears entirely -
the objective branch sums the per-sequence values directly, because
only their total is wanted. The same slicing carries the numeric
forward-backward and Viterbi passes.

Validation. depmixS4 on gaussian emissions with constant transitions
(3.2e-9) and with `trans = ~x` (1.3e-8, every transition coefficient to
five decimals, which pins the state-1 reference and the "covariate at t
drives t -> t+1" convention); poisson (1.3e-9) and categorical
(multinomial, 6.2e-9, emission probabilities to five decimals). hmmTMB
on the stationary fixed-effect model, 1.1e-12 on logLik and seven
decimals on every parameter. The probe's own reference numbers
reproduce exactly: probeA1 -235.204362981 (2.9e-9), probeB1
-891.019360499 (8.3e-9), F3's free/stationary pair -1608.76549264 (df
7) and -1609.41007570 (df 6) to the digit, F4's masked-NA
-691.400711096 to the digit, the stationary solve 6.8e-13 against the
numeric stationary forward, and forward-backward 6.7e-16 against a
brute-force per-sequence reference with Viterbi identical.

Residue, all of it deliberate and documented in ?hmm.

- The T = 1 collapse is refused rather than allowed: with every
  sequence a singleton the transition logits are flat and the reported
  df counts them (probe F2's 7 against a mixture's 5). Holding them at
  constants lifts the refusal, and the fit then reproduces `mixture()`
  to 6.8e-13 with df 5 on both sides, which is the shipped test.
- Multimodality is unfixed, as the probe said it would be. The cold
  start on probe D4's model still lands at -1096.09575602 with
  convergence 0; restarted at the global point it reaches
  -1087.99646521, hmmTMB's value to the last digit. `frm_allfit()`
  varies optimizers, not starts, so a multi-START helper is the
  obvious follow-on.
- K is capped at 9 because the `tr{i}{j}` dpar names concatenate the
  state indices. A separator would lift the cap at the cost of the
  naming that reads naturally.
- `multinomial(K = C)` emissions inherit its dpar names, so the
  per-state copies are `mu21`, `mu31`, ... (category, then state).
  Legible enough at C, K <= 9, ugly beyond.
- `predict(type = "response")` on newdata, `se.fit` on the response
  scale, `conditional_effects()`, OSA and deviance residuals, and
  `emmeans` are refused or untested: each needs the occupancy
  probability, which conditions on the observed responses of a whole
  sequence.
- Continuous-time transitions, higher-order chains, and hmm() inside
  hmm() stay out, as the memo's section 8 said.
- No vignette yet (the memo budgeted a day for one).

## Pharmacometrics tier (rxode2/nlmixr2 comparison, user request 2026-09-02)

What replicates today: smooth-ODE population models via frm_ode()
(dynamics as plain R), REs/covariates on parameters, bolus/replace/
multiply/reset dosing, constant-rate infusions, repeated doses,
ii/addl compact repetition, steady-state records, piecewise-constant
time-varying covariates, per-row output selection, event_scale
bioavailability, frm_simulate population simulation. Analytic solved
systems (linCmt analogs) are writable by hand as nl bodies.

Status of the six items (closed 2026-09-02, branch wt-odegaps; see
dev/ode-feasibility.md section 10 for numbers):
1. TIME-VARYING COVARIATES - **DONE**, piecewise-constant only.
   frm_ode(tv = , tv_break = ) rides the segmented solve; LOCF, the
   rxode2 convention for covariates. Estimated VALUES allowed
   (tv_break carries the change points as data); rxode2's linear
   interpolation is NOT available and cannot be (the value inside a
   segment would have to depend on t, and t is an advector).
2. STEADY-STATE dosing - **DONE** as an approximation.
   events$ss with events$ii, run-in of n_ss cycles (default 20) from
   a zeroed system. Iteration-until-convergence is impossible on the
   tape (it branches on a value), so n_ss is data and the shortfall
   is exp(-n_ss * k * ii) for linear kinetics; the numeric path
   compares the last two cycles and warns past ss_tol.
   Reset events - **DONE**. method = "reset" sets every state to
   value, which is NONMEM/rxode2 evid 3 at value = 0; reset + add at
   one instant is evid 4, reset first.
3. Lag times (alag) - still refused by design (estimated event times
   change tape structure); revisit only with a new design.
4. Combined error models - **DONE** as a documented example, no
   machinery: nlf(sigma ~ 0.5 * log(exp(2*ladd) + exp(2*lprop)*mu^2))
   + lf(ladd ~ 1, lprop ~ 1). In ?frm_ode and vignette("ode").
5. Multiple-endpoint ODE models - **DONE** via output-by-row.
   frm_ode(output = <column>) shares one solve and gathers; the
   two-call indicator spelling gives an identical logLik. mvbf was
   not needed and was not tried.
6. ii/addl sugar - **DONE**, pure preprocessing, expansion is
   bit-identical to the hand-written table.

Left open after this round:
- Estimated event times / lag times / estimated ii (as above).
- rxode2-style linear covariate interpolation.
- A NONMEM-shaped reader (one table with evid/amt/cmt/rate split
  into data + events). Deliberate: it is a data-reshaping helper,
  not model machinery, and the two-table split is the honest shape.
- More than one ss row per group (a second run-in would discard the
  first); write the later doses out with ii/addl.

## Next-round small items (queued 2026-09-02, DELIVERED same day, branch wt-smallitems)

All three shipped in v0.37: plot.frmtmb_influence (Cook's + dfbetas
panels, 2/sqrt(n) band, car::influencePlot refusal documented),
the one-time shadowing message with leading-dot aliases (.sigma,
registered only on collision), and bare-nlpar bound aliases. One
probe claim below was corrected in the process: nlpar-vs-dpar names
do NOT collide-and-error - a name in an nl body that matches a dpar
is a dpar REFERENCE by construction (nlpars is setdiff'd against
fam$dpars), and vignette("inputs") now records the actual precedence
(nlpar, then data column, then dpar). Original queue entry kept for
the record:

- plot.frmtmb_influence: index plot of Cook's distances with top
  cases labeled + dfbetas panels, from the refit-based quantities
  influence() already computes exactly. Deliberately NOT
  car::influencePlot compatibility: hatvalues/rstudent are
  OLS-geometry approximations that are ill-defined for a
  marginalized Laplace fit, and our refit-based deletion is the
  exact version of what those approximate. Afternoon-scale.
- Reserved-name shadowing disambiguation: a covariate literally
  named sigma fits fine (names namespace apart: mu.sigma vs
  sigma.(Intercept)), but hypothesis(fit, "sigma = 0") silently
  resolves to the COVARIATE coefficient via the v0.21 shadowing
  guard. Add a one-time message when shadowing occurs in the
  hypothesis env, and a naming-collisions paragraph in
  vignette("inputs") (probe results in the 2026-09-02 session:
  nlpar/dpar collisions error; undeclared dpars error; .eta_
  reserved names are unreachable from data; nlf bodies use
  column-wins deliberately).

## frm_sample on ODE fits: BROKEN at the tmbstan/RTMBode boundary
## (found 2026-09-02, deterministic repro)

Sampling a fit whose tape contains frm_ode() adjoint nodes fails at
warmup iteration 1 EVEN AT THE FITTED OPTIMUM, where direct
obj$fn/obj$gr succeed: DLSODA reports TOLSF = nan / T-illegal-at-
TCUR=0 signatures, then RTMB's ADjoint raises "Wrong output length".
Reproduced with raw tmbstan(fit$obj, init = list(mode)) - NOT an
frmtmb init-handling issue. Repro: scratchpad lv-sample-debug.R /
lv-init-probe.R (the Lotka-Volterra two-solve model; whether
single-node ODE fits also fail is untested). Prime suspect: RTMBode's
C shim keeps STATIC global tape pointers (set_pointers X/Y/F, set by
setTape before each solve) and tmbstan's evaluation context breaks
the sequencing that direct evaluation happens to satisfy. Secondary
damage: the R error crossing Stan's C++ boundary corrupts rstan's
nested AD arena for the WHOLE SESSION (later calls fail with
"empty_nested() must be true before calling recover_memory()") -
retry needs a fresh session.

frmtmb side (done): frm_sample now raises an informative error when
the sampler returns no draws, naming this case, instead of a
seq_len(NA) error from the extraction code.

ROOT CAUSE FOUND (2026-09-02, branch wt-rtmbode; dev/upstream/ holds
the full analysis, a three-patch series against kaskr/RTMB@5242257,
13 repro scripts and issue/PR skeletons). NOT static pointers, NOT
tmbstan, NOT gc - all three eliminated experimentally. RTMBode's
updateSolution() calls deSolve unguarded; deSolve fails two ways
that both escape as R errors (a hard "illegal input" error when a
rate reaches Inf through exp(), and an EARLY RETURN with fewer rows
than requested, which the fixed-length ADjoint node reports as
"Wrong output length"). Stan reaches such parameters by construction
(unconstrained initial step size 1 against a gradient of order 1e2+
at the mode), so warmup iteration 1 fails deterministically; both
signatures reproduce from a plain obj$fn() call with no sampler.
Patch 01 (NaN fill + mapping the solution onto the requested times)
makes the LV model sample end to end, 0 divergences. Bonus finding:
the >8-state ceiling is a deSolve integer overflow (lsoda's lrw
formula overflows R integers at neq = 46337; order-3 sensitivities
for a Laplace gradient get there fast); separate deSolve issue
skeleton in the write-up.

Our side, future: auto-derived containment bounds (mode +/- k*SE)
for sampling ODE fits once a patched RTMBode makes failure
rejectable. Bare-nlpar bound aliases: DONE in v0.37 (wt-smallitems).

## Functional data analysis: probed 2026-09-02 (lane wt-vignettes)

Probed while writing sections 10 to 12 of `vignette("case-studies")`.
Verdicts, with the numbers that produced them.

**Function-on-scalar regression: WORKS, no gap in the model.** Long
format plus `s(t) + s(t, by = x)` is the mgcv reduction, and frmtmb
reproduces `mgcv::gam(method = "ML")` on the same data: the fitted
coefficient functions agree to 2e-6, the smoothing parameters to 9e-5
relative, and frmtmb's marginal `logLik()` equals mgcv's ML
smoothness-selection score (`-gam$gcv.ubre`) to 1e-8. Note for anyone
comparing again: `logLik.gam` is NOT that quantity (it is the
unpenalized likelihood at the fit), so comparing `logLik()` to
`logLik()` across the two packages looks like a 30-point disagreement
and is not one.

**`bs = "fs"` factor-smooth interactions: WORK, undocumented until
now.** `s(t, subject, bs = "fs", k = 5)` survives
`smoothCon()` + `smooth2random()` and comes out as three variance
components (the wiggly part plus the two null-space directions).
Agreement with `mgcv::gam` on the same term: coefficient functions to
2.4e-6, marginal logLik to 1.4e-10. This is the honest FoSR spelling
for a per-subject curve, and it beats `(1 + t | subject)` on the same
parameter count (AIC 899.86 against 901.65 on the vignette's data).
Nothing in `?frm` or `vignette("inputs")` said `fs` was available.

**Scalar-on-function regression: WORKS.** `s(Smat, by = LX)` with
MATRIX columns in the data frame survives frame assembly, fitting and
`predict(newdata =)` with a matrix column. Fitted values match
`mgcv::gam` to 9e-8 and the recovered coefficient function to four
decimals. The earlier expectation that this would fail was wrong.

### The gaps that are real

- **CLOSED 2026-09-02 (lane wt-smoothgaps).**
  `conditional_effects()` on a matrix covariate said the wrong thing.
  The refusal now names the matrix columns and the smooth term that
  carries them, says a matrix column is a whole function per row and so
  has no one-dimensional axis to vary along, and points at
  `predict(newdata = )` over a grid the user builds. The generic "No
  plottable predictors found" stays for the genuinely empty case, so
  the two faults read apart.
- **CLOSED 2026-09-02 (lane wt-smoothgaps).** `re.form = NA` is now
  population-level for smooths too. Probed first: the old behavior
  dropped `(1 | subject)` and kept BOTH the population smooth and the
  `fs` per-subject curves, which is
  `mgcv::predict.gam(exclude = "s(subject)")`, not the population
  curve. The rule now implemented classifies each smooth from its
  smoothCon object: group-indexed (`fs.interaction`, `random.effect`
  over a factor, a tensor product with an `re` margin) is dropped,
  everything else (`s()`, `s(by = )`, `te()`, `t2()`, `sz`, `gp()`,
  `hsgp()`) is kept. It agrees with
  `predict(gam, exclude = c("s(t,subject)", "s(subject)"))` to 3e-9 in
  sample and on newdata. `conditional_effects()` draws at the
  population level, so it inherits the rule, and the grouping factor of
  a factor-smooth is no longer offered as an effect to plot.
- **CLOSED 2026-09-02 (lane wt-smoothgaps).** `predict(newdata = )` no
  longer needs the grouping column for an `fs` term under
  `re.form = NA`: the dropped block is never built, so `PredictMat()`
  is never called for it. When the column IS needed (a conditional
  prediction), a frmtmb error names the column and the term before
  mgcv's internal message. An unseen level of an `fs` term now errors
  by default, naming the term, and `allow_new_levels = TRUE` predicts
  it at the population level, which is what mgcv's basis already
  produced: it matches levels by label and returns a zero row for one
  it does not know.
- **No functional-response object.** There is no `refund::pffr`
  equivalent: no `ff()` term, no automatic long-format reshape, and no
  penalty that couples `b0(t)` with `b1(t)`. This is a genuine
  difference in scope, not a defect, and the long-format spelling
  covers the models the short course teaches.

### A custom link with a data-dependent bound: works, thinly documented

`frmtmb_family(links = )` accepts a link OBJECT as well as a name
(`get_link()` passes a list through), so a scaled logit onto `(0, U)`
with `U` read off the data at family-construction time is expressible
and is what the vignette's Wiener family uses for `ndt`, whose support
constraint is `rt > ndt`. The contract a custom link must satisfy is
`name`, `linkfun`, `linkinv`, `mu_eta`, and it is stated only in the
`@param links` line of `?frmtmb_family`. Worth an example in the
custom-family docs. Note also that such a link saturates in double
precision past a linear predictor near 37, so the bound is exact in
algebra and only nearly exact in arithmetic.

## Prior-machinery edge cases (review findings, 2026-09-02, logged not fixed)

Three low-severity behaviors in R/priors.R, confirmed by reading
during the v0.39 review, left as-is deliberately (the machinery was
just reviewed end to end; each is an edge case with a safe default):

1. An explicit `class = "cor"` prior silently skips a refused block
   (for example toep) when another block matched; the refusal text
   surfaces only when nothing matched. A per-block note would be
   kinder.
2. A bounds-only `class = "theta"` spec on a position already covered
   by a joint LKJ entry adds a hard bound without retiring or
   renormalizing the LKJ term: an undocumented truncation. Document
   or refuse.
3. The cs marginal map with eta < 1 evaluates (eta - 1) * log(0) at
   logistic saturation (|t| >= ~745) and returns NaN on the tape;
   unreachable in normal sampling, but a clamp would close it.
