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

## t2() newdata prediction defect (found 2026-09-02, refused not fixed)

`predict(newdata = )` on any `t2()` smooth has failed since the
prediction path was written: `smooth2random(type = 2)` returns no
`trans.U` for t2 (the s() path's `PredictMat %*% U * D` multiplies by
NULL). Now an informative refusal; in-sample fitted()/predict() are
unaffected (they use the fit-time design matrices).

What the investigation established (probe scripts in scratchpad,
findings preserved here):
- The correct column mapping is NOT identity: `smooth2random.t2`
  returns `pen.ind` (penalty index per original basis column, 0 =
  fixed), and `sweep(sm$X, 2, trans.D, "*")` followed by grouping
  columns as (pen 1, pen 2, ..., pen 0) reproduces
  `cbind(rand[[1]], rand[[2]], ..., Xf)` EXACTLY (diff 0).
- BUT `mgcv::PredictMat(sm, data)` does not reproduce `sm$X` on the
  TRAINING rows for a t2 smooth built with
  `smoothCon(absorb.cons = TRUE)` (max diff 5.9 on the probe), while
  it does for s(). So the newdata basis needs more than the pen.ind
  mapping; the discrepancy (likely a t2-specific reparameterization
  smoothCon applies to `sm$X` that PredictMat does not, or a
  constraint-absorption ordering issue) must be understood before a
  fix ships. Compare how gamm4 predicts t2 terms.
- Fix sketch: store `ord = order-by-pen.ind` and `trans.D` on sm_info
  at frame time, resolve the PredictMat basis question, then newdata
  M = PredictMat-equivalent basis, scaled and permuted. Validate
  newdata == in-sample rows at 1e-10 and against
  mgcv::predict.gam on the same model.

## Robustness items (user discussion 2026-09-02, HELD for later)

Four distinct items under the word "robust", none queued yet:

1. NUMERICAL robustness audit (the user's original point): TMB/RTMB
   ship dbinom_robust(x, size, logit_p) and dnbinom_robust, which
   evaluate the log-density directly from the linear-predictor scale
   and stay finite and differentiable at extreme eta, where
   dbinom(y, n, plogis(eta)) saturates to -Inf with garbage
   gradients. glmmTMB uses the robust forms internally, so our
   at-the-optimum agreement tests would NOT have caught off-mode
   fragility. Audit every family lpdf whose parameterization
   round-trips through an inverse link (bernoulli/binomial via
   plogis, negbinomial small-mu, zi gates at extreme zi, ordinal
   threshold differences) and switch to the robust RTMB forms where
   exposed; probe with extreme-eta gradient checks (the regime:
   separation, GK quadrature nodes, frm_sample tails, mixture
   gating). Check what RTMB exports before assuming.
2. t-distributed random effects, brms spelling gr(g, dist =
   "student") (brms#1876 - brms has it, so this is grammar-matching
   not invention). TMB latents need not be Gaussian; registry blocks
   would swap dnorm for dt. Open question: Laplace accuracy over the
   non-log-concave t latent density - feasibility probe vs adaptive
   quadrature and frm_sample BEFORE shipping.
3. huber() family: Huber's least-favorable distribution is a real
   density (gaussian center, Laplace tails), so ML with it is
   legitimate, same working-likelihood caveats as asym_laplace
   (document, point at frm_bootstrap). Validate point estimates vs
   MASS::rlm(psi = psi.huber). Afternoon-scale.
4. Cluster-robust (sandwich) vcov: per-OBSERVATION scores fight the
   marginalized objective (why sandwich::estfun was skipped), but
   per-CLUSTER scores are computable - the marginal likelihood
   factors over independent groups. CR0/CR2 with clubSandwich as
   reference. Highest practical demand of the four; changes
   inference for every clustered model.
robustlmm-style bounded-influence estimating equations (DAStau) are
NOT a likelihood and stay out of scope; item 2 is our answer to the
same concern.
## `nlf()`: recorded, not implemented (v0.34.x, wt-polish)

`lf()` shipped with the brms-portability batch: it is sugar over the
dpar formulas `bf()` already takes, so `+.frmtmb_formula` merges its
formulas into `pforms`. `nlf()` is not the same shape. In brms it
declares a NONLINEAR formula for one parameter,
`bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1)`, which needs a per-dpar
`nl` flag. frmtmb's `nl` is one flag on the whole `bf()`, and the
nonlinear branch of `parse_one_response()` reads the main formula's
right-hand side as the body with every `pforms` entry as a linear
nonlinear-parameter formula. Supporting `nlf()` means a nonlinearity
flag per dpar and a nested body in the objective, which is real work in
`parse.R` and `objective.R`, not sugar. The existing spelling stays
`bf(y ~ exp(b * x), b ~ 1, nl = TRUE)`.
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
