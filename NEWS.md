# frmtmb 0.27.0

Fixes from a full review of the v0.22-v0.25 waves.

## Corrected behavior

* `cens()` combined with `trunc()` produced a silently wrong
  likelihood: the censoring contribution was not restricted to the
  truncation window (right-censoring now contributes
  (F(ub) - F(y))/Z), which inflated the dispersion by ~14% on the
  reference problem. The corrected objective matches a hand-rolled
  windowed likelihood to 1e-8, collapses exactly to the old form
  without truncation, and the cens+trunc OSA residuals recalibrate
  to the analytic PIT (5e-14).
* `residuals(type = "deviance")` under `se()` now weights each
  row's unit deviance by its own variance (the glm prior-weight
  form) instead of treating a common dispersion as shared.
* `predict(se.fit = TRUE)` on quadrature fits no longer dies with a
  non-conformable error: it reports mode-conditional standard errors
  with a warning. Singular joint precisions in the predict path now
  degrade like `vcov()` does, and the shared solver uses
  `Matrix::solve` (strictly more robust for ill-conditioned but
  invertible GLMM precisions).
* `frm_sample()` chain inits are clamped inside any bounds passed to
  Stan, with a warning naming the parameters when the ML mode itself
  violates a bound (previously an incomprehensible rstan error).
* `vint()`/`vreal()` variables are required in `newdata`; a custom
  family's prediction previously returned a length-0 vector with no
  message when they were missing.
* Character censoring codes like "0"/"1"/"-1" decode (the error
  message had advertised them).

## Compatibility registry and fuzzer

* The registry's precedence is redesigned (lexicographic comparison
  of sorted side specificities, with a validator that forbids
  file-order ties except documented overrides). 488 of 3750
  resolutions were corrected against probed reality, including 429
  pairs whose covstruct conditions a family-level rule had erased,
  and the newly recorded caveat that spatial structures over a plain
  factor silently use level order as coordinates. `frm_compat()`
  accepts vectors and errors on empty input.
* The fuzz tier's invariants were strengthened (a previously
  unfailable interval check replaced by the Wald identity plus
  parameter-coverage tallies; per-row simulate agreement; grammar
  divergences asserted; convergence demotion keyed to actual
  convergence verdicts after fixing a truncation that discarded
  them). The strengthened tier reproduces only known findings.

# frmtmb 0.26.0

Pooled model comparison across imputations and the diagnostics/UX
backlog.

## New

* `anova()` on `frm_multiple()` fits pools nested-model tests across
  imputations: D1 (multivariate Wald), D2 (chi-square combining),
  and D3 (Meng-Rubin likelihood pooling, with the plug-in leg
  evaluated by re-taping each imputation's objective at the pooled
  parameters - no refits). D1/D2 match `mice` to 1e-7 on an
  exactly-shared reference; D3 is validated against the Meng-Rubin
  formula directly, since `mice::D3`'s `fix.coef` variant is not
  the plug-in statistic. Includes an ARIV clamp and a Reiter-df
  fallback for a small-m case where mice returns NaN.
* `cbind(successes, failures)` binomial responses are accepted (the
  glm/lme4/glmmTMB spelling), rewritten internally to
  `successes | trials(successes + failures)`; bit-identical to the
  `trials()` form.
* `simulate()` returns ordinal draws as ordered factors with the
  original levels and multinomial draws as count matrices (both
  families previously had no simulator), and respects `na.exclude`
  padding.
* `frmtmb_control(check_nlev_1 =, check_olre =)`: lme4-style
  warning/ignore/stop vocabulary for one-level grouping factors and
  gaussian observation-level random effects (the `se()`-based
  meta-analysis idiom is recognized and not flagged).
* `diagnose()` adds a complete-separation heuristic, predictor-scale
  warnings pointing at `autoscale`, and an isSingular-style verdict
  independent of the Hessian; it also no longer errors on fits
  without random effects.

## Corrected behavior

* Models with zero free outer parameters fit degenerately instead of
  dying inside nlminb (fixing latent empty-sdreport and
  empty-gradient bugs found along the way).
* A nonlinear-parameter name colliding with a data column errors
  instead of silently shadowing the column; nonlinear fits that fail
  from default zero starts name `start =` in the error.
* REML `anova()` compares fits whose fixed-effect designs span the
  same column space (term reordering included) instead of refusing
  every REML pair; genuinely different designs and REML/ML mixes
  still refuse with the reason.
* Offsets in distributional-parameter formulas were verified correct
  (to 1e-13) and are now regression-tested against the silent-drop
  failure mode reported upstream (glmmTMB#625).

# frmtmb 0.25.0

Simulation-workflow ergonomics, the last deferred method-surface
items, and CRAN readiness.

## New

* `frm_simulate()` accepts natural-scale parameter names - the same
  vocabulary as `hypothesis()`/`variables()`: coefficients by name,
  `sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`, response-scale
  `sigma` - inverted to the internal parameterization through the
  covariance registry (us/diag/homdiag/smooth/gr structures;
  others refuse by block name). The internal spelling still works
  and is byte-identical.
* Prior-predictive simulation, the `sample_prior = "only"` analog:
  `frm_simulate(..., priors = set_prior(...))` draws a fresh
  parameter vector per replicate on each prior class's documented
  scale (natural SDs for class `sd`, truncation by rejection) and
  attaches the drawn parameters for prior-predictive checks. Every
  coefficient and SD must be pinned by a prior or `newparams`;
  omissions error instead of silently becoming zero effects.
* `simulate(censored = TRUE)` applies the fitted type-I censoring
  mechanism to the draws; the default remains the latent uncensored
  response, which is also brms's `posterior_predict` convention
  (verified against brms source) and is now documented.
* `residuals(type = "deviance")` across the GLM family set
  (gaussian, poisson, binomial/bernoulli, Gamma, exponential,
  inverse.gaussian, nbinom1/2, geometric, beta, tweedie), exact
  against `stats::glm` where glm offers them; families without a
  standard unit deviance refuse by name. `deviance()` itself stays
  `-2 logLik` (lme4 convention).
* `predict(type = "response", se.fit = TRUE)` now works for
  non-identity-mean families (zi/hurdle, lognormal, trials-binomial,
  truncated responses) through a joint delta method across all dpar
  linear predictors, including cross-dpar covariance (agrees with
  glmmTMB's response-scale standard errors to 1.2e-5 on a
  zero-inflated poisson mixed model).

## Corrected behavior

* `conditional_effects(method = "predict")` evaluates addition terms
  on the effect grid, so its bands respect `trials()` and `trunc()`
  (binomial bands are counts in [0, n], not Bernoulli 0/1); aterm
  variables must be pinned in `conditions`, and the error names the
  variable.

## CRAN and infrastructure

* Heavy reference-validation test files are `skip_on_cran()`-gated:
  the CRAN-condition suite drops from ~128s to ~72s on the
  development machine while CI (NOT_CRAN=true) keeps full coverage.
* nlme added to Suggests (a test uses `nlme::Soybean` as reference
  data; CI's `--as-cran` unstated-dependency check halts without the
  declaration).
* Benchmark verdict recorded in dev/benchmarks.md: optimParallel is
  not adopted - with exact AD gradients its concurrency caps at two
  evaluations (measured 1.03-1.22x end to end, slower with cold
  clusters), RTMB tapes cannot ship to PSOCK workers, and 100% of
  InstEval's optimization time is inside the taped objective and the
  inner sparse Cholesky. `frmtmb_control(profile = TRUE)` remains
  the measured lever for many-coefficient models (1.6x there).

# frmtmb 0.24.0

The quadrature and OSA defect clusters surfaced by the fuzzer and
compatibility registry.

## Corrected behavior

* `quadrature = TRUE` is rebuilt around one root cause: TMB's
  Gauss-Kronrod marginalization calibrates each integrand's rescaling
  once, at the parameter values in hand when the tape is built, and
  frmtmb taped at the cold start. A plain Laplace fit now runs first
  and the quadrature tape is built at its optimum. This fixes:
  conditional modes returned as NA (`ranef()`/`fitted()`/`predict()`
  were silently NA for all groups but the first, whose slot held a
  wrong value) - modes now come from the Laplace inner solve and
  match `glmer(nAGQ = 25)`'s `ranef()` to 3e-5; and bare "NA/NaN
  gradient evaluation" crashes for poisson/Gamma/Beta with nested or
  even single scalar blocks - all now fit with gradients < 3e-4.
* `quadrature` combined with `trunc()` produced logLik = +Inf as a
  successful fit; the combination is refused (the CDF difference
  underflows at quadrature nodes; a stable fix needs log-CDF forms).
  `mixture()` with `REML` or `profile = TRUE` is refused: both
  Laplace-expand the mu coefficients around a single mode, and a
  mixture likelihood is permutation-multimodal in exactly those
  coefficients. `mixture()` with quadrature remains supported.
* `residuals(type = "osa")` on censored fits was broken (raw LAPACK
  singularity or NaN for every censored row). Uncensored rows now get
  calibrated one-step residuals by conditioning on the censoring
  window and renormalizing the one-step CDF to it (matches the
  analytic conditional PIT to 7e-9); censored rows are NA (an event
  has no one-step CDF), and row-varying censoring points refuse.
  Note `simulate()` draws the latent uncensored response under
  `cens()`, so DHARMa is not a substitute there (documented).
* `residuals(type = "osa")` on ordinal fits crashed inside the OSA
  machinery; the ordinal log-densities now handle `oneStepPredict`'s
  taped-observation objects (exact Lagrange-basis category
  indicator), and all four ordinal families produce calibrated
  residuals (cumulative matches the analytic randomized-quantile
  residual to 4e-14).
* A singular joint precision under REML or `profile = TRUE` no longer
  throws a raw LAPACK error from `vcov()`/`summary()`/`confint()`/
  `hypothesis()`: it degrades to NaN standard errors with one warning
  pointing at `diagnose()`, matching the ML branch.
* The brms-migration vignette documents that `binomial()` without
  `trials()` means Bernoulli here (glm convention) where brms
  requires `trials()` or `bernoulli()`.

The fuzz tier's open findings drop from 28 to 16 on the identical
plan; most of the remainder are now informative refusals on thin-data
edge cases rather than defects.

# frmtmb 0.23.0

Defect wave driven by an open-issue sweep of brms/lme4/glmmTMB, plus
a feature-compatibility registry.

## Corrected behavior

* Truncation now reaches the whole post-fit surface, not just the
  likelihood: `fitted()`, `predict(type = "response")`, and
  `residuals()` report the truncated mean E[Y | lb <= Y <= ub]
  (closed forms for all six CDF families, validated to 1e-15);
  `simulate()`, `posterior_predict()`, and `frm_simulate()`
  rejection-sample within bounds; newdata re-evaluates variable
  bounds. `residuals(type = "osa")` on truncated (and untruncated
  gaussian) responses was also miscalibrated and now integrates over
  the truncated support (KS uniformity restored from p ~ 1e-15 to
  0.76-0.97). Dpar-scale predictions stay untruncated by design.
  DHARMa/pp_check on truncated models are meaningful again.
* `(f || g)` with a factor now yields independent per-level effects
  (the diag structure the syntax promises) instead of a silently
  fully correlated block; identical to an explicit `diag(f | g)`.
  Numeric double-bars are unchanged. Existing factor-double-bar fits
  change, because the old ones were wrong (lme4 has the same open
  bug, lme4#818).
* `ar1()`/`hetar1()` warn when the ordering factor's integer levels
  have gaps: levels correlate by position (the glmmTMB reading,
  unchanged), so a gap counts as one step; the warning points at
  `ou()` over `num_factor()` for irregular spacing (glmmTMB#1278).
* Predictions at non-estimable points of a rank-deficient design
  return NA with one warning naming the aliased columns, instead of
  silently returning the partial sum (predict.lm semantics;
  lme4#303). Collinear restatements of kept columns stay exact.
* Grouping factors written as calls work: `(1 | factor(x))`,
  `(1 | interaction(a, b))` (lme4#464).
* A random-effect term crossed with `*` or `:` errors instead of
  being silently refit additively (lme4#196); `mo()`/`mi()`
  interaction multipliers must be numeric (brms#1828); `anova()`
  requires equal `nobs` (lme4#622).
* `frm_sample()` chain initialization: chain 1 anchors at the ML
  mode, further chains are jittered on the unconstrained scale
  (`init_jitter`), restoring the overdispersion Rhat needs; a
  boundary-mode warning fires for singular fits; mixture posteriors
  are documented as needing `init = "random"`.

## New

* Feature compatibility registry: `frm_compat()` answers what plays
  with what (works / conditional / refused / broken / untested)
  across 3750 feature pairs, and the new "Feature compatibility"
  article is generated from the registry at build time so it cannot
  drift. Known-broken pairs are listed there; the registry probing
  itself surfaced trunc x quadrature and cens x OSA as broken
  (queued for the next fix wave).
* An audit of 559 currently open brms/lme4/glmmTMB issues is
  recorded in dev/test-backlog.md, including the pathologies frmtmb
  is structurally immune to.
* A pairwise grammar fuzzer (env-gated: `FRMTMB_FUZZ=true`): a
  310-spec covering array over the grammar with metamorphic
  invariants (predict/fitted identity, permutation invariance,
  simulator support membership, vcov sanity) and a brms
  `make_standata` structural oracle. Its first run surfaced the
  quadrature defect cluster now queued for fixing (conditional
  modes left NA, crashes with nested groups and Beta, +Inf logLik
  with trunc()); the tier reports those as failures until they are
  fixed.

# frmtmb 0.22.0

Cross-validation against brms itself, closer brms compatibility, and
fit ergonomics.

## brms agreement suite

* New test tier validating the elaborate grammar structurally against
  `brms::make_standata()`/`brmsterms()` without any Stan compilation
  (designs, RE structures, `|ID|` merging, ordinal thresholds, `cs()`,
  `mo()` codes, gp bases, smooths, addition terms, mvbf, nl, mixture
  naming; most exact to 1e-10 or better), plus an opt-in numeric tier
  (`FRMTMB_BRMS_FIT_TESTS=true`) showing fixed-effects estimates match
  the mode of brms's own generated Stan program to 1e-4 and that our
  `mo()` parameterization is brms's likelihood exactly (via
  `rstan::log_prob`). brms is in Suggests only.
* Found by the suite and fixed: `cs()` predictor variables never
  reached the combined model frame, so `bf(y ~ x + cs(z))` errored
  unless `z` already appeared elsewhere.

## brms compatibility

* Hilbert-space `gp(x, k =)` now uses brms's input convention exactly:
  coordinates are rescaled by the largest pairwise distance (one
  shared factor across dimensions) and centered, so the boundary is
  `L = c` and the same `gp()` call is the same approximation in both
  packages (basis agreement with brms to 1e-16, tied coordinates and
  vector-valued `c =` included). This roughly doubles the effective
  boundary at a given `c`: accuracy against the exact gp improves
  about 25x at the same k (k = 40 now within 2e-4 logLik), and the
  default `c = 1.25` is the right choice in 2-D as well. Reported
  `range(gp)` values stay in data units (exact log-shift
  back-transform); fitted lengthscales differ from pre-0.22 fits by
  the data-dependent scale factor.
* `cens()` accepts brms's character and factor codes
  (`"left"`/`"none"`/`"right"`/`"interval"`, prefix-matched,
  case-insensitively); unknown labels and out-of-range numeric codes
  error informatively instead of coercing to NA.
* `gp()`'s `k`/`c`/`iso`, `rr()`'s `d`, and `se()`'s `sigma` arguments
  now evaluate in the formula environment (brms behavior), so
  variables and expressions work; invalid values error with the
  offending expression named.

## Fit ergonomics

* `frmtmb_control(verbose =)` (and the `frm(verbose =)` shortcut):
  level 1 prints timed stage lines (parse, frame, tape, optimize,
  restarts, sdreport) through `message()` so a slow fit shows where
  it is slow; level 2 adds the optimizer's own iteration trace.
  Bootstrap, influence, and allfit refit loops stay quiet.
* The large-gradient warning now points to `diagnose()` and the new
  "Convergence problems" remedy ladder in the diagnostics vignette
  (scaling, restarts, optimizer comparison, profiles, boundary fits,
  and judging marginal gradients).
* `frm_allfit()`'s nloptr optimizer produced an empty row: nloptr
  rejects callbacks that declare `...` (RTMB's `obj$fn`/`obj$gr` do),
  and the error was swallowed; separately, missing `ftol_rel` made
  NLopt report failure at the converged optimum. Both fixed; the
  agreement test now requires every optimizer to converge and match.
* The brms-migration vignette documents how to recover the model
  function (the `stancode()` analog): the objective closure via
  `build_objective(fit$frame)`, its joint-vs-marginal relationship to
  `fit$obj$fn`, and how to reproduce `fit$obj`.

# frmtmb 0.21.0

Method-surface audit against lme4/glmmTMB/brms, plus fixes from a
full code review.

## Corrected behavior

* `predict(type = "response")` now returns the expected response for
  every family (equal to `fitted()` in-sample), not the primary
  dpar's natural scale; zi/hurdle, lognormal, and trials-binomial
  fits were affected. The glmmTMB aliases `"conditional"`, `"zprob"`,
  `"zlink"`, and `"disp"` are accepted (validated against glmmTMB on
  a zero-inflated Poisson to ~8e-7). Note `"disp"` returns `sigma` on
  its natural scale, where glmmTMB returns the variance-scale
  dispersion for gaussian. `se.fit` for the mean of a
  non-identity-mean family is not yet available and errors with
  guidance.
* `rescor = TRUE` now refuses `cens()`, `trunc()`, and `se()`
  addition terms, which the joint-gaussian likelihood silently
  ignored (wrong likelihood with no warning).
* Bounds (`lower`/`upper`, `set_prior` lb/ub) under
  `frmtmb_control(profile = TRUE)` were positionally misaligned:
  a bound on a fixed coefficient could silently pin a covariance
  parameter instead. Bounding a profiled coefficient now errors;
  covariance-parameter bounds land on the right slot.
* `confint()`, `vcov(full = TRUE)`, and `diagnose()` labels were
  broken or misaligned on `mi()` fits: the latent imputation
  component was counted as an outer parameter.
* `frm_sample(priors = )` failed on fixed-effects-only fits of
  single-dpar families (a `$` partial-match bug), and
  `frm_sample(laplace = TRUE)` was broken twice: the default init
  had the wrong length, and draw columns were mislabeled (a column
  named `b[1]` actually held a covariance parameter).
* `influence()` tables and `cooks.distance()` now align with
  `vcov()` (fits with a constant dpar errored; column names are now
  the coefficient labels).
* `simulate()` follows the stats seed contract: a `"seed"` attribute
  and a restored RNG state instead of clobbering the global stream.
* A covariate literally named `sigma` is no longer shadowed by the
  residual sigma in `hypothesis()`/`variables()`.
* In `predict(newdata = )`, an NA in a variable used only in a
  random-effect design now propagates to the prediction instead of
  being silently zeroed (only genuinely new levels predict at the
  population value).

## New methods

* `drop1()` (AIC/LRT, marginality-aware; matches
  `lme4::drop1.merMod` to 1e-4), `cooks.distance()` directly on a
  fit, `dfbeta()`/`dfbetas()` on `influence()` results (stats sign
  convention; lme4 returns the negation), `na.action()`, and an
  lme4-style "Groups:" line in `summary()`.
* `confint()` accepts lme4's `"Wald"` spelling and gains
  `method = "boot"` (percentile intervals via [frm_bootstrap()]).
* `vcov(full = TRUE)` rows carry per-parameter names (glmmTMB
  convention).
* `predict()` accepts the lme4/glmmTMB `allow.new.levels` dot
  spelling.
* On draws objects: `posterior_linpred()` and `ranef()`
  (brms-shaped arrays).
* `emmeans::recover_data` had been silently missing from NAMESPACE
  since v0.13 (orphaned export tag); emmeans support was broken.

## Internals

* One authoritative outer-parameter map shared by `confint()`,
  bounds resolution, `vcov(full = TRUE)`, and sampling; `graphics`
  and `grDevices` added to Imports; duplicated formula/family
  coercion and gp position-key helpers unified; REML `logLik` df
  counts only outer parameters (documented; lme4 counts the
  integrated fixed effects too, so REML AICs are not comparable
  across packages).

# frmtmb 0.20.0

* `mixture_mvn(K, D)`: multivariate gaussian mixture components for
  model-based clustering of an n x D matrix response (the mclust /
  clustTMB use case). Every class mean is a full linear predictor, so
  cluster means may depend on covariates and random effects; mixing
  weights are multinomial-logit dpars with their own formulas.
  Validated against a hand-rolled ML fit to 1.3e-10 and against
  faithful-data cluster recovery. Limitations (documented in
  `?mixture_mvn`): class covariances are unstructured (`us`) and
  covariate-free - no mclust-style constrained covariance taxonomy
  (EII..VEV); `cens()`/`trunc()`, `mvbf()`, and `simulate()` are not
  supported.
* Multi-dimensional Gaussian processes: `gp(x1, x2, ...)` takes up to
  three variables, with a separate lengthscale per dimension by
  default (the brms convention) or one shared lengthscale via
  `iso = TRUE`, in both the exact and the Hilbert-space (`k =`) form.
  The tensor HSGP basis is capped at `k^D <= 1000` columns.
  `confint_varcorr()` reports one range row per dimension. Validated
  against a direct 2-D ML fit to 3.8e-11.
* Kriging prediction for the exact `gp()`: `predict(newdata = )` at
  positions not seen at fit time now returns the GP conditional mean,
  and `se.fit = TRUE` adds the conditional variance, instead of
  erroring. Kriging mean validated to 1.7e-15 and the standard error
  to the exact decomposition identity. The Hilbert-space form already
  interpolated through its basis.
* `frmtmb_control(sparse_x = TRUE)`: sparse fixed-effect design
  matrices (`Matrix::sparse.model.matrix`), the analog of
  `glmmTMB(sparseX = )`, for models where a many-level fixed factor
  dominates memory. Estimates are identical to the dense path
  (gradient agreement 3e-14; 13.8% lower tape memory in the
  validation model); `model.matrix()` on the fit returns a sparse
  matrix. Prediction falls back to the dense builder for newdata
  containing NA factor rows, which `sparse.model.matrix` would
  silently zero.
* `frmtmb_control(autoscale = TRUE)`: internal standardization of
  badly scaled continuous predictors (the lme4 >= 1.1.37 feature): a
  scaled pre-fit is mapped back exactly and warm-starts the reported
  fit, parameter magnitudes feed `nlminb`'s scale hook, and
  convergence and standard-error Hessian checks run in scaled units.
  Turns 12 scaling warnings into 0 on the validation problem with a
  0.0 log-likelihood difference against a manually standardized
  reference. Reported results are always on the original scale.
* `frm_multiple()` pooling extensions: `$pooled_varcorr` pools
  random-effect SDs and correlations across imputations on
  transformed scales (log for SDs and GP ranges, Fisher z for
  correlations) with Barnard-Rubin degrees of freedom, and
  `hypothesis()` on a `frmtmb_multiple` object pools delta-method
  hypothesis estimates by Rubin's rules with t-based intervals.
  Fixed-effect pooling matches `mice::pool()` to 2.1e-6.
  Likelihood-ratio pooling across imputations (mice's D3) is not
  implemented; `anova()` on pooled fits remains per-fit.

# frmtmb 0.19.0

* Latent classes now combine with continuous random effects
  (growth-mixture models): `mixture(..., groups = ~g)` accepts
  random effects, smooths, and gp() terms in the component formulas.
  The class sum is computed conditional on the latent effects, so
  one Laplace approximation integrates them - the sum-and-integral
  swap makes this exact as an integrand identity, and even crossed
  random effects are structurally fine. Random effects written in a
  component formula are class-specific by construction. Accuracy:
  the Laplace approximation of the class-mixture integrand carries a
  small documented bias (about 0.1 log-likelihood units in the
  validation problem, with parameter estimates matching the exact
  closed-form marginal to 0.01); `quadrature = TRUE` is numerically
  exact when the per-group integrand is univariate (validated to
  5e-3) and approximate when class-specific intercepts couple.
  `simulate()` now draws one class per group and then simulates each
  observation from its group's component; `mixture_probs()` gives
  empirical-Bayes classification conditional on the modes.
* `variables()` (brms spelling): lists every parameter name usable in
  `hypothesis()` expressions (coefficients, `sd_`/`cor_` summaries,
  `sigma`); on `frm_sample()` output it lists the draw columns.
* `get_prior()` (brms spelling): enumerates every slot `set_prior()`
  can target (class/coef/dpar/group rows), from a formula plus data
  or from a fit; the default in every slot is flat.
* README rewritten for the current scope, with related-work pointers
  (glmmTMB, BayesRTMB); getting-started vignette extended to the full
  grammar; SPEC.md status and deviations brought current.

# frmtmb 0.18.0

The last three roadmap features.

* `gp(x)` Gaussian-process terms: exact (a dense squared-exponential
  block over the unique positions, with a standard 1e-6 nugget) or
  the Hilbert-space approximation via `gp(x, k = 30)` (sine basis
  with spectral-density prior SDs), in any linear predictor. Exact
  fits match direct GP marginal ML; the approximation converges to
  the exact answer and predicts at arbitrary positions (the exact
  form predicts at observed positions only). `confint_varcorr()`
  reports `sd(gp)` and `range(gp)`.
* `mo()` and `mi()` interactions: `mo(x):z`, `mo(x)*z`, `mi(x):z`,
  `mi(x)*z` with numeric multipliers; `mo()` interactions share
  their variable's simplex (brms convention). Both validated exact
  against direct ML (the `mi()` interaction stays linear in the
  latent value, so the closed-form marginal still applies).
* Group-level latent-class mixtures: `mixture(..., groups = ~g)`
  sums the class assignment per group (growth-mixture /
  latent-class regression; the tractable nested case of brms#1905).
  Exact against direct ML; `mixture_probs()` returns posterior class
  probabilities per group, or per observation for ordinary mixtures.
  Restrictions: no random effects or smooths alongside, mixing
  predictors evaluated at each group's first row, no `simulate()`
  yet. Crossed-design group mixtures remain out of the Laplace
  class.

# frmtmb 0.17.0

Clearing the roadmap: mixtures, measurement error, and the remaining
grammar gaps.

* `mixture(fam1, fam2, ...)` finite-mixture families: each component
  keeps its own suffixed dpars (`mu1`, `sigma1`, ...), mixing
  proportions are multinomial-logit dpars with their own linear
  predictors (mixture-of-experts comes free), and the likelihood is a
  branch-free logsumexp - so random effects and smooths work in any
  component formula. Exact against direct ML for gaussian and poisson
  mixtures. Component means initialize on spread response quantiles;
  the usual finite-mixture multimodality caveats apply
  (`frm_allfit()`, or order intercepts via bounds).
* `mi(sdx)` measurement error (brms `me()`): known per-observation
  measurement SDs make every true value latent, with the observed
  values entering through a measurement model. Exact against the
  closed-form bivariate-normal marginal; attenuation bias is
  corrected and NAs (missing + mismeasured) combine.
* `cs()` category-specific effects for `sratio`, `cratio`, and `acat`
  (refused for `cumulative`, as in brms). Exact against direct ML.
* `equalto(x + 0 | g, V)` covariance structure: fully fixed known V,
  zero parameters (meta-analytic sampling covariances).
* rr() fits now support `se.fit`: the delta method routes the factor
  columns through the loadings Jacobian and adds loading-parameter
  columns; verified parameterization-invariant against the
  full-rank-vs-`us()` equivalence.
* Registered insight methods (`find_formula`, `find_random`,
  `get_parameters`, `get_varcov`, `find_statistic`, link accessors):
  `insight::is_mixed_model()` and the easystats accessor layer now
  see the random-effect structure.

# frmtmb 0.16.0

The two deferred architecture features.

* `rr(x | g, d = k)` reduced-rank / factor-analytic covariance
  (glmmTMB's GLVM structure): per-level coefficients are loadings
  times iid standard-normal factors, so `b` holds the factors while
  the design matrices span the coefficient space. Matches glmmTMB
  exactly and reduces to `us()` at full rank. `se.fit` on rr models
  is deferred (the loadings enter the coefficients nonlinearly);
  `predict`, `simulate`, `ranef`, and `VarCorr` (the rank-k
  covariance) all work.
* `mi()` one-step imputation for continuous predictors, reversing the
  original exclusion: `bf(y ~ mi(x) + z) + bf(x | mi() ~ z)` turns
  missing `x` values into latent parameters that Laplace integrates
  (exactly, in the linear-gaussian case - validated against the
  closed-form marginal likelihood). Imputation uncertainty
  propagates into the coefficient standard errors automatically.
  Rows with NAs in non-`mi()` variables still drop; imputation
  models must be gaussian or student; discrete predictors remain
  impossible (as in Stan). `me()` measurement error is the natural
  follow-on.

# frmtmb 0.15.0

Tier-2 sweep of dev/feature-gaps.md plus a method surface for sampled
fits.

* `mo()` monotonic effects (brms syntax): scale coefficient times an
  estimated simplex, exact against direct ML; works in prediction,
  `conditional_effects()`, and with ordered factors. No frequentist
  package offers this. Standalone additive terms only for now.
* Sequential ordinal families `sratio`, `cratio`, `acat`, validated
  against direct ML and collapsing to logistic regression at K = 2.
* Covariance structures: `hetar1`, `homcs`, `homtoep` (glmmTMB
  parameterizations, exact where glmmTMB converges) and spatial
  `exp()`, `gau()`, `mat()` over `num_factor(x, y)` coordinates
  (Matern uses bounded internal transforms for stability where both
  we and glmmTMB otherwise diverge).
* Fixed a latent parser bug: `exp(x)` (or any covariance-structure
  name) used as a plain function in a fixed formula was silently
  stripped by the formula splitter - since v0.1. Such calls are now
  protected automatically.
* `influence(fit, groups = )` and `cooks.distance()`: case-deletion
  diagnostics over warm-started refits.
* `frm_multiple()`: fits across multiply-imputed datasets (list or
  `mice::mids`) pooled by Rubin's rules with Barnard-Rubin df.
* `conditional_effects()` gains `method = "predict"` (prediction
  intervals) and data-frame `conditions` (one condition set per row,
  brms style).
* `vint()`/`vreal()` addition terms pass arbitrary data vectors to
  custom families (brms custom-family convention). The test suite
  includes a full Wiener drift-diffusion model written as a
  `custom_family()` in plain R - the workflow brms needs raw Stan
  code for.
* Method surface for `frm_sample()` draws: `summary()` (with
  Rhat/ESS), `fixef()`, `VarCorr()` (natural scale), `hypothesis()`
  (exact posterior version), `posterior_epred()` /
  `posterior_predict()` (each draw runs the full prediction
  machinery), `pp_check()`, and `posterior::as_draws()`.
* Examples added across the reference (run by R CMD check).

# frmtmb 0.14.0

Gap-closing milestone from a sweep of the lme4 / glmmTMB / brms
vignettes (dev/feature-gaps.md).

* `se()` addition term: known per-observation sampling SDs
  (meta-analysis), gaussian and student. `se(x)` fixes the residual SD
  (sigma is mapped out); `se(x, sigma = TRUE)` adds an estimated sigma
  in quadrature. With `(1 | obs)` this is random-effects
  meta-analysis, and with `gr(g, cov = A)` the phylogenetic version.
* Binomial (and beta-binomial) responses may be proportions of
  `trials()`, the glm/glmer idiom, converted to counts internally.
* New families: `bernoulli`, `geometric`, `exponential`, `weibull`
  (mean-parameterized, matches survreg's likelihood),
  `shifted_lognormal`, `hurdle_gamma`, `hurdle_lognormal`,
  `zero_inflated_binomial`, `zero_inflated_beta` (validated against
  glmmTMB), and `asym_laplace` for quantile regression (fixed
  `quantile` dpar reproduces `quantreg::rq` estimates).
* `ranef(condVar = TRUE)`: conditional SDs of the modes (TMB/glmmTMB
  convention: fixed-effect uncertainty propagates), plus tidy
  `as.data.frame()` methods for `ranef` and `VarCorr` output.
* `frm_allfit()`: refit under every available optimizer (nlminb,
  optim, bobyqa, NLopt) and compare, the lme4 `allFit()` analog.
* `frm_simulate()`: de novo simulation from a `bf()` formula and
  supplied parameters with no fitted model (power analysis; the
  glmmTMB `simulate_new()` analog).
* `frmtmb_control(profile = TRUE)` (experimental): profile the
  primary-coefficient vector into the inner Laplace problem, the
  `glmmTMBControl(profile = TRUE)` / `glmer(nAGQ = 0)` speed
  approximation for many-coefficient models; covariance comes from
  the joint precision.
* `anova()` no longer fails on models sharing a primary formula
  (distributional vs plain fits); row labels now include dpar
  formulas.

# frmtmb 0.13.0

Method-surface ("sugar") milestone: the conventional S3 methods that
downstream packages and user muscle memory dispatch on.

* New accessors: `sigma()`, `terms()`, `weights()`, `model.matrix()`,
  `deviance()`, `extractAIC()` (enables `step()` / `drop1()`),
  `ngrps()`, `prior_summary()`, and a `profile()` method wrapping
  `TMB::tmbprofile()`.
* **Breaking**: `coef()` now follows the lme4 / glmmTMB / brms
  convention (per-group coefficients: fixed effects broadcast over
  levels plus the conditional modes). Use `fixef()` for the fixed
  effects alone. Fits without random effects still return the
  coefficient vector.
* `conditional_effects()` with a base-graphics `plot()` method: grid
  predictions per predictor (or `"x:z"` pair) with Wald bands, smooths
  included, matrix covariates held at column means.
* `plot(fit)`: Pearson-residual diagnostics (residuals vs fitted,
  normal QQ).
* `pp_check()` (registered on the bayesplot generic): predictive checks
  from `simulate()` draws through any `ppc_*` function.
* `hypothesis()`: tests of arbitrary expressions of the parameters,
  brms spelling, with three methods: delta-method Wald (default),
  `"profile"` (profile-likelihood intervals via `TMB::tmbroot()`
  lincombs; linear hypotheses, ML fits), and `"boot"` (parametric
  bootstrap percentile intervals; any expression). The expression
  environment includes natural-scale random-effect names
  (`sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`, `sigma`), so ICC
  inference is `hypothesis(fit, "sd_g__Intercept^2 /
  (sd_g__Intercept^2 + sigma^2)", method = "boot")`. The result is a
  `frmtmb_hypothesis` data frame carrying the bootstrap draws or
  profile curves in attributes, with a `plot()` method (bootstrap
  histogram, profile curve, or the implied Wald normal, per
  hypothesis).
* `frm_bootstrap(fit, FUN, nsim)`: parametric bootstrap over
  `refit()`s (the `bootMer` analog), with `confint()` percentile
  intervals.
* `refit(fit, newresp)`: refit to a new response reusing the assembled
  design with a warm start - one re-tape plus optimization; the engine
  for parametric bootstrap over `simulate()` draws.
* The accessor set makes insight's default methods work
  (`find_response()`, `get_data()`, `get_sigma()`, `n_obs()`, ...);
  dedicated insight methods (random-effect formula parts,
  `get_parameters()`) are a roadmap item.

# frmtmb 0.12.0

Pre-release consolidation toward CRAN. Highlights by milestone:

## Model grammar (brms-compatible spelling)

* `bf()` with distributional-parameter formulas carrying the full
  predictor grammar (random effects and smooths in any dpar), constant
  dpars, and addition terms `weights()`, `trials()`, `cens()` (left /
  right / interval), `trunc()` (with the inclusive discrete correction).
* Nonlinear formulas (`nl = TRUE`) with a full linear predictor per
  parameter.
* Multivariate models: `mvbf()` / `bf() + bf()` / `mvbind()`, a family
  per response, gaussian residual correlation (`set_rescor()`), and
  `|ID|` random-effect correlation across formulas.
* mgcv smooths `s()` / `t2()` in any linear predictor, including
  matrix-covariate summation-convention terms: scalar-on-function,
  function-on-scalar, and function-on-function regression validated
  against mgcv.
* Covariance structures: `us`, `diag`, `homdiag`, `cs`, `ar1`, `toep`,
  `ou` (positions via `num_factor()`), and known-structure terms
  `gr(g, cov = A)` (dense, with correlated slopes via a Kronecker
  product) and `gr(g, prec = Q)` (sparse GMRF).

## Families

gaussian, poisson, binomial, Gamma, lognormal, student, negbinomial,
nbinom1, beta, tweedie, compois, beta_binomial, skew_normal,
inverse.gaussian, exgaussian, zero-inflated and hurdle counts,
cumulative ordinal, matrix-response multinomial, and `custom_family()`
as a plain R log-density with `check_custom_family()` AD verification.

## Estimation and inference

* ML and REML by Laplace approximation; `quadrature = TRUE` upgrades
  scalar random effects to adaptive Gauss-Kronrod marginalization
  (matches `glmer(nAGQ = 25)` and GLMMadaptive).
* MAP / regularized ML via brms-style `set_prior()`; hard bounds
  (`lower` / `upper`, and per-prior `lb` / `ub`).
* `sdreport` is deferred until standard errors are first needed
  (roughly a quarter off fit time; `se = TRUE` restores eager mode).
* Pluggable optimizers (`frmtmb_control(optimizer = )`), including
  arbitrary user functions.
* `confint()` (Wald / profile / likelihood-root), natural-scale
  `confint_varcorr()`, `anova()` likelihood-ratio tests, `diagnose()`.

## Diagnostics and ecosystem

* Residuals: response, pearson, and one-step-ahead (`type = "osa"`).
* DHARMa simulation-based residuals (`dharma_residuals()`).
* NUTS on the fitted objective: `frm_sample()` (ML-mode
  initialization, named draws, priors, bounds) and `check_laplace()`;
  bayesplot works on the draws.
* emmeans and marginaleffects support, both verified against glmmTMB.

## Verification

Roughly 500 tests compare fits against glmmTMB, lme4, mgcv, MASS,
survival, nnet, GLMMadaptive, and hand-written RTMB references, plus
an edge-case suite mined from the lme4 / glmmTMB / brms issue
trackers (dev/test-backlog.md).
