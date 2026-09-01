# Test backlog mined from lme4 / glmmTMB / brms

Edge cases harvested from the reference packages' test suites and
issue/PR history (2026-08-31). Status: DONE items have tests in
tests/testthat/ (mostly test-edgecases.R); the rest are open work.

## Done

- Frozen data-dependent bases (poly/ns/scale) via predvars stored from
  the fit-time model frame, applied to every linear predictor and RE
  component at prediction; single-row newdata. [glmmTMB#402, #512, #853;
  lme4 predict_basis; brms#494]
- Rank-deficient X: drop aliased columns per linear predictor with a
  message naming the component; frozen column set reused for newdata.
  [lme4#144; glmmTMB test-checkRank.R]
- scale()-style n x 1 matrix responses dropped to vectors. [glmmTMB#937]
- Row-permutation invariance for us/ar1; relevel invariance for us.
  [brms#1747 - silent bias with sorted rows; glmmTMB test-varstruc.R]
- Numeric/character/factor grouping equivalence. [lme4 test-factors.R]
- Duplicate multivariate responses rejected. [brms tests.standata.R]
- trials() validation: y > trials, non-integer y, constant literal.
  [brms data-response.R taxonomy]
- Gaussian rescaling equivariance of coefficients and logLik.
- NA rows dropped across response/predictors/grouping jointly.

## Addressed in v0.6 (moved up from Open)

- Non-default contrasts in prediction (global option and per-factor).
- re.form = NA without grouping columns in newdata.
- Formula-environment robustness: combined model frame stored on the
  fit; model.frame/predict/emmeans survive the calling env vanishing.
- na.exclude: fitted/residuals/predict padded via napredict; na.pass
  errors informatively; Inf responses rejected.
- Slash grouping (1|a/b) equivalence with (1|a) + (1|a:b).
- cens() x dpar-formula cross-product (vs hand-rolled reference).
- Frequency-weight equivalence for aggregated poisson data.
- dpar/nlpar names with dots or underscores rejected.
- predict() warns on unknown arguments.

## Open - high priority

1. Interval censoring (cens code 2 + y2), with NA y2 allowed on
   non-interval rows. [brms#1070]
2. Discrete-family truncation off-by-one: P(lb <= Y) needs F(lb-1), not
   F(lb), once count families get lcdf. [brms#1903]
3. |ID| new-level prediction must draw from the joint block, not the
   marginals; dpar-formula reordering must not change logLik.
   [brms#779, #674]
4. Level ordering of (1|a:b) interaction factors matching
   droplevels(a:b) label order (count verified; label order and star
   syntax (1|a*b) untested). [lme4#635/#636/#945]
5. Constant-weight (non-)invariance for gaussian documented and tested
   (frequency semantics now verified for poisson). [lme4 priorWeights.R]

## Open - medium

- One-level grouping factors / gaussian OLRE: three-way
  ignore/warn/stop control vocabulary. [lme4 lmerControl checks]
- Complete separation in binomial: diagnose() thresholds for big
  coefficients (|b|>10) and flat directions. [glmmTMB diagnose()]
- Singular-fit detection independent of the Hessian (isSingular
  equivalent); profile CIs on boundary parameters must warn not hang.
  [lme4 test-isSingular.R, #660]
- Predictor-scale warnings (|log10 sd| > 3). [glmmTMB diagnose()]
- Offset in dpar formulas; offset-only models (zero-column X in every
  predictor - partially covered by ordinal); offset argument form.
  [glmmTMB test-offset.R, #625, #286]
- predict type grid for zi families: response = (1-zprob)*conditional;
  zprob on non-zi model returns 0 not garbage; truncated conditional
  means. [glmmTMB#798, #873, #634]
- gam-style exclude= for zeroing individual smooths in prediction;
  re.form = NA keeps smooths (implemented - keep regression test).
  [glmmTMB test-smooths.R; mgcv semantics]
- simulate() returns original response type (ordered factor for
  ordinal, matrix for multinomial) and respects na.action. [glmmTMB
  test-simulate.R; lme4#737]
- nlpar name colliding with a data column: error or documented
  precedence (currently the nlpar value wins silently - decide).
  [brms#391]
- nl start values: fail informatively when optimizer diverges from
  0-starts, or require start for nl fits. [brms#734 doctrine]
- REML logLik invariance to fixed-effect term order; anova() refusing
  REML fits with different X column spaces (currently refuses all REML).
  [glmmTMB#776]
- confint/vcov excluding mapped parameters everywhere (constant dpars
  covered; extend to any future map use). [glmmTMB#1120]
- dpar/nlpar names with dots or underscores rejected (collision with
  coefficient naming). [brms tests.brmsformula.R]

## Open-issue sweep (2026-09-01)

Mined from the *currently open* trackers of brms (145), lme4 (191) and
glmmTMB (223); 34 shortlisted and run against frmtmb. Fixed items have
regression tests in tests/testthat/test-open-issues.R, except lme4#303
and lme4#464/#156, which live in
tests/testthat/test-aliased-grouping.R.

### Fixed

- Random-effect terms crossed with `*` or `:` (`y ~ x * (1 | g)`) now
  error instead of being silently refit as `+`. [lme4#196]
- `mo()`/`mi()` interaction multipliers: the numeric type gate never
  fired for character vectors (`is.numeric(as.numeric("a"))` is TRUE),
  so the column went all-NA and the fit died at "NA/NaN gradient
  evaluation". [brms#1828]
- `anova()` rejected fits with different `nobs`; previously it compared
  likelihoods across data sets and returned a negative Chisq.
  [lme4#622]
- `||` over a factor produced a fully correlated `us` block (|cor| up to
  0.99), the opposite of what the syntax promises. `parse_linpred()` now
  expands each `||` term itself and tags every piece `diag()`, so a
  factor's levels get independent variances; the fit is identical to an
  explicit `diag(f | g)`. Numeric double bars keep lme4's block split
  and their old estimates, because `diag` and `us` coincide at dimension
  one. [lme4#818]
- Prediction at an aliased cell of a rank-deficient design returned the
  random-effect contribution alone. The fit now freezes the null space
  of the fixed-effect design, and `predict(newdata =)` returns NA (with
  NA standard errors) plus one warning naming the dropped columns for
  the rows that load on it. Rows that merely restate a kept column
  (`x2 = 2 * x`) stay exact, and the in-sample paths are untouched.
  lme4 still returns the partial sum silently. [lme4#303]
- Grouping factors written as calls - `(1 | factor(x))`,
  `(1 | interaction(a, b))` - now fit. reformulas re-evaluates the bar
  RHS inside the model frame, where the call's own arguments are not
  columns, and died with an error raised several frames down
  ("unique() applies only to vectors"). Frame assembly now points the
  bar at the frame column the expression already produced; the original
  expression is kept for labels and for newdata prediction, and `:`
  / `/` groupings keep their reformulas expansion. [lme4#464, #156]
- The quadrature defect cluster (dev/fuzz-findings.md N1-N4). TMBad's
  `marginal_gk` rescales each integrand ONCE, at whatever parameter
  values the template holds when `MakeADFun` tapes it, and freezes that
  `(mu, sigma)` pair. One cold calibration explained all of it: every
  conditional mode but the first came back `NA` (the marginalized
  objective carries none, and `parList()` slid an outer value into the
  first slot), and poisson, Gamma and Beta over nested scalar blocks -
  Beta over a single one - died at a bare `NA/NaN gradient evaluation`.
  `frm()` now fits the plain Laplace objective first, tapes the
  marginalized one at that optimum, and reads the modes back from the
  inner Newton solve there. They match `glmer(nAGQ = 25)`'s `ranef()`
  to 3e-05. `quadrature` crossed with `trunc()` is refused instead:
  the normalizer `log(F(ub) - F(lb))` underflows at the Gauss-Kronrod
  nodes, so the objective is `-Inf` even at the Laplace optimum, and
  the fit used to report `logLik = +Inf` as converged. What survives is
  a runtime limitation on hard likelihoods (singular variance
  components, nested blocks on thin data), reported as an error naming
  `quadrature` rather than an RTMB string. Regression tests in
  tests/testthat/test-quadrature-defects.R.
- `mixture()` under `REML = TRUE` or `frmtmb_control(profile = TRUE)`
  is refused. Both integrate the fixed effects out with a Laplace
  approximation about a single inner mode, and a mixture likelihood is
  invariant to permuting its components, so it is multimodal in exactly
  those coefficients. The fits used to stop at `NA/NaN gradient
  evaluation` or report a gradient near 1e9 with no guard.
  `quadrature = TRUE` stays allowed - it marginalizes the random
  effects, not the coefficients - and test-v19.R pins down that it is
  exact when the per-group integrand is univariate.
- Truncation reached only the likelihood, never the post-fit surface.
  `fitted()`, `predict(type = "response")` and `residuals()` now report
  E[Y | lb <= Y <= ub] (closed forms for every family with an `lcdf`:
  gaussian, lognormal, poisson, exponential, weibull, inverse.gaussian;
  the poisson form reuses the objective's F(lb-1) convention), and
  `simulate()`/`posterior_predict()` reject out-of-bounds draws instead
  of sampling the untruncated distribution (23% of draws used to land
  outside [lb, ub] on a Poisson trunc(2, 6) fit, which silently
  invalidated DHARMa and pp_check). Dpar-scale predictions stay
  untruncated on purpose. `residuals(type = "osa")` was wrong too: the
  taped density integrates to 1 only over [lb, ub], so the conditional
  CDF is now built on that domain (and never with `fullGaussian`, since
  a truncated gaussian is not gaussian). [brms#1923, #1903]

### Mitigated

- `ar1()`/`hetar1()` over an ordering factor with gaps still treat level
  POSITION as time - dropping times 7-9 from a 1..10 series makes
  cor(t6, t10) come back as rho, not rho^4, and biases rho itself. That
  reading is glmmTMB's, and the glmmTMB agreement tests pin it down, so
  the likelihood is unchanged and frame assembly warns instead: when the
  ordering levels are whole numbers but not consecutive, the warning
  names the gap and points at `ou()` over `num_factor()`, which is the
  correct spelling for irregular spacing. Non-integer labels stay
  silent, since position is then the only available meaning.
  [glmmTMB#1278]

### Open - medium

- `quadrature = TRUE` still breaks down on hard likelihoods: a variance
  component near zero (the fuzzer found `sd = 8.8e-05`) defeats the
  finite-difference curvature estimate `marginal_gk` calibrates with
  (`dx = 1`), and nested scalar blocks on thin data make the outer
  integrand the output of a frozen inner rescaling. `quad_fit()` tries
  three calibration points and then errors, naming `quadrature`. A real
  fix wants an integrator that recalibrates per evaluation; TMBad's
  `adaptive = TRUE` is meant to be that and is measurably worse, so it
  would have to be built rather than switched on.
- `cbind(successes, failures)` responses are rejected with a message
  about `trials` that does not name the real problem. Every reference
  package accepts the spelling, so ported code lands here first.
  [glmmTMB#1319, #1325]
- A model with zero free parameters (`y | trials(n) ~ 0`) fails inside
  nlminb with "'d' must be a nonempty numeric (double) vector" rather
  than reporting the degenerate model. [glmmTMB#1325, #1317]
- No `link_zi` / `link_hu`: the zero-inflation and hurdle parts are
  logit-only, as in glmmTMB. [glmmTMB#847]

### Verified immune (regression tests added)

- Zero prior weights are exactly equivalent to subsetting. [lme4#880]
- `(1 | a * b)` expands to a + b + a:b; lme4 still cannot. [lme4#234]
- Nonlinear fixed-effect SEs match nlme; nlmer is ~100x too small on
  the same fit, and `predict(newdata =)` works. [lme4#819, #164]
- REML predictions agree with `fixef()`. [glmmTMB#1143, #983]
- Numeric vs character grouping levels in newdata. [lme4#616]
- Non-integer binomial/poisson responses rejected. [lme4#682, #180]
- `t2()` matches mgcv; `te()`/`ti()` refused clearly. [glmmTMB#1082]
- Discrete truncation normalizes with F(lb-1). [brms#1903, #1923]

## OSA and inference surface (2026-09-01)

From the compatibility-registry probes and the grammar fuzzer.
Regression tests in tests/testthat/test-osa-inference.R.

### Fixed

- `cens()` x `residuals(type = "osa")`. A censored row contributes a
  probability MASS, and in the tape that contribution does not depend
  on the observation at all, so `fullGaussian` inverted an exactly
  singular Hessian block and `oneStepGeneric` integrated a flat slice
  to infinity and returned NaN. What is well defined is the CDF of the
  UNCENSORED rows conditional on the censoring events, so those rows go
  in `subset` and the censored rows in `conditional`; an uncensored row
  is a draw that landed inside the censoring window, so its integration
  domain is that window, exactly as a `trunc()` fit's is. Censored rows
  return NA. Verified against the analytic conditional PIT
  `qnorm(F(y) / F(c))` to 7e-09 and by KS uniformity. Row-varying
  censoring points (the distribution of an uncensored response is then
  not identified without a model for the censoring process) and
  interval censoring are refused with a message that also says why
  `dharma_residuals()` is not the fallback: `simulate()` draws the
  LATENT uncensored response, so its draws are not comparable with the
  observed censored values.
- `cens()` x `simulate()` semantics. Drawing the latent response is
  CORRECT and is what brms does: no `posterior_predict_*` method in
  brms 2.23.0 reads the censoring column (only `log_lik_censor` does),
  and brms's `pp_check()` therefore drops the censored rows outright
  ("Censored responses are not included"). The model describes the
  latent distribution; censoring belongs to the observation process.
  Documented as the default, with `simulate(censored = TRUE)` as the
  opt-in that applies the mechanism: every draw is recorded at the
  edge of the observation window, so the draws become comparable with
  the observed data. That needs one censoring point per side (type-I
  censoring), because an uncensored row's censoring point is unknown
  when the times vary by row; row-varying times and interval censoring
  are refused.
- Ordinal x `residuals(type = "osa")`. `oneStepPredict` re-tapes with
  the response promoted to a parameter, and the ordinal lpdfs index and
  compare with it (`y == K`, `tau[pmin(y, K1)]`), which no advector
  supports. The four ordinal families now carry an OSA branch that
  selects the category with a Lagrange basis over 1..K - exact in
  floating point at integer y - and applies the `@keep` data-term
  indicator that RTMB's own densities get from `dGenericOSA`. The
  residuals match the analytic randomized quantile residual to 4e-14
  and are uniform under KS for cumulative, sratio, cratio and acat,
  with or without random effects.
- Raw LAPACK error from `solve(jointPrecision)` under REML or
  `profile = TRUE`. See fuzz finding N4.

### Open - medium

- OSA under `cens()` covers the uncensored rows only. A residual for a
  censored row would have to be a randomized quantile inside the
  censoring interval, which `oneStepPredict` cannot produce; doing it
  would mean computing the conditional CDF at the censoring point
  outside TMB.
- `dharma_residuals()` on a censored fit compares latent draws with
  observed censored values and is not valid. It neither warns nor
  refuses. `simulate(censored = TRUE)` makes the draws comparable, but
  the point mass it puts at each censoring point is not a distribution
  DHARMa's rank transform can use, so it is not a drop-in fix.

## Reference

Full agent report with per-item repro sketches and issue links:
sourced from glmmTMB/lme4/brms GitHub test suites and NEWS, 2026-08-31.
Key files to mirror: glmmTMB test-predict.R, test-formulas.R,
test-varstruc.R, test-offset.R, test-weight.R, test-NAhandling.R (lme4),
tests.brmsformula.R (brms).
