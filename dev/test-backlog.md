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
regression tests in tests/testthat/test-open-issues.R.

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

### Open - high priority

1. Truncation is ignored by every post-fit mean: the likelihood
   normalizes correctly with F(lb-1), but `fitted()`, `predict(type =
   "response")` and `residuals()` return the *untruncated* family mean,
   and `simulate()` draws from the untruncated distribution (23% of
   draws land outside [lb, ub] on a Poisson trunc(2, 6) fit). brms is
   wrong here too but only by an off-by-one; we omit the correction
   entirely. Needs E[Y | lb <= Y <= ub] per family (only gaussian,
   lognormal and poisson accept trunc()) plus rejection sampling in
   `simulate()`. [brms#1923, #1903]

### Open - medium

- Prediction at an aliased cell of a rank-deficient design returns the
  random-effect contribution alone, i.e. the dropped column silently
  contributes 0. Upstream wants NA or an error; new *fixed* levels
  already error correctly. [lme4#303]
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

## Reference

Full agent report with per-item repro sketches and issue links:
sourced from glmmTMB/lme4/brms GitHub test suites and NEWS, 2026-08-31.
Key files to mirror: glmmTMB test-predict.R, test-formulas.R,
test-varstruc.R, test-offset.R, test-weight.R, test-NAhandling.R (lme4),
tests.brmsformula.R (brms).
