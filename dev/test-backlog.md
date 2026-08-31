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

## Reference

Full agent report with per-item repro sketches and issue links:
sourced from glmmTMB/lme4/brms GitHub test suites and NEWS, 2026-08-31.
Key files to mirror: glmmTMB test-predict.R, test-formulas.R,
test-varstruc.R, test-offset.R, test-weight.R, test-NAhandling.R (lme4),
tests.brmsformula.R (brms).
