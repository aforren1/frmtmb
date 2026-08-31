# Package index

## Model specification

- [`bf()`](bf.md) : Set up a model formula
- [`mvbf()`](mvbf.md) [`set_rescor()`](mvbf.md) : Combine formulas into
  a multivariate model
- [`num_factor()`](num_factor.md) : Factor with numeric-coded levels for
  coordinate covariance structures

## Fitting

- [`frm()`](frm.md) : Fit a model
- [`frmtmb_control()`](frmtmb_control.md) : Control parameters for
  frmtmb fits

## Families

- [`frmtmb_family()`](frmtmb_family.md)
  [`custom_family()`](frmtmb_family.md) : Define a model family
- [`student()`](frmtmb-families.md) [`lognormal()`](frmtmb-families.md)
  [`negbinomial()`](frmtmb-families.md)
  [`nbinom1()`](frmtmb-families.md) [`Beta()`](frmtmb-families.md)
  [`tweedie()`](frmtmb-families.md) [`compois()`](frmtmb-families.md)
  [`zero_inflated_poisson()`](frmtmb-families.md)
  [`zero_inflated_negbinomial()`](frmtmb-families.md)
  [`hurdle_poisson()`](frmtmb-families.md)
  [`multinomial()`](frmtmb-families.md)
  [`cumulative()`](frmtmb-families.md)
  [`beta_binomial()`](frmtmb-families.md)
  [`skew_normal()`](frmtmb-families.md)
  [`exgaussian()`](frmtmb-families.md) : Additional response families
- [`check_custom_family()`](check_custom_family.md) : Check a custom
  family's log-density for AD safety

## Inference and diagnostics

- [`confint(`*`<frmtmb_fit>`*`)`](confint.frmtmb_fit.md) : Confidence
  intervals for frmtmb fits
- [`confint_varcorr()`](confint_varcorr.md) : Natural-scale confidence
  intervals for covariance parameters
- [`diagnose()`](diagnose.md) : Convergence diagnostics for a frmtmb fit
- [`dharma_residuals()`](dharma_residuals.md) : DHARMa residual
  diagnostics
- [`residuals(`*`<frmtmb_fit>`*`)`](residuals.frmtmb_fit.md) : Residuals
  from a frmtmb fit
- [`anova(`*`<frmtmb_fit>`*`)`](anova.frmtmb_fit.md) : Likelihood-ratio
  tests between nested frmtmb fits
- [`as_tmbstan()`](as_tmbstan.md) : Sample from a frmtmb fit with
  tmbstan (NUTS)
- [`frm_sample()`](frm_sample.md) : Sample the fitted model with NUTS
- [`check_laplace()`](check_laplace.md) : Check the Laplace/Wald
  approximation against NUTS
- [`set_prior()`](set_prior.md) : Set up priors brms-style
- [`prior_normal()`](frmtmb-priors.md) [`prior_t()`](frmtmb-priors.md) :
  Prior specifications for frm_sample

## Extractors and prediction

- [`predict(`*`<frmtmb_fit>`*`)`](predict.frmtmb_fit.md) : Predictions
  from a frmtmb fit
- [`residuals(`*`<frmtmb_fit>`*`)`](residuals.frmtmb_fit.md) : Residuals
  from a frmtmb fit
- [`simulate(`*`<frmtmb_fit>`*`)`](simulate.frmtmb_fit.md) : Simulate
  responses from a frmtmb fit
- [`fixef()`](fixef.md) : Extract fixed effects
- [`ranef()`](ranef.md) : Extract random-effect modes
- [`VarCorr()`](VarCorr.md) : Extract random-effect covariance matrices
- [`rescor_matrix()`](rescor_matrix.md) : Estimated residual correlation
  matrix (rescor fits), else NULL
- [`vcov(`*`<frmtmb_fit>`*`)`](vcov.frmtmb_fit.md) : Covariance matrix
  of the fixed-effect estimates
