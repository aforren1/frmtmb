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
- [`frm_allfit()`](frm_allfit.md) : Refit a model with every available
  optimizer
- [`frm_simulate()`](frm_simulate.md) : Simulate responses from a
  formula and parameters
- [`frm_multiple()`](frm_multiple.md) : Fit a model across multiply
  imputed datasets

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
  [`exgaussian()`](frmtmb-families.md)
  [`bernoulli()`](frmtmb-families.md)
  [`geometric()`](frmtmb-families.md)
  [`exponential()`](frmtmb-families.md)
  [`weibull()`](frmtmb-families.md)
  [`shifted_lognormal()`](frmtmb-families.md)
  [`hurdle_gamma()`](frmtmb-families.md)
  [`hurdle_lognormal()`](frmtmb-families.md)
  [`zero_inflated_binomial()`](frmtmb-families.md)
  [`zero_inflated_beta()`](frmtmb-families.md)
  [`asym_laplace()`](frmtmb-families.md)
  [`sratio()`](frmtmb-families.md) [`cratio()`](frmtmb-families.md)
  [`acat()`](frmtmb-families.md) : Additional response families
- [`mixture()`](mixture.md) : Finite mixture families
- [`mixture_mvn()`](mixture_mvn.md) : Multivariate gaussian mixture
  family
- [`mixture_probs()`](mixture_probs.md) : Posterior class probabilities
  of a mixture fit
- [`check_custom_family()`](check_custom_family.md) : Check a custom
  family's log-density for AD safety

## Inference and diagnostics

- [`confint(`*`<frmtmb_fit>`*`)`](confint.frmtmb_fit.md) : Confidence
  intervals for frmtmb fits
- [`confint_varcorr()`](confint_varcorr.md) : Natural-scale confidence
  intervals for covariance parameters
- [`profile(`*`<frmtmb_fit>`*`)`](profile.frmtmb_fit.md) : Likelihood
  profiles
- [`hypothesis()`](hypothesis.md) : Hypothesis tests on parameter
  expressions
- [`variables()`](variables.md) : Usable parameter names
- [`frm_bootstrap()`](frm_bootstrap.md) : Parametric bootstrap
- [`diagnose()`](diagnose.md) : Convergence diagnostics for a frmtmb fit
- [`dharma_residuals()`](dharma_residuals.md) : DHARMa residual
  diagnostics
- [`residuals(`*`<frmtmb_fit>`*`)`](residuals.frmtmb_fit.md) : Residuals
  from a frmtmb fit
- [`plot(`*`<frmtmb_fit>`*`)`](plot.frmtmb_fit.md) : Diagnostic plots
  for a fit
- [`pp_check(`*`<frmtmb_fit>`*`)`](pp_check.frmtmb_fit.md) : Predictive
  check against simulated responses
- [`anova(`*`<frmtmb_fit>`*`)`](anova.frmtmb_fit.md) : Likelihood-ratio
  tests between nested frmtmb fits
- [`as_tmbstan()`](as_tmbstan.md) : Sample from a frmtmb fit with
  tmbstan (NUTS)
- [`frm_sample()`](frm_sample.md) : Sample the fitted model with NUTS
- [`posterior_epred()`](posterior_epred.md)
  [`posterior_predict()`](posterior_epred.md) : Expected-value and
  predictive draws from sampled parameters
- [`check_laplace()`](check_laplace.md) : Check the Laplace/Wald
  approximation against NUTS
- [`influence(`*`<frmtmb_fit>`*`)`](influence.frmtmb_fit.md) : Influence
  measures by case deletion
- [`set_prior()`](set_prior.md) : Set up priors brms-style
- [`get_prior()`](get_prior.md) : Enumerate the targetable prior slots
- [`prior_normal()`](frmtmb-priors.md) [`prior_t()`](frmtmb-priors.md) :
  Prior specifications for frm_sample
- [`prior_summary()`](prior_summary.md) : Priors used in a fit

## Extractors and prediction

- [`predict(`*`<frmtmb_fit>`*`)`](predict.frmtmb_fit.md) : Predictions
  from a frmtmb fit
- [`conditional_effects()`](conditional_effects.md) : Conditional
  effects of predictors
- [`residuals(`*`<frmtmb_fit>`*`)`](residuals.frmtmb_fit.md) : Residuals
  from a frmtmb fit
- [`simulate(`*`<frmtmb_fit>`*`)`](simulate.frmtmb_fit.md) : Simulate
  responses from a frmtmb fit
- [`refit()`](refit.md) : Refit a model to a new response
- [`coef(`*`<frmtmb_fit>`*`)`](coef.frmtmb_fit.md) : Per-group
  coefficients (fixed effects plus conditional modes)
- [`fixef()`](fixef.md) : Extract fixed effects
- [`ranef()`](ranef.md) : Extract random-effect modes
- [`ngrps()`](ngrps.md) : Number of levels per random-effect grouping
  factor
- [`VarCorr()`](VarCorr.md) : Extract random-effect covariance matrices
- [`rescor_matrix()`](rescor_matrix.md) : Estimated residual correlation
  matrix (rescor fits), else NULL
- [`vcov(`*`<frmtmb_fit>`*`)`](vcov.frmtmb_fit.md) : Covariance matrix
  of the fixed-effect estimates
- [`sigma(`*`<frmtmb_fit>`*`)`](sigma.frmtmb_fit.md) : Residual standard
  deviation
