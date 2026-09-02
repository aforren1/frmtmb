# Package index

## Model specification

- [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) : Set up a
  model formula

- [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md)
  [`set_rescor()`](https://aforren1.github.io/frmtmb/reference/mvbf.md)
  : Combine formulas into a multivariate model

- [`num_factor()`](https://aforren1.github.io/frmtmb/reference/num_factor.md)
  : Factor with numeric-coded levels for coordinate covariance
  structures

- [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
  : Solve an ODE once per group inside a nonlinear predictor

- [`frm_ode_failures()`](https://aforren1.github.io/frmtmb/reference/frm_ode_failures.md)
  :

  Groups whose ODE solve failed in the last
  [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
  call

## Fitting

- [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) : Fit a
  model
- [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md)
  : Control parameters for frmtmb fits
- [`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)
  : Refit a model with every available optimizer
- [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
  : Simulate responses from a formula and parameters
- [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  : Fit a model across multiply imputed datasets

## Families

- [`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  : Define a model family
- [`student()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`lognormal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`negbinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`nbinom1()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`Beta()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`tweedie()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`compois()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`zero_inflated_poisson()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`zero_inflated_negbinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`hurdle_poisson()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`cumulative()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`beta_binomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`skew_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`exgaussian()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`bernoulli()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`geometric()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`exponential()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`weibull()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`shifted_lognormal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`hurdle_gamma()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`hurdle_lognormal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`zero_inflated_binomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`zero_inflated_beta()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`asym_laplace()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`sratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`cratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  : Additional response families
- [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  : Finite mixture families
- [`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
  : Multivariate gaussian mixture family
- [`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
  : Posterior class probabilities of a mixture fit
- [`check_custom_family()`](https://aforren1.github.io/frmtmb/reference/check_custom_family.md)
  : Check a custom family's log-density for AD safety

## Inference and diagnostics

- [`confint(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/confint.frmtmb_fit.md)
  : Confidence intervals for frmtmb fits
- [`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
  : Natural-scale confidence intervals for covariance parameters
- [`profile(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/profile.frmtmb_fit.md)
  : Likelihood profiles
- [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  : Hypothesis tests on parameter expressions
- [`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
  : Usable parameter names
- [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  : Parametric bootstrap
- [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
  : Convergence diagnostics for a frmtmb fit
- [`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
  : DHARMa residual diagnostics
- [`residuals(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/residuals.frmtmb_fit.md)
  : Residuals from a frmtmb fit
- [`plot(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/plot.frmtmb_fit.md)
  : Diagnostic plots for a fit
- [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
  : Predictive check against simulated responses
- [`anova(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/anova.frmtmb_fit.md)
  : Likelihood-ratio tests between nested frmtmb fits
- [`anova(`*`<frmtmb_multiple>`*`)`](https://aforren1.github.io/frmtmb/reference/anova.frmtmb_multiple.md)
  : Pooled model comparison across imputations (D1, D2, D3)
- [`drop1(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/drop1.frmtmb_fit.md)
  : Single-term deletions
- [`as_tmbstan()`](https://aforren1.github.io/frmtmb/reference/as_tmbstan.md)
  : Sample from a frmtmb fit with tmbstan (NUTS)
- [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  : Sample the fitted model with NUTS
- [`posterior_epred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  : Expected-value and predictive draws from sampled parameters
- [`as_draws()`](https://aforren1.github.io/frmtmb/reference/as_draws.md)
  : Convert draws to a posterior draws object
- [`getME(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/getME.frmtmb_fit.md)
  : Extract components of a fit, lme4 style
- [`check_laplace()`](https://aforren1.github.io/frmtmb/reference/check_laplace.md)
  : Check the Laplace/Wald approximation against NUTS
- [`influence(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/influence.frmtmb_fit.md)
  [`cooks.distance(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/influence.frmtmb_fit.md)
  [`dfbeta(`*`<frmtmb_influence>`*`)`](https://aforren1.github.io/frmtmb/reference/influence.frmtmb_fit.md)
  [`dfbetas(`*`<frmtmb_influence>`*`)`](https://aforren1.github.io/frmtmb/reference/influence.frmtmb_fit.md)
  : Influence measures by case deletion
- [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  : Set up priors brms-style
- [`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md)
  : Enumerate the targetable prior slots
- [`prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
  [`prior_t()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
  : Prior specifications for frm_sample
- [`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.md)
  : Priors used in a fit

## Extractors and prediction

- [`predict(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
  : Predictions from a frmtmb fit
- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  : Conditional effects of predictors
- [`simulate(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.md)
  : Simulate responses from a frmtmb fit
- [`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md) :
  Refit a model to a new response
- [`coef(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/coef.frmtmb_fit.md)
  : Per-group coefficients (fixed effects plus conditional modes)
- [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) :
  Extract fixed effects
- [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md) :
  Extract random-effect modes
- [`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.md) :
  Number of levels per random-effect grouping factor
- [`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
  : Extract random-effect covariance matrices
- [`rescor_matrix()`](https://aforren1.github.io/frmtmb/reference/rescor_matrix.md)
  : Estimated residual correlation matrix (rescor fits), else NULL
- [`vcov(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/vcov.frmtmb_fit.md)
  : Covariance matrix of the fixed-effect estimates
- [`sigma(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/sigma.frmtmb_fit.md)
  : Residual standard deviation

## Feature compatibility

- [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
  : Query the feature compatibility registry
- [`frm_compat_features()`](https://aforren1.github.io/frmtmb/reference/frm_compat_features.md)
  : Feature metadata for the compatibility registry
- [`frm_compat_rules()`](https://aforren1.github.io/frmtmb/reference/frm_compat_rules.md)
  : Compatibility rules, before resolution

## Package

- [`frmtmb`](https://aforren1.github.io/frmtmb/reference/frmtmb-package.md)
  [`frmtmb-package`](https://aforren1.github.io/frmtmb/reference/frmtmb-package.md)
  : frmtmb: Formula-Based Regression Models via 'RTMB'
