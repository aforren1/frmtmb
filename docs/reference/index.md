# Package index

## Model specification

- [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) : Set up a
  model formula
- [`lf()`](https://aforren1.github.io/frmtmb/reference/lf.md) : Add
  parameter formulas to a model formula
- [`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) : Add a
  nonlinear parameter formula to a model formula
- [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md)
  [`set_rescor()`](https://aforren1.github.io/frmtmb/reference/mvbf.md)
  : Combine formulas into a multivariate model
- [`num_factor()`](https://aforren1.github.io/frmtmb/reference/num_factor.md)
  : Factor with numeric-coded levels for coordinate covariance
  structures
- [`frmtmb-autocor`](https://aforren1.github.io/frmtmb/reference/frmtmb-autocor.md)
  : Within-group residual correlation (R-side autocorrelation)
- [`frmtmb-multimembership`](https://aforren1.github.io/frmtmb/reference/frmtmb-multimembership.md)
  : Multi-membership random effects

## Fitting

- [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) : Fit a
  model
- [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md)
  : Control parameters for frmtmb fits
- [`par_template()`](https://aforren1.github.io/frmtmb/reference/par_template.md)
  : Parameter names and starting values
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

- [`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
  : Declare a non-rowwise likelihood to the core

- [`frmtmb_register_frame_check()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_frame_check.md)
  : Check an assembled model frame from another package

- [`frmtmb_ad_overload()`](https://aforren1.github.io/frmtmb/reference/frmtmb_ad_overload.md)
  : Make a user-written function tape-safe

- [`response_mean()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`mixture_multimodal_refusals()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`mixture_posterior()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`as_frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`eval_dpars()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`single_response()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`fit_extras()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`dpar_linpred()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`structure_supports_all()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  [`frame_block_of()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
  : The two fitting options every mixture-type family refuses, in its
  own name.

- [`latent_probs()`](https://aforren1.github.io/frmtmb/reference/latent_probs.md)
  : Posterior probabilities of the latent states

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
  [`zero_inflated_asym_laplace()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`huber()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`sratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`cratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`von_mises()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`categorical()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  [`cox()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  : Additional response families

- [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  : Finite mixture families

- [`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
  : Multivariate gaussian mixture family

- [`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
  : Posterior class probabilities of a mixture fit

- [`cox_baseline()`](https://aforren1.github.io/frmtmb/reference/cox_baseline.md)
  :

  The fitted baseline-hazard simplex of a
  [`cox()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  fit.

- [`frmtmb-student-re`](https://aforren1.github.io/frmtmb/reference/frmtmb-student-re.md)
  : Student-t distributed random effects

- [`check_custom_family()`](https://aforren1.github.io/frmtmb/reference/check_custom_family.md)
  : Check a custom family's log-density for AD safety

## Inference and diagnostics

- [`confint(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/confint.frmtmb_fit.md)
  : Confidence intervals for frmtmb fits
- [`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
  : Natural-scale confidence intervals for covariance parameters
- [`vcov_cluster()`](https://aforren1.github.io/frmtmb/reference/vcov_cluster.md)
  : Cluster-robust (sandwich) covariance
- [`cluster_scores()`](https://aforren1.github.io/frmtmb/reference/cluster_scores.md)
  : Per-cluster score matrix
- [`autocor_matrix()`](https://aforren1.github.io/frmtmb/reference/autocor_matrix.md)
  : Estimated within-group residual correlation matrix
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
- [`as_draws()`](https://aforren1.github.io/frmtmb/reference/as_draws.md)
  [`as_draws_matrix()`](https://aforren1.github.io/frmtmb/reference/as_draws.md)
  [`as_draws_array()`](https://aforren1.github.io/frmtmb/reference/as_draws.md)
  [`as_draws_df()`](https://aforren1.github.io/frmtmb/reference/as_draws.md)
  [`as_draws_list()`](https://aforren1.github.io/frmtmb/reference/as_draws.md)
  [`as_draws_rvars()`](https://aforren1.github.io/frmtmb/reference/as_draws.md)
  : Convert to a posterior draws object
- [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md)
  [`waic()`](https://aforren1.github.io/frmtmb/reference/loo.md)
  [`loo_compare()`](https://aforren1.github.io/frmtmb/reference/loo.md)
  [`LOO()`](https://aforren1.github.io/frmtmb/reference/loo.md)
  [`WAIC()`](https://aforren1.github.io/frmtmb/reference/loo.md) :
  Approximate leave-one-out cross-validation
- [`bayes_R2()`](https://aforren1.github.io/frmtmb/reference/bayes_R2.md)
  : Bayesian R-squared
- [`posterior_summary()`](https://aforren1.github.io/frmtmb/reference/posterior_summary.md)
  : Summaries and intervals of draws
- [`expose_functions()`](https://aforren1.github.io/frmtmb/reference/expose_functions.md)
  : Expose a model's compiled functions
- [`ndraws()`](https://aforren1.github.io/frmtmb/reference/draws-dimensions.md)
  [`nchains()`](https://aforren1.github.io/frmtmb/reference/draws-dimensions.md)
  [`niterations()`](https://aforren1.github.io/frmtmb/reference/draws-dimensions.md)
  [`nvariables()`](https://aforren1.github.io/frmtmb/reference/draws-dimensions.md)
  : Size of a draws object
- [`getME(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/getME.frmtmb_fit.md)
  : Extract components of a fit, lme4 style
- [`frmtmb-sampling-api`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`build_objective`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`row_lpdf`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`with_cs_offsets`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`us_chol_cor`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`aterms_for_newdata`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`has_trunc`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`as_priorlist`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`resolve_prior_input`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`neg_log_prior_fn`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`resolve_bounds`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`spec_target`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`frmtmb_register_prior_defaults`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ncp_eligible`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ncp_scale_b`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ncp_unscale_b`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`covstruct_has_chol`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`block_sd_idx`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`block_cor_prior`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`block_n_cor`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`is_student_block`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`sim_can`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`sim_note`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`sim_context`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`sim_draw`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`sim_is_structured`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`par_name_bare`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`outer_par_names`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`estimated_coef_names`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`log_sd_theta_index`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`sdr_of`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`require_fitted`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`hyp_parse_all`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`hyp_vals_only`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`hyp_env_vals`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`hyp_eval`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`hyp_tail_p`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ce_grids_build`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ce_boot_one`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ce_finalize`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ce_cats_display`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ce_structure_check`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`ce_re_formula`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`find_linpred`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`arg_unset`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`re_form_arg`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`frmtmb_register_compat`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  [`compat_rule_builder`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
  : Core internals for a sampling extension
- [`influence(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/influence.frmtmb_fit.md)
  [`cooks.distance(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/influence.frmtmb_fit.md)
  [`dfbeta(`*`<frmtmb_influence>`*`)`](https://aforren1.github.io/frmtmb/reference/influence.frmtmb_fit.md)
  [`dfbetas(`*`<frmtmb_influence>`*`)`](https://aforren1.github.io/frmtmb/reference/influence.frmtmb_fit.md)
  : Influence measures by case deletion
- [`plot(`*`<frmtmb_influence>`*`)`](https://aforren1.github.io/frmtmb/reference/plot.frmtmb_influence.md)
  : Plot case-deletion influence
- [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  : Set up priors brms-style
- [`prior()`](https://aforren1.github.io/frmtmb/reference/prior.md)
  [`prior_()`](https://aforren1.github.io/frmtmb/reference/prior.md)
  [`prior_string()`](https://aforren1.github.io/frmtmb/reference/prior.md)
  : Set up priors with brms's quoting spelling
- [`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md)
  : Enumerate the targetable prior slots
- [`prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
  [`prior_t()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
  [`prior_lkj()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
  : Prior objects, addressed by internal parameter name
- [`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.md)
  : Priors used in a fit

## Extractors and prediction

- [`predict(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
  : Predictions from a frmtmb fit
- [`fitted(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/fitted.frmtmb_fit.md)
  : Fitted values
- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  : Conditional effects of predictors
- [`simulate(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.md)
  : Simulate responses from a frmtmb fit
- [`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md) :
  Refit a model to a new response
- [`update(`*`<frmtmb_fit>`*`)`](https://aforren1.github.io/frmtmb/reference/update.frmtmb_fit.md)
  : Update and refit a model
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
