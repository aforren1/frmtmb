#' srr_stats
#'
#' rOpenSci statistical software standards for frmtmb.
#'
#' Categories: **General** and **Regression and Supervised Learning**.
#'
#' The Bayesian and Monte Carlo category does not apply. frmtmb is a
#' maximum-likelihood package: parameters are point estimates found by
#' numerical optimization, and random effects are integrated out with the
#' Laplace approximation, not sampled. Inference is Wald, profile
#' likelihood, likelihood-root, or nonparametric bootstrap. `set_prior()`
#' produces penalized (MAP) point estimates, which is regularization,
#' not Bayesian inference.
#'
#' The optional `frm_sample()` sampler runs NUTS through tmbstan. On a
#' fitted model it explores the same objective the fit maximized, with
#' flat priors, which is a diagnostic for the Laplace and Wald
#' approximations (see `check_laplace()`) rather than the inferential
#' path. From a formula it does sample a posterior, under weakly
#' informative default priors matching brms, and it reports `n_eff` and
#' `Rhat` for every parameter; it is nonetheless an adjunct rather than
#' the package's inferential surface, which is why the Bayesian
#' standards' requirements on prior specification, posterior output and
#' chain-convergence contracts are not claimed here. A model whose
#' primary answer is a posterior is better served by brms and Stan.
#'
#' The Time Series, Spatial, Dimensionality Reduction, Machine Learning,
#' Probability Distributions, Unsupervised Learning, and Exploratory Data
#' Analysis categories also do not apply. The package fits spatial
#' covariance structures (`car()`, `spde()`) and correlated-error
#' structures (`ar1()`, `ou()`) as random-effect covariance vocabulary
#' inside a regression model, not as a spatial or time-series analysis
#' front end. There is no train/test split, resampling protocol, or
#' predictive-performance workflow of the kind the Machine Learning
#' category assumes.
#'
#' Standards that are met are tagged with `@srrstats` at the code that
#' implements them. Standards that a document rather than code satisfies
#' are tagged in that document instead, in a dedicated non-evaluated
#' chunk: `vignettes/inputs.Rmd` carries the input, terminology,
#' attribute, scaling and benchmark-reproduction standards it closes.
#' Standards that do not apply are collected below with a reason.
#' Standards that are not yet met carry no tag; those are the compliance
#' gap list.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G5.2a} Every condition message raised in `R/` is unique.
#'   The 646 `stop()` calls that carry literal text produce 646 distinct
#'   messages (mechanically re-counted by an AST walk; 529 at v0.35.0,
#'   plus nine net from the structured simulator contract and the
#'   formula interface of `frm_sample()` with its default priors, plus
#'   32 from the `brmsfit` method surface on draws (the leave-one-out
#'   cluster and the principled refusals that replace "could not find
#'   function" for a ported brms script), plus 19 from the `frm_ode()`
#'   pharmacometrics tier, plus two from the influence plot and the
#'   ambiguous bare nonlinear parameter, plus four from
#'   `conditional_effects()` on draws and the brms::bf() masking
#'   detection, plus one for the non-centered sampling switch
#'   (`frm_sample(reparameterize = )` rejecting a non-logical value),
#'   plus six for the LKJ correlation prior (its shape argument, the
#'   two directions of the class-and-distribution check in
#'   `set_prior()`, a bound offered to a whole correlation matrix, a
#'   class `"cor"` that matches no block of the model, and the by-name
#'   prior spelling refusing a density that belongs to a whole block)),
#'   so a
#'   message a user reports names one line of source. Two calls that
#'   would otherwise read the
#'   same are separated by the context that tells the two faults apart,
#'   for example `car()` on the left of a bar term versus `car()` as its
#'   grouping factor, `posterior_predict()` versus `simulate()` versus
#'   `frm_simulate()` meeting a family with no simulator - each names
#'   itself and then repeats the family's own reason, so the three read
#'   consistently and stay distinguishable - and `confint(parm =)`
#'   versus `profile(parm =)` rejecting an unknown parameter name.
#'   Where one refusal serves many callers it is written once and told
#'   which caller it speaks for rather than copied: `require_fitted()`
#'   raises a single message naming the method that called it, which
#'   keeps the source line unique while the user still reads their own
#'   call back. The property is mechanically checkable: parse `R/`,
#'   collect the literal argument text of every `stop()`, and require
#'   no duplicates. The count is a whole-package property, so a
#'   development lane that adds or removes a refusal must re-run the
#'   walk and update the number here; concurrent lanes editing this
#'   file therefore conflict on this paragraph by design.
#'
#' @noRd
NULL

#' NA_standards
#'
#' Standards that do not apply to frmtmb, with reasons.
#'
#' @srrstatsNA {G2.4d} The package never converts an input to `factor`.
#'   Factor structure belongs to the user's data and is resolved by
#'   `stats::model.frame()` and `stats::model.matrix()` with the stored
#'   contrasts. Creating factors silently would change the model.
#' @srrstatsNA {G3.1} The package estimates no covariance matrix from
#'   data with `stats::cov()`, which it never calls. The covariance of
#'   the parameters is the inverse observed information taken from the
#'   automatic-differentiation Hessian, and residual correlations
#'   (`rescor`) and random-effect covariances are model parameters
#'   estimated by maximum likelihood. There is therefore no covariance
#'   algorithm to choose between.
#' @srrstatsNA {G3.1a} Follows from G3.1: there is no covariance method
#'   argument to document.
#' @srrstatsNA {G4.0} The package writes no files. No exported function
#'   takes a file name or path, and `R/` contains no `write*()`,
#'   `saveRDS()`, or `sink()` call.
#' @srrstatsNA {G5.1} The package creates and exports no data sets. Tests
#'   use data sets from other packages (`lme4::sleepstudy`,
#'   `lme4::cbpp`, `datasets::faithful`, `nlme::Soybean`,
#'   `brms::inhaler`, `brms::kidney`) or seeded simulators defined in
#'   `tests/testthat/helper-reference.R`, which ships in the package
#'   tarball, so every test data set can be reproduced by a user.
#' @srrstatsNA {G5.4c} No method in the package comes from a published
#'   numerical example whose implementation is unavailable. Every
#'   reference implementation is an installable R package and is called
#'   live in the tests, so stored published values are not needed.
#' @srrstatsNA {G5.11} The extended tests need no additional data. They
#'   simulate their own data from a seed, so there is nothing to
#'   download.
#' @srrstatsNA {G5.11a} Follows from G5.11: there are no downloads that
#'   could fail.
#' @srrstatsNA {RE4.15} frmtmb is not forecasting software. It fits
#'   regression models to observed data; there is no time index, no
#'   forecast horizon, and no extrapolation beyond the predictor space
#'   that a horizon could parameterize. Prediction uncertainty away from
#'   the observed data is still reported: `predict(se.fit = TRUE)` adds
#'   the marginal variance of unseen grouping levels and the kriging
#'   variance of a Gaussian process.
#' @srrstatsNA {RE6.3} Follows from RE4.15. There is no forecast, so a
#'   plot method has no interpolated/extrapolated distinction to draw.
#'   Non-estimable prediction rows under a rank-deficient fit are marked
#'   instead, by returning `NA` with a warning.
#' @srrstatsNA {RE7.4} Follows from RE4.15: there is no forecast horizon
#'   along which to demonstrate growing intervals.
#' @noRd
NULL
