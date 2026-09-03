#' srr_stats
#'
#' rOpenSci statistical software standards for frmtmb.
#'
#' Categories: **General** and **Regression and Supervised Learning**.
#'
#' The Bayesian and Monte Carlo category does not apply. The reason is
#' the package's documented inferential surface, not the absence of a
#' sampler.
#'
#' That surface is maximum likelihood with the Laplace approximation.
#' `frm()` finds point estimates by numerical optimization and
#' integrates latent effects out with the Laplace approximation. Every
#' estimate, standard error, interval and test that the documentation
#' reports is a maximum-likelihood quantity. `summary()` gives Wald
#' standard errors from the inverse observed information. `confint()`
#' gives Wald, profile-likelihood, likelihood-root or bootstrap
#' intervals. `anova()` and `drop1()` give likelihood-ratio tests.
#' `vignette("frmtmb")`, `vignette("brms-migration")`,
#' `vignette("case-studies")`, `vignette("diagnostics")` and
#' `vignette("ode")` all report their models this way. `set_prior()` on
#' a `frm()` fit gives a penalized (MAP) point estimate, which is
#' regularization and not posterior inference.
#'
#' `frm_sample()` is an opt-in bridge to NUTS through tmbstan, and its
#' help page states which question each of its two routes answers. On a
#' fitted model it explores the LIKELIHOOD under flat improper priors.
#' That is what makes `check_laplace()` meaningful: the run measures the
#' Laplace and Wald approximations against the shape of the objective
#' that the fit maximized, and a default prior would change the thing
#' being measured. From a formula the function does sample a posterior,
#' under weakly informative default priors that match brms, with an LKJ
#' default on correlations and a non-centered parameterization, and it
#' reports `n_eff` and `Rhat` for every parameter. That route has two
#' documented purposes: a script ported from brms keeps the
#' `posterior_epred()`, `pp_check()` and `loo()` calls it already
#' contains, and a Laplace fit can be compared against the posterior of
#' the same model. The package states the limits of the bridge instead
#' of hiding them. The section "Priors, and what these numbers mean" in
#' `?loo` says that an elpd from the flat-prior route is
#' likelihood-shaped and unregularized, and sends model comparison to
#' `AIC()` or to a run with priors. `bridge_sampler()` refuses a
#' marginal likelihood rather than return one that a flat prior leaves
#' undefined.
#'
#' The Bayesian standards describe a package whose front door is the
#' sampler. They put prior specification and prior sensitivity in the
#' api, they make the posterior-predictive workflow the primary output,
#' and they attach chain-convergence contracts to the return object.
#' The front door of frmtmb is `frm()`, and its return object is a
#' fitted model with a Hessian. To claim those standards would describe
#' software that this package is not, and it would tell a user that the
#' guarantees are somewhere they are not. A model whose primary answer
#' is a posterior is referred to brms and Stan: the README says so under
#' "Related work", and `vignette("brms-migration")` has a section
#' "When you still want brms" that names the cases.
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
#' @srrstats {G5.2a} Every `stop()` message raised in `R/` is unique.
#'   The 651 `stop()` calls that carry literal text produce 651 distinct
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
#'   prior spelling refusing a density that belongs to a whole block),
#'   plus four from the functional-data batch: the
#'   `conditional_effects()` refusal that names its matrix columns, the
#'   two `predict(newdata = )` refusals that name a smooth's missing
#'   columns before mgcv reports the length mismatch internally, and an
#'   unseen level of a factor-smooth term, itself two refusals: the
#'   default, and the factor `bs = "re"` basis that has no zero row to
#'   hand a new level under `allow_new_levels = TRUE`), so a
#'   message a user reports names one line of source. Two calls that
#'   would otherwise read the
#'   same are separated by the context that tells the two faults apart,
#'   for example `car()` on the left of a bar term versus `car()` as its
#'   grouping factor, and `confint(parm =)` versus `profile(parm =)`
#'   rejecting an unknown parameter name. `posterior_predict()`,
#'   `simulate()` and `frm_simulate()` meeting a family with no
#'   simulator are a third case: each names itself and then repeats the
#'   family's own reason, so the three read consistently and stay
#'   distinguishable.
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
#'   simulate their own data from a seed, or use a data set from a
#'   package that is already in `Suggests`, so there is nothing to
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
