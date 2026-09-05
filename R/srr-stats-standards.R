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
#' `vignette("case-studies")` and `vignette("diagnostics")` all report
#' their models this way. `set_prior()` on
#' a `frm()` fit gives a penalized (MAP) point estimate, which is
#' regularization and not posterior inference.
#'
#' Sampling is not in this package. The companion package
#' frmtmb.sample provides `frm_sample()`, an opt-in bridge to NUTS
#' through tmbstan, together with the posterior method surface it
#' needs (`posterior_epred()`, `pp_check()`, `loo()`) and
#' `check_laplace()`, which measures the Laplace and Wald
#' approximations against a posterior of the objective the fit
#' maximized. On either of its routes, a fitted model or a formula,
#' `frm_sample()` samples under weakly informative default priors that
#' match brms, with an LKJ default on correlations and a non-centered
#' parameterization, and it reports `n_eff` and `Rhat` for every
#' parameter. This package keeps only the generics and the registration
#' seam those methods attach to (`?"frmtmb-sampling-api"`), so that a
#' script ported from brms keeps the calls it already contains once the
#' companion is attached. The companion states the limits of the bridge
#' instead of hiding them: its `loo()` help page says that an elpd from
#' a `prior = "flat"` run is likelihood-shaped and unregularized, and
#' sends model comparison to `AIC()` or to a run with priors, and its
#' `bridge_sampler()` refuses a marginal likelihood that a flat prior
#' leaves undefined.
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
#' @srrstats {G5.2a} Every condition message raised in `R/` is unique,
#'   for `stop()`, `warning()` and `message()` alike, and the property
#'   is ENFORCED rather than recorded:
#'   `tests/testthat/test-message-uniqueness.R` parses `R/` at test
#'   time, collects the literal string fragments of every condition
#'   call, and fails on the first duplicate template. No count is kept
#'   in prose on purpose; a number here was stale the commit after it
#'   was written, and the test fails at the moment of drift instead.
#'
#'   What the walk certifies is TEMPLATE uniqueness: each message's
#'   literal text resolves to one line of source, so a reported message
#'   is findable. Two templates that interpolate different runtime
#'   values can still render similar final text; the shared
#'   input-validation helpers in `R/utils.R` are the known case, one
#'   template each with the offending argument's name filled in at run
#'   time, so `nsim` misuse reads the same from `simulate()` and
#'   `frm_bootstrap()`. Where the caller matters to the user, the
#'   template carries it: `require_fitted()` names the method that
#'   called it, `car()` on the left of a bar term reads differently
#'   from `car()` as a grouping factor, and `posterior_predict()`,
#'   `simulate()` and `frm_simulate()` meeting a family with no
#'   simulator each name themselves before repeating the family's own
#'   reason.
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
