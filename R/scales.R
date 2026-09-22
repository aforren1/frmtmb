#' Which scale each method reports
#'
#' Every quantity a fitted model reports lives on one of three scales,
#' and nothing in the output says which. This page says it once, per
#' method, so that a ported script can be read without guessing.
#'
#' The three scales are:
#'
#' * **Link.** The linear predictor of one distributional parameter,
#'   before its inverse link is applied. This is where the coefficients
#'   live, where the standard errors are symmetric, and where a random
#'   effect has its variance.
#' * **Response.** The units of the outcome itself.
#' * **Unitless.** A standardized quantity, such as a Pearson or
#'   quantile residual, which has no units by construction.
#'
#' @section What each method returns:
#'
#' | method | returns | scale |
#' |---|---|---|
#' | `predict()` | a summary of the PREDICTIVE distribution | response |
#' | `frm_linpred()`, default | the `mu` linear predictor | link |
#' | `frm_linpred(type = "response")` | the conditional mean | response |
#' | `frm_linpred(dpar = "sigma")` | that predictor | link; `"response"` gives response |
#' | `fitted()` | a summary of the conditional mean | response |
#' | `conditional_effects()` | `estimate__` and its band | response |
#' | `residuals()`, default | observed minus `fitted()` | response |
#' | `residuals(type = "pearson")` | that, over the conditional SD | unitless |
#' | `residuals(type = "deviance")` | signed root unit deviance | unitless |
#' | `residuals(type = "osa")` | one-step-ahead quantile residual | unitless |
#' | `simulate()` | draws of the outcome | response |
#' | `coef()`, `fixef()`, `ranef()` | coefficients, modes | link, per dpar |
#' | `confint()`, `vcov()`, `summary()` | the same | link, per dpar |
#' | `VarCorr()` | SDs and correlations of a predictor | link, per dpar |
#' | `sigma()` | the residual SD, `NA` if it varies by row | response |
#' | `hypothesis()` | `Estimate` and its interval | link, per dpar |
#' | `posterior_summary()` | a summary of a draws MATRIX | that matrix's own |
#' | `pp_check()` | the outcome against replicates | response |
#' | `bayes_R2()` | refuses on a `frmtmb_fit` | neither |
#' | `logLik()`, `AIC()`, `BIC()` | the fitted likelihood | the data's own |
#'
#' @section Where this differs from brms:
#'
#' `predict()` used to be the first, and is no longer: since item 2.6d
#' it is brms's, a summary of the predictive distribution on the
#' response scale. The LINEAR PREDICTOR, which is what it returned
#' before and what the glmmTMB `type` vocabulary reaches, moved to
#' [frm_linpred()] unchanged. On a lognormal fit the two differ by
#' about 8 against about 6,700, which is why they could not share a
#' name. Measured with `dev/generics-scale.R` and
#' `dev/generics-scale-brms.R`.
#'
#' `sigma` is the one that remains. `summary()`'s
#' `Regression Coefficients` block and `fixef()` report every
#' coefficient on its own link, so a `sigma` WITH A FORMULA and the
#' default log link is reported as `log(sigma)` and can be negative.
#' A `sigma` nobody wrote a formula for is not a coefficient at all: it
#' is in `summary()`'s `Further Distributional Parameters` block on its
#' own response scale, as in brms, and [sigma()] returns that number.
#'
#' @section What agrees with brms without any conversion:
#'
#' `fitted()`, `residuals()`, `simulate()`, `ngrps()`, `ranef()`,
#' `conditional_effects()` and the coefficients of a linear predictor
#' whose link is the identity.
#'
#' For a lognormal fit with a CONSTANT sigma and no truncation,
#' `fitted()` is `exp(mu + sigma^2 / 2)`, the mean of the distribution
#' rather than its median `exp(mu)`; the ratio between the two is
#' `exp(sigma^2 / 2)`, which is the whole of the discrepancy a ported
#' script sees. **Both conditions are load-bearing**, and neither is
#' obvious from the formula:
#'
#' * With a distributional sigma (`bf(y ~ x, sigma ~ x)`), `sigma()`
#'   returns `NA` with a warning, because there is no one number to
#'   return, and the formula gives `NA` with it. Use
#'   `frm_linpred(dpar = "sigma", type = "response")`, which reproduces
#'   `fitted(dpar = "sigma")[, "Estimate"]` exactly (`identical()` is
#'   `TRUE`).
#' * Under truncation, `fitted()` is the TRUNCATED mean and the
#'   formula answers a different question. With a lower bound `lb`,
#'   the truncated mean is `exp(mu + sigma^2 / 2) * pnorm(sigma - a) /
#'   pnorm(-a)` with `a = (log(lb) - mu) / sigma`, so the naive
#'   formula falls short by a relative `1 - pnorm(-a) / pnorm(sigma -
#'   a)` in each row. That is an identity, not an estimate, and there
#'   is no single number for it: it depends on where the bound sits in
#'   each row's distribution, near zero for a row far above the bound
#'   and approaching one for a row whose mean is below it
#'   (`dev/generics-trunc.R` sweeps three bounds on one design).
#'
#' Measured on one 400-row lognormal fit against a two-chain brms fit
#' of the same model (`dev/generics-scale-brms.R`, seed 2026): the
#' largest relative disagreement in `fitted()` over the 400 rows is
#' 0.0106, which is 0.113 of brms's own posterior standard deviation
#' for that row and 0.058 of it at the median. `predict()` disagrees
#' by roughly 440 at the median row and 765 at row 1, because it is
#' a different scale rather than a different answer. Those two are
#' given loosely on purpose: brms's `predict()` summarizes fresh
#' predictive draws, so it moves between calls on one fit, where
#' the `fitted()` figures above are deterministic and repeat to the
#' last digit.
#'
#' @seealso [predict.frmtmb_fit()], [frm_linpred()],
#'   [fitted.frmtmb_fit()],
#'   [residuals.frmtmb_fit()], [sigma.frmtmb_fit()], [fixef()],
#'   [VarCorr()], [confint.frmtmb_fit()]
#' @examples
#' set.seed(2026)
#' dd <- data.frame(x = rnorm(200))
#' dd$y <- exp(rnorm(200, 8 + 0.4 * dd$x, 0.4))
#' fit <- frm(bf(y ~ x) + lognormal(), data = dd)
#'
#' # predict() summarizes the predictive distribution; frm_linpred()
#' # is the linear predictor and fitted() the expected response
#' head(predict(fit, ndraws = 200))
#' head(frm_linpred(fit))
#' head(fitted(fit))
#'
#' # the last two are related by the family's own mean
#' head(exp(frm_linpred(fit) + sigma(fit)^2 / 2) -
#'        fitted(fit)[, "Estimate"])
#'
#' # a sigma nobody wrote a formula for is a distributional parameter
#' # on its own scale, not a coefficient on its link
#' summary(fit)$spec_pars
#' sigma(fit)
#' @return This page documents a convention. It is not a function, so it
#'   returns no value.
#' @name frmtmb-scales
NULL
