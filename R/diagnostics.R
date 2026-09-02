# Simulation-based residual diagnostics through DHARMa.

#' DHARMa residual diagnostics
#'
#' Builds a DHARMa object from the fit's family simulator: scaled
#' quantile residuals that are uniform under a correctly specified
#' model. Plot the result with `plot()`, or run
#' [DHARMa::testUniformity()], [DHARMa::testDispersion()],
#' [DHARMa::testZeroInflation()] and friends on it.
#'
#' @section Ordinal responses:
#' An ordinal fit is supported: the draws are the simulated categories
#' and the response is its own integer codes, so the rank transform runs
#' on the order alone (`integerResponse = TRUE`, which is what makes
#' DHARMa randomize within ties). The `fittedPredictedResponse` DHARMa
#' plots against is the expected category index `sum_k k * P(y = k)`,
#' the same scalar [residuals()] scores an ordinal fit by; it sets the
#' horizontal axis of the display and nothing else, since the residuals
#' themselves come from the ranks of the draws.
#'
#' @param fit A `frmtmb_fit` (univariate; the family needs a simulator).
#' @param nsim Number of simulated response vectors.
#' @param re.form Passed to [simulate.frmtmb_fit()]: `NULL` (default)
#'   conditions on the estimated random effects; `NA` redraws them.
#' @param seed Optional RNG seed for the simulations.
#' @param ... Passed to [DHARMa::createDHARMa()].
#' @return A `DHARMa` object.
#' @examples
#' if (requireNamespace("DHARMa", quietly = TRUE)) {
#'   set.seed(1)
#'   dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#'   dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
#'   fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#'
#'   # scaled quantile residuals: uniform under a correct model, which
#'   # is what Pearson residuals cannot give for a discrete family
#'   res <- dharma_residuals(fit, nsim = 100, seed = 1)
#'   plot(res)
#'   DHARMa::testUniformity(res, plot = FALSE)
#'   DHARMa::testDispersion(res, plot = FALSE)
#' }
#' @export
dharma_residuals <- function(fit, nsim = 250, re.form = NULL,
                             seed = NULL, ...) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop("dharma_residuals() needs the 'DHARMa' package", call. = FALSE)
  }
  rspec <- uni_resp(fit, "dharma_residuals()")
  ordinal <- identical(rspec$family$type, "ordinal")
  # DHARMa works in fitted-row space, so the na.exclude padding
  # simulate() adds has to come back off
  sims <- na_unpad(fit, simulate(fit, nsim = nsim, seed = seed,
                                 re.form = re.form))
  sims <- if (ordinal) {
    # simulate() hands ordinal draws back as ordered factors; the rank
    # transform needs the integer codes the response itself carries
    matrix(unlist(lapply(sims, as.integer), use.names = FALSE),
           nrow = nrow(sims))
  } else {
    as.matrix(sims)
  }
  fpr <- if (ordinal) {
    # DHARMa wants ONE number per observation and an ordinal response
    # has no mean. The expected category index sum_k k * p_k is used:
    # it is the quantity residuals(type = "response") is taken against,
    # it is monotone in the latent predictor (so a trend against it
    # reads the same way as a trend against eta), and unlike the latent
    # predictor it lives on the scale of the observed response, which is
    # what DHARMa's residual-versus-fitted display assumes. The
    # calibration itself never sees it: the quantile residuals come from
    # the ranks of the simulated draws alone.
    ord_cat_moments(fit, rspec)$mean
  } else {
    as.vector(stats::na.omit(fitted(fit)))
  }
  DHARMa::createDHARMa(
    simulatedResponse = sims,
    observedResponse = fit$frame$y[[rspec$resp_name]],
    fittedPredictedResponse = fpr,
    integerResponse = identical(rspec$family$type, "discrete") || ordinal,
    ...
  )
}