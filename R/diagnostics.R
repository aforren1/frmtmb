# Simulation-based residual diagnostics through DHARMa.

#' DHARMa residual diagnostics
#'
#' Builds a DHARMa object from the fit's family simulator: scaled
#' quantile residuals that are uniform under a correctly specified
#' model. Plot the result with `plot()`, or run
#' [DHARMa::testUniformity()], [DHARMa::testDispersion()],
#' [DHARMa::testZeroInflation()] and friends on it.
#'
#' @param fit A `frmtmb_fit` (univariate; the family needs a simulator).
#' @param nsim Number of simulated response vectors.
#' @param re.form Passed to [simulate.frmtmb_fit()]: `NULL` (default)
#'   conditions on the estimated random effects; `NA` redraws them.
#' @param seed Optional RNG seed for the simulations.
#' @param ... Passed to [DHARMa::createDHARMa()].
#' @return A `DHARMa` object.
#' @export
dharma_residuals <- function(fit, nsim = 250, re.form = NULL,
                             seed = NULL, ...) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop("dharma_residuals() needs the 'DHARMa' package", call. = FALSE)
  }
  rspec <- uni_resp(fit, "dharma_residuals()")
  if (identical(rspec$family$type, "ordinal")) {
    stop("dharma_residuals() has no ordinal support: DHARMa's rank ",
         "transform needs a response on a numeric scale, and an ordinal ",
         "response has only an order. Use residuals(type = \"osa\") ",
         "instead", call. = FALSE)
  }
  # DHARMa works in fitted-row space, so the na.exclude padding
  # simulate() adds has to come back off
  sims <- as.matrix(na_unpad(fit, simulate(fit, nsim = nsim, seed = seed,
                                           re.form = re.form)))
  DHARMa::createDHARMa(
    simulatedResponse = sims,
    observedResponse = fit$frame$y[[rspec$resp_name]],
    fittedPredictedResponse = as.vector(stats::na.omit(fitted(fit))),
    integerResponse = identical(rspec$family$type, "discrete"),
    ...
  )
}