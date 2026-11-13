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
#' A NOMINAL response (the `categorical()` family) is refused. A scaled
#' quantile residual is the predictive CDF evaluated at the observation,
#' and a CDF needs an ordered support. Ordinal categories have one, so
#' their codes carry real information; nominal ones do not, and their
#' 1..K codes are an arbitrary labeling - relabel the levels and every
#' residual moves, which is not a diagnostic. DHARMa says the same of
#' multinomial responses in its own vignette. Check a nominal fit with
#' [pp_check()] instead (`type = "bars"` compares the observed category
#' counts with the simulated ones), or read the per-category
#' probabilities through `conditional_effects()`.
#'
#' A composition (the [dirichlet()] and [logistic_normal()] families) is
#' refused for the same reason: a row of proportions has no scalar
#' distribution function to evaluate.
#'
#' @section Censored responses:
#' A `cens()` row observes an event, not a value: `Y > y` when right
#' censored, `Y <= y` when left censored, `y < Y <= y2` for an interval.
#' Its residual is a uniform draw between the predictive distribution
#' function at the two ends of that event, the randomized quantile
#' residual of Dunn and Smyth (1996) extended to censoring, and an
#' uncensored row keeps `F(y)`. To get it, the draws of a censored row
#' that fall inside its event are set to one point of the event, and
#' the observation to the same point; DHARMa's default PIT residual
#' (`method = "PIT"`) then randomizes over the tie, as it does for a
#' count. Under a correct model the residuals are uniform whatever each
#' row's censoring time is, so the censoring may vary by row and needs
#' no model of its own. On a censored row, `observedResponse` and the
#' row of `simulatedResponse` hold that point, not the recorded value
#' and the latent draws. On a fit with a censored row,
#' `method = "traditional"` and `rotation`, passed on to
#' [DHARMa::createDHARMa()], are refused: the first adds noise that moves
#' a censored row off its event, and the second mixes rows before the
#' event is read, so neither gives residuals that are uniform under a
#' correct model.
#'
#' @param fit A `frmtmb_fit` (univariate; the family needs a simulator).
#' @param nsim Number of simulated response vectors.
#' @param re_formula Passed to [simulate.frmtmb_fit()]: `NULL` (default)
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
dharma_residuals <- function(fit, nsim = 250, re_formula = NULL,
                             seed = NULL, ...) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    frm_stop("dharma_residuals() needs the 'DHARMa' package", call. = FALSE)
  }
  rspec <- single_response(fit, "dharma_residuals()")
  # A quantile residual is the predictive CDF at the observation, so it
  # needs an ordered support. Nominal category codes are an arbitrary
  # labeling: relabel the levels and every residual changes. Refusing is
  # the honest answer, not a rank transform over the level order.
  if (identical(rspec$family[["type"]], "categorical")) {
    frm_stop("dharma_residuals() has no meaning for a nominal response: a ",
             "scaled quantile residual is the predictive distribution ",
             "function evaluated at the observation, and an unordered ",
             "category has no distribution function - relabeling the levels ",
             "would change every residual. Use pp_check(type = \"bars\") to ",
             "compare observed and simulated category counts",
             call. = FALSE)
  }
  # A row of proportions has no scalar predictive CDF either: the draws
  # are compositions, and a rank of one row among them is not defined
  if (fam_is_simplex(rspec$family)) {
    frm_stop("dharma_residuals() has no meaning for the '",
             rspec$family[["family"]], "' family: a scaled quantile ",
             "residual is the predictive distribution function evaluated ",
             "at the observation, and a row of proportions has no such ",
             "function. Compare predict(fit), a summary of simulated ",
             "compositions, with the observed proportions instead",
             call. = FALSE)
  }
  # a row of counts over categories has no scalar CDF either; 0.71.0
  # stopped on the missing mean first, and once the family had one the
  # rank code failed on the matrix with a subscript error
  if (fam_is_multinomial(rspec$family)) {
    frm_stop("dharma_residuals() has no meaning for the '",
             rspec$family[["family"]], "' ",
             "family: a scaled quantile residual is the predictive ",
             "distribution function evaluated at the observation, and a ",
             "row of counts over categories has no such function. Compare ",
             "predict(fit), a summary of simulated counts, or fitted(fit), ",
             "the expected counts, with the observed counts instead",
             call. = FALSE)
  }
  ordinal <-identical(rspec$family[["type"]], "ordinal")
  # DHARMa works in fitted-row space, so the na.exclude padding
  # simulate() adds has to come back off
  sims <- na_unpad(fit, simulate(fit, nsim = nsim, seed = seed,
                                 re_formula = re_formula))
  sims <- if (ordinal) {
    # simulate() hands ordinal draws back as ordered factors; the rank
    # transform needs the integer codes the response itself carries,
    # which start at 0 for a hurdle family
    shift <- ord_code0(rspec$family) - 1L
    matrix(unlist(lapply(sims, function(s) {
      if (is.factor(s)) as.integer(s) + shift else as.integer(s)
    }), use.names = FALSE), nrow = nrow(sims))
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
    as.vector(stats::na.omit(fitted_point(fit)))
  }
  # a missing mi() response is no observation: the frame's placeholder
  # would be scored against the draws as if it had been seen
  yobs <- frame_observed_y(fit, rspec$resp_name)
  cen <- fit$frame[["aterm_values"]][[rspec$resp_name]][["cens"]]
  if (!is.null(cen) && any(cen != 0)) {
    dharma_censor_args(list(...))
    cc <- dharma_censor(sims, yobs, cen,
                        fit$frame[["aterm_values"]][[rspec$resp_name]][[
                          "cens_y2"]],
                        identical(rspec$family[["type"]], "discrete"))
    sims <- cc$sims
    yobs <- cc$y
  }
  keep <- !is.na(yobs)
  DHARMa::createDHARMa(
    simulatedResponse = sims[keep, , drop = FALSE],
    observedResponse = yobs[keep],
    fittedPredictedResponse = fpr[keep],
    integerResponse = identical(rspec$family[["type"]], "discrete") || ordinal,
    ...
  )
}

#' Refuse the createDHARMa() arguments under which the censored residual
#' is not the one dharma_censor() sets up. "traditional" adds uniform
#' noise to the observation and the draws, which spreads a censored
#' row's tie onto values outside its event, and `rotation` mixes rows
#' before the tie is read. On a correct right-censored gaussian fit (94
#' of 300 rows censored) they gave KS p-values of 2.3e-25 and 2.3e-14,
#' where the default gives 0.86 (dev/silent-findings.md, "Fix round").
#'
#' @noRd
dharma_censor_args <- function(dots) {
  m <- dots[["method"]]
  if (is.character(m) && length(m) >= 1L &&
        identical(pmatch(m[1L], c("PIT", "traditional")), 2L)) {
    frm_stop("dharma_residuals(method = \"traditional\") is refused on a ",
             "cens() fit with a censored row: it adds noise to the ",
             "observation and the draws, which moves a censored row off ",
             "the event it observes, so the residuals are not uniform ",
             "under a correct model. Use the default method = \"PIT\"",
             call. = FALSE)
  }
  if (!is.null(dots[["rotation"]])) {
    frm_stop("dharma_residuals(rotation = ) is refused on a cens() fit ",
             "with a censored row: the rotation mixes rows before each ",
             "row's censored event is read, so the residuals are not ",
             "uniform under a correct model", call. = FALSE)
  }
  invisible(NULL)
}

#' Put a censored row's draws on the event its observation is, so that
#' DHARMa's PIT residual is the randomized quantile residual of a
#' censored observation.
#'
#' A row observes the event `a < Y <= b`: `(y, Inf)` when right
#' censored, `(-Inf, y]` when left censored, `(y, y2]` for an interval.
#' Each of the row's draws inside the event is set to one point `s` of
#' it, and the observation to `s` as well. DHARMa's PIT residual is then
#' uniform on `(mean(draws <= a), mean(draws <= b))`, the empirical
#' `(F(a), F(b))`, which is the residual of Dunn and Smyth (1996)
#' extended to censoring. It is uniform under the model whatever the
#' row's censoring time, so the censoring needs no model and may vary by
#' row. An uncensored row keeps its draws: its residual is `F(y)`.
#' Through 0.81.0 a censored row was scored as the exact value `y`, so a
#' right-censored row read `F(y)` instead of a draw from `(F(y), 1)`,
#' and nothing was said (dev/silent-findings.md, item 650).
#'
#' `s` must sit above every draw at or below `a` and below every draw
#' above `b`. On the integer lattice `a < Y` is `Y >= a + 1`, so a
#' right-censored count takes `s = a + 1`; a continuous one takes a
#' double just above `a`, so that a draw AT `a` (an atom of the family
#' there) still counts below. Left censoring and an interval take
#' `s = b`. Draws come from simulate(), so they already lie inside a
#' `trunc()` window, and the event is read inside it.
#'
#' @noRd
dharma_censor <- function(sims, y, cen, y2, discrete) {
  n <- length(y)
  for (i in which(cen != 0)) {
    d <- sims[i, ]
    if (cen[i] == 1) {
      s <- if (discrete) {
        y[i] + 1
      } else {
        y[i] + 4 * .Machine$double.eps * max(abs(y[i]), 1)
      }
      d[d > y[i]] <- s
    } else if (cen[i] == -1) {
      s <- y[i]
      d[d <= y[i]] <- s
    } else {
      s <- y2[i]
      d[d > y[i] & d <= y2[i]] <- s
    }
    sims[i, ] <- d
    y[i] <- s
  }
  stopifnot(length(y) == n)
  list(sims = sims, y = y)
}
