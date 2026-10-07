# brms post-processing calls on a fit that a ported script meets:
# stancode(), standata(), pp_mixture() and add_criterion().
#
# Their generics are brms's (frm_generic_owners) and lived in
# frmtmb.sample until the methods for a fit were added here, as
# parnames() moved before them: a generic defined in the extension had
# no method for a fit, so `stancode(fit)` said "no applicable method",
# and with brms loaded it reached brms's default, which said "Data must
# be specified using the 'data' argument" about a fit that has data.

#' Stan code and Stan data of a model
#'
#' brms generates a Stan program and a Stan data list for every model,
#' and `stancode()` and `standata()` return them. frmtmb generates
#' neither: the model is an R closure that `build_objective()` builds
#' from the assembled frame, and RTMB differentiates it. So on a fit
#' both stop with that reason and name what holds the same content.
#'
#' `standata()` does not return the frame instead. The frame's names
#' are not brms's (`X`, `Z_1_1`, `J_1`, ...), and a list that answers
#' some of brms's names and not others would be read as brms's by a
#' ported script.
#'
#' @param object A `frmtmb_fit`.
#' @param ... Ignored; these methods always stop.
#' @return These methods never return; they signal an error.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(40))
#' dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
#' fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
#' try(stancode(fit))
#' try(standata(fit))
#'
#' # the closure that IS the model evaluates the joint negative log
#' # density; with no random effects that is the negative log likelihood
#' nll <- build_objective(fit$frame)
#' all.equal(nll(fit$estimates), -as.numeric(logLik(fit)))
#' @name stancode
#' @export
stancode <- function(object, ...) UseMethod("stancode")

#' @rdname stancode
#' @exportS3Method brms::stancode
#' @export
stancode.frmtmb_fit <- function(object, ...) {
  frm_stop("stancode() has no meaning for a frmtmb fit: there is no Stan ",
           "program. The model is an R closure built by build_objective() ",
           "from the assembled frame and differentiated by RTMB, and the ",
           "closure IS the source: build_objective(fit$frame) returns it, ",
           "fit$obj$fn is the taped marginal likelihood, and fit$frame ",
           "holds everything baked into it. See the 'Seeing the model ",
           "code' section of vignette(\"brms-migration\")", call. = FALSE)
}

#' @rdname stancode
#' @export
standata <- function(object, ...) UseMethod("standata")

#' @rdname stancode
#' @exportS3Method brms::standata
#' @export
standata.frmtmb_fit <- function(object, ...) {
  frm_stop("standata() has no meaning for a frmtmb fit: nothing is ",
           "exported to a Stan data list. The assembled frame fit$frame ",
           "holds the same content (the response, the design matrices, ",
           "the sparse Z, the addition terms) under frmtmb's names, and ",
           "model.matrix(), getME() and model.frame() read the pieces of ",
           "it individually", call. = FALSE)
}

#' Posterior mixture-component probabilities at the estimates
#'
#' brms's `pp_mixture()` on a maximum-likelihood fit. For a [mixture()]
#' or [mixture_mvn()] fit, including a mixture of ordinal families, it
#' gives the probability that each observation came from each
#' component, given its response, with the parameters at their
#' estimates. `Estimate` is [mixture_probs()]. The quantile columns are
#' a Wald interval on the logit scale, so they stay inside (0, 1): the
#' logit of the estimate plus a normal quantile times its delta-method
#' standard error, taken over the estimated parameters and the
#' group-level effects. `Est.Error` is the standard deviation of that
#' same law, the probability `plogis(Z)` with `Z` normal around the
#' logit, integrated numerically (100-point Gauss-Hermite, and adaptive
#' quadrature where the law is wide). brms reports
#' the mean, the SD and the quantiles of the same probability over the
#' posterior draws; `pp_mixture()` on `frmtmb.sample::frm_sample()`
#' draws does that here.
#'
#' Why not the delta-method standard error of the probability itself: a
#' probability near 0 or 1 is far from linear in the parameters, and
#' that first-order error understates the spread there by orders of
#' magnitude. On a two-component gaussian mixture of 300 rows, against
#' brms 2.23.0's `pp_mixture()` on the same data (4 chains of 1500
#' draws, brms's default priors), the delta-method error was a median
#' 0.36 of brms's posterior SD over all 300 rows, and 0.008 within 1e-4
#' of an edge. `Est.Error` as reported here was a median 1.11 of brms's
#' SD over all rows, and 0.84 to 1.27 (medians of bins by distance from
#' the edge) for probabilities between 0.001 and 0.999. Within 0.001 of
#' an edge it overstates: 2.4 and 4.8 times brms's SD (bin medians),
#' where brms's is a median 0.015 and 0.007. The interval's width was a
#' median 1.04 of brms's, and each of the 300 estimates was inside
#' brms's 95% interval.
#'
#' The array has brms's shape and names: observations by statistics by
#' components, with the components named `P(K = k | Y)`. For a
#' group-level mixture (`mixture(groups = )`) the rows are groups, as in
#' [mixture_probs()].
#'
#' @param x A `frmtmb_fit` with a mixture family.
#' @param newdata,re_formula brms's arguments, in brms's positions.
#'   `newdata` is refused: the probabilities classify the rows the model
#'   was fitted on. `re_formula = NULL` (the default) conditions on the
#'   group-level effects, which is what is computed; any other value is
#'   refused.
#' @param resp The response, for brms's signature. A multivariate fit is
#'   refused.
#' @param ndraws,draw_ids,summary,robust brms's draws arguments. A
#'   maximum-likelihood fit has no draws, so a value other than the
#'   default is refused by name with the reason.
#' @param log If `TRUE`, log probabilities: the estimate and the
#'   interval ends are logged, and `Est.Error` is the standard deviation
#'   of the log under the same law.
#' @param probs Probabilities of the quantile columns, as in brms.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it.
#' @return An `observations x statistics x components` array.
#' @seealso [mixture_probs()] for the bare matrix.
#' @examples
#' set.seed(4)
#' dd <- data.frame(y = c(rnorm(60, -2), rnorm(60, 3)))
#' fit <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()),
#'            data = dd)
#' pm <- pp_mixture(fit)
#' dim(pm)
#' head(pm[, , "P(K = 1 | Y)"])
#' all.equal(unname(pm[, "Estimate", ]), unname(mixture_probs(fit)))
#' @export
pp_mixture <- function(x, ...) UseMethod("pp_mixture")

#' @rdname pp_mixture
#' @exportS3Method brms::pp_mixture
#' @export
pp_mixture.frmtmb_fit <- function(x, newdata = NULL, re_formula = NULL,
                                  resp = NULL, ndraws = NULL,
                                  draw_ids = NULL, log = FALSE,
                                  summary = TRUE, robust = FALSE,
                                  probs = c(0.025, 0.975), ...) {
  frm_check_dots(...)
  fit_refuse_draws_args("pp_mixture()", summary = summary, robust = robust)
  check_flag(log, "log")
  require_fitted(x, "pp_mixture()")
  if (!is.null(ndraws) || !is.null(draw_ids)) {
    frm_stop("pp_mixture() cannot honor `",
             if (!is.null(ndraws)) "ndraws" else "draw_ids",
             "`: brms chooses posterior draws with it, and a ",
             "maximum-likelihood fit has none. Sample with ",
             "frmtmb.sample::frm_sample() and call pp_mixture() on the ",
             "draws for that", call. = FALSE)
  }
  if (!is.null(newdata)) {
    frm_stop("pp_mixture() does not take newdata on a fit: the component ",
             "probability of a row is a statement about that row's own ",
             "response, and it is computed for the rows the model was ",
             "fitted on", call. = FALSE)
  }
  if (!is.null(re_formula)) {
    frm_stop("pp_mixture() does not take re_formula: the probabilities ",
             "condition on the group-level effects at their modes, which ",
             "is brms's default, re_formula = NULL", call. = FALSE)
  }
  if (length(x$spec$responses) > 1L) {
    frm_stop("pp_mixture() is not supported for a multivariate fit yet: ",
             "mixture_probs() reads one response", call. = FALSE)
  }
  rn <- names(x$spec$responses)
  if (!is.null(resp) && !identical(resp, rn)) {
    frm_stop("pp_mixture(resp = \"", paste(resp, collapse = ", "),
             "\") names no response of this model; it has ", rn,
             call. = FALSE)
  }
  if (is.null(x$spec$responses[[1L]]$family[["mix"]])) {
    # brms's words
    frm_stop("Method 'pp_mixture' can only be applied to mixture models.",
             call. = FALSE)
  }
  qn <- brms_prob_cols(probs)
  P <- mixture_probs(x)
  # the group effects move a row's probability through its predictors,
  # so they join the differenced vector, as fitted()'s matrix route
  # takes them
  b_idx <- if (length(x$estimates[["b"]])) {
    ids <- sort(unique(c(smooth_b_idx(x), re_governed_b(x))))
    if (length(ids)) ids
  }
  K <- ncol(P)
  # each probability beside its complement, the sum of the other
  # components: near 1 a probability has lost the digits of 1 - p, and
  # its finite-difference error is rounding, while the complement keeps
  # both. So the logit and its standard error are read off whichever of
  # the two is the smaller
  pq <- function(fit) {
    Pf <- mixture_probs(fit)
    cbind(Pf, vapply(seq_len(K), function(k) {
      rowSums(Pf[, -k, drop = FALSE])
    }, numeric(nrow(Pf))))
  }
  Q <- pq(x)[, K + seq_len(K), drop = FALSE]
  se2 <- fit_fd_se(x, pq, b_idx = b_idx)
  se2 <- if (is.null(se2)) {
    matrix(NA_real_, nrow(P), 2L * K)
  } else {
    matrix(as.numeric(se2), nrow(P))
  }
  out <- array(NA_real_, c(nrow(P), 2L + length(qn), K),
               dimnames = list(rownames(P) %||% seq_len(nrow(P)),
                               c("Estimate", "Est.Error", qn),
                               paste0("P(K = ", seq_len(K), " | Y)")))
  for (k in seq_len(K)) {
    p <- P[, k]
    q <- Q[, k]
    se <- ifelse(p <= q, se2[, k], se2[, K + k])
    ends <- mixture_prob_band(p, q, se, probs)
    err <- logitnormal_sd(p, q, se / (p * q), log = log)
    # a probability that rounded to 0 or 1 has no logit; its delta
    # error is 0 there, and so is the spread of the law
    flat <- !is.finite(err) & is.finite(se)
    err[flat] <- 0
    if (log) {
      out[, "Estimate", k] <- base::log(p)
      for (j in seq_along(qn)) out[, qn[j], k] <- base::log(ends[[j]])
    } else {
      out[, "Estimate", k] <- p
      for (j in seq_along(qn)) out[, qn[j], k] <- ends[[j]]
    }
    out[, "Est.Error", k] <- err
  }
  out
}

#' Nodes and weights of the n-point Gauss-Hermite rule for an
#' expectation under the standard normal, by Golub-Welsch: the
#' eigenvalues of the probabilists' Hermite Jacobi matrix are the nodes
#' and the squared first components of its eigenvectors the weights.
#'
#' @noRd
gh_std_normal <- function(n = 40L) {
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L))
  J[cbind(seq_len(n - 1L), 2:n)] <- off
  J[cbind(2:n, seq_len(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1L, ]^2)
}

#' The standard deviation of `plogis(Z)`, or of `log(plogis(Z))`, for
#' `Z ~ N(log(p) - log(q), lse^2)` row by row, `q` the complement of
#' `p`: the law the logit Wald interval of `pp_mixture()` is built from.
#'
#' Integrated, not linearized: near 0 or 1 the probability is far from
#' linear in `Z`, which is why the delta-method error understated the
#' spread. Two passes (mean, then squared deviations) and, above one
#' half, the complement `plogis(-Z)`, so that a spread far below the
#' value is not lost to cancellation against it.
#'
#' @noRd
logitnormal_sd <- function(p, q, lse, log = FALSE) {
  n <- length(p)
  out <- rep(NA_real_, n)
  ok <- which(is.finite(p) & is.finite(q) & is.finite(lse) & p > 0 &
                q > 0)
  if (!length(ok)) return(out)
  hi <- p[ok] > q[ok]
  # the logit of the smaller of p and its complement q, so the values
  # integrated are small where the probability is near 1
  lg <- log(p[ok]) - log(q[ok])
  mu <- ifelse(hi, -lg, lg)
  s <- abs(lse[ok])
  # where p is the larger, z is the complement's logit and log(p) is
  # log(1 - plogis(z)), which plogis(-z, log.p = TRUE) keeps finite
  tr <- function(z, h) {
    if (!log) return(stats::plogis(z))
    stats::plogis(if (h) -z else z, log.p = TRUE)
  }
  # Gauss-Hermite converges slowly once the normal is wide against the
  # logistic's curvature: 100 nodes hold 4e-6 relative up to a logit SD
  # of 4 and lose 0.6% at 10 (dev/surface-p1-gh.R), so a wider law is
  # integrated adaptively, row by row
  wide <- s > 4
  if (any(!wide)) {
    g <- gh_std_normal(100L)
    i <- which(!wide)
    Z <- outer(mu[i], rep(1, length(g$x))) + outer(s[i], g$x)
    V <- if (log) {
      # the test has Z's shape, or ifelse() returns hi's length
      ifelse(matrix(hi[i], nrow(Z), ncol(Z)),
             stats::plogis(-Z, log.p = TRUE), stats::plogis(Z, log.p = TRUE))
    } else {
      stats::plogis(Z)
    }
    dim(V) <- dim(Z)
    m <- as.vector(V %*% g$w)
    out[ok[i]] <- sqrt(as.vector(((V - m)^2) %*% g$w))
  }
  for (i in which(wide)) {
    dn <- function(z) stats::dnorm(z, mu[i], s[i])
    lim <- mu[i] + c(-14, 14) * s[i]
    int <- function(f) {
      tryCatch(
        stats::integrate(f, lim[1L], lim[2L], rel.tol = 1e-10,
                         subdivisions = 1000L)$value,
        error = function(e) {
          logitnormal_pieces(f, mu[i], s[i])
        })
    }
    m <- int(function(z) tr(z, hi[i]) * dn(z))
    out[ok[i]] <- sqrt(int(function(z) (tr(z, hi[i]) - m)^2 * dn(z)))
  }
  out
}

#' An integral over the logit `z` of a law `N(mu, s^2)`, in pieces.
#'
#' The fallback of `logitnormal_sd()`'s adaptive integral. At logit SDs
#' of thousands one `integrate()` over `mu +- 14 s` can stop with "the
#' integral is probably divergent": the logistic's whole transition is
#' a sliver of the range (the re-check of 2026-10-07, n1, at p 1e-3 and
#' 1e-5 with SD 3000 and p 1e-30 with SD 1e4). Breaking the range at
#' the normal's centre, a few SDs either side, and at the transition
#' `z = 0` with a few logistic widths either side keeps each piece
#' smooth, and no piece may stop the call: a piece that still fails
#' reports its estimate instead of an error.
#'
#' @noRd
logitnormal_pieces <- function(f, mu, s) {
  lo <- mu - 14 * s
  hi <- mu + 14 * s
  br <- c(mu + c(-14, -6, -3, -1, 0, 1, 3, 6, 14) * s,
          c(-40, -10, -3, 0, 3, 10, 40))
  br <- sort(unique(br[br >= lo & br <= hi]))
  sum(vapply(seq_len(length(br) - 1L), function(j) {
    stats::integrate(f, br[j], br[j + 1L], rel.tol = 1e-10,
                     subdivisions = 1000L, stop.on.error = FALSE)$value
  }, 0))
}

#' Wald interval ends of a probability on the logit scale, one vector
#' per entry of `probs`.
#'
#' A symmetric band around 0.99 would leave the unit interval, so the
#' band is built on the logit, `log(p) - log(q)` with `q` the
#' complement, where the standard error is `se / (p q)`. A probability
#' that rounded to 0 or 1 has no logit and a standard error of 0 there,
#' so its ends are the point.
#'
#' @noRd
mixture_prob_band <- function(p, q, se, probs) {
  lse <- se / (p * q)
  flat <- !is.finite(lse) & is.finite(se)
  lg <- log(p) - log(q)
  lapply(probs %||% numeric(0), function(pr) {
    v <- stats::plogis(lg + stats::qnorm(pr) * lse)
    v[flat] <- p[flat]
    v
  })
}

#' Add a model comparison criterion to a fit
#'
#' brms's `add_criterion()` computes `loo()`, `waic()`, `kfold()`,
#' `loo_subsample()`, `bayes_R2()`, `loo_R2()` or the marginal
#' likelihood once and stores it in the fit, so that a later `loo()` or
#' `loo_compare()` reads it instead of computing it again. Every one of
#' those criteria is a posterior quantity: an elpd, a pointwise
#' predictive density or an R2 averaged over draws, or a marginal
#' likelihood integrated over the prior. A maximum-likelihood fit has
#' none of them, and frmtmb refuses `loo()`, `waic()` and `bayes_R2()`
#' on a fit for that reason. So `add_criterion()` on a fit refuses each
#' criterion by name with the same reason. `AIC()` and `BIC()` are the
#' maximum-likelihood comparison, and they need nothing stored: they
#' read the log likelihood the fit already holds. On
#' `frmtmb.sample::frm_sample()` draws, `add_criterion()` stores the
#' criteria as brms does.
#'
#' @param x A `frmtmb_fit`.
#' @param criterion brms's criterion names: a subset of `"loo"`,
#'   `"waic"`, `"kfold"`, `"loo_subsample"`, `"bayes_R2"`, `"loo_R2"`
#'   and `"marglik"`.
#' @param model_name,overwrite,file,force_save brms's other arguments,
#'   accepted so that a ported call reaches the refusal. This method
#'   always stops.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it.
#' @return This method never returns; it signals an error.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(40))
#' dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
#' fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
#' try(add_criterion(fit, "loo"))
#'
#' # the maximum-likelihood comparison
#' fit0 <- frm(bf(y ~ 1), family = gaussian(), data = dd)
#' AIC(fit, fit0)
#' @export
add_criterion <- function(x, ...) UseMethod("add_criterion")

#' brms's criterion names, in brms's order, and why each needs draws.
#' Read by the fit method here and by frmtmb.sample's draws method,
#' which validates against the same names.
#'
#' @noRd
add_criterion_reasons <- c(
  loo = paste("an elpd averages the pointwise likelihood over posterior",
              "draws"),
  waic = paste("WAIC averages the pointwise likelihood over posterior",
               "draws"),
  kfold = paste("brms's kfold() refits the model and averages the",
                "held-out predictive density over posterior draws"),
  loo_subsample = paste("an elpd averages the pointwise likelihood over",
                        "posterior draws"),
  bayes_R2 = "Bayesian R2 is computed per posterior draw",
  loo_R2 = "LOO R2 weights posterior draws by their PSIS-LOO weights",
  marglik = paste("the marginal likelihood integrates the likelihood over",
                  "the prior, which bridge sampling estimates from draws")
)

#' brms's validation of `criterion`, with its deprecated `"R2"` alias.
#'
#' @noRd
add_criterion_names <- function(criterion) {
  ok <- !missing(criterion) && is.character(criterion) &&
    length(criterion) && !anyNA(criterion)
  if (ok && any(criterion == "R2")) {
    # brms's warning and its mapping
    frm_warning("Criterion 'R2' is deprecated. Please use 'bayes_R2' ",
                "instead.", call. = FALSE)
    criterion[criterion == "R2"] <- "bayes_R2"
  }
  criterion <- if (ok) unique(criterion)
  if (!ok || !all(criterion %in% names(add_criterion_reasons))) {
    # brms's words
    frm_stop("Argument 'criterion' should be a subset of ",
             paste(names(add_criterion_reasons), collapse = ", "),
             call. = FALSE)
  }
  criterion
}

#' @rdname add_criterion
#' @exportS3Method brms::add_criterion
#' @export
add_criterion.frmtmb_fit <- function(x, criterion, model_name = NULL,
                                     overwrite = FALSE, file = NULL,
                                     force_save = FALSE, ...) {
  frm_check_dots(...)
  criterion <- add_criterion_names(criterion)
  c1 <- criterion[1L]
  frm_stop("add_criterion(criterion = \"", c1, "\") has no meaning on a ",
           "maximum-likelihood fit: ", add_criterion_reasons[[c1]],
           ", and a fit has one estimate and no draws. Compare ",
           "maximum-likelihood fits with AIC() or BIC(), which read the ",
           "log likelihood the fit already holds, or with anova() for ",
           "nested ones. For the criterion itself, sample first: ",
           "add_criterion(frmtmb.sample::frm_sample(fit), \"", c1, "\")",
           call. = FALSE)
}
