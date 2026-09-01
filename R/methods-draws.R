# Method surface for frm_sample() output: every draw is a complete
# parameter vector (tmbstan samples the random effects too), so the
# fitted-model machinery runs per draw exactly.

# Column positions of each template component inside the draws matrix
# (which is in template order, mapped betad entries absent, lp__ last).
draws_par_index <- function(fit) {
  tpl <- fit$frame$par_template
  idx <- list()
  pos <- 0L
  for (cp in names(tpl)) {
    len <- length(tpl[[cp]])
    if (cp == "betad" && length(fx <- fit$frame$betad_fixed_idx)) {
      len <- len - length(fx)
    }
    idx[[cp]] <- pos + seq_len(len)
    pos <- pos + len
  }
  idx
}

# The originating fit with its estimates replaced by one draw.
draws_fit_at <- function(x, i, idx = draws_par_index(x$fit)) {
  fit <- x$fit
  est <- fit$frame$par_template   # mapped betad entries keep link(const)
  row <- x$draws[i, ]
  for (cp in names(idx)) {
    if (cp == "betad" && length(fx <- fit$frame$betad_fixed_idx)) {
      pos <- setdiff(seq_along(est[[cp]]), fx)
      est[[cp]][pos] <- row[idx[[cp]]]
    } else {
      est[[cp]][] <- row[idx[[cp]]]
    }
  }
  for (cp in names(fit$estimates)) {
    names(est[[cp]]) <- names(fit$estimates[[cp]])
  }
  fit$estimates <- est
  fit$cache <- new.env(parent = emptyenv())   # no stale sdreport
  fit
}

draws_subsample <- function(x, ndraws) {
  n <- nrow(x$draws)
  if (is.null(ndraws) || ndraws >= n) return(seq_len(n))
  round(seq(1, n, length.out = ndraws))
}

#' @export
summary.frmtmb_draws <- function(object, ...) {
  m <- object$draws
  keep <- setdiff(colnames(m),
                  c("lp__", grep("^b\\[", colnames(m), value = TRUE)))
  tab <- t(vapply(keep, function(nm) {
    c(mean = mean(m[, nm]), sd = stats::sd(m[, nm]),
      `2.5%` = unname(stats::quantile(m[, nm], 0.025)),
      `97.5%` = unname(stats::quantile(m[, nm], 0.975)))
  }, numeric(4)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    ss <- rstan::summary(object$stanfit)$summary
    n_keep <- min(nrow(ss), ncol(m))
    conv <- ss[seq_len(n_keep), c("n_eff", "Rhat"), drop = FALSE]
    rownames(conv) <- colnames(m)[seq_len(n_keep)]
    ok <- intersect(keep, rownames(conv))
    tab <- cbind(tab, n_eff = NA_real_, Rhat = NA_real_)
    tab[ok, c("n_eff", "Rhat")] <- conv[ok, ]
  }
  tab
}

#' @export
fixef.frmtmb_draws <- function(object, ...) {
  nm <- gsub("[()]", "", estimated_coef_names(object$fit))
  nm0 <- estimated_coef_names(object$fit)
  idx <- draws_par_index(object$fit)
  cols <- c(idx$beta, idx$betad)
  m <- object$draws[, cols, drop = FALSE]
  out <- cbind(
    Estimate = colMeans(m),
    Est.Error = apply(m, 2, stats::sd),
    Q2.5 = apply(m, 2, stats::quantile, 0.025),
    Q97.5 = apply(m, 2, stats::quantile, 0.975)
  )
  rownames(out) <- nm0
  out
}

#' @export
VarCorr.frmtmb_draws <- function(x, ...) {
  fit <- x$fit
  idx <- draws_par_index(fit)
  if (is.null(idx$theta)) return(NULL)
  th_draws <- x$draws[, idx$theta, drop = FALSE]
  per_draw <- lapply(seq_len(nrow(th_draws)), function(i) {
    fit$estimates$theta <- th_draws[i, ]
    as.data.frame(VarCorr.frmtmb_fit(fit))$sdcor
  })
  base <- as.data.frame(VarCorr(fit))
  M <- do.call(rbind, per_draw)
  base$estimate <- colMeans(M)
  base$lwr <- apply(M, 2, stats::quantile, 0.025)
  base$upr <- apply(M, 2, stats::quantile, 0.975)
  base$vcov <- NULL
  base$sdcor <- NULL
  base
}

#' @rdname ranef
#' @export
ranef.frmtmb_draws <- function(object, ...) {
  fit <- object$fit
  if (!length(fit$frame$re_blocks)) return(list())
  idx <- draws_par_index(fit)
  n <- nrow(object$draws)
  per <- vector("list", n)
  for (i in seq_len(n)) {
    per[[i]] <- ranef(draws_fit_at(object, i, idx))
  }
  # brms shape: per term, a levels x statistics x coefficients array
  out <- list()
  for (tn in names(per[[1]])) {
    M0 <- per[[1]][[tn]]
    A <- vapply(per, function(r) r[[tn]], M0)
    st <- array(NA_real_, c(nrow(M0), 4L, ncol(M0)),
                dimnames = list(rownames(M0),
                                c("Estimate", "Est.Error",
                                  "Q2.5", "Q97.5"),
                                colnames(M0)))
    st[, "Estimate", ] <- apply(A, c(1, 2), mean)
    st[, "Est.Error", ] <- apply(A, c(1, 2), stats::sd)
    st[, "Q2.5", ] <- apply(A, c(1, 2), stats::quantile, 0.025)
    st[, "Q97.5", ] <- apply(A, c(1, 2), stats::quantile, 0.975)
    out[[tn]] <- st
  }
  out
}

#' @rdname hypothesis
#' @export
hypothesis.frmtmb_draws <- function(x, hypothesis, alpha = 0.05, ...) {
  fit <- x$fit
  exs <- lapply(hypothesis, hyp_parse)
  idx <- draws_par_index(fit)
  n <- nrow(x$draws)
  draws <- matrix(NA_real_, n, length(exs),
                  dimnames = list(NULL, hypothesis))
  for (i in seq_len(n)) {
    sh <- draws_fit_at(x, i, idx)
    w <- hyp_vals_only(sh)
    draws[i, ] <- vapply(exs, function(ex) {
      hyp_eval(sh, ex, w$vals, w$comp)
    }, numeric(1))
  }
  rows <- lapply(seq_along(exs), function(k) {
    t_k <- draws[, k]
    p <- min(1, 2 * min((1 + sum(t_k <= 0)) / (1 + n),
                        (1 + sum(t_k >= 0)) / (1 + n)))
    data.frame(hypothesis = hypothesis[k], estimate = mean(t_k),
               se = stats::sd(t_k),
               lwr = unname(stats::quantile(t_k, alpha / 2)),
               upr = unname(stats::quantile(t_k, 1 - alpha / 2)),
               z = mean(t_k) / stats::sd(t_k), p = p)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "method") <- "posterior"
  attr(out, "alpha") <- alpha
  attr(out, "draws") <- draws
  attr(out, "nsim") <- n
  attr(out, "converged") <- rep(TRUE, n)
  class(out) <- c("frmtmb_hypothesis", "data.frame")
  out
}

#' Expected-value and predictive draws from sampled parameters
#'
#' `posterior_epred()` evaluates the response-scale expectation per
#' draw; `posterior_predict()` additionally simulates responses from
#' the family, giving the posterior predictive distribution. Both
#' condition on each draw's own random effects (`re.form = NA` drops
#' them).
#'
#' @param object A `frmtmb_draws` from [frm_sample()].
#' @param newdata,resp,re.form As in [predict.frmtmb_fit()].
#' @param ndraws Number of draws to use (default: all).
#' @param ... Unused.
#' @return A draws-by-observations matrix.
#' @export
posterior_epred <- function(object, ...) UseMethod("posterior_epred")

#' @rdname posterior_epred
#' @export
posterior_epred.frmtmb_draws <- function(object, newdata = NULL,
                                         resp = NULL, re.form = NULL,
                                         ndraws = NULL, ...) {
  idx <- draws_par_index(object$fit)
  rows <- draws_subsample(object, ndraws)
  out <- NULL
  for (k in seq_along(rows)) {
    sh <- draws_fit_at(object, rows[k], idx)
    p <- predict(sh, newdata = newdata, resp = resp, re.form = re.form,
                 type = "response")
    if (is.null(out)) out <- matrix(NA_real_, length(rows), length(p))
    out[k, ] <- p
  }
  out
}

#' @rdname posterior_epred
#' @export
posterior_linpred <- function(object, ...) UseMethod("posterior_linpred")

#' @rdname posterior_epred
#' @param transform For `posterior_linpred()`: if `TRUE`, apply the
#'   inverse link (the value of the `mu` dpar on its natural scale,
#'   brms's convention; unlike `posterior_epred()` this is not the
#'   response mean for zero-inflated and similar families).
#' @param dpar For `posterior_linpred()`: which distributional
#'   parameter's linear predictor to evaluate.
#' @export
posterior_linpred.frmtmb_draws <- function(object, transform = FALSE,
                                           newdata = NULL, resp = NULL,
                                           re.form = NULL, dpar = NULL,
                                           ndraws = NULL, ...) {
  idx <- draws_par_index(object$fit)
  rows <- draws_subsample(object, ndraws)
  out <- NULL
  for (k in seq_along(rows)) {
    sh <- draws_fit_at(object, rows[k], idx)
    p <- predict(sh, newdata = newdata, resp = resp, dpar = dpar,
                 re.form = re.form,
                 type = if (transform) "conditional" else "link")
    if (is.null(out)) out <- matrix(NA_real_, length(rows), length(p))
    out[k, ] <- p
  }
  out
}

#' @rdname posterior_epred
#' @export
posterior_predict <- function(object, ...) UseMethod("posterior_predict")

#' @rdname posterior_epred
#' @export
posterior_predict.frmtmb_draws <- function(object, newdata = NULL,
                                           resp = NULL, re.form = NULL,
                                           ndraws = NULL, ...) {
  fit <- object$fit
  resp <- resp %||% names(fit$spec$responses)[1L]
  rspec <- fit$spec$responses[[resp]]
  if (is.null(rspec$family$sim)) {
    stop("Family '", rspec$family$family, "' has no simulator yet",
         call. = FALSE)
  }
  idx <- draws_par_index(object$fit)
  rows <- draws_subsample(object, ndraws)
  av <- if (is.null(newdata)) {
    fit$frame$aterm_values[[resp]]
  } else if (has_trunc(rspec)) {
    # truncation bounds must follow the newdata rows, or the draws land
    # outside the support the likelihood was normalized on
    aterms_for_newdata(rspec, newdata)
  } else {
    list()
  }
  out <- NULL
  for (k in seq_along(rows)) {
    sh <- draws_fit_at(object, rows[k], idx)
    dp <- if (is.null(newdata)) {
      eval_dpars(sh)[[resp]]
    } else {
      dpv <- list()
      for (dnm in names(rspec$dpars)) {
        dpv[[dnm]] <- as.vector(predict(sh, newdata = newdata,
                                        dpar = dnm, resp = resp,
                                        re.form = re.form,
                                        type = "response"))
      }
      dpv
    }
    ys <- sim_response(rspec$family, dp, av, length(dp[[1]]),
                       extra = fit_extras(sh))
    if (is.null(out)) out <- matrix(NA_real_, length(rows), length(ys))
    out[k, ] <- ys
  }
  out
}

#' @rdname pp_check
#' @exportS3Method bayesplot::pp_check
#' @export
pp_check.frmtmb_draws <- function(object, type = "dens_overlay",
                                  ndraws = 50, ...) {
  fit <- object$fit
  rspec <- uni_resp(fit, "pp_check()")
  y <- fit$frame$y[[1L]]
  if (is.matrix(y)) {
    stop("pp_check() supports vector responses", call. = FALSE)
  }
  yrep <- posterior_predict(object, ndraws = ndraws)
  fun <- get(paste0("ppc_", type), envir = asNamespace("bayesplot"))
  fun(as.numeric(y), yrep, ...)
}

#' Convert draws to a posterior draws object
#'
#' @param x A `frmtmb_draws` object.
#' @param ... Unused.
#' @return A `posterior::draws_matrix`: one column per sampled variable
#'   and one row per draw.
#' @export
as_draws <- function(x, ...) {
  # own generic for the same reason as pp_check: posterior stays in
  # Suggests, and as_draws(ds) must work without attaching it
  UseMethod("as_draws")
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws
#' @export
as_draws.frmtmb_draws <- function(x, ...) {
  posterior::as_draws_matrix(x$draws)
}

#' @rdname variables
#' @export
variables.frmtmb_draws <- function(x, ...) {
  colnames(x$draws)
}

#' @export
as.data.frame.frmtmb_draws <- function(x, ...) {
  as.data.frame(x$draws)
}