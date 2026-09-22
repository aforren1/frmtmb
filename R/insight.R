# Registered methods for the insight package (the accessor layer under
# performance, parameters, and the rest of easystats). The default
# methods already work off our standard accessors; these fill the gaps
# the defaults cannot see: the random-effect formula split (which gates
# is_mixed_model and r2_nakagawa) and clean parameter extraction.

#' @exportS3Method insight::find_formula
find_formula.frmtmb_fit <- function(x, verbose = TRUE, ...) {
  rspec <- x$spec$responses[[1L]]
  dp <- rspec$dpars[["mu"]] %||% rspec$dpars[[1L]]
  out <- list(
    conditional = stats::as.formula(
      call("~", rspec$resp_expr, reformulas::RHSForm(dp[["fixed"]])),
      env = rspec$formula_env
    )
  )
  bars <- lapply(dp[["re"]] %||% list(), function(rt) {
    stats::as.formula(call("~", rt$bar), env = rspec$formula_env)
  })
  if (length(bars)) {
    out$random <- if (length(bars) == 1L) bars[[1L]] else unname(bars)
  }
  structure(out, class = c("insight_formula", "list"))
}

#' @exportS3Method insight::find_random
find_random.frmtmb_fit <- function(x, split_nested = FALSE,
                                   flatten = FALSE, ...) {
  bks <- Filter(function(bk) bk[["covstruct"]] != "smooth",
                x$frame[["re_blocks"]])
  grps <- unique(vapply(bks, `[[`, "", "group_name"))
  if (!length(grps)) return(NULL)
  if (flatten) grps else list(random = grps)
}

#' @exportS3Method insight::get_parameters
get_parameters.frmtmb_fit <- function(x, ...) {
  bd <- x$estimates[["betad"]]
  if (length(fx <- x$frame[["betad_fixed_idx"]])) bd <- bd[-fx]
  # every estimated coefficient plus an ordinal fit's thresholds and
  # cs() coefficients, the thresholds on the MODEL'S scale and under
  # fixef()'s names (`Intercept[1]`), not as the internal (tau_1, log
  # increments) vector marginaleffects perturbs; get_varcov() below is
  # their covariance, so the two line up row for row
  ts <- insight_param_scale(x)
  n_ex <- length(ts$est) - length(x$estimates[["beta"]]) - length(bd)
  data.frame(
    Parameter = names(ts$est),
    Estimate = unname(ts$est),
    Component = c(rep("conditional", length(x$estimates[["beta"]])),
                  rep("dispersion", length(bd)),
                  rep("conditional", n_ex))
  )
}

#' @exportS3Method insight::get_varcov
get_varcov.frmtmb_fit <- function(x, ...) {
  # the covariance of the parameters get_parameters() reports, which is
  # every estimated coefficient plus an ordinal fit's thresholds on the
  # model's scale; vcov() is brms's population-level block and leaves
  # out an intercept-only dispersion parameter that get_parameters()
  # does report
  insight_param_scale(x)$V
}

#' The parameter vector insight's seams report, and its covariance.
#'
#' `interop_coef_vector()` holds an ordinal fit's extras on the
#' INTERNAL scale, which is what marginaleffects needs: it perturbs the
#' vector and writes it back. insight only READS, and a reader expects
#' the thresholds fixef() reports. So the extras are mapped through the
#' family's threshold map, under fixef()'s names, and the covariance
#' follows by the delta method. Every other fit gets
#' `interop_coef_vector()` and `interop_vcov()` unchanged.
#'
#' @noRd
insight_param_scale <- function(x) {
  est <- interop_coef_vector(x)
  V <- interop_vcov(x)
  ex <- brms_extra_fixef(x)
  # a model whose thresholds brms_extra_fixef() cannot name (more than
  # one ordinal response shares one component) keeps the internal scale
  # rather than a guessed alignment
  if (!length(ex) ||
        !identical(vapply(ex, `[[`, "", "comp"), ord_extra_comps(x))) {
    return(list(est = est, V = V))
  }
  n0 <- length(estimated_coef_names(x))
  A <- diag(length(est))
  nm <- names(est)
  off <- n0
  for (e in ex) {
    k <- length(e$raw)
    pos <- off + seq_len(k)
    for (i in seq_len(k)) {
      A[pos[i], pos] <- fd_gradient_row(e$map, e$raw, i)
    }
    est[pos] <- e$values
    nm[pos] <- e$names
    off <- off + k
  }
  names(est) <- nm
  V <- A %*% V %*% t(A)
  dimnames(V) <- list(nm, nm)
  list(est = est, V = V)
}

# insight::get_residuals() has no method for this class either, so it
# fell through to its default, which returns stats::residuals(x)
# unchanged. residuals() is brms's n x 4 summary matrix since item
# 2.6f, so that default started handing every caller a matrix where
# insight's own contract is one number per observation: mean(), sd()
# and qqnorm() on the result all changed value in silence. insight has
# the same problem with brms and solves it in get_residuals.brmsfit,
# whose body this is: take the Estimate column, keep the whole matrix
# on a "full" attribute, and carry insight's own class.

#' @exportS3Method insight::get_residuals
get_residuals.frmtmb_fit <- function(x, ...) {
  r_attr <- stats::residuals(x, ...)
  r <- as.vector(r_attr[, "Estimate"])
  attr(r, "full") <- r_attr
  class(r) <- c("insight_residuals", class(r))
  r
}

# insight::get_predicted() has no method for this class, so it fell
# through to its default, which calls predict(x, type = ) and then
# fitted(x, type = ). Neither takes `type` any more (item 2.6d), and
# insight's default WARNS and returns NULL rather than erroring, which
# is a silent loss of every prediction. A method of its own says which
# quantity each of insight's `predict` values is here, and it can carry
# the standard errors that the default never had: the base build
# returned NA in every CI column.

#' @exportS3Method insight::get_predicted
get_predicted.frmtmb_fit <- function(x, data = NULL,
                                     predict = "expectation",
                                     ci = 0.95, verbose = TRUE, ...) {
  pr <- predict %||% "expectation"
  if (!pr %in% c("expectation", "response", "link", "prediction")) {
    frm_stop("insight::get_predicted(predict = \"", pr, "\") has no ",
             "counterpart here. \"expectation\" and \"response\" are ",
             "fitted(), \"link\" is the linear predictor, and ",
             "\"prediction\" is predict()'s predictive summary",
             call. = FALSE)
  }
  level <- if (is.null(ci) || !is.numeric(ci)) 0.95 else ci
  a <- (1 - level) / 2
  if (identical(pr, "prediction")) {
    sm <- stats::predict(x, newdata = data, probs = c(a, 1 - a))
    est <- sm[, "Estimate"]
    se <- sm[, "Est.Error"]
    lo <- sm[, 3L]
    hi <- sm[, 4L]
  } else {
    # a category-valued response has no single expectation, so the
    # latent linear predictor is what one number per row can be, which
    # is what insight's default produced here before
    type <- if (identical(pr, "link")) "link" else "response"
    fam <- x$spec$responses[[1L]]$family
    if (fam[["type"]] %in% c("ordinal", "categorical")) type <- "link"
    p <- tryCatch(frm_linpred(x, newdata = data, type = type,
                              se.fit = TRUE),
                  error = function(e) NULL)
    if (is.null(p)) {
      est <- as.numeric(frm_linpred(x, newdata = data, type = type))
      se <- rep(NA_real_, length(est))
    } else {
      est <- as.numeric(p$fit)
      se <- as.numeric(p$se.fit)
    }
    q <- stats::qnorm(1 - a)
    lo <- est - q * se
    hi <- est + q * se
  }
  out <- as.numeric(est)
  attr(out, "ci_data") <- data.frame(SE = as.numeric(se),
                                     CI_low = as.numeric(lo),
                                     CI_high = as.numeric(hi))
  attr(out, "data") <- data
  attr(out, "predict") <- pr
  class(out) <- c("get_predicted", "numeric")
  out
}

#' @exportS3Method insight::find_statistic
find_statistic.frmtmb_fit <- function(x, ...) {
  # No measured call site in insight or marginaleffects passes a
  # name this method lacks (dev/argspell-exempt2.R), so the
  # exemption was a hole rather than a contract.
  frm_check_dots(...)
  "z-statistic"
}

#' @exportS3Method insight::link_inverse
link_inverse.frmtmb_fit <- function(x, ...) {
  # No measured call site in insight or marginaleffects passes a
  # name this method lacks (dev/argspell-exempt2.R), so the
  # exemption was a hole rather than a contract.
  frm_check_dots(...)
  rspec <- x$spec$responses[[1L]]
  dp <- rspec$dpars[["mu"]] %||% rspec$dpars[[1L]]
  dp[["link"]]$linkinv
}

#' @exportS3Method insight::link_function
link_function.frmtmb_fit <- function(x, ...) {
  # No measured call site in insight or marginaleffects passes a
  # name this method lacks (dev/argspell-exempt2.R), so the
  # exemption was a hole rather than a contract.
  frm_check_dots(...)
  rspec <- x$spec$responses[[1L]]
  dp <- rspec$dpars[["mu"]] %||% rspec$dpars[[1L]]
  dp[["link"]]$linkfun
}
