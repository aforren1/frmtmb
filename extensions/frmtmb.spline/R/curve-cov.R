# The joint covariance of a grid prediction, taken from core's seam.
#
# A simultaneous band, a delta-method derivative and an implicit-function
# feature all need the same object: the covariance of the WHOLE grid
# prediction, Sigma = A V A', where A maps the coefficient vector to the
# curve and V is the joint covariance of the fixed and random
# coefficients. A smooth's wiggly part is a random-effect block even when
# the smooth is a population term, so V must carry the b block.
# `vcov(fit, full = TRUE)` does not carry it under either of its
# branches, and `predict(se.fit = TRUE)` forms Sigma internally and
# returns only sqrt(diag(Sigma)).
#
# Until frmtmb 0.52.0 this package rebuilt A by unit perturbation, one
# predict() call per contributing coefficient, and read V out of
# `fit$cache$Vjoint`, which was an internal. Core now exports both
# halves: `frm_lp_basis()` returns A, its coefficient positions, the
# covariance at those positions and the variance that is NOT coefficient
# uncertainty, and `frm_joint_cov()` returns V on its own. Everything
# below is a consumer of the first of those.
#
# The check survives the change and is still not decoration:
# `diag(A V A') + extra_var` must equal `predict(se.fit = TRUE)^2`
# wherever core has an independent route to that number. What the check
# now verifies is that this package reads the seam correctly, rather
# than that a reconstruction reproduced it.

#' Run `expr`, holding back any `ps()` knot-span warning frmtmb raises
#' inside it.
#'
#' The seam is a CLASS, `frmtmb_ps_span_warning`, not the text: the
#' three exported functions here evaluate the curve at three, five or
#' one grid position per point the user asked about, so core's warning
#' would arrive once per internal evaluation and count stencil rows
#' rather than grid rows. Held back and returned, it can be surfaced
#' once under the name of the function the user actually called, or
#' turned into a refusal where a warning is the wrong answer.
#'
#' @noRd
sp_catch_span <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    expr,
    frmtmb_ps_span_warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = val, span = unique(msgs))
}

#' Refuse a feature search whose bracket leaves a `ps()` knot span.
#'
#' A warning is the right answer for a band, which the reader can see
#' bending to the intercept, and the wrong one for a feature: a peak or
#' a crossing located past the span is a peak or a crossing OF THE
#' DECAYING PARTIAL SUM, the search reports it as a number with a
#' standard error, and nothing in the output says which curve it belongs
#' to. One template, called from both places the bracket is checked.
#'
#' @noRd
sp_span_stop <- function(span) {
  if (!length(span)) return(invisible(NULL))
  stop("frm_curve_feature(): the search bracket leaves a ps() term's ",
       "knot span, so a root located in it would be a root of the ",
       "decaying partial sum rather than of the fitted curve, and the ",
       "implicit-function standard error beside it would describe ",
       "neither. Narrow newdata to the span. ",
       paste(span, collapse = " "), call. = FALSE)
}

#' The `ps()` span messages for the grid the USER passed, rather than
#' for the stencil the caller widened it into.
#'
#' [frm_curve_deriv()] evaluates the design at `c(x - e, x, x + e)` and
#' [frm_curve_feature()]'s stationary-point scan at `c(x - e1, x + e1)`,
#' with `e` a millionth of the grid's range. Core is handed those points
#' and counts those rows, so a grid whose endpoint sits exactly on a
#' knot is outside the span by `e` and gets a warning, or in the feature
#' case a refusal, for a bracket the user chose entirely inside it.
#' Following the refusal's own advice ("Narrow newdata to the span")
#' reproduced the refusal.
#'
#' One `predict()` on the grid itself settles it. Callers only reach
#' here when the widened call already reported something, which for the
#' derivative stencil is exact (its point set CONTAINS the grid, so a
#' clean stencil implies a clean grid) and for the feature scan holds
#' whenever the spline argument is monotone in `var`, which is every
#' shape this package documents.
#'
#' @noRd
sp_span_on_grid <- function(fit, nd, dpar, resp, re.form) {
  sp_catch_span(sp_predict_eta(fit, nd, dpar, resp, re.form))$span
}

#' One prediction on the link scale, as a plain numeric vector.
#'
#' @noRd
sp_predict_eta <- function(fit, newdata, dpar, resp, re.form) {
  as.numeric(stats::predict(fit, newdata = newdata, type = "link",
                            dpar = dpar, resp = resp, re.form = re.form))
}

#' Everything the three exported functions share: the grid, the design,
#' the covariance, and the check that the covariance is the right one.
#'
#' `frm_lp_basis()` (frmtmb >= 0.52.0) returns `A`, `V` and the variance
#' that is not coefficient uncertainty. For a linear predictor it is the
#' same object `predict(se.fit = TRUE)` reduces to a diagonal, so the
#' two are compared on every call. For a NONLINEAR body `A` is a
#' Jacobian, `predict(se.fit = TRUE)` refuses outright, and there is
#' nothing to compare against: the check is skipped and `rel` is `NA`,
#' which `print()` reports rather than hides.
#'
#' @noRd
sp_curve_parts <- function(fit, newdata, dpar, resp, re.form, tol) {
  if (!inherits(fit, "frmtmb_fit")) {
    stop("frm_curve(): `object` must be a frmtmb fit, the model a curve ",
         "is read off, not an object of class ", class(fit)[1L],
         call. = FALSE)
  }
  if (!is.data.frame(newdata) || !nrow(newdata)) {
    stop("`newdata` must be a data frame with at least one row: it is ",
         "the grid the curve is evaluated on", call. = FALSE)
  }
  lbc <- sp_catch_span(frmtmb::frm_lp_basis(fit, newdata = newdata,
                                            dpar = dpar, resp = resp,
                                            re.form = re.form))
  lb <- lbc$value
  C <- as.matrix(lb$A)
  Sigma <- unname(C %*% lb$V %*% t(C))
  se <- unname(sqrt(pmax(diag(Sigma) + lb$extra_var, 0)))

  # A nonlinear body is the case core refuses se.fit for, so there is no
  # second number to check against. Everything else is checked.
  if (sp_is_nl(fit, dpar, resp)) {
    return(list(eta = lb$eta, C = C, V = lb$V, Sigma = Sigma, se = se,
                se_ref = rep(NA_real_, length(se)), rel = NA_real_,
                n_predict = 0L, newdata = newdata, dpar = dpar,
                resp = resp, re.form = re.form, fit = fit,
                span = lbc$span))
  }
  ref <- stats::predict(fit, newdata = newdata, type = "link", dpar = dpar,
                        resp = resp, re.form = re.form, se.fit = TRUE)
  se_ref <- as.numeric(ref$se.fit)
  rel <- max(abs(se / pmax(se_ref, .Machine$double.eps) - 1))
  if (!is.finite(rel) || rel > tol) {
    stop("frm_curve(): the assembled curve covariance disagrees with ",
         "predict(se.fit = TRUE) by ", format(rel, digits = 3),
         " relative, which is above the tolerance ", format(tol),
         ". Both come from frm_lp_basis(); a disagreement means this ",
         "package is reading the seam wrongly, and a fit where the two ",
         "disagree is one it must not report a band for", call. = FALSE)
  }
  list(eta = lb$eta, C = C, V = lb$V, Sigma = Sigma, se = se,
       se_ref = se_ref, rel = rel, n_predict = 1L,
       newdata = newdata, dpar = dpar, resp = resp, re.form = re.form,
       fit = fit, span = lbc$span)
}

#' Is this linear predictor computed by a nonlinear body?
#'
#' @noRd
sp_is_nl <- function(fit, dpar, resp) {
  rn <- resp %||% names(fit$spec$responses)[1L]
  rspec <- fit$spec$responses[[rn]]
  dp <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1L]
  lp <- fit$frame[["linpreds"]][[paste0(rn, ".", dp)]]
  !is.null(lp) && !is.null(lp[["nl_body"]])
}

#' The max-deviation simultaneous critical value of Ruppert, Wand and
#' Carroll (2003, ch. 6).
#'
#' Draw the curve's own deviation process, standardize it by `div`, take
#' the largest absolute value over the grid, and report the `level`
#' quantile. `div` is a separate argument rather than `sqrt(diag(S))`
#' because the critical value is comparable between two packages only
#' when the divisor is: gratia standardizes a smooth-only deviation by
#' the FULL predictor's standard error, and the critical value that
#' comes back is 8 percent smaller than the self-standardized one on the
#' same fit. The BAND is right either way, because the same divisor
#' calibrates it and scales it.
#'
#' @section A divisor of exactly zero:
#' A grid point whose standard error is exactly zero has a
#' DETERMINISTIC deviation, not a small one. `S` is positive
#' semi-definite, so a zero diagonal entry forces that whole row to zero,
#' and `(L z)` is exactly 0 there: the standardized deviation is `0 / 0`.
#' Such a point is covered with probability one, so it cannot be the
#' argmax and it leaves the maximization rather than turning it into
#' `NaN` and killing `quantile()` one line later.
#'
#' This is not a corner case. Past a [frmtmb::ps()] term's knot span the
#' basis is exactly zero, so its DERIVATIVE design is exactly zero and
#' every such row of a [frm_curve_deriv()] grid has a standard error of
#' exactly zero. Before this guard, `frm_curve_deriv()` on its own
#' defaults raised the span warning and then died in `quantile.default`
#' with a message naming neither the span nor the function.
#'
#' Dropping nothing is bit-identical to not having the guard, so a grid
#' with no zero divisor is unaffected.
#'
#' @noRd
sp_sim_crit <- function(S, div, nsim, level, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- nrow(S)
  keep <- is.finite(div) & div > 0
  if (!any(keep)) {
    stop("simultaneous = TRUE: every point on this grid has a standard ",
         "error of exactly zero, so the deviation process is degenerate ",
         "and there is no maximum to take a quantile of. Past a ps() ",
         "term's knot span the basis is exactly zero and so is the ",
         "derivative design, which is the usual way to arrive here. Use ",
         "simultaneous = FALSE, or move the grid inside the span",
         call. = FALSE)
  }
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = m)
  z <- matrix(stats::rnorm(m * nsim), m, nsim)
  mx <- apply(abs((L[keep, , drop = FALSE] %*% z) / div[keep]), 2L, max)
  crit <- unname(stats::quantile(mx, level, type = 8))
  # The standard error of a sample quantile, sqrt(p(1-p)/n) / f(q). A
  # critical value reported without it invites a comparison across
  # packages that simulation noise alone would fail.
  d <- stats::density(mx)
  f <- stats::approx(d$x, d$y, xout = crit)$y
  list(crit = crit,
       mcse = if (isTRUE(is.finite(f) && f > 0)) {
         sqrt(level * (1 - level) / nsim) / f
       } else NA_real_)
}
