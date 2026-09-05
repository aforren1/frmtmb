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
  lb <- frmtmb::frm_lp_basis(fit, newdata = newdata, dpar = dpar,
                             resp = resp, re.form = re.form)
  C <- as.matrix(lb$A)
  Sigma <- unname(C %*% lb$V %*% t(C))
  se <- unname(sqrt(pmax(diag(Sigma) + lb$extra_var, 0)))

  # A nonlinear body is the case core refuses se.fit for, so there is no
  # second number to check against. Everything else is checked.
  if (sp_is_nl(fit, dpar, resp)) {
    return(list(eta = lb$eta, C = C, V = lb$V, Sigma = Sigma, se = se,
                se_ref = rep(NA_real_, length(se)), rel = NA_real_,
                n_predict = 0L, newdata = newdata, dpar = dpar,
                resp = resp, re.form = re.form, fit = fit))
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
       fit = fit)
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
#' @noRd
sp_sim_crit <- function(S, div, nsim, level, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- nrow(S)
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = m)
  z <- matrix(stats::rnorm(m * nsim), m, nsim)
  mx <- apply(abs((L %*% z) / div), 2L, max)
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
