# Derivatives of a fitted curve, and the step size that is not a
# constant.
#
# The derivative of a curve is a LINEAR functional of the same
# coefficients, so its standard error is an ordinary delta method once
# the derivative of the design is in hand. This package takes that
# derivative by central differences of the design rows, which is what
# gratia does, but not with gratia's step size: a fixed eps = 1e-7 is
# fine at order one and catastrophic at order two, because the second
# difference divides by eps^2 and turns floating-point cancellation into
# the dominant error. Measured against the exact derivative of
# f(t) = 2 sin(pi t) + 0.6 t on [0.05, 0.95]:
#
#   eps      max |d1 error|   max |d2 error|
#   1e-04    1.021e-07        2.234e-07
#   1e-06    4.033e-10        8.219e-04
#   1e-07    5.660e-09        8.165e-02
#   1e-08    4.652e-08        1.007e+01
#
# So the default step is scaled by the covariate's range and set per
# order: 1e-6 of the range for the first derivative, 1e-4 for the
# second, each within a factor of about 15 of its own optimum.
#
# Exact basis derivatives were the alternative and are not reachable:
# the design is rebuilt through predict(), which evaluates a basis and
# never differentiates one. Differentiating the basis itself would mean
# reading mgcv smooth objects out of the fitted frame, which is a deeper
# reach into core than this package makes anywhere else, for an error
# already at the tenth significant figure.

#' Derivatives of a fitted curve
#'
#' The first or second derivative of a fitted linear predictor with
#' respect to one covariate, on a grid, with delta-method standard errors
#' and the same pair of intervals [frm_curve()] gives: pointwise, and
#' simultaneous over the whole grid.
#'
#' The derivative is taken by central differences of the DESIGN, not of
#' the fitted values, so the estimate and its standard error describe one
#' function rather than two. Writing `D` for the differenced design and
#' `V` for the joint coefficient covariance, the reported curve is `D c`
#' and its covariance is `D V D'`, which is exact for the differenced
#' basis; the only approximation is the difference itself, and `eps`
#' controls that.
#'
#' @param object A `frmtmb_fit`, or a `frmtmb_curve` from [frm_curve()],
#'   in which case its grid and its `dpar`, `resp` and `re.form` are
#'   reused.
#' @param var Name of the covariate to differentiate with respect to. It
#'   must be a numeric column of the grid.
#' @param order 1 or 2.
#' @param newdata The grid. Required when `object` is a fit; taken from
#'   the curve otherwise.
#' @param eps Step size. `NULL`, the default, is `1e-6` of the grid's
#'   range at `order = 1` and `1e-4` of it at `order = 2`.
#' @inheritParams frm_curve
#'
#' @return A `frmtmb_curve` data frame, as [frm_curve()] returns, whose
#'   `.estimate` is the derivative.
#'
#' @seealso [frm_curve()], [frm_curve_feature()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200)))
#' dd$y <- 2 * sin(pi * dd$x) + rnorm(200, 0, 0.4)
#' fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
#'                    family = stats::gaussian(), data = dd)
#' g <- data.frame(x = seq(0.05, 0.95, length.out = 19))
#' d1 <- frm_curve_deriv(fit, var = "x", order = 1, newdata = g,
#'                       simultaneous = FALSE)
#' # the derivative of 2 sin(pi x) is 2 pi cos(pi x)
#' head(cbind(g, fitted = d1$.estimate, truth = 2 * pi * cos(pi * g$x)))
#' @export
frm_curve_deriv <- function(object, var, order = 1L, newdata = NULL,
                            dpar = NULL, resp = NULL, re.form = NA,
                            level = 0.95, simultaneous = TRUE,
                            nsim = 10000L, eps = NULL, seed = NULL,
                            tol = 1e-6) {
  sp_check_level(level)
  sp_check_flag(simultaneous, "simultaneous")
  if (!identical(order, 1L) && !identical(order, 2L) &&
      !identical(order, 1) && !identical(order, 2)) {
    stop("`order` must be 1 (the slope) or 2 (the curvature). Higher ",
         "orders are not offered: a third central difference divides by ",
         "eps^3 and there is no step size at which it is accurate",
         call. = FALSE)
  }
  order <- as.integer(order)
  sp <- sp_spec(object, newdata, dpar, resp, re.form)
  sp_rp_gate(sp$fit)
  nd <- sp$newdata
  sp_check_var(nd, var)
  x <- as.numeric(nd[[var]])
  e <- sp_eps(eps, x, order)

  # One perturbation pass over the stacked grid: the design at the three
  # stencil positions costs the same number of predict() calls as the
  # design at one, because a call already returns every row.
  stack <- rbind(nd, nd, nd)
  stack[[var]] <- c(x - e, x, x + e)
  parts <- sp_curve_parts(sp$fit, stack, sp$dpar, sp$resp, sp$re.form, tol)
  m <- nrow(nd)
  lo <- seq_len(m)
  mid <- m + lo
  hi <- 2L * m + lo
  if (order == 1L) {
    D <- (parts$C[hi, , drop = FALSE] - parts$C[lo, , drop = FALSE]) / (2 * e)
    est <- (parts$eta[hi] - parts$eta[lo]) / (2 * e)
  } else {
    D <- (parts$C[hi, , drop = FALSE] - 2 * parts$C[mid, , drop = FALSE] +
            parts$C[lo, , drop = FALSE]) / e^2
    est <- (parts$eta[hi] - 2 * parts$eta[mid] + parts$eta[lo]) / e^2
  }
  Sigma <- D %*% parts$V %*% t(D)
  se <- sqrt(pmax(diag(Sigma), 0))
  parts$newdata <- nd
  out <- sp_assemble(parts, est, se, Sigma, level, simultaneous, nsim,
                     FALSE, seed, nd,
                     what = paste0("derivative of order ", order,
                                   " in ", var))
  attr(out, "eps") <- e
  out
}

#' Resolve the (fit, grid, predictor) triple from either input form.
#'
#' @noRd
sp_spec <- function(object, newdata, dpar, resp, re.form) {
  if (inherits(object, "frmtmb_curve")) {
    s <- attr(object, "spec")
    nd <- if (is.null(newdata)) s$newdata else newdata
    return(list(fit = attr(object, "fit"), newdata = nd, dpar = s$dpar,
                resp = s$resp, re.form = s$re.form))
  }
  if (is.null(newdata)) {
    stop("`newdata` is required when the first argument is a fit: it is ",
         "the grid the curve is evaluated on. Pass a frmtmb_curve from ",
         "frm_curve() to reuse a grid instead", call. = FALSE)
  }
  list(fit = object, newdata = newdata, dpar = dpar, resp = resp,
       re.form = re.form)
}

#' @noRd
sp_check_var <- function(nd, var) {
  if (!is.character(var) || length(var) != 1L) {
    stop("`var` must name one covariate of the grid, as a string",
         call. = FALSE)
  }
  if (is.null(nd[[var]])) {
    stop("`var` names '", var, "', which is not a column of the grid. ",
         "The grid has: ", paste(names(nd), collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(nd[[var]])) {
    stop("A curve is differentiated with respect to a NUMERIC covariate ",
         "and '", var, "' is ", class(nd[[var]])[1L],
         ". A factor has no derivative; take a contrast instead",
         call. = FALSE)
  }
  invisible(NULL)
}

#' The measured default step: a fraction of the covariate's own range,
#' per order. See the note at the top of this file for the sweep.
#'
#' @noRd
sp_eps <- function(eps, x, order) {
  if (!is.null(eps)) {
    if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) ||
        eps <= 0) {
      stop("`eps` must be one positive finite number, or NULL for the ",
           "measured default", call. = FALSE)
    }
    return(eps)
  }
  r <- diff(range(x))
  if (!is.finite(r) || r <= 0) r <- max(abs(x), 1)
  r * if (order == 1L) 1e-6 else 1e-4
}
