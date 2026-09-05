# Where a curve peaks, and where it crosses a level, with standard
# errors.
#
# These are the numbers a movement paper reports: the time of peak
# velocity, movement onset as the crossing of a speed threshold, the time
# of peak acceleration. Each is a ROOT of some derivative of the fitted
# curve, so each is a random variable with a distribution, and the
# implicit-function theorem turns the curve's own covariance into that
# distribution's variance.
#
# For a peak, t* solves f'(t*) = 0. Perturbing the coefficients moves
# f', and t* moves to keep the root: dt*/dc = -f''(t*)^-1 df'(t*)/dc, so
#
#     var(t*) = var(f'(t*)) / f''(t*)^2.
#
# For a crossing of a level a, t* solves f(t*) = a and the same argument
# gives var(t*) = var(f(t*)) / f'(t*)^2.
#
# Both are first-order and both degrade in the same place: a flat peak
# (f'' near zero) or a shallow crossing (f' near zero) has an
# ill-determined location, and the delta method says so by returning a
# large standard error rather than by failing. What it CANNOT say is
# that the root has stopped being unique, so the returned interval is
# reported per root and the search reports every root it found.

#' Features of a fitted curve: a peak, a trough, or a level crossing
#'
#' Locates a stationary point or a level crossing of a fitted curve and
#' gives its position a standard error, by the implicit-function delta
#' method.
#'
#' The search is a scan of the grid for a sign change, then Newton
#' refinement on the fitted curve, then one design pass at the located
#' root for the variance. Every root the grid brackets is returned, so a
#' curve with two peaks gives two rows; a curve with none gives a
#' zero-row result rather than an error, because "this curve does not
#' peak in this window" is an answer.
#'
#' At a stationary point the curve's own value has a delta-method
#' simplification worth knowing: the derivative of `f(t*)` with respect
#' to the coefficients is `df/dc + f'(t*) dt*/dc`, and `f'(t*)` is zero
#' there, so the standard error of the PEAK HEIGHT is just the pointwise
#' standard error of the curve at `t*`. It is reported as `.value_se`,
#' and it is not inflated by the uncertainty in the peak's location.
#'
#' @param object A `frmtmb_fit`, or a `frmtmb_curve` from [frm_curve()].
#' @param var Name of the covariate the feature is located along.
#' @param type `"maximum"`, `"minimum"` or `"extremum"` for a stationary
#'   point of the given kind; `"crossing"` for the points where the curve
#'   passes `at`.
#' @param at The level to cross, for `type = "crossing"`.
#' @param newdata The grid the search scans, and the values every other
#'   covariate is held at. Required when `object` is a fit.
#' @param eps Step size for the differences. `NULL` is the measured
#'   default of [frm_curve_deriv()].
#' @param maxit Newton iterations allowed per root.
#' @inheritParams frm_curve
#'
#' @return A data frame with one row per root: `.feature`, `.var`,
#'   `.estimate` (the located position), `.se`, `.lower_ci`, `.upper_ci`,
#'   `.value` (the curve there) and `.value_se`. The `"check"` attribute
#'   carries the covariance agreement, as [frm_curve()]'s does.
#'
#' @seealso [frm_curve()], [frm_curve_deriv()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(300)))
#' dd$y <- 2 * sin(pi * dd$x) + rnorm(300, 0, 0.3)
#' fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
#'                    family = stats::gaussian(), data = dd)
#' # 2 sin(pi x) peaks at x = 0.5
#' frm_curve_feature(fit, var = "x", type = "maximum",
#'                   newdata = data.frame(x = seq(0.05, 0.95, length.out = 41)))
#' @export
frm_curve_feature <- function(object, var,
                              type = c("maximum", "minimum", "extremum",
                                       "crossing"),
                              at = 0, newdata = NULL, dpar = NULL,
                              resp = NULL, re.form = NA, level = 0.95,
                              eps = NULL, maxit = 50L, tol = 1e-6) {
  type <- match.arg(type)
  sp_check_level(level)
  sp_check_count(maxit, "maxit")
  if (!is.numeric(at) || length(at) != 1L || !is.finite(at)) {
    stop("`at` must be one finite number: the level the crossing is of",
         call. = FALSE)
  }
  sp <- sp_spec(object, newdata, dpar, resp, re.form)
  # before the root scan, not after: a search that finds no root returns
  # early and would otherwise never reach the covariance
  sp_rp_gate(sp$fit)
  nd <- sp$newdata
  sp_check_var(nd, var)
  x <- as.numeric(nd[[var]])
  if (length(x) < 3L) {
    stop("A feature search needs a grid of at least three points along '",
         var, "': it scans for a sign change before it refines one",
         call. = FALSE)
  }
  ord <- order(x)
  nd <- nd[ord, , drop = FALSE]
  x <- x[ord]
  e1 <- sp_eps(eps, x, 1L)
  e2 <- sp_eps(eps, x, 2L)
  row1 <- nd[1L, , drop = FALSE]

  # g() is the function whose root is wanted, evaluated with one
  # predict() call per point set: the curve itself for a crossing, its
  # first derivative for a stationary point.
  eta_at <- function(tv) {
    d <- row1[rep(1L, length(tv)), , drop = FALSE]
    d[[var]] <- tv
    sp_predict_eta(sp$fit, d, sp$dpar, sp$resp, sp$re.form)
  }
  gfun <- if (type == "crossing") {
    function(tv) eta_at(tv) - at
  } else {
    function(tv) {
      v <- eta_at(c(tv - e1, tv + e1))
      k <- length(tv)
      (v[k + seq_len(k)] - v[seq_len(k)]) / (2 * e1)
    }
  }
  gp <- if (type == "crossing") {
    function(tv) {
      v <- eta_at(c(tv - e1, tv + e1))
      k <- length(tv)
      (v[k + seq_len(k)] - v[seq_len(k)]) / (2 * e1)
    }
  } else {
    function(tv) {
      v <- eta_at(c(tv - e2, tv, tv + e2))
      k <- length(tv)
      (v[2L * k + seq_len(k)] - 2 * v[k + seq_len(k)] + v[seq_len(k)]) / e2^2
    }
  }

  gv <- gfun(x)
  cross <- which(gv[-length(gv)] * gv[-1L] < 0)
  if (type %in% c("maximum", "minimum")) {
    want <- if (type == "maximum") 1 else -1
    # a maximum is where the slope falls through zero from above
    cross <- cross[sign(gv[cross] - gv[cross + 1L]) == want]
  }
  roots <- numeric(0)
  for (i in cross) {
    lo <- x[i]
    hi <- x[i + 1L]
    # linear interpolation starts Newton inside the bracket; Newton is
    # kept inside it, because a near-zero second derivative would
    # otherwise throw the step out of the window the user asked about
    t <- lo + (hi - lo) * gv[i] / (gv[i] - gv[i + 1L])
    for (k in seq_len(maxit)) {
      g <- gfun(t)
      d <- gp(t)
      if (!is.finite(d) || d == 0) break
      tn <- t - g / d
      if (!is.finite(tn) || tn < lo || tn > hi) {
        tn <- (lo + hi) / 2
      }
      if (g * gv[i] < 0) hi <- t else lo <- t
      if (abs(tn - t) <= 1e-12 * max(abs(t), diff(range(x)))) {
        t <- tn
        break
      }
      t <- tn
    }
    roots <- c(roots, t)
  }
  out <- data.frame(.feature = character(0), .var = character(0),
                    .estimate = numeric(0), .se = numeric(0),
                    .lower_ci = numeric(0), .upper_ci = numeric(0),
                    .value = numeric(0), .value_se = numeric(0))
  ck <- list(cov_rel_error = NA_real_, n_predict = NA_integer_,
             crit_mcse = NA_real_)
  if (!length(roots)) {
    return(structure(out, class = c("frmtmb_feature", "data.frame"),
                     check = ck, level = level, type = type,
                     row.names = integer(0)))
  }

  # One design pass over the whole five-point stencil at every root: the
  # variances below are exact functionals of the same C and V the bands
  # use, so a feature and a band on one fit cannot disagree about the
  # covariance.
  stk <- row1[rep(1L, 5L * length(roots)), , drop = FALSE]
  stk[[var]] <- c(roots - e2, roots - e1, roots, roots + e1, roots + e2)
  parts <- sp_curve_parts(sp$fit, stk, sp$dpar, sp$resp, sp$re.form, tol)
  nr <- length(roots)
  blk <- function(k) parts$C[(k - 1L) * nr + seq_len(nr), , drop = FALSE]
  eta <- function(k) parts$eta[(k - 1L) * nr + seq_len(nr)]
  D1 <- (blk(4L) - blk(2L)) / (2 * e1)
  D2 <- (blk(5L) - 2 * blk(3L) + blk(1L)) / e2^2
  f1 <- (eta(4L) - eta(2L)) / (2 * e1)
  f2 <- (eta(5L) - 2 * eta(3L) + eta(1L)) / e2^2
  f0 <- eta(3L)
  C0 <- blk(3L)
  qf <- function(A) sqrt(pmax(rowSums((A %*% parts$V) * A), 0))
  if (type == "crossing") {
    denom <- f1
    se_t <- qf(C0) / abs(denom)
  } else {
    denom <- f2
    se_t <- qf(D1) / abs(denom)
  }
  crit <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(
    .feature = type, .var = var, .estimate = roots, .se = se_t,
    .lower_ci = roots - crit * se_t, .upper_ci = roots + crit * se_t,
    .value = f0, .value_se = qf(C0), stringsAsFactors = FALSE)
  structure(out, class = c("frmtmb_feature", "data.frame"),
            check = list(cov_rel_error = parts$rel,
                         n_predict = parts$n_predict,
                         crit_mcse = NA_real_),
            level = level, type = type,
            slope = f1, curvature = f2,
            row.names = seq_len(nrow(out)))
}

#' @export
print.frmtmb_feature <- function(x, ...) {
  ck <- attr(x, "check")
  cat("<frmtmb curve feature> ", attr(x, "type"), ", ", nrow(x),
      " found, level ", attr(x, "level"), "\n", sep = "")
  if (nrow(x)) {
    cat("  covariance checked against predict(se.fit = TRUE) to ",
        format(ck$cov_rel_error, digits = 3), " relative\n", sep = "")
    print(as.data.frame(x))
  } else {
    cat("  the grid brackets no sign change, so the curve has no ",
        attr(x, "type"), " in this window\n", sep = "")
  }
  invisible(x)
}
