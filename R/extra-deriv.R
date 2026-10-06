# The covariance of a DERIVATIVE of the variance that is not coefficient
# uncertainty, built from its sources rather than by differencing
# extra_cov. A central difference of a covariance divides its rounding
# by e^2 at order 1 and by e^4 at order 2; at frmtmb.spline's design
# step (1e-4 of the range at order 2) that left the second derivative's
# extra variance at the size of its own rounding, a standard error of 0
# past an exact gp()'s data and up to 1.5 at an unseen level of a term
# linear in the variable, whose true value is 0
# (dev/reviews/2026-10-05-gpby.md, B1).

#' The covariance of a curve's derivative that is not coefficient
#' uncertainty
#'
#' The derivative counterpart of `extra_cov` from [frm_lp_basis()]: the
#' `n x n` covariance, over the rows of `newdata`, of the `order`-th
#' derivative in `var` of the part of the linear predictor that is not
#' a function of the coefficients. Two sources contribute, and each is
#' differentiated before its covariance is formed:
#'
#' * An exact [gp()] term's kriging residual. Its derivative's
#'   covariance is the kernel's closed form,
#'   `d^o_x d^o_x' k(x, x') - d^o_x k(x, X) K^-1 d^o_x' k(X, x')`, from
#'   the squared exponential kernel's Hermite-polynomial derivatives. It
#'   is not zero at an observed position: the field's slope there is
#'   not fixed by its value. The nugget is white noise and has no
#'   derivative, so it is not in this covariance.
#' * A grouping level the fit never saw (`allow_new_levels = TRUE`).
#'   The design rows of its draw are differenced with the central
#'   stencil of step `eps`, as a design is, and the covariance is then
#'   `(L M) S (L M)'`, which carries no differenced rounding.
#'
#' A Hilbert-space `gp(k = )` and a smooth have no such part: their
#' uncertainty is all in their coefficients.
#'
#' @inheritParams frm_lp_basis
#' @param newdata The grid, a data frame.
#' @param var The variable to differentiate in, a column of `newdata`.
#' @param order 1 or 2.
#' @param eps The stencil's step for a new level's design rows.
#' @return An `n x n` matrix, as `extra_cov` is: a sparse
#'   `Matrix::sparseMatrix()` where the sources cover few of its entries
#'   (an empty one when there is no source), dense otherwise.
#' @seealso [frm_lp_basis()], whose `extra_cov` this differentiates.
#' @examples
#' set.seed(5)
#' d <- data.frame(x = sort(runif(40, 0, 4)))
#' d$y <- sin(d$x) + rnorm(40, 0, 0.2)
#' fit <- frm(bf(y ~ gp(x)), data = d)
#' nd <- data.frame(x = c(4.5, 5, 6))
#' # the slope's and the curvature's variance past the data
#' diag(frm_extra_cov_deriv(fit, nd, var = "x", order = 1))
#' diag(frm_extra_cov_deriv(fit, nd, var = "x", order = 2))
#' @export
frm_extra_cov_deriv <- function(object, newdata, var, order = 1L,
                                eps = NULL, dpar = NULL, resp = NULL,
                                re_formula = NULL,
                                allow_new_levels = FALSE) {
  require_frmtmb_fit(object, "frm_extra_cov_deriv()")
  require_fitted(object, "frm_extra_cov_deriv()")
  check_flag(allow_new_levels, "allow_new_levels")
  if (!is.data.frame(newdata) || !nrow(newdata)) {
    frm_stop("frm_extra_cov_deriv(): `newdata` must be a data frame with ",
             "at least one row, not ", arg_desc(newdata), call. = FALSE)
  }
  if (!is.character(var) || length(var) != 1L || !var %in% names(newdata)) {
    frm_stop("frm_extra_cov_deriv(): `var` must name one column of ",
             "`newdata`", call. = FALSE)
  }
  if (!isTRUE(order %in% c(1, 2))) {
    frm_stop("frm_extra_cov_deriv(): `order` must be 1 or 2",
             call. = FALSE)
  }
  order <- as.integer(order)
  x <- newdata[[var]]
  if (!is.numeric(x)) {
    frm_stop("frm_extra_cov_deriv(): `var` must be numeric", call. = FALSE)
  }
  if (is.null(eps)) {
    r <- diff(range(x))
    if (!is.finite(r) || r <= 0) r <- max(abs(x), 1)
    eps <- r * if (order == 1L) 1e-6 else 1e-4
  }
  check_re_form(re_formula)
  rr <- re_resolve(object, re_formula, "frm_extra_cov_deriv()")
  object <- rr$fit
  use_re <- re_form_keeps(rr$re_formula)
  resp <- resp %||% names(object$spec$responses)[1]
  rspec <- object$spec$responses[[resp]]
  if (is.null(rspec)) stop_unknown_response(object, resp)
  dpar <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1]
  lp <- object$frame[["linpreds"]][[linpred_key(resp, dpar)]]
  if (is.null(lp)) {
    frm_stop("frm_extra_cov_deriv(): unknown dpar '", dpar,
             "' for response '", resp, "'", call. = FALSE)
  }
  n <- nrow(newdata)
  blocks <- list()
  # a nonlinear body refuses both sources (lp_basis_nl())
  if (is.null(lp[["nl_body"]])) {
    env <- object$spec$responses[[lp[["resp"]]]]$formula_env
    for (gi in lp[["gps"]] %||% list()) {
      if (!identical(gi$type, "exact")) next
      b <- gp_krig_deriv_cov(object, gi, newdata, env, var, order)
      if (!is.null(b)) blocks[[length(blocks) + 1L]] <- b
    }
    if (use_re) {
      blocks <- c(blocks, new_level_deriv_cov(object, lp, newdata, use_re,
                                              allow_new_levels, var, order,
                                              eps))
    }
  }
  extra_cov_assemble(blocks, n)
}

#' Probabilists' Hermite polynomials, the derivatives of the squared
#' exponential kernel: `d^n/du^n exp(-u^2 / 2) = (-1)^n He_n(u) exp(...)`.
#'
#' @noRd
hermite_he <- function(n, z) {
  switch(n + 1L, rep(1, length(z)), z, z^2 - 1, z^3 - 3 * z,
         z^4 - 6 * z^2 + 3)
}

#' The covariance of the `o`-th derivative in `var` of one exact gp()
#' sub-GP's kriging residual at the rows of `newdata`, in closed form, as
#' a block `list(rows, C)` over the rows the sub-GP reaches, or `NULL`
#' when it reaches none or does not depend on `var`.
#'
#' With `u = (x - x') / l` along `var`'s dimension, the kernel's
#' derivatives are `d^o_x c = (-1)^o He_o(u) / l^o c` and
#' `d^o_x d^o_x' c = (-1)^o He_2o(u) / l^2o c`, times the kernel along
#' the other dimensions. The residual's derivative has covariance
#' `sd^2 (c^(o,o) - g (C + nugget I)^-1 g')`, `g` the cross derivative
#' against the fitted positions, times the rows' `by` multipliers.
#'
#' @noRd
gp_krig_deriv_cov <- function(fit, gi, newdata, env, var, o) {
  n <- nrow(newdata)
  vars <- vapply(gi$exprs, deparse1, "")
  d <- match(var, vars)
  uses <- vapply(gi$exprs, function(ex) var %in% all.vars(ex), NA)
  if (is.na(d)) {
    if (any(uses)) {
      frm_stop("frm_extra_cov_deriv(): `", var, "` enters a gp() term ",
               "through a transformed covariate (",
               paste(vars[uses], collapse = ", "), "), and its derivative ",
               "would need the chain rule through that transform. Write ",
               "the covariate as a column of its own", call. = FALSE)
    }
    return(NULL)
  }
  if (identical(gi$by$type, "numeric") &&
      identical(gi$by$label, var)) {
    frm_stop("frm_extra_cov_deriv(): `", var, "` is the numeric by ",
             "variable of a gp() term as well as its covariate, so the ",
             "derivative would need the product rule through the ",
             "multiplier", call. = FALSE)
  }
  bk <- fit$frame[["re_blocks"]][[gi$block_id]]
  w <- gp_by_mult(gi$by, newdata, env, n)
  on <- which(w != 0)
  if (!length(on)) return(NULL)
  Xc <- do.call(cbind, lapply(gi$exprs, function(ex) {
    as.numeric(eval(ex, newdata, env))
  }))[on, , drop = FALSE]
  pos <- gi$positions
  th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
  ell <- function(j) if (isTRUE(bk[["gp_iso"]])) exp(th[2]) else
    exp(th[1 + j])
  # the smooth kernel (no nugget) between two coordinate sets
  kern <- function(A, B) {
    Q <- 0
    for (j in seq_len(ncol(pos))) {
      Q <- Q + outer(A[, j], B[, j], "-")^2 / (2 * ell(j)^2)
    }
    exp(-Q)
  }
  l <- ell(d)
  g <- kern(Xc, pos) * ((-1)^o / l^o) *
    hermite_he(o, outer(Xc[, d], pos[, d], "-") / l)
  # He_2o is even and the kernel symmetric, so this is symmetric as
  # computed, and so is the A A' form below: no symmetrizing copy
  S <- kern(Xc, Xc) * ((-1)^o / l^(2L * o)) *
    hermite_he(2L * o, outer(Xc[, d], Xc[, d], "-") / l)
  thc <- th
  thc[1] <- 0
  K <- unname(covstruct_registry[["gp"]]$vcov(thc, bk))
  R <- tryCatch(chol(K), error = function(e) NULL)
  if (!is.null(R)) {
    S <- S - tcrossprod(t(backsolve(R, t(g), transpose = TRUE)))
  } else {
    S <- S - g %*% solve(K, t(g))
    S <- (S + t(S)) / 2
  }
  diag(S) <- pmax(diag(S), 0)
  S <- exp(2 * th[1]) * S
  if (any(w[on] != 1)) S <- S * outer(w[on], w[on])
  list(rows = on, C = S)
}

#' The covariance of the `o`-th derivative in `var` of every unseen
#' grouping level's draw at the rows of `newdata`: the central stencil
#' of step `eps` on the draw's design rows `M`, then `(L M) S (L M)'`.
#' One block `list(rows, C)` per draw.
#'
#' @noRd
new_level_deriv_cov <- function(fit, lp, newdata, use_re, allow_new_levels,
                                var, o, eps) {
  n <- nrow(newdata)
  lo <- hi <- newdata
  lo[[var]] <- newdata[[var]] - eps
  hi[[var]] <- newdata[[var]] + eps
  st <- rbind(lo, newdata, hi)
  ed <- lp_eta_design(fit, lp, st, use_re, allow_new_levels)
  ev <- lp_extra_var(fit, ed, use_re)
  out <- list()
  for (B in extra_var_blocks(ev$new_levels, 3L * n)) {
    # the grid rows whose stencil touches this draw
    r <- sort(unique((B$rows - 1L) %% n + 1L))
    M <- B$M
    Md <- if (o == 1L) {
      (M[2L * n + r, , drop = FALSE] - M[r, , drop = FALSE]) / (2 * eps)
    } else {
      (M[2L * n + r, , drop = FALSE] - 2 * M[n + r, , drop = FALSE] +
         M[r, , drop = FALSE]) / eps^2
    }
    C <- Md %*% B$S %*% t(Md)
    out[[length(out) + 1L]] <- list(rows = r, C = (C + t(C)) / 2)
  }
  out
}
