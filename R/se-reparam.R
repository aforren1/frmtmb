# The outer Hessian in coordinates where the fixed designs are well
# conditioned.
#
# On a raw polynomial of degree 8 on [1, 2] the Hessian scaled to unit
# diagonal has a condition number near 4e16, above 1 / machine epsilon.
# Formed in the coefficients' own coordinates its smallest eigenvalue is
# rounding, so the SE check could only call those coefficients
# "illcond" and give them no SE, while lm() reports them through its QR
# to 3.6e-9 of the exact covariance (dev/illse-floor.R). The Hessian
# formed along the columns of `M`, the inverse of each design's R
# factor, is close to the identity, and `M (M' H M)^-1 M'` is the same
# covariance as `H^-1` for any invertible `M`: the coordinates change
# only what rounding costs.

#' A basis of the outer parameters in which each linear predictor's
#' fixed design has orthogonal columns.
#'
#' For each linear predictor whose coefficients are outer parameters and
#' enter through a dense design `X` (not a nonlinear body, not a `mo()`
#' term, whose column moves with its simplex), the block of `M` at its
#' positions is `T` with `X T = Q sqrt(n)`: orthogonal columns of unit
#' root mean square, so a unit step in each new coordinate moves the
#' predictor by about one unit per row. Every other parameter keeps its
#' own coordinate. The columns are scaled before the QR so that the
#' factor sees the design's shape and not its units, as lm() pivots on
#' them. A design whose QR finds it rank deficient keeps its own
#' coordinates, since its alias is an exact flat direction judged
#' elsewhere (se_design_null()), and so does one whose condition number
#' is below the square root of `se_reparam_kappa`, which is not
#' ill-conditioned (design_whitener()). Rows are thinned evenly above
#' `frmtmb.unit_ls_cells`, as predictor_design_kappa() does; `T` only
#' has to make the Hessian well conditioned, not equal to the identity.
#'
#' Returns `M` (`np x np`) with attribute `pos`, the positions its
#' blocks cover, or `NULL` when no design gets a block.
#'
#' @noRd
reparam_basis <- function(frame, pn,
                          cells = getOption("frmtmb.unit_ls_cells", 5e6)) {
  np <- length(pn)
  M <- diag(1, np)
  covered <- integer(0)
  for (lp in frame[["linpreds"]] %||% list()) {
    if (!is.null(lp[["constant"]]) || !is.null(lp[["nl_body"]]) ||
          length(lp[["mo"]] %||% list()) ||
          !lp[["par"]] %in% c("beta", "betad")) {
      next
    }
    X <- lp[["X"]]
    if (!design_ok(X)) next
    pos <- outer_index_of(lp, pn, frame)
    if (length(pos) != ncol(X)) next
    # an ordinal predictor has no intercept: a common shift of the
    # thresholds is one, and a raw polynomial without a constant column
    # nearly contains the constant (on a cumulative degree-8 fit the
    # designs' coordinates alone left a curvature at 1.3e-11 of the
    # largest, dev/illse-glm.R), so the shift joins the design
    sh <- ord_shift(frame, lp, pn)
    if (!is.null(sh)) {
      Ta <- tryCatch(design_whitener(cbind(-1, as.matrix(X)), cells),
                     error = function(e) NULL)
      if (is.null(Ta)) next
      tp <- sh$pos
      j0 <- which.max(abs(sh$s))
      blk <- c(tp, pos)
      B <- matrix(0, length(blk), length(blk))
      k <- ncol(X) + 1L
      B[seq_along(tp), seq_len(k)] <- outer(sh$s, Ta[1L, ])
      B[length(tp) + seq_along(pos), seq_len(k)] <- Ta[-1L, ]
      rest <- setdiff(seq_along(tp), j0)
      B[cbind(rest, k + seq_along(rest))] <- 1
      M[blk, blk] <- B
      covered <- c(covered, blk)
      next
    }
    if (ncol(X) < 2L) next
    Tb <- tryCatch(design_whitener(X, cells), error = function(e) NULL)
    if (is.null(Tb)) next
    M[pos, pos] <- Tb
    covered <- c(covered, pos)
  }
  if (length(covered)) structure(M, pos = sort(unique(covered)))
}

#' Where an ordinal response's thresholds sit among the outer parameters
#' `pn`, and the direction `s` in them that shifts every threshold by
#' one, read off the family's own maps (`post$ord_thresholds` and its
#' inverse): the first coordinate alone for ordered or equidistant
#' thresholds, all of them for free ones. NULL for any other predictor,
#' and for thresholds that vary by group, which have no single shift.
#'
#' @noRd
ord_shift <- function(frame, lp, pn) {
  if (!identical(lp[["dpar"]], "mu")) return(NULL)
  rs <- frame[["spec"]]$responses
  if (length(rs) != 1L) return(NULL)
  post <- rs[[1L]]$family[["post"]]
  fwd <- post[["ord_thresholds"]]
  inv <- post[["ord_thresholds_raw"]]
  tp <- which(pn == "tau_raw")
  if (!is.function(fwd) || !is.function(inv) || !length(tp)) return(NULL)
  raw0 <- numeric(length(tp))
  s <- tryCatch(inv(fwd(raw0) + 1) - raw0, error = function(e) NULL)
  if (length(s) != length(tp) || !all(is.finite(s)) || !any(s != 0)) {
    return(NULL)
  }
  list(pos = tp, s = s)
}

#' `T` with `X T` orthogonal, of unit root mean square per column, from
#' the QR of the column-scaled design; NULL when the QR is rank
#' deficient, the design's condition number is below the square root of
#' `se_reparam_kappa` (it is not ill-conditioned), or `T` is not finite.
#'
#' @noRd
design_whitener <- function(X, cells) {
  p <- ncol(X)
  n <- nrow(X)
  m <- max(2L * p, floor(cells / p))
  rows <- if (n > m) unique(round(seq(1, n, length.out = m))) else
    seq_len(n)
  Xk <- as.matrix(X[rows, , drop = FALSE])
  if (!all(is.finite(Xk))) return(NULL)
  s <- sqrt(colSums(Xk^2) / nrow(Xk))
  if (!all(s > 0)) return(NULL)
  q <- base::qr(Xk / rep(s, each = nrow(Xk)))
  if (q$rank < p) return(NULL)
  R <- qr.R(q)
  # A design that is not itself ill-conditioned keeps its coordinates:
  # the ill-conditioning is then not its, and new coordinates for it
  # would only rescale a direction no design makes weak. Under OpenBLAS
  # the nonlinear alias a_k + b of test-se-tier.R (a ~ 0 + f, a design
  # with condition number 1) measured just above se_reparam_floor in
  # such coordinates and kept SEs on a flat ridge
  if (1 / rcond(R, triangular = TRUE) < sqrt(se_reparam_kappa)) {
    return(NULL)
  }
  # X[, piv] / s[piv] = Q R, so X T = Q sqrt(rows) with T[piv, ] = R^-1
  # scaled back by s
  Ri <- backsolve(R, diag(p)) * sqrt(nrow(Xk))
  Tb <- matrix(0, p, p)
  Tb[q$pivot, ] <- Ri / s[q$pivot]
  if (all(is.finite(Tb))) Tb
}

#' The outer Hessian in the coordinates of `M`, `M' H M`, by central
#' differences of the gradient along each column of `M`, with the step
#' optimHess() takes (1e-3), here in units that move every predictor by
#' about one unit per row. The differences are never formed in the
#' natural coordinates and then rotated: there the curvature of the
#' weak directions is below the rounding of the larger ones.
#'
#' Returns `list(H, E)`, `E` half the asymmetry before symmetrizing (the
#' differences' own noise), or NULL when a gradient is not finite.
#'
#' @noRd
reparam_hessian <- function(obj, p, M, step = 1e-3, richardson = TRUE) {
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  k <- ncol(M)
  H <- matrix(0, k, k)
  cd <- function(j, h) {
    gp <- tryCatch(se_quiet(as.numeric(obj$gr(p + h * M[, j]))),
                   error = function(e) NULL)
    gm <- tryCatch(se_quiet(as.numeric(obj$gr(p - h * M[, j]))),
                   error = function(e) NULL)
    if (length(gp) != length(p) || length(gm) != length(p) ||
          !all(is.finite(gp)) || !all(is.finite(gm))) {
      return(NULL)
    }
    crossprod(M, gp - gm) / (2 * h)
  }
  for (j in seq_len(k)) {
    c1 <- cd(j, step)
    if (is.null(c1)) return(NULL)
    if (richardson) {
      c2 <- cd(j, 2 * step)
      if (is.null(c2)) return(NULL)
      c1 <- (4 * c1 - c2) / 3
    }
    H[, j] <- c1
  }
  list(H = (H + t(H)) / 2, E = abs(H - t(H)) / 2)
}

#' reparam_hessian(), and once more along the eigenvectors of its result
#' when that is still ill-conditioned (above `se_reparam_refine` at unit
#' diagonal). The designs' coordinates do not see a collinearity between
#' a design and a parameter outside it: a cumulative fit's columns
#' `x, ..., x^8` without an intercept nearly contain the constant that
#' its first threshold carries, which left one curvature at 1.3e-11 of
#' the largest (dev/illse-glm.R). That direction is resolved to the
#' differences' accuracy, so the second pass, scaled to unit curvature
#' along it, is well conditioned.
#'
#' Returns `list(H, E, M)` with `M` the coordinates `H` is in, or NULL.
#'
#' The last result is kept on the objective (`obj$env`), keyed on `p`
#' and `M`: newton_polish(), conv_headroom(), grad_headroom() and
#' se_reparam_cov() each ask for it, mostly at the fit's own point. On
#' n = 1e5 with 100 parameters one took 100 to 200 s, and six of them
#' made the fit 2.5 times as long as on 0.70.0 (review of lane illse,
#' item 7).
#'
#' @noRd
reparam_hessian_refined <- function(obj, p, M) {
  env <- obj$env
  key <- unname(as.numeric(p))
  hit <- if (is.environment(env)) env$.frmtmb_reparam_h
  if (!is.null(hit) && identical(hit$p, key) && identical(hit$M, M)) {
    return(hit$h)
  }
  h <- reparam_hessian_refined_at(obj, p, M)
  if (is.environment(env)) {
    env$.frmtmb_reparam_h <- list(p = key, M = M, h = h)
  }
  h
}

#' `obj$he(p)`, the last one kept on the objective keyed on `p`: the
#' Newton steps, conv_headroom(), grad_headroom() and the SE check each
#' ask for it at the fit's own point, and on a wide design one costs as
#' much as several gradients of the whole data. An AD Hessian is a
#' function of `p` alone, so a kept one is the one a new call returns.
#' An error is not kept.
#'
#' @noRd
obj_he_cached <- function(obj, p) {
  env <- obj$env
  key <- unname(as.numeric(p))
  hit <- if (is.environment(env)) env$.frmtmb_he
  if (!is.null(hit) && identical(hit$p, key)) return(hit$H)
  H <- obj$he(p)
  if (is.environment(env)) env$.frmtmb_he <- list(p = key, H = H)
  H
}

#' reparam_hessian_refined() without its cache.
#'
#' @noRd
reparam_hessian_refined_at <- function(obj, p, M) {
  h <- reparam_hessian(obj, p, M)
  if (is.null(h) || !all(is.finite(h$H)) || !all(diag(h$H) > 0)) {
    return(NULL)
  }
  h$M <- M
  D <- sqrt(diag(h$H))
  e <- eigen(h$H / outer(D, D), symmetric = TRUE)
  big <- max(e$values)
  if (min(e$values) > big / se_reparam_refine) return(h)
  # coordinates W with W' H W = I: H = D U L U' D, W = D^-1 U L^-1/2. An
  # eigenvalue at the differences' own noise, or rounded to zero or
  # below, is scaled by that noise: its eigenvector is still accurate to
  # the noise over the gap, and the second pass measures the curvature
  # along it (one step from the stop of the cumulative fit above the
  # first pass gave the threshold direction -2.1e-12 of the largest,
  # against a noise of 3.3e-7; dev/illse-refine.R)
  nk <- sqrt(colSums(((h$E / outer(D, D)) %*% e$vectors)^2))
  lam <- pmax(e$values, nk, big * .Machine$double.eps)
  W <- sweep(e$vectors / D, 2L, sqrt(lam), "/")
  M2 <- M %*% W
  h2 <- reparam_hessian(obj, p, M2)
  if (is.null(h2) || !all(is.finite(h2$H)) || !all(diag(h2$H) > 0)) {
    return(h)
  }
  h2$M <- M2
  # the directions the second pass rescaled, in natural units: along a
  # flat one it reads the gradient's noise over a step that rescaling
  # made long as curvature (se_reparam_cov() asks the objective)
  h2$weak <- M %*% (e$vectors[, e$values <= big / se_reparam_refine,
                              drop = FALSE] / D)
  h2
}

#' The condition number at unit diagonal of the Hessian in the designs'
#' coordinates above which reparam_hessian_refined() measures it once
#' more: products in those coordinates keep about twelve digits below it.
#'
#' @noRd
se_reparam_refine <- 1e4

#' The covariance of the outer parameters in the coordinates of
#' reparam_basis(), where the fixed design makes the Hessian in natural
#' coordinates too ill-conditioned to invert.
#'
#' Asked by cov_from_hessian() before any of its tiers, so that the
#' answer does not depend on which tier would have accepted the natural
#' inverse: on degree 6 with (1 | g) solve() of the finite-difference
#' Hessian kept every SE at 0.3 to 1.8 percent from the exact ones,
#' and at degree 8 the SE check could only call nine coefficients
#' "illcond" (dev/illse-proto.R, dev/illse-proto3.R). Each condition
#' below is one under which the new coordinates are known to hold the
#' answer; where any fails the tiers run as before:
#'
#' - the fit has at most `se_reparam_np` outer parameters, so that the
#'   decomposition below and the design's QR cost little beside it;
#' - the Hessian scaled to unit diagonal has a direction at or below
#'   `1 / se_reparam_kappa` of its largest curvature, and every such
#'   direction lies mostly (`se_reparam_load`) in the coefficients the
#'   new coordinates cover: a weak direction of a variance component or
#'   a nonlinear parameter is not the design's, and the tiers judge it;
#' - in the new coordinates every curvature is above `se_reparam_floor`
#'   of the largest and above `se_noise_mult` times its own noise, so a
#'   flat direction, a ridge or an ascent the design does not cause
#'   still goes to the tiers and their verdicts;
#' - a Newton step from the estimates gains at most `se_ill_gain`: a fit
#'   that stopped short keeps the verdict that says so.
#'
#' Returns `list(V, shown, lost, tier = 4L, root)` with `root` a square
#' factor, `V = root root'`, or NULL. A product `a' V a` formed from `V`
#' in double precision loses the digits the coordinates kept (on
#' degree 8 the exact covariance rounded to double gave 2.8 times the
#' true variance of a fitted value, dev/illse-floor.R); `root' a` does
#' not, and the propagating routes read it (jc_rowvar()).
#'
#' `free`, the positions the covariance is over, is all of them when
#' cov_from_hessian() asks. se_tier3() asks with the parameters it
#' removed left out: a group sd at its boundary, a row it emptied, a
#' bound. The covariance of the others is then conditional on them, as
#' tier 3's own is, with zero rows and columns for them; it needs every
#' removed parameter to be outside the designs' blocks. On a degree-8
#' raw polynomial with `(1 | g)` and the sd at zero, tier 3 in natural
#' coordinates called x3 to x6 flat and kept the other coefficients'
#' SEs at 5 to 10 percent of lm()'s (review of lane illse, item 4).
#'
#' @noRd
se_reparam_cov <- function(fit, H, p, free = seq_len(NROW(H))) {
  np <- NROW(H)
  k <- length(free)
  if (k < 2L || np > se_reparam_np || length(p) != np ||
        !all(is.finite(H[free, free])) || !is.null(fit$importance) ||
        isTRUE(fit$quadrature)) {
    return(NULL)
  }
  Hf <- H[free, free, drop = FALSE]
  D <- sqrt(abs(diag(Hf)))
  if (!all(D > 0)) return(NULL)
  e <- eigen(Hf / outer(D, D), symmetric = TRUE)
  big <- max(abs(e$values))
  weak <- which(e$values <= big / se_reparam_kappa)
  if (!length(weak)) return(NULL)
  M <- tryCatch(reparam_basis(fit$frame, names(fit$opt$par)),
                error = function(e) NULL)
  if (is.null(M) || !all(dim(M) == np)) return(NULL)
  pos <- attr(M, "pos")
  out <- setdiff(seq_len(np), free)
  if (length(intersect(out, pos))) return(NULL)
  posf <- which(free %in% pos)
  load <- colSums(e$vectors[posf, weak, drop = FALSE]^2)
  if (any(load < se_reparam_load)) return(NULL)
  if (se_reparam_defer(fit, H, c(pos, out))) return(NULL)
  if (length(out)) M <- M[, free, drop = FALSE]
  saved <- obj_state_save(fit$obj)
  on.exit(obj_state_restore(fit$obj, saved), add = TRUE)
  h <- reparam_hessian_refined(fit$obj, p, M)
  if (is.null(h)) return(NULL)
  M <- h$M
  Dt <- sqrt(diag(h$H))
  et <- eigen(h$H / outer(Dt, Dt), symmetric = TRUE)
  nk <- sqrt(colSums(((h$E / outer(Dt, Dt)) %*% et$vectors)^2))
  if (min(et$values) <= se_reparam_floor * max(et$values) ||
        any(et$values <= se_noise_mult * nk)) {
    return(NULL)
  }
  # a direction the second pass rescaled must curve the objective itself
  # (se_obj_curv(), the step test's second differences): a curvature
  # the differences of the gradient made from their own noise does not
  Dn <- sqrt(abs(diag(H)))
  for (j in seq_len(NCOL(h$weak))) {
    dj <- h$weak[, j]
    dj <- dj / sqrt(sum((dj * Dn)^2))
    cj <- se_obj_curv(fit, p, dj)
    if (!is.finite(cj) || cj <= 0) return(NULL)
  }
  R <- tryCatch(chol(h$H), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  root <- M %*% backsolve(R, diag(k))
  g <- tryCatch(se_quiet(as.numeric(fit$obj$gr(p))), error = function(e) NULL)
  if (length(g) != np || !all(is.finite(g))) return(NULL)
  # what a Newton step over `free` gains: g' V g / 2, through the factor
  gain <- sum(crossprod(root, g)^2) / 2
  if (!is.finite(gain) || gain > se_ill_gain) return(NULL)
  V <- tcrossprod(root)
  dimnames(V) <- dimnames(H)
  rownames(root) <- rownames(H)
  # whether a' V a formed from the whole matrix keeps a few digits:
  # it loses about kappa * epsilon (dev/illse-floor.R)
  attr(root, "lossy") <- min(e$values) <= big / se_product_kappa
  none <- stats::setNames(character(0), character(0))
  list(V = V, shown = V, lost = none, tier = 4L, root = root,
       ev = et$values)
}

#' Must the tiers judge a parameter outside the designs' blocks (`pos`)?
#' The unit-diagonal test above sees only how a direction compares with
#' the others. A log standard deviation run to zero has a curvature near
#' zero in absolute terms, but it is decoupled from the rest, so at unit
#' diagonal its eigenvalue is near 1 and the route took it: on raw
#' polynomials with `(1 | g)` and no group effect, 2 of 80 fits gave
#' SE(log sd) of 4.75e5 and 3.5e5 with no warning, where 0.70.0 named
#' the sd (review of lane illse, item 4). So the route is refused when
#' a parameter outside `pos` has a row that tier 3 empties
#' (`se_empty_row`, in the optimizer's units), or is a log standard
#' deviation whose curvature is below `se_edge_screen`, where
#' se_edge_sd() and tier 3's "short" and "boundary" verdicts look.
#'
#' @noRd
se_reparam_defer <- function(fit, H, pos) {
  np <- NROW(H)
  out <- setdiff(seq_len(np), pos)
  if (!length(out)) return(FALSE)
  u <- rep_len(fit$par_units %||% 1, np)
  Hu <- abs(H * outer(u, u))
  if (any(apply(Hu[out, , drop = FALSE], 1L, max) <= se_empty_row)) {
    return(TRUE)
  }
  th <- fit$estimates[["theta"]]
  if (!length(th)) return(FALSE)
  om <- outer_par_map(fit)
  thpos <- which(om$comp == "theta")
  if (length(thpos) != length(th) || length(om$comp) != np) return(TRUE)
  sd_i <- tryCatch(log_sd_theta_index(fit), error = function(e) integer(0))
  j <- intersect(thpos[unique(sd_i)], out)
  any(diag(Hu)[j] <= se_edge_screen)
}

#' Does a Hessian (natural units, symmetric) have a direction at or
#' below `1 / se_reparam_kappa` of its largest curvature at unit
#' diagonal, so that the design's coordinates are worth trying? Only for
#' at most `se_reparam_np` parameters: the decomposition and the
#' design's QR then cost little beside the fit.
#'
#' @noRd
hess_needs_reparam <- function(H) {
  np <- NROW(H)
  if (np < 2L || np > se_reparam_np || !all(is.finite(H))) return(FALSE)
  D <- sqrt(abs(diag(H)))
  if (!all(D > 0)) return(FALSE)
  ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  # A direction that curves downward by more than se_ridge_screen of
  # the largest is a saddle, not ill-conditioning: test-grad-verdict.R
  # stops a fit one iteration from a start far from the maximum (an
  # eigenvalue of -0.42 against 2.46) and asks for its verdict there.
  # A finite-difference Hessian's own noise stays below that bound, and
  # must not keep a fit with (1 | g) from the steps it needs (poisson
  # degree 8, seed 2 of dev/illse-glm.R)
  big <- max(abs(ev))
  any(abs(ev) <= big / se_reparam_kappa) && all(ev >= -se_ridge_screen * big)
}

#' The most outer parameters se_reparam_cov() and the Newton steps take
#' the design's coordinates for, as newton_polish() caps its own steps.
#'
#' @noRd
se_reparam_np <- 100L

#' The condition number, at unit diagonal, of an outer Hessian above
#' which the covariance is taken in the design's coordinates
#' (se_reparam_cov()). solve() of an exact Hessian loses about
#' `kappa * epsilon` of each entry, 1e-8 here, and a finite-difference
#' one more; in the design's coordinates the loss on the raw
#' polynomials of dev/illse-proto.R was 2e-12 to 3e-6 at every degree
#' from 3 to 8.
#'
#' @noRd
se_reparam_kappa <- 1e8

#' The smallest curvature, as a fraction of the largest at unit
#' diagonal, that the design's coordinates must leave (se_reparam_cov()).
#' On the raw polynomials those coordinates put every curvature within
#' a factor of 1.5 of the largest (gaussian, with and without (1 | g));
#' a direction a millionth of the largest is not the design's.
#'
#' @noRd
se_reparam_floor <- 1e-6

#' The share of a weak direction (squared loading at unit diagonal) that
#' must lie in the coefficients reparam_basis() covers.
#'
#' @noRd
se_reparam_load <- 0.5

#' `diag(G V G')` over positions `pos` of a joint covariance `jc`
#' (get_joint_cov()). Where the fit carries a factor of its covariance
#' (se_reparam_cov()), through the factor: `V` itself cannot be
#' multiplied by a row of an ill-conditioned design in double
#' precision without losing every digit. Otherwise exactly the product
#' each caller formed before.
#'
#' @noRd
jc_rowvar <- function(jc, G, pos) {
  if (is.null(jc$root)) {
    V <- jc$V[pos, pos, drop = FALSE]
    return(rowSums((G %*% V) * G))
  }
  v <- rowSums((G %*% jc$root[pos, , drop = FALSE])^2)
  if (!is.null(jc$Vrest)) {
    Vr <- jc$Vrest[pos, pos, drop = FALSE]
    v <- v + rowSums((G %*% Vr) * G)
  }
  v
}

#' `G V G'` over positions `pos` of `jc`, the whole matrix; see
#' jc_rowvar().
#'
#' @noRd
jc_quad <- function(jc, G, pos) {
  if (is.null(jc$root)) {
    return(G %*% jc$V[pos, pos, drop = FALSE] %*% t(G))
  }
  out <- tcrossprod(G %*% jc$root[pos, , drop = FALSE])
  if (!is.null(jc$Vrest)) {
    out <- out + G %*% jc$Vrest[pos, pos, drop = FALSE] %*% t(G)
  }
  out
}

#' The condition number, at unit diagonal, above which a product
#' `a' V a` formed from the whole covariance in double precision keeps
#' fewer than about four digits: kappa times epsilon is 2e-4 here. On
#' the raw polynomials of dev/illse-floor.R the loss was 8e-5 at
#' degree 6 (kappa 1.6e13), 0.2 percent at degree 7 and a factor of 2.8
#' at degree 8.
#'
#' @noRd
se_product_kappa <- 1e12

#' Is the covariance of `fit` held as a factor whose whole matrix cannot
#' be multiplied without losing its digits (se_reparam_cov())?
#'
#' @noRd
cov_products_lossy <- function(fit) {
  root <- tryCatch(sdr_of(fit)$cov_fixed_root, error = function(e) NULL)
  isTRUE(attr(root, "lossy"))
}

#' Say once per fit and call (`what`), at a call that hands out the
#' whole covariance, that a product formed from it loses its digits,
#' which routes form such products from the factor and are right, and
#' which read the whole matrix and are not.
#'
#' @noRd
warn_cov_products <- function(fit, what) {
  if (!cov_products_lossy(fit)) return(invisible(FALSE))
  cache <- fit$cache
  seen <- if (is.environment(cache)) cache$warned_cov_products
  if (what %in% seen) return(invisible(TRUE))
  if (is.environment(cache)) cache$warned_cov_products <- c(seen, what)
  frm_warning(
    what, ": the fixed design of this fit is ill-conditioned. The ",
    "covariance is computed in coordinates where it is not, so the ",
    "entries of this matrix are accurate. But a product such as a' V a ",
    "formed from the whole matrix in double precision is not accurate: ",
    "it can lose every digit of a fitted value's variance. ",
    cov_routes_right(), ". These read the whole matrix: a product you ",
    "form yourself, and car and lmtest, are not right; ",
    "cooks.distance() and the D1 and D2 tests of anova() on ",
    "frm_multiple() fits invert it and can stop with an error; ",
    "marginaleffects and insight::get_varcov() get NaN. Centering the ",
    "covariate, or poly() for its powers, removes the ill-conditioning",
    call. = FALSE, class = "frmtmb_cov_products")
  invisible(TRUE)
}
