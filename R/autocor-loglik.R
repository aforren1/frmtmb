# The pointwise log-density of a response under a residual covariance
# term: brms 2.23.0's leave-one-out conditional density of each row
# given the other rows of its group (log_lik_gaussian_time(),
# log_lik_gaussian_fcor(), log_lik_gaussian_lagsar(),
# log_lik_gaussian_errorsar() and their student versions, after the loo
# vignette on non-factorized models).
#
# Everything is read off the precision matrix Q of the response, never
# off the covariance: with e = y - E(y),
#
#   g = Q e,  cbar = diag(Q),  E(y_i | y_-i) = y_i - g_i / cbar_i,
#   Var(y_i | y_-i) = 1 / cbar_i,
#
# and for the multivariate t (Shah et al. 2014, eq. 6, brms's
# student_t_cov_factor()) the conditional is a t with nu + k - 1 degrees
# of freedom and scale^2 (nu + e'Qe - g_i^2 / cbar_i) / (nu + k - 1) /
# cbar_i, k the size of the group. brms forms Q by inverting a dense
# n x n covariance per draw. Here the ar() family inverts one k x k
# correlation per pattern of present time points, shared by every group
# of that pattern, and a registered term gives Q = K'K from the Jacobian
# K of its whitening map, which for sar() is the sparse (I - rho W) /
# sigma: Q itself is never formed.

#' brms's pointwise log-density for one response with a residual
#' covariance term, at the dpar values `dpv` (the whole list
#' `eval_dpars()` returns) and the fit's own `thetaac`. One value per
#' row of the response, in the response's row order.
#'
#' @noRd
autocor_row_loglik <- function(fit, dpv, resp) {
  ac <- fit$frame[["autocor"]][[resp]]
  if (is.null(ac) || autocor_is_cond(ac)) {
    frm_stop("internal: autocor_row_loglik() on response '", resp,
             "', which has no residual covariance term", call. = FALSE)
  }
  y <- as.numeric(fit$frame[["y"]][[resp]])
  n <- length(y)
  dp <- dpv[[resp]]
  mu <- rep_len(as.numeric(dp[["mu"]]), n)
  sigma <- as.numeric(dp[["sigma"]])
  th <- fit$estimates[["thetaac"]][ac[["theta_idx"]]]
  pt <- if (autocor_is_seam(ac)) {
    autocor_seam_loo_parts(ac, y, mu, sigma, th)
  } else {
    autocor_group_loo_parts(ac, y, mu, rep_len(sigma, n), th)
  }
  nu <- if (isTRUE(ac[["student"]])) as.numeric(dp[["nu"]])[1L] else NULL
  autocor_loo_density(pt, nu)
}

#' The leave-one-out pieces of the ar() family, one correlation pattern
#' at a time: `g = Q e`, `cbar = diag(Q)`, each row's group quadratic
#' form `e'Qe` and group size `k`. With `D = diag(sigma)` and `R` the
#' pattern's correlation, `Q = D^-1 R^-1 D^-1`, so on the standardized
#' residuals `z` the pattern needs one `k x k` inverse for all its
#' groups.
#'
#' @noRd
autocor_group_loo_parts <- function(ac, y, mu, sigma, th) {
  n <- length(y)
  R <- as.matrix(autocor_cor(th, ac))
  z <- (y - mu) / sigma
  g <- numeric(n)
  cbar <- numeric(n)
  quad <- numeric(n)
  k <- integer(n)
  for (pt in ac[["patterns"]]) {
    kk <- pt$k
    Rp <- if (kk == ac[["d"]]) R else {
      matrix(as.vector(R)[pt$gather], kk, kk)
    }
    Ri <- chol2inv(chol(Rp))
    Z <- matrix(z[pt$rows], kk, pt$G)
    GZ <- Ri %*% Z
    rows <- pt$rows
    g[rows] <- as.vector(GZ) / sigma[rows]
    cbar[rows] <- rep(diag(Ri), pt$G) / sigma[rows]^2
    quad[rows] <- rep(colSums(Z * GZ), each = kk)
    k[rows] <- kk
  }
  list(g = g, cbar = cbar, quad = quad, k = k)
}

#' The leave-one-out pieces of a registered term, from its whitening map
#' `u = K (y - m)`: `Q = K'K`, so `g = K'u`, `cbar` the squared column
#' norms of `K` and `e'Qe = u'u`, all at the cost of `K`'s nonzeros.
#' `K` comes from the term's `jacobian` slot, or, for a term without
#' one, as the inverse of the `B` that `unwhiten()` of the identity
#' gives: exact, and a dense `n x n` solve.
#'
#' @noRd
autocor_seam_loo_parts <- function(ac, y, mu, sigma, th) {
  n <- length(y)
  e <- ac_block_entry(ac)
  sg <- autocor_seam_sigma(e, sigma, n)
  # a block built over other rows reads its parameters through its
  # theta_map, as autocor_seam_B() and the simulator do; the map is
  # applied here once, and the fallback gets the raw th because
  # autocor_seam_B() applies it itself
  thm <- autocor_seam_theta_on(ac, th)
  u <- as.numeric(e[["whiten"]](y, mu, sg, thm, ac[["aux"]])[["u"]])
  K <- autocor_seam_jacobian(ac, e, sg, th, thm, mu)
  if (inherits(K, "Matrix")) {
    g <- as.numeric(Matrix::crossprod(K, u))
    cbar <- as.numeric(Matrix::colSums(K * K))
  } else {
    g <- as.numeric(crossprod(K, u))
    cbar <- colSums(K * K)
  }
  list(g = g, cbar = cbar, quad = rep(sum(u * u), n), k = rep(n, n))
}

#' The Jacobian `K` of a registered term's whitening map in `y`, at
#' `sigma` as the map receives it.
#'
#' @noRd
autocor_seam_jacobian <- function(ac, e, sg, th, thm, mu) {
  if (!is.null(e[["jacobian"]])) {
    return(e[["jacobian"]](sg, thm, ac[["aux"]]))
  }
  n <- ac[["d"]]
  if (n > autocor_loo_dense_max) {
    frm_warning("The '", ac[["struct"]], "' term of ", ac[["label"]],
                " has no jacobian slot, so its pointwise log-density ",
                "inverts a dense ", n, " x ", n, " matrix per draw, a cost ",
                "that grows as n^3 (?frmtmb_register_autocor, \"The ",
                "Jacobian\", has the measured times). The values are ",
                "exact. The term's package can give ",
                "frmtmb_register_autocor(jacobian =), a sparse matrix ",
                "where the map is sparse", call. = FALSE,
                class = "frmtmb_autocor_dense_fallback")
  }
  B <- autocor_seam_B(ac, th, mu, sg)
  solve(B)
}

# Rows above which a registered term without a jacobian slot warns that
# its pointwise log-density takes the dense route: measured on this
# machine (dev/loglik-fallback-cost.R) at 0.5 s and 48 Mb per draw at
# 1000 rows, 0.56 h for 4000 draws, and 1.9 h at 2000 rows. A warning,
# not a refusal: the values are exact, and 345 Mb per draw at 3000 rows
# is not a memory hazard on a machine that fits such a model.
autocor_loo_dense_max <- 1000L

#' brms's leave-one-out conditional density from the pieces: the normal
#' of `y_i` around `y_i - g_i / cbar_i` with variance `1 / cbar_i`, or
#' for `student()` the t of `student_t_cov_factor()`. Both are written in
#' the residual `g / cbar` itself, which is what `y - yloo` is, so no
#' `y` is subtracted from a nearly equal `yloo`.
#'
#' @noRd
autocor_loo_density <- function(pt, nu = NULL) {
  r <- pt$g / pt$cbar
  if (is.null(nu)) {
    return(-0.5 * log(2 * pi) + 0.5 * log(pt$cbar) - 0.5 * r * pt$g)
  }
  df <- nu + pt$k - 1
  # nu + e'Qe - g^2 / cbar: the group's quadratic form without row i
  s2 <- (nu + pt$quad - r * pt$g) / df / pt$cbar
  stats::dt(r / sqrt(s2), df, log = TRUE) - 0.5 * log(s2)
}
