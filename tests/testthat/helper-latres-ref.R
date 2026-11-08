# An independent reference for the standard errors of a fit with a latent
# residual: the joint covariance formed densely here from joint_precision()
# alone, so a test of the sparse route (R/joint-sparse.R) does not compare it
# with itself. Where the standard-error check kept a propagation covariance,
# the reference is joint_cov_repair()'s block formula with it, written out
# again here; otherwise it is solve(Q).

latres_ref_V <- function(fit) {
  Q <- as.matrix(joint_precision(fit))
  sdr <- sdr_of(fit)
  r <- fit$obj$env$random
  f <- setdiff(seq_len(nrow(Q)), r)
  Vp <- sdr$cov_fixed_prop
  if (is.null(Vp) || nrow(Vp) != length(f)) return(solve(Q))
  W <- solve(Q[r, r])
  J <- -W %*% Q[r, f]
  V <- matrix(0, nrow(Q), nrow(Q))
  V[r, r] <- W + J %*% Vp %*% t(J)
  V[r, f] <- J %*% Vp
  V[f, r] <- t(V[r, f])
  V[f, f] <- Vp
  V
}

# the link-scale standard error of each row of frm_lp_basis()'s design
latres_ref_se <- function(fit, newdata = NULL, ...) {
  lb <- frm_lp_basis(fit, newdata = newdata, ...)
  V <- latres_ref_V(fit)
  A <- as.matrix(lb$A)
  pos <- lb$coef_pos
  sqrt(pmax(rowSums((A %*% V[pos, pos]) * A), 0) + lb$extra_var)
}

# max |a / b - 1| over the rows finite on both, after checking that the same
# rows are finite
latres_ref_ratio <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  testthat::expect_identical(is.finite(a), is.finite(b))
  ok <- is.finite(b)
  if (!any(ok)) return(0)
  max(abs(a[ok] / b[ok] - 1))
}
