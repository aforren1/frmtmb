# Lane gpfield: the first-order Laplace correction of a fit with an exact
# gp(), prototype (no frmtmb change). For the linear predictor u that holds
# the field (the mean, or a nonlinear body's parameter), each row's
# derivatives of -log density in u by a 5-point stencil of the family's
# own density (row_lpdf(), shifting the predictor's offset), the latent
# covariance v = A H^-1 A' with A = Z dc/db (rr_jacobians()) and H the inner
# Hessian at the mode, and
#   tk = -1/8 sum(h4 v_ii^2) + 1/8 sum(h3_i h3_j v_ii v_jj v_ij)
#        + 1/12 sum(h3_i h3_j v_ij^3)
# (latres_tk()'s formula). The marginal is about laplace + tk.
gf_tk <- function(fit) {
  ns <- asNamespace("frmtmb")
  fr <- fit$frame
  bks <- Filter(function(b) identical(b$covstruct, "gp"), fr$re_blocks)
  if (!length(bks)) return(NA_real_)
  key <- ns$gp_block_lp(bks[[1]], fr)
  lp <- fr$linpreds[[key]]
  r <- lp$resp
  # every random effect must sit in that predictor
  Z <- ns$lp_Z(fit, lp)
  fam <- fit$spec$responses[[r]]$family
  y <- fr$y[[r]]
  av <- fr$aterm_values[[r]]
  ex <- ns$fit_extras(fit)
  off0 <- lp$offset %||% 0
  h <- 1e-2
  ll <- function(k) {
    f2 <- fit
    f2$frame$linpreds[[key]]$offset <- off0 + k * h
    dp <- ns$eval_dpars(f2)[[r]]
    as.numeric(ns$row_lpdf(fam, y, y, dp, av, ex))
  }
  l <- lapply(-2:2, ll)
  h3 <- -(l[[5]] - 2 * l[[4]] + 2 * l[[2]] - l[[1]]) / (2 * h^3)
  h4 <- -(l[[5]] - 4 * l[[4]] + 6 * l[[3]] - 4 * l[[2]] + l[[1]]) / h^4
  Jb <- ns$rr_jacobians(fit)$Jb
  A <- as.matrix(Z %*% Jb)
  env <- fit$obj$env
  fit$obj$fn(fit$opt$par)
  H <- as.matrix(env$spHess(env$last.par, random = TRUE))
  H <- 0.5 * (H + t(H))
  v <- A %*% solve(H, t(A))
  dv <- diag(v)
  c(tk = -sum(h4 * dv^2) / 8 + sum(outer(h3 * dv, h3 * dv) * v) / 8 +
      sum(outer(h3, h3) * v^3) / 12,
    t4 = -sum(h4 * dv^2) / 8, vmax = max(dv))
}
