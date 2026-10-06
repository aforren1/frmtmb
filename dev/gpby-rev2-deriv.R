# Reviewer round 2: frm_extra_cov_deriv() and frm_curve_deriv() against
# references written here independently of R/extra-deriv.R (explicit
# kernel derivatives, not Hermite polynomials), on forms the punch did
# not test: an anisotropic 2-D exact gp() differentiated in each
# covariate, an isotropic 2-D one, gp(x, by = f) per sub-GP (a numeric
# by too), an exact gp() beside an unseen level of (1 + x | g) on one
# grid, and difference curves of derivatives across levels and
# positions.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
relmax <- function(a, b) max(abs(as.matrix(a) - b)) / max(abs(b))
# the conditional covariance of the o-th derivative in dimension d of a
# zero-mean SE field with sd^2 s2, length scales ls, fitted positions P
# (rows), nugget 1e-6 on K, at new rows X
ref_deriv <- function(X, P, s2, ls, d, o, nug = 1e-6) {
  kern <- function(A, B) {
    Q <- 0
    for (j in seq_len(ncol(A))) Q <- Q + outer(A[, j], B[, j], "-")^2 /
      (2 * ls[j]^2)
    exp(-Q)
  }
  l <- ls[d]
  D <- outer(X[, d], P[, d], "-")
  g <- kern(X, P) * if (o == 1) -D / l^2 else D^2 / l^4 - 1 / l^2
  DD <- outer(X[, d], X[, d], "-")
  c2 <- kern(X, X) * if (o == 1) 1 / l^2 - DD^2 / l^4 else
    DD^4 / l^8 - 6 * DD^2 / l^6 + 3 / l^4
  K <- kern(P, P) + diag(nug, nrow(P))
  s2 * (c2 - g %*% solve(K, t(g)))
}
blk_par <- function(fit, j = 1L) {
  gis <- fit$frame$linpreds[["y.mu"]]$gps
  bk <- fit$frame$re_blocks[[gis[[j]]$block_id]]
  th <- fit$estimates$theta[bk$theta_idx]
  D <- ncol(gis[[j]]$positions)
  ls <- if (isTRUE(bk$gp_iso)) rep(exp(th[2]), D) else exp(th[-1])
  list(P = gis[[j]]$positions, s2 = exp(2 * th[1]), ls = ls)
}
set.seed(3)
n <- 70
d <- data.frame(x = round(runif(n, 0, 4), 2), z = round(runif(n, 0, 3), 2))
d$y <- sin(1.2 * d$x) + 0.5 * cos(2 * d$z) + rnorm(n, 0, 0.15)
nd <- data.frame(x = c(1.05, 3.7, 4.6, 5.2), z = c(0.4, 2.2, 3.4, 1.5))
for (iso in c(FALSE, TRUE)) {
  fit <- frm(bf(y ~ gp(x, z, iso = iso)), data = d)
  p <- blk_par(fit)
  cat(sprintf("== gp(x, z, iso = %s): ls = %s\n", iso,
              paste(sprintf("%.4f", p$ls), collapse = ", ")))
  for (v in c("x", "z")) for (o in 1:2) {
    E <- frm_extra_cov_deriv(fit, nd, var = v, order = o)
    R <- ref_deriv(as.matrix(nd[, c("x", "z")]), p$P, p$s2, p$ls,
                   match(v, c("x", "z")), o)
    cat(sprintf("  var %s order %d: max rel |E - ref| %.3e | diag %s\n",
                v, o, relmax(E, R),
                paste(sprintf("%.4g", diag(as.matrix(E))), collapse = " ")))
  }
}
# gp(x, by = f): each sub-GP its own theta and positions
set.seed(4)
d2 <- data.frame(x = round(runif(90, 0, 5), 2),
                 f = factor(rep(c("a", "b", "c"), 30)),
                 w = runif(90, 0.5, 2))
d2$y <- ifelse(d2$f == "a", sin(d2$x), ifelse(d2$f == "b", cos(2 * d2$x),
               0.3 * d2$x)) + rnorm(90, 0, 0.2)
fb <- frm(bf(y ~ gp(x, by = f)), data = d2)
nd2 <- data.frame(x = c(2.333, 5.4, 6.1, 2.333, 5.4, 6.1),
                  f = factor(c("a", "a", "b", "b", "c", "c"),
                             levels = c("a", "b", "c")))
for (o in 1:2) {
  E <- as.matrix(frm_extra_cov_deriv(fb, nd2, var = "x", order = o))
  R <- matrix(0, 6, 6)
  for (j in 1:3) {
    r <- which(nd2$f == levels(nd2$f)[j])
    p <- blk_par(fb, j)
    R[r, r] <- ref_deriv(as.matrix(nd2$x[r]), p$P, p$s2, p$ls, 1, o)
  }
  cat(sprintf("== gp(x, by = f) order %d: max rel |E - ref| %.3e | cross-level block max |E| %.1e\n",
              o, relmax(E, R), max(abs(E[nd2$f == "a", nd2$f != "a"]))))
}
fw <- frm(bf(y ~ gp(x, by = w)), data = d2)
nd3 <- data.frame(x = c(2.333, 5.4, 6.1), w = c(0.5, 1, 2))
p <- blk_par(fw)
for (o in 1:2) {
  E <- as.matrix(frm_extra_cov_deriv(fw, nd3, var = "x", order = o))
  R <- ref_deriv(as.matrix(nd3$x), p$P, p$s2, p$ls, 1, o) *
    outer(nd3$w, nd3$w)
  cat(sprintf("== gp(x, by = w) order %d: max rel |E - ref| %.3e\n", o,
              relmax(E, R)))
}
# an exact gp() and an unseen level of (1 + x | g) on one grid
set.seed(6)
ng <- 20
d4 <- data.frame(g = factor(rep(seq_len(ng), each = 8)),
                 x = round(runif(8 * ng, 0, 4), 2))
u0 <- rnorm(ng, 0, 0.6); u1 <- rnorm(ng, 0, 0.3)
d4$y <- sin(d4$x) + u0[d4$g] + u1[d4$g] * d4$x + rnorm(nrow(d4), 0, 0.2)
f4 <- frm(bf(y ~ gp(x) + (1 + x | g)), data = d4)
Sg <- frmtmb::VarCorr(f4)$g$cov[, "Estimate", ]
lev <- c(levels(d4$g), "n1", "n2")
xs <- c(1.111, 4.5, 5.3)
g1 <- data.frame(x = xs, g = factor("n1", levels = lev))
g2 <- data.frame(x = xs, g = factor("n2", levels = lev))
p <- blk_par(f4)
for (o in 1:2) {
  E <- as.matrix(frm_extra_cov_deriv(f4, g1, var = "x", order = o,
                                     allow_new_levels = TRUE))
  R <- ref_deriv(as.matrix(xs), p$P, p$s2, p$ls, 1, o) +
    if (o == 1) Sg[2, 2] else 0
  cat(sprintf("== gp(x) + (1 + x | g), unseen level, order %d: max rel |E - ref| %.3e\n",
              o, relmax(E, R)))
}
# difference curves of derivatives: two unseen levels at the same x; one
# level at two different x grids; the full .se against D V D' + ref
Dref <- function(fit, a, b, e, o, ranl) {
  st <- function(g) {
    lo <- g; hi <- g; lo$x <- g$x - e; hi$x <- g$x + e
    A <- function(h) as.matrix(frm_lp_basis(fit, newdata = h,
                                            allow_new_levels = ranl)$A)
    if (o == 1) (A(hi) - A(lo)) / (2 * e) else
      (A(hi) - 2 * A(g) + A(lo)) / e^2
  }
  lb <- frm_lp_basis(fit, newdata = a, allow_new_levels = ranl)
  D <- st(a) - if (is.null(b)) 0 else st(b)
  D %*% lb$V %*% t(D)
}
for (o in 1:2) {
  e <- if (o == 1) 4.2e-6 else 4.2e-4
  gp_part <- ref_deriv(as.matrix(xs), p$P, p$s2, p$ls, 1, o)
  dv <- frm_curve_deriv(f4, var = "x", order = o, newdata = g1,
                        contrast = g2, re_formula = NULL,
                        allow_new_levels = TRUE, simultaneous = FALSE,
                        eps = e)
  want <- diag(Dref(f4, g1, g2, e, o, TRUE)) +
    2 * (if (o == 1) Sg[2, 2] else 0)
  cat(sprintf(paste0("== deriv difference n1 - n2 (same x), order %d: ",
                     ".se^2 / (D V D' + 2 Sg22) %s\n"), o,
              paste(sprintf("%.8f", dv$.se^2 / want), collapse = " ")))
  ga <- data.frame(x = xs, g = factor("n1", levels = lev))
  gb <- data.frame(x = xs + 0.5, g = factor("n1", levels = lev))
  dv2 <- frm_curve_deriv(f4, var = "x", order = o, newdata = ga,
                         contrast = gb, re_formula = NULL,
                         allow_new_levels = TRUE, simultaneous = FALSE,
                         eps = e)
  Rab <- ref_deriv(as.matrix(c(xs, xs + 0.5)), p$P, p$s2, p$ls, 1, o)
  i1 <- 1:3; i2 <- 4:6
  gpd <- diag(Rab[i1, i1] + Rab[i2, i2] - Rab[i1, i2] - Rab[i2, i1])
  want2 <- diag(Dref(f4, ga, gb, e, o, TRUE)) + gpd
  cat(sprintf(paste0("== deriv difference x vs x + 0.5 (one unseen level),",
                     " order %d: .se^2 / (D V D' + gp diff) %s\n"), o,
              paste(sprintf("%.8f", dv2$.se^2 / want2), collapse = " ")))
}
cat("DONE\n")
