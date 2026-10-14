# Reviewer of lane illse, attack 4: the three thresholds of the new route
# (natural kappa 1e8: route or tiers; natural kappa 1e12: lossy, which
# decides vcov()'s warning and marginaleffects' NaN; first-pass kappa
# 1e4 in the design's coordinates: refine or not). Raw polynomials of
# degree d on [a, a + 1], a swept so that the natural kappa crosses a
# threshold; gaussian, exact Hessian in 256 bits at the fit's point.
# Per fit: natural kappa, route, lossy, first-pass kappa, refined, the
# largest coefficient SE error against exact, and the refined against
# the unrefined covariance (cor scale).
#   Rscript dev/illse-rev-thresh.R <deg> <a values> [seed]
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
a <- commandArgs(TRUE)
deg <- as.integer(a[1])
avals <- eval(parse(text = a[2]))
seed <- if (length(a) > 2) as.integer(a[3]) else 1L
for (a0 in avals) {
  set.seed(seed)
  d <- data.frame(x = a0 + stats::runif(200))
  d$y <- sin(3 * (d$x - a0)) + stats::rnorm(200, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + ")))
  r <- cap(frm(fo, data = d))
  f <- r$v
  if (inherits(f, "rev_err")) { cat("a", a0, "ERROR", f, "\n"); next }
  p <- f$opt$par; pn <- names(p); pb <- which(pn == "beta"); ps <- which(pn == "betad")
  X <- stats::model.matrix(stats::as.formula(paste("~", paste(xn, collapse = " + "))), d)
  Xf <- as.matrix(f$frame$linpreds[[1]]$X)
  if (!identical(dim(X), dim(Xf))) { cat("a", a0, "frame dropped a column\n"); next }
  Xm <- mM(X)
  rr <- m(d$y) - Xm %*% mpfr2array(m(p[pb]), dim = c(length(pb), 1L))
  s2 <- exp(2 * m(p[ps]))
  H <- mM(matrix(0, length(p), length(p)))
  H[pb, pb] <- crossprod(Xm) / s2
  hb <- 2 * crossprod(Xm, rr) / s2
  H[pb, ps] <- hb; H[ps, pb] <- hb
  H[ps, ps] <- 2 * sum(rr^2) / s2
  Ve <- mp_inv(H)
  Hn <- num(H); D <- sqrt(diag(Hn))
  ev <- eigen(Hn / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  kap <- max(ev) / min(ev)
  # what the fit's own AD Hessian says (the route's gate reads it)
  Ha <- f$obj$he(p); Da <- sqrt(abs(diag(Ha)))
  eva <- eigen((Ha + t(Ha)) / 2 / outer(Da, Da), symmetric = TRUE, only.values = TRUE)$values
  kap_ad <- max(abs(eva)) / min(abs(eva))
  root <- ns$sdr_of(f)$cov_fixed_root
  V <- suppressWarnings(vcov(f, full = TRUE))
  e <- cov_err(V, Ve)
  # first pass in the design's coordinates, and refined vs unrefined
  M <- ns$reparam_basis(f$frame, names(p))
  k1 <- NA; dref <- NA; refined <- NA
  if (!is.null(M)) {
    h1 <- ns$reparam_hessian(f$obj, p, M)
    D1 <- sqrt(diag(h1$H))
    e1 <- eigen(h1$H / outer(D1, D1), symmetric = TRUE, only.values = TRUE)$values
    k1 <- max(e1) / min(e1)
    h2 <- ns$reparam_hessian_refined(f$obj, p, M)
    refined <- !identical(h2$M, M)
    V1 <- M %*% solve(h1$H) %*% t(M)
    V2 <- h2$M %*% solve(h2$H) %*% t(h2$M)
    dg <- sqrt(diag(V2))
    dref <- max(abs(V1 - V2) / outer(dg, dg))
    f$obj$fn(p)
  }
  vw <- cap(vcov(f))
  cat(sprintf(paste0("deg %d a %-7g s%d | kappa exact %.3e AD %.3e | route %s lossy %s | ",
                     "vcov warn %-8s | cov err cor %.1e diag %.1e | first-pass kappa %.2e ",
                     "refined %s | refined vs unrefined %.1e\n"),
              deg, a0, seed, kap, kap_ad, !is.null(root), isTRUE(attr(root, "lossy")),
              wtag(vw$w), e[["cor"]], e[["diag"]], k1, refined, dref))
}
