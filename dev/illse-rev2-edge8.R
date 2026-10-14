# Reviewer of lane illse, round 2, item 4b: degree-8 raw polynomials with
# (1 | g) on data with no group effect (seeds 1 to 20). Where the sd runs
# to zero, the reference is the exact 256-bit covariance of (beta, log
# sigma) at the fit's own point with the sd at zero (the gaussian
# likelihood of lm(), ML scale): every entry of that block of vcov(full =
# TRUE), and fitted(re_formula = NA) at x = 1.1, 1.5, 1.9.
# Also the cross-curvature between log sd and the rest, to show that
# conditioning on the sd is the full-information answer there.
#   REV_ARM=lane|base Rscript dev/illse-rev2-edge8.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
deg <- 8
xn <- paste0("x", seq_len(deg))
for (seed in 1:20) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  p <- f$opt$par; pn <- names(p)
  pb <- which(pn == "beta"); ps <- which(pn == "betad"); pt <- which(pn == "theta")
  sd <- exp(p[pt])
  l <- ns$sdr_of(f)$se_lost
  th <- if ("theta_1" %in% names(l)) l[["theta_1"]] else "kept"
  if (sd > 1e-4) {
    cat(sprintf("s%2d | sd %.1e away from zero | theta_1 %s | lost %d | warn %s\n",
                seed, sd, th, length(l), wtag(r$w)))
    next
  }
  X <- as.matrix(f$frame$linpreds[[1]]$X)
  Xm <- mM(X)
  rr <- m(d$y) - Xm %*% mpfr2array(m(p[pb]), dim = c(length(pb), 1L))
  s2 <- exp(2 * m(p[ps]))
  q <- c(pb, ps)
  H <- mM(matrix(0, length(q), length(q)))
  ib <- seq_along(pb); is <- length(q)
  H[ib, ib] <- crossprod(Xm) / s2
  hb <- 2 * crossprod(Xm, rr) / s2
  H[ib, is] <- hb; H[is, ib] <- hb
  H[is, is] <- 2 * sum(rr^2) / s2
  Ve <- mp_inv(H)
  V <- suppressWarnings(vcov(f, full = TRUE))
  e <- cov_err(V[q, q], Ve)
  nd <- data.frame(x = c(1.1, 1.5, 1.9), g = d$g[1])
  for (k in seq_len(deg)) nd[[xn[k]]] <- nd$x^k
  A <- cbind(1, as.matrix(nd[, xn]))
  ve <- mp_rowvar(A, Ve[ib, ib])
  fr <- cap(fitted(f, newdata = nd, re_formula = NA))
  lmf <- stats::lm(stats::as.formula(paste("y ~", paste(xn, collapse = " + "))), data = d)
  lms <- sqrt(diag(stats::vcov(lmf)) * (n - deg - 1) / n)
  short <- as.numeric(stats::logLik(lmf)) - as.numeric(stats::logLik(f))
  # cross-curvature of log sd with the rest, relative (natural units)
  Hn <- tryCatch(stats::optimHess(p, f$obj$fn, f$obj$gr), error = function(e) NULL)
  f$obj$fn(p)
  cross <- if (is.null(Hn)) NA else max(abs(Hn[pt, -pt]) / sqrt(abs(Hn[pt, pt]) * abs(diag(Hn)[-pt])))
  cat(sprintf(paste0("s%2d | sd %.1e | theta_1 %s | lost %d | warn %s | ll short of lm %.1e | ",
                     "(beta, log sigma) vs exact: cor %.1e diag %.1e | coef SE / lm %s | ",
                     "fitted SE / exact %s | H row log sd: diag %.1e, cross/sqrt %.1e\n"),
              seed, sd, th, length(l), wtag(r$w), short, e[["cor"]], e[["diag"]],
              ratio_str(sqrt(diag(V))[pb] / lms),
              if (inherits(fr$v, "rev_err")) "ERR" else ratio_str(fr$v[, "Est.Error"] / sqrt(ve)),
              if (is.null(Hn)) NA else Hn[pt, pt], cross))
}
