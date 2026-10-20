# Review of lane illse2, item 5: reproduce the lane's accuracy of
# cooks.distance() and D1 / D2 (Wald) at degree 8 without random effects
# (dev/illse2-p1.R's data), with an independent reference: the precision
# is the closed-form gaussian Hessian at each fit's point (no inversion
# for Cook's), and the Wald precision of a subset is its Schur
# complement, all in 256 bits.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-p1.R <deg> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
deg <- as.integer(a[1])
seeds <- eval(parse(text = a[2]))
xn <- paste0("x", seq_len(deg))
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + ")))
fo0 <- stats::as.formula(paste("~", paste(xn, collapse = " + ")))
mk <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x)
  d$y <- d$y + stats::rnorm(200, 0, 0.3)
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  d
}
# the precision of (beta, log sigma) at th, closed form, 256 bits
prec <- function(X, y, th) {
  p <- ncol(X)
  Xm <- mM(X)
  r <- m(y) - Xm %*% col1(m(th[seq_len(p)]))
  s2 <- exp(2 * m(th[p + 1L]))
  H <- mM(matrix(0, p + 1L, p + 1L))
  H[seq_len(p), seq_len(p)] <- crossprod(Xm) / s2
  hb <- 2 * crossprod(Xm, r) / s2
  H[seq_len(p), p + 1L] <- hb
  H[p + 1L, seq_len(p)] <- hb
  H[p + 1L, p + 1L] <- 2 * sum(r^2) / s2
  H
}
# d' H d for each row of D, in 256 bits
qf_rows <- function(H, D) {
  vapply(seq_len(nrow(D)), function(i) {
    v <- col1(m(D[i, ]))
    num(sum(v * (H %*% v)))
  }, 0)
}
# Schur complement of H on positions t (the precision of the subset)
schur <- function(H, t) {
  o <- setdiff(seq_len(nrow(H)), t)
  Hoo <- H[o, o]
  # solve Hoo Z = Hot by 256-bit Gauss-Jordan through mp_inv
  H[t, t] - H[t, o] %*% mp_inv(Hoo) %*% H[o, t]
}
for (seed in seeds) {
  d <- mk(seed)
  r <- cap(frm(fo, data = d))
  fit <- r$v
  X <- stats::model.matrix(fo0, d)
  H <- prec(X, d$y, fit$opt$par)
  inf <- cap(stats::influence(fit, groups = "g", data = d))
  dd <- sweep(inf$v$fixed, 2, inf$v$fixed_full)
  cde <- qf_rows(H, dd) / ncol(dd)
  cd <- cap(stats::cooks.distance(inf$v))
  cat(sprintf("seed %d | code %d | cooks: %s | warn %s\n", seed,
              fit$opt$convergence,
              if (is_err(cd$v)) paste("ERROR", substr(cd$v, 1, 60)) else
                sprintf("max |ratio - 1| %.2e", max(abs(cd$v / cde - 1))),
              wtag(c(r$w, inf$w, cd$w))))
  dl <- lapply(1:3, function(j) {
    set.seed(1000 * seed + j)
    dj <- d
    dj$y <- dj$y + stats::rnorm(nrow(dj), 0, 0.05)
    dj
  })
  mm <- cap(frm_multiple(fo, data = dl))
  fits <- mm$v$fits
  cn <- ns$estimated_coef_names(fits[[1]])
  pk <- match(xn, cn)
  k <- length(pk)
  Q <- vapply(fits, function(f) f$opt$par[pk], numeric(k))
  P <- lapply(1:3, function(j) schur(prec(X, dl[[j]]$y, fits[[j]]$opt$par),
                                     pk))
  # D2 (Wald): per-imputation chi-squares through each precision
  dj <- vapply(1:3, function(j) {
    v <- col1(m(Q[, j])); num(sum(v * (P[[j]] %*% v)))
  }, 0)
  re2 <- max(0, (1 + 1 / 3) * stats::var(sqrt(dj)))
  st2 <- (mean(dj) / k - (3 + 1) / (3 - 1) * re2) / (1 + re2)
  # D1: ubar = mean of the covariances, which needs the inverses
  ub <- (mp_inv(P[[1]]) + mp_inv(P[[2]]) + mp_inv(P[[3]])) / 3
  ubi <- mp_inv(ub)
  qb <- rowMeans(Q)
  C <- (Q - qb) / sqrt(2)
  tr <- sum(vapply(1:3, function(j) {
    v <- col1(m(C[, j])); num(sum(v * (ubi %*% v)))
  }, 0))
  re1 <- max(0, (1 + 1 / 3) * tr / k)
  v <- col1(m(qb))
  st1 <- num(sum(v * (ubi %*% v))) / ((1 + re1) * k)
  a1 <- cap(anova(mm$v, method = "D1", constraint = xn))
  a2 <- cap(anova(mm$v, method = "D2", use = "wald", constraint = xn))
  sh <- function(an, st) {
    if (is_err(an$v)) paste("ERROR", substr(an$v, 1, 50)) else
      sprintf("|stat ratio - 1| %.2e", abs(an$v$statistic / st - 1))
  }
  cat(sprintf("        D1 all %s | D2 all %s | warn %s\n", sh(a1, st1),
              sh(a2, st2), wtag(c(a1$w, a2$w))))
}
