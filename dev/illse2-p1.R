# Lane illse2, item 1: cooks.distance() and the D1 and D2 tests of
# anova() on frm_multiple() fits, against exact references at each
# fit's own point (256-bit Hessian, inverted in 256 bits).
#   ILLSE2_ARM=lane|base Rscript dev/illse2-p1.R <re 0|1> <deg> <seeds>
# Raw polynomial of degree `deg` on [1, 2], n 200, 20 groups of 10.
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
re <- as.integer(a[1]) == 1L
deg <- as.integer(a[2])
seeds <- eval(parse(text = a[3]))
xn <- paste0("x", seq_len(deg))
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                              if (re) "+ (1 | g)"))
fo0 <- stats::as.formula(paste("~", paste(xn, collapse = " + ")))
mk <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + if (re) stats::rnorm(20, 0, 0.5)[d$g] else 0
  d$y <- d$y + stats::rnorm(200, 0, 0.3)
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  d
}
# exact covariance of the outer parameters at a fit's own point
exact_V <- function(fit, d) {
  p <- fit$opt$par
  X <- stats::model.matrix(fo0, d)
  stopifnot(max(abs(X - as.matrix(fit$frame$linpreds[[1]]$X))) == 0)
  H <- if (re) mp_hessian(lmm_nll(X, d$y, d$g), p) else
    gauss_hess(X, d$y, p)
  mp_inv(H)
}
cat(sprintf("design: degree %d, %s\n", deg, if (re) "(1 | g)" else "no RE"))
for (seed in seeds) {
  d <- mk(seed)
  r <- cap(frm(fo, data = d))
  if (is_err(r$v)) { cat("seed", seed, "fit ERROR", r$v, "\n"); next }
  fit <- r$v
  pn <- names(fit$opt$par)
  posv <- which(pn %in% c("beta", "betad"))
  Ve <- exact_V(fit, d)
  cat(sprintf("seed %d | code %d | fit warn %s | lost %s | root %s\n",
              seed, fit$opt$convergence, wtag(r$w), lost_of(fit),
              !is.null(ns$sdr_of(fit)$cov_fixed_root)))
  # cooks.distance() over the groups
  inf <- cap(stats::influence(fit, groups = "g", data = d))
  if (!is_err(inf$v)) {
    Pe <- mp_inv(Ve[posv, posv])
    dd <- sweep(inf$v$fixed, 2, inf$v$fixed_full)
    cde <- vapply(seq_len(nrow(dd)), function(i) mp_quad(Pe, dd[i, ]),
                  0) / ncol(dd)
    cd <- cap(stats::cooks.distance(inf$v))
    cat(sprintf("  cooks       %s | warn %s | exact range %.3g-%.3g\n",
                if (is_err(cd$v)) paste("ERROR", substr(cd$v, 1, 70)) else
                  ratio_str(cd$v / cde), wtag(c(inf$w, cd$w)),
                min(cde), max(cde)))
  } else cat("  influence ERROR", substr(inf$v, 1, 80), "\n")
  # frm_multiple() over three imputations that differ
  dl <- lapply(1:3, function(j) {
    set.seed(1000 * seed + j)
    dj <- d
    dj$y <- dj$y + stats::rnorm(nrow(dj), 0, 0.05)
    dj
  })
  mm <- cap(frm_multiple(fo, data = dl))
  if (is_err(mm$v)) { cat("  multiple ERROR", mm$v, "\n"); next }
  fits <- mm$v$fits
  Ves <- lapply(1:3, function(j) exact_V(fits[[j]], dl[[j]]))
  cn <- ns$estimated_coef_names(fits[[1]])
  for (set in list(top2 = xn[deg - 1:0], all = xn)) {
    k <- length(set)
    lab <- if (k == 2L) "top2" else "all"
    pk <- posv[match(set, cn)]
    Q <- vapply(fits, function(f) f$opt$par[pk], numeric(k))
    Ue <- lapply(Ves, function(V) V[pk, pk])
    ub <- (Ue[[1]] + Ue[[2]] + Ue[[3]]) / 3
    ubi <- mp_inv(ub)
    qb <- rowMeans(Q)
    C <- (Q - qb) / sqrt(2)
    tr <- sum(vapply(1:3, function(j) mp_quad(ubi, C[, j]), 0))
    re1 <- max(0, (1 + 1 / 3) * tr / k)
    st1 <- mp_quad(ubi, qb) / ((1 + re1) * k)
    dj <- vapply(1:3, function(j) mp_quad(mp_inv(Ue[[j]]), Q[, j]), 0)
    re2 <- max(0, (1 + 1 / 3) * stats::var(sqrt(dj)))
    st2 <- (mean(dj) / k - (3 + 1) / (3 - 1) * re2) / (1 + re2)
    a1 <- cap(anova(mm$v, method = "D1", constraint = set))
    a2 <- cap(anova(mm$v, method = "D2", use = "wald", constraint = set))
    show <- function(an, st, rv) {
      if (is_err(an$v)) return(paste("ERROR", substr(an$v, 1, 60)))
      sprintf("stat %.9f riv %.9f", an$v$statistic / st, an$v$riv / rv)
    }
    cat(sprintf("  D1 %-4s %s | warn %s | exact stat %.4g riv %.4g\n",
                lab, show(a1, st1, re1), wtag(a1$w), st1, re1))
    cat(sprintf("  D2 %-4s %s | warn %s | exact stat %.4g riv %.4g\n",
                lab, show(a2, st2, re2), wtag(a2$w), st2, re2))
  }
}
