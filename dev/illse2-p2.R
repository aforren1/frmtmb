# Lane illse2, item 2: degree 7, seed 4 of the 80 sd-at-zero fits of
# dev/illse-rev-edge4.R. Is the fit short? frmtmb against lme4 (ML) and
# against an exact 256-bit profile of the log-likelihood in the sd.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-p2.R [deg] [seed]
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
deg <- if (length(a)) as.integer(a[1]) else 7L
seed <- if (length(a) > 1) as.integer(a[2]) else 4L
set.seed(seed)
n <- 300
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(1:20, length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
xn <- paste0("x", seq_len(deg))
for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                              "+ (1 | g)"))
X <- stats::model.matrix(stats::as.formula(paste("~", paste(xn,
                                                 collapse = " + "))), d)
r <- cap(frm(fo, data = d))
f <- r$v
p <- f$opt$par
cat(sprintf("frmtmb: code %d | logLik %.10f | sd %.4g | sigma %.6g\n",
            f$opt$convergence, as.numeric(logLik(f)),
            exp(p[names(p) == "theta"]), exp(p[names(p) == "betad"])))
cat("  lost:", lost_of(f), "\n")
for (w in r$w) cat("  WARN:", substr(w, 1, 400), "\n")
V <- suppressWarnings(vcov(f))
cat("  coefficient SEs:", paste(signif(sqrt(diag(V)), 4), collapse = " "),
    "\n")
cat("  se_sd_gain_up:", format(ns$se_sd_gain_up(f, 1L)), "\n")
# lme4, ML; on a raw polynomial lme4 rescales and may warn
l <- cap(lme4::lmer(fo, data = d, REML = FALSE))
if (!is_err(l$v)) {
  cat(sprintf("lme4: logLik %.10f | sd %.4g | sigma %.6g | warn %s\n",
              as.numeric(logLik(l$v)),
              attr(lme4::VarCorr(l$v)$g, "stddev"), sigma(l$v),
              wtag(l$w)))
} else cat("lme4 ERROR", l$v, "\n")
# exact profile: for tau = sd / sigma, beta is GLS and sigma^2 = Q / n,
# in 256 bits; the per-group inverse is I - c J with c = t2 / (1 + n_g t2)
gi <- split(seq_len(n), d$g)
Xm <- mM(X); ym <- m(d$y)
XX <- crossprod(Xm); Xy <- crossprod(Xm, ym); yy <- sum(ym^2)
Sg <- lapply(gi, function(ix) list(s = colSums(Xm[ix, , drop = FALSE]),
                                   sy = sum(ym[ix]), n = length(ix)))
prof <- function(tau) {
  t2 <- m(tau)^2
  A <- XX; bv <- Xy; c0 <- yy; ld <- m(0)
  for (st in Sg) {
    cg <- t2 / (1 + st$n * t2)
    sm <- col1(st$s)
    A <- A - cg * (sm %*% t(sm))
    bv <- bv - cg * st$sy * sm
    c0 <- c0 - cg * st$sy^2
    ld <- ld + log(1 + st$n * t2)
  }
  bh <- mp_inv(A) %*% bv
  Q <- c0 - sum(bv * bh)
  s2 <- Q / n
  # log-likelihood
  -(n / 2) * log(2 * Rmpfr::Const("pi", PREC) * s2) - n / 2 - ld / 2
}
taus <- c(0, 1e-4, 1e-3, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2)
ll <- vapply(taus, function(t) num(prof(t)), 0)
for (i in seq_along(taus)) {
  cat(sprintf("profile tau %-7g logLik %.10f\n", taus[i], ll[i]))
}
opt <- stats::optimize(function(t) -num(prof(t)), c(0, 0.3), tol = 1e-8)
cat(sprintf("profile max: tau %.6g logLik %.10f (frmtmb %+.3e)\n",
            opt$minimum, -opt$objective,
            as.numeric(logLik(f)) + opt$objective))
# the covariance at the fit's own point against the exact one
H <- mp_hessian(lmm_nll(X, d$y, d$g), p)
Ve <- mp_inv(H)
Vf <- suppressWarnings(vcov(f, full = TRUE))
if (all(is.finite(Vf))) {
  dg <- sqrt(num(diag(Ve)))
  E <- abs(num(mM(Vf) - Ve)) / outer(dg, dg)
  cat(sprintf("vcov(full) against exact: cor %.2e | diag %.2e | root %s\n",
              max(E), max(abs(diag(Vf) / num(diag(Ve)) - 1)),
              !is.null(ns$sdr_of(f)$cov_fixed_root)))
  cat("  exact SE(log sd)", signif(dg[length(dg)], 4), "frmtmb",
      signif(sqrt(Vf[nrow(Vf), nrow(Vf)]), 4), "\n")
} else cat("vcov(full) not finite; exact SE(log sd)",
           signif(sqrt(num(diag(Ve)))[length(p)], 4), "\n")
