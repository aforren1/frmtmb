# Lane adhess: y ~ a * exp(b * x), gaussian, without random effects.
# frmtmb's covariance of (a, b, log sigma) against the exact one: the
# observed information at frmtmb's estimate, in closed form and inverted
# in 256 bits. nls() beside it, against its own formula evaluated in 256
# bits at its estimate, (RSS / (n - 2)) (J' J)^-1. The two definitions
# differ (nls drops the residual curvature term and divides by n - 2),
# so the SE ratio frmtmb / nls is shown only as context.
#   ADHESS_ARM=lane|base Rscript dev/adhess-nls.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
nll_hess_mp <- function(x, y, a, b, ls) {
  x <- m(x); y <- m(y); a <- m(a); b <- m(b)
  s2 <- exp(2 * m(ls))
  e <- exp(b * x)
  r <- y - a * e
  fa <- e; fb <- a * x * e; fab <- x * e; fbb <- a * x^2 * e
  H <- mM(matrix(0, 3, 3))
  H[1, 1] <- sum(fa^2) / s2
  H[1, 2] <- H[2, 1] <- sum(fa * fb - r * fab) / s2
  H[2, 2] <- sum(fb^2 - r * fbb) / s2
  H[1, 3] <- H[3, 1] <- -2 * sum(r * fa) / s2
  H[2, 3] <- H[3, 2] <- -2 * sum(r * fb) / s2
  H[3, 3] <- 2 * sum(r^2) / s2
  H
}
nls_cov_mp <- function(x, y, a, b) {
  x <- m(x); y <- m(y); a <- m(a); b <- m(b)
  e <- exp(b * x)
  r <- y - a * e
  J <- mM(matrix(0, length(x), 2))
  J[, 1] <- e
  J[, 2] <- a * x * e
  s2 <- sum(r^2) / (length(x) - 2)
  num(mp_inv(crossprod(J)) * s2)
}
cases <- list(x20 = 20, x01 = 0, x5 = 5)
for (cn in names(cases)) {
  for (seed in 1:3) {
    set.seed(seed)
    x0 <- cases[[cn]]
    d <- data.frame(x = x0 + runif(200))
    d$y <- 2 * exp(0.8 * (d$x - x0)) + rnorm(200, 0, 0.3)
    st <- c(2 * exp(-0.8 * x0), 0.8)
    r <- cap(frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d,
                 start = list(beta = st)))
    f <- r$v
    p <- f$opt$par
    Vf <- suppressWarnings(vcov(f, full = TRUE))
    Vx <- num(mp_inv(nll_hess_mp(d$x, d$y, p[1], p[2], p[3])))
    ef <- cov_err(Vf, Vx)
    g <- nls(y ~ a * exp(b * x), data = d, start = list(a = p[[1]], b = p[[2]]),
             control = nls.control(tol = 1e-8, maxiter = 200,
                                   minFactor = 1e-10, warnOnly = TRUE))
    cg <- coef(g)
    Vn <- nls_cov_mp(d$x, d$y, cg[1], cg[2])
    en <- cov_err(unname(vcov(g)), Vn)
    cat(sprintf(paste0("%-4s seed %d | frmtmb SE err %.2e cor err %.2e | ",
                       "nls SE err %.2e cor err %.2e | ",
                       "SE frmtmb / nls %s | |b - b_nls| / SE %.1e | ",
                       "warn %s\n"),
                cn, seed, ef[1], ef[2], en[1], en[2],
                paste(sprintf("%.4f", sqrt(diag(Vf))[1:2] /
                                sqrt(diag(vcov(g)))), collapse = ","),
                max(abs(p[1:2] - cg) / sqrt(diag(Vx))[1:2]), wtag(r$w)))
  }
}
