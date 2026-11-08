# Designs with a scaled condition number of 1e4 or more that the base
# build may already fit right: each arm against the same model on the
# centered column (or orthogonal polynomial). REV_ARM selects the build.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
kap <- function(X) {
  X <- as.matrix(X)
  R <- qr.R(qr(X / rep(sqrt(colSums(X^2) / nrow(X)), each = nrow(X))))
  1 / rcond(R, triangular = TRUE)
}
cat(sprintf("%-34s %8s %6s %5s %12s %10s %s\n", "case", "kappa", "basis", "conv",
            "dll", "maxshift", "warn"))
row <- function(lbl, X, rw, fz, map) {
  if (rw$err) { cat(lbl, "ERROR", conditionMessage(rw$f), "\n"); return() }
  bw <- quiet(fixef(rw$f)); bz <- quiet(fixef(fz))
  sh <- if (is.null(map)) NA else max(abs(map(bw[, 1]) - bz[, 1]) / bz[, 2])
  cat(sprintf("%-34s %8.2g %6s %5s %+12.3e %10.3g %s\n", lbl, kap(X),
              basis_on(rw$f), rw$f$opt$convergence, ll_of(rw$f) - ll_of(fz),
              sh, short_w(rw)))
}
for (s in 1:5) {
  set.seed(800 + s)
  n <- 200
  z <- rnorm(n)
  g <- factor(rep(1:10, length.out = n))
  u <- rnorm(10, 0, 0.5)[g]
  for (m in c(1e4, 1e5)) {
    d <- data.frame(z = z, w = z + m, g = g,
                    yg = 0.5 + 0.4 * z + rnorm(n),
                    yp = rpois(n, exp(0.5 + 0.4 * z)),
                    yb = rbinom(n, 1, plogis(0.3 + 0.6 * z)),
                    ygr = 0.5 + 0.4 * z + u + rnorm(n))
    map <- function(b) c(b[1] + b[2] * m, b[2])
    for (cc in list(c("yg", "gaussian"), c("yp", "poisson"),
                    c("yb", "bernoulli"))) {
      fam <- switch(cc[2], gaussian = gaussian(), poisson = poisson(),
                    bernoulli = bernoulli())
      fz <- quiet(frm(bf(as.formula(paste(cc[1], "~ z"))), family = fam, data = d))
      rw <- rfit(bf(as.formula(paste(cc[1], "~ w"))), family = fam, data = d)
      row(sprintf("%s ~ w, m=%g, seed %d", cc[1], m, 800 + s), cbind(1, d$w),
          rw, fz, map)
    }
    fz <- quiet(frm(bf(ygr ~ z + (1 | g)), family = gaussian(), data = d))
    rw <- rfit(bf(ygr ~ w + (1 | g)), family = gaussian(), data = d)
    row(sprintf("lmm ygr ~ w, m=%g, seed %d", m, 800 + s), cbind(1, d$w),
        rw, fz, map)
  }
  # raw polynomials on [1, 2], against the orthogonal spelling (ll only)
  x <- runif(n, 1, 2)
  d <- data.frame(x = x, g = g, y = sin(3 * x) + u + rnorm(n, 0, 0.2))
  for (deg in 3:6) {
    fr <- as.formula(sprintf("y ~ poly(x, %d, raw = TRUE)", deg))
    fo <- as.formula(sprintf("y ~ poly(x, %d)", deg))
    fz <- quiet(frm(bf(fo), family = gaussian(), data = d))
    rw <- rfit(bf(fr), family = gaussian(), data = d)
    row(sprintf("raw poly %d, seed %d", deg, 800 + s),
        model.matrix(fr, d), rw, fz, NULL)
    fz <- quiet(frm(bf(update(fo, . ~ . + (1 | g))), family = gaussian(), data = d))
    rw <- rfit(bf(update(fr, . ~ . + (1 | g))), family = gaussian(), data = d)
    row(sprintf("raw poly %d + (1|g), seed %d", deg, 800 + s),
        model.matrix(fr, d), rw, fz, NULL)
  }
}
