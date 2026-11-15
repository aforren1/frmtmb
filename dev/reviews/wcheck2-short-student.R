# Short check of lane wcheck2's fix round: the exact-mode Laplace
# approximation of asym_laplace y ~ x + (1 | g) + (1 | gr(h, dist =
# "student", dist_nu = nu)) with h a copy of g (one t effect per group),
# by hand, along g's log sd. No frmtmb: the same data as
# dev/reviews/wcheck2-short-ald.R part student (seed 20261116).
# l is piecewise linear in (b_g, b_t), so TMB's inner Hessian is
# diag(1 / C, H_t(b_t)) with C = sd_g^2 and H_t = -d2 log t(b_t); the
# gaussian block's constant and curvature cancel, the t block's do not:
#   Laplace_g = max_{b_g, b_t} [l(b_g + b_t) + log t(b_t; s, nu)
#               - b_g^2 / (2 C)] + log(2 pi) / 2 - log(H_t(b_t*)) / 2.
#   Rscript dev/reviews/wcheck2-short-student.R [nu ...]
nus <- as.numeric(commandArgs(TRUE))
if (!length(nus)) nus <- c(2.5, 5)
set.seed(20261116)
G <- 20; n <- 8 * G
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = 8)))
d$h <- factor(paste0("h", d$g))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 1.5)[as.integer(d$g)] +
  rt(nlevels(d$h), 2.5)[as.integer(d$h)] + rnorm(n, 0, 0.5)
ald <- function(y, mu, sigma, p) {
  u <- (y - mu) / sigma
  log(p) + log(1 - p) - log(sigma) - 0.5 * (abs(u) + (2 * p - 1) * u)
}
gi <- split(seq_len(n), d$g)
lap <- function(C, s, nu, b0, b1, sigma, p = 0.5) {
  mu <- b0 + b1 * d$x
  tot <- 0
  for (rows in gi) {
    yr <- d$y[rows]; mr <- mu[rows]
    lt <- function(bt) dt(bt / s, nu, log = TRUE) - log(s)
    # max over m = b_g + b_t of l(m) - (m - b_t)^2 / (2 C), concave in m:
    # on the segment with k residuals below m the slope of l is
    # (n p - k) / sigma, so the max is a segment's stationary point
    # clipped to the segment (exact)
    r <- sort(yr - mr); nr <- length(r)
    lo <- c(-Inf, r); hi <- c(r, Inf)
    sl <- (nr * p - 0:nr) / sigma
    inner1 <- function(bt) {
      m <- pmin(pmax(bt + C * sl, lo), hi)
      max(vapply(m, function(mm) sum(ald(yr, mr + mm, sigma, p)), 0) -
            (m - bt)^2 / (2 * C)) + lt(bt)
    }
    inner <- function(bt) vapply(bt, inner1, 0)
    R <- 25 + 10 * s
    bg <- seq(-R, R, length.out = 401)
    v <- inner(bg)
    k <- which.max(v)
    o <- optimize(inner, c(bg[max(1, k - 1)], bg[min(401, k + 1)]),
                  maximum = TRUE, tol = 1e-10)
    if (o$objective < v[k]) o <- list(maximum = bg[k], objective = v[k])
    z2 <- (o$maximum / s)^2
    Ht <- (nu + 1) / (nu * s^2) * (1 - z2 / nu) / (1 + z2 / nu)^2
    if (Ht <= 0) return(NA_real_)
    tot <- tot + o$objective + 0.5 * log(2 * pi) - 0.5 * log(Ht)
  }
  tot
}
lsd <- c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8)
co <- coef(lm(y ~ x, d))
for (nu in nus) {
  cat(sprintf("--- nu %g, h a copy of g, b0 %.4g, b1 %.4g, sigma 0.3, quantile 0.5 fixed\n",
              nu, co[1], co[2]))
  # profile over the t block's log scale at each C
  prof <- vapply(lsd, function(l) {
    o <- optimize(function(ls) {
      v <- lap(exp(2 * l), exp(ls), nu, co[1], co[2], 0.3)
      if (is.na(v)) -1e10 else v
    }, c(-4, 3), maximum = TRUE, tol = 1e-6)
    c(o$objective, o$maximum)
  }, c(0, 0))
  cat(sprintf("  log sd_g     %s\n  max over s   %s\n  at log s     %s\n",
              paste(format(lsd), collapse = " "),
              paste(format(prof[1, ], digits = 8), collapse = " "),
              paste(format(prof[2, ], digits = 4), collapse = " ")))
  cat(sprintf("  never falls: %s; argmax at log sd_g %g (value %.6g), at +8 %.6g\n",
              all(diff(prof[1, ]) >= -1e-6), lsd[which.max(prof[1, ])],
              max(prof[1, ]), prof[1, length(lsd)]))
  # also with s held at its value at the argmax
  s0 <- exp(prof[2, which.max(prof[1, ])])
  v <- vapply(lsd, function(l) lap(exp(2 * l), s0, nu, co[1], co[2], 0.3), 0)
  cat(sprintf("  s held at %.4g: %s\n", s0,
              paste(format(v, digits = 8), collapse = " ")))
}
cat("DONE\n")
