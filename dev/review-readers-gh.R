# Review of lane readers: accuracy of the product Gauss-Hermite rule
# (me_gh_moments() at me_gh_size[r] nodes per direction) at 1 to 4
# directions, with a large sd on a log link.
# Part 1: the rule on integrands with a closed form (exact references).
# Part 2: real fits, references by 1-D and 2-D stats::integrate() over the
# fit's own draw (the sources are independent, so the mean factors), and a
# Monte Carlo check at seed 31, 4e6 draws per row.
source("dev/review-readers-setup.R")
gh <- ns$me_gh_moments
K <- ns$me_gh_size
cat("me_gh_size:", K, "\n\n== Part 1: exp(s z1 + 0.5 (z2 + ... + zr))\n")
cat(sprintf("%2s %4s %4s %12s %12s\n", "r", "K", "s", "rel.mean", "rel.sd"))
for (r in 1:4) {
  for (s in c(1, 1.5, 2, 2.5, 3)) {
    v <- s^2 + 0.25 * (r - 1)
    f <- function(z) exp(s * z[[1L]] + 0.5 * Reduce(`+`, z[-1L], 0))
    m <- gh(f, r, K[r])
    ex1 <- exp(v / 2)
    exsd <- sqrt(exp(2 * v) - exp(v))
    cat(sprintf("%2d %4d %4.1f %12.3e %12.3e %s\n", r, K[r], s,
                m$m1 / ex1 - 1, sqrt(m$m2 - m$m1^2) / exsd - 1,
                if (isTRUE(m$inf)) "inf" else ""))
  }
}
cat("\n== Part 1b: plogis(eta + s z1) at r directions (others flat)\n")
for (r in 1:4) {
  for (s in c(1.5, 3)) {
    for (eta in c(-4, 0)) {
      ref <- stats::integrate(function(z) plogis(eta + s * z) * dnorm(z),
                              -Inf, Inf, rel.tol = 1e-12)$value
      ref2 <- stats::integrate(function(z) plogis(eta + s * z)^2 * dnorm(z),
                               -Inf, Inf, rel.tol = 1e-12)$value
      m <- gh(function(z) plogis(eta + s * z[[1L]]), r, K[r])
      cat(sprintf("r %d K %2d s %.1f eta %3.0f rel.mean %10.3e rel.sd %10.3e\n",
                  r, K[r], s, eta, m$m1 / ref - 1,
                  sqrt(m$m2 - m$m1^2) / sqrt(ref2 - ref^2) - 1))
    }
  }
}

cat("\n== Part 2a: poisson new level, data sd 2.5, seed 6101\n")
set.seed(6101)
J <- 40; n <- 800
d <- data.frame(x = runif(n, -1, 1), g = factor(sample(J, n, TRUE)))
d$y <- rpois(n, exp(-1 + 0.3 * d$x + rnorm(J, 0, 2.5)[d$g]))
fa <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = d)
s <- VarCorr(fa)$g$sd[1, 1]
nd <- data.frame(x = c(-0.5, 0.5), g = "new")
eta <- as.numeric(frm_linpred(fa, newdata = nd, re_formula = NA))
ff <- fitted(fa, newdata = nd, allow_new_levels = TRUE)
cat("sd", s, " rel.mean", format(ff[, "Estimate"] / exp(eta + s^2 / 2) - 1,
                                 digits = 3),
    " Est.Error / u-only sd", format(ff[, "Est.Error"] /
      sqrt((exp(s^2) - 1) * exp(2 * eta + s^2)), digits = 4), "\n")

dm <- function(f, nd, dp) as.numeric(frm_linpred(f, newdata = nd, dpar = dp))
slope <- function(f, nd, dp, v) {
  nd1 <- nd
  nd1[[v]] <- nd1[[v]] + 1
  dm(f, nd1, dp) - dm(f, nd, dp)
}
i1 <- function(fn) stats::integrate(function(z) fn(z) * dnorm(z), -14, 14,
                                    rel.tol = 1e-12, subdivisions = 1000L)$value

cat("\n== Part 2b: hurdle_poisson, me() on mu and hu: 2 directions, seed 6102\n")
set.seed(6102)
n <- 1000
d <- data.frame(l1 = rnorm(n), l2 = rnorm(n), s1 = 0.3, s2 = 0.3)
d$x1 <- d$l1 + rnorm(n, 0, 0.3)
d$x2 <- d$l2 + rnorm(n, 0, 0.3)
lam <- exp(0.5 + 1.0 * d$l1)
cnt <- rpois(n, lam)
while (any(cnt == 0)) cnt[cnt == 0] <- rpois(sum(cnt == 0), lam[cnt == 0])
d$y <- ifelse(runif(n) < plogis(-0.3 + 0.8 * d$l2), 0L, cnt)
fb <- frm(bf(y ~ me(x1, s1), hu ~ me(x2, s2)), family = hurdle_poisson(),
          data = d)
for (s1 in c(1, 2, 2.5, 3)) {
  nd <- data.frame(x1 = 0, x2 = 0.5, s1 = s1, s2 = 1)
  em <- dm(fb, nd, "mu"); bm <- slope(fb, nd, "mu", "x1")
  eh <- dm(fb, nd, "hu"); bh <- slope(fb, nd, "hu", "x2")
  sm <- bm * s1; sh <- bh * 1
  a <- i1(function(z) 1 - plogis(eh + sh * z))
  a2 <- i1(function(z) (1 - plogis(eh + sh * z))^2)
  b <- i1(function(z) { l <- exp(em + sm * z); l / (-expm1(-l)) })
  b2 <- i1(function(z) { l <- exp(em + sm * z); (l / (-expm1(-l)))^2 })
  ref <- a * b
  refsd <- sqrt(a2 * b2 - ref^2)
  mm <- ns$fitted_me_moments(fb, nd, NULL, NULL, NULL, FALSE)
  cat(sprintf("link sd on mu %.2f, on hu %.2f: rel.mean %10.3e rel.sd %10.3e\n",
              sm, sh, mm$mean / ref - 1, sqrt(mm$var) / refsd - 1))
}

cat("\n== Part 2c: hurdle_negbinomial, me() on mu, shape, hu: 3 dirs, seed 6103\n")
set.seed(6103)
n <- 1500
d <- data.frame(l1 = rnorm(n), l2 = rnorm(n), l3 = rnorm(n), s1 = 0.3,
                s2 = 0.3, s3 = 0.3)
d$x1 <- d$l1 + rnorm(n, 0, 0.3)
d$x2 <- d$l2 + rnorm(n, 0, 0.3)
d$x3 <- d$l3 + rnorm(n, 0, 0.3)
mu <- exp(0.7 + 1.0 * d$l1)
sh <- exp(1 + 0.4 * d$l2)
cnt <- rnbinom(n, mu = mu, size = sh)
while (any(cnt == 0)) {
  i <- cnt == 0
  cnt[i] <- rnbinom(sum(i), mu = mu[i], size = sh[i])
}
d$y <- ifelse(runif(n) < plogis(-0.3 + 0.8 * d$l3), 0L, cnt)
fc <- frm(bf(y ~ me(x1, s1), shape ~ me(x2, s2), hu ~ me(x3, s3)),
          family = hurdle_negbinomial(), data = d)
for (s1 in c(1, 2, 2.5, 3)) {
  nd <- data.frame(x1 = 0, x2 = 0, x3 = 0.5, s1 = s1, s2 = 2, s3 = 1)
  em <- dm(fc, nd, "mu"); bm <- slope(fc, nd, "mu", "x1")
  es <- dm(fc, nd, "shape"); bs <- slope(fc, nd, "shape", "x2")
  eh <- dm(fc, nd, "hu"); bh <- slope(fc, nd, "hu", "x3")
  sm <- bm * s1; ss <- bs * 2; shu <- bh * 1
  a <- i1(function(z) 1 - plogis(eh + shu * z))
  a2 <- i1(function(z) (1 - plogis(eh + shu * z))^2)
  inner <- function(p) function(z1) vapply(z1, function(zz) {
    m <- exp(em + sm * zz)
    i1(function(z2) {
      k <- exp(es + ss * z2)
      (m / (-expm1(-k * log1p(m / k))))^p
    })
  }, 0)
  b <- i1(inner(1)); b2 <- i1(inner(2))
  ref <- a * b
  refsd <- sqrt(a2 * b2 - ref^2)
  mm <- ns$fitted_me_moments(fc, nd, NULL, NULL, NULL, FALSE)
  cat(sprintf(paste("link sd mu %.2f shape %.2f hu %.2f: rel.mean %10.3e",
                    "rel.sd %10.3e\n"), sm, ss, shu, mm$mean / ref - 1,
              sqrt(mm$var) / refsd - 1))
}
# a Monte Carlo check of the 3-D reference at one row, seed 31, 4e6 draws
set.seed(31)
nd <- data.frame(x1 = 0, x2 = 0, x3 = 0.5, s1 = 1, s2 = 2, s3 = 1)
em <- dm(fc, nd, "mu"); bm <- slope(fc, nd, "mu", "x1")
es <- dm(fc, nd, "shape"); bs <- slope(fc, nd, "shape", "x2")
eh <- dm(fc, nd, "hu"); bh <- slope(fc, nd, "hu", "x3")
N <- 4e6
m <- exp(em + bm * 1 * rnorm(N)); k <- exp(es + bs * 2 * rnorm(N))
v <- (1 - plogis(eh + bh * rnorm(N))) * m / (-expm1(-k * log1p(m / k)))
mm <- ns$fitted_me_moments(fc, nd, NULL, NULL, NULL, FALSE)
cat("MC row s1 = 1: mean", mean(v), "+-", sd(v) / sqrt(N), " fitted",
    mm$mean, " sd", sd(v), " fitted sd", sqrt(mm$var), "\n")
