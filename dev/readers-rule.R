# Part 2c of dev/review-readers-gh.R on this lane's build, with the warnings
# fitted_me_moments() gives (the rule check at three directions). Seed 6103.
source("dev/review-readers-setup.R")
dm <- function(f, nd, dp) as.numeric(frm_linpred(f, newdata = nd, dpar = dp))
slope <- function(f, nd, dp, v) {
  nd1 <- nd
  nd1[[v]] <- nd1[[v]] + 1
  dm(f, nd1, dp) - dm(f, nd, dp)
}
i1 <- function(fn) stats::integrate(function(z) fn(z) * dnorm(z), -14, 14,
                                    rel.tol = 1e-12, subdivisions = 1000L)$value
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
  r <- with_warn(ns$fitted_me_moments(fc, nd, NULL, NULL, NULL, FALSE))
  mm <- r$value
  cat(sprintf(paste("link sd mu %.2f shape %.2f hu %.2f: rel.mean %10.3e",
                    "rel.sd %10.3e warn %s\n"), sm, ss, shu,
              mm$mean / ref - 1,
              sqrt(mm$var) / refsd - 1, paste(wclass(r$warnings), collapse = ",")))
}
