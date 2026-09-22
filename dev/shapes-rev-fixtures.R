# Reviewer fixtures. Sourced by every shapes-rev-* script so that base
# and lane see byte-identical data. Seed and construction live here so
# that a number recorded anywhere in the review has one place to come
# from.

rev_data <- function() {
  set.seed(20260918)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = rnorm(n),
                   f = factor(rep(c("a", "b", "c"), length.out = n)),
                   g = factor(rep(1:15, each = 10)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + c(a = 0, b = 0.4, c = -0.3)[dd$f], 1)
  dd$bin <- rbinom(n, 1, plogis(0.2 + 0.7 * dd$x + 0.5 * dd$z))
  dd$ymix <- rnorm(n, 1 + 0.5 * dd$x + rnorm(15, 0, 0.6)[dd$g], 1)
  dd$ord <- factor(cut(1 + 0.8 * dd$x + rnorm(n),
                       c(-Inf, -0.3, 0.8, Inf), labels = 1:3),
                   ordered = TRUE)
  dd$cnt <- rpois(n, exp(0.5 + 0.4 * dd$x))
  dd
}

rev_fits <- function(dd = rev_data()) {
  list(
    gaussian = frm(bf(y ~ x + f) + gaussian(), data = dd),
    binomial = frm(bf(bin ~ x + z) + bernoulli(), data = dd),
    mixed    = frm(bf(ymix ~ x + (1 | g)) + gaussian(), data = dd),
    ordinal  = frm(bf(ord ~ x) + cumulative(), data = dd),
    distreg  = frm(bf(y ~ x + f, sigma ~ z) + gaussian(), data = dd),
    poisson  = frm(bf(cnt ~ x) + poisson(), data = dd)
  )
}
