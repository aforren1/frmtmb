# v0.19: latent classes combined with continuous random effects
# (growth-mixture models). The class sum happens conditional on the
# latent effects, so one Laplace integrates them; validation is
# against the exact closed-form gaussian marginal.

sim_lcm <- function(seed = 23, ng = 50, m = 8) {
  set.seed(seed)
  cls <- rbinom(ng, 1, 0.4)
  bg <- rnorm(ng, 0, 0.5)
  g <- rep(seq_len(ng), each = m)
  y <- rnorm(ng * m, (c(-1, 2)[cls + 1] + bg)[g], 0.8)
  list(dd = data.frame(y = y, g = factor(g)), cls = cls, m = m,
       ng = ng)
}

# exact per-group-per-class gaussian marginal (Sherman-Morrison)
lcm_llk <- function(Ym, mu, tau, sig) {
  m <- nrow(Ym)
  s2 <- sig^2
  t2 <- tau^2
  ld <- m * log(s2) + log1p(m * t2 / s2)
  r <- Ym - mu
  q <- colSums(r^2) / s2 -
    (t2 / (s2 * (s2 + m * t2))) * colSums(r)^2
  -0.5 * (m * log(2 * pi) + ld + q)
}

test_that("latent classes with class-specific intercept REs work", {
  s <- sim_lcm()
  dd <- s$dd
  Ym <- matrix(dd$y, nrow = s$m)
  fit <- frm(bf(y ~ 1 + (1 | g)) +
               mixture(gaussian(), gaussian(), groups = ~g), data = dd)
  expect_length(fit$frame$re_blocks, 2L)   # class-specific effects

  nll <- function(p) {
    pi1 <- stats::plogis(p[7])
    M <- cbind(log(pi1) + lcm_llk(Ym, p[1], exp(p[3]), exp(p[5])),
               log(1 - pi1) + lcm_llk(Ym, p[2], exp(p[4]), exp(p[6])))
    mx <- pmax(M[, 1], M[, 2])
    -sum(mx + log(exp(M[, 1] - mx) + exp(M[, 2] - mx)))
  }
  op <- stats::optim(c(-1, 2, log(0.5), log(0.5), log(0.8), log(0.8),
                       stats::qlogis(0.4)), nll, method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  # Laplace approximates the class-mixture integrand: a small,
  # documented bias against the exact marginal
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 0.5)
  mus <- sort(c(fixef(fit)$mu1[[1]], fixef(fit)$mu2[[1]]))
  expect_equal(mus, sort(op$par[1:2]), tolerance = 0.05)

  # empirical-Bayes classification and simulation both work
  P <- mixture_probs(fit)
  acc <- mean((P[, 1] > 0.5) == (s$cls == 0))
  expect_gt(max(acc, 1 - acc), 0.95)
  sm <- simulate(fit, nsim = 1, re.form = NA, seed = 1)
  expect_equal(stats::sd(tapply(sm[[1]], dd$g, mean)),
               stats::sd(tapply(dd$y, dd$g, mean)), tolerance = 0.25)
})

test_that("quadrature is exact for a univariate class-RE integrand", {
  s <- sim_lcm(seed = 29)
  dd <- s$dd
  Ym <- matrix(dd$y, nrow = s$m)
  # the RE lives in class 1 only: per-group integrand is univariate,
  # so adaptive GK integrates it exactly
  fq <- frm(bf(y ~ 1 + (1 | g), mu2 ~ 1) +
              mixture(gaussian(), gaussian(), groups = ~g),
            data = dd, quadrature = TRUE)
  nll <- function(p) {
    pi1 <- stats::plogis(p[6])
    l2 <- colSums(stats::dnorm(Ym, p[2], exp(p[5]), log = TRUE))
    M <- cbind(log(pi1) + lcm_llk(Ym, p[1], exp(p[3]), exp(p[4])),
               log(1 - pi1) + l2)
    mx <- pmax(M[, 1], M[, 2])
    -sum(mx + log(exp(M[, 1] - mx) + exp(M[, 2] - mx)))
  }
  op <- stats::optim(c(-1, 2, log(0.5), log(0.8), log(0.8),
                       stats::qlogis(0.4)), nll, method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fq)) + op$value), 5e-3)
})
