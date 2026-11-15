# A weight multiplies a row's log-likelihood, as brms's weights() does.
# So in a maximum-likelihood fit with no random effect and no prior, one
# constant weight c on every row multiplies the objective by c (a prior,
# a random effect's density and REML's fixed-effect term are not
# multiplied, so there it moves the estimates): every estimate stays and
# the log-likelihood is c times the unweighted one. lme4's weights are
# precision weights instead, under which the coefficients stay and sigma
# moves; the lme4 test this item comes from (priorWeights.R) pins that
# reading, so the difference is written down here.

constant_weight_data <- function(seed = 31, n = 80L) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n))
  d$y <- 1 + 0.6 * d$x - 0.4 * d$z + stats::rnorm(n, 0, 0.7)
  d
}

test_that("a constant weight leaves a gaussian fit and multiplies logLik", {
  d <- constant_weight_data()
  d$w2 <- 2
  d$w3 <- 3
  f1 <- frm(bf(y ~ x + z), family = gaussian(), data = d)
  f2 <- frm(bf(y | weights(w2) ~ x + z), family = gaussian(), data = d)
  f3 <- frm(bf(y | weights(w3) ~ x + z, sigma ~ x), family = gaussian(),
            data = d)
  f1s <- frm(bf(y ~ x + z, sigma ~ x), family = gaussian(), data = d)
  # the estimates agree to a small fraction of their own standard errors,
  # sigma's included: the optimizer stops at a gradient tolerance, so
  # exact equality is not owed
  for (p in list(list(f1, f2), list(f1s, f3))) {
    a <- p[[1L]]
    b <- p[[2L]]
    est_a <- a$opt$par
    est_b <- b$opt$par
    se <- sqrt(diag(vcov(a, full = TRUE)))
    expect_identical(names(est_b), names(est_a))
    expect_length(se, length(est_a))
    expect_true(all(abs(est_b - est_a) <= 1e-3 * se),
                info = paste(format(abs(est_b - est_a) / se), collapse = " "))
  }
  # and the maximized log-likelihood is c times the unweighted one; the
  # mismatch of the two optima is second order in the estimates, so a
  # relative 1e-8 is loose
  l1 <- as.numeric(logLik(f1))
  expect_lt(abs(as.numeric(logLik(f2)) - 2 * l1), 1e-8 * abs(l1))
  l1s <- as.numeric(logLik(f1s))
  expect_lt(abs(as.numeric(logLik(f3)) - 3 * l1s), 1e-8 * abs(l1s))
})
