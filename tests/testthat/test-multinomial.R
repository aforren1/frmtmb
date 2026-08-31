sim_multinom_data <- function(seed = 81, n = 300, size = 10) {
  set.seed(seed)
  x <- rnorm(n)
  e2 <- 0.5 + 0.8 * x
  e3 <- -0.3 + 0.4 * x
  P <- cbind(1, exp(e2), exp(e3))
  P <- P / rowSums(P)
  Y <- t(vapply(seq_len(n),
                function(i) stats::rmultinom(1, size, P[i, ])[, 1],
                numeric(3)))
  dd <- data.frame(x = x)
  dd$Y <- Y
  dd
}

test_that("multinomial matrix response matches a hand-rolled reference", {
  dd <- sim_multinom_data()
  fit <- frm(bf(Y ~ x) + multinomial(K = 3), data = dd)

  Y <- dd$Y; x <- dd$x
  nll_ref <- function(p) {
    e2 <- p$b2[1] + p$b2[2] * x
    e3 <- p$b3[1] + p$b3[2] * x
    denom <- 1 + exp(e2) + exp(e3)
    ll <- -rowSums(Y) * log(denom) + Y[, 2] * e2 + Y[, 3] * e3
    -sum(ll) - sum(lgamma(rowSums(Y) + 1) - rowSums(lgamma(Y + 1)))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b2 = c(0, 0), b3 = c(0, 0)),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  est <- opt$par
  expect_vector_equal(fixef(fit)$mu2, est[1:2], tol = 1e-4)
  expect_vector_equal(fixef(fit)$mu3, est[3:4], tol = 1e-4)
})

test_that("multinomial matches nnet::multinom", {
  skip_if_not_installed("nnet")
  dd <- sim_multinom_data(seed = 82)
  fit <- frm(bf(Y ~ x) + multinomial(K = 3), data = dd)
  ref <- nnet::multinom(Y ~ x, data = dd, trace = FALSE,
                        reltol = 1e-14, maxit = 500)
  # multinom's loglik omits the multinomial coefficient
  const <- sum(lgamma(rowSums(dd$Y) + 1)) - sum(lgamma(dd$Y + 1))
  expect_lt(abs((as.numeric(logLik(fit)) - const) -
                  as.numeric(logLik(ref))), 1e-3)
  cf <- coef(ref)
  expect_vector_equal(fixef(fit)$mu2, cf[1, ], tol = 1e-2)
  expect_vector_equal(fixef(fit)$mu3, cf[2, ], tol = 1e-2)
})

test_that("per-category dpar formulas can be overridden", {
  dd <- sim_multinom_data(seed = 83)
  dd$z <- rnorm(nrow(dd))
  fit <- frm(bf(Y ~ x, mu3 ~ z) + multinomial(K = 3), data = dd)
  expect_named(fixef(fit)$mu2, c("(Intercept)", "x"))
  expect_named(fixef(fit)$mu3, c("(Intercept)", "z"))
})

test_that("multinomial validation", {
  dd <- sim_multinom_data()
  expect_error(frm(bf(Y ~ x) + multinomial(K = 4), data = dd),
               "n x 4")
  expect_error(multinomial(), "number of categories")
  expect_error(fitted(frm(bf(Y ~ x) + multinomial(K = 3), data = dd)),
               "not defined")
})
