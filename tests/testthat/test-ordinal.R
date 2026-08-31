sim_ord_data <- function(seed = 91, n = 600, K = 4) {
  set.seed(seed)
  x <- rnorm(n)
  eta <- 1.2 * x
  tau <- c(-1, 0.3, 1.4)
  u <- runif(n)
  p <- plogis(outer(rep(1, n), tau) - eta)
  y <- rowSums(u > cbind(p, 1)) + 1L
  data.frame(y = y, x = x)
}

# raw thresholds -> ordered thresholds
tau_of <- function(fit) {
  raw <- fit$estimates$tau_raw
  cumsum(c(raw[1], exp(raw[-1])))
}

test_that("cumulative logit matches MASS::polr", {
  skip_if_not_installed("MASS")
  dd <- sim_ord_data()
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  ref <- MASS::polr(factor(y) ~ x, data = dd, method = "logistic",
                    Hess = TRUE)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-5)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-3)
  expect_vector_equal(tau_of(fit), ref$zeta, tol = 1e-3)
})

test_that("cumulative probit matches MASS::polr", {
  skip_if_not_installed("MASS")
  dd <- sim_ord_data(seed = 92)
  fit <- frm(bf(y ~ x) + cumulative(link = "probit"), data = dd)
  ref <- MASS::polr(factor(y) ~ x, data = dd, method = "probit")
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-4)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-3)
})

test_that("ordered factor responses and threshold-only models work", {
  dd <- sim_ord_data(seed = 93, n = 300)
  dd$yf <- factor(dd$y, ordered = TRUE)
  f1 <- frm(bf(yf ~ x) + cumulative(), data = dd)
  f2 <- frm(bf(y ~ x) + cumulative(), data = dd)
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f2))), 1e-8)

  f0 <- frm(bf(y ~ 1) + cumulative(), data = dd)
  # threshold-only: thresholds are the sample cumulative logits
  p <- cumsum(prop.table(table(dd$y)))[-4]
  expect_vector_equal(tau_of(f0), qlogis(p), tol = 1e-4)
})

test_that("ordinal random-intercept model matches a hand-rolled reference", {
  set.seed(94)
  n <- 800
  g <- factor(rep(1:40, 20))
  x <- rnorm(n)
  eta <- 0.8 * x + rnorm(40, 0, 0.9)[g]
  tau <- c(-0.8, 0.8)
  u <- runif(n)
  p <- plogis(outer(rep(1, n), tau) - eta)
  dd <- data.frame(y = rowSums(u > cbind(p, 1)) + 1L, x = x, g = g)
  fit <- frm(bf(y ~ x + (1 | g)) + cumulative(), data = dd)

  yv <- dd$y; xv <- dd$x; gi <- as.integer(dd$g)
  nll_ref <- function(pp) {
    "[<-" <- RTMB::ADoverload("[<-")
    "c" <- RTMB::ADoverload("c")
    nll <- -sum(RTMB::dnorm(pp$u, 0, exp(pp$lsd), log = TRUE))
    eta <- pp$beta * xv + pp$u[gi]
    t2 <- c(pp$t1, pp$t1 + exp(pp$ld))
    Fv <- function(z) 1 / (1 + exp(-z))
    i1 <- as.numeric(yv == 1); i3 <- as.numeric(yv == 3)
    up <- Fv(t2[pmin(yv, 2)] - eta) * (1 - i3) + i3
    lo <- Fv(t2[pmax(yv - 1, 1)] - eta) * (1 - i1)
    nll - sum(log(up - lo))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(beta = 0, t1 = -0.5, ld = 0, lsd = 0,
                              u = numeric(40)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-5)
})
