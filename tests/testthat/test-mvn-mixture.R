# mixture_mvn(): multivariate gaussian mixture components
# (mclust-style model-based clustering with linear-predictor means).

sim_mvnmix_data <- function(seed = 42, n = 400, p1 = 0.4,
                            m1 = c(0, 0), m2 = c(3, 4)) {
  set.seed(seed)
  cl <- rbinom(n, 1, p1)
  L1 <- t(chol(matrix(c(1, 0.5, 0.5, 1), 2)))
  L2 <- t(chol(matrix(c(0.5, -0.2, -0.2, 0.8), 2)))
  E <- matrix(rnorm(2 * n), 2)
  Y <- t(ifelse(matrix(cl == 1, 2, n, byrow = TRUE),
                m1 + L1 %*% E, m2 + L2 %*% E))
  dd <- data.frame(row = seq_len(n))
  dd$Y <- Y
  list(dd = dd, cl = cl)
}

test_that("mixture_mvn matches direct ML", {
  sim <- sim_mvnmix_data()
  fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = sim$dd)

  Y <- sim$dd$Y
  # hand-coded bivariate normal over Cholesky-parameterized covariances:
  # p = (m1, m2, log diag / off chol 1, same for 2, logit pi1)
  ldmvn2 <- function(Y, m, S) {
    det_ <- S[1, 1] * S[2, 2] - S[1, 2]^2
    d1 <- Y[, 1] - m[1]
    d2 <- Y[, 2] - m[2]
    q <- (S[2, 2] * d1^2 - 2 * S[1, 2] * d1 * d2 + S[1, 1] * d2^2) / det_
    -log(2 * pi) - 0.5 * log(det_) - 0.5 * q
  }
  nll <- function(p) {
    L1 <- matrix(c(exp(p[5]), p[7], 0, exp(p[6])), 2)
    L2 <- matrix(c(exp(p[8]), p[10], 0, exp(p[9])), 2)
    pi1 <- stats::plogis(p[11])
    -sum(log(pi1 * exp(ldmvn2(Y, p[1:2], L1 %*% t(L1))) +
               (1 - pi1) * exp(ldmvn2(Y, p[3:4], L2 %*% t(L2)))))
  }
  op <- stats::optim(c(0, 0, 3, 4, 0, 0, 0.5, log(0.7), log(0.8), -0.3,
                       stats::qlogis(0.4)),
                     nll, method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-5)

  # posterior class probabilities recover the simulated classes (up to
  # label swap)
  P <- mixture_probs(fit)
  expect_equal(dim(P), c(nrow(Y), 2L))
  expect_equal(unname(rowSums(P)), rep(1, nrow(Y)), tolerance = 1e-10)
  acc <- mean((P[, 1] > 0.5) == (sim$cl == 1))
  expect_gt(max(acc, 1 - acc), 0.95)

  # fitted() is the n x D mixture-mean matrix; its column means equal
  # the response column means at the gaussian-mixture ML optimum
  fv <- fitted(fit)
  expect_equal(dim(fv), dim(Y))
  expect_vector_equal(colMeans(fv), colMeans(Y), tol = 1e-4)
  expect_equal(dim(residuals(fit)), dim(Y))
})

test_that("mixture_mvn recovers the faithful clusters", {
  Y <- data.matrix(datasets::faithful)   # eruptions, waiting
  dd <- data.frame(row = seq_len(nrow(Y)))
  dd$Y <- Y
  fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)

  fx <- fixef(fit)
  m <- rbind(c(fx$mu1d1, fx$mu1d2), c(fx$mu2d1, fx$mu2d2))
  m <- m[order(m[, 1]), ]   # label-swap invariant: order by eruptions
  # reference means from the standard 2-cluster EM solution
  expect_lt(abs(m[1, 1] - 2.04), 0.3)
  expect_lt(abs(m[1, 2] - 54.5), 2)
  expect_lt(abs(m[2, 1] - 4.29), 0.3)
  expect_lt(abs(m[2, 2] - 80.0), 2)

  # class assignments agree with the eruptions < 3 threshold split
  P <- mixture_probs(fit)
  short <- Y[, 1] < 3
  agree <- mean((P[, 1] > 0.5) == short)
  expect_gt(max(agree, 1 - agree), 0.95)
})

test_that("mixture_mvn class means can depend on covariates", {
  set.seed(43)
  n <- 300
  x <- rnorm(n)
  cl <- rbinom(n, 1, 0.5)
  L <- t(chol(matrix(c(0.49, 0.15, 0.15, 0.49), 2)))
  E <- L %*% matrix(rnorm(2 * n), 2)
  Y <- cbind(ifelse(cl == 1, 0, 4) + 1.2 * x + E[1, ],
             ifelse(cl == 1, 1, 5) - 0.8 * x + E[2, ])
  dd <- data.frame(x = x)
  dd$Y <- Y

  f_x <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = dd)
  f_0 <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
  expect_gt(as.numeric(logLik(f_x)), as.numeric(logLik(f_0)) + 20)
  # both classes share the true slopes (loose recovery)
  fx <- fixef(f_x)
  for (nm in c("mu1d1", "mu2d1")) {
    expect_lt(abs(fx[[nm]][["x"]] - 1.2), 0.3)
  }
  for (nm in c("mu1d2", "mu2d2")) {
    expect_lt(abs(fx[[nm]][["x"]] + 0.8), 0.3)
  }
})

test_that("mixture_mvn dpar overrides and covariate gating work", {
  sim <- sim_mvnmix_data(seed = 44, n = 150)
  dd <- sim$dd
  dd$x <- rnorm(nrow(dd))

  # per-class per-dimension override
  f_ov <- frm(bf(Y ~ x, mu2d1 ~ 1) + mixture_mvn(K = 2, D = 2),
              data = dd)
  fx <- fixef(f_ov)
  expect_named(fx$mu1d1, c("(Intercept)", "x"))
  expect_named(fx$mu2d1, "(Intercept)")

  # gating: mixing weights take a full linear predictor
  f_gate <- frm(bf(Y ~ 1, theta1 ~ x) + mixture_mvn(K = 2, D = 2),
                data = dd)
  f_flat <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
  expect_named(fixef(f_gate)$theta1, c("(Intercept)", "x"))
  expect_gte(as.numeric(logLik(f_gate)),
             as.numeric(logLik(f_flat)) - 1e-6)
})

test_that("mixture_mvn validation and guards", {
  expect_error(mixture_mvn(2), "K >= 2")
  expect_error(mixture_mvn(1, 2), "K >= 2")
  expect_error(mixture_mvn(2, 1), "mixture\\(gaussian")

  sim <- sim_mvnmix_data(seed = 45, n = 80)
  dd <- sim$dd
  expect_error(frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 3), data = dd),
               "n x 3")
  expect_error(frm(bf(row ~ 1) + mixture_mvn(K = 2, D = 2), data = dd),
               "numeric matrix")
  # extra-parameter families refuse multivariate specs
  dd$z <- rnorm(nrow(dd))
  expect_error(
    frm(mvbf(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2),
             bf(z ~ 1) + gaussian()), data = dd),
    "multivariate"
  )
  # cens()/trunc() need a CDF, which a mixture density does not carry
  dd$cc <- 0
  expect_error(frm(bf(Y | cens(cc) ~ 1) + mixture_mvn(K = 2, D = 2),
                   data = dd),
               "CDF")

  fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
  expect_error(residuals(fit, type = "pearson"), "variance function")
  expect_error(simulate(fit), "no simulator")
})
