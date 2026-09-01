# v0.17: mixture families, me() measurement error, equalto(), cs()
# category-specific effects, rr se.fit, insight methods.

test_that("gaussian and poisson mixtures match direct ML", {
  set.seed(81)
  n <- 400
  y <- ifelse(rbinom(n, 1, 0.35) == 1, rnorm(n, 3, 0.6), rnorm(n, 0, 1))
  dd <- data.frame(y = y)
  fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
  nll <- function(p) {
    pi1 <- stats::plogis(p[5])
    -sum(log(pi1 * stats::dnorm(y, p[1], exp(p[3])) +
               (1 - pi1) * stats::dnorm(y, p[2], exp(p[4]))))
  }
  op <- stats::optim(c(0, 3, 0, log(0.6), 0.6), nll, method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)
  expect_length(fitted(fit), n)
  expect_lt(abs(mean(fitted(fit)) - mean(y)), 1e-6)
  s <- simulate(fit, nsim = 1, seed = 1)
  expect_lt(abs(stats::sd(s[[1]]) - stats::sd(y)), 0.3)

  y2 <- ifelse(rbinom(n, 1, 0.4) == 1, rpois(n, 2), rpois(n, 12))
  f2 <- frm(bf(y ~ 1) + mixture(poisson(), poisson()),
            data = data.frame(y = y2))
  nll2 <- function(p) {
    pi1 <- stats::plogis(p[3])
    -sum(log(pi1 * stats::dpois(y2, exp(p[1])) +
               (1 - pi1) * stats::dpois(y2, exp(p[2]))))
  }
  op2 <- stats::optim(c(log(2), log(12), 0), nll2, method = "BFGS",
                      control = list(reltol = 1e-13))
  expect_lt(abs(as.numeric(logLik(f2)) + op2$value), 1e-6)
})

test_that("mixing weights can depend on covariates (experts)", {
  set.seed(82)
  n <- 500
  x <- rnorm(n)
  cl <- rbinom(n, 1, stats::plogis(-0.5 + 1.5 * x))
  y <- ifelse(cl == 1, rnorm(n, 2, 0.5), rnorm(n, -1, 0.8))
  dd <- data.frame(y = y, x = x)
  f_moe <- frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()),
               data = dd)
  f_flat <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
  expect_gt(as.numeric(logLik(f_moe)), as.numeric(logLik(f_flat)) + 10)
  expect_gt(abs(fixef(f_moe)$theta1[["x"]]), 0.5)
})

test_that("mi(sd) measurement error matches the closed form", {
  set.seed(91)
  n <- 250
  z <- rnorm(n)
  t_true <- rnorm(n, 0.5 + 0.8 * z, 0.7)
  su <- 0.4
  x_obs <- t_true + rnorm(n, 0, su)
  y <- rnorm(n, 1 + 0.6 * t_true + 0.3 * z, 0.9)
  dd <- data.frame(y = y, x = x_obs, z = z, su = su)

  fit <- frm(bf(y ~ mi(x) + z) + gaussian() +
               bf(x | mi(su) ~ z) + gaussian(), data = dd)

  # joint (y, x_obs) | z is bivariate normal
  nll <- function(p) {
    sy <- exp(p[6]); sx <- exp(p[7])
    mux <- p[4] + p[5] * z
    muy <- p[1] + p[2] * mux + p[3] * z
    Sy <- p[2]^2 * sx^2 + sy^2
    Sx <- sx^2 + su^2
    Sxy <- p[2] * sx^2
    det_ <- Sy * Sx - Sxy^2
    dy <- y - muy
    dx <- x_obs - mux
    q <- (Sx * dy^2 - 2 * Sxy * dy * dx + Sy * dx^2) / det_
    -sum(-log(2 * pi) - 0.5 * log(det_) - 0.5 * q)
  }
  op <- stats::optim(c(1, 0.6, 0.3, 0.5, 0.8, log(0.9), log(0.7)), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-5)
  expect_equal(unname(fixef(fit)$y_mu[c("(Intercept)", "mix", "z")]),
               op$par[1:3], tolerance = 1e-3)
  # attenuation corrected: naive regression on x_obs shrinks b1
  naive <- stats::coef(stats::lm(y ~ x_obs + z))[["x_obs"]]
  expect_gt(fixef(fit)$y_mu[["mix"]], naive + 0.02)
  # latent truths tracked
  expect_gt(cor(fit$estimates$miss, t_true), 0.8)
})

test_that("equalto() fixes the covariance exactly", {
  set.seed(93)
  ng <- 30
  V <- matrix(c(1.2, 0.5, 0.5, 0.8), 2, 2)
  b <- t(chol(V)) %*% matrix(rnorm(2 * ng), 2)
  f <- factor(rep(c("a", "b"), ng))
  g <- factor(rep(seq_len(ng), each = 2))
  y <- 1 + as.vector(b) + rnorm(2 * ng, 0, 0.5)
  dd <- data.frame(y = y, f = f, g = g)

  fit <- frm(bf(y ~ 1 + equalto(f + 0 | g, V)) + gaussian(), data = dd)
  expect_length(fit$estimates$theta, 0L)
  expect_equal(unname(VarCorr(fit)[[1]]), unname(V))

  # direct marginal-gaussian reference over (mu, log sigma)
  Z <- stats::model.matrix(~ 0 + g:f)
  ord <- order(rep(seq_len(ng), each = 2), rep(1:2, ng))
  Sig_b <- kronecker(diag(ng), V)
  nll <- function(p) {
    S <- Z[, ord] %*% Sig_b %*% t(Z[, ord]) + exp(2 * p[2]) * diag(2 * ng)
    -stats::dnorm(0)   # placeholder replaced below
  }
  # simpler: per-group 2x2 marginal (groups independent)
  nll <- function(p) {
    S <- V + exp(2 * p[2]) * diag(2)
    Si <- solve(S)
    ld <- determinant(S)$modulus
    r <- matrix(y - p[1], 2)
    -sum(-log(2 * pi) - 0.5 * ld - 0.5 * colSums(r * (Si %*% r)))
  }
  op <- stats::optim(c(1, log(0.5)), nll, method = "BFGS",
                     control = list(reltol = 1e-13))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)
})

test_that("cs() category-specific effects match direct ML", {
  set.seed(95)
  n <- 500
  x <- rnorm(n)
  tau_t <- c(-1, 0.3, 1.5)
  g_t <- c(0.8, 0.2, -0.4)   # category-specific slopes
  P <- local({
    Fm <- stats::plogis(outer(-0.3 * x, tau_t, `+`) - outer(x, g_t))
    Pm <- matrix(0, n, 4)
    surv <- rep(1, n)
    for (k in 1:3) {
      Pm[, k] <- Fm[, k] * surv
      surv <- surv * (1 - Fm[, k])
    }
    Pm[, 4] <- surv
    Pm
  })
  y <- apply(P, 1, function(p) sample.int(4, 1, prob = p))
  dd <- data.frame(y = y, x = x)

  fit <- frm(bf(y ~ x + cs(x)) + sratio(), data = dd)
  nll <- function(p) {
    Fm <- stats::plogis(outer(-p[1] * x, p[2:4], `+`) -
                          outer(x, p[5:7]))
    Pm <- matrix(0, n, 4)
    surv <- rep(1, n)
    for (k in 1:3) {
      Pm[, k] <- Fm[, k] * surv
      surv <- surv * (1 - Fm[, k])
    }
    Pm[, 4] <- surv
    -sum(log(Pm[cbind(seq_len(n), y)]))
  }
  op <- stats::optim(c(0.3, tau_t, g_t), nll, method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)
  # a plain sratio fit is strictly worse (cs is real here)
  f0 <- frm(bf(y ~ x) + sratio(), data = dd)
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(f0)) + 5)
  # cs works for cratio and acat, refuses cumulative
  expect_no_error(frm(bf(y ~ x + cs(x)) + cratio(), data = dd))
  expect_no_error(frm(bf(y ~ x + cs(x)) + acat(), data = dd))
  expect_error(frm(bf(y ~ x + cs(x)) + cumulative(), data = dd),
               "sratio, cratio, or acat")
})

test_that("cs() variables reach the model frame on their own", {
  # cs(z) with z absent from every other term: the combined model frame
  # has to collect cs() variables too, or assembly fails on lookup
  set.seed(96)
  dd <- data.frame(y = sample(1:4, 300, TRUE), x = rnorm(300),
                   z = rnorm(300))
  fr <- frm(bf(y ~ x + cs(z)) + sratio(), data = dd, dry_run = "frame")
  cs <- fr$linpreds[["y.mu"]]$cs
  expect_length(cs, 1)
  expect_vector_equal(cs[[1]]$vals, dd$z, tol = 1e-12)
})

test_that("rr() se.fit runs and is parameterization-invariant", {
  set.seed(97)
  n_site <- 50
  L <- matrix(c(1, 0.6, 0.3, 0, 0.8, -0.5), 3, 2)
  f <- matrix(rnorm(n_site * 2), n_site, 2)
  dd <- data.frame(
    y = rpois(n_site * 3, exp(0.4 + as.vector(t(f %*% t(L))))),
    spp = factor(rep(1:3, n_site)),
    site = factor(rep(seq_len(n_site), each = 3))
  )
  # full-rank rr and us are the same model in different
  # parameterizations: var(eta-hat) must agree
  fr <- suppressWarnings(
    frm(bf(y ~ 1 + rr(spp + 0 | site, d = 3)) + poisson(), data = dd)
  )
  fu <- suppressWarnings(
    frm(bf(y ~ 1 + us(spp + 0 | site)) + poisson(), data = dd)
  )
  pr <- predict(fr, se.fit = TRUE)
  pu <- predict(fu, se.fit = TRUE)
  expect_vector_equal(pr$fit, pu$fit, tol = 1e-3)
  expect_vector_equal(pr$se.fit, pu$se.fit, tol = 0.02)
  # newdata path too
  nr <- predict(fr, newdata = dd[1:3, ], se.fit = TRUE)
  expect_vector_equal(nr$se.fit, pr$se.fit[1:3], tol = 1e-6)
})

test_that("insight sees the mixed-model structure", {
  skip_if_not_installed("insight")
  set.seed(99)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.7)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  expect_true(insight::is_mixed_model(fit))
  expect_equal(insight::find_random(fit)$random, "g")
  fl <- insight::find_formula(fit)
  expect_equal(deparse1(fl$conditional), "y ~ x")
  gp <- insight::get_parameters(fit)
  expect_true(all(c("Parameter", "Estimate") %in% names(gp)))
  expect_equal(nrow(gp), 3L)
  expect_equal(insight::find_statistic(fit), "z-statistic")
  expect_equal(insight::link_inverse(fit)(0.5), 0.5)
})
