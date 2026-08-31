# v0.16: rr() reduced-rank covariance and mi() one-step imputation.

test_that("rr() matches glmmTMB and nests us() at full rank", {
  set.seed(61)
  n_site <- 60
  S <- 6
  L <- matrix(0, S, 2)
  L[, 1] <- c(1.2, 0.8, 0.5, -0.4, 0.9, 0.2)
  L[2:S, 2] <- c(0.7, -0.6, 0.4, 0.3, -0.5)
  f <- matrix(rnorm(n_site * 2), n_site, 2)
  u <- f %*% t(L)
  dd <- data.frame(
    y = rpois(n_site * S, exp(0.5 + as.vector(t(u)))),
    spp = factor(rep(seq_len(S), n_site)),
    site = factor(rep(seq_len(n_site), each = S))
  )
  fit <- suppressWarnings(
    frm(bf(y ~ 1 + rr(spp + 0 | site, d = 2)) + poisson(), data = dd)
  )
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    gt <- try(suppressWarnings(
      glmmTMB::glmmTMB(y ~ 1 + rr(spp + 0 | site, d = 2), data = dd,
                       family = poisson())
    ), silent = TRUE)
    if (!inherits(gt, "try-error") &&
        is.finite(as.numeric(logLik(gt)))) {
      expect_loglik_equal(fit, gt, tol = 1e-4)
    }
  }
  # rank-2 covariance really is rank 2
  V <- VarCorr(fit)[[1]]
  expect_equal(sum(eigen(V, only.values = TRUE)$values > 1e-8), 2L)

  # full-rank rr must equal us exactly
  d2 <- dd[dd$spp %in% c("1", "2", "3"), , drop = FALSE]
  d2$spp <- droplevels(d2$spp)
  fr <- suppressWarnings(
    frm(bf(y ~ 1 + rr(spp + 0 | site, d = 3)) + poisson(), data = d2)
  )
  fu <- suppressWarnings(
    frm(bf(y ~ 1 + us(spp + 0 | site)) + poisson(), data = d2)
  )
  expect_loglik_equal(fr, fu, tol = 1e-4)

  # methods run through the coefficient-space expansion
  expect_length(fitted(fit), nrow(dd))
  expect_equal(unname(predict(fit, newdata = dd[1:6, ],
                              type = "response")),
               unname(fitted(fit)[1:6]), tolerance = 1e-10)
  r <- ranef(fit)
  expect_equal(dim(r[[1]]), c(60L, 6L))
  expect_equal(nrow(as.data.frame(VarCorr(fit))), 6L + 15L)
  s <- simulate(fit, nsim = 2, re.form = NA)
  expect_equal(nrow(s), nrow(dd))
  # se.fit works since v0.17 (loadings Jacobian)
  ps <- predict(fit, se.fit = TRUE)
  expect_true(all(is.finite(ps$se.fit)))
  expect_error(frm(bf(y ~ 1 + rr(spp + 0 | site, d = 9)) + poisson(),
                   data = dd), "must not exceed")
})

test_that("mi() matches the closed-form marginal likelihood", {
  set.seed(71)
  n <- 300
  z <- rnorm(n)
  x <- rnorm(n, 0.5 + 0.8 * z, 0.7)
  y <- rnorm(n, 1 + 0.6 * x + 0.3 * z, 0.9)
  x_mis <- x
  mis <- sort(sample(n, 90))
  x_mis[mis] <- NA
  dd <- data.frame(y = y, x = x_mis, z = z)

  fit <- frm(bf(y ~ mi(x) + z) + gaussian() +
               bf(x | mi() ~ z) + gaussian(), data = dd)

  # linear-gaussian case: missing x integrates out in closed form, and
  # the Laplace approximation is exact
  obs <- !is.na(x_mis)
  nll <- function(p) {
    sy <- exp(p[6]); sx <- exp(p[7])
    mux <- p[4] + p[5] * z
    -(sum(dnorm(y[obs], p[1] + p[2] * x_mis[obs] + p[3] * z[obs], sy,
                log = TRUE)) +
        sum(dnorm(x_mis[obs], mux[obs], sx, log = TRUE)) +
        sum(dnorm(y[!obs], p[1] + p[2] * mux[!obs] + p[3] * z[!obs],
                  sqrt(sy^2 + p[2]^2 * sx^2), log = TRUE)))
  }
  op <- stats::optim(c(1, 0.6, 0.3, 0.5, 0.8, log(0.9), log(0.7)), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)
  expect_equal(unname(fixef(fit)$y_mu[c("(Intercept)", "mix", "z")]),
               op$par[1:3], tolerance = 1e-4)
  expect_equal(unname(fixef(fit)$x_mu), op$par[4:5], tolerance = 1e-4)

  # imputations track the truth (x | y, z posterior modes)
  expect_gt(cor(fit$estimates$miss, x[mis]), 0.5)

  # imputation uncertainty flows into the coefficient SEs: the mi(x)
  # SE exceeds the complete-data SE
  s_mi <- summary(fit)$coefficients$y_mu["mix", "Std. Error"]
  f_full <- frm(bf(y ~ x + z) + gaussian(),
                data = data.frame(y = y, x = x, z = z))
  s_full <- summary(f_full)$coefficients$mu["x", "Std. Error"]
  expect_gt(s_mi, s_full)

  # predict: in-sample and newdata (which must be complete)
  expect_length(predict(fit, resp = "y"), n)
  expect_length(predict(fit, newdata = data.frame(x = 0.5, z = 0),
                        resp = "y"), 1L)
  expect_error(predict(fit, newdata = data.frame(x = NA_real_, z = 0),
                       resp = "y"), "complete")
})

test_that("mi() degenerates and guards correctly", {
  set.seed(73)
  n <- 150
  z <- rnorm(n)
  x <- rnorm(n, 0.5 + 0.8 * z, 0.7)
  y <- rnorm(n, 1 + 0.6 * x + 0.3 * z, 0.9)
  d2 <- data.frame(y = y, x = x, z = z)

  # no missing values: the joint fit equals the two separate fits
  f2 <- frm(bf(y ~ mi(x) + z) + gaussian() +
              bf(x | mi() ~ z) + gaussian(), data = d2)
  fy <- frm(bf(y ~ x + z) + gaussian(), data = d2)
  fx <- frm(bf(x ~ z) + gaussian(), data = d2)
  expect_lt(abs(as.numeric(logLik(f2)) -
                  as.numeric(logLik(fy)) - as.numeric(logLik(fx))),
            1e-6)

  # rows with NA in a non-mi variable still drop
  d3 <- d2
  d3$x[1:10] <- NA
  d3$z[1:5] <- NA
  f3 <- frm(bf(y ~ mi(x) + z) + gaussian() +
              bf(x | mi() ~ z) + gaussian(), data = d3)
  expect_equal(nobs(f3), n - 5L)
  expect_length(f3$estimates$miss, 5L)   # rows 6:10 latent

  # guards
  expect_error(frm(bf(y ~ mi(x) + z) + gaussian(), data = d2),
               "matching imputation model")
  expect_error(frm(bf(y ~ mi(x) + z) + gaussian() +
                     bf(x | mi() ~ z) + poisson(), data = d3),
               "gaussian or student")
  expect_error(frm(bf(y ~ mi(x) * z) + gaussian() +
                     bf(x | mi() ~ z) + gaussian(), data = d2),
               "standalone")
})
