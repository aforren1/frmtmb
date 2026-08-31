# mgcv::gam(method = "ML") minimizes the negative log marginal likelihood
# (stored in $gcv.ubre), the same criterion frmtmb maximizes with the
# smoothing variance as a variance component; for gaussian models the
# Laplace approximation is exact, so agreement is tight.
gam_ml_loglik <- function(g) -as.numeric(g$gcv.ubre)

test_that("s(x) gaussian matches mgcv::gam ML", {
  set.seed(41)
  n <- 300
  dd <- data.frame(x = runif(n))
  dd$y <- sin(3 * dd$x) + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ s(x)) + gaussian(), data = dd)
  ref <- mgcv::gam(y ~ s(x), data = dd, method = "ML")
  expect_lt(abs(as.numeric(logLik(fit)) - gam_ml_loglik(ref)), 1e-4)
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 1e-2)
})

test_that("s(x) poisson matches mgcv::gam ML", {
  set.seed(42)
  n <- 400
  dd <- data.frame(x = runif(n))
  dd$y <- rpois(n, exp(1 + sin(2 * dd$x)))
  fit <- frm(bf(y ~ s(x)) + poisson(), data = dd)
  ref <- mgcv::gam(y ~ s(x), data = dd, family = poisson, method = "ML")
  expect_lt(abs(as.numeric(logLik(fit)) - gam_ml_loglik(ref)), 1e-2)
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 0.05)
})

test_that("s(x) + (1|g) matches gam with a re smooth", {
  set.seed(43)
  n <- 400
  dd <- data.frame(x = runif(n), g = factor(rep(1:20, 20)))
  dd$y <- sin(3 * dd$x) + rnorm(20, 0, 0.5)[dd$g] + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ s(x) + (1 | g)) + gaussian(), data = dd)
  ref <- mgcv::gam(y ~ s(x) + s(g, bs = "re"), data = dd, method = "ML")
  expect_lt(abs(as.numeric(logLik(fit)) - gam_ml_loglik(ref)), 1e-4)
  # group effect present and fitted values agree with gam
  vc <- VarCorr(fit)
  sd_frm <- sqrt(vc[["1 | g"]][1, 1])
  expect_gt(sd_frm, 0.2)
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 1e-2)
})

test_that("smooths work in dpar formulas (sigma ~ s(z))", {
  set.seed(44)
  n <- 500
  dd <- data.frame(x = runif(n), z = runif(n))
  dd$y <- rnorm(n, sin(3 * dd$x), exp(0.3 * cos(3 * dd$z) - 0.5))
  fit <- frm(bf(y ~ s(x), sigma ~ s(z)) + gaussian(), data = dd)
  # mgcv's gaulss fits the same model space (its second predictor is
  # log(1/sigma)); its ML criterion drops constants so compare the fitted
  # mean and sigma curves, not the criterion values
  ref <- mgcv::gam(list(y ~ s(x), ~ s(z)), data = dd,
                   family = mgcv::gaulss(b = 0), method = "ML")
  mu_frm <- predict(fit, type = "response")
  sig_frm <- predict(fit, dpar = "sigma", type = "response")
  sig_gam <- 1 / fitted(ref)[, 2]
  expect_lt(max(abs(mu_frm - fitted(ref)[, 1])), 0.05)
  expect_lt(max(abs(sig_frm - sig_gam)), 0.05)
})

test_that("t2 tensor smooths fit and match gam ML", {
  set.seed(45)
  n <- 400
  dd <- data.frame(x = runif(n), w = runif(n))
  dd$y <- sin(2 * dd$x) * cos(2 * dd$w) + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ t2(x, w)) + gaussian(), data = dd)
  ref <- mgcv::gam(y ~ t2(x, w), data = dd, method = "ML")
  expect_lt(abs(as.numeric(logLik(fit)) - gam_ml_loglik(ref)), 0.05)
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 0.1)
})

test_that("te() errors with guidance; smooth predictions round-trip", {
  expect_error(frm(bf(y ~ te(x, w)) + gaussian(),
                   data = NULL, dry_run = "spec"),
               "t2")

  set.seed(46)
  n <- 300
  dd <- data.frame(x = runif(n))
  dd$y <- sin(3 * dd$x) + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ s(x)) + gaussian(), data = dd)
  # newdata equal to training reproduces in-sample predictions
  expect_equal(predict(fit, newdata = dd), predict(fit), tolerance = 1e-8)
  # smooth curve survives population-level predictions
  expect_equal(predict(fit, newdata = dd, re.form = NA), predict(fit),
               tolerance = 1e-8)
  # interpolation with standard errors
  nd <- data.frame(x = seq(0.1, 0.9, length.out = 11))
  p <- predict(fit, newdata = nd, se.fit = TRUE)
  expect_true(all(is.finite(p$fit)) && all(p$se.fit > 0))
  # se blows up under extrapolation
  expect_gt(predict(fit, newdata = data.frame(x = 1.3),
                    se.fit = TRUE)$se.fit,
            max(p$se.fit))
})
