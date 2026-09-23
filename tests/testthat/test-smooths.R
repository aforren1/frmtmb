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
  expect_lt(max(abs(fitted(fit)[, "Estimate"] - fitted(ref))), 1e-2)
})

test_that("s(x) poisson matches mgcv::gam ML", {
  set.seed(42)
  n <- 400
  dd <- data.frame(x = runif(n))
  dd$y <- rpois(n, exp(1 + sin(2 * dd$x)))
  fit <- frm(bf(y ~ s(x)) + poisson(), data = dd)
  ref <- mgcv::gam(y ~ s(x), data = dd, family = poisson, method = "ML")
  expect_lt(abs(as.numeric(logLik(fit)) - gam_ml_loglik(ref)), 1e-2)
  expect_lt(max(abs(fitted(fit)[, "Estimate"] - fitted(ref))), 0.05)
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
  vc <- varcorr_matrices(fit)
  sd_frm <- sqrt(vc[["1 | g"]][1, 1])
  expect_gt(sd_frm, 0.2)
  expect_lt(max(abs(fitted(fit)[, "Estimate"] - fitted(ref))), 1e-2)
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
  mu_frm <- frm_linpred(fit, type = "response")
  sig_frm <- frm_linpred(fit, dpar = "sigma", type = "response")
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
  expect_lt(max(abs(fitted(fit)[, "Estimate"] - fitted(ref))), 0.1)
})

test_that("t2 smooths predict on newdata", {
  set.seed(45)
  n <- 400
  dd <- data.frame(x = runif(n), w = runif(n))
  dd$y <- sin(2 * dd$x) * cos(2 * dd$w) + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ t2(x, w)) + gaussian(), data = dd)

  # a t2 basis carries a separate prediction constraint (mgcv `Cp`) that
  # PredictMat honors; the frame drops it (smoothCon modCon = 3) so the
  # newdata basis IS the fit basis and the round-trip is exact
  expect_equal(frm_linpred(fit, newdata = dd), frm_linpred(fit),
               tolerance = 1e-10)
  expect_equal(frm_linpred(fit, newdata = dd, re_formula = NA),
               frm_linpred(fit),
               tolerance = 1e-10)
  s_in <- frm_linpred(fit, se.fit = TRUE)
  s_nd <- frm_linpred(fit, newdata = dd, se.fit = TRUE)
  expect_equal(s_nd$se.fit, s_in$se.fit, tolerance = 1e-10)

  # off the training grid, against the same model fitted by mgcv
  ref <- mgcv::gam(y ~ t2(x, w), data = dd, method = "ML")
  nd <- expand.grid(x = seq(0.05, 0.95, length.out = 11),
                    w = seq(0.05, 0.95, length.out = 11))
  p <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
  expect_lt(max(abs(p$fit - as.numeric(predict(ref, newdata = nd)))), 1e-3)
  expect_true(all(is.finite(p$se.fit)) && all(p$se.fit > 0))
  # se blows up under extrapolation
  expect_gt(frm_linpred(fit, newdata = data.frame(x = 1.6, w = 1.6),
                    se.fit = TRUE)$se.fit,
            max(p$se.fit))
})

test_that("t2 newdata prediction survives by=, dpars and random effects", {
  set.seed(7)
  n <- 500
  dd <- data.frame(x = runif(n), z = runif(n),
                   f = factor(sample(c("a", "b"), n, TRUE)),
                   g = factor(sample(letters[1:8], n, TRUE)),
                   u = rnorm(n))
  dd$y <- sin(3 * dd$x) * dd$z + 0.5 * dd$u +
    as.numeric(dd$f) + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ u + t2(x, z, by = f) + (1 | g)) + gaussian(), data = dd)
  expect_equal(frm_linpred(fit, newdata = dd), frm_linpred(fit),
               tolerance = 1e-10)
  nd <- dd[1:20, ]
  nd$x <- nd$x * 0.5
  expect_true(all(is.finite(frm_linpred(fit, newdata = nd,
                                        se.fit = TRUE)$se.fit)))
  # a two-variable conditional effect goes through the same newdata path
  ce <- conditional_effects(fit, effects = "x:z")
  expect_true(all(is.finite(ce[[1]]$estimate__)))

  set.seed(8)
  d3 <- data.frame(x = runif(400), z = runif(400))
  d3$y <- rnorm(400, sin(3 * d3$x), exp(-1 + 0.5 * d3$z))
  f3 <- frm(bf(y ~ t2(x, z), sigma ~ s(z)) + gaussian(), data = d3)
  expect_equal(frm_linpred(f3, newdata = d3), frm_linpred(f3),
               tolerance = 1e-10)
  expect_equal(frm_linpred(f3, dpar = "sigma", newdata = d3),
               frm_linpred(f3, dpar = "sigma"), tolerance = 1e-10)
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
  expect_equal(frm_linpred(fit, newdata = dd), frm_linpred(fit),
               tolerance = 1e-8)
  # smooth curve survives population-level predictions
  expect_equal(frm_linpred(fit, newdata = dd, re_formula = NA),
               frm_linpred(fit),
               tolerance = 1e-8)
  # interpolation with standard errors
  nd <- data.frame(x = seq(0.1, 0.9, length.out = 11))
  p <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
  expect_true(all(is.finite(p$fit)) && all(p$se.fit > 0))
  # se blows up under extrapolation
  expect_gt(frm_linpred(fit, newdata = data.frame(x = 1.3),
                    se.fit = TRUE)$se.fit,
            max(p$se.fit))
})

test_that("REML on a location-scale smooth differs from mgcv's by design", {
  # A MEASUREMENT, not an agreement. frmtmb's REML = TRUE integrates the
  # mu coefficients and keeps sigma's outer (the double-GLM REML), while
  # mgcv's REML for gaulss() integrates the coefficients of both
  # predictors, so the two criteria put the sigma smooth in different
  # places. ?frm says so under REML. Seeds 44 to 48 of
  # dev/correct-reml-mgcv.R (dev/correct-log/reml-mgcv-lane.txt), as
  # log(frmtmb lambda) - log(mgcv sp), lambda = 1 / sd^2:
  #   mu smooth:    ML 0.129 to 0.262, REML -0.011 to -0.005
  #   sigma smooth: ML 0.138 to 0.165, REML 0.142 to 0.170
  # REML closes the mu smooth's gap and leaves the sigma smooth's where
  # ML left it, which is the difference in what each integrates.
  set.seed(44)
  n <- 500
  dd <- data.frame(x = runif(n), z = runif(n))
  dd$y <- rnorm(n, sin(3 * dd$x), exp(0.8 * cos(5 * dd$z) - 0.5))
  gap <- function(method) {
    fit <- frm(bf(y ~ s(x), sigma ~ s(z)) + gaussian(), data = dd,
               REML = identical(method, "REML"))
    ref <- mgcv::gam(list(y ~ s(x), ~ s(z)), data = dd,
                     family = mgcv::gaulss(b = 0), method = method)
    lam <- vapply(fit$frame$re_blocks, function(b) {
      exp(-2 * fit$estimates$theta[b$theta_idx])
    }, 1)
    dp <- vapply(fit$frame$re_blocks, `[[`, "", "dpar")
    mu_g <- fitted(ref)[, 1]
    sg_g <- 1 / fitted(ref)[, 2]
    lam <- unname(lam)
    c(sp_mu = abs(log(lam[dp == "mu"]) - log(unname(ref$sp[1]))),
      sp_sigma = abs(log(lam[dp == "sigma"]) - log(unname(ref$sp[2]))),
      curve_mu = max(abs(as.numeric(frm_linpred(fit, type = "response")) -
                           mu_g)) / diff(range(mu_g)),
      curve_sigma = max(abs(as.numeric(frm_linpred(
        fit, dpar = "sigma", type = "response")) - sg_g)) /
        diff(range(sg_g)))
  }
  ml <- gap("ML")
  reml <- gap("REML")
  # the sigma smooth's gap under REML is the ML one: measured ratio
  # 1.026 to 1.037 over the five seeds
  expect_gt(reml[["sp_sigma"]] / ml[["sp_sigma"]], 0.8)
  expect_lt(reml[["sp_sigma"]] / ml[["sp_sigma"]], 1.25)
  # and it is the larger one by far: measured 13.2 to 28.6 times the mu
  # smooth's REML gap
  expect_gt(reml[["sp_sigma"]] / reml[["sp_mu"]], 10)
  # the mu smooth is where the criteria agree: REML closes most of its
  # ML gap (measured 22 to 28 times smaller, and the curve gap 6.7 to
  # 21 times smaller) ...
  expect_lt(reml[["sp_mu"]] / ml[["sp_mu"]], 0.2)
  expect_lt(reml[["curve_mu"]] / ml[["curve_mu"]], 0.2)
  # ... and the sigma curve under REML is further from mgcv's than the mu
  # curve, measured 23.5 to 65.5 times further relative to each curve's
  # range
  expect_gt(reml[["curve_sigma"]] / reml[["curve_mu"]], 10)
})
