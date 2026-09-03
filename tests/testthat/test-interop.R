test_that("emmeans marginal means match glmmTMB", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("glmmTMB")
  set.seed(131)
  n <- 300
  dd <- data.frame(
    f = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    x = rnorm(n),
    g = factor(rep(1:15, length.out = n))
  )
  dd$y <- rnorm(n, 1 + c(a = 0, b = 1, c = -1)[dd$f] + 0.5 * dd$x +
                  rnorm(15, 0, 0.7)[dd$g], 1)

  fit <- frm(bf(y ~ f + x + (1 | g)) + gaussian(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ f + x + (1 | g), data = dd, REML = FALSE)

  em_f <- as.data.frame(emmeans::emmeans(fit, "f"))
  em_g <- as.data.frame(emmeans::emmeans(ref, "f"))
  expect_vector_equal(em_f$emmean, em_g$emmean, tol = 1e-4)
  expect_vector_equal(em_f$SE, em_g$SE, tol = 1e-4)

  ct_f <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit, "f"),
                                          "pairwise"))
  expect_length(ct_f$estimate, 3)
})

test_that("as_tmbstan hands the objective to NUTS", {
  skip_if_not_installed("tmbstan")
  set.seed(132)
  dd <- data.frame(x = rnorm(80))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  sf <- suppressWarnings(as_tmbstan(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 1))
  expect_s4_class(sf, "stanfit")
  skip_if_not_installed("rstan")
  dr <- rstan::extract(sf, "beta")$beta
  # judged against the chain's own spread: a seeded chain is not
  # platform-deterministic, and this asserts wiring, not mixing
  if (sampler_gates_on()) {
    expect_lt(abs(mean(dr[, 1]) - fixef(fit)$mu[[1]]),
              5 * stats::sd(dr[, 1]) + 1e-8)
  }
})

# --- lme4::getME ------------------------------------------------------

getme_fit <- local({
  data(sleepstudy, package = "lme4")
  frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
      data = sleepstudy)
})

test_that("getME extracts the designs and agrees with the frame", {
  fit <- getme_fit
  lp <- fit$frame$linpreds[["Reaction.mu"]]
  X <- frmtmb:::getME.frmtmb_fit(fit, "X")
  expect_identical(dim(X), c(180L, 2L))
  expect_identical(colnames(X), c("(Intercept)", "Days"))
  expect_equal(unname(as.matrix(X)), unname(as.matrix(lp$X)))

  Z <- frmtmb:::getME.frmtmb_fit(fit, "Z")
  expect_s4_class(Z, "Matrix")
  expect_identical(dim(Z), c(180L, 36L))
  expect_equal(unname(as.matrix(Z)), unname(as.matrix(lp$Z)))
  # level-major within the block, one column per (level, coefficient)
  expect_identical(colnames(Z)[1:3],
                   c("308.(Intercept)", "308.Days", "309.(Intercept)"))
  Zt <- frmtmb:::getME.frmtmb_fit(fit, "Zt")
  expect_identical(dim(Zt), c(36L, 180L))
  expect_equal(unname(as.matrix(Zt)), t(unname(as.matrix(Z))))

  # b lines up with the Z columns, so Z %*% b is the RE contribution
  b <- frmtmb:::getME.frmtmb_fit(fit, "b")
  expect_identical(names(b), colnames(Z))
  re <- t(ranef(fit)[[1]])
  expect_equal(unname(b), as.numeric(re))
})

test_that("getME serves the scalar vocabulary", {
  fit <- getme_fit
  expect_equal(frmtmb:::getME.frmtmb_fit(fit, "beta"),
               fit$estimates$beta)
  expect_identical(frmtmb:::getME.frmtmb_fit(fit, "fixef"),
                   frmtmb:::getME.frmtmb_fit(fit, "beta"))
  th <- frmtmb:::getME.frmtmb_fit(fit, "theta")
  expect_length(th, 3L)
  expect_identical(names(th), c("theta_1", "theta_2", "theta_3"))
  # the internal covariance scale is unbounded, unlike lme4's theta
  expect_identical(frmtmb:::getME.frmtmb_fit(fit, "lower"),
                   stats::setNames(rep(-Inf, 3), names(th)))
  expect_equal(frmtmb:::getME.frmtmb_fit(fit, "sigma"), sigma(fit))
  fl <- frmtmb:::getME.frmtmb_fit(fit, "flist")
  expect_named(fl, "Subject")
  expect_s3_class(fl$Subject, "factor")
  expect_identical(levels(fl$Subject), levels(fit$frame$data_frame$Subject))
  expect_length(fl$Subject, 180L)
  expect_identical(frmtmb:::getME.frmtmb_fit(fit, "n_rtrms"), 1L)
  expect_identical(frmtmb:::getME.frmtmb_fit(fit, "n_rfacs"), 1L)
  # several names come back as a named list
  both <- frmtmb:::getME.frmtmb_fit(fit, c("theta", "sigma"))
  expect_named(both, c("theta", "sigma"))
})

test_that("getME refuses names outside the vocabulary", {
  expect_error(frmtmb:::getME.frmtmb_fit(getme_fit, "Lambdat"),
               "unknown name")
  expect_error(frmtmb:::getME.frmtmb_fit(getme_fit, "Lambdat"),
               "n_rfacs")
})

test_that("getME dispatches from lme4's generic and matches lmer", {
  skip_if_not_installed("lme4")
  fit <- getme_fit
  expect_identical(dim(lme4::getME(fit, "X")), c(180L, 2L))
  ref <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                    data = fit$frame$data_frame, REML = FALSE)
  # lme4 hangs a msgScaleX attribute on its own X, so compare contents
  Xf <- as.matrix(lme4::getME(fit, "X"))
  Xr <- as.matrix(lme4::getME(ref, "X"))
  expect_identical(dim(Xf), dim(Xr))
  expect_equal(as.vector(Xf), as.vector(Xr))
  Zf <- as.matrix(lme4::getME(fit, "Z"))
  Zr <- as.matrix(lme4::getME(ref, "Z"))
  expect_identical(dim(Zf), dim(Zr))
  expect_equal(as.vector(Zf), as.vector(Zr))
  expect_identical(lme4::getME(fit, "n_rtrms"),
                   as.integer(lme4::getME(ref, "n_rtrms")))
  expect_identical(lme4::getME(fit, "n_rfacs"),
                   as.integer(lme4::getME(ref, "n_rfacs")))
})

test_that("getME needs resp= for the designs of a multivariate fit", {
  set.seed(414)
  n <- 60
  dd <- data.frame(x = rnorm(n), g = factor(rep(1:6, 10)))
  dd$y1 <- rnorm(n, 1 + 0.5 * dd$x + rnorm(6, 0, 0.4)[dd$g], 1)
  dd$y2 <- rnorm(n, -1 + 0.3 * dd$x, 1)
  mv <- frm(bf(y1 ~ x + (1 | g)) + gaussian() +
              bf(y2 ~ x) + gaussian(), data = dd)
  expect_error(frmtmb:::getME.frmtmb_fit(mv, "X"), "resp")
  expect_error(frmtmb:::getME.frmtmb_fit(mv, "Zt"), "resp")
  expect_identical(dim(frmtmb:::getME.frmtmb_fit(mv, "X", resp = "y1")),
                   c(60L, 2L))
  # the scalar names answer without one
  expect_length(frmtmb:::getME.frmtmb_fit(mv, "beta"), 4L)
  expect_identical(frmtmb:::getME.frmtmb_fit(mv, "n_rfacs"), 1L)
  expect_length(frmtmb:::getME.frmtmb_fit(mv, "sigma"), 2L)
})
