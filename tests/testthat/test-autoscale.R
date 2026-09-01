# frmtmb_control(autoscale = TRUE): internal predictor standardization
# via the two-stage warm start (scaled pre-fit, exact back-transform,
# ordinary unscaled fit). ML is invariant under affine column
# reparameterization, so hand-standardized references are exact.

# Collect a fit together with every warning it raised.
fit_catching <- function(expr) {
  w <- character(0)
  fit <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warnings = w)
}

sim_badscale_glmm <- function(seed = 202, n_g = 30, n_per = 20) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_g), each = n_per))
  x1 <- stats::rnorm(n_g * n_per) * 1e6
  x2 <- stats::rnorm(n_g * n_per) * 1e-6
  b0 <- stats::rnorm(n_g, 0, 0.5)
  eta <- 0.3 + 4e-7 * x1 + 3e5 * x2 + b0[g]
  data.frame(y = stats::rpois(length(eta), exp(eta)),
             x1 = x1, x2 = x2, g = g)
}

test_that("autoscale reproduces the plain fit on a well-scaled model", {
  dd <- sim_pois_glmm()
  form <- bf(y ~ x + (1 | g)) + poisson()
  f0 <- frm(form, data = dd)
  fa <- frm(form, data = dd,
            control = frmtmb_control(autoscale = TRUE))
  expect_loglik_equal(fa, f0, tol = 1e-8)
  expect_vector_equal(fixef(fa)$mu, fixef(f0)$mu, tol = 1e-6)
  expect_vector_equal(fa$estimates$theta, f0$estimates$theta, tol = 1e-5)
})

test_that("autoscale is a silent no-op when nothing qualifies", {
  set.seed(303)
  dd <- data.frame(f = factor(rep(letters[1:3], 100)),
                   g = factor(rep(1:30, each = 10)))
  b0 <- stats::rnorm(30, 0, 0.5)
  dd$y <- stats::rpois(300, exp(0.5 + c(0, 0.3, -0.2)[dd$f] + b0[dd$g]))
  form <- bf(y ~ f + (1 | g)) + poisson()
  f0 <- frm(form, data = dd)
  expect_null(autoscale_plan(f0$frame))
  fa <- NULL
  expect_silent(fa <- frm(form, data = dd,
                          control = frmtmb_control(autoscale = TRUE)))
  # no qualifying column: exactly the single ordinary fit
  expect_identical(fa$opt$par, f0$opt$par)
  expect_null(fa$par_units)
})

test_that("autoscale rescues a badly scaled poisson GLMM", {
  db <- sim_badscale_glmm()
  form <- bf(y ~ x1 + x2 + (1 | g)) + poisson()

  auto <- fit_catching(frm(form, data = db,
                           control = frmtmb_control(autoscale = TRUE)))
  expect_s3_class(auto$fit, "frmtmb_fit")
  expect_identical(auto$warnings, character(0))
  expect_equal(auto$fit$opt$convergence, 0L)
  fa <- auto$fit

  # hand-standardized reference with the standardization undone by
  # hand: exact under ML reparameterization invariance
  m1 <- mean(db$x1); s1 <- stats::sd(db$x1)
  m2 <- mean(db$x2); s2 <- stats::sd(db$x2)
  db2 <- db
  db2$z1 <- (db$x1 - m1) / s1
  db2$z2 <- (db$x2 - m2) / s2
  fr <- frm(bf(y ~ z1 + z2 + (1 | g)) + poisson(), data = db2)
  expect_loglik_equal(fa, fr, tol = 1e-6)
  cz <- fixef(fr)$mu
  ref <- c(cz[["(Intercept)"]] - cz[["z1"]] * m1 / s1 -
             cz[["z2"]] * m2 / s2,
           cz[["z1"]] / s1, cz[["z2"]] / s2)
  expect_equal(unname(fixef(fa)$mu), ref, tolerance = 1e-6)
  expect_vector_equal(fa$estimates$theta, fr$estimates$theta, tol = 1e-4)

  # SEs transform by the same column scales (finite-difference Hessians
  # on both sides, so agreement is at step-size accuracy)
  se_a <- summary(fa)$coefficients$mu[, "Std. Error"]
  se_r <- summary(fr)$coefficients$mu[, "Std. Error"]
  expect_equal(unname(se_a[c("x1", "x2")]),
               unname(se_r[c("z1", "z2")]) / c(s1, s2),
               tolerance = 1e-3)
  # prediction SEs ride on the joint precision and match too
  pa <- predict(fa, se.fit = TRUE)
  pr <- predict(fr, se.fit = TRUE)
  expect_true(all(is.finite(pa$se.fit)))
  expect_vector_equal(pa$se.fit, pr$se.fit, tol = 1e-3)

  # the plain fit either matches or (as here, typically) fails where
  # autoscale did not
  plain <- fit_catching(frm(form, data = db))
  if (!inherits(plain$fit, "error") && !length(plain$warnings)) {
    expect_loglik_equal(plain$fit, fa, tol = 1e-6)
  } else {
    expect_gt(length(plain$warnings), 0)
  }

  # warm-started refit skips the pre-fit (template short-circuit)
  rf <- refit(fa, db$y)
  expect_loglik_equal(rf, fa, tol = 1e-6)
})

test_that("autoscale leaves smooth and mo() columns untouched", {
  skip_if_not_installed("mgcv")
  set.seed(404)
  n <- 300
  g <- factor(rep(1:30, each = 10))
  xbig <- stats::rnorm(n) * 1e6
  xs <- stats::runif(n)
  m <- factor(sample(1:4, n, replace = TRUE), ordered = TRUE)
  b0 <- stats::rnorm(30, 0, 0.4)
  mo_eff <- c(0, 0.3, 0.5, 0.9)[as.integer(m)]
  dd <- data.frame(xbig = xbig, xs = xs, m = m, g = g)
  dd$y <- 1 + 4e-7 * xbig + sin(2 * pi * xs) + mo_eff + b0[g] +
    stats::rnorm(n, 0, 0.5)

  form <- bf(y ~ xbig + mo(m) + s(xs) + (1 | g)) + gaussian()
  auto <- fit_catching(frm(form, data = dd,
                           control = frmtmb_control(autoscale = TRUE)))
  expect_identical(auto$warnings, character(0))
  fa <- auto$fit

  # only the parametric xbig column is planned; the smooth basis and
  # the mo() placeholder stay untouched in the scaled frame
  plan <- autoscale_plan(fa$frame)
  expect_named(plan, "y.mu")
  lp <- fa$frame$linpreds$y.mu
  expect_identical(plan$y.mu$cols, match("xbig", colnames(lp$X)))
  sf <- autoscale_frame(fa$frame, plan)
  other <- setdiff(seq_len(ncol(lp$X)), plan$y.mu$cols)
  expect_identical(sf$linpreds$y.mu$X[, other], lp$X[, other])

  # hand-standardized reference
  mb <- mean(dd$xbig); sb <- stats::sd(dd$xbig)
  dd2 <- dd
  dd2$xz <- (dd$xbig - mb) / sb
  fr <- frm(bf(y ~ xz + mo(m) + s(xs) + (1 | g)) + gaussian(),
            data = dd2)
  expect_loglik_equal(fa, fr, tol = 1e-6)
  ca <- fixef(fa)$mu
  cr <- fixef(fr)$mu
  expect_equal(ca[["xbig"]], cr[["xz"]] / sb, tolerance = 1e-6)
  expect_equal(ca[["(Intercept)"]],
               cr[["(Intercept)"]] - cr[["xz"]] * mb / sb,
               tolerance = 1e-6)
  expect_equal(unname(ca[["mom"]]), unname(cr[["mom"]]),
               tolerance = 1e-4)
  sm <- grep("s(xs)", names(ca), fixed = TRUE)
  expect_vector_equal(ca[sm], cr[sm], tol = 1e-4)
  expect_vector_equal(fa$estimates$theta, fr$estimates$theta, tol = 1e-4)
})

test_that("autoscale scales without centering when there is no intercept", {
  db <- sim_badscale_glmm(seed = 505)
  form <- bf(y ~ 0 + x1 + x2 + (1 | g)) + poisson()
  auto <- fit_catching(frm(form, data = db,
                           control = frmtmb_control(autoscale = TRUE)))
  expect_identical(auto$warnings, character(0))
  fa <- auto$fit
  plan <- autoscale_plan(fa$frame)
  expect_identical(plan$y.mu$center, c(0, 0))

  s1 <- stats::sd(db$x1); s2 <- stats::sd(db$x2)
  db2 <- db
  db2$z1 <- db$x1 / s1
  db2$z2 <- db$x2 / s2
  fr <- frm(bf(y ~ 0 + z1 + z2 + (1 | g)) + poisson(), data = db2)
  expect_loglik_equal(fa, fr, tol = 1e-6)
  cz <- fixef(fr)$mu
  expect_equal(unname(fixef(fa)$mu),
               c(cz[["z1"]] / s1, cz[["z2"]] / s2),
               tolerance = 1e-6)
})

test_that("autoscale covers dpar formulas and combines with profile", {
  set.seed(606)
  n <- 300
  g <- factor(rep(1:30, each = 10))
  xbig <- stats::rnorm(n) * 1e5
  b0 <- stats::rnorm(30, 0, 0.4)
  dd <- data.frame(xbig = xbig, g = g)
  dd$y <- stats::rnorm(n, 1 + 5e-6 * xbig + b0[g],
                       exp(-0.5 + 2e-6 * xbig))
  form <- bf(y ~ xbig + (1 | g), sigma ~ xbig) + gaussian()
  auto <- fit_catching(frm(form, data = dd,
                           control = frmtmb_control(autoscale = TRUE)))
  expect_identical(auto$warnings, character(0))
  fa <- auto$fit
  # both the beta (mu) and betad (sigma) columns are planned
  plan <- autoscale_plan(fa$frame)
  expect_setequal(vapply(plan, `[[`, "", "par"), c("beta", "betad"))

  mb <- mean(dd$xbig); sb <- stats::sd(dd$xbig)
  dd2 <- dd
  dd2$xz <- (dd$xbig - mb) / sb
  fr <- frm(bf(y ~ xz + (1 | g), sigma ~ xz) + gaussian(), data = dd2)
  expect_loglik_equal(fa, fr, tol = 1e-6)
  expect_equal(fixef(fa)$mu[["xbig"]], fixef(fr)$mu[["xz"]] / sb,
               tolerance = 1e-6)
  expect_equal(fixef(fa)$sigma[["xbig"]], fixef(fr)$sigma[["xz"]] / sb,
               tolerance = 1e-6)

  # profile = TRUE moves beta into the inner problem; autoscale rides
  # on top (both stages profiled)
  db <- sim_badscale_glmm(seed = 707)
  fp <- fit_catching(frm(bf(y ~ x1 + x2 + (1 | g)) + poisson(),
                         data = db,
                         control = frmtmb_control(autoscale = TRUE,
                                                  profile = TRUE)))
  expect_identical(fp$warnings, character(0))
  m1 <- mean(db$x1); s1 <- stats::sd(db$x1)
  m2 <- mean(db$x2); s2 <- stats::sd(db$x2)
  db2 <- db
  db2$z1 <- (db$x1 - m1) / s1
  db2$z2 <- (db$x2 - m2) / s2
  frp <- frm(bf(y ~ z1 + z2 + (1 | g)) + poisson(), data = db2,
             control = frmtmb_control(profile = TRUE))
  expect_loglik_equal(fp$fit, frp, tol = 1e-6)
  expect_equal(fixef(fp$fit)$mu[["x1"]], fixef(frp)$mu[["z1"]] / s1,
               tolerance = 1e-5)
})
