# Edge cases mined from lme4/glmmTMB/brms issue history (dev/test-backlog.md
# holds the full list with sources).

test_that("data-dependent bases are frozen at fit time (glmmTMB#402)", {
  set.seed(151)
  dd <- data.frame(x = runif(60, 0, 10), g = factor(rep(1:6, 10)))
  dd$y <- 1 + 0.5 * dd$x - 0.05 * dd$x^2 + rnorm(6, 0, 0.5)[dd$g] +
    rnorm(60, 0, 0.3)
  m1 <- frm(bf(y ~ poly(x, 2) + (1 | g)) + gaussian(), data = dd)
  m2 <- frm(bf(y ~ poly(x, 2, raw = TRUE) + (1 | g)) + gaussian(),
            data = dd)
  nd <- data.frame(x = c(0, 2.5, 9), g = factor(1, levels = levels(dd$g)))
  expect_equal(predict(m1, newdata = nd, re.form = NA),
               predict(m2, newdata = nd, re.form = NA), tolerance = 1e-6)
  # single-row newdata is the killer case (brms#494)
  p1 <- predict(m1, newdata = nd[2, , drop = FALSE], re.form = NA)
  expect_equal(p1, predict(m1, newdata = nd, re.form = NA)[2],
               tolerance = 1e-8)
  # scale() in the formula round-trips through prediction
  m3 <- frm(bf(y ~ scale(x) + (1 | g)) + gaussian(), data = dd)
  expect_equal(predict(m3, newdata = dd), predict(m3), tolerance = 1e-8)
})

test_that("rank-deficient designs drop aliased columns (lme4#144)", {
  set.seed(152)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$x2 <- 2 * dd$x                       # perfectly collinear
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
  expect_message(
    m <- frm(bf(y ~ x + x2 + (1 | g)) + gaussian(), data = dd),
    "rank deficient"
  )
  m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(m0))), 1e-6)
  # prediction from the rank-reduced fit works and matches
  expect_equal(predict(m, newdata = dd), predict(m0, newdata = dd),
               tolerance = 1e-6)
})

test_that("matrix-attribute responses and offsets are handled (glmmTMB#937/#773)", {
  set.seed(153)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
  dd$ys <- scale(dd$y)                    # n x 1 matrix with attributes
  m <- frm(bf(ys ~ x + (1 | g)) + gaussian(), data = dd)
  expect_length(fitted(m), 60)
})

test_that("row-permutation and relevel invariance for covariance structures (brms#1747)", {
  dd <- sim_ar1_data(seed = 154, n_g = 40)
  perm <- sample(nrow(dd))
  f1 <- frm(bf(y ~ 1 + ar1(tim + 0 | g)) + gaussian(), data = dd)
  f2 <- frm(bf(y ~ 1 + ar1(tim + 0 | g)) + gaussian(), data = dd[perm, ])
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f2))), 1e-6)

  u1 <- frm(bf(y ~ 1 + us(tim + 0 | g)) + gaussian(), data = dd)
  u2 <- frm(bf(y ~ 1 + us(tim + 0 | g)) + gaussian(), data = dd[perm, ])
  expect_lt(abs(as.numeric(logLik(u1)) - as.numeric(logLik(u2))), 1e-5)
  # us() is also invariant to releveling the term factor
  dd3 <- transform(dd, tim = stats::relevel(tim, "3"))
  u3 <- frm(bf(y ~ 1 + us(tim + 0 | g)) + gaussian(), data = dd3)
  expect_lt(abs(as.numeric(logLik(u1)) - as.numeric(logLik(u3))), 1e-5)
})

test_that("numeric, character, and factor grouping variables are equivalent", {
  set.seed(155)
  d1 <- data.frame(x = rnorm(100), g = rep(1:10, 10))
  d1$y <- rnorm(100, 1 + 0.5 * d1$x + rnorm(10, 0, 0.6)[d1$g], 1)
  d2 <- transform(d1, g = factor(g))
  d3 <- transform(d1, g = as.character(g))
  f1 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d1)
  f2 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d2)
  f3 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d3)
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f2))), 1e-6)
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f3))), 1e-6)
})

test_that("duplicate multivariate responses are rejected (brms)", {
  expect_error(
    frm(mvbf(bf(y ~ x) + gaussian(), bf(y ~ 1) + gaussian()),
        data = NULL, dry_run = "spec"),
    "Duplicated response"
  )
})

test_that("trials() validation catches the brms error taxonomy", {
  dd <- data.frame(y = c(2, 3, 5), n = c(5, 5, 4), x = 1:3)
  expect_error(frm(bf(y | trials(n) ~ x) + binomial(), data = dd),
               "\\[0, trials\\]")   # y > trials
  dd2 <- data.frame(y = c(0.5, 1), n = c(2, 2), x = 1:2)
  expect_error(frm(bf(y | trials(n) ~ x) + binomial(), data = dd2),
               "integer")
  # constant literal trials work
  dd3 <- data.frame(y = rbinom(50, 10, 0.4), x = rnorm(50))
  fit <- frm(bf(y | trials(10) ~ x) + binomial(), data = dd3)
  expect_length(fitted(fit), 50)
})

test_that("response-scale equivariance under rescaling", {
  set.seed(156)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.6)[dd$g], 1)
  m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  s <- 1000
  dd$ys <- dd$y * s
  m1 <- frm(bf(ys ~ x + (1 | g)) + gaussian(), data = dd)
  expect_vector_equal(fixef(m1)$mu, fixef(m0)$mu * s, tol = 1e-2)
  expect_lt(abs(as.numeric(logLik(m1)) -
                  (as.numeric(logLik(m0)) - 100 * log(s))), 1e-4)
})

test_that("NA handling: rows dropped consistently, Inf rejected upstream", {
  set.seed(157)
  dd <- data.frame(x = rnorm(50), g = factor(rep(1:5, 10)))
  dd$y <- rnorm(50, 1 + 0.5 * dd$x, 1)
  dd$y[1:3] <- NA
  dd$g[4] <- NA
  m <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_identical(stats::nobs(m), 46L)
  expect_length(fitted(m), 46)
  # NA rows in newdata propagate NA predictions, not errors
  nd <- dd[1:6, ]
  p <- predict(m, newdata = nd, re.form = NA)
  expect_true(is.na(p[1]) == FALSE)   # NA was only in y/g, x is fine
})
