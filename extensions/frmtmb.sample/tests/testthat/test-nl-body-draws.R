# frm_sample() on a nonlinear body: the shapes frm() dropped or refused
# through core 0.78.0. The sampler reads core's frame and objective, so the
# fix is core's; these hold the sampling route to it. Through 0.78.0 an
# nl body under an ordinal family sampled the thresholds alone
# (dev/nlfix-sample.R, seed 49: variables b_Intercept[1:3] and lp__).

test_that("an nl body under an ordinal family is sampled", {
  skip_on_cran()
  set.seed(49)
  n <- 150
  d <- data.frame(x = rnorm(n))
  d$y <- ordered(cut(0.8 * d$x + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                     labels = 1:4))
  s <- suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d, family = cumulative(),
    chains = 1, iter = 300, warmup = 150, seed = 49, refresh = 0)))
  expect_true("b_b1_Intercept" %in% variables(s))
  m <- as.matrix(s)
  # the slope the data carry, against the spread of its own draws
  expect_gt(mean(m[, "b_b1_Intercept"]) / stats::sd(m[, "b_b1_Intercept"]), 2)
  e <- posterior_epred(s)
  expect_equal(dim(e), c(150L, n, 4L))
  # the predictor moves the category probabilities within a draw
  expect_gt(stats::sd(e[1L, , 1L]), 0)
  expect_equal(dim(posterior_predict(s)), c(150L, n))
  pr <- default_prior(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d,
                      family = cumulative(), route = "sample")
  expect_identical(pr$prior[pr$class == "Intercept" & pr$coef == ""],
                   "student_t(3, 0, 2.5)")
})

test_that("an nl body with a residual term is sampled", {
  skip_on_cran()
  set.seed(158)
  n <- 120
  d <- data.frame(x = rnorm(n), g = factor(sort(rep(1:6, length.out = n))))
  d$t <- stats::ave(seq_len(n), d$g, FUN = seq_along)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(n, 0, 0.5)
  pr <- c(set_prior("normal(1, 2)", nlpar = "b1"),
          set_prior("normal(0, 1)", nlpar = "b2"))
  s <- suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
       autocor = ~ ar(time = t, gr = g, cov = TRUE)),
    data = d, family = gaussian(), prior = pr, chains = 1, iter = 300,
    warmup = 150, seed = 158, refresh = 0)))
  expect_true(all(c("b_b1_Intercept", "b_b2_Intercept", "ar[1]") %in%
                    variables(s)))
  expect_equal(dim(posterior_predict(s)), c(150L, n))
  # cov = FALSE reads the same body, with the shift on earlier residuals
  s2 <- suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
       autocor = ~ ar(time = t, gr = g, cov = FALSE)),
    data = d, family = gaussian(), prior = pr, chains = 1, iter = 300,
    warmup = 150, seed = 158, refresh = 0)))
  expect_true("ar[1]" %in% variables(s2))
  e <- posterior_epred(s2)
  expect_equal(dim(e), c(150L, n))
  # the first row of a group has no earlier residual; the second has one
  m <- as.matrix(s2)
  body1 <- m[, "b_b1_Intercept"] * exp(m[, "b_b2_Intercept"] * d$x[1L])
  expect_lt(max(abs(e[, 1L] - body1)) / max(abs(body1)), 1e-8)
  body2 <- m[, "b_b1_Intercept"] * exp(m[, "b_b2_Intercept"] * d$x[2L])
  by_hand <- body2 + m[, "ar[1]"] * (d$y[1L] - body1)
  expect_lt(max(abs(e[, 2L] - by_hand)) / max(abs(by_hand)), 1e-8)
})

test_that("the guard leaves a linear combination of parameters alone", {
  # its first form stepped along the null direction of 2 * a - b and refused
  # this on frm_sample() too (dev/reviews/2026-11-02-nlfix.md, major 1)
  skip_on_cran()
  set.seed(6201)
  n <- 100
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 + rnorm(n)
  s <- suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ 2 * a - b, a + b ~ 1, nl = TRUE), data = d, family = gaussian(),
    prior = c(set_prior("normal(0, 1)", nlpar = "a"),
              set_prior("normal(0, 1)", nlpar = "b")),
    chains = 1, iter = 200, warmup = 100, seed = 1, refresh = 0)))
  expect_true(all(c("b_a_Intercept", "b_b_Intercept") %in% variables(s)))
})

test_that("nl = TRUE on a family with several locations is sampled", {
  skip_on_cran()
  set.seed(5001)
  n <- 120
  d <- data.frame(x = rnorm(n))
  d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
  s <- suppressWarnings(suppressMessages(frm_sample(
    bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), data = d,
    family = categorical(), chains = 1, iter = 300, warmup = 150, seed = 1,
    refresh = 0,
    prior = c(set_prior("normal(0, 2)", nlpar = "b0"),
              set_prior("normal(0, 2)", nlpar = "b1")))))
  expect_true(all(c("b_b0_Intercept", "b_b1_Intercept") %in% variables(s)))
  expect_equal(dim(posterior_epred(s)), c(150L, n, 3L))
})
