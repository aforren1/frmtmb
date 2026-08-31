# v0.13 sugar, part 2: conditional_effects, plot, pp_check, profile,
# hypothesis.

test_that("conditional_effects builds sensible grids with Wald bands", {
  set.seed(21)
  n <- 200
  dd <- data.frame(x = rnorm(n), f = factor(rep(c("a", "b"), n / 2)),
                   g = factor(rep(1:10, n / 10)))
  dd$y <- rnorm(n, 1 + 0.8 * dd$x + (dd$f == "b") * 0.5 +
                  rnorm(10, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + f + (1 | g)) + gaussian(), data = dd)

  ce <- conditional_effects(fit)
  expect_s3_class(ce, "frmtmb_conditional_effects")
  expect_setequal(names(ce), c("x", "f"))

  dx <- ce$x
  expect_equal(nrow(dx), 100L)
  expect_true(all(dx$lower__ <= dx$estimate__ &
                    dx$estimate__ <= dx$upper__))
  # the varied grid spans the data; other predictors held at reference
  expect_equal(range(dx$x), range(dd$x))
  expect_true(all(dx$f == "a"))
  # a population prediction at the same point matches
  p <- predict(fit, newdata = data.frame(x = dx$x[1], f = "a"),
               re.form = NA, se.fit = TRUE)
  expect_equal(dx$estimate__[1], unname(p$fit), tolerance = 1e-8)
  expect_equal(dx$se__[1], unname(p$se.fit), tolerance = 1e-8)

  df <- ce$f
  expect_equal(nrow(df), 2L)
  expect_equal(as.character(df$f), c("a", "b"))

  # two-variable effect: second variable at mean and mean +/- sd
  ce2 <- conditional_effects(fit, effects = "x:f")
  expect_equal(nrow(ce2$`x:f`), 200L)
  expect_equal(sort(unique(as.character(ce2$`x:f`$f))), c("a", "b"))

  # conditions override the reference value
  ce3 <- conditional_effects(fit, effects = "x",
                             conditions = list(f = "b"))
  expect_true(all(ce3$x$f == "b"))
  expect_true(all(ce3$x$estimate__ > ce$x$estimate__))
})

test_that("conditional_effects respects the link scale", {
  dd <- sim_pois_glmm()
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  ce <- conditional_effects(fit, effects = "x")
  expect_true(all(ce$x$lower__ > 0))
  eta <- predict(fit, newdata = data.frame(x = ce$x$x), re.form = NA)
  expect_equal(ce$x$estimate__, unname(exp(eta)), tolerance = 1e-8)
})

test_that("plot methods draw without error", {
  dd <- sim_pois_glmm(n_g = 10, n_per = 10)
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  tmp <- file.path(tempdir(), "frmtmb-plots.pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  expect_no_error(plot(fit))
  expect_no_error(plot(conditional_effects(fit)))
  expect_no_error(print(conditional_effects(fit, effects = "x")))
})

test_that("pp_check hands simulations to bayesplot", {
  skip_if_not_installed("bayesplot")
  dd <- sim_pois_glmm(n_g = 10, n_per = 10)
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  p <- bayesplot::pp_check(fit, ndraws = 3)
  expect_s3_class(p, "ggplot")
  p2 <- bayesplot::pp_check(fit, type = "stat", ndraws = 3, stat = "mean")
  expect_s3_class(p2, "ggplot")
})

test_that("profile wraps tmbprofile and agrees with confint", {
  dd <- data.frame(y = rnorm(80, 1, 2), x = rnorm(80))
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  pr <- profile(fit, parm = "x")
  expect_s3_class(pr, "data.frame")
  ci_p <- confint(pr)
  ci_c <- confint(fit, parm = "x", method = "profile")
  expect_equal(unname(ci_p[1, ]), unname(ci_c[1, 1:2]), tolerance = 1e-6)
  prs <- profile(fit, parm = c("(Intercept)", "x"))
  expect_named(prs, c("(Intercept)", "x"))
})

test_that("hypothesis reproduces Wald results and the delta method", {
  set.seed(5)
  dd <- data.frame(x1 = rnorm(150), x2 = rnorm(150))
  dd$y <- rnorm(150, 1 + 0.6 * dd$x1 + 0.4 * dd$x2, 1)
  fit <- frm(bf(y ~ x1 + x2) + gaussian(), data = dd)

  # single coefficient: matches summary's z-test
  h1 <- hypothesis(fit, "x1")
  sm <- summary(fit)$coefficients$mu
  expect_equal(h1$estimate, sm["x1", "Estimate"], tolerance = 1e-10)
  expect_equal(h1$se, sm["x1", "Std. Error"], tolerance = 1e-8)

  # linear contrast: matches direct computation from vcov
  h2 <- hypothesis(fit, "x1 - x2 = 0")
  V <- vcov(fit)
  se_ref <- sqrt(V["x1", "x1"] + V["x2", "x2"] - 2 * V["x1", "x2"])
  expect_equal(h2$se, se_ref, tolerance = 1e-6)

  # nonlinear expression with a dpar coefficient: exp(log sigma) = sigma
  h3 <- hypothesis(fit, "exp(sigma_Intercept)")
  expect_equal(h3$estimate, sigma(fit), tolerance = 1e-8)

  # multiple hypotheses come back as rows
  hh <- hypothesis(fit, c("x1", "x2", "x1 + x2 = 1"))
  expect_equal(nrow(hh), 3L)
  expect_error(hypothesis(fit, "x1 = 0 = 1"), "at most one")
})
