test_that("insight::get_sigma() keeps insight's default on a current fit", {
  skip_if_not_installed("insight")
  # the 0.77.0 consolidation: the method refuses only a fit saved with
  # the old sigma and hands every other fit to insight's default. Under
  # sigma ~ x, sigma() is NA, and insight falls back on
  # sqrt(|deviance| / residual df), as for any model (1.7285282 at
  # frmtmb 0.76.0 on these data; NA when the method returned sigma())
  set.seed(3)
  d <- data.frame(x = stats::rnorm(120), x2 = stats::rnorm(120))
  d$y <- 1 + d$x + stats::rnorm(120, sd = exp(0.2 * d$x2))
  fit <- frm(bf(y ~ x, sigma ~ x2), data = d)
  allow_warnings({
    gs <- as.numeric(insight::get_sigma(fit))
    ref <- sqrt(abs(insight::get_deviance(fit, verbose = FALSE)) /
                  insight::get_df(fit, type = "residual", verbose = FALSE))
  }, "sigma varies by observation")
  expect_true(is.finite(gs))
  expect_identical(gs, as.numeric(ref))
  # where sigma() answers, insight's default is sigma() itself
  fg <- frm(y ~ x, data = d)
  expect_identical(as.numeric(insight::get_sigma(fg)), unname(sigma(fg)))
})
