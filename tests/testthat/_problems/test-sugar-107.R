# Extracted from test-sugar.R:107

# test -------------------------------------------------------------------------
skip_if_not_installed("insight")
set.seed(11)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.7)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
expect_equal(insight::find_response(fit), "y")
expect_equal(insight::n_obs(fit), 100)
expect_equal(nrow(insight::get_data(fit)), 100)
expect_equal(as.numeric(insight::get_sigma(fit)), sigma(fit))
expect_setequal(insight::find_predictors(fit)$conditional, c("x", "g"))
