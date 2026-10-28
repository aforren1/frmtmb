test_that("marginaleffects reads sigma as fitted() does under ar(cov = TRUE)", {
  skip_if_not_installed("marginaleffects")
  # the 0.77.0 consolidation: marginaleffects' dpar = (one lane) on the
  # innovation sd (the other): the same estimate and standard error as
  # fitted(dpar = "sigma"); a saved 0.76.0 fit answers on the mean and
  # is refused on sigma with frmtmb's reason
  set.seed(1)
  d <- data.frame(x = stats::rnorm(60), h = factor(rep(1:5, each = 12)),
                  t = rep(1:12, 5))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(60)
  nd <- data.frame(x = c(-1, 0, 1), t = 1, h = factor(1, levels = 1:5))
  fs <- frm(bf(y ~ x + ar(time = t, gr = h, cov = TRUE), sigma ~ x),
            data = d)
  me <- allow_warnings(marginaleffects::predictions(fs, newdata = nd,
                                                    dpar = "sigma"),
                       "not known to be supported")
  ft <- fitted(fs, newdata = nd, dpar = "sigma")
  expect_equal(me$estimate, unname(ft[, "Estimate"]), tolerance = 1e-10)
  expect_equal(me$std.error, unname(ft[, "Est.Error"]), tolerance = 1e-6)
  fit <- frm(y ~ x + ar(time = t, gr = h, cov = TRUE), data = d)
  m1 <- allow_warnings(marginaleffects::predictions(fit, newdata = nd[1, ],
                                                    dpar = "sigma"),
                       "not known to be supported")
  expect_equal(m1$estimate, unname(sigma(fit)), tolerance = 1e-12)
  old <- fit
  old$frame$autocor[[1L]]$sigma_is <- NULL
  expect_identical(
    marginaleffects::predictions(old, newdata = nd)$estimate,
    marginaleffects::predictions(fit, newdata = nd)$estimate)
  expect_error(marginaleffects::get_predict(old, newdata = nd,
                                            dpar = "sigma"),
               "frmtmb 0.76.0 or earlier", class = "frmtmb_error")
})
