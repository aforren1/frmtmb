test_that("log_lik(newdata =) names only an alternative that answers", {
  # the 0.76.0 release review, finding 5: the refusal named
  # posterior_predict(newdata =), which a correlated residual refuses
  d <- lla_data()
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective")
  ds <- lla_draws(fit, list(b_Intercept = b0, b_x = b1, sigma = sg,
                            thetaac_1 = lla_pinv(c(0.5, -0.3, 0.8))))
  e <- tryCatch(log_lik(ds, newdata = d), error = identity)
  expect_s3_class(e, "frmtmb_error")
  expect_match(conditionMessage(e), "does not take newdata", fixed = TRUE)
  expect_match(conditionMessage(e), "Use posterior_epred(newdata =) for",
               fixed = TRUE)
  expect_match(conditionMessage(e), "is refused for this model too",
               fixed = TRUE)
  # the alternative it rules out is refused indeed
  pe <- tryCatch(posterior_predict(ds, newdata = d), error = identity)
  expect_s3_class(pe, "error")
  expect_match(conditionMessage(pe), "its draws are structured",
               fixed = TRUE)
  # the guard absent: without the residual term both are named
  f2 <- frm(y ~ x, data = d, dry_run = "objective")
  ds2 <- lla_draws(f2, list(b_Intercept = b0, b_x = b1, sigma = sg))
  e2 <- tryCatch(log_lik(ds2, newdata = d), error = identity)
  expect_s3_class(e2, "frmtmb_error")
  expect_match(conditionMessage(e2),
               "Use posterior_epred(newdata =) or posterior_predict(",
               fixed = TRUE)
})
