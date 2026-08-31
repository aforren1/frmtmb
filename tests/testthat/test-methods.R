fit_sleep <- local({
  if (requireNamespace("lme4", quietly = TRUE)) {
    data(sleepstudy, package = "lme4")
    frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
           data = sleepstudy)
  } else {
    NULL
  }
})

test_that("accessor methods are consistent", {
  skip_if(is.null(fit_sleep))
  fit <- fit_sleep

  expect_identical(stats::nobs(fit), 180L)
  expect_s3_class(logLik(fit), "logLik")
  expect_identical(attr(logLik(fit), "df"), length(fit$opt$par))

  fe <- fixef(fit)
  expect_named(fe, c("mu", "sigma"))
  expect_named(fe$mu, c("(Intercept)", "Days"))

  V <- vcov(fit)
  expect_identical(dim(V), c(3L, 3L))
  expect_true(isSymmetric(V, tol = 1e-8))
  expect_true(all(diag(V) > 0))

  re <- ranef(fit)
  expect_length(re, 1)
  expect_identical(dim(re[[1]]), c(18L, 2L))

  vc <- VarCorr(fit)
  expect_length(vc, 1)
  expect_identical(dim(vc[[1]]), c(2L, 2L))

  expect_identical(family(fit)$family, "gaussian")
  expect_s3_class(formula(fit), "formula")
})

test_that("predict/fitted/residuals invariants hold", {
  skip_if(is.null(fit_sleep))
  fit <- fit_sleep

  mu <- predict(fit, type = "response")
  expect_length(mu, stats::nobs(fit))
  expect_identical(mu, fitted(fit))
  expect_identical(predict(fit, type = "link"), mu)  # identity link

  sig <- predict(fit, type = "response", dpar = "sigma")
  expect_true(all(sig > 0))
  expect_lt(stats::sd(sig), 1e-10)  # intercept-only sigma is constant

  r <- residuals(fit)
  expect_equal(r, fit$frame$y$Reaction - fitted(fit))
  rp <- residuals(fit, type = "pearson")
  expect_equal(rp, r / sig, tolerance = 1e-10)
})

test_that("print and summary run without error", {
  skip_if(is.null(fit_sleep))
  expect_output(print(fit_sleep), "frmtmb fit")
  s <- summary(fit_sleep)
  expect_s3_class(s, "summary.frmtmb_fit")
  expect_output(print(s), "Coefficients")
  expect_named(s$coefficients, c("mu", "sigma"))
})

test_that("start values are validated", {
  skip_if(is.null(fit_sleep))
  data(sleepstudy, package = "lme4")
  expect_error(frm(bf(Reaction ~ Days) + gaussian(), sleepstudy,
                      start = list(bogus = 1)),
               "Unknown start component")
  expect_error(frm(bf(Reaction ~ Days) + gaussian(), sleepstudy,
                      start = list(beta = 1)),
               "length")
})
