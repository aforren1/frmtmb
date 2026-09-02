# Regression pins for brms vignette code that ports on the mechanical
# brm -> frm transform. Each block is the vignette's own call with the
# MCMC arguments removed and nothing else changed, so a future grammar or
# default change that breaks a ported vignette breaks a test here.
#
# The full audit is dev/brms-vignette-port.md; this file is the subset
# that is cheap enough to run every time.

test_that("brms_overview: kidney lognormal with censoring ports verbatim", {
  skip_on_cran()
  skip_if_not_installed("brms")
  utils::data("kidney", package = "brms", envir = environment())

  # the patient variance sits on the boundary under ML, so nlminb chatters
  # about NA/NaN function evaluations on its way there
  fit1 <- suppressWarnings(
    frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
        data = kidney, family = lognormal())
  )

  fx <- fixef(fit1)$mu
  expect_equal(unname(fx[["(Intercept)"]]), 2.6807, tolerance = 1e-3)
  expect_equal(unname(fx[["sexfemale"]]), 2.4730, tolerance = 1e-3)
  expect_equal(as.numeric(sigma(fit1)), 1.1209, tolerance = 1e-3)
  # the vignette's brms posterior means, for orientation: 2.73 and 2.42
  expect_lt(abs(fx[["(Intercept)"]] - 2.73), 0.2)
  expect_lt(abs(fx[["sexfemale"]] - 2.42), 0.2)
})

test_that("brms_overview: inhaler ordinal ports once the family is called", {
  skip_on_cran()
  skip_if_not_installed("brms")
  utils::data("inhaler", package = "brms", envir = environment())

  # brms writes family = cumulative; frm() needs the constructor called
  fit3 <- frm(rating ~ treat + period + carry + (1 | subject),
              data = inhaler, family = cumulative())

  fx <- fixef(fit3)$mu
  expect_named(fx, c("treat", "period", "carry"))
  expect_equal(unname(fx[["treat"]]), -0.9978, tolerance = 1e-3)
  expect_equal(as.numeric(logLik(fit3)), -448.29, tolerance = 1e-3)
})

test_that("brms_nonlinear: the loss growth curve ports with start values", {
  skip_on_cran()
  skip_if_not_installed("brms")
  utils::data("loss", package = "brms", envir = environment())

  # brms gets the starting region from its priors; frmtmb needs `start`
  fit_loss <- frm(bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
                     ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE),
                  data = loss, family = gaussian(),
                  start = list(beta = c(5000, 1, 45)))

  expect_equal(unname(fixef(fit_loss)$ult[["(Intercept)"]]), 5292.5,
               tolerance = 1e-3)
  expect_equal(unname(fixef(fit_loss)$omega[["(Intercept)"]]), 1.3370,
               tolerance = 1e-3)
  expect_equal(unname(fixef(fit_loss)$theta[["(Intercept)"]]), 45.899,
               tolerance = 1e-3)
})

test_that("brms_customfamilies: cbpp binomial with trials() ports verbatim", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  utils::data("cbpp", package = "lme4", envir = environment())

  fit1 <- frm(incidence | trials(size) ~ period + (1 | herd),
              data = cbpp, family = binomial())

  fx <- fixef(fit1)$mu
  expect_equal(unname(fx[["(Intercept)"]]), -1.3985, tolerance = 1e-3)
  expect_true(all(fx[c("period2", "period3", "period4")] < 0))
})

test_that("brms_distreg: a sigma submodel and its hypothesis port verbatim", {
  skip_on_cran()
  set.seed(1234)
  group <- rep(c("treat", "placebo"), each = 30)
  symptom_post <- c(rnorm(30, mean = 1, sd = 2), rnorm(30, mean = 0, sd = 1))
  dat1 <- data.frame(group, symptom_post)

  fit1 <- frm(bf(symptom_post ~ group, sigma ~ group), data = dat1,
              family = gaussian())

  expect_named(fixef(fit1), c("mu", "sigma"))
  # the vignette's own two-sided hypotheses, unchanged
  hyp <- hypothesis(fit1, c("exp(sigma_Intercept) = 0",
                            "exp(sigma_Intercept + sigma_grouptreat) = 0"))
  expect_equal(nrow(as.data.frame(hyp)), 2L)
  # the treated group is the more variable one, as simulated
  expect_gt(unname(fixef(fit1)$sigma[["grouptreat"]]), 0)
})

test_that("brms_monotonic: mo() and its interaction port once a family is given", {
  skip_on_cran()
  set.seed(1234)
  income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(income_options, 100, TRUE),
                   levels = income_options, ordered = TRUE)
  mean_ls <- c(30, 60, 70, 75)
  ls <- mean_ls[income] + rnorm(100, sd = 7)
  dat <- data.frame(income, ls)
  dat$age <- rnorm(100, mean = 40, sd = 10)

  fit1 <- frm(ls ~ mo(income), data = dat, family = gaussian())
  fit5 <- frm(ls ~ mo(income) * age, data = dat, family = gaussian())

  # brms's mo() coefficient is the average step, so three steps span the
  # simulated 30 -> 75 range; the vignette prints moincome = 15.73
  b <- unname(fixef(fit1)$mu[["moincome"]])
  expect_equal(3 * b, 45, tolerance = 0.1)
  expect_s3_class(conditional_effects(fit5, "income:age"),
                  "frmtmb_conditional_effects")
})
