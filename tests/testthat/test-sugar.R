# v0.13 sugar: conventional accessor methods, coef convention, refit.

test_that("accessors match lme4 conventions on sleepstudy", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
             data = sleepstudy)
  ref <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                    REML = FALSE)

  expect_equal(sigma(fit), sigma(ref), tolerance = 1e-4)
  # brms's named list since item 2.6f, not lme4's integer vector
  expect_equal(ngrps(fit), list(Subject = 18L))
  expect_equal(weights(fit), rep(1, 180))
  expect_equal(deviance(fit), -2 * as.numeric(logLik(fit)))

  ea <- extractAIC(fit)
  expect_equal(ea[1], attr(logLik(fit), "df"))
  expect_equal(ea[2], stats::AIC(fit))

  X <- model.matrix(fit)
  expect_equal(dim(X), c(180L, 2L))
  expect_equal(colnames(X), c("(Intercept)", "Days"))
  expect_s3_class(terms(fit), "terms")

  # conventional coef: fixef broadcast + conditional modes, per level
  cc <- coef(fit)
  rc <- coef(ref)$Subject
  expect_equal(dim(cc$Subject), dim(rc))
  expect_lt(max(abs(as.matrix(cc$Subject) - as.matrix(rc))), 0.5)
  expect_equal(rownames(cc$Subject), rownames(rc))
  # coef = fixef + ranef exactly, within our own fit
  expect_equal(cc$Subject$Days,
               unname(fixef_by_dpar(fit)$mu["Days"] + ranef(fit)[[1]][,
                                                                      "Days"]),
               tolerance = 1e-10)
})

test_that("coef falls back to fixef without random effects", {
  dd <- data.frame(y = rnorm(50), x = rnorm(50))
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  cc <- coef(fit)
  expect_type(cc, "double")
  # the VALUES are fixef_by_dpar()'s; the NAMES are brms's, the same
  # ones fixef() and vcov() put on their rows, so that anything pairing
  # coef() with vcov() by name keeps every row
  expect_equal(unname(cc), unname(fixef_by_dpar(fit)$mu))
  expect_named(cc, rownames(vcov(fit)))
  expect_named(cc, c("Intercept", "x"))
})

test_that("two terms on one factor share a coef frame", {
  dd <- sim_pois_glmm()
  fit <- frm(bf(y ~ x + (1 | g) + (0 + x | g)) + poisson(), data = dd)
  cc <- coef(fit)
  expect_named(cc, "g")
  expect_equal(colnames(cc$g), c("Intercept", "x"))
  expect_equal(cc$g$x,
               unname(fixef_by_dpar(fit)$mu["x"] + ranef(fit)[[2]][, "x"]),
               tolerance = 1e-10)
})

test_that("sigma handles non-gaussian and modeled-sigma fits", {
  dd <- sim_pois_glmm()
  fp <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  expect_equal(sigma(fp), 1)

  set.seed(7)
  dg <- data.frame(x = rnorm(200))
  dg$y <- rnorm(200, 1 + dg$x, exp(0.2 + 0.3 * dg$x))
  fg <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dg)
  expect_warning(s <- sigma(fg), "varies by observation")
  expect_true(is.na(s))
})

test_that("prior_summary returns the fit priors", {
  dd <- data.frame(y = rnorm(60), x = rnorm(60))
  f0 <- frm(bf(y ~ x) + gaussian(), data = dd)
  expect_output(expect_null(prior_summary(f0)), "No priors")
  pr <- set_prior("normal(0, 1)", class = "b")
  f1 <- frm(bf(y ~ x) + gaussian(), data = dd, prior = pr)
  expect_s3_class(prior_summary(f1), "frmtmb_priorlist")
})

test_that("refit matches a fresh fit on the new response", {
  dd <- sim_pois_glmm(n_g = 20, n_per = 10)
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  set.seed(99)
  ysim <- simulate(fit, nsim = 1, re_formula = NA)[[1L]]

  rf <- refit(fit, ysim)
  dd2 <- dd
  dd2$y <- ysim
  fresh <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd2)

  expect_loglik_equal(rf, fresh, tol = 1e-6)
  expect_vector_equal(fixef_by_dpar(rf)$mu, fixef_by_dpar(fresh)$mu, tol = 1e-5)
  # the refit is a full frmtmb_fit: methods work
  expect_s3_class(summary(rf), "summary.frmtmb_fit")
  expect_error(refit(fit, ysim[-1]), "length")
})

test_that("insight defaults work off the standard accessors", {
  skip_if_not_installed("insight")
  set.seed(11)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.7)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  expect_equal(insight::find_response(fit), "y")
  expect_equal(insight::n_obs(fit), 100)
  expect_equal(nrow(insight::get_data(fit)), 100)
  expect_equal(as.numeric(insight::get_sigma(fit)), sigma(fit))
  # registered insight methods (v0.17) split fixed from random parts
  expect_equal(insight::find_predictors(fit)$conditional, "x")
  expect_equal(insight::find_random(fit)$random, "g")
})

test_that("refit powers a small parametric bootstrap", {
  set.seed(31)
  dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
  dd$y <- rnorm(120, 1 + 0.5 * dd$x + rnorm(12, 0, 0.7)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  sims <- simulate(fit, nsim = 5, re_formula = NA, seed = 1)
  boots <- vapply(sims, function(ys) fixef_by_dpar(refit(fit, ys))$mu["x"],
                  numeric(1))
  expect_length(boots, 5)
  expect_true(all(is.finite(boots)))
  expect_lt(abs(mean(boots) - fixef_by_dpar(fit)$mu["x"]), 0.3)
})

test_that("insight still falls back to NA when model.matrix() refuses", {
  skip_if_not_installed("insight")
  # insight's .get_predicted_ci_modelmatrix() tests the error it catches
  # for "simpleError"; a multivariate fit is where model.matrix() refuses
  set.seed(20260917)
  d <- data.frame(y = rnorm(60), z = rnorm(60), x = rnorm(60))
  mv <- frm(mvbf(bf(y ~ x), bf(z ~ x)), data = d)
  e <- tryCatch(model.matrix(mv), error = identity)
  expect_identical(class(e), c("frmtmb_error", "simpleError", "error",
                               "condition"))
  expect_match(conditionMessage(e), "disambiguate with resp")
  ci <- NULL
  expect_warning(
    ci <- insight::get_predicted_ci(mv, predictions = rep(0, 60), data = d),
    "Something went wrong")
  expect_identical(dim(ci), c(60L, 3L))
  expect_true(all(is.na(ci$SE)))
})
