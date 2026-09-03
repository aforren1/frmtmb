# v0.13 sugar: conventional accessor methods, coef convention, refit.

test_that("accessors match lme4 conventions on sleepstudy", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
             data = sleepstudy)
  ref <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                    REML = FALSE)

  expect_equal(sigma(fit), sigma(ref), tolerance = 1e-4)
  expect_equal(ngrps(fit), c(Subject = 18L))
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
               unname(fixef(fit)$mu["Days"] + ranef(fit)[[1]][, "Days"]),
               tolerance = 1e-10)
})

test_that("coef falls back to fixef without random effects", {
  dd <- data.frame(y = rnorm(50), x = rnorm(50))
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  cc <- coef(fit)
  expect_type(cc, "double")
  expect_equal(cc, fixef(fit)$mu)
})

test_that("two terms on one factor share a coef frame", {
  dd <- sim_pois_glmm()
  fit <- frm(bf(y ~ x + (1 | g) + (0 + x | g)) + poisson(), data = dd)
  cc <- coef(fit)
  expect_named(cc, "g")
  expect_equal(colnames(cc$g), c("(Intercept)", "x"))
  expect_equal(cc$g$x,
               unname(fixef(fit)$mu["x"] + ranef(fit)[[2]][, "x"]),
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
  ysim <- simulate(fit, nsim = 1, re.form = NA)[[1L]]

  rf <- refit(fit, ysim)
  dd2 <- dd
  dd2$y <- ysim
  fresh <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd2)

  expect_loglik_equal(rf, fresh, tol = 1e-6)
  expect_vector_equal(fixef(rf)$mu, fixef(fresh)$mu, tol = 1e-5)
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

  sims <- simulate(fit, nsim = 5, re.form = NA, seed = 1)
  boots <- vapply(sims, function(ys) fixef(refit(fit, ys))$mu["x"],
                  numeric(1))
  expect_length(boots, 5)
  expect_true(all(is.finite(boots)))
  expect_lt(abs(mean(boots) - fixef(fit)$mu["x"]), 0.3)
})
