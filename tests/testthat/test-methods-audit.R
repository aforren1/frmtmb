# Method-surface audit fixes: stats/lme4/glmmTMB argument conventions
# and the expected-response predict() semantics.

meth_env <- local({
  set.seed(31)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.5)[dd$g], 1)
  list(dd = dd,
       fit = frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd))
})

test_that("vcov(full = TRUE) is labeled like confint rows", {
  fit <- meth_env$fit
  V <- vcov(fit, full = TRUE)
  expect_identical(rownames(V), rownames(confint(fit)))
  # the coefficient block is vcov() itself
  nm <- rownames(vcov(fit))
  expect_equal(V[nm, nm], vcov(fit), tolerance = 1e-12)
})

test_that("confint accepts the Wald spelling and a boot method", {
  fit <- meth_env$fit
  expect_identical(confint(fit, method = "Wald"),
                   confint(fit, method = "wald"))
  ci <- confint(fit, parm = "x", method = "boot", nsim = 30, seed = 1)
  expect_identical(colnames(ci), c("lwr", "upr", "est"))
  expect_identical(rownames(ci), "x")
  expect_true(ci[, "lwr"] < ci[, "est"] && ci[, "est"] < ci[, "upr"])
  # percentile interval roughly agrees with Wald on a clean gaussian fit
  cw <- confint(fit, parm = "x", method = "wald")
  expect_lt(abs(ci[, "lwr"] - cw[, "lwr"]), 0.15)
  expect_lt(abs(ci[, "upr"] - cw[, "upr"]), 0.15)
})

test_that("simulate follows the stats seed contract", {
  fit <- meth_env$fit
  s1 <- simulate(fit, nsim = 2, seed = 42)
  s2 <- simulate(fit, nsim = 2, seed = 42)
  expect_identical(names(s1), c("sim_1", "sim_2"))
  expect_equal(s1$sim_1, s2$sim_1)
  expect_identical(attr(attr(s1, "seed"), "kind")[[1]], RNGkind()[1])
  # a seeded call must not disturb the global RNG stream
  set.seed(99)
  before <- .Random.seed
  simulate(fit, nsim = 1, seed = 7)
  expect_identical(before, .Random.seed)
  # unseeded calls store the pre-call RNG state, as stats::simulate does
  s3 <- simulate(fit, nsim = 1)
  expect_identical(attr(s3, "seed"), before)
})

test_that("na.action returns the fit's na.action", {
  fit <- meth_env$fit
  expect_null(stats::na.action(fit))
  dd <- meth_env$dd
  dd$x[3] <- NA
  fit_na <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
                na.action = stats::na.exclude)
  expect_s3_class(stats::na.action(fit_na), "exclude")
  expect_identical(as.integer(stats::na.action(fit_na)), 3L)
})

test_that("summary reports the grouping-factor sizes", {
  expect_output(print(summary(meth_env$fit)), "Groups: g, 10")
})

test_that("drop1 matches lme4", {
  skip_if_not_installed("lme4")
  dd <- meth_env$dd
  fit <- frm(bf(y ~ x + I(x^2) + (1 | g)) + gaussian(), data = dd)
  d1 <- drop1(fit, test = "Chisq")
  expect_s3_class(d1, "anova")
  expect_identical(rownames(d1), c("<none>", "x", "I(x^2)"))

  ref <- lme4::lmer(y ~ x + I(x^2) + (1 | g), data = dd, REML = FALSE)
  dref <- drop1(ref, test = "Chisq")
  expect_vector_equal(d1$AIC, dref$AIC, tol = 1e-4)
  expect_vector_equal(d1$LRT[-1], dref$LRT[-1], tol = 1e-4)
  expect_vector_equal(d1$Df[-1], dref$npar[-1], tol = 1e-12)

  # marginality: interactions protect their main effects
  fit2 <- frm(bf(y ~ x * g) + gaussian(), data = dd)
  expect_identical(rownames(drop1(fit2)), c("<none>", "x:g"))
  fit_reml <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
                  REML = TRUE)
  expect_error(drop1(fit_reml), "REML")
})

test_that("cooks.distance on the fit and dfbeta/dfbetas match lme4", {
  fit <- meth_env$fit
  infl <- influence(fit, groups = "g")
  expect_equal(cooks.distance(fit, groups = "g"),
               cooks.distance(infl))
  db <- dfbeta(infl)
  dbs <- dfbetas(infl)
  expect_identical(dim(db), dim(infl$fixed))
  # dfbetas is dfbeta scaled by the coefficient SEs
  expect_equal(dbs, sweep(db, 2, sqrt(diag(vcov(fit))), `/`),
               tolerance = 1e-12)
  # sign convention: full-data estimate minus leave-one-out estimate
  # (influence columns carry the vcov labels since the v0.21 alignment)
  expect_equal(db[1, "x"],
               unname(fixef(fit)$mu["x"] - infl$fixed[1, "x"]))

  skip_if_not_installed("lme4")
  ref <- lme4::lmer(y ~ x + (1 | g), data = meth_env$dd, REML = FALSE)
  iref <- stats::influence(ref, groups = "g")
  # lme4's influence dfbeta is deleted-minus-full; ours follows the
  # stats::dfbeta.lm sign (full-minus-deleted), hence the negation
  dref <- -stats::dfbeta(iref)
  expect_vector_equal(db[rownames(dref), "x"], dref[, "x"],
                      tol = 1e-3)
})

test_that("predict type = 'response' is the expected response (zi)", {
  skip_if_not_installed("glmmTMB")
  set.seed(71)
  n <- 400
  x <- rnorm(n)
  z <- rnorm(n)
  mu <- exp(0.5 + 0.4 * x)
  zi <- plogis(-0.5 + 0.8 * z)
  dd <- data.frame(y = rbinom(n, 1, 1 - zi) * rpois(n, mu), x = x, z = z)
  fit <- suppressWarnings(
    frm(bf(y ~ x, zi ~ z) + zero_inflated_poisson(), data = dd)
  )
  ref <- glmmTMB::glmmTMB(y ~ x, ziformula = ~z, family = poisson,
                          data = dd)

  # the fitted() invariant every reference package satisfies
  expect_equal(predict(fit, type = "response"), fitted(fit),
               tolerance = 1e-12)
  # glmmTMB type values agree numerically
  expect_vector_equal(predict(fit, type = "response"),
                      predict(ref, type = "response"), tol = 1e-4)
  expect_vector_equal(predict(fit, type = "conditional"),
                      predict(ref, type = "conditional"), tol = 1e-4)
  expect_vector_equal(predict(fit, type = "zprob"),
                      predict(ref, type = "zprob"), tol = 1e-4)
  expect_vector_equal(predict(fit, type = "zlink"),
                      predict(ref, type = "zlink"), tol = 1e-3)
  # an explicit dpar keeps the per-dpar meaning
  expect_equal(predict(fit, dpar = "mu", type = "response"),
               predict(fit, type = "conditional"), tolerance = 1e-12)
  # newdata goes through the same mean
  nd <- dd[1:5, ]
  expect_equal(predict(fit, newdata = nd, type = "response"),
               unname(fitted(fit)[1:5]), tolerance = 1e-8)
  # the expected response now carries joint delta-method SEs, and they
  # agree with glmmTMB's own response-scale delta method
  ps <- predict(fit, type = "response", se.fit = TRUE)
  rs <- predict(ref, type = "response", se.fit = TRUE)
  expect_vector_equal(ps$fit, rs$fit, tol = 1e-3)
  expect_vector_equal(ps$se.fit, rs$se.fit, tol = 1e-3)
  expect_error(predict(fit, type = "disp"), "dispersion")
})

test_that("predict type aliases and spellings on gaussian fits", {
  fit <- meth_env$fit
  # identity-mean family: response semantics are unchanged
  expect_equal(predict(fit, type = "response"), fitted(fit),
               tolerance = 1e-12)
  expect_equal(predict(fit, type = "conditional"),
               predict(fit, type = "response"), tolerance = 1e-12)
  expect_equal(unique(round(predict(fit, type = "disp"), 10)),
               round(sigma(fit), 10))
  # the lme4/glmmTMB dot spelling of allow_new_levels is accepted
  nd <- data.frame(x = 0, g = factor("99"))
  expect_identical(predict(fit, newdata = nd, allow.new.levels = TRUE),
                   predict(fit, newdata = nd, allow_new_levels = TRUE))
  expect_warning(predict(fit, bogus_arg = 1), "bogus_arg")
})

test_that("predict type = 'response' scales binomial means by trials", {
  set.seed(32)
  n <- 200
  x <- rnorm(n)
  sz <- sample(5:10, n, replace = TRUE)
  p <- plogis(-0.3 + 0.6 * x)
  # the trials variable lives only in the data, so a newdata call that
  # omits it cannot fall back to the formula environment
  dd <- data.frame(y = rbinom(n, sz, p), size = sz, x = x)
  fit <- frm(bf(y | trials(size) ~ x) + binomial(), data = dd)
  # brms epred convention (and the fitted() invariant): counts scale
  expect_equal(predict(fit, type = "response"), fitted(fit),
               tolerance = 1e-12)
  nd <- data.frame(x = c(0, 0), size = c(1, 10))
  pr <- predict(fit, newdata = nd, type = "response")
  expect_equal(pr[2], 10 * pr[1], tolerance = 1e-10)
  # missing trials in newdata is an error, not a silent 1
  expect_error(predict(fit, newdata = data.frame(x = 0),
                       type = "response"), "trials")
})
