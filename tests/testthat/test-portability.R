# brms portability: the spellings a ported vignette uses. Every block
# names the finding it closes in dev/brms-vignette-port.md.

port_data <- function(seed = 1, n = 60) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(n / 10, 10), y = 0)
  # portability is about argument spellings reaching the same model, so
  # the response is drawn from the model itself
  # The draw takes its own seed, away from the fixture's: reusing
  # the fixture seed restarts the same random stream that made the
  # covariates, and the residuals come out equal to x.
  d$y <- frm_simulate(bf(y ~ x) + gaussian(), d,
                      newparams = list(Intercept = 1, x = 2, sigma = 1),
                      nsim = 1, seed = seed + 1000L)[[1]]
  d$k <- as.integer(cut(d$y, 3))
  d
}

test_that("FN-1: an unnamed family is gaussian, as in brms", {
  d <- port_data()
  f_def <- frm(y ~ x, data = d)
  f_exp <- frm(bf(y ~ x), family = gaussian(), data = d)
  expect_identical(as.numeric(logLik(f_def)), as.numeric(logLik(f_exp)))
  expect_identical(attr(logLik(f_def), "df"), attr(logLik(f_exp), "df"))
  expect_identical(f_def$spec$responses[[1]]$family$family, "gaussian")
  # silent: a message here would fire on every linear mixed model
  expect_silent(as_bform(y ~ x))
})

test_that("FN-1: an attached family still wins over the default", {
  d <- port_data()
  d$cnt <- rpois(nrow(d), 3)
  b <- as_bform(bf(cnt ~ x) + poisson())
  expect_identical(b$family$family, "poisson")
  # the family argument overrides a univariate attachment, and fills
  # the empty slots of a multivariate one
  b2 <- as_bform(bf(cnt ~ x) + poisson(), family = gaussian())
  expect_identical(b2$family$family, "gaussian")
  mv <- as_bform(mvbf(bf(cnt ~ x) + poisson(), bf(y ~ x)))
  expect_identical(vapply(mv$forms, function(f) f$family$family, ""),
                   c("poisson", "gaussian"))
  mv2 <- as_bform(mvbf(bf(cnt ~ x) + poisson(), bf(y ~ x)),
                  family = Gamma())
  expect_identical(vapply(mv2$forms, function(f) f$family$family, ""),
                   c("poisson", "Gamma"))
})

test_that("FN-6: a family given as a bare constructor is accepted", {
  d <- port_data()
  fit <- frm(bf(k ~ x), data = d, family = cumulative)
  expect_identical(fit$spec$responses[[1]]$family$family, "cumulative")
  expect_identical(as_frmtmb_family(gaussian)$family, "gaussian")
  expect_identical(as_frmtmb_family(poisson)$family, "poisson")
  expect_error(as_frmtmb_family(1:3), "Cannot interpret")
})

test_that("FN-5: hypothesis() takes brms's directional form", {
  d <- port_data()
  fit <- frm(bf(y ~ x), family = gaussian(), data = d)
  two <- hypothesis(fit, "x")
  gt <- hypothesis(fit, "x > 0")
  lt <- hypothesis(fit, "x < 0")
  expect_identical(gt$estimate, two$estimate)
  # one-sided p is exactly half the two-sided one, on the right side
  expect_equal(gt$p, two$p / 2, tolerance = 1e-10)
  expect_equal(lt$p, 1 - two$p / 2, tolerance = 1e-10)
  # one-sided interval: unbounded on the side the alternative points at
  expect_identical(gt$upr, Inf)
  expect_identical(lt$lwr, -Inf)
  expect_equal(gt$lwr, two$estimate - qnorm(0.95) * two$se,
               tolerance = 1e-10)
  expect_equal(lt$upr, two$estimate + qnorm(0.95) * two$se,
               tolerance = 1e-10)
  expect_identical(attr(gt, "direction"), "greater")
  # "lhs > rhs" is the difference of the two sides
  expect_equal(hypothesis(fit, "x > Intercept")$estimate,
               hypothesis(fit, "x - Intercept")$estimate,
               tolerance = 1e-10)
  expect_output(print(gt), "one-sided")
  expect_error(hypothesis(fit, "x > 0 > 1"), "at most one")
  expect_error(hypothesis(fit, "x > 0 = 1"), "not both")
})

test_that("FN-5: class/group name the natural-scale summaries", {
  d <- port_data()
  # the correlated-slope block is singular on this data, which is
  # irrelevant here: the test is about how the names are built
  fit <- suppressWarnings(
    frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d)
  )
  a <- hypothesis(fit, "Intercept - x > 0", class = "sd", group = "g")
  b <- hypothesis(fit, "sd_g__Intercept - sd_g__x > 0")
  expect_equal(a$estimate, b$estimate, tolerance = 1e-10)
  expect_equal(a$p, b$p, tolerance = 1e-10)
  # class "b" and no class are the plain coefficient names
  expect_equal(hypothesis(fit, "x", class = "b")$estimate,
               hypothesis(fit, "x")$estimate, tolerance = 1e-10)
  # a class without a group is a dpar prefix
  fd <- frm(bf(y ~ x, sigma ~ x), family = gaussian(), data = d)
  expect_equal(hypothesis(fd, "Intercept", class = "sigma")$estimate,
               hypothesis(fd, "sigma_Intercept")$estimate,
               tolerance = 1e-10)
  # a class/group that names nothing is refused, not silently dropped
  flin <- frm(bf(y ~ x), family = gaussian(), data = d)
  expect_error(hypothesis(flin, "x > 0", class = "sd", group = "g"),
               "not a parameter of this model")
  expect_error(hypothesis(fit, "x > 0", class = "sd", group = "typo"),
               "not a parameter of this model")
  # a name already written in full keeps its spelling
  expect_equal(
    hypothesis(fit, "sd_g__Intercept > 0", class = "sd", group = "g")$estimate,
    hypothesis(fit, "sd_g__Intercept")$estimate, tolerance = 1e-10)
})

test_that("FN-5: the bootstrap and profile bounds are one-sided too", {
  d <- port_data()
  fit <- frm(bf(y ~ x), family = gaussian(), data = d)
  bo <- hypothesis(fit, c("x > 0", "x"), method = "boot", nsim = 40,
                   seed = 3)
  expect_identical(bo$upr[1], Inf)
  expect_true(bo$lwr[1] > bo$lwr[2])
  # tail proportion with the (1 + k) / (1 + n) correction
  dr <- attr(bo, "draws")[, 1]
  expect_equal(bo$p[1], (1 + sum(dr <= 0)) / (1 + length(dr)),
               tolerance = 1e-10)
  pr <- hypothesis(fit, c("x > 0", "x"), method = "profile")
  expect_identical(pr$upr[1], Inf)
  expect_true(pr$lwr[1] > pr$lwr[2])
})

test_that("FN-9: update() takes the brms argument spellings", {
  d <- port_data()
  fit <- frm(bf(y ~ x), family = gaussian(), data = d)
  u1 <- update(fit, formula. = ~ . + z)
  expect_identical(deparse1(formula(u1)), "y ~ x + z")
  u2 <- update(fit, formula = y ~ x, newdata = d[1:40, ])
  expect_identical(nobs(u2), 40L)
  u3 <- update(fit, ~ . - x)
  expect_identical(deparse1(formula(u3)), "y ~ 1")
  # the stats::update() spelling with a dotted left-hand side
  u3b <- update(fit, . ~ . + z)
  expect_identical(deparse1(formula(u3b)), "y ~ x + z")
  expect_equal(as.numeric(logLik(u3b)),
               as.numeric(logLik(update(fit, ~ . + z))),
               tolerance = 1e-10)
  # a delta may change the response as well
  d$w <- d$y + 1
  u3c <- update(fit, w ~ . + z, data = d)
  expect_identical(deparse1(formula(u3c)), "w ~ x + z")
  # a two-sided formula with no dot replaces the stored one
  u3d <- update(fit, y ~ z)
  expect_identical(deparse1(formula(u3d)), "y ~ z")
  # the delta keeps the dpar formulas and the family
  fd <- frm(bf(y ~ x, sigma ~ x), family = gaussian(), data = d)
  u4 <- update(fd, ~ . + z)
  expect_identical(deparse1(formula(u4)), "y ~ x + z")
  expect_named(u4$bform$pforms, "sigma")
  u5 <- update(fd, . ~ . + z)
  expect_identical(deparse1(formula(u5)), "y ~ x + z")
  expect_named(u5$bform$pforms, "sigma")
  expect_error(update(fit, ~ . + z, data = d, newdata = d), "not both")
  # a delta is refused where it would be ambiguous
  fnl <- frm(bf(y ~ a * x + b, a ~ 1, b ~ 1, nl = TRUE), data = d,
             family = gaussian(), start = list(beta = c(1, 1)))
  expect_error(update(fnl, . ~ . + z), "nonlinear formula")
  fmv <- frm(bf(y ~ x) + bf(z ~ x), family = gaussian(), data = d)
  expect_error(update(fmv, . ~ . + x), "which response")
})

test_that("FN-10: lf() adds parameter formulas to a bf()", {
  d <- port_data()
  a <- frm(bf(y ~ x) + lf(sigma ~ x), family = gaussian(), data = d)
  b <- frm(bf(y ~ x, sigma ~ x), family = gaussian(), data = d)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-6)
  expect_s3_class(lf(sigma ~ x), "frmtmb_lf")
  expect_error(lf(~x), "two-sided")
  expect_error(lf(), "at least one")
  expect_error(lf(sigma ~ x, sigma ~ z), "Duplicated")
  expect_error(bf(y ~ x, sigma ~ z) + lf(sigma ~ x), "already sets")
  expect_error(mvbf(bf(y ~ x), bf(z ~ x)) + lf(sigma ~ x),
               "which response")
})

test_that("FN-8: one formula can name several parameters", {
  d <- port_data()
  a <- frm(bf(y ~ a * x + b, a + b ~ 1, nl = TRUE), data = d,
           family = gaussian(), start = list(beta = c(1, 1)))
  b <- frm(bf(y ~ a * x + b, a ~ 1, b ~ 1, nl = TRUE), data = d,
           family = gaussian(), start = list(beta = c(1, 1)))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-8)
  # the right-hand side is copied to each name, dpars included
  f <- bf(y ~ x, sigma + nu ~ x)
  expect_named(f$pforms, c("sigma", "nu"))
  expect_identical(deparse1(f$pforms$nu), "nu ~ x")
  expect_error(bf(y ~ x, sigma + sigma ~ x), "Duplicated")
  expect_error(bf(y ~ x, sigma + a.b ~ x), "Invalid parameter name")
})

test_that("parameter names: one vocabulary across the methods", {
  set.seed(4)
  ng <- 24
  n <- ng * 10
  d <- data.frame(x = rnorm(n), fosternest = factor(rep(1:ng, 10)))
  z1 <- rnorm(ng)
  z2 <- rnorm(ng)
  u1 <- 0.7 * z1
  u2 <- 0.5 * (0.6 * z1 + sqrt(0.64) * z2)
  d$tarsus <- 1 + 0.5 * d$x + u1[d$fosternest] + rnorm(n)
  d$back <- -1 + 0.3 * d$x + u2[d$fosternest] + rnorm(n)
  # the user's repro: an mvbf fit, where the two vocabularies differ
  fit <- frm(bf(mvbind(tarsus, back) ~ x + (1 | p | fosternest)),
             family = gaussian(), data = d)
  expect_true("tarsus_(Intercept)" %in% rownames(confint(fit)))
  expect_true("tarsus_Intercept" %in% variables(fit))
  # each entry point now takes the other's spelling
  expect_identical(profile(fit, "tarsus_Intercept"),
                   profile(fit, "tarsus_(Intercept)"))
  expect_identical(confint(fit, "tarsus_Intercept"),
                   confint(fit, "tarsus_(Intercept)"))
  expect_equal(hypothesis(fit, "`tarsus_(Intercept)`")$estimate,
               hypothesis(fit, "tarsus_Intercept")$estimate,
               tolerance = 1e-12)
  # variables() still lists one spelling only
  expect_false(any(grepl("[()]", variables(fit))))
  # an sd_ alias resolves to the theta it names, and says so
  al <- frmtmb:::par_alias_index(fit)
  sd_nm <- "sd_fosternest__tarsus.muIntercept"
  expect_true(sd_nm %in% names(al))
  th <- outer_par_names(fit)[al[[sd_nm]]]
  expect_message(a <- confint(fit, sd_nm), "internal")
  expect_identical(a, suppressMessages(confint(fit, th)))
  expect_identical(rownames(a), th)
  # a 2x2 us block carries one internal correlation parameter, so its
  # cor name is one-to-one
  expect_true(paste0("cor_fosternest__tarsus.muIntercept__",
                     "back.muIntercept") %in% names(al))
  # unknown names advertise both vocabularies
  expect_error(confint(fit, "nope"), "Parentheses may be dropped")
  expect_error(profile(fit, "nope"), "natural-scale")
  # bounds take the parenthesis-free spelling too
  expect_identical(
    frmtmb:::resolve_bounds(fit, c(tarsus_Intercept = 0), NULL)$lower,
    frmtmb:::resolve_bounds(fit, c(`tarsus_(Intercept)` = 0),
                            NULL)$lower)
})

test_that("an autocorrelation parameter is addressable by its name", {
  set.seed(11)
  ns <- 20
  nt <- 6
  d <- expand.grid(t = 1:nt, subj = factor(1:ns))
  d$x <- rnorm(nrow(d))
  e <- unlist(lapply(seq_len(ns), function(i) {
    as.numeric(stats::arima.sim(list(ar = 0.6), nt, sd = 0.8))
  }))
  d$y <- 1 + 0.5 * d$x + e
  fit <- frm(bf(y ~ x + ar(t, subj, cov = TRUE)), family = gaussian(),
             data = d)
  al <- frmtmb:::par_alias_index(fit)
  # one thetaac entry, so the natural name addresses it
  expect_true("ar1" %in% names(al))
  th <- outer_par_names(fit)[al[["ar1"]]]
  expect_identical(th, "thetaac_1")
  expect_message(a <- confint(fit, "ar1"), "internal")
  expect_identical(a, suppressMessages(confint(fit, th)))
})

test_that("a name that is not one internal parameter is refused", {
  set.seed(6)
  n <- 200
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:20, 10)))
  b <- matrix(rnorm(60, 0, 0.6), 20, 3)
  d$y <- 1 + b[d$g, 1] + (0.5 + b[d$g, 2]) * d$x +
    (0.2 + b[d$g, 3]) * d$z + rnorm(n)
  fit <- suppressWarnings(
    frm(bf(y ~ x + z + (1 + x + z | g)), family = gaussian(), data = d)
  )
  al <- frmtmb:::par_alias_index(fit)
  # a 3x3 us block: the sds are one-to-one
  expect_true("sd_g__Intercept" %in% names(al))
  # its correlations are Cholesky mixtures, so they are not
  expect_false("cor_g__Intercept__x" %in% names(al))
  expect_error(profile(fit, "cor_g__Intercept__x"),
               "does not stand for a single internal one")
  expect_error(confint(fit, "cor_g__Intercept__x", method = "profile"),
               "method = 'profile'")
  # hypothesis() is the route named, and it handles the combination
  expect_true(is.finite(hypothesis(fit, "cor_g__Intercept__x")$estimate))
})

test_that("FN-11: frm_multiple refuses what it cannot pool", {
  d <- port_data()
  imps <- lapply(1:3, function(i) {
    z <- d
    z$y <- z$y + rnorm(nrow(d), 0, 0.1)
    z
  })
  fm <- frm_multiple(y ~ x, data = imps, family = gaussian())
  expect_error(plot(fm), "no pooled display")
  expect_error(conditional_effects(fm, "x"), "no pooled version")
  skip_if_not_installed("posterior")
  expect_error(posterior::as_draws_array(fm), "needs draws")
  expect_error(posterior::nchains(fm), "needs draws")
  # pooled hypothesis tests keep working, directionally too
  h <- hypothesis(fm, "x > 0")
  expect_identical(h$upr, Inf)
  expect_equal(h$p, hypothesis(fm, "x")$p / 2, tolerance = 1e-10)
})
