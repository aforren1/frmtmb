# brms 2.23.0's model surface on items a fit used to refuse or read
# differently: s() of an expression, attr(x, "min") of a mo() variable,
# a later lf() on a parameter bf() sets, the REML df, unused grouping
# levels, the place of a non-centered intercept, default_prior() of
# rr(), frm_multiple() displays and insight::get_predicted() on a
# multinomial fit. brms's values are recorded in dev/parity2-findings.md.

test_that("s() of an expression fits, as s() of the column it computes", {
  set.seed(4)
  d <- data.frame(x = stats::runif(80, 0.1, 2))
  d$y <- sin(2 * d$x) + stats::rnorm(80, sd = 0.3)
  d$ex <- exp(d$x)
  f1 <- frm(bf(y ~ s(exp(x))), family = gaussian(), data = d)
  f2 <- frm(bf(y ~ s(ex)), family = gaussian(), data = d)
  ll1 <- as.numeric(logLik(f1))
  ll2 <- as.numeric(logLik(f2))
  expect_lt(abs(ll1 - ll2), 1e-6 * abs(ll2))
  # brms's names: bs_sexpx_1, sds_sexpx_1
  expect_true(all(c("bs_sexpx_1", "sds_sexpx_1") %in% variables(f1)))
  nd <- data.frame(x = c(0.5, 1.5))
  p1 <- fitted(f1, newdata = nd)[, "Estimate"]
  p2 <- fitted(f2, newdata = data.frame(ex = exp(nd$x)))[, "Estimate"]
  expect_lt(max(abs(p1 - p2)), 1e-5 * max(abs(p2)))
  ce <- suppressMessages(conditional_effects(f1, resolution = 5))
  expect_identical(names(ce), "x")
  # a non-syntactic column has no name brms gives; brms 2.23.0 stops
  # parsing s(`my x`) ("unexpected symbol"), and frmtmb refuses by name
  d[["my x"]] <- d$x
  expect_error(frm(bf(y ~ s(`my x`)), family = gaussian(), data = d),
               "not a syntactic name")
})

test_that("attr(x, \"min\") of a mo() variable sets its lowest category", {
  set.seed(9)
  d <- data.frame(k = rep(1:4, 10))
  d$y <- 0.5 * d$k + stats::rnorm(40)
  attr(d$k, "min") <- 0
  # the declared category 0 holds no row, so its simplex cell and the
  # intercept trade off: the SE check names them, as it should
  flat <- c("Standard errors are not available",
            "move along a direction the fit does not determine")
  fit <- allow_warnings(frm(bf(y ~ mo(k)), family = gaussian(), data = d),
                        flat, require = flat[1L])
  m <- fit$frame$linpreds[[1L]]$mo[[1L]]
  # brms 2.23.0: Jmo = 4 and Xmo = 1 2 3 4 for the values 1..4
  expect_identical(m$min, 0)
  expect_identical(as.integer(m$D), 4L)
  expect_identical(as.integer(m$codes[1:4]), 1:4)
  # new data are coded from the declared minimum, so a value of the fit
  # predicts what the fit has at it, and the declared category answers
  nd <- data.frame(k = 0:4)
  p <- allow_warnings(fitted(fit, newdata = nd)[, "Estimate"], flat)
  p_in <- allow_warnings(fitted(fit)[1:4, "Estimate"], flat)
  expect_equal(p[2:5], unname(p_in), tolerance = 64 * .Machine$double.eps)
  expect_length(p, 5L)
  d2 <- d
  attr(d2$k, "min") <- 2
  expect_error(frm(bf(y ~ mo(k)), family = gaussian(), data = d2),
               "no larger than the smallest value")
})

test_that("a later lf() replaces the formula bf() set, as in brms", {
  f <- NULL
  expect_message(f <- bf(y ~ x, sigma ~ z) + lf(sigma ~ w),
                 "Replacing initial definitions of parameters 'sigma'")
  expect_identical(deparse1(f$pforms$sigma), "sigma ~ w")
  f2 <- NULL
  expect_message(f2 <- bf(y ~ x, sigma = 1) + lf(sigma ~ w), "Replacing")
  expect_null(f2$pfix$sigma)
  expect_identical(deparse1(f2$pforms$sigma), "sigma ~ w")
  set.seed(2)
  d <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60),
                  w = stats::rnorm(60))
  d$y <- stats::rnorm(60, d$x, exp(0.3 * d$w))
  fa <- frm(f, family = gaussian(), data = d)
  fb <- frm(bf(y ~ x, sigma ~ w), family = gaussian(), data = d)
  expect_identical(as.numeric(logLik(fa)), as.numeric(logLik(fb)))
})

test_that("the df of a REML fit counts the fixed effects, as lme4's", {
  set.seed(5)
  d <- data.frame(x = stats::rnorm(90), g = factor(rep(1:15, 6)))
  d$y <- 1 + 0.4 * d$x + stats::rnorm(15, sd = 0.7)[d$g] +
    stats::rnorm(90)
  ml <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  re <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d,
            REML = TRUE)
  # lme4: 2 coefficients, 1 sd, sigma, under ML and REML alike
  expect_identical(attr(logLik(re), "df"), attr(logLik(ml), "df"))
  expect_identical(as.integer(attr(logLik(re), "df")), 4L)
  expect_equal(AIC(re), -2 * as.numeric(logLik(re)) + 2 * 4,
               tolerance = 64 * .Machine$double.eps)
})

test_that("drop_unused_levels = FALSE keeps an unused level in the term", {
  set.seed(3)
  du <- data.frame(g = factor(sample(1:6, 60, TRUE), levels = 1:7),
                   x = stats::rnorm(60))
  du$y <- stats::rnorm(60, 1 + 0.4 * du$x) + c(-0.5, 0.5)[(du$g %in% 1:3) +
                                                           1]
  fk <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = du,
            drop_unused_levels = FALSE)
  fd <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = du)
  # brms 2.23.0 make_standata(drop_unused_levels = FALSE): N_1 = 7
  re <- ranef(fk)$g
  expect_identical(rownames(re), as.character(1:7))
  # no row reads level 7, so its mode is its prior's, 0
  expect_lte(abs(re["7", 1]), 1e-8 * max(abs(re[, 1])))
  expect_identical(rownames(ranef(fd)$g), as.character(1:6))
  # a level no row reads leaves the likelihood where it was
  expect_lt(abs(as.numeric(logLik(fk)) - as.numeric(logLik(fd))),
            1e-6 * abs(as.numeric(logLik(fd))))
})

test_that("a non-centered intercept keeps its place among the b_", {
  set.seed(6)
  d <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60))
  d$y <- stats::rnorm(60, 1 + d$x, exp(0.2 * d$z))
  v <- function(f) {
    grep("^b_", variables(frm(f, family = gaussian(), data = d)),
         value = TRUE)
  }
  # brms 2.23.0's variables(), in order
  expect_identical(v(bf(y ~ 0 + Intercept + x, sigma ~ z)),
                   c("b_sigma_Intercept", "b_Intercept", "b_x",
                     "b_sigma_z"))
  expect_identical(v(bf(y ~ x, sigma ~ 0 + Intercept + z)),
                   c("b_Intercept", "b_x", "b_sigma_Intercept",
                     "b_sigma_z"))
  expect_identical(v(bf(y ~ x, sigma ~ z)),
                   c("b_Intercept", "b_sigma_Intercept", "b_x",
                     "b_sigma_z"))
})

test_that("default_prior() lists a model with an rr() term", {
  set.seed(2)
  d <- data.frame(g = factor(rep(1:10, each = 6)),
                  tf = factor(rep(c("a", "b", "c"), 20)),
                  x = stats::rnorm(60))
  d$y <- stats::rnorm(60)
  dp <- default_prior(bf(y ~ x + rr(0 + tf | g, d = 1)),
                      family = gaussian(), data = d)
  expect_true("b" %in% dp$class)
  expect_false(any(dp$class == "sd" & dp$group == "g"))
})

test_that("frm_multiple() displays take int_conditions", {
  set.seed(7)
  d <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60))
  d$y <- 1 + d$x * d$z + stats::rnorm(60)
  fm <- frm_multiple(bf(y ~ x * z), family = gaussian(),
                     data = list(d, d[c(2:60, 1), ]))
  ce <- conditional_effects(fm, effects = "x:z",
                            int_conditions = list(z = c(-1, 1)))
  expect_identical(sort(unique(as.numeric(as.character(ce[[1L]]$z)))),
                   c(-1, 1))
})

test_that("insight::get_predicted() refuses a multinomial fit by name", {
  skip_if_not_installed("insight")
  set.seed(5)
  d <- data.frame(x = stats::rnorm(60))
  pm <- cbind(1, exp(0.5 * d$x), exp(-0.3 * d$x))
  d$Y <- t(apply(pm / rowSums(pm), 1, function(p) stats::rmultinom(1, 10, p)))
  colnames(d$Y) <- c("a", "b", "c")
  d$size <- 10L
  fit <- frm(bf(Y | trials(size) ~ x), family = multinomial(), data = d)
  expect_error(insight::get_predicted(fit), "one per row and category")
})
