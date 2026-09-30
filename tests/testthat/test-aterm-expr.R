# An addition term given an expression, `weights(wt * 2)`, is evaluated
# on the data rows, as brms's get_ad_values() evaluates it. Before, the
# expression went into the model-frame formula as it stood, where
# `wt * 2` is an interaction and `s / 2` a nesting, and the frame died
# with "invalid model formula in ExtractVars" (dev/formrobust-findings.md
# section 1, dev/formrobust-repro1.R and -repro1b.R).

ae_data <- local({
  set.seed(11)
  n <- 60
  d <- data.frame(x = stats::rnorm(n), wt = stats::runif(n, 0.5, 2),
                  time = stats::runif(n, 1, 3), n = stats::rpois(n, 5) + 2L,
                  s = stats::runif(n, 0.2, 0.6))
  d$yc <- stats::rpois(n, exp(0.2 + 0.3 * d$x) * 2 * d$time)
  d$yb <- stats::rbinom(n, d$n + 1, stats::plogis(0.3 * d$x))
  d$y <- stats::rnorm(n, 0.5 * d$x, 1)
  d$y2 <- stats::rnorm(n, -0.5 * d$x, 1)
  d$lb <- -3
  # the same values, precomputed as columns
  d$w2 <- d$wt * 2
  d$t2 <- d$time * 2
  d$n1 <- d$n + 1
  d$s2 <- d$s / 2
  d$lbm <- d$lb - 1
  d
})

# The expression and its precomputed column are the same data, so the
# two fits are the same fit: an identity, compared bitwise
ae_same <- function(fa, fb, resp = NULL) {
  expect_identical(logLik(fa), logLik(fb))
  nd <- ae_data[1:6, ]
  expect_identical(unname(fitted(fa, newdata = nd, resp = resp)),
                   unname(fitted(fb, newdata = nd, resp = resp)))
  set.seed(3)
  pa <- predict(fa, newdata = nd, resp = resp)
  set.seed(3)
  pb <- predict(fb, newdata = nd, resp = resp)
  expect_identical(unname(pa), unname(pb))
  if (is.null(resp)) {
    expect_identical(unname(as.matrix(simulate(fa, nsim = 2, seed = 1,
                                               newdata = nd))),
                     unname(as.matrix(simulate(fb, nsim = 2, seed = 1,
                                               newdata = nd))))
  }
}

test_that("weights(), rate(), se() and trunc() take an expression", {
  ae_same(frm(bf(y | weights(wt * 2) ~ x), data = ae_data),
          frm(bf(y | weights(w2) ~ x), data = ae_data))
  ae_same(frm(bf(yc | rate(time * 2) ~ x), data = ae_data,
              family = poisson()),
          frm(bf(yc | rate(t2) ~ x), data = ae_data, family = poisson()))
  ae_same(frm(bf(y | se(s / 2) ~ x), data = ae_data),
          frm(bf(y | se(s2) ~ x), data = ae_data))
  ae_same(frm(bf(y | se(s / 2, sigma = TRUE) ~ x), data = ae_data),
          frm(bf(y | se(s2, sigma = TRUE) ~ x), data = ae_data))
  ae_same(frm(bf(y | trunc(lb = lb - 1) ~ x), data = ae_data),
          frm(bf(y | trunc(lb = lbm) ~ x), data = ae_data))
  ae_same(frm(bf(yb | trials(n + 1) ~ x), data = ae_data,
              family = binomial()),
          frm(bf(yb | trials(n1) ~ x), data = ae_data, family = binomial()))
})

test_that("an expression term works in a multivariate model", {
  fa <- frm(bf(y | weights(wt * 2) ~ x) + bf(y2 | se(s / 2) ~ x),
            data = ae_data)
  fb <- frm(bf(y | weights(w2) ~ x) + bf(y2 | se(s2) ~ x), data = ae_data)
  ae_same(fa, fb, resp = "y")
  ae_same(fa, fb, resp = "y2")
  d <- ae_data
  d$sx <- d$x > 0
  fa <- frm(bf(yc | rate(time * 2) + subset(x > 0) ~ x, family = poisson()) +
              bf(y ~ x), data = d)
  fb <- frm(bf(yc | rate(t2) + subset(sx) ~ x, family = poisson()) +
              bf(y ~ x), data = d)
  expect_identical(logLik(fa), logLik(fb))
})

test_that("a term reads its variables from the data alone, as brms", {
  # an outside value would be read again by predictions and refits, and
  # after saveRDS()/readRDS() in another session it is missing or
  # changed; on the response side, as trials(k), the change is silent
  k <- 2
  kt <- 12L
  d <- ae_data
  d$yt <- pmin(d$yb, 12L)
  msg <- "reads `k`, which is not a column of `data`"
  expect_error(frm(bf(y | weights(wt * k) ~ x), data = ae_data), msg,
               fixed = TRUE)
  expect_error(frm(bf(yc | rate(time * k) ~ x), data = ae_data,
                   family = poisson()), msg, fixed = TRUE)
  expect_error(frm(bf(yt | trials(kt) ~ x), data = d, family = binomial()),
               "reads `kt`, which is not a column of `data`", fixed = TRUE)
  expect_error(frm(bf(yt | trials(kt + 0L) ~ x), data = d,
                   family = binomial()),
               "Put `kt` in the data", fixed = TRUE)
  # a full-length vector outside the data too, which 0.66.0 read
  w_out <- ae_data$wt
  expect_error(frm(bf(y | weights(w_out) ~ x), data = ae_data),
               "reads `w_out`, which is not a column of `data`",
               fixed = TRUE)
  # the guard-absent cases: the same constant as a data column
  d$kt <- 12L
  f_lit <- frm(bf(yt | trials(12) ~ x), data = d, family = binomial())
  expect_identical(logLik(frm(bf(yt | trials(kt) ~ x), data = d,
                              family = binomial())), logLik(f_lit))
  expect_identical(logLik(frm(bf(yt | trials(kt + 0L) ~ x), data = d,
                              family = binomial())), logLik(f_lit))
  d$k <- 2
  ae_same(frm(bf(y | weights(wt * k) ~ x), data = d),
          frm(bf(y | weights(w2) ~ x), data = d))
  # function calls are not variables
  d$lw <- log(d$wt) + 1
  expect_identical(logLik(frm(bf(y | weights(log(wt) + 1) ~ x),
                              data = d)),
                   logLik(frm(bf(y | weights(lw) ~ x), data = d)))
  # and a constant in a predictor keeps the R convention
  kx <- 3
  d2 <- ae_data
  d2$x3 <- d2$x * 3
  expect_identical(logLik(frm(bf(y ~ I(x * kx)), data = ae_data)),
                   logLik(frm(bf(y ~ x3), data = d2)))
})

test_that("a variable R finds only as a function is refused by name", {
  # without a column t, t is base::t() and the product failed with
  # "non-numeric argument to binary operator"
  expect_error(frm(bf(y | weights(t * 2) ~ x), data = ae_data),
               "R finds only the function t()", fixed = TRUE)
  # the guard-absent case: a column named t is read
  d <- ae_data
  d$t <- d$wt
  ae_same(frm(bf(y | weights(t * 2) ~ x), data = d),
          frm(bf(y | weights(w2) ~ x), data = d))
})

test_that("an expression that gives NA is refused, as brms refuses it", {
  expect_error(frm(bf(y | trunc(lb = ifelse(x > 1, NA, -5)) ~ x),
                   data = ae_data),
               "trunc(lb = ifelse(x > 1, NA, -5)) of response 'y' is NA",
               fixed = TRUE)
  expect_error(frm(bf(yb | trials(ifelse(x > 1, NA, n + 1)) ~ x),
                   data = ae_data, family = binomial()),
               "is NA on", fixed = TRUE)
  expect_error(frm(bf(y | weights(ifelse(x > 1, NA, wt)) ~ x),
                   data = ae_data), "is NA on", fixed = TRUE)
  # the guard-absent case: an NA in the variable itself drops the row,
  # through na.action, as it does for any model variable
  d <- ae_data
  d$wt[1] <- NA
  f <- suppressMessages(frm(bf(y | weights(wt * 2) ~ x), data = d))
  expect_identical(nobs(f), nrow(d) - 1L)
})

test_that("the interval bound of cens() takes one value per row", {
  d <- ae_data
  d$cc <- ifelse(d$x > 1, "interval", "none")
  expect_error(frm(bf(y | cens(cc, max(y) + 10) ~ x), data = d),
               "the interval upper bound has 1 value(s) where the data have",
               fixed = TRUE)
  d$y2 <- max(d$y) + 10
  expect_identical(
    logLik(frm(bf(y | cens(cc, y2) ~ x), data = d)),
    logLik(frm(bf(y | cens(cc, y2 + 0) ~ x), data = d)))
})

test_that("a single value is recycled, as brms recycles it", {
  # min(y) is one value; the frame used to take it as a column and die
  # with "variable lengths differ"
  d <- ae_data
  d$lbmin <- min(d$y) - 1
  fa <- frm(bf(y | trunc(lb = min(y) - 1) ~ x), data = ae_data)
  fb <- frm(bf(y | trunc(lb = lbmin) ~ x), data = d)
  expect_identical(logLik(fa), logLik(fb))
  expect_identical(fitted(fa), fitted(fb))
  # on newdata the summary is taken over newdata's rows, as brms takes it
  expect_true(all(is.finite(fitted(fa, newdata = ae_data[1:6, ]))))
})

test_that("conditional_effects() draws an expression term", {
  fit <- frm(bf(y | weights(wt * 2) ~ x), data = ae_data)
  ce <- conditional_effects(fit, effects = "x")
  expect_true(all(is.finite(ce[[1]]$estimate__)))
  fr <- frm(bf(yc | rate(time * 2) ~ x), data = ae_data, family = poisson())
  ref <- frm(bf(yc | rate(t2) ~ x), data = ae_data, family = poisson())
  expect_identical(conditional_effects(fr, effects = "x")[[1]]$estimate__,
                   conditional_effects(ref, effects = "x")[[1]]$estimate__)
})

test_that("weights(scale = TRUE) is brms's scaling", {
  d <- ae_data
  d$wsc <- d$wt / sum(d$wt) * nrow(d)
  expect_identical(
    logLik(frm(bf(y | weights(wt, scale = TRUE) ~ x), data = ae_data)),
    logLik(frm(bf(y | weights(wsc) ~ x), data = d)))
  # the guard-absent case: without scale the weights are the raw ones
  expect_false(identical(
    logLik(frm(bf(y | weights(wt) ~ x), data = ae_data)),
    logLik(frm(bf(y | weights(wsc) ~ x), data = d))))
  expect_error(frm(bf(y | weights(wt, scale = "yes") ~ x), data = ae_data),
               "must be TRUE or FALSE", fixed = TRUE)
})

test_that("a term whose value has the wrong length is refused by name", {
  expect_error(frm(bf(y | weights(wt[1:3] * 2) ~ x), data = ae_data),
               "has 3 values where the data have 60 rows", fixed = TRUE)
})

test_that("the values are brms's own", {
  skip_unless_brms()
  sd <- brms_standata(y | weights(wt * 2) + se(s / 2, sigma = TRUE) ~ x,
                      data = ae_data)
  fr <- frm(bf(y | weights(wt * 2) + se(s / 2, sigma = TRUE) ~ x),
            data = ae_data, dry_run = "frame")
  av <- fr$aterm_values[[1L]]
  expect_identical(as.numeric(sd$weights), av[["weights"]])
  expect_identical(as.numeric(sd$se), av[["se"]])
  sd <- brms_standata(yc | rate(time * 2) ~ x, data = ae_data,
                      family = poisson())
  fr <- frm(bf(yc | rate(time * 2) ~ x), data = ae_data, family = poisson(),
            dry_run = "frame")
  expect_identical(as.numeric(sd$denom), fr$aterm_values[[1L]][["rate"]])
  sd <- brms_standata(y | trunc(lb = min(y) - 1) ~ x, data = ae_data)
  fr <- frm(bf(y | trunc(lb = min(y) - 1) ~ x), data = ae_data,
            dry_run = "frame")
  expect_identical(as.numeric(sd$lb), fr$aterm_values[[1L]][["trunc_lb"]])
  sd <- brms_standata(y | weights(wt, scale = TRUE) ~ x, data = ae_data)
  fr <- frm(bf(y | weights(wt, scale = TRUE) ~ x), data = ae_data,
            dry_run = "frame")
  expect_identical(as.numeric(sd$weights),
                   fr$aterm_values[[1L]][["weights"]])
})
