# Date, POSIXct and difftime predictors: accepted, but fitted as the
# number underneath, in an origin the user did not choose. srr G2.5.

date_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- data.frame(day = as.Date("2020-01-01") + seq_len(n) - 1L)
  d$days <- as.numeric(d$day - min(d$day))
  d$y <- stats::rnorm(n, 1 + 0.02 * d$days, 0.5)
  d
}

test_that("a Date predictor is reported, with its unit and origin", {
  d <- date_data()
  expect_message(
    frm(bf(y ~ day) + gaussian(), data = d, dry_run = "frame"),
    "day \\(Date, days since 1970-01-01\\)"
  )
  expect_message(
    frm(bf(y ~ day) + gaussian(), data = d, dry_run = "frame"),
    "Center the column"
  )
})

test_that("POSIXct and difftime are reported with their own units", {
  d <- date_data(30)
  d$ts <- as.POSIXct("2020-01-01", tz = "UTC") + d$days * 86400
  expect_message(
    frm(bf(y ~ ts) + gaussian(), data = d, dry_run = "frame"),
    "ts \\(POSIXct, seconds since 1970-01-01\\)"
  )
  d$dt <- difftime(d$day, min(d$day), units = "days")
  expect_message(
    frm(bf(y ~ dt) + gaussian(), data = d, dry_run = "frame"),
    "dt \\(difftime, days\\)"
  )
})

test_that("an ordinary numeric predictor says nothing", {
  d <- date_data(30)
  expect_no_message(
    frm(bf(y ~ days) + gaussian(), data = d, dry_run = "frame")
  )
})

test_that("the coercion keeps the slope and moves the intercept", {
  d <- date_data()
  raw <- suppressWarnings(suppressMessages(
    frm(bf(y ~ day) + gaussian(), data = d)))
  ctr <- frm(bf(y ~ days) + gaussian(), data = d)

  # same model, so the same slope
  expect_equal(unname(fixef(raw)$mu[["day"]]),
               unname(fixef(ctr)$mu[["days"]]), tolerance = 1e-5)
  # but the intercept is the value at 1970-01-01, thousands of days out,
  # which is why the message points at centering
  expect_lt(fixef(raw)$mu[["(Intercept)"]], -100)
  expect_equal(
    fixef(raw)$mu[["(Intercept)"]] +
      fixef(raw)$mu[["day"]] * as.numeric(min(d$day)),
    unname(fixef(ctr)$mu[["(Intercept)"]]), tolerance = 1e-4
  )
})

test_that("a Date grouping variable says nothing", {
  # it is used for its distinct levels, where a Date behaves exactly
  # like a factor, so advice about coefficients and intercepts would be
  # wrong
  d <- date_data(40)
  d$g <- as.Date("2020-01-01") + rep(0:3, each = 10)
  expect_no_message(
    fit <- frm(bf(y ~ days + (1 | g)) + gaussian(), data = d))
  expect_identical(unname(ngrps(fit)[["g"]]), 4L)
})

test_that("a Date response says nothing", {
  # an explicit as.numeric() converts the response, and for a gaussian
  # a location shift of it is absorbed by the intercept
  set.seed(2)
  d <- data.frame(x = stats::rnorm(40))
  d$yd <- as.Date("2020-01-01") + round(10 * d$x) + stats::rpois(40, 20)
  expect_no_message(fit <- frm(bf(yd ~ x) + gaussian(), data = d))
  expect_equal(unname(fixef(fit)$mu[["x"]]), 10, tolerance = 1)
})

test_that("a difftime response says nothing either", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(40))
  d$yd <- as.difftime(10 * d$x + stats::rnorm(40), units = "days")
  expect_no_message(frm(bf(yd ~ x) + gaussian(), data = d,
                        dry_run = "frame"))
})
