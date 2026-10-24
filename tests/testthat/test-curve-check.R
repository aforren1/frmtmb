# post$curve_check: the hook a family uses to refuse a curve off a fit
# whose curve it knows to be an artifact. frmtmb.survival's
# royston_parmar() refuses a fit with a non-monotone cumulative hazard
# through it (its test-rp-floored.R); here it is checked on a family of
# this file's own, with the case where the hook is absent beside it.

curve_check_fit <- function(hook) {
  set.seed(3)
  d <- data.frame(x = sort(stats::runif(120)))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(120, 0, 0.3)
  fam <- as_frmtmb_family(gaussian())
  seen <- new.env()
  if (!is.null(hook)) {
    fam[["post"]][["curve_check"]] <- function(fit, resp) {
      seen$resp <- c(seen$resp, resp)
      hook(fit, resp)
    }
  }
  list(fit = frm(y ~ s(x, k = 6), family = fam, data = d), seen = seen,
       g = data.frame(x = seq(0.1, 0.9, length.out = 9)))
}

test_that("each curve function asks the family before it assembles", {
  o <- curve_check_fit(function(fit, resp) {
    frm_stop("tcheck: this fit has no curve to report", call. = FALSE)
  })
  expect_error(frm_curve(o$fit, newdata = o$g, simultaneous = FALSE),
               "tcheck: this fit has no curve to report")
  expect_error(frm_curve_deriv(o$fit, var = "x", newdata = o$g,
                               simultaneous = FALSE),
               "tcheck: this fit has no curve to report")
  expect_error(frm_curve_feature(o$fit, var = "x", newdata = o$g),
               "tcheck: this fit has no curve to report")
  # the hook is handed the response it is asked about, once per call
  expect_identical(o$seen$resp, rep("y", 3L))
})

test_that("a hook that returns, and a family without one, change nothing", {
  ref <- curve_check_fit(NULL)
  ok <- curve_check_fit(function(fit, resp) invisible(NULL))
  a <- frm_curve(ref$fit, newdata = ref$g, simultaneous = FALSE)
  b <- frm_curve(ok$fit, newdata = ok$g, simultaneous = FALSE)
  expect_identical(a$.estimate, b$.estimate)
  expect_identical(a$.se, b$.se)
  expect_identical(ok$seen$resp, "y")
  expect_null(ref$seen$resp)
})
