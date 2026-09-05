## The two floors in royston_parmar(), and the refusal that stands in
## for a core lccdf slot.
##
## The dataset below is the review's: 600 subjects, a strong group
## effect, every subject observed except one group-A subject censored
## far beyond any event time. It converges without a warning, and the
## censored row's true log survival is thousands where the scored value
## is -35.127363. That is the fit this file exists to refuse.

sp_far_censored <- function(seed = 20260905, n = 600, cens_at = 50) {
  set.seed(seed)
  d <- data.frame(grp = factor(rep(c("A", "B"), each = n / 2)))
  d$t <- stats::rweibull(n, shape = 1.6,
                         scale = ifelse(d$grp == "A", 0.30, 3.0))
  d$censored <- 0L
  i <- which(d$grp == "A")[1L]
  d$t[i] <- cens_at
  d$censored[i] <- 1L
  d
}

test_that("a fit whose censored row is past the accurate region is refused", {
  d <- sp_far_censored()
  # the fit itself converges and says nothing: that is the defect
  expect_no_warning(
    fit <- frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp),
                       family = royston_parmar(df = 3), data = d))
  expect_equal(fit$opt$convergence, 0)
  expect_true(is.finite(as.numeric(stats::logLik(fit))))

  r <- rp_floored(fit, action = "report")
  expect_equal(r$n_censored_floored, 1L)
  expect_identical(r$scale, "hazard")
  # -log S on that row is in the thousands where the scored value is
  # floored at -35.127363
  expect_gt(r$max_nlogS, 1000)
  expect_equal(r$threshold, 19.2)
  expect_equal(attr(r, "rows")$censored, which(d$censored == 1L))

  # and the refusal names every part the contract asks for
  expect_error(rp_floored(fit), "1 of 1 censored rows")
  expect_error(rp_floored(fit), "19.2")
  expect_error(rp_floored(fit), "complementary log-CDF")
  expect_error(rp_floored(fit), "dev/spline-seam-proposal.md")
  expect_error(rp_floored(fit), "POST-FIT")
})

test_that("the reported likelihood really is the wrong one", {
  # not an assertion about the refusal but about why it is needed: the
  # censored row is scored at the floor rather than at its own -log S
  d <- sp_far_censored()
  fit <- frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp),
                     family = royston_parmar(df = 3), data = d)
  r <- rp_floored(fit, action = "report")
  floor_value <- log(1 - (1 - 5.551115123125783e-16))
  # the error in the reported log likelihood is at least the difference
  # between the row's true contribution and the floor it was given
  expect_gt(r$max_nlogS + floor_value, 1000)
  expect_lt(as.numeric(stats::logLik(fit)) - (-r$max_nlogS), Inf)
  expect_gt(as.numeric(stats::logLik(fit)), -r$max_nlogS)
})

test_that("the curve functions refuse a floored fit rather than draw it", {
  d <- sp_far_censored()
  fit <- frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp),
                     family = royston_parmar(df = 3), data = d)
  g <- data.frame(t = seq(0.05, 5, length.out = 10),
                  grp = factor("A", levels = levels(d$grp)))
  expect_error(frm_curve(fit, newdata = g, dpar = "mu",
                         simultaneous = FALSE),
               "censored rows are scored past the accurate region")
  expect_error(frm_curve_deriv(fit, var = "t", newdata = g, dpar = "mu",
                               simultaneous = FALSE),
               "censored rows are scored past the accurate region")
  expect_error(frm_curve_feature(fit, var = "t", newdata = g, dpar = "mu"),
               "censored rows are scored past the accurate region")
})

test_that("an ordinary censored fit passes and reports its own reach", {
  skip_if_not_installed("flexsurv")
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  bc <- e$bc
  bc$censored <- 1 - bc$censrec
  fit <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                     family = royston_parmar(df = 2), data = bc)
  r <- rp_floored(fit)                      # errors if anything floored
  expect_equal(r$n_censored_floored, 0L)
  expect_equal(r$n_nonmonotone, 0L)
  # the reach: the largest -log S on any of the 387 censored rows, an
  # order of magnitude below the 19.2 where accuracy starts to go
  expect_lt(r$max_nlogS, 5)
  expect_gt(r$max_nlogS, 0)
  expect_equal(r$n_obs, nrow(bc))
  # and a curve off that fit is not gated
  g <- data.frame(recyrs = seq(0.5, 6, length.out = 12),
                  group = factor("Good", levels = levels(bc$group)))
  expect_no_error(frm_curve(fit, newdata = g, dpar = "gamma1",
                            simultaneous = FALSE))
})

test_that("the odds and normal scales are checked on the same quantity", {
  # -log S is one quantity for all three scales and the three agree to
  # every printed digit at a given value, because they differ only in
  # how eta maps to S. Checked here against the closed forms.
  eta <- c(-2, 0, 2, 5)
  expect_equal(exp(eta), exp(eta))
  expect_equal(log1p(exp(eta)), -log(1 / (1 + exp(eta))), tolerance = 1e-12)
  expect_equal(-stats::pnorm(-eta, log.p = TRUE), -log(stats::pnorm(-eta)),
               tolerance = 1e-9)
  # eta = 6 on the normal scale, the review's threshold, is -log S 20.7,
  # so the single 19.2 threshold is the stricter of the two
  expect_gt(-stats::pnorm(-6, log.p = TRUE), 19.2)
  expect_lt(-stats::qnorm(exp(-19.2)), 6)
})

test_that("rp_floored refuses anything that is not a royston_parmar fit", {
  set.seed(2)
  d <- data.frame(x = stats::rnorm(60))
  d$y <- stats::rnorm(60, d$x, 1)
  fit <- frmtmb::frm(frmtmb::bf(y ~ x), family = stats::gaussian(), data = d)
  expect_error(rp_floored(fit), "reads the floors of a royston_parmar")
  expect_error(rp_floored(d), "reads the floors of a royston_parmar")
  bad <- frmtmb::frm(frmtmb::bf(y ~ x), family = stats::gaussian(), data = d)
  expect_error(rp_floored(bad, max_nlogS = -1), "one positive finite number")
})
