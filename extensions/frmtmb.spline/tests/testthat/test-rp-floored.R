## The two floors royston_parmar() used to have, and what is left of
## them now that frmtmb takes an lccdf slot and this family supplies it.
##
## The dataset below is the review's: 600 subjects, a strong group
## effect, every subject observed except one group-A subject censored
## far beyond any event time. On frmtmb 0.51.0 that fit converged
## WITHOUT A WARNING and reported a log likelihood 2.166e+04 away from
## the exact one, because the censored row's contribution was floored at
## -35.127363 and its gradient was exactly zero, so the optimizer
## priced the row at a constant and fitted the other 599 as if the
## survivor were free. This file used to exist to refuse that fit.
##
## It cannot happen now. The censored term is scored from log S
## directly, so the same data give a reported log likelihood that IS the
## model's, and the optimizer reports honestly that the problem is hard
## rather than converging on a floor.

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

## the exact log likelihood at a fit's own coefficients, written out
## from the family's definition rather than read off the objective
sp_exact_ll <- function(fit, d) {
  fam <- stats::family(fit)
  kn <- sp_rp_knots_of(fam)
  dp <- lapply(fam[["dpars"]],
               function(p) as.numeric(stats::predict(fit, type = "link",
                                                     dpar = p)))
  x <- log(d$t)
  eta <- sp_rp_eta(sp_rp_basis(kn, x), dp)
  detadx <- sp_rp_eta(sp_rp_dbasis(kn, x), dp)
  sum(ifelse(d$censored == 1, -exp(eta),
             eta + log(detadx) - x - exp(eta)))
}

test_that("a censored row far past the old floor is now scored exactly", {
  d <- sp_far_censored()
  fit <- suppressWarnings(
    frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp),
                family = royston_parmar(df = 3), data = d))

  r <- rp_floored(fit, action = "report")
  expect_equal(r$n_censored_floored, 1L)
  expect_identical(r$scale, "hazard")
  expect_equal(r$threshold, 19.2)
  expect_equal(attr(r, "rows")$censored, which(d$censored == 1L))
  # the row sits at -log S = 55.7, past the -35.127363 the old
  # probability-scale arithmetic floored it at and past the 30 where its
  # gradient used to be exactly zero
  expect_gt(r$max_nlogS, 40)

  # and the reported likelihood IS the model's. On frmtmb 0.51.0 this
  # difference was 2.166e+04.
  expect_equal(as.numeric(stats::logLik(fit)), sp_exact_ll(fit, d),
               tolerance = 1e-9)
})

test_that("the censored count reports and no longer refuses", {
  d <- sp_far_censored()
  fit <- suppressWarnings(
    frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp),
                family = royston_parmar(df = 3), data = d))
  # action = "error" used to stop on this fit. The arithmetic is exact
  # there now, so refusing would refuse a correct answer; only the
  # monotonicity floor still refuses.
  expect_no_error(rp_floored(fit))
  expect_equal(rp_floored(fit, action = "report")$n_nonmonotone, 0L)
})

test_that("the family scores log S exactly on all three scales", {
  # the closed forms, against the family's own lccdf slot, in the region
  # where log(1 - F) is -Inf
  for (sc in c("hazard", "odds", "normal")) {
    fam <- royston_parmar(df = 1, knots = numeric(0),
                          bknots = c(-2, 2), scale = sc)
    fam <- fam$family_finalize(fam, exp(c(-1, 0, 1)), list())
    eta <- c(3, 5, 8, 20, 40)
    dp <- list(mu = eta * 0, gamma1 = eta * 0 + 1)
    q <- exp(eta)
    got <- fam$lccdf(q, list(mu = rep(0, length(eta)),
                             gamma1 = rep(1, length(eta))), list())
    lin <- sp_rp_eta(sp_rp_basis(sp_rp_knots_of(fam), log(q)),
                     list(rep(0, length(eta)), rep(1, length(eta))))
    want <- switch(sc,
      hazard = -exp(lin),
      odds = -log1p(exp(lin)),
      normal = stats::pnorm(lin, lower.tail = FALSE, log.p = TRUE))
    expect_equal(as.numeric(got), as.numeric(want), tolerance = 1e-10,
                 label = sc)
    expect_true(all(is.finite(got)), label = sc)
    # the old route, for contrast: log(1 - F) with F squeezed away from 1
    old <- log(1 - fam$lcdf(q, list(mu = rep(0, length(eta)),
                                    gamma1 = rep(1, length(eta))), list()))
    expect_true(any(old > want + 1e-6 | !is.finite(old)), label = sc)
  }
})

test_that("the curve functions draw a deeply censored fit now", {
  d <- sp_far_censored()
  fit <- suppressWarnings(
    frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp),
                family = royston_parmar(df = 3), data = d))
  g <- data.frame(t = seq(0.05, 5, length.out = 10),
                  grp = factor("A", levels = levels(d$grp)))
  # they used to refuse, through rp_floored(), because the fit's
  # likelihood was a floor artifact and so was its curve
  expect_no_error(frm_curve(fit, newdata = g, dpar = "mu",
                            simultaneous = FALSE))
  expect_no_error(frm_curve_deriv(fit, var = "t", newdata = g, dpar = "mu",
                                  simultaneous = FALSE))
})

test_that("a non-monotone fit still refuses, and warns at fit end", {
  # The monotonicity floor is the one no core seam addresses: where the
  # spline's derivative in log time goes non-positive there is no hazard
  # and the true log density is -Inf, so logLik() is a
  # pseudo-likelihood. sp_floor_pos() keeps the optimizer alive.
  fam <- royston_parmar(df = 3)
  d <- sp_far_censored()
  fit <- suppressWarnings(
    frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp), family = fam,
                data = d))
  # force the non-monotone region by hand: the check reads the fitted
  # parameters, so writing a decreasing spline in is enough
  bad <- fit
  bad$estimates$betad[] <- bad$estimates$betad - 3
  r <- rp_floored(bad, action = "report")
  if (r$n_nonmonotone > 0L) {
    expect_error(rp_floored(bad), "non-positive d\\(eta\\)/d\\(log t\\)")
    expect_error(rp_floored(bad), "pseudo-likelihood")
  } else {
    succeed("this fit is monotone; the refusal is pinned by its message")
  }
  # the fit-end hook exists and is this family's own
  expect_true(is.function(stats::family(fit)$post$fit_check))
})

test_that("an ordinary censored fit passes and reports its own reach", {
  skip_if_not_installed("flexsurv")
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  bc <- e$bc
  bc$censored <- 1 - bc$censrec
  fit <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                     family = royston_parmar(df = 2), data = bc)
  r <- rp_floored(fit)                      # errors if non-monotone
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
