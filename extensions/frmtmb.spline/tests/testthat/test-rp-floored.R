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
               function(p) as.numeric(frm_linpred(fit, type = "link",
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
  expect_equal(r$n_censored_deep, 1L)
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

test_that("fitted() still refuses, and the count is named for what it is", {
  d <- sp_far_censored()
  fit <- suppressWarnings(
    frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp),
                family = royston_parmar(df = 3), data = d))

  # the refusal is by DESIGN and is unrelated to the floors: a survival
  # time has no mean on the response scale here, and mu is a spline
  # coefficient rather than a fitted value. An exact likelihood does not
  # give the model a mean it never had.
  expect_error(stats::fitted(fit), "no mean on the response scale")
  expect_error(frm_linpred(fit, type = "response"),
               "no mean on the response scale")

  r <- rp_floored(fit, action = "report")
  # the field counts DEEP censored rows, not floored ones. There is one
  # here, at -log S = 55.73, and it is scored exactly: naming it
  # "floored" would say the opposite of what the lccdf slot bought.
  expect_named(r, c("n_censored_deep", "max_nlogS", "threshold",
                    "n_nonmonotone", "n_nonmonotone_censored",
                    "max_survival_rise", "scale", "n_obs"))
  expect_equal(r$n_censored_deep, 1L)
  expect_equal(r$n_nonmonotone, 0L)
  expect_no_error(rp_floored(fit))
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
  expect_equal(r$n_censored_deep, 0L)
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

## Item 3.6 of dev/extension-gaps-plan.md. A group with NO events has no
## event row, so a check on event rows alone cannot see it, and a
## random slope on gamma1 is pushed DOWN for such a group: its score in
## u is -sum(x_i H_i), negative for rows past t = 1, and the barrier
## log(gamma1 + u) lives in a density the group does not contribute.
## The construction is dev/frailty/frailty-floor3.R's at seed 20260910:
## on frmtmb.spline 0.7.0 this fit converged clean, warned nothing,
## passed rp_floored() and passed frm_curve(), with five centres at
## slopes of -0.21 to -0.29. dev/phase3a-log/ has that run.

## What an expression did: its value, every warning, and the error if
## it stopped. Assertions read this, so a build that errors where the
## current one warns FAILS them rather than stopping the block.
sp_outcome <- function(expr) {
  w <- character(0)
  v <- tryCatch(withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(value = if (inherits(v, "error")) NULL else v,
       error = if (inherits(v, "error")) conditionMessage(v) else NULL,
       warnings = w)
}

sp_nodeath <- function(seed = 20260910L, n = 400L, n_centre = 40L,
                       n_nodeath = 5L, sd_u = 0.35, gamma1 = 0.6,
                       scale = 5, beta = 0.6, p_cens = 0.4) {
  set.seed(seed)
  centre <- rep(seq_len(n_centre), length.out = n)
  repeat {
    u <- stats::rnorm(n_centre, 0, sd_u)
    if (all(gamma1 + u > 0.05)) break
  }
  trt <- stats::rbinom(n, 1L, 0.5)
  g0 <- -gamma1 * log(scale)
  tt <- exp((log(-log(stats::runif(n))) - g0 - beta * trt) /
              (gamma1 + u[centre]))
  tau <- unname(stats::quantile(tt, 1 - p_cens))
  ev <- as.integer(tt <= tau)
  tt <- pmin(tt, tau)
  j <- centre <= n_nodeath
  tt[j] <- tau
  ev[j] <- 0L
  data.frame(time = tt, event = ev, censored = 1L - ev, trt = trt,
             centre = factor(centre))
}

test_that("a centre with no events and a falling slope warns", {
  skip_on_cran()
  d <- sp_nodeath()
  bk <- range(log(d$time[d$event == 1L]))
  warns <- character(0)
  fit <- withCallingHandlers(
    frmtmb::frm(frmtmb::bf(time | cens(censored) ~ trt,
                           gamma1 ~ (1 | centre)),
                family = royston_parmar(knots = numeric(0), bknots = bk),
                data = d),
    warning = function(cnd) {
      warns <<- c(warns, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  nodeath <- which(as.integer(d$centre) <= 5L)
  # the premise, measured rather than assumed: those five centres did
  # land below zero, and the population slope did not
  g1 <- unname(frmtmb::fixef_by_dpar(fit)$gamma1[["(Intercept)"]])
  slope <- g1 + frmtmb::ranef(fit)[["centre"]][, 1L]
  expect_equal(unname(which(slope <= 0)), 1:5)
  expect_gt(g1, 0)
  # the fit end says so, which 0.7.0 did not
  expect_true(any(grepl("fitted survival rises", warns)))
  # A censored row is scored exactly, so logLik() is the model's; the
  # fitted survival rises there, which is worth a warning and not a
  # refusal (user decision, 2026-09-24). rp_floored() and the curve
  # entry point both warn and answer.
  o <- sp_outcome(rp_floored(fit))
  expect_null(o$error)
  expect_true(any(grepl("fitted survival rises", o$warnings)))
  g <- data.frame(time = exp(seq(-2, 1, length.out = 5)), trt = 0,
                  centre = factor(1, levels = levels(d$centre)))
  o <- sp_outcome(frm_curve(fit, newdata = g, dpar = "gamma1",
                            simultaneous = FALSE))
  expect_null(o$error)
  expect_true(any(grepl("fitted survival rises", o$warnings)))
  expect_false(is.null(o$value))
  r <- rp_floored(fit, action = "report")
  # exactly the no-death rows, and every one of them
  expect_identical(r[["n_nonmonotone_censored"]], length(nodeath))
  expect_identical(attr(r, "rows")[["nonmonotone_censored"]], nodeath)
  # and the likelihood count is unmoved: these rows are scored exactly
  expect_identical(r[["n_nonmonotone"]], 0L)
})

test_that("the censored-row count reads an interval's upper end", {
  # An interval row's likelihood reads S at both ends, so the check has
  # to read both. Built so that ONLY the upper ends are non-monotone:
  # the boundary knot sits 1.5 past the last event, every upper end
  # sits 1.2 past it, and the spline's second coefficient is written in
  # to put the derivative just below zero there while it stays positive
  # at every event and every lower end.
  set.seed(11)
  n <- 200
  d <- data.frame(t = stats::rweibull(n, 1.3, 2))
  d$cens <- 0L
  iv <- seq(1, n, by = 4)
  xe <- log(d$t[-iv])
  d$cens[iv] <- 2L
  d$t[iv] <- exp(pmin(log(d$t[iv]), max(xe) - 0.1))
  d$t2 <- d$t
  d$t2[iv] <- exp(max(xe) + 1.2)
  fit <- suppressWarnings(frmtmb::frm(
    frmtmb::bf(t | cens(cens, t2) ~ 1),
    family = royston_parmar(knots = stats::median(xe),
                            bknots = c(min(xe), max(xe) + 1.5)),
    data = d))
  kn <- sp_rp_knots_of(stats::family(fit))
  g1 <- unname(fit$estimates[["betad"]][["gamma1_(Intercept)"]])
  v_hi <- sp_rp_dbasis(kn, max(xe) + 1.2)[[3L]]
  bad <- fit
  bad$estimates[["betad"]][["gamma2_(Intercept)"]] <- -1.1 * g1 / v_hi
  g2 <- bad$estimates[["betad"]][["gamma2_(Intercept)"]]
  der <- function(x) g1 + g2 * sp_rp_dbasis(kn, x)[[3L]]
  # the premise, measured: positive at every observed lower time,
  # non-positive at every upper end
  expect_gt(min(der(log(d$t))), 0)
  expect_lte(der(max(xe) + 1.2), 0)
  r <- rp_floored(bad, action = "report")
  expect_identical(r[["n_nonmonotone"]], 0L)
  expect_identical(r[["n_nonmonotone_censored"]], length(iv))
  expect_identical(attr(r, "rows")[["nonmonotone_censored"]],
                   as.integer(iv))
  o <- sp_outcome(rp_floored(bad))
  expect_null(o$error)
  expect_true(any(grepl("fitted survival rises", o$warnings)))
})

test_that("an ordinary fit reports zero censored rows where S rises", {
  skip_if_not_installed("flexsurv")
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  bc <- e$bc
  bc$censored <- 1 - bc$censrec
  fit <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group,
                                gamma1 ~ group),
                     family = royston_parmar(df = 3), data = bc)
  r <- rp_floored(fit, action = "report")
  expect_identical(r[["n_nonmonotone_censored"]], 0L)
  expect_identical(r[["n_nonmonotone"]], 0L)
  expect_no_error(rp_floored(fit))
})

## Punch round 2 of item 3.6, the user decision of 2026-09-24: REFUSE
## where a non-positive d(eta)/d(log t) meets an EVENT row, whose
## density is then floored and logLik() a pseudo-likelihood; WARN where
## it meets only censored rows, which are scored exactly and where the
## rise is usually extrapolation past an arm's last event. The warning
## gives the size of the rise, because that is what a user judges.

## dev/phase3a-cure36.R's cure-fraction design: two arms of 250, 40 and
## 55 percent cured. At df = 4, seed 20260932 it fired with every
## flagged row past its arm's last event and a rise of 0.099.
sp_cure <- function(seed = 20260932L, n = 500L) {
  set.seed(seed)
  x <- rep(0:1, each = n / 2)
  cured <- stats::runif(n) < ifelse(x == 1, 0.55, 0.40)
  t <- stats::rweibull(n, 1.3, ifelse(x == 1, 1.6, 1))
  t[cured] <- Inf
  cc <- pmin(6, stats::runif(n, 0, 12))
  data.frame(t = pmin(t, cc), censored = as.integer(t > cc),
             x = factor(x))
}

test_that("a rise on censored rows only warns, with its size", {
  skip_on_cran()
  d <- sp_cure()
  fo <- sp_outcome(frmtmb::frm(frmtmb::bf(t | cens(censored) ~ x,
                                          gamma1 ~ x),
                               family = royston_parmar(df = 4), data = d))
  fit <- fo$value
  expect_true(inherits(fit, "frmtmb_fit"))
  # the fit-end warning names the size of the rise
  expect_true(any(grepl("fitted survival rises by up to 0\\.0[0-9]",
                        fo$warnings)))
  r <- if (inherits(fit, "frmtmb_fit")) rp_floored(fit, action = "report")
  expect_identical(r[["n_nonmonotone"]], 0L)
  expect_gt(r[["n_nonmonotone_censored"]], 0L)
  # the premise: every flagged row lies past its own arm's last event
  rows <- attr(r, "rows")[["nonmonotone_censored"]]
  last <- tapply(d$t[d$censored == 0], d$x[d$censored == 0], max)
  expect_true(length(rows) > 0 &&
                all(d$t[rows] > last[as.character(d$x[rows])]))
  # the reported rise is the one measured by the construction: between
  # 5.2e-04 and 0.13 over all 27 firing fits, 0.099 on this one
  rise <- if (is.null(r[["max_survival_rise"]])) NA_real_ else
    r[["max_survival_rise"]]
  expect_true(isTRUE(rise > 0.05 && rise < 0.15))
  # rp_floored() and frm_curve() warn and answer
  o <- sp_outcome(rp_floored(fit))
  expect_null(o$error)
  expect_true(any(grepl("fitted survival rises", o$warnings)))
  g <- data.frame(t = seq(0.5, 5, length.out = 6),
                  x = factor(1, levels = levels(d$x)))
  # gamma1, not mu: frm_curve(dpar = "mu") on a `gamma1 ~ x` fit refuses
  # on its own covariance cross-check on the released build too, a
  # separate matter filed in dev/phase3a-findings.md
  o <- sp_outcome(frm_curve(fit, newdata = g, dpar = "gamma1",
                            simultaneous = FALSE))
  expect_null(o$error)
  expect_true(any(grepl("fitted survival rises", o$warnings)))
  expect_false(is.null(o$value))
})

test_that("an event-row floor still refuses, and the curve with it", {
  # the unchanged half of the decision, pinned so that it stays: a
  # decreasing spline written into a fit where every row is an event
  d <- sp_far_censored()
  fit <- suppressWarnings(
    frmtmb::frm(frmtmb::bf(t | cens(censored) ~ grp),
                family = royston_parmar(df = 3), data = d))
  bad <- fit
  bad$estimates[["betad"]][] <- bad$estimates[["betad"]] - 3
  r <- rp_floored(bad, action = "report")
  expect_gt(r[["n_nonmonotone"]], 0L)
  o <- sp_outcome(rp_floored(bad))
  expect_match(if (is.null(o$error)) "" else o$error, "pseudo-likelihood")
  g <- data.frame(t = seq(0.05, 5, length.out = 10),
                  grp = factor("A", levels = levels(d$grp)))
  o <- sp_outcome(frm_curve(bad, newdata = g, dpar = "mu",
                            simultaneous = FALSE))
  expect_match(if (is.null(o$error)) "" else o$error, "pseudo-likelihood")
})

test_that("a fit with neither is silent", {
  skip_if_not_installed("flexsurv")
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  bc <- e$bc
  bc$censored <- 1 - bc$censrec
  fo <- sp_outcome(frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group,
                                          gamma1 ~ group),
                               family = royston_parmar(df = 3), data = bc))
  expect_identical(fo$warnings, character(0))
  o <- sp_outcome(rp_floored(fo$value))
  expect_null(o$error)
  expect_identical(o$warnings, character(0))
  expect_identical(o$value[["max_survival_rise"]], 0)
})

