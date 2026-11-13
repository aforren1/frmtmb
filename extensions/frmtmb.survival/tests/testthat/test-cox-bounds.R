# Every time the cox() likelihood reads the cumulative baseline hazard at
# must lie inside the basis. Past a boundary knot the I-spline design is
# clamped to the knot (cox_cb_entry()), so a bound outside the knots
# reads the total cumulative hazard there, as though no hazard accrued
# beyond it. Through 0.81.0 the default knots and the check of given
# knots both read range(y) alone, and an interval-censoring upper bound
# past the last event was clamped without a word.

bounds_data <- function(seed = 11, n = 200) {
  set.seed(seed)
  x <- stats::rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  d <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
  d$hi <- d$time
  # the two latest rows become intervals that end past every response
  top <- order(d$time, decreasing = TRUE)[1:2]
  d$cens[top] <- 2
  d$hi[top] <- max(d$time) * c(1.5, 2)
  # and some ordinary intervals inside the range
  mid <- setdiff(seq(5, n, by = 10), top)
  d$cens[mid] <- 2
  d$hi[mid] <- d$time[mid] * 1.2
  d
}

# brms's default boundary from the response alone
brms_boundary <- function(y) {
  r <- range(y)
  c(max(r[1L] - diff(r) / 50, 0), r[2L] + diff(r) / 50)
}

cox_frame <- function(fo, d, fam = cox()) {
  frm(bf(fo), family = fam, data = d, dry_run = "frame")
}

test_that("the default knots cover an interval bound past the data", {
  d <- bounds_data()
  fr <- cox_frame(time | cens(cens, hi) ~ x, d)
  av <- fr$aterm_values$time
  bnd <- av$bhaz_spec$boundary
  i2 <- d$cens == 2
  expect_gte(bnd[2L], max(d$hi[i2]))
  # the cumulative design at the upper bounds: a bound past brms's knot
  # but below the last bound is not at the top of the basis, so its
  # I-spline row is not all ones as a clamped row is
  old_hi <- brms_boundary(d$time)[2L]
  e <- Filter(function(e) identical(e$q, as.numeric(av$cens_y2)),
              av$cox_cb)[[1L]]
  past <- which(i2 & d$hi > old_hi & d$hi < max(d$hi[i2]))
  expect_gte(length(past), 1L)
  eps <- .Machine$double.eps
  below_top <- apply(e$Z[past, , drop = FALSE], 1L,
                     function(r) any(r < 1 - 64 * eps))
  expect_true(all(below_top))
  # the lower knot did not move: no bound lies below the data
  expect_identical(bnd[1L], brms_boundary(d$time)[1L])
})

test_that("without such a bound the default is brms's, unchanged", {
  d <- bounds_data()
  d$cens[d$cens == 2] <- 0
  fr <- cox_frame(time | cens(cens) ~ x, d)
  expect_identical(fr$aterm_values$time$bhaz_spec$boundary,
                   brms_boundary(d$time))
  # an interval that ends inside the default knots leaves them alone
  d2 <- bounds_data()
  d2$cens[d2$hi > max(d2$time)] <- 0
  fr2 <- cox_frame(time | cens(cens, hi) ~ x, d2)
  expect_identical(fr2$aterm_values$time$bhaz_spec$boundary,
                   brms_boundary(d2$time))
})

test_that("given knots that miss an interval bound are refused", {
  d <- bounds_data()
  # they cover every response, so the old check passed them
  bk <- c(0, max(d$time) * 1.1)
  expect_lt(max(d$time), bk[2L])
  expect_error(cox_frame(time | cens(cens, hi) ~ x, d,
                         cox(Boundary.knots = bk)),
               "must cover every time the likelihood reads", fixed = TRUE)
  expect_error(cox_frame(time | cens(cens, hi) + bhaz(Boundary.knots = bk)
                         ~ x, d),
               "must cover every time the likelihood reads", fixed = TRUE)
  # the absent case: knots out to the last bound, and exactly at it
  for (top in c(max(d$hi) * 1.01, max(d$hi))) {
    fr <- cox_frame(time | cens(cens, hi) ~ x, d,
                    cox(Boundary.knots = c(0, top)))
    expect_identical(fr$aterm_values$time$bhaz_spec$boundary, c(0, top))
  }
  # a response outside the knots keeps its own message
  expect_error(cox_frame(time | cens(cens, hi) ~ x, d,
                         cox(Boundary.knots = c(0, max(d$time) / 2))),
               "must cover every observed time", fixed = TRUE)
})

test_that("a truncation upper bound is read like an interval bound", {
  d <- bounds_data()
  d$cens[d$cens == 2] <- 0
  d$ub <- max(d$time) * 1.5
  old <- brms_boundary(d$time)
  expect_gt(d$ub[1L], old[2L])
  fr <- cox_frame(time | cens(cens) + trunc(ub = ub) ~ x, d)
  expect_identical(fr$aterm_values$time$bhaz_spec$boundary,
                   c(old[1L], d$ub[1L]))
  expect_error(cox_frame(time | cens(cens) + trunc(ub = ub) ~ x, d,
                         cox(Boundary.knots = old)),
               "must cover every time the likelihood reads", fixed = TRUE)
})

test_that("delayed entry keeps brms's knots, and the fit is exact on them", {
  skip_if_not_installed("splines2")
  d <- bounds_data()
  d$cens[d$cens == 2] <- 0
  # times away from zero, so that brms's lower knot is above zero, and
  # entry times below it: each row was at risk from its entry time on
  d$time <- d$time + 2
  d$entry <- d$time * 0.5
  old <- brms_boundary(d$time)
  expect_lt(min(d$entry), old[1L])
  fo <- bf(time | cens(cens) + trunc(lb = entry) ~ x)
  fit <- frm(fo, family = cox(), data = d)
  sp <- fit$frame$aterm_values$time$bhaz_spec
  # the lower bounds do not move the knots: below the lower knot the
  # model has no hazard, so H0(entry) = 0 there is the model's own value.
  # Widening the knot to the entry times fitted 26.5 units worse here.
  expect_identical(sp$boundary, old)
  # brms's knots given by hand are the same fit, not a refusal
  given <- frm(fo, family = cox(Boundary.knots = old), data = d)
  expect_identical(as.numeric(logLik(given)), as.numeric(logLik(fit)))
  # the log-likelihood is the hand-written one at the fit's estimates,
  # through splines2's bases: an event adds log h(t), every row loses
  # H(t) and gets back H(entry), which is 0 below the lower knot
  bas <- function(f, t) {
    f(pmin(pmax(t, sp$boundary[1L]), sp$boundary[2L]), knots = sp$internal,
      Boundary.knots = sp$boundary, degree = sp$degree,
      intercept = sp$intercept)
  }
  s <- unname(cox_baseline(fit))
  eta <- as.vector(frm_linpred(fit, type = "link"))
  h0 <- drop(bas(splines2::mSpline, d$time) %*% s)
  H0 <- drop(bas(splines2::iSpline, d$time) %*% s)
  He <- ifelse(d$entry < sp$boundary[1L], 0,
               drop(bas(splines2::iSpline, d$entry) %*% s))
  ll <- sum(ifelse(d$cens == 0, log(h0) + eta, 0) - exp(eta) * (H0 - He))
  expect_lt(abs(ll - as.numeric(logLik(fit))), 1e-8 * abs(ll))
})
