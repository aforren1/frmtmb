## The family as frmtmb sees it: does it tape, does it fit, does it
## recover, and does it refuse the things it cannot do.

test_that("check_custom_family() passes the wiener density", {
  set.seed(9)
  y <- stats::runif(60, 0.5, 2.0)
  expect_true(check_custom_family(
    wiener(max_ndt = 0.4), y = y,
    dpars = list(mu = rep(0.8, 60), bs = rep(1.4, 60),
                 ndt = rep(0.25, 60), bias = rep(0.45, 60)),
    aterms = list(vint1 = rep(0:1, 30))))
})

test_that("the AD gradient agrees with finite differences", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")
  tt <- c(0.35, 0.6, 0.9, 1.4, 2.6, 5.0)
  upv <- c(0, 1, 0, 1, 1, 0)
  lp <- function(p) sum(ddm_lpdf_both(tt - p[4], p[1], p[2], p[3], upv))
  # the points span normalized times on both sides of the blend band,
  # so the smooth switch itself is differentiated here
  pts <- list(c(0.8, 1.4, 0.5, 0.2), c(-1.2, 0.9, 0.35, 0.15),
              c(0.0, 2.4, 0.7, 0.30), c(2.5, 0.6, 0.45, 0.25),
              c(-0.3, 3.2, 0.2, 0.05))
  for (p in pts) {
    tp <- RTMB::MakeTape(lp, p)
    g_ad <- as.numeric(tp$jacobian(p))
    g_fd <- numDeriv::grad(lp, p)
    expect_equal(g_ad, g_fd, tolerance = 1e-7)
  }
})

test_that("a two-condition fit recovers the generating parameters", {
  skip_if_not_installed("RWiener")
  set.seed(404)
  cond <- rep(c(0, 1), each = 350)
  dat <- ddm_simulate(700, mu = 0.4 + 0.9 * cond, bs = 1.4, ndt = 0.28)
  dat$cond <- cond
  fit <- frm(bf(rt | vint(upper) ~ cond, bias = 0.5),
             family = wiener(), data = dat)
  e <- unlist(fixef(fit))
  ub <- min(dat$rt)
  ndt_hat <- ub / (1 + exp(-e[["ndt.(Intercept)"]]))

  # within the sampling noise of 700 trials: the boundary and the
  # non-decision time are the sharp ones, the drift intercept is loose
  expect_equal(exp(e[["bs.(Intercept)"]]), 1.4, tolerance = 0.06)
  expect_equal(ndt_hat, 0.28, tolerance = 0.02)
  expect_equal(e[["mu.cond"]], 0.9, tolerance = 0.15)
  expect_equal(e[["mu.(Intercept)"]], 0.4, tolerance = 0.25)

  # the data-bounded link keeps ndt strictly inside the support
  expect_lt(ndt_hat, ub)
  expect_gt(ndt_hat, 0)
  # fixing bias on the RESPONSE scale puts its linear predictor at 0
  expect_equal(unname(e[["bias.(Intercept)"]]), 0)

  # and the fitted log likelihood is RWiener's at the same parameters
  drift <- e[["mu.(Intercept)"]] + e[["mu.cond"]] * dat$cond
  ll_ref <- sum(mapply(function(q, up, v) {
    log(RWiener::dwiener(q, exp(e[["bs.(Intercept)"]]), ndt_hat, 0.5, v,
                         resp = if (up == 1) "upper" else "lower"))
  }, dat$rt, dat$upper, drift))
  expect_equal(as.numeric(logLik(fit)), ll_ref, tolerance = 1e-8)
})

test_that("the non-decision-time bound is found from the data", {
  skip_if_not_installed("RWiener")
  set.seed(12)
  dat <- ddm_simulate(200, mu = 0.9, bs = 1.3, ndt = 0.25)
  # no max_ndt argument: valid_y fills it in from min(rt)
  fit <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(),
             data = dat)
  e <- unlist(fixef(fit))
  ndt_hat <- min(dat$rt) / (1 + exp(-e[["ndt.(Intercept)"]]))
  expect_lt(ndt_hat, min(dat$rt))
  expect_equal(ndt_hat, 0.25, tolerance = 0.05)

  # and an explicit bound gives the same answer when it is the same bound
  fit2 <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5),
              family = wiener(max_ndt = min(dat$rt)), data = dat)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(fit2)),
               tolerance = 1e-6)
})

test_that("the scaled-logit ndt link round-trips and keeps the support", {
  fam <- wiener(max_ndt = 0.5)
  lk <- fam$links[["ndt"]]
  v <- c(0.001, 0.05, 0.25, 0.45, 0.4999)
  expect_equal(lk$linkinv(lk$linkfun(v)), v, tolerance = 1e-10)
  eta <- c(-50, -5, 0, 5, 30)
  expect_true(all(lk$linkinv(eta) > 0))
  expect_true(all(lk$linkinv(eta) < 0.5))
  # the derivative agrees with a central difference
  h <- 1e-6
  expect_equal(lk$mu_eta(eta),
               (lk$linkinv(eta + h) - lk$linkinv(eta - h)) / (2 * h),
               tolerance = 1e-6)
  # a link used before frm() has seen data says so rather than
  # returning a silent NA
  expect_error(wiener()$links[["ndt"]]$linkinv(0),
               "bound is not set yet")
})

test_that("the missing decision indicator is refused, not ignored", {
  skip_if_not_installed("RWiener")
  set.seed(5)
  dat <- ddm_simulate(80, mu = 0.6, bs = 1.3, ndt = 0.2)
  # frmtmb has no way for a family to declare a required addition term,
  # so this refusal is the package's own. Without it the density reads
  # a NULL and the log likelihood collapses to a sum over no rows.
  expect_error(frm(bf(rt ~ 1, bias = 0.5), family = wiener(), data = dat),
               "decision indicator is missing")
  # the message points at the brms spelling the user probably tried
  expect_error(frm(bf(rt ~ 1, bias = 0.5), family = wiener(), data = dat),
               "dec\\(decision\\)")
})

test_that("a decision indicator that is not 0/1 is refused", {
  skip_if_not_installed("RWiener")
  set.seed(5)
  dat <- ddm_simulate(80, mu = 0.6, bs = 1.3, ndt = 0.2)
  dat$bad <- dat$upper + 1
  expect_error(frm(bf(rt | vint(bad) ~ 1, bias = 0.5), family = wiener(),
                   data = dat),
               "must be 0 \\(lower boundary\\) or 1")
})

test_that("a non-decision-time bound at or above min(rt) is refused", {
  skip_if_not_installed("RWiener")
  set.seed(5)
  dat <- ddm_simulate(80, mu = 0.6, bs = 1.3, ndt = 0.2)
  expect_error(frm(bf(rt | vint(upper) ~ 1, bias = 0.5),
                   family = wiener(max_ndt = 5), data = dat),
               "above the smallest response time")
  expect_error(wiener(max_ndt = -1), "one positive finite number")
  # min(rt) itself is fine, and is what the default picks: the scaled
  # logit never reaches its bound at a finite linear predictor
  expect_no_error(frm(bf(rt | vint(upper) ~ 1, bias = 0.5),
                      family = wiener(max_ndt = min(dat$rt)), data = dat))
})

test_that("a non-positive response is refused", {
  set.seed(5)
  dat <- data.frame(rt = c(0.5, 0.7, -0.1), upper = c(1, 0, 1))
  expect_error(frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(),
                   data = dat),
               "strictly positive, finite response time")
})

test_that("the compat rows this package registers are present", {
  tb <- frm_compat("wiener")
  expect_true(nrow(tb) > 0)
  expect_true("vint()" %in% tb$feature_b)
  expect_equal(tb$status[tb$feature_b == "vint()"], "works")
  expect_equal(tb$status[tb$feature_b == "trunc()"], "refused")
  expect_true("wiener" %in% frm_compat_features()$name)
})
