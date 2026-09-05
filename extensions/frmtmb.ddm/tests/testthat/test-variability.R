## Ratcliff's across-trial variability: the analytic Wiener density
## averaged over a normal drift rate, a uniform start point and a
## uniform non-decision time.
##
## The reference below is written out here rather than taken from a
## package, because no package computes this density in a way that is
## independent of the way this one does. It is the plain Wiener density
## from RWiener, fed to R's own adaptive quadrature three times over.
## Slow, and that is the point: the thing being tested is whether fixed
## nodes reach the answer adaptive nodes give.

ref_full <- function(y, v, a, w0, t0, sv, sz, st, up = 0, tol = 1e-11) {
  dw <- function(dt, om, nu) {
    if (dt <= 0) return(0)
    RWiener::dwiener(dt + 1e-12, a, 1e-12, om, nu,
                     resp = if (up == 1) "upper" else "lower")
  }
  over_t <- function(om, nu) {
    if (st <= 0) return(dw(y - t0, om, nu))
    hi <- min(t0 + st / 2, y)
    if (hi <= t0 - st / 2) return(0)
    stats::integrate(
      function(tau) vapply(tau, function(z) dw(y - z, om, nu), 0),
      t0 - st / 2, hi, rel.tol = tol, subdivisions = 500L)$value / st
  }
  over_w <- function(nu) {
    if (sz <= 0) return(over_t(w0, nu))
    stats::integrate(
      function(om) vapply(om, function(z) over_t(z, nu), 0),
      w0 - sz / 2, w0 + sz / 2, rel.tol = tol,
      subdivisions = 500L)$value / sz
  }
  if (sv <= 0) return(over_w(v))
  stats::integrate(
    function(nu) vapply(nu, function(z) stats::dnorm(z, v, sv) * over_w(z), 0),
    v - 10 * sv, v + 10 * sv, rel.tol = tol, subdivisions = 500L)$value
}

## the density under test, with the node counts named explicitly
var_lpdf <- function(y, v, a, w, t0, sv, sz, st, up, nsz = 7L, nst = 21L) {
  nd <- ddm_nodes(c("sz", "st"), c(sz = nsz, st = nst))
  ddm_lpdf_var(y, v, a, w, t0, sv, sz, st, up, nd, TRUE, 1e-9)
}

## Values in the range the literature reports for the full model: a
## drift standard deviation of order the drift, a start-point range of
## a fifth of the boundary separation, and a non-decision-time range of
## about a third of the non-decision time.
lit <- list(
  list(y = 0.55, v = 1.2, a = 1.5, w = 0.50, t0 = 0.30,
       sv = 1.0, sz = 0.20, st = 0.10, up = 1, tag = "typical"),
  list(y = 1.80, v = 2.0, a = 1.2, w = 0.60, t0 = 0.25,
       sv = 0.5, sz = 0.10, st = 0.05, up = 1, tag = "slow"),
  list(y = 0.42, v = 0.8, a = 2.0, w = 0.50, t0 = 0.35,
       sv = 1.5, sz = 0.30, st = 0.20, up = 0, tag = "wide, range cut"),
  list(y = 0.36, v = 1.0, a = 1.5, w = 0.45, t0 = 0.30,
       sv = 1.0, sz = 0.20, st = 0.20, up = 0, tag = "range cut hard"))

test_that("variability at zero reproduces the plain Wiener density", {
  # The degenerate check. It needs no external reference and no
  # tolerance argument worth arguing about: the drift integral's closed
  # form reduces term for term at sv = 0, and a one-node rule on a
  # zero-width interval is an evaluation at its midpoint. What is left
  # is floating-point rounding on a differently associated sum.
  nd <- ddm_nodes(character(0), c(sz = 1L, st = 1L))
  set.seed(2)
  worst <- 0
  for (i in 1:200) {
    y <- stats::runif(1, 0.4, 4)
    v <- stats::runif(1, -3, 3)
    a <- stats::runif(1, 0.4, 3)
    w <- stats::runif(1, 0.15, 0.85)
    t0 <- stats::runif(1, 0.05, 0.3)
    up <- sample(0:1, 1)
    worst <- max(worst, abs(
      ddm_lpdf_var(y, v, a, w, t0, 0, 0, 0, up, nd, FALSE, 1e-9) -
        ddm_lpdf_both(y - t0, v, a, w, up)))
  }
  expect_lt(worst, 1e-13)

  # and with the quadratures switched ON but their widths driven almost
  # to nothing, so that the node loops run and still collapse to the
  # same number. The tolerance here is not slack: an average over a
  # uniform of width h differs from the value at its midpoint by
  # h^2 f'' / 24 f, which is a property of averaging and not an error in
  # the quadrature, and at the other end the surviving share of the
  # range is a difference of two endpoints either side of the
  # non-decision time and loses digits to cancellation as the width
  # shrinks. 1e-6 sits between the two.
  nd2 <- ddm_nodes(c("sz", "st"), c(sz = 7L, st = 21L))
  expect_equal(
    ddm_lpdf_var(0.9, 1.1, 1.4, 0.45, 0.25, 0, 0, 1e-6, 1, nd2, TRUE, 1e-9),
    ddm_lpdf_both(0.9 - 0.25, 1.1, 1.4, 0.45, 1),
    tolerance = 1e-9)
})

test_that("the drift integral's closed form is exact, not a quadrature", {
  skip_if_not_installed("RWiener")
  # This is the claim that saves a whole set of nodes: the drift enters
  # the log density as -v a w - v^2 t / 2 and nowhere else, so averaging
  # over a normal drift is a Gaussian integral of an
  # exponential-quadratic and completes the square. Checked against
  # adaptive quadrature of the same thing.
  worst <- 0
  for (p in list(c(0.40, 1.2, 1.4, 0.45, 0.5), c(1.10, -2.0, 0.8, 0.30, 1.5),
                 c(2.50, 0.5, 2.2, 0.60, 1.0), c(0.15, 3.0, 1.0, 0.50, 0.8),
                 c(3.00, 0.0, 1.6, 0.70, 2.0))) {
    t <- p[1]; v <- p[2]; a <- p[3]; w <- p[4]; sv <- p[5]
    num <- stats::integrate(
      function(nu) vapply(nu, function(z) {
        stats::dnorm(z, v, sv) *
          RWiener::dwiener(t + 1e-10, a, 1e-10, w, z, resp = "lower")
      }, 0), v - 12 * sv, v + 12 * sv, rel.tol = 1e-12)$value
    worst <- max(worst, abs(exp(ddm_lpdf_lower_sv(t, v, a, w, sv)) - num) / num)
  }
  # machine precision, which is what "closed form" has to mean
  expect_lt(worst, 1e-13)

  # and it reduces to the plain density at sv = 0
  expect_equal(ddm_lpdf_lower_sv(0.7, 1.3, 1.4, 0.4, 0),
               ddm_lpdf_lower(0.7, 1.3, 1.4, 0.4), tolerance = 1e-13)
})

test_that("the fixed-node density matches adaptive quadrature", {
  skip_if_not_installed("RWiener")
  for (cs in lit) {
    r <- log(with(cs, ref_full(y, v, a, w, t0, sv, sz, st, up)))
    o <- with(cs, var_lpdf(y, v, a, w, t0, sv, sz, st, up))
    expect_lt(abs(o - r), 1e-8)
  }
})

test_that("the node counts are where the defaults say they are", {
  skip_if_not_installed("RWiener")
  # The measurement behind wiener(nodes =). Both integrands are
  # analytic, so Gauss-Legendre converges geometrically, but not at the
  # same rate: the start-point integrand is smooth across its whole
  # range, while the non-decision-time integrand is cut by the response
  # time on a fast trial and turns on sharply at the cut. That is the
  # entire reason the two defaults differ.
  cs <- lit[[4]]                       # the range-cut case, the hard one
  r <- log(with(cs, ref_full(y, v, a, w, t0, sv, sz, st, up)))
  err <- function(nsz, nst) {
    abs(with(cs, var_lpdf(y, v, a, w, t0, sv, sz, st, up, nsz, nst)) - r)
  }
  # start point: 7 nodes is already at the level the non-decision-time
  # rule limits the answer to, and more buys nothing
  expect_lt(err(7L, 41L), 1e-11)
  expect_gt(err(3L, 41L), 1e-8)
  # non-decision time: the default is 21 and each step up gains orders
  expect_gt(err(15L, 7L), 1e-5)
  expect_lt(err(15L, 21L), 1e-8)
  expect_lt(err(15L, 41L), 1e-11)

  # on a trial whose range is NOT cut, the same rule is exact much
  # sooner, which is why the default is a worst case and not a typical
  # one
  cs2 <- lit[[2]]
  r2 <- log(with(cs2, ref_full(y, v, a, w, t0, sv, sz, st, up)))
  expect_lt(abs(with(cs2, var_lpdf(y, v, a, w, t0, sv, sz, st, up,
                                   5L, 7L)) - r2), 1e-12)
})

test_that("the Gauss-Legendre rule is the rule it claims to be", {
  # A quadrature nobody checks is a quadrature nobody can trust. An
  # n-point Gauss-Legendre rule is exact for polynomials of degree
  # 2n - 1 and no higher, and this one is normalized to average rather
  # than integrate, so it must reproduce the moments of the uniform.
  for (n in c(1L, 2L, 3L, 7L, 21L)) {
    g <- ddm_gauss_legendre(n)
    expect_length(g$x, n)
    expect_equal(sum(g$w), 1)
    expect_true(all(g$x > 0 & g$x < 1))
    for (k in 0:(2 * n - 1)) {
      expect_equal(sum(g$w * g$x^k), 1 / (k + 1), tolerance = 1e-12)
    }
  }
  # and it is symmetric about the midpoint, so a symmetric integrand
  # gets no spurious skew
  g <- ddm_gauss_legendre(9L)
  expect_equal(g$x, rev(1 - g$x), tolerance = 1e-14)
})

test_that("the variability density differentiates correctly", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")
  nd <- ddm_nodes(c("sz", "st"), c(sz = 7L, st = 21L))
  yv <- c(0.36, 0.55, 0.9, 1.6, 3.0)
  upv <- c(0, 1, 1, 0, 1)
  # p = (v, a, w, t0, sv, sz, st)
  lp <- function(p) {
    sum(ddm_lpdf_var(yv, p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                     upv, nd, TRUE, 1e-9))
  }
  pts <- list(c(1.2, 1.5, 0.50, 0.30, 1.0, 0.20, 0.10),
              c(-0.8, 2.0, 0.40, 0.20, 0.5, 0.30, 0.15),
              # a point where the fastest row's range IS cut, so the
              # smooth minimum that does the cutting is differentiated
              c(1.0, 1.5, 0.45, 0.30, 1.0, 0.20, 0.20))
  for (p in pts) {
    tp <- RTMB::MakeTape(lp, p)
    expect_equal(as.numeric(tp$jacobian(p)), numDeriv::grad(lp, p),
                 tolerance = 1e-6)
  }
})

test_that("check_custom_family() passes the full DDM", {
  n <- 60
  set.seed(9)
  y <- stats::runif(n, 0.5, 2.0)
  expect_true(check_custom_family(
    wiener(max_ndt = 0.4, variability = c("sv", "sz", "st")), y = y,
    dpars = list(mu = rep(0.8, n), bs = rep(1.4, n), ndt = rep(0.25, n),
                 bias = rep(0.45, n), sv = rep(0.8, n), sz = rep(0.2, n),
                 st = rep(0.08, n)),
    aterms = list(vint1 = rep(0:1, n / 2))))
})

test_that("the full DDM recovers what it was simulated from", {
  skip_if_not_installed("RWiener")
  skip_on_cran()
  # sz is left out of the estimated set on purpose and pinned by the
  # separate test below: the start-point range is the parameter the
  # literature reports as the hardest to identify, and a recovery study
  # that needs tens of thousands of trials to see it is not a unit test.
  set.seed(707)
  reps <- 8
  truth <- c(mu = 1.2, bs = 1.5, ndt = 0.30, sv = 1.0, st = 0.10)
  est <- matrix(NA_real_, reps, length(truth),
                dimnames = list(NULL, names(truth)))
  for (r in seq_len(reps)) {
    dat <- ddm_simulate(1500, mu = 1.2, bs = 1.5, ndt = 0.30, bias = 0.5,
                        sv = 1.0, st = 0.10)
    fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
               family = wiener(variability = c("sv", "st")), data = dat)
    e <- unlist(fixef(fit))
    ub <- min(dat$rt)
    est[r, ] <- c(e[["mu.(Intercept)"]],
                  exp(e[["bs.(Intercept)"]]),
                  ub / (1 + exp(-e[["ndt.(Intercept)"]])),
                  exp(e[["sv.(Intercept)"]]),
                  2 * ub / (1 + exp(-e[["st.(Intercept)"]])))
  }
  m <- colMeans(est)
  mcse <- apply(est, 2, stats::sd) / sqrt(reps)
  # the claim: every parameter's Monte Carlo mean is within three of
  # its own standard errors of the truth. Three rather than two because
  # eight replicates give the standard error itself only eight degrees
  # of freedom.
  z <- (m - truth) / mcse
  expect_true(all(abs(z) < 3),
              info = paste(names(truth), round(m, 4), round(mcse, 4),
                           round(z, 2), collapse = " | "))
})

test_that("the start-point range is recovered when there is enough data", {
  skip_if_not_installed("RWiener")
  skip_on_cran()
  # Not a fit: a likelihood profile, which separates "the parameter is
  # not identified" from "the optimizer did not find it". The start
  # point range shows up in the difference between the two boundaries'
  # fast tails, so it takes a weak drift and many trials to see, which
  # is what the literature reports and what a fit of a thousand trials
  # will not deliver.
  set.seed(13)
  dat <- ddm_simulate(12000, mu = 0.3, bs = 1.5, ndt = 0.30, bias = 0.5,
                      sv = 1.0, sz = 0.3, st = 0.10)
  nd <- ddm_nodes(c("sz", "st"), c(sz = 9L, st = 21L))
  ll <- function(sz) {
    sum(ddm_lpdf_var(dat$rt, 0.3, 1.5, 0.5, 0.30, 1.0, sz, 0.10,
                     dat$upper, nd, TRUE, 1e-9 * min(dat$rt)))
  }
  grid <- c(1e-6, 0.1, 0.2, 0.3, 0.4, 0.5)
  v <- vapply(grid, ll, 0)
  expect_equal(grid[which.max(v)], 0.3)
  expect_gt(v[grid == 0.3] - v[grid == 1e-6], 10)
})

test_that("wiener() refuses a variability spec it cannot act on", {
  expect_error(wiener(variability = "sw"), "names the across-trial")
  expect_error(wiener(variability = c("sv", "sv")), "names the across-trial")
  expect_error(wiener(variability = 1), "names the across-trial")
  expect_error(wiener(nodes = c(sz = 0)), "node counts")
  expect_error(wiener(nodes = c(sz = 2.5)), "node counts")
  expect_error(wiener(nodes = c(bogus = 7)), "node counts")
  expect_error(wiener(nodes = 7), "node counts")
  # a partial node spec keeps the other default
  f <- wiener(max_ndt = 0.4, variability = "st", nodes = c(st = 31))
  expect_setequal(f$dpars, c("mu", "bs", "ndt", "bias", "st"))
  # and the dpar order is canonical, not the order the user wrote
  g <- wiener(max_ndt = 0.4, variability = c("st", "sv"))
  expect_equal(g$dpars, c("mu", "bs", "ndt", "bias", "sv", "st"))
})

test_that("the variability parameters carry the links they claim", {
  f <- wiener(max_ndt = 0.4, variability = c("sv", "sz", "st"))
  expect_equal(f$links$sv$name, "log")
  expect_equal(f$links$sz$name, "logit")
  # st is bounded by twice the non-decision time bound, so the width
  # cannot run away to a value no response time could have produced
  lk <- f$links$st
  expect_equal(lk$name, "scaled_logit")
  expect_equal(lk$linkinv(0), 0.4)               # half of 2 * max_ndt
  expect_true(all(lk$linkinv(c(-30, -5, 0, 5, 30)) < 0.8))
  expect_equal(lk$linkinv(lk$linkfun(0.3)), 0.3, tolerance = 1e-12)

  # the logit on sz is what keeps the uniform start point inside the
  # boundaries at an unbiased start: the width is below 1, so the range
  # is inside (0, 1) whenever bias is 0.5
  expect_lt(f$links$sz$linkinv(30), 1)
  expect_gt(f$links$sz$linkinv(-30), 0)
  # both saturate past a linear predictor of about 37, as the plain
  # ndt link does and for the same reason; nothing guards it, because
  # the likelihood is long unreachable by then
  expect_equal(f$links$sz$linkinv(40), 1)
})

test_that("a fixed variability parameter is a plain fixed dpar", {
  skip_if_not_installed("RWiener")
  set.seed(88)
  dat <- ddm_simulate(300, mu = 1.0, bs = 1.4, ndt = 0.28, sv = 0.8)
  # holding sv on the response scale is the ordinary dpar mechanism,
  # not something the variability code has to know about
  fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5, sv = 0.8),
             family = wiener(variability = "sv"), data = dat)
  e <- unlist(fixef(fit))
  expect_equal(unname(e[["sv.(Intercept)"]]), log(0.8), tolerance = 1e-10)
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("the conditional mean response time accounts for variability", {
  skip_if_not_installed("RWiener")
  skip_on_cran()
  # The closed-form conditional mean is the mean of the plain model, and
  # it is the wrong number once the parameters vary between trials: the
  # trials that reached a boundary are not a fair sample of the drift
  # rates that could have produced them. The reference is the empirical
  # mean of data simulated from the model, which knows nothing about how
  # the mean is computed.
  set.seed(55)
  n <- 40000
  p <- list(mu = 0.8, bs = 1.5, ndt = 0.30, bias = 0.5,
            sv = 1.2, sz = 0.3, st = 0.10)
  dat <- ddm_simulate(n, mu = p$mu, bs = p$bs, ndt = p$ndt, bias = p$bias,
                      sv = p$sv, sz = p$sz, st = p$st)
  nd <- ddm_nodes(c("sz", "st"), c(sz = 15L, st = 21L))
  gh <- ddm_gauss_hermite(21L)

  mc <- c(NA_real_, NA_real_)
  ours <- c(NA_real_, NA_real_)
  plain <- c(NA_real_, NA_real_)
  for (b in c(0, 1)) {
    got <- dat$rt[dat$upper == b]
    mcse <- stats::sd(got) / sqrt(length(got))
    mc[b + 1] <- mean(got)
    ours[b + 1] <- ddm_mean_rt_var(p, list(dec = b), nd, gh)
    plain[b + 1] <- ddm_mean_rt(p, list(dec = b))
    expect_lt(abs(ours[b + 1] - mc[b + 1]), 4 * mcse)
  }

  # And the sharpest statement of why the extra work is not decoration.
  # At an unbiased start point the plain conditional mean decision time
  # is the SAME for both boundaries, by a symmetry of the driftless
  # process that the drift does not break. The data do not agree: with
  # a positive drift and a drift range, the trials that reached the
  # lower boundary are the ones whose own drift was small or negative,
  # and they are slower. The variability mean reproduces that gap and
  # the plain one cannot represent it at all.
  expect_equal(plain[1], plain[2], tolerance = 1e-10)
  expect_gt(mc[1] - mc[2], 0.04)
  expect_equal(ours[1] - ours[2], mc[1] - mc[2], tolerance = 0.1)
})

test_that("the simulator draws from the model the density scores", {
  skip_if_not_installed("RWiener")
  skip_on_cran()
  # sim() has to condition on the boundary the same way the mean does,
  # by reweighting which per-trial parameters the row could have had.
  # Drawing them from their own distributions and ignoring the boundary
  # would give a different distribution, and this compares against data
  # generated by an entirely separate route.
  set.seed(56)
  p <- list(mu = 0.8, bs = 1.5, ndt = 0.30, bias = 0.5,
            sv = 1.2, sz = 0.3, st = 0.10)
  ref <- ddm_simulate(20000, mu = p$mu, bs = p$bs, ndt = p$ndt,
                      bias = p$bias, sv = p$sv, sz = p$sz, st = p$st)
  nd <- ddm_nodes(c("sz", "st"), c(sz = 15L, st = 21L))
  for (b in c(0, 1)) {
    m <- 6000
    dp <- lapply(p, function(z) rep(z, m))
    drew <- ddm_sim_rt_var(dp, list(dec = rep(b, m)), m, nd)
    got <- ref$rt[ref$upper == b]
    expect_true(all(is.finite(drew)))
    se <- sqrt(stats::var(drew) / m + stats::var(got) / length(got))
    expect_lt(abs(mean(drew) - mean(got)), 4 * se)
    # the shape too, not only the location
    expect_equal(stats::quantile(drew, c(0.1, 0.5, 0.9), names = FALSE),
                 stats::quantile(got, c(0.1, 0.5, 0.9), names = FALSE),
                 tolerance = 0.05)
  }
})

test_that("fitted() on a variability fit uses the variability mean", {
  skip_if_not_installed("RWiener")
  skip_on_cran()
  set.seed(57)
  dat <- ddm_simulate(500, mu = 1.0, bs = 1.4, ndt = 0.28, sv = 1.0)
  fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
             family = wiener(variability = "sv"), data = dat)
  ft <- fitted(fit)
  expect_true(all(is.finite(ft)))
  expect_true(all(ft > 0))
  # the two boundaries get different means, which is the conditioning
  expect_false(isTRUE(all.equal(ft[dat$upper == 1][1],
                                ft[dat$upper == 0][1])))
  # and simulate() runs through the reweighted draw
  s <- stats::simulate(fit, nsim = 1)
  expect_true(all(is.finite(unlist(s))))
})
