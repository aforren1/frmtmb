## The drift-diffusion custom family formerly worked in
## vignette("case-studies") (that section now points at frmtmb.ddm,
## which supersedes the fixed-truncation density kept here). The test
## stays in core as a custom_family() regression: it exercises the
## constructor, check_custom_family() and the fit surface on a density
## defined entirely in user code.
##
## The vignette validated a research-grade density it wrote itself, so
## the density is pinned here as well: a vignette that only re-derives
## its own arithmetic proves nothing, and a change in RTMB's arithmetic
## or in the family contract has to break a test, not a page build.
##
## The reference is RWiener's C implementation, which chooses its term
## count adaptively. The version here uses a FIXED truncation, because
## the tape cannot branch on a parameter, so the pin covers both the
## agreement and the range of normalized times over which it holds.

## log first-passage density at the LOWER boundary with drift v,
## boundary separation a and relative start point w, Navarro and Fuss
## (2009) small-time series with 2K + 1 terms. Plain arithmetic only:
## no c(), no [<-, so it needs no ADoverload bindings of its own and it
## runs on advectors and on doubles alike.
wiener_lpdf1 <- function(t, v, a, w, K = 10L) {
  u <- t / a^2
  s <- 0
  for (k in -K:K) s <- s + (w + 2 * k) * exp(-(w + 2 * k)^2 / (2 * u))
  -v * a * w - v^2 * t / 2 - 2 * log(a) -
    0.5 * (log(2 * pi) + 3 * log(u)) + log(s)
}

wiener_family <- function(max_ndt) {
  force(max_ndt)
  ndt_link <- list(
    name = "scaled_logit",
    linkfun = function(mu) log(mu / (max_ndt - mu)),
    linkinv = function(eta) max_ndt / (1 + exp(-eta)),
    mu_eta = function(eta) {
      p <- 1 / (1 + exp(-eta)); max_ndt * p * (1 - p)
    })
  custom_family(
    "wiener",
    dpars = c("mu", "bs", "ndt", "bias"),
    links = list(mu = "identity", bs = "log", ndt = ndt_link,
                 bias = "logit"),
    lpdf = function(y, dpars, aterms) {
      up <- aterms$vint1
      wiener_lpdf1(y - dpars$ndt, dpars$mu * (1 - 2 * up), dpars$bs,
                   dpars$bias + up * (1 - 2 * dpars$bias))
    },
    init_dpars = list(mu = function(y, aterms) 0.5,
                      bs = function(y, aterms) 1.5,
                      ndt = function(y, aterms) 0.5 * min(y),
                      bias = function(y, aterms) 0.5),
    type = "continuous")
}

## the vignette's simulator, verbatim
r_ddm <- function(n, v, a, w, ndt, dt = 1e-4, tmax = 5) {
  v <- rep_len(v, n); x <- rep(a * w, n); live <- rep(TRUE, n)
  tt <- numeric(n); up <- integer(n); step <- 0L
  while (any(live) && step * dt < tmax) {
    step <- step + 1L
    x[live] <- x[live] + v[live] * dt +
      stats::rnorm(sum(live), 0, sqrt(dt))
    hi <- live & x >= a; lo <- live & x <= 0
    tt[hi | lo] <- step * dt; up[hi] <- 1L
    live <- live & !hi & !lo
  }
  data.frame(rt = tt + ndt, upper = up)[!live, ]
}

test_that("the fixed-truncation density matches RWiener pointwise", {
  skip_if_not_installed("RWiener")
  gr <- expand.grid(t = c(0.05, 0.15, 0.4, 0.8, 1.5, 2.5),
                    a = c(0.8, 1.4, 2.2), w = c(0.35, 0.5, 0.7),
                    v = c(-1.5, 0, 1.5))
  rel <- mapply(function(t, a, w, v) {
    ref <- log(RWiener::dwiener(t + 0.2, a, 0.2, w, v, resp = "lower"))
    abs(wiener_lpdf1(t, v, a, w) - ref) / abs(ref)
  }, gr$t, gr$a, gr$w, gr$v)
  expect_equal(length(rel), nrow(gr))
  expect_lt(max(rel), 1e-8)

  # the UPPER boundary comes from the reflection, not a second series
  up <- mapply(function(t, a, w, v) {
    ref <- log(RWiener::dwiener(t + 0.2, a, 0.2, w, v, resp = "upper"))
    abs(wiener_lpdf1(t, -v, a, 1 - w) - ref) / abs(ref)
  }, gr$t, gr$a, gr$w, gr$v)
  expect_lt(max(up), 1e-8)
})

test_that("the truncation bound is where the vignette says it is", {
  skip_if_not_installed("RWiener")
  # exact to 1e-9 out to a normalized time of 4, and past 6 it is
  # CANCELLATION rather than truncation: more terms do not help
  err <- function(u, K) {
    ref <- log(RWiener::dwiener(u + 1e-9, 1, 1e-9, 0.45, 1,
                                resp = "lower"))
    abs(wiener_lpdf1(u, 1, 1, 0.45, K) - ref) / abs(ref)
  }
  expect_lt(max(vapply(c(0.05, 0.5, 1, 2, 3, 4), err, 0, K = 10L)), 1e-9)
  expect_gt(err(8, 10L), 1e-3)
  expect_gt(err(8, 20L), 1e-3)
})

test_that("check_custom_family() passes the wiener lpdf", {
  set.seed(9)
  y <- stats::runif(50, 0.4, 1.5)
  expect_true(check_custom_family(
    wiener_family(max_ndt = 0.4), y = y,
    dpars = list(mu = rep(1, 50), bs = rep(1.6, 50), ndt = rep(0.2, 50),
                 bias = rep(0.5, 50)),
    aterms = list(vint1 = rep(0:1, 25))))
})

test_that("the fitted log likelihood is RWiener's at the same parameters", {
  skip_if_not_installed("RWiener")
  set.seed(404)
  cond <- rep(c(-1, 1), each = 350)
  dat <- r_ddm(700, v = 0.4 + 0.9 * cond, a = 1.4, w = 0.5, ndt = 0.28)
  dat$x <- cond
  fit <- frm(bf(rt | vint(upper) ~ x, bias = 0.5),
             family = wiener_family(max_ndt = min(dat$rt)), data = dat)
  e <- unlist(fixef(fit))
  ndt_hat <- min(dat$rt) / (1 + exp(-e[["ndt.(Intercept)"]]))
  drift <- e[["mu.(Intercept)"]] + e[["mu.x"]] * dat$x
  ll_ref <- sum(mapply(function(q, up, v) {
    log(RWiener::dwiener(q, exp(e[["bs.(Intercept)"]]), ndt_hat, 0.5, v,
                         resp = if (up == 1) "upper" else "lower"))
  }, dat$rt, dat$upper, drift))
  expect_equal(as.numeric(logLik(fit)), ll_ref, tolerance = 1e-8)

  # the generating parameters come back, within the sampling noise of
  # 700 trials (the boundary and the non-decision time are the sharp
  # ones; the drift intercept is the loose one)
  expect_equal(exp(e[["bs.(Intercept)"]]), 1.4, tolerance = 0.06)
  expect_equal(ndt_hat, 0.28, tolerance = 0.02)
  expect_equal(e[["mu.x"]], 0.9, tolerance = 0.15)
  # the data-bounded link keeps ndt strictly inside the support
  expect_lt(ndt_hat, min(dat$rt))
  expect_gt(ndt_hat, 0)
  # fixing bias on the RESPONSE scale puts its linear predictor at 0
  expect_equal(unname(e[["bias.(Intercept)"]]), 0)
})

test_that("the scaled-logit ndt link round-trips and keeps the support", {
  fam <- wiener_family(max_ndt = 0.5)
  lk <- fam$links[["ndt"]]
  v <- c(0.001, 0.05, 0.25, 0.45, 0.4999)
  expect_equal(lk$linkinv(lk$linkfun(v)), v, tolerance = 1e-10)
  # over the working range the value stays strictly inside (0, max_ndt)
  eta <- c(-50, -5, 0, 5, 30)
  expect_true(all(lk$linkinv(eta) > 0))
  expect_true(all(lk$linkinv(eta) < 0.5))
  # past about 37 the logit saturates in double precision and the value
  # rounds to max_ndt exactly. The density itself is what keeps the
  # optimizer away from there: it goes to zero faster than any power of
  # the decision time, so the log likelihood falls off a cliff long
  # before the link does
  expect_equal(lk$linkinv(50), 0.5)
  expect_equal(wiener_lpdf1(1e-9, 1, 1.4, 0.5), -Inf)
  expect_lt(wiener_lpdf1(1e-3, 1, 1.4, 0.5),
            wiener_lpdf1(0.5, 1, 1.4, 0.5) - 100)
  # the derivative agrees with a finite difference
  h <- 1e-6
  expect_equal(lk$mu_eta(eta),
               (lk$linkinv(eta + h) - lk$linkinv(eta - h)) / (2 * h),
               tolerance = 1e-6)
})
