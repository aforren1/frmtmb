# predict() draws the parameters from their asymptotic normal law, with
# the side law of test-predict-draw-edge.R where one parameter at a time
# shows a flat or steep side. A weakly identified mixture has no
# parameter at an end, yet its log-likelihood falls far faster than the
# normal law says along directions no single parameter shows. The draw
# check compares, at 20 test draws, the objective's drop with the drop
# the law implies there; where the median excess is above 1 the draws
# are reweighted by likelihood over proposal density (Pareto smoothed)
# and resampled. dev/predraw2-findings.md has the measurements and the
# output of these tests on 0.82.0.

# One-component poisson data fitted with a two-component mixture (the
# generator of lane ppdraw's sweep, seed 101): the check's median excess
# is 9.0, far above its threshold of 1.
pdc_mix_fit <- function() {
  set.seed(101)
  x <- stats::rnorm(300)
  invisible(stats::runif(300))
  d <- data.frame(x = x, y = stats::rpois(300, exp(1 + 0.3 * x)))
  frm(bf(y ~ x), family = mixture(poisson(), poisson()), data = d)
}

# The objective's drop from its minimum at each parameter draw predict()
# takes, under the same seed (the drawer is the first thing predict()
# takes from the stream).
pdc_drops <- function(fit, ndraws, seed = 1) {
  set.seed(seed)
  dr <- allow_warnings(
    suppressMessages(predict_par_drawer(fit, TRUE, ndraws)),
    "Pareto")
  p <- fit$opt$par
  f0 <- fit$obj$fn(p)
  vapply(seq_len(ndraws), function(s) {
    fit$obj$fn(fit_outer_vector(dr(s))) - f0
  }, 0)
}

test_that("a weakly identified mixture's draws are reweighted", {
  skip_on_cran()
  skip_if_not_installed("loo")
  fit <- pdc_mix_fit()
  said <- character(0)
  set.seed(1)
  p <- withCallingHandlers(
    predict(fit, ndraws = 1000, summary = FALSE),
    frmtmb_draw_reweighted = function(m) {
      said <<- c(said, "message")
      invokeRestart("muffleMessage")
    },
    frmtmb_draw_pareto = function(w) {
      said <<- c(said, "warning")
      invokeRestart("muffleWarning")
    })
  # the draws predict() returns sit where the likelihood is: their median
  # drop is a fraction of that at draws from the normal law of vcov(),
  # which 0.82.0 returned (no side of this fit is flat or steep, so its
  # draws were the normal law's)
  drops <- pdc_drops(fit, 1000L)
  V <- vcov(fit, full = TRUE)
  set.seed(2)
  Z <- matrix(stats::rnorm(nrow(V) * 1000), nrow(V), 1000)
  D <- fit$opt$par + crossprod(chol(V), Z)
  f0 <- fit$obj$fn(fit$opt$par)
  raw <- apply(D, 2L, function(v) fit$obj$fn(v) - f0)
  expect_lt(stats::median(drops), stats::median(raw) / 2)
  expect_false(anyNA(p))
  # either way is valid: the weights' Pareto k decides which is said
  expect_length(said, 1L)
  chk <- fit$cache[["draw_check"]]
  expect_true(chk$fire)
  expect_gt(chk$excess, 3 * draw_check_tol)
})

test_that("too few draws to judge the weights is said", {
  skip_on_cran()
  skip_if_not_installed("loo")
  fit <- pdc_mix_fit()
  set.seed(1)
  allow_warnings(suppressMessages(predict(fit, ndraws = 10)),
                 "too few to judge the weights",
                 require = "too few to judge the weights")
})

test_that("a fit the normal law describes keeps its draws and stream", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(100))
  d$y <- stats::rpois(100, exp(0.5 + 0.4 * d$x))
  fit <- frm(bf(y ~ x), family = poisson(), data = d)
  set.seed(1)
  p1 <- predict(fit, ndraws = 200, summary = FALSE)
  s1 <- .Random.seed
  chk <- fit$cache[["draw_check"]]
  expect_false(chk$fire)
  # far below the threshold, not near it
  expect_lt(abs(chk$excess), draw_check_tol / 10)
  local_mocked_bindings(draw_reweight = function(object, dw) dw$D)
  set.seed(1)
  p0 <- predict(fit, ndraws = 200, summary = FALSE)
  expect_identical(p1, p0)
  expect_identical(s1, .Random.seed)
})

test_that("many well-identified coefficients do not fire the check", {
  skip_on_cran()
  # y ~ 0 + f, poisson, 200 levels of about 10 counts each: each level's
  # log-likelihood is a little skewed, and the summed excess passes 1.
  # The first build of the check fired here, and its reweighting in 200
  # dimensions cut the coefficients' 90% coverage from 0.91 to 0.83. The
  # excess relative to the law's own drop stays small whatever the
  # number of coefficients
  set.seed(1)
  d <- data.frame(f = factor(rep(seq_len(200), each = 5)))
  d$y <- stats::rpois(1000, exp(log(2) + stats::rnorm(200, 0, 0.3))[d$f])
  fit <- frm(bf(y ~ 0 + f), family = poisson(), data = d)
  set.seed(1)
  p1 <- predict(fit, ndraws = 50, summary = FALSE)
  chk <- fit$cache[["draw_check"]]
  expect_gt(chk$excess, draw_check_tol)
  expect_lt(chk$rel, draw_check_rel / 3)
  expect_false(chk$fire)
  local_mocked_bindings(draw_reweight = function(object, dw) dw$D)
  set.seed(1)
  expect_identical(predict(fit, ndraws = 50, summary = FALSE), p1)
})

test_that("the weights carry the side law's own density", {
  skip_if_not_installed("loo")
  # a 1-D objective t^2 / 2 above 0 and (t / a)^2 / 2 below, whose lower
  # side the side law narrows exactly to its shape: the weighted share
  # below 0 is the target's a / (1 + a), up to the deviates' own split
  # between the sides. Weights anchored at each side's slope at 0 kept
  # the proposal's half (0.5065 at a = 0.3 and 0.1)
  for (a in c(0.3, 0.1)) {
    obj <- list(fn = function(p) {
      t <- p[1]
      if (t >= 0) t^2 / 2 else (t / a)^2 / 2
    }, env = new.env())
    object <- list(obj = obj, opt = list(par = 0))
    law <- list(lo = list(a * (1:4)), hi = list(1), s = 1, Lc = NULL,
                edge = integer(0), narrowed = 1L)
    set.seed(1)
    N <- 20000
    Z <- matrix(stats::rnorm(N), 1, N)
    L <- matrix(1, 1, 1)
    D <- draw_edge_apply(law, L, Z)
    r <- draw_logratio(object, list(L = L, law = law, v0 = 0,
                                    held = integer(0), off = 0L), D, Z)
    neg <- Z[1, ] < 0
    # the log ratio is log(a) on every draw below 0 and 0 above: the
    # target's density over the proposal's, a on the narrowed side
    expect_lt(max(abs(r$lw[neg] - log(a))), 1e-8 * abs(log(a)))
    expect_lt(max(abs(r$lw[!neg])), 1e-8)
    # resampled, the share below 0 is the target's a / (1 + a) within the
    # Monte Carlo error of the deviates' split between the sides and of
    # the resampling
    rs <- draw_resample(r$lw, r$stratum)
    share <- mean(D[1, rs$idx] < 0)
    nl <- sum(neg)
    q <- N * a / (nl * a + N - nl)
    se <- sqrt((2 * a / (1 + a)^2)^2 * 0.25 / N + nl * q * (1 - q) / N^2)
    expect_lt(abs(share - a / (1 + a)), 4 * se)
    # the check reads none of this as an excess: the side law follows the
    # slice exactly
    expect_lt(max(abs(r$e[1:20])), 1e-8)
  }
})

test_that("the excess is 0 where the objective is the law's quadratic", {
  # se(s) with sigma fixed: the objective is exactly quadratic in the
  # coefficients, so the drop at every draw is half its Mahalanobis
  # distance, and the excess is rounding
  set.seed(4)
  d <- data.frame(x = stats::rnorm(80), s = stats::runif(80, 0.5, 1.5))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(80, 0, d$s)
  fit <- frm(bf(y | se(s) ~ x), family = gaussian(), data = d)
  ds <- fit_draw_space(fit)
  L <- chol(ds$V)
  v0 <- fit_outer_vector(fit, ds$map)
  Z <- matrix(stats::rnorm(length(v0) * 20), length(v0), 20)
  D <- v0 + crossprod(L, Z)
  r <- draw_logratio(fit, list(L = L, law = NULL, v0 = v0, held = integer(0),
                               off = 0L), D, Z)
  drop <- r$e + colSums(Z^2) / 2
  expect_lt(max(abs(r$e)), 1e-6 * max(drop))
})

test_that("a side only just steep is narrowed only a little", {
  # 0.82.0 narrowed a side in full up to its threshold (r = v[2] / 2 at
  # 0.5) and not at all past it: the 98% width of a new level's
  # prediction stepped from 5.36 to 5.76 between neighboring data sets
  z <- c(0.5, 1, 2, 3, 5)
  full <- function(v) {
    g <- stats::approx(0:4, c(0, v), xout = pmin(z, 4))$y
    g[z > 4] <- v[4] + (z[z > 4] - 4) * max(0, v[4] - v[3])
    pmin(z, z * v[2] / 2, g)
  }
  # r = 0.35: narrowed in full, as before
  v <- c(0.4, 0.7, 1.1, 1.4)
  expect_identical(draw_edge_map(z, 1, v), full(v))
  # r = 0.499: the map is within 1% of the normal law's distance from
  # the full narrowing, so it meets the normal law at the threshold
  v <- c(0.6, 0.998, 1.7, 2.3)
  h <- full(v)
  expect_true(all(h < z))
  m <- draw_edge_map(z, 1, v)
  expect_lt(max(abs(m - z) / abs(h - z)), 0.011)
  # and halfway, r = 0.45, halfway
  v <- c(0.55, 0.9, 1.5, 2)
  m <- draw_edge_map(z, 1, v)
  expect_equal((m - z) / (full(v) - z), rep(0.5, length(z)))
})

test_that("the Jacobian the check reads is the side law map's slope", {
  # a steep side's law: knots where the log-likelihood falls by z^2 / 2
  for (v in list(c(0.3, 0.5, 0.6, 0.7), c(0.05, 0.08, 0.1, 0.11),
                 c(0.9, 1.6, 2.1, 2.4))) {
    z <- c(0.1, 0.7, 1.3, 2.5, 3.2, 4.5, 6)
    h <- 1e-6
    num <- (draw_edge_map(z + h, 1, v) - draw_edge_map(z - h, 1, v)) / (2 * h)
    expect_equal(draw_edge_slope(z, v), num, tolerance = 1e-6)
  }
  expect_identical(draw_edge_slope(c(0.5, 2), 1), c(1, 1))
})
