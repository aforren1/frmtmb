# The warm start of an exact gp() beside brms's cov = FALSE ar(): the model's
# marginal with the AR term read exactly (R/gp-start.R, gp_marginal_ar()).
# The reference density is written from the formula with no frmtmb code in it.
# Design: dev/arstep-findings.md.

# y - beta = f + e: f the field, sdgp^2 times the kernel plus the relative
# nugget; e the AR(p) residual of brms's cov = FALSE form, started at zero in
# each group, so A e is white noise for A = I - sum_k ar_k (lag k). Rows are in
# time order within each group.
ga_ref_nll <- function(y, x, g, beta, sigma, sdgp, lscale, ar) {
  n <- length(y)
  A <- diag(n)
  for (grp in unique(g)) {
    r <- which(g == grp)
    for (k in seq_along(ar)) {
      for (i in seq_along(r)) if (i > k) A[r[i], r[i - k]] <- -ar[k]
    }
  }
  u <- sort(unique(x))
  J <- match(x, u)
  K <- exp(-outer(u, u, "-")^2 / (2 * lscale^2))
  diag(K) <- diag(K) + 1e-12
  Ai <- solve(A)
  V <- sdgp^2 * K[J, J] + sigma^2 * Ai %*% t(Ai)
  R <- chol(V)
  e <- backsolve(R, y - beta, transpose = TRUE)
  sum(log(diag(R))) + 0.5 * sum(e^2) + 0.5 * n * log(2 * pi)
}

# the review's step-shaped field beside ar(): seed 74000 + 31 seed + n
ga_step <- function(n, seed, p = 1L, groups = 1L) {
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  ar <- if (p == 1L) 0.4 else c(0.4, 0.2)
  data.frame(x = x, t = rep(seq_len(n / groups), groups),
             g = rep(seq_len(groups), each = n / groups),
             y = 2 * (x > 5) - 1 +
               as.numeric(stats::arima.sim(list(ar = ar), n, sd = 0.03)))
}

ga_natural <- function(fit, p) {
  pl <- fit$obj$env$parList(p)
  list(beta = pl$beta[[1L]], sigma = exp(pl$betad[[1L]]),
       sdgp = exp(pl$theta[[1L]]), lscale = exp(pl$theta[[2L]]),
       ar = pl$thetaac)
}

test_that("the AR marginal reads the model's own density", {
  skip_on_cran()
  cases <- list(
    list(d = ga_step(40L, 2L), f = y ~ gp(x) + ar(time = t), g = FALSE),
    list(d = ga_step(40L, 1L, p = 2L), f = y ~ gp(x) + ar(time = t, p = 2),
         g = FALSE),
    list(d = ga_step(48L, 1L, groups = 3L),
         f = y ~ gp(x) + ar(time = t, gr = g), g = TRUE))
  for (cs in cases) {
    fit <- frm(cs$f, data = cs$d, family = gaussian())
    s <- frmtmb:::gp_marginal_setup(fit$frame)
    expect_false(is.null(s))
    nll <- frmtmb:::gp_marginal_nll(s)
    p0 <- fit$opt$par
    iac <- which(names(p0) == "thetaac")
    for (k in 0:2) {
      p <- p0 + k * 0.2 * cos(seq_along(p0) + k)
      # away from the optimum, coefficients whose inverse of A stays small:
      # past a unit root A^-1 grows like ar^n and the formula's solve() loses
      # the digits the comparison needs (1e-3 at ar 1.74, n 40)
      if (k > 0L) p[iac] <- c(0.3, -0.4)[k] * 0.5^(seq_along(iac) - 1L)
      v <- fit$obj$fn(p)
      nat <- ga_natural(fit, p)
      grp <- if (cs$g) cs$d$g else rep(1L, nrow(cs$d))
      ref <- ga_ref_nll(cs$d$y, cs$d$x, grp,
                        nat$beta, nat$sigma, nat$sdgp, nat$lscale, nat$ar)
      pl <- fit$obj$env$parList(p)
      m <- nll(list(beta = pl$beta, betad = pl$betad, theta = pl$theta,
                    thetaac = pl$thetaac))
      # three readings of one number: the tape, the formula, the marginal
      expect_lt(abs(v - ref), 1e-8 * abs(ref))
      expect_lt(abs(m - ref), 1e-8 * abs(ref))
    }
  }
})

test_that("a step beside ar() reaches the maximum, not the field's basin", {
  skip_on_cran()
  # 0.79.1 stopped at ll 130.31 with conv 0 and no warning, ar 0.38; the
  # maximum is 156.17 with ar 1.008, a unit root that takes the step
  d <- ga_step(160L, 1L)
  fit <- frm(y ~ gp(x) + ar(time = t), data = d, family = gaussian())
  ref <- frm(y ~ gp(x) + ar(time = t), data = d, family = gaussian(),
             start = list(theta = c(log(1.5), log(0.1))))
  ll <- as.numeric(logLik(fit))
  lr <- as.numeric(logLik(ref))
  expect_gte(ll, lr - 1e-6 * abs(lr))
  # the start came from the marginal, which the objective agreed with
  expect_false(is.null(fit$obj$env$.frmtmb_gp_marginal))
  # a row that needs the scan's unit-root cell: from AR coefficient 0
  # alone the start reached 15.47 and the maximum is 28.79
  d <- ga_step(80L, 22L)
  fit <- frm(y ~ gp(x) + ar(time = t), data = d, family = gaussian())
  ref <- frm(y ~ gp(x) + ar(time = t), data = d, family = gaussian(),
             start = list(theta = c(log(0.4782), log(0.1701)),
                          thetaac = 1.022))
  ll <- as.numeric(logLik(fit))
  lr <- as.numeric(logLik(ref))
  expect_gte(ll, lr - 1e-6 * abs(lr))
})

test_that("the marginal reads only the AR forms it states", {
  skip_on_cran()
  d <- ga_step(40L, 2L)
  setup <- function(f, ...) {
    fit <- frm(f, data = d, family = gaussian(), ...)
    frmtmb:::gp_marginal_setup(fit$frame)
  }
  # read: pure AR, cov = FALSE
  expect_false(is.null(setup(y ~ gp(x) + ar(time = t))))
  # not read: a moving-average part, the stationary covariance form, a row's
  # own weight
  expect_null(setup(y ~ gp(x) + ma(time = t)))
  expect_null(setup(y ~ gp(x) + arma(time = t)))
  expect_null(setup(y ~ gp(x) + ar(time = t, cov = TRUE)))
  d$w <- rep(c(1, 2), length.out = nrow(d))
  expect_null(setup(y | weights(w) ~ gp(x) + ar(time = t)))
  # and with no residual term the marginal is the one it was
  s <- setup(y ~ gp(x))
  expect_false(is.null(s))
  expect_null(s$ar)
})

test_that("an AR marginal the objective does not agree with is not used", {
  skip_on_cran()
  d <- ga_step(40L, 2L)
  honest <- frm(y ~ gp(x) + ar(time = t), data = d, family = gaussian())
  expect_false(is.null(honest$obj$env$.frmtmb_gp_marginal))
  real <- frmtmb:::gp_marginal_max
  local_mocked_bindings(gp_marginal_max = function(frame, ...) {
    m <- real(frame, ...)
    if (!is.null(m)) m$value <- m$value + 1
    m
  }, .package = "frmtmb")
  off <- frm(y ~ gp(x) + ar(time = t), data = d, family = gaussian())
  # the working start ran instead, and nothing claims the marginal's proof
  expect_null(off$obj$env$.frmtmb_gp_marginal)
  expect_true(is.finite(as.numeric(logLik(off))))
})

test_that("off the gaussian, the second maximum of the marginal is a start", {
  skip_on_cran()
  # the review's bernoulli row (dev/review-ncgp-r3-ng.R, seed 1, n 40): 0.79.1
  # stopped 1.07 short at a length scale of 0.012, conv 0 and no warning
  n <- 40L
  set.seed(74000 + 31 + n)
  x <- sort(round(stats::runif(n, 0, 10), 4))
  eta <- 2 * sin(3 * x) + 0.5 * x
  d <- data.frame(x = x,
                  y = stats::rbinom(n, 1, stats::plogis(eta - mean(eta))))
  fit <- frm(y ~ gp(x), data = d, family = bernoulli())
  ref <- frm(y ~ gp(x), data = d, family = bernoulli(),
             start = list(theta = c(0.18101, -0.30611)))
  w <- frmtmb:::gp_start_working(fit$frame)
  expect_gt(length(w$alts), 0L)
  ll <- as.numeric(logLik(fit))
  lr <- as.numeric(logLik(ref))
  expect_gte(ll, lr - 1e-6 * abs(lr))
  # the objective is left at the optimum the fit reports
  env <- fit$obj$env
  outer <- env$last.par.best[-env$random]
  expect_equal(unname(as.numeric(outer)), unname(as.numeric(fit$opt$par)))
  expect_equal(as.numeric(env$value.best), as.numeric(fit$opt$objective))
  # a gaussian response takes the model's own marginal, whose start returns
  # before any second start is read
  g <- data.frame(x = x, y = eta + 0.1 * stats::rnorm(n))
  fg <- frm(y ~ gp(x), data = g, family = gaussian())
  expect_false(is.null(fg$obj$env$.frmtmb_gp_marginal))
})

test_that("a second descent that loses leaves the modes at the kept optimum", {
  skip_on_cran()
  # the ncgp review's Gamma row (seed 1, n 40, curve sin(x) + 0.4 sin(10x)):
  # the second maximum of the working marginal is 18.8 below its best, and
  # the descent from it ends about 20 below the first. Without obj_settle()
  # the modes and fitted values read that descent's end
  n <- 40L
  set.seed(74000 + 31 + n)
  x <- sort(round(stats::runif(n, 0, 10), 4))
  eta <- sin(x) + 0.4 * sin(10 * x)
  d <- data.frame(x = x, y = stats::rgamma(n, 400, 400 / exp(eta)))
  fit <- frm(y ~ gp(x), data = d, family = Gamma(link = "log"))
  # the modes at the reported optimum, solved afresh from the same point
  ref <- frm(y ~ gp(x), data = d, family = Gamma(link = "log"),
             start = list(theta = fit$estimates$theta,
                          beta = fit$estimates$beta))
  b <- fit$estimates$b
  br <- ref$estimates$b
  # whether the second descent ran, won or lost here, the modes and the
  # fitted values must be those of the optimum the fit reports
  expect_lt(max(abs(b - br)), 1e-4 * max(abs(br)))
  fe <- function(f) {
    v <- fitted(f)
    as.numeric(if (is.matrix(v) || is.data.frame(v)) v[, 1] else v)
  }
  f1 <- fe(fit)
  f2 <- fe(ref)
  expect_lt(max(abs(f1 - f2)), 1e-4 * max(abs(f2)))
  env <- fit$obj$env
  expect_equal(unname(env$last.par[-env$random]), unname(fit$opt$par))
})

test_that("a constant sigma is a known noise in the marginal", {
  skip_on_cran()
  set.seed(12)
  n <- 50L
  x <- sort(stats::runif(n, 0, 10))
  d <- data.frame(x = x, t = seq_len(n), y = sin(x) + 0.3 * stats::rnorm(n))
  # without ar(): the marginal reads the objective at the fit's point. Read
  # as no noise, it stood 1.6e12 from it, and its maximum was the start and
  # the number the convergence check compares the stop with
  fit <- frm(bf(y ~ gp(x), sigma = 0.3), data = d)
  s <- frmtmb:::gp_marginal_setup(fit$frame)
  expect_false(is.null(s))
  pl <- fit$obj$env$parList(fit$opt$par)
  m <- frmtmb:::gp_marginal_nll(s)(list(beta = pl$beta, theta = pl$theta))
  v <- fit$obj$fn(fit$opt$par)
  expect_lt(abs(m - v), 1e-8 * abs(v))
  # beside ar() the AR marginal reads a free sigma only: no scan, and no
  # warning from one (an escaped warning fails this test)
  fa <- frm(bf(y ~ gp(x) + ar(time = t), sigma = 0.3), data = d)
  expect_null(frmtmb:::gp_marginal_setup(fa$frame))
  expect_true(is.finite(as.numeric(logLik(fa))))
})

test_that("on a panel the AR marginal factors positions, not rows", {
  skip_on_cran()
  # 25 series of 12 times on 12 shared positions: 300 rows. Factoring the
  # rows took 180 s on 1000 rows of 40 positions, against 0.4 s for the
  # working start
  set.seed(9602)
  S <- 25L
  Tn <- 12L
  d <- expand.grid(t = seq_len(Tn), s = factor(seq_len(S)))
  d$x <- d$t / 2
  d$y <- sin(d$x) + unlist(lapply(seq_len(S), function(i) {
    as.numeric(stats::arima.sim(list(ar = 0.5), Tn, sd = 0.2))
  }))
  fit <- frm(y ~ gp(x) + ar(time = t, gr = s), data = d, family = gaussian())
  s <- frmtmb:::gp_marginal_setup(fit$frame)
  expect_identical(s$u, Tn)
  # nothing the marginal holds is rows by rows
  dims <- function(o) {
    if (is.matrix(o)) return(list(dim(o)))
    if (is.list(o) && !is.data.frame(o)) return(do.call(c, lapply(o, dims)))
    list()
  }
  big <- vapply(dims(s[setdiff(names(s), "frame")]), function(dm) {
    all(dm >= nrow(d))
  }, NA)
  expect_false(any(big))
  expect_false(is.null(fit$obj$env$.frmtmb_gp_marginal))
})
