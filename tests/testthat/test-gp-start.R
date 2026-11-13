# Where the descent of a fit with an exact gp() begins (R/gp-start.R), and what
# follows from it: the maximum on a trend plus a wiggle, y's and x's units, the
# warning at the nugget, and standard errors in units. References are written
# from the formula with no frmtmb code in them. Designs: dev/ncgp-findings.md,
# fix round 1.

# the profile of the marginal over length scale and noise ratio
gs_prof <- function(y, x, q) {
  u <- sort(unique(x))
  J <- match(x, u)
  n <- length(y)
  C <- exp(-outer(u, u, "-")^2 / (2 * exp(2 * q[1])))
  diag(C) <- diag(C) + 1e-12
  V <- C[J, J]
  diag(V) <- diag(V) + exp(q[2])
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(1e12)
  wy <- backsolve(R, y, transpose = TRUE)
  w1 <- backsolve(R, rep(1, n), transpose = TRUE)
  e <- wy - w1 * sum(w1 * wy) / sum(w1^2)
  s2 <- sum(e^2) / n
  sum(log(diag(R))) + 0.5 * n * (log(2 * pi * s2) + 1)
}

# the marginal's maximum from a grid and from the fit's own point, so a fit
# cannot pass by a reference stopping short of it
gs_max <- function(y, x, fit = NULL) {
  d <- diff(range(x))
  grid <- expand.grid(l = log(d / 2^(0:9)), r = log(10^(1:-9)))
  v <- mapply(function(a, b) gs_prof(y, x, c(a, b)), grid$l, grid$r)
  st <- lapply(order(v)[1:3], function(i) unlist(grid[i, ]))
  if (!is.null(fit)) {
    th <- fit$estimates$theta
    st <- c(st, list(c(th[2L], 2 * (fit$estimates$betad[[1L]] - th[1L]))))
  }
  min(vapply(st, function(q0) {
    stats::nlminb(q0, function(q) gs_prof(y, x, q),
                  control = list(rel.tol = 1e-12))$objective
  }, 0))
}

gs_fit <- function(d, ...) {
  allow_warnings(frm(bf(y ~ gp(x)), data = d, family = gaussian(), ...),
                 c("Optimizer did not report", "gradient", "noise has reached",
                   "Standard errors are not"))
}

gs_curves <- list(
  mix = function(x) 2 * sin(3 * x) + x,
  mix2 = function(x) sin(4 * x) + 0.5 * x,
  mix3 = function(x) 3 * cos(2 * x) - 2 * x,
  mix4 = function(x) sin(5 * x) + x^2 / 10,
  two = function(x) sin(x) + 0.3 * sin(10 * x),
  sin10 = function(x) sin(10 * x),
  sin3 = function(x) sin(3 * x)
)

test_that("frm() reaches the maximum on a trend plus a wiggle", {
  # 95 of 168 rows of this set stopped short by more than 1e-3 nats, worst 499,
  # conv 0 with no warning; 0.78.0 stopped short on 57
  for (nm in names(gs_curves)) for (s in c(0.3, 0.05, 1e-2)) {
    set.seed(31047)
    x <- sort(round(stats::runif(40L, 0, 10), 4))
    y <- gs_curves[[nm]](x) + s * stats::rnorm(40L)
    fit <- gs_fit(data.frame(x = x, y = y))
    gap <- as.numeric(logLik(fit)) + gs_max(y, x, fit)
    expect_gt(gap, -1e-3)
  }
})

test_that("a trend the field also explains reaches the end of its ridge", {
  # exp(0.8 x): the intercept and the sd move along a ridge, and a descent that
  # carries them stops short: the last two rows read 0.2 and 1.1 nats short
  # without the profile's descent, and 0.002 with it, 0.006 under OpenBLAS. The
  # profile's simplex from the fit's own point is the reference; the value rounds
  # at 1e-3. The bound is 3 times the largest reading and a tenth of the smallest
  # shortfall
  rows <- list(c(5701, 100, 0.3), c(5811, 100, 0.3), c(5751, 40, 0.01))
  for (r in rows) {
    set.seed(r[[1L]])
    x <- sort(stats::runif(r[[2L]], 0, 10))
    y <- exp(0.8 * x) + r[[3L]] * stats::rnorm(r[[2L]])
    fit <- gs_fit(data.frame(x = x, y = y))
    th <- fit$estimates$theta
    q0 <- c(th[2L], 2 * (fit$estimates$betad[[1L]] - th[1L]))
    ref <- stats::optim(q0, function(q) gs_prof(y, x, q), method = "Nelder-Mead",
                        control = list(maxit = 600L, reltol = 1e-14))$value
    expect_gt(as.numeric(logLik(fit)) + ref, -0.02)
  }
})

test_that("at low noise the field holds the wiggle and sigma the noise", {
  # sigma read 1.19 at a true 1e-3, ll -133.95 against 282.84, conv 0, silent
  set.seed(1)
  x <- sort(stats::runif(80L, 0, 10))
  y <- 2 * sin(3 * x) + x + 1e-3 * stats::rnorm(80L)
  fit <- gs_fit(data.frame(x = x, y = y))
  expect_lt(abs(log(exp(fit$estimates$betad[[1L]]) / 1e-3)), log(1.5))
})

test_that("an exact gp() fit is y's and x's up to the optimizer's tolerance", {
  # y * 1000 on sin(6 x) stopped 66 to 103 nats short with an empty field
  designs <- list(
    list(seed = 4140, n = 40L, to = 6, f = function(x) sin(6 * x), s = 1e-3),
    list(seed = 1, n = 80L, to = 10, s = 1e-3,
         f = function(x) 2 * sin(3 * x) + x),
    list(seed = 2, n = 60L, to = 10, s = 0.1,
         f = function(x) sin(4 * x) + 0.5 * x),
    list(seed = 20261102, n = 60L, to = 6, f = function(x) sin(2 * x), s = 1e-3)
  )
  for (ds in designs) {
    set.seed(ds$seed)
    x <- round(stats::runif(ds$n, 0, ds$to), 3)
    y <- ds$f(x) + ds$s * stats::rnorm(ds$n)
    nd <- seq(min(x) + 0.05, max(x) - 0.07, length.out = 9L)
    read <- function(ys, xs) {
      f <- gs_fit(data.frame(x = x * xs, y = y * ys))
      ft <- fitted(f, newdata = data.frame(x = nd * xs))
      list(ll = as.numeric(logLik(f)) + ds$n * log(ys),
           sigma = exp(f$estimates$betad[[1L]]) / ys,
           sd = exp(f$estimates$theta[[1L]]) / ys,
           rho = exp(f$estimates$theta[[2L]]) / xs,
           fit = ft[, 1L] / ys, se = ft[, 2L] / ys)
    }
    b <- read(1, 1)
    for (sc in list(c(1000, 1), c(1e-3, 1), c(1, 1000), c(1, 1e-3))) {
      a <- read(sc[1L], sc[2L])
      expect_lt(abs(a$ll - b$ll), 1e-5 * abs(b$ll))
      expect_lt(abs(a$sigma / b$sigma - 1), 1e-3)
      expect_lt(abs(a$sd / b$sd - 1), 1e-2)
      expect_lt(abs(a$rho / b$rho - 1), 1e-3)
      expect_lt(max(abs(a$fit - b$fit) / b$se), 1e-2)
      expect_true(all(is.finite(a$se)))
      expect_lt(max(abs(a$se / b$se - 1)), 1e-2)
    }
  }
})

test_that("a reader's draw of a far length scale reads a field", {
  # the warm start reached a long, flat field on one level of by = f, whose
  # length scale has a large standard error: predict()'s draw underflowed, and
  # the error took every row
  set.seed(20261106)
  n <- 60L
  d <- data.frame(x = round(stats::runif(n, 0, 6), 3),
                  f = factor(rep(c("a", "b", "c"), length.out = n)))
  stats::runif(2L * n)
  d$y <- sin(2 * d$x) + 0.2 * stats::rnorm(n)
  fit <- allow_warnings(frm(bf(y ~ gp(x, by = f, cmc = FALSE)), data = d),
                        c("Optimizer did not report", "gradient"))
  bk <- fit$frame$re_blocks[[1L]]
  L <- frmtmb:::gp_factor(c(0, -800), bk)
  expect_true(all(is.finite(L)))
  expect_lt(max(abs(tcrossprod(L) - diag(1 + 1e-12, nrow(L)))), 1e-10)
  expect_no_error(allow_warnings(stats::predict(fit, ndraws = 20L),
                                 "standard error"))
})

# the best of a few starts on the fit's own tape, for a model with no closed form
gs_best <- function(args, fit) {
  best <- as.numeric(logLik(fit))
  th0 <- fit$estimates$theta
  for (rho0 in c(0.15, 0.5, 1.5)) for (sd0 in c(0.5, 3)) {
    th <- th0
    for (b in fit$frame$re_blocks) {
      if (identical(b$covstruct, "gp")) {
        th[b$theta_idx[1L]] <- log(sd0)
        th[b$theta_idx[-1L]] <- log(rho0)
      }
    }
    a <- args
    a$start <- c(a$start, list(theta = th))
    v <- tryCatch(as.numeric(logLik(suppressWarnings(do.call(frm, a)))),
                  error = function(e) NA_real_)
    if (isTRUE(v > best)) best <- v
  }
  best
}

gs_edge <- local({
  set.seed(58001)
  x <- sort(stats::runif(60L, 0, 10))
  z <- stats::rnorm(60L)
  zm <- z
  zm[seq(5L, 60L, 12L)] <- NA
  d <- data.frame(x, z, zm, zo = z + stats::rnorm(60L, sd = 0.1), s = 0.1,
                  y = 2 * sin(3 * x) + x + 0.3 * z + 0.1 * stats::rnorm(60L))
  d$o <- ordered(sample(1:4, 60L, TRUE))
  d
})

test_that("what the marginal's gate refuses reaches the maximum all the same", {
  # conv 0 and silent, on a curve 0.78.0 reached: 50, 52, 94, 102 and 52 nats
  # short (dev/ncgp-findings.md, round 2)
  specs <- list(
    list(y ~ 0 + gp(x)),
    list(y ~ gp(x), start = list(beta = 5)),
    list(bf(y ~ mi(zm) + gp(x)) + bf(zm | mi() ~ 1) + set_rescor(FALSE)),
    list(y ~ me(zo, s) + gp(x)),
    list(y ~ mo(o) + gp(x)),
    # a field in a parameter of a nonlinear body: 100.4 nats short on 0.79.1
    # (dev/fitedge-findings.md, item 1)
    list(bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE))
  )
  for (sp in specs) {
    a <- c(sp[1L], list(data = gs_edge), sp[-1L])
    fit <- suppressWarnings(do.call(frm, a))
    expect_gt(as.numeric(logLik(fit)) - gs_best(a, fit), -1e-3)
  }
})

test_that("a link-scale response starts where its working response reads", {
  # 6 of 6 negbinomial rows short, worst 167, and the like on these
  rows <- list(
    list(y = function(e, n) stats::rnbinom(n, size = 20, mu = exp(4 + e)),
         fam = negbinomial()),
    list(y = function(e, n) stats::rgamma(n, 400, 400 / exp(e)),
         fam = Gamma(link = "log")),
    list(y = function(e, n) exp(e + stats::rnorm(n, sd = 0.02)), fam = lognormal()),
    list(y = function(e, n) e + 0.02 * stats::rt(n, 5), fam = student())
  )
  for (r in rows) {
    set.seed(56073)
    x <- sort(round(stats::runif(60L, 0, 10), 4))
    d <- data.frame(x = x, y = r$y(2 * sin(3 * x) + 0.5 * x, 60L))
    a <- list(y ~ gp(x), data = d, family = r$fam)
    fit <- suppressWarnings(do.call(frm, a))
    expect_gt(as.numeric(logLik(fit)) - gs_best(a, fit), -1e-3)
  }
})

test_that("the working response takes each family's edge", {
  fr <- function(...) suppressWarnings(frm(...))$frame
  set.seed(3)
  x <- stats::runif(30L, 0, 6)
  lp1 <- function(f) f$linpreds[[1L]]
  d <- data.frame(x = x, y = stats::rpois(30L, 0.3), t = 5L)
  f <- fr(y ~ gp(x), data = d, family = poisson())
  expect_true(all(is.finite(frmtmb:::gp_working_y(f, lp1(f)))))
  d$b <- stats::rbinom(30L, 5L, 0.5)
  f <- fr(b | trials(t) ~ gp(x), data = d, family = binomial())
  expect_true(all(is.finite(frmtmb:::gp_working_y(f, lp1(f)))))
  d$ber <- stats::rbinom(30L, 1L, 0.5)
  f <- fr(ber ~ gp(x), data = d, family = bernoulli())
  expect_true(all(is.finite(frmtmb:::gp_working_y(f, lp1(f)))))
  # an ordinal response has no such reading, and takes the scan
  d$o <- ordered(sample(1:3, 30L, TRUE))
  f <- fr(o ~ gp(x), data = d, family = cumulative())
  expect_null(frmtmb:::gp_working_y(f, lp1(f)))
  expect_null(frmtmb:::gp_start_working(f))
})

test_that("fresh curves reach the maximum, on a bound and on the review's worst", {
  # the best 3 cells overall share a basin: 14 of 192 short, 13 of them silent,
  # worst 40.5 at sigma 0.14 against 0.0011
  fns <- list(
    saw = function(x) (x %% 2) - 1 + 0.2 * x,
    abs3 = function(x) abs(sin(x)) * 3 - 0.5 * x,
    dbl = function(x) exp(-(x - 3)^2) + 2 * exp(-(x - 7)^2 / 0.1),
    gauss4 = function(x) sin(7 * x) * exp(-x / 5) + 0.3 * x,
    poly = function(x) (x - 2) * (x - 5) * (x - 8) / 10 + sin(4 * x),
    sq = function(x) sign(sin(1.5 * x)) + 0.1 * x,
    tan1 = function(x) atan(5 * (x - 5)) + cos(9 * x) / 2,
    big = function(x) 50 + 10 * sin(2.5 * x) + 3 * x)
  gaps <- numeric(0)
  for (nm in names(fns)) for (s in c(0.1, 1e-3)) for (n in c(30L, 70L)) {
    set.seed(55011 + n)
    x <- sort(round(stats::runif(n, 0, 10), 4))
    y <- fns[[nm]](x) + s * stats::rnorm(n)
    fit <- gs_fit(data.frame(x = x, y = y))
    gaps <- c(gaps, as.numeric(logLik(fit)) + gs_max(y, x, fit))
  }
  expect_lte(mean(gaps < -1e-3), 1 / 16)
  # the review's worst row by name
  set.seed(55011 + 150L)
  x <- sort(round(stats::runif(150L, 0, 10), 4))
  y <- fns$sq(x) + 1e-3 * stats::rnorm(150L)
  fit <- gs_fit(data.frame(x = x, y = y))
  expect_lt(exp(fit$estimates$betad[[1L]]), 10 * 1e-3)
})

gs_false <- function(mock_shift = NULL) {
  # the rows of noise 3e-4, where nlminb from the maximum stops at false
  # convergence on most: returns, per row, its code and the warnings it raised
  lapply(c(1L, 3L, 5L), function(seed) {
    set.seed(seed)
    n <- 100L
    x <- round(stats::runif(n, 0, 6), 3)
    y <- sin(2 * x) + 3e-4 * stats::rnorm(n)
    w <- character()
    fit <- withCallingHandlers(
      frm(bf(y ~ gp(x)), data = data.frame(x = x, y = y)),
      warning = function(cn) {
        w <<- c(w, conditionMessage(cn))
        invokeRestart("muffleWarning")
      })
    list(conv = fit$opt$convergence, w = w, fit = fit)
  })
}

test_that("false convergence is quiet only where the marginal proves the stop", {
  # 13 of 15 fits at 3e-4 warned, all at the maximum, because the warm start starts
  # there. A Newton step alone does not see a ridge, so its clause is not enough
  rows <- gs_false()
  said <- function(r) any(grepl("did not report convergence", r$w, fixed = TRUE))
  stopped <- Filter(function(r) r$conv != 0L, rows)
  expect_gt(length(stopped), 0L)
  for (r in stopped) {
    expect_false(said(r))
    # the code stays
    expect_match(r$fit$opt$message, "false convergence", fixed = TRUE)
  }
  # the same stop, against a marginal that read another number, still warns
  real <- frmtmb:::gp_marginal_max
  local_mocked_bindings(gp_marginal_max = function(frame, keep = 3L, distinct = 6L) {
    m <- real(frame, keep, distinct)
    m$value <- m$value + 1
    m
  }, .package = "frmtmb")
  rows2 <- gs_false()
  stopped2 <- Filter(function(r) r$conv != 0L, rows2)
  expect_gt(length(stopped2), 0L)
  for (r in stopped2) {
    expect_true(said(r))
    expect_match(paste(r$w, collapse = " "), "Newton step", fixed = TRUE)
  }
})

test_that("no warm start, no proof: a theta from the caller still warns", {
  set.seed(3)
  n <- 100L
  x <- round(stats::runif(n, 0, 6), 3)
  y <- sin(2 * x) + 3e-4 * stats::rnorm(n)
  d <- data.frame(x = x, y = y)
  f0 <- suppressWarnings(frm(bf(y ~ gp(x)), data = d))
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(y ~ gp(x)), data = d, start = list(theta = f0$estimates$theta,
                                              betad = f0$estimates$betad,
                                              beta = f0$estimates$beta)),
    warning = function(cn) {
      w <<- c(w, conditionMessage(cn))
      invokeRestart("muffleWarning")
    })
  expect_null(fit$obj$env$.frmtmb_gp_marginal)
  if (fit$opt$convergence != 0L) {
    expect_true(any(grepl("did not report convergence", w, fixed = TRUE)))
  }
})

test_that("a start from the caller is the caller's", {
  set.seed(2)
  x <- sort(stats::runif(30L, 0, 6))
  d <- data.frame(x = x, y = sin(2 * x) + 0.1 * stats::rnorm(30L))
  local_mocked_bindings(gp_marginal_max = function(...) stop("read it"),
                        .package = "frmtmb")
  expect_error(gs_fit(d), "read it")
  expect_no_error(gs_fit(d, start = list(theta = c(log(0.8), log(1.1)))))
  # a beta replaces only itself: the marginal is still read
  expect_error(gs_fit(d, start = list(beta = 0.2)), "read it")
})

test_that("the marginal reads the tape's objective at one point", {
  # every spelling it serves, so a warm start is a start on this model
  set.seed(6101)
  n <- 60L
  d <- data.frame(x = stats::runif(n, 0, 8), z = stats::runif(n, 0, 8),
                  g = factor(sample(letters[1:5], n, TRUE)),
                  f = factor(sample(c("a", "b"), n, TRUE)),
                  w = stats::runif(n, 0.5, 2), s = 0.05,
                  xr = rep(sort(stats::runif(15L, 0, 8)), each = 4L))
  d$y <- sin(2 * d$x) + 0.5 * sin(d$z) + sin(2 * d$xr) + 0.1 * stats::rnorm(n) +
    0.2 * stats::rnorm(5L)[as.integer(d$g)]
  specs <- list(y ~ gp(x), y ~ z + gp(x), y ~ gp(x, by = f), y ~ gp(x, by = w),
                y ~ gp(x) + gp(z), y ~ z + gp(x) + (1 | g), y ~ gp(xr),
                y ~ gp(x, z, iso = FALSE), y ~ gp(x) + s(z),
                y | se(s) ~ gp(x))
  for (sp in specs) {
    fit <- allow_warnings(frm(sp, data = d),
                          c("Optimizer did not report", "gradient"))
    s <- frmtmb:::gp_marginal_setup(fit$frame)
    expect_false(is.null(s))
    p <- list(beta = fit$estimates$beta, theta = fit$estimates$theta)
    if (s$has_sigma) p$betad <- fit$estimates$betad
    v <- frmtmb:::gp_marginal_nll(s)(p)
    expect_lt(abs(v - fit$opt$objective), 1e-7 * abs(fit$opt$objective))
  }
})

test_that("the marginal is not used where it is not the model", {
  set.seed(6102)
  n <- 40L
  d <- data.frame(x = stats::runif(n, 0, 6), z = stats::runif(n))
  d$y <- sin(2 * d$x) + 0.2 * stats::rnorm(n)
  d$k <- stats::rpois(n, exp(1 + sin(2 * d$x)))
  gate <- function(...) {
    f <- allow_warnings(frm(...), c("Optimizer did not report", "gradient",
                                    "Standard errors are not"))
    frmtmb:::gp_marginal_setup(f$frame)
  }
  expect_null(gate(bf(k ~ gp(x)), data = d, family = poisson()))
  expect_null(gate(bf(y ~ gp(x), sigma ~ z), data = d))
  expect_null(gate(bf(y ~ gp(x)), data = d, family = student()))
  expect_false(is.null(gate(bf(y ~ gp(x)), data = d)))
})

test_that("a fit out of the marginal's reach reads the scan in y's units", {
  # sigma ~ z: no marginal, so the scan decides, and at unit sd it read every
  # length scale alike at y * 1000 and stopped with an empty field
  set.seed(4141)
  n <- 40L
  x <- round(stats::runif(n, 0, 6), 3)
  e <- stats::rnorm(n)
  z <- stats::runif(n)
  d <- data.frame(x = x, z = z, y = 1000 * (sin(4 * x) + 0.01 * e))
  fit <- allow_warnings(frm(bf(y ~ gp(x), sigma ~ z), data = d),
                        c("Optimizer did not report", "gradient"))
  sdgp <- exp(fit$estimates$theta[[1L]])
  expect_gt(sdgp, 100)
})

test_that("the warning at the nugget names it, and only there", {
  set.seed(8901)
  n <- 60L
  x <- round(stats::runif(n, 0, 6), 3)
  e <- stats::rnorm(n)
  warned <- function(s) {
    seen <- NULL
    withCallingHandlers(
      frm(bf(y ~ gp(x)), data = data.frame(x = x, y = sin(2 * x) + s * e)),
      warning = function(w) {
        if (inherits(w, "frmtmb_gp_nugget")) seen <<- w
        invokeRestart("muffleWarning")
      })
    seen
  }
  w <- warned(1e-7)
  expect_s3_class(w, c("frmtmb_gp_nugget", "frmtmb_warning"))
  expect_match(conditionMessage(w), "1e-6 * sdgp", fixed = TRUE)
  expect_null(warned(1e-3))
  expect_null(warned(1e-4))
  # with sigma known and positions repeating, sdgp is what moves
  u <- sort(stats::runif(20L, 0, 6))
  xr <- rep(u, each = 5L)
  dr <- data.frame(x = xr, s = 1e-6)
  dr$y <- sin(2 * xr) + 1e-6 * stats::rnorm(100L)
  seen <- NULL
  withCallingHandlers(
    frm(bf(y | se(s) ~ gp(x)), data = dr),
    warning = function(w) {
      if (inherits(w, "frmtmb_gp_nugget")) seen <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(seen, "frmtmb_gp_nugget")
  # a large mean's information is far from it
  dp <- data.frame(x = x)
  dp$y <- stats::rpois(n, exp(16 + sin(2 * x)))
  fp <- allow_warnings(frm(bf(y ~ gp(x)), data = dp, family = poisson()),
                       c("Optimizer did not report", "gradient"))
  expect_null(frmtmb::gp_nugget_message(fp))
  # no gp, no question
  f0 <- frm(bf(y ~ x), data = data.frame(x = x, y = e))
  expect_null(frmtmb::gp_nugget_message(f0))
})

test_that("standard errors are read in the response's units", {
  # (Intercept) read NaN, "flat", at y * 1000, and every fitted standard error
  # with it, on a raw difference of 1e-3 against sigma 0.1
  set.seed(14)
  x <- sort(stats::runif(40L, 0, 8))
  y <- sin(x) + 1e-4 * stats::rnorm(40L)
  se_of <- function(c0) {
    f <- gs_fit(data.frame(x = x, y = y * c0))
    sqrt(diag(as.matrix(vcov(f, full = TRUE))))[1L] / c0
  }
  s1 <- se_of(1)
  for (c0 in c(300, 1000, 1e4, 1e6)) {
    expect_lt(abs(se_of(c0) / s1 - 1), 1e-2)
  }
  # and at unit scale outside a gp, a large unit reads the small one's
  set.seed(15)
  g <- gl(12, 8)
  z <- stats::rnorm(96L)
  v <- 2 + 0.5 * z + stats::rnorm(12L, 0, 0.7)[g] + 0.5 * stats::rnorm(96L)
  lme_se <- function(c0) {
    f <- frm(bf(v ~ z + (1 | g)), data = data.frame(v = v * c0, z = z, g = g))
    sqrt(diag(as.matrix(vcov(f))))[1:2] / c0
  }
  a <- lme_se(1)
  for (c0 in c(1e3, 1e6)) {
    expect_lt(max(abs(lme_se(c0) / a - 1)), 1e-5)
  }
})

test_that("a sigma the field took stays flat at the row check", {
  # at noise 1 the field can take the noise, and sigma runs to e^-10 with next to
  # no information: the row check read that as curved and gave it 6341
  set.seed(60002)
  n <- 50L
  x <- sort(round(stats::runif(n, 0, 8), 4))
  z <- stats::rnorm(n)
  stats::runif(0L)
  g <- factor(sample(letters[1:5], n, TRUE))
  y <- sin(2 * x) + 0.3 * z + stats::rnorm(n)
  w <- character()
  m <- character()
  fit <- withCallingHandlers(
    frm(y ~ z + gp(x), data = data.frame(x, z, g, y)),
    warning = function(cn) {
      w <<- c(w, conditionMessage(cn))
      invokeRestart("muffleWarning")
    },
    message = function(cn) {
      m <<- c(m, conditionMessage(cn))
      invokeRestart("muffleMessage")
    })
  expect_lt(exp(fit$estimates$betad[[1L]]), 1e-2)
  se <- sqrt(diag(as.matrix(vcov(fit, full = TRUE))))
  expect_true(is.nan(se[["sigma_(Intercept)"]]))
  # either verdict says so: "flat", or sigma at the end of its range,
  # where the field interpolates the rows (se_edge_limit())
  expect_true(any(grepl("Standard errors are not available", w,
                        fixed = TRUE)) ||
                any(grepl("sigma_(Intercept) is at the end of its range", m,
                          fixed = TRUE)))
})

test_that("a row under its noise is flat only if it reads no curvature", {
  # a parameter nothing depends on stays flat through the second reading
  set.seed(16)
  n <- 30L
  d <- data.frame(x = sort(stats::runif(n, 0, 6)))
  d$y <- sin(d$x) + 0.1 * stats::rnorm(n)
  fit <- gs_fit(d)
  p <- fit$opt$par
  u <- fit$par_units %||% rep(1, length(p))
  expect_true(frmtmb:::se_row_curved(1L, fit, p, u))
})

test_that("retaping a fit outside any scope keeps the tape's settings", {
  # a worker of parallel chains retapes on its own: 11 GB at 300 positions
  set.seed(17)
  d <- data.frame(x = sort(stats::runif(60L, 0, 6)))
  d$y <- sin(2 * d$x) + 0.1 * stats::rnorm(60L)
  fit <- gs_fit(d)
  env <- fit$obj$env
  v0 <- fit$obj$fn(fit$opt$par)
  g0 <- fit$obj$gr(fit$opt$par)
  flag <- function() TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  before <- flag()
  env$retape(set.defaults = FALSE)
  expect_identical(flag(), before)
  expect_equal(fit$obj$fn(fit$opt$par), v0, tolerance = 1e-10)
  expect_equal(as.numeric(fit$obj$gr(fit$opt$par)), as.numeric(g0),
               tolerance = 1e-6)
})

test_that("the outer Hessian of an exact gp() reads the marginal's gradient", {
  # each gradient of the z tape is an inner solve over a dense Hessian: 11.3 s
  # of a 25 s fit at 500 positions (dev/fitedge-findings.md, item 3)
  hess <- function(fit) {
    g0 <- fit$obj$gr
    n_gr <- 0L
    fit$obj$gr <- function(...) {
      n_gr <<- n_gr + 1L
      g0(...)
    }
    env <- fit$obj$env
    rm(list = intersect(".frmtmb_outer_h", ls(env, all.names = TRUE)),
       envir = env)
    h <- frmtmb:::fit_outer_hessian(fit)
    # the z tape's own differences at the same steps, for reference
    pn <- names(fit$opt$par)
    u <- fit$par_units
    ref <- if (is.null(u) || all(u == 1)) {
      frmtmb:::fd_hessian(h$at, g0,
                          eps = frmtmb:::outer_fd_steps(fit$frame, pn))$H
    } else {
      frmtmb:::fd_hessian(h$at / u, function(q) g0(q * u) * u,
                          eps = frmtmb:::outer_fd_steps(fit$frame, pn, u))$H /
        outer(u, u)
    }
    list(n = n_gr, H = h$H, ref = ref)
  }
  ok <- c("Optimizer did not report", "gradient", "Standard errors are not")
  fit <- allow_warnings(frm(y ~ z + gp(x), data = gs_edge), ok)
  h <- hess(fit)
  expect_identical(h$n, 0L)
  expect_lt(max(abs(h$H - h$ref)), 1e-6 * max(abs(h$ref)))
  # the guard's absent case: a prior is another objective, and so is a family
  # with no closed form; both keep the z tape
  fp <- allow_warnings(
    frm(y ~ z + gp(x), data = gs_edge,
        prior = set_prior("normal(0, 1)", class = "b")), ok)
  expect_null(frmtmb:::gp_marginal_outer_gr(fp, fp$opt$par))
  expect_gt(hess(fp)$n, 0L)
  d <- gs_edge
  d$k <- stats::rpois(nrow(d), exp(1 + 0.3 * sin(d$x)))
  fk <- allow_warnings(frm(k ~ gp(x), data = d, family = poisson()), ok)
  expect_null(frmtmb:::gp_marginal_outer_gr(fk, fk$opt$par))
  expect_gt(hess(fk)$n, 0L)
  # and a covariance too ill conditioned for the marginal's rounding: at noise
  # 1e-3 on 2 sin(3x) + x its standard errors were 35 percent off an exact
  # Hessian, and the z tape's 3e-5
  set.seed(1)
  x <- round(stats::runif(80L, 0, 10), 3)
  f1 <- gs_fit(data.frame(x = x, y = 2 * sin(3 * x) + x +
                            1e-3 * stats::rnorm(80L)))
  expect_null(frmtmb:::gp_marginal_outer_gr(f1, f1$opt$par))
  expect_gt(hess(f1)$n, 0L)
})

# rows of the body sweeps: dev/fitedge-gpbody-sweep.R (`fe`) and
# dev/review-fitedge-designs.R (`rv`), same seeds and draws
gs_fe_row <- function(seed, n, curve, body) {
  set.seed(18000 + 100 * seed + n + match(curve, c("sin", "wig", "exp")))
  x <- sort(stats::runif(n, 0, 10))
  z <- stats::rnorm(n)
  w <- stats::runif(n, 1, 2)
  f <- switch(curve, sin = sin(x), wig = 2 * sin(3 * x) + x)
  y <- switch(body,
              A = f + 0.3 * z + 0.1 * stats::rnorm(n),
              C = exp(0.5 + 0.5 * f / max(abs(f))) + 0.3 * z +
                0.1 * stats::rnorm(n))
  data.frame(x, z, w, y)
}
gs_rv_row <- function(design, seed, n, ci) {
  di <- c(E = 1L, F = 2L, K = 7L, L = 8L)[[design]]
  set.seed(52000 + 100 * seed + n + 10 * di + ci)
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  g <- factor(sample(1:8, n, TRUE))
  fby <- factor(sample(c("p", "q"), n, TRUE))
  f <- switch(ci, sin(x), 2 * sin(3 * x) + x, exp(0.3 * x))
  fs <- f / max(abs(f))
  ge <- stats::rnorm(8, 0, 0.5)[g]
  y <- switch(design,
              E = (2 + fs) * exp(0.4 * z) + 0.1 * stats::rnorm(n),
              F = (3 + fs) / (1 + exp(-1.5 * (z - 0.2))) +
                0.1 * stats::rnorm(n),
              K = stats::rgamma(n, 5, 5 / exp(0.5 + 0.5 * fs + 0.3 * z)),
              L = 2 * exp((0.3 + 0.2 * fs) * z) + 0.1 * stats::rnorm(n))
  data.frame(x, z, g, fby, y)
}
gs_ll <- function(...) {
  as.numeric(logLik(allow_warnings(frm(...), c("Optimizer did not report",
                                               "gradient",
                                               "Standard errors are not"))))
}
# 0.79.1's start for a body's field: the body start switched off
gs_ll_off <- function(...) {
  local_mocked_bindings(gp_body_start = function(...) NULL, .package = "frmtmb")
  gs_ll(...)
}

test_that("a field in a body never ends below the start it replaced", {
  # Gamma(log) `a + b * z` read 3.30 nats below 0.79.1 with the body start
  # alone; the scan's start now descends too and the better fit is kept
  # (dev/reviews/2026-11-06-fitedge.md, B1)
  dK <- gs_rv_row("K", 1L, 60L, 3L)
  fK <- bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE)
  l0 <- gs_ll_off(fK, data = dK, family = Gamma(link = "log"))
  l1 <- gs_ll(fK, data = dK, family = Gamma(link = "log"))
  expect_gt(l1, l0 - 1e-6 * abs(l0))
  # a caller's start is the caller's: from c = 8 the body start read 125 nats
  # below 0.79.1
  dF <- gs_rv_row("F", 1L, 60L, 2L)
  fF <- bf(y ~ a / (1 + exp(-b * (z - c))), a ~ gp(x), b ~ 1, c ~ 1,
           nl = TRUE)
  st <- list(beta = c("c_(Intercept)" = 8))
  l0 <- gs_ll_off(fF, data = dF, start = st)
  l1 <- gs_ll(fF, data = dF, start = st)
  expect_lt(abs(l1 - l0), 1e-8 * abs(l0))
})

test_that("the body start reads every coefficient, in one step", {
  # from y alone, `a + b * z` stopped 0.44 nats below `y ~ z + gp(x)`, which
  # is the same model, and `exp(a) + b * z` 17 nats below the best of 12
  # starts; a second step left `a * exp(b * z)` 44 nats below it
  # (dev/fitedge-findings.md, item 1 and fix round)
  dA <- gs_fe_row(2L, 60L, "sin", "A")
  lin <- gs_ll(y ~ z + gp(x), data = dA)
  body <- gs_ll(bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE), data = dA)
  expect_lt(abs(body - lin), 1e-6 * abs(lin))
  dC <- gs_fe_row(2L, 60L, "wig", "C")
  fC <- bf(y ~ exp(a) + b * z, a ~ gp(x), b ~ 1, nl = TRUE)
  best <- -Inf
  for (rho in c(0.15, 0.5, 1.5, 5)) for (sd0 in c(0.3, 1, 3)) {
    v <- tryCatch(gs_ll(fC, data = dC,
                        start = list(theta = c(log(sd0), log(rho)))),
                  error = function(e) -Inf)
    best <- max(best, v)
  }
  expect_gt(gs_ll(fC, data = dC), best - 1e-5 * abs(best))
  dE <- gs_rv_row("E", 1L, 120L, 2L)
  fE <- bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE)
  best <- -Inf
  for (rho in c(0.15, 0.5, 1.5, 5)) for (sd0 in c(0.3, 1, 3)) {
    v <- tryCatch(gs_ll(fE, data = dE,
                        start = list(theta = c(log(sd0), log(rho)))),
                  error = function(e) -Inf)
    best <- max(best, v)
  }
  expect_gt(gs_ll(fE, data = dE), best - 1e-5 * abs(best))
})

test_that("a multiplicative body is linearized again at the fit's optimum", {
  skip_on_cran()
  # at the template's zeros `a * exp(b * z)` scales a's field by 1, and
  # `b ~ gp(x)` by 0: on 0.81.0 these read 36.50 and 5.24 nats below one
  # theta start, conv 0 and no warning (dev/gpstop-findings.md, item 1)
  dE <- gs_rv_row("E", 2L, 120L, 2L)
  fE <- bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE)
  ref <- gs_ll(fE, data = dE, start = list(theta = c(log(1.5), log(0.3))))
  expect_gt(gs_ll(fE, data = dE), ref - 1e-5 * abs(ref))
  dL <- gs_rv_row("L", 2L, 120L, 2L)
  fL <- bf(y ~ a * exp(b * z), a ~ 1, b ~ gp(x), nl = TRUE)
  ref <- gs_ll(fL, data = dL, start = list(theta = c(log(0.3), log(0.1))))
  expect_gt(gs_ll(fL, data = dL), ref - 1e-5 * abs(ref))
})

# rows of the review's sweep (dev/reviews/gpstop-review-sweep2.R): a
# logistic body with the field in its asymptote
gs_logis_row <- function(fam, n, seed) {
  set.seed(81100 + 1000 * match(fam, c("gaussian", "poisson", "bernoulli",
                                       "student")) +
             100 * seed + n + 30)
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  w <- round(stats::runif(n, 0.2, 5), 3)
  f <- if (seed %% 2) sin(x) else 2 * sin(3 * x) + x
  eta <- (3 + f / max(abs(f))) / (1 + exp(-1.5 * z))
  y <- switch(fam,
              student = eta + 0.1 * stats::rt(n, 3),
              bernoulli = stats::rbinom(n, 1, stats::plogis(
                1.5 * (eta - stats::median(eta)) / stats::sd(eta))))
  data.frame(x, z, w, y)
}
# the smallest eigenvalue of the fit's inner Hessian at its optimum
gs_inner_ev <- function(fit) {
  env <- fit$obj$env
  fit$obj$fn(fit$opt$par)
  H <- as.matrix(env$spHess(env$last.par, random = TRUE))
  min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
}

test_that("a second reading of the body keeps only a regular maximum", {
  skip_on_cran()
  # the descent from it rose to a student() spike (code 1, gradient 1.2e13,
  # inner eigenvalue 1.97e-7) and to a bernoulli field of length scale
  # e^-323 at code 0 with no warning; importance sampling read both below
  # 0.81.0's point (dev/reviews/2026-11-11-gpstop.md, MAJOR 1 and 2)
  fl <- bf(y ~ a / (1 + exp(-b * z)), a ~ gp(x), b ~ 1, nl = TRUE)
  ok <- c("Optimizer did not report", "gradient", "Standard errors are not")
  for (fam in c("student", "bernoulli")) {
    d <- gs_logis_row(fam, 40L, 2L)
    family <- if (fam == "student") student() else bernoulli()
    fit <- allow_warnings(frm(fl, data = d, family = family), ok)
    ref <- local({
      local_mocked_bindings(gp_body_relin_applies = function(...) FALSE,
                            .package = "frmtmb")
      allow_warnings(frm(fl, data = d, family = family), ok)
    })
    ll <- as.numeric(logLik(fit))
    l0 <- as.numeric(logLik(ref))
    # never below the fit without the second reading ...
    expect_gt(ll, l0 - 1e-8 * abs(l0))
    # ... and, where above it, at a regular maximum of the field
    if (ll > l0 + 1e-8 * abs(l0)) {
      expect_identical(fit$opt$convergence, 0L)
      gap <- min(diff(sort(unique(d$x))))
      expect_gt(exp(fit$estimates$theta[2L]), gap / 10)
      expect_gt(gs_inner_ev(fit), 0.01 * gs_inner_ev(ref))
    }
  }
})

test_that("a white-noise field is refused at any length scale below the gaps", {
  skip_on_cran()
  # the bernoulli artifact polished from a length scale near the smallest
  # gap is a stationary point at code 0 with a healthy inner Hessian, at
  # 0.67 of the smallest gap, and at 52 times it beside one near tie; a
  # rule on the smallest gap kept both (review, short check S1)
  fl <- bf(y ~ a / (1 + exp(-b * z)), a ~ gp(x), b ~ 1, nl = TRUE)
  ok <- c("Optimizer did not report", "gradient", "Standard errors are not")
  regular <- frmtmb:::gp_relin_regular
  for (v in c("plain", "tie")) {
    d <- gs_logis_row("bernoulli", 40L, 2L)
    if (v == "tie") d$x[20] <- d$x[19] + 1e-4
    fit <- allow_warnings(frm(fl, data = d, family = bernoulli()), ok)
    obj <- fit$obj
    # the absent case: the fit's own regular maximum passes
    expect_true(regular(obj, fit$frame, fit$opt, 1e-3))
    mg <- min(diff(sort(unique(d$x))))
    pn <- names(obj$par)
    p <- obj$par
    p[pn == "beta"] <- c(-11.046, -135.46)
    p[pn == "theta"] <- c(3.5988, log((if (v == "tie") 20 else 5) * mg / 10))
    o <- suppressWarnings(stats::nlminb(p, obj$fn, obj$gr,
                                        control = list(iter.max = 200)))
    o$convergence <- as.integer(o$convergence)
    expect_false(regular(obj, fit$frame, o, 1e-3))
    obj$fn(fit$opt$par)
  }
})

test_that("the body is linearized again only for its field, with no start", {
  n_relin <- 0L
  real <- frmtmb:::gp_body_relin
  local_mocked_bindings(gp_body_relin = function(...) {
    n_relin <<- n_relin + 1L
    real(...)
  }, .package = "frmtmb")
  ok <- c("Optimizer did not report", "gradient", "Standard errors are not")
  dE <- gs_rv_row("E", 1L, 60L, 2L)
  fE <- bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE)
  # the absent case: a field in a body, no start of the caller's
  fit <- allow_warnings(frm(fE, data = dE), ok)
  expect_gt(n_relin, 0L)
  # update() refits from the cold start, as frm() does
  n_relin <- 0L
  allow_warnings(update(fit, data = dE), ok)
  expect_gt(n_relin, 0L)
  # refit() is handed the fit's estimates, so it has a start
  n_relin <- 0L
  allow_warnings(refit(fit, newresp = dE$y * 1.01), ok)
  expect_identical(n_relin, 0L)
  n_relin <- 0L
  allow_warnings(frm(fE, data = dE,
                     start = list(beta = c("b_(Intercept)" = 0.3))), ok)
  expect_identical(n_relin, 0L)
  # a field in the mean: the working response reads no body
  set.seed(9)
  dP <- data.frame(x = dE$x, y = stats::rpois(nrow(dE), exp(1 + sin(dE$x))))
  allow_warnings(frm(y ~ gp(x), data = dP, family = poisson()), ok)
  expect_identical(n_relin, 0L)
})

test_that("a caller's start turns the body start off on every refit path", {
  # frm_allfit(), the autoscale pre-fit and refit(newresp =) forwarded a
  # caller's `start` with the body start still on
  # (dev/reviews/2026-11-06-fitedge.md, short check, m1)
  dF <- gs_rv_row("F", 1L, 60L, 2L)
  fF <- bf(y ~ a / (1 + exp(-b * (z - c))), a ~ gp(x), b ~ 1, c ~ 1,
           nl = TRUE)
  n_body <- 0L
  real <- frmtmb:::gp_body_start
  local_mocked_bindings(gp_body_start = function(...) {
    n_body <<- n_body + 1L
    real(...)
  }, .package = "frmtmb")
  ok <- c("Optimizer did not report", "gradient", "Standard errors are not")
  # the absent case: no start of the caller's, and the body start runs
  allow_warnings(frm(fF, data = dF), ok)
  expect_gt(n_body, 0L)
  n_body <- 0L
  fit <- allow_warnings(frm(fF, data = dF,
                            start = list(beta = c("c_(Intercept)" = 8))), ok)
  expect_identical(n_body, 0L)
  af <- frm_allfit(fit, optimizers = "nlminb", start_from_mle = FALSE)
  expect_identical(n_body, 0L)
  expect_s3_class(af, "frmtmb_allfit")
})
