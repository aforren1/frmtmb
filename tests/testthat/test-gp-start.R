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
    list(y ~ mo(o) + gp(x))
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
  local_mocked_bindings(gp_marginal_max = function(frame, keep = 3L, distinct = 4L) {
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
  fit <- withCallingHandlers(
    frm(y ~ z + gp(x), data = data.frame(x, z, g, y)),
    warning = function(cn) {
      w <<- c(w, conditionMessage(cn))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  expect_lt(exp(fit$estimates$betad[[1L]]), 1e-2)
  se <- sqrt(diag(as.matrix(vcov(fit, full = TRUE))))
  expect_true(is.nan(se[["sigma_(Intercept)"]]))
  expect_true(any(grepl("Standard errors are not available", w, fixed = TRUE)))
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
