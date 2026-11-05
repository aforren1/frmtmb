# The exact gp() is non-centered at a nugget relative to sdgp: b is z, a
# standard normal, and the field is sd L(lscale) z. On a gaussian response the
# field integrates in closed form, so every number here has a reference written
# from the formula, with no frmtmb code in it. The designs are the ones that
# failed on the centered field at a nugget of 1e-6
# (dev/ncgp-findings.md has each failure as measured on 0.78.0).

gpnc_marginal <- function(p, y, x, se = NULL, nug = 1e-12) {
  u <- sort(unique(x))
  J <- match(x, u)
  s2 <- exp(2 * p[2])
  K <- s2 * exp(-outer(u, u, "-")^2 / (2 * exp(2 * p[3])))
  diag(K) <- diag(K) * (1 + nug)
  V <- K[J, J, drop = FALSE]
  diag(V) <- diag(V) + if (is.null(se)) exp(2 * p[4]) else se^2
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(1e12)
  r <- backsolve(R, y - p[1], transpose = TRUE)
  0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * length(y) * log(2 * pi)
}

# the marginal's maximum from a grid of starts and from frm()'s own point, so
# that frm() cannot pass by a reference stopping short of it
gpnc_max <- function(y, x, fit = NULL, se = NULL) {
  ry <- stats::sd(y)
  rx <- diff(range(x))
  starts <- list()
  for (r in rx / c(16, 4, 1)) {
    for (ns in ry * c(1e-1, 1e-3)) {
      starts[[length(starts) + 1L]] <- c(mean(y), log(ry), log(r),
                                         if (is.null(se)) log(ns))
    }
  }
  if (!is.null(fit)) {
    bk <- fit$frame$re_blocks[[1L]]
    th <- fit$estimates$theta[bk$theta_idx]
    starts <- c(starts, list(c(fit$estimates$beta[[1L]], th,
                               if (is.null(se)) log(sigma(fit)[1L]))))
  }
  best <- NULL
  for (st in starts) {
    sc <- c(ry, rep(1, length(st) - 1L))
    o <- tryCatch(
      stats::nlminb(st / sc, function(q) gpnc_marginal(q * sc, y, x, se),
                    control = list(iter.max = 500, eval.max = 800,
                                   rel.tol = 1e-12)),
      error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$objective < best$objective)) {
      o$par <- o$par * sc
      best <- o
    }
  }
  best
}

gpnc_read <- function(fit) {
  bk <- fit$frame$re_blocks[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  list(ll = -fit$opt$objective, sdgp = exp(th[[1L]]), rho = exp(th[[2L]]))
}

gpnc_low <- function(s, n = 60L) {
  set.seed(20261102)
  x <- round(stats::runif(n, 0, 6), 3)
  e0 <- stats::rnorm(n)
  data.frame(x = x, y = sin(2 * x) + s * e0)
}

gpnc_fit <- function(d, ...) {
  allow_warnings(frm(bf(y ~ gp(x)), data = d, family = gaussian(), ...),
                 c("Optimizer did not report", "gradient"))
}

test_that("frm() reaches the closed form's maximum at low noise", {
  # On 0.78.0: sigma 1.27e-7 at a true 1e-3, and the fit 1.77 nats short
  for (s in c(1e-2, 3e-3, 1e-3, 3e-4)) {
    d <- gpnc_low(s)
    fit <- gpnc_fit(d)
    cf <- gpnc_max(d$y, d$x, fit)
    r <- gpnc_read(fit)
    expect_lt(abs(r$ll + cf$objective), 1e-5 * abs(cf$objective))
    expect_lt(abs(log(sigma(fit)[1L]) - cf$par[4L]), 0.01)
    expect_lt(abs(log(r$sdgp) - cf$par[2L]), 0.02)
  }
})

test_that("frm() reaches it at ordinary noise on a field of large range", {
  # On 0.78.0, x^2 with sigma 0.3: sigma 6.1e-5, against the marginal's 0.348
  set.seed(5701)
  x <- sort(stats::runif(100, 0, 10))
  d <- data.frame(x = x, y = x^2 + 0.3 * stats::rnorm(100))
  fit <- gpnc_fit(d)
  cf <- gpnc_max(d$y, d$x, fit)
  expect_lt(abs(-fit$opt$objective + cf$objective),
            1e-5 * abs(cf$objective))
  expect_lt(abs(sigma(fit)[1L] / exp(cf$par[4L]) - 1), 0.01)
})

test_that("frm() finds a rough field's maximum, not the empty one", {
  # On 0.78.0 and without the start scan: ll -116.11, sigma 0.77, the field
  # near zero, against the marginal's -72.37
  set.seed(5702)
  x <- sort(stats::runif(100, 0, 10))
  e0 <- stats::rnorm(100)
  d <- data.frame(x = x, y = sin(8 * x) + 0.3 * e0)
  fit <- gpnc_fit(d)
  cf <- gpnc_max(d$y, d$x, fit)
  expect_lt(abs(-fit$opt$objective + cf$objective),
            1e-5 * abs(cf$objective))
  expect_lt(sigma(fit)[1L], 0.4)
})

test_that("a known se() meets no floor on the field", {
  # On 0.78.0: 286.5 nats against the marginal's 318.3 at se 3e-4
  d <- gpnc_low(3e-4)
  d$s <- 3e-4
  fit <- allow_warnings(frm(bf(y | se(s) ~ gp(x)), data = d,
                            family = gaussian()),
                        c("Optimizer did not report", "gradient"))
  cf <- gpnc_max(d$y, d$x, fit, se = d$s)
  expect_lt(abs(-fit$opt$objective + cf$objective),
            1e-5 * abs(cf$objective))
})

test_that("replicated positions leave sdgp at the marginal's", {
  # On 0.78.0, 20 positions x 30 at noise 3e-4: sdgp 0.85 against 2.11
  set.seed(131)
  u <- sort(stats::runif(20, 0, 6))
  x <- rep(u, each = 30L)
  e0 <- stats::rnorm(length(x))
  d <- data.frame(x = x, y = sin(2 * x) + 3e-4 * e0)
  fit <- gpnc_fit(d)
  cf <- gpnc_max(d$y, d$x, fit)
  expect_lt(abs(log(gpnc_read(fit)$sdgp) - cf$par[2L]), 0.02)
  expect_lt(abs(-fit$opt$objective + cf$objective),
            1e-5 * abs(cf$objective))
})

test_that("the fit does not depend on y's units or on x's", {
  # An absolute nugget fails this: 0 of 18 fits converged at y times 1000 with
  # 1e-12 absolute, ll -1118 against -377
  d <- gpnc_low(1e-3)
  n <- nrow(d)
  f0 <- gpnc_read(fit0 <- gpnc_fit(d))
  dy <- transform(d, y = y * 1000)
  f1 <- gpnc_read(fit1 <- gpnc_fit(dy))
  expect_lt(abs(f1$ll - (f0$ll - n * log(1000))), 1e-5 * abs(f0$ll))
  expect_lt(abs(log(sigma(fit1)[1L] / sigma(fit0)[1L]) - log(1000)), 0.01)
  expect_lt(abs(log(f1$sdgp / f0$sdgp) - log(1000)), 0.02)
  expect_lt(abs(log(f1$rho / f0$rho)), 0.02)
  dy <- transform(d, y = y / 1000)
  f2 <- gpnc_read(gpnc_fit(dy))
  expect_lt(abs(f2$ll - (f0$ll + n * log(1000))), 1e-5 * abs(f0$ll))
  dx <- transform(d, x = x * 1000)
  f3 <- gpnc_read(gpnc_fit(dx))
  expect_lt(abs(f3$ll - f0$ll), 1e-5 * abs(f0$ll))
  expect_lt(abs(log(f3$rho / f0$rho) - log(1000)), 0.02)
})

test_that("b is z, and the coefficients are the field", {
  d <- gpnc_low(1e-1)
  fit <- gpnc_fit(d)
  bk <- fit$frame$re_blocks[[1L]]
  r <- gpnc_read(fit)
  th <- fit$estimates$theta[bk$theta_idx]
  u <- sort(unique(d$x))
  K <- exp(-outer(u, u, "-")^2 / (2 * r$rho^2))
  diag(K) <- diag(K) + 1e-12
  z <- fit$estimates$b[bk$b_idx]
  cvec <- frmtmb:::coef_b(fit)
  p <- bk$gp_order
  Lh <- t(chol(K[p, p]))[order(p), , drop = FALSE]
  expect_equal(cvec[bk$c_idx], r$sdgp * as.vector(Lh %*% z),
               tolerance = 1e-10)
  # the density on b is the standard normal's, whatever theta is
  nll <- frmtmb:::cs_entry("gp")$nll
  expect_identical(nll(z, th, bk), sum(stats::dnorm(z, 0, 1, log = TRUE)))
  expect_identical(nll(z, th + 1, bk), nll(z, th, bk))
  # a field is not a positionwise copy of b, and every caller that perturbs
  # or batches b one entry at a time has to know
  expect_false(frmtmb:::block_b_positionwise(bk))
  expect_true(frmtmb:::frame_needs_expand(fit$frame))
})

test_that("d field / d theta at z agrees with the tape's", {
  d <- gpnc_low(1e-1)
  d2 <- transform(d, w = rep(c(1, 2), length.out = nrow(d)))
  for (ff in list(bf(y ~ gp(x)), bf(y ~ gp(x, w, iso = FALSE)))) {
    fit <- allow_warnings(frm(ff, data = d2, family = gaussian()),
                          c("Optimizer did not report", "gradient"))
    bk <- Filter(function(b) identical(b$covstruct, "gp"),
                 fit$frame$re_blocks)[[1L]]
    th <- fit$estimates$theta[bk$theta_idx]
    z <- fit$estimates$b[bk$b_idx]
    u <- frmtmb:::gp_field_dtheta(z, th, bk)
    ana <- exp(th[1L]) * frmtmb:::gp_factor(th, bk) %*% u
    tp <- RTMB::MakeTape(function(p) frmtmb:::gp_field(z, p, bk), th)
    expect_lt(max(abs(ana - tp$jacobian(th))),
              1e-4 * max(abs(ana)))
  }
})

test_that("fitted(), predict() and new positions read the conditional mean", {
  d <- gpnc_low(1e-2)
  fit <- gpnc_fit(d)
  r <- gpnc_read(fit)
  sg <- sigma(fit)[1L]
  mu <- fit$estimates$beta[[1L]]
  u <- sort(unique(d$x))
  J <- match(d$x, u)
  k <- function(a, b) r$sdgp^2 * exp(-outer(a, b, "-")^2 / (2 * r$rho^2))
  V <- k(u, u)[J, J] + diag(sg^2, nrow(d))
  al <- solve(V, d$y - mu)
  # at the rows, in their own order (x is not sorted here)
  expect_gt(sum(diff(d$x) < 0), 0L)
  mean_rows <- mu + as.vector(k(d$x, d$x) %*% al)
  fv <- fitted(fit)[, "Estimate"]
  expect_lt(max(abs(fv - mean_rows)), 1e-7 * max(abs(mean_rows)))
  # at unseen positions
  nd <- data.frame(x = c(0.123, 2.345, 5.555))
  kr <- mu + as.vector(k(nd$x, d$x) %*% al)
  expect_lt(max(abs(frm_linpred(fit, newdata = nd) - kr)),
            1e-6 * max(abs(kr)))
  # ranef() names the field at the sorted positions
  re <- ranef(fit)[[1L]]
  expect_lt(max(abs(as.vector(re) - as.vector(k(u, d$x) %*% al))),
            1e-6 * max(abs(as.vector(re))))
})

test_that("the factor runs farthest-first, and its field is the kernel's", {
  d <- gpnc_low(1e-1, n = 30L)
  fit <- gpnc_fit(d)
  bk <- fit$frame$re_blocks[[1L]]
  p <- bk$gp_order
  expect_identical(sort(p), seq_along(bk$b_idx))
  # by hand, on a small one: the medoid, then each farthest from those taken
  u <- c(0, 1, 2, 3.5, 7)
  hand <- list(aux_D2 = list(outer(u, u, "-")^2))
  expect_identical(frmtmb:::gp_farthest_first(hand), c(3L, 5L, 1L, 4L, 2L))
  # rows at one position take their turn at the end
  tied <- list(aux_D2 = list(outer(c(0, 0, 5), c(0, 0, 5), "-")^2))
  expect_identical(frmtmb:::gp_farthest_first(tied)[3L], 2L)
  th <- fit$estimates$theta[bk$theta_idx]
  L <- frmtmb:::gp_factor(th, bk)
  K <- frmtmb:::gp_corr(th, bk)
  expect_lt(max(abs(tcrossprod(L) - K)), 64 * .Machine$double.eps)
  Lo <- frmtmb:::gp_factor_ordered(th, bk)
  expect_true(all(Lo[upper.tri(Lo)] == 0))
  # the reason for the order: over sorted positions the second pivot is already
  # small, and here it is a position far from the first
  expect_gt(diag(Lo)[2L], 2 * diag(t(chol(K)))[2L])
})

test_that("a kriging row's standard error holds at low noise", {
  # the delta method's columns go through the factor, whose rounding is largest
  # here: 0.992 to 1.001 of the closed form on this design (dev/ncgp-perm3.R)
  set.seed(6601)
  n <- 40L
  x <- round(stats::runif(n, 0, 6), 3)
  e0 <- stats::rnorm(n)
  y <- sin(2 * x) + 3e-4 * e0
  d <- data.frame(x = x, y = y)
  fit <- gpnc_fit(d)
  m <- gpnc_max(y, x, fit)
  H <- stats::optimHess(m$par, gpnc_marginal, y = y, x = x)
  Vo <- solve(H)
  nd <- data.frame(x = seq(min(x) + 0.13, max(x) - 0.17, length.out = 7L))
  cond <- function(p) {
    Ko <- exp(2 * p[2]) * exp(-outer(nd$x, x, "-")^2 / (2 * exp(2 * p[3])))
    V <- exp(2 * p[2]) * exp(-outer(x, x, "-")^2 / (2 * exp(2 * p[3])))
    diag(V) <- diag(V) * (1 + 1e-12) + exp(2 * p[4])
    R <- chol(V)
    W <- backsolve(R, t(Ko), transpose = TRUE)
    r <- backsolve(R, y - p[1], transpose = TRUE)
    list(m = p[1] + as.vector(t(W) %*% r), v = exp(2 * p[2]) - colSums(W^2))
  }
  J <- vapply(seq_along(m$par), function(i) {
    e <- replace(numeric(4L), i, 1e-5)
    (cond(m$par + e)$m - cond(m$par - e)$m) / 2e-5
  }, numeric(nrow(nd)))
  ref <- sqrt(pmax(cond(m$par)$v, 0) + rowSums((J %*% Vo) * J))
  se <- suppressWarnings(stats::fitted(fit, newdata = nd))[, 2L]
  expect_lt(max(abs(se / ref - 1)), 0.05)
})

test_that("the delta method reads the same variance as the closed form", {
  d <- gpnc_low(1e-1, n = 40L)
  fit <- gpnc_fit(d)
  r <- gpnc_read(fit)
  y <- d$y
  x <- d$x
  # eta at the rows as a function of the outer parameters on the closed form
  # alone, with its Hessian for the outer parameters' covariance
  p0 <- c(fit$estimates$beta[[1L]], log(r$sdgp), log(r$rho),
          log(sigma(fit)[1L]))
  cond <- function(p, xo) {
    Ko <- exp(2 * p[2]) * exp(-outer(xo, x, "-")^2 / (2 * exp(2 * p[3])))
    V <- exp(2 * p[2]) * exp(-outer(x, x, "-")^2 / (2 * exp(2 * p[3])))
    diag(V) <- diag(V) + exp(2 * p[4])
    p[1] + as.vector(Ko %*% solve(V, y - p[1]))
  }
  H <- stats::optimHess(p0, gpnc_marginal, y = y, x = x)
  Vo <- solve(H)
  grid <- list(rows = x, new = c(0.5, 3.3))
  for (g in names(grid)) {
    xo <- grid[[g]]
    J <- vapply(seq_along(p0), function(i) {
      e <- replace(numeric(4L), i, 1e-5)
      (cond(p0 + e, xo) - cond(p0 - e, xo)) / 2e-5
    }, numeric(length(xo)))
    K <- exp(2 * p0[2]) * exp(-outer(x, x, "-")^2 / (2 * exp(2 * p0[3])))
    Ko <- exp(2 * p0[2]) * exp(-outer(xo, x, "-")^2 / (2 * exp(2 * p0[3])))
    V <- K
    diag(V) <- diag(V) + exp(2 * p0[4])
    cv <- exp(2 * p0[2]) - rowSums((Ko %*% solve(V)) * Ko)
    ref <- sqrt(pmax(cv, 0) + rowSums((J %*% Vo) * J))
    got <- if (g == "rows") {
      frm_linpred(fit, se.fit = TRUE)$se.fit
    } else {
      frm_linpred(fit, newdata = data.frame(x = xo), se.fit = TRUE)$se.fit
    }
    expect_lt(max(abs(got / ref - 1)), 0.02, label = g)
  }
  # ranef()'s condSD is the field's alone, without the intercept
  re <- ranef(fit, condVar = TRUE)[[1L]]
  sdv <- as.vector(attr(re, "condSD"))
  u <- sort(unique(x))
  field <- function(p) cond(p, u) - p[1]
  J <- vapply(seq_along(p0), function(i) {
    e <- replace(numeric(4L), i, 1e-5)
    (field(p0 + e) - field(p0 - e)) / 2e-5
  }, numeric(length(u)))
  Ko <- exp(2 * p0[2]) * exp(-outer(u, x, "-")^2 / (2 * exp(2 * p0[3])))
  V <- exp(2 * p0[2]) * exp(-outer(x, x, "-")^2 / (2 * exp(2 * p0[3])))
  diag(V) <- diag(V) + exp(2 * p0[4])
  cv <- exp(2 * p0[2]) - rowSums((Ko %*% solve(V)) * Ko)
  ref <- sqrt(pmax(cv, 0) + rowSums((J %*% Vo) * J))
  expect_lt(max(abs(sdv / ref - 1)), 0.02)
})

test_that("a factor that fails is NaN on the tape, never a number", {
  # At 1e-12 no length scale failed on 18 designs up to 500 positions
  # (dev/ncgp-cholscan.R), so the failure is forced by a smaller nugget
  d <- gpnc_low(1e-1, n = 60L)
  fit <- gpnc_fit(d)
  local_mocked_bindings(gp_nugget = 1e-17, .package = "frmtmb")
  fr <- fit$frame
  tpl <- fit$estimates
  # a length scale short enough to factor, for the tape to be built at
  tpl$theta[2L] <- log(1e-2)
  obj <- frmtmb:::gp_tape_scope(fr, RTMB::MakeADFun(
    frmtmb:::build_objective(fr), tpl, random = "b", silent = TRUE))
  good <- bad <- obj$par
  bad[which(names(bad) == "theta")[2L]] <- log(5)
  expect_false(is.nan(obj$fn(good)))
  expect_true(is.nan(suppressWarnings(obj$fn(bad))))
  expect_true(all(is.nan(suppressWarnings(obj$gr(bad)))))
  # and the next good point reads the same as before it
  expect_identical(obj$fn(good), obj$fn(good))
})

test_that("the sparse Hessian setting does not leak out of the scope", {
  loadNamespace("RTMB")
  before <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  tc0 <- RTMB::TapeConfig()
  d <- gpnc_low(1e-1, n = 30L)
  fit <- gpnc_fit(d)
  expect_identical(TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress,
                   before)
  expect_identical(RTMB::TapeConfig()$matmul, tc0$matmul)
  # a frame with no exact field takes neither, so its tape is the same
  seen <- NULL
  f0 <- frm(bf(y ~ x), data = d, family = gaussian())
  frmtmb:::gp_tape_scope(f0$frame, {
    seen <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  })
  expect_identical(seen, before)
  frmtmb:::gp_tape_scope(fit$frame, {
    seen <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  })
  expect_identical(seen, 1L)
})

test_that("REML, two fields, by and sigma ~ gp() reach the closed form", {
  d <- gpnc_low(1e-2, n = 40L)
  d$z <- round(stats::runif(40, 0, 5), 2)
  n <- nrow(d)
  # two fields: y ~ N(mu, K1 + K2 + sigma^2)
  fit <- allow_warnings(frm(bf(y ~ gp(x) + gp(z)), data = d,
                            family = gaussian()),
                        c("Optimizer did not report", "gradient"))
  th <- fit$estimates$theta
  bks <- fit$frame$re_blocks
  V <- diag(sigma(fit)[1L]^2, n)
  for (i in 1:2) {
    v <- d[[c("x", "z")[i]]]
    tt <- th[bks[[i]]$theta_idx]
    V <- V + exp(2 * tt[1]) * (exp(-outer(v, v, "-")^2 / (2 * exp(2 * tt[2]))) +
                                 1e-12 * outer(v, v, "==") )
  }
  R <- chol(V)
  r <- backsolve(R, d$y - fit$estimates$beta[[1L]], transpose = TRUE)
  ref <- 0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
  expect_lt(abs(fit$opt$objective - ref), 1e-6 * abs(ref))
  # REML integrates the intercept out on top
  fr <- allow_warnings(frm(bf(y ~ gp(x)), data = d, family = gaussian(),
                           REML = TRUE),
                       c("Optimizer did not report", "gradient"))
  bk <- fr$frame$re_blocks[[1L]]
  tt <- fr$estimates$theta[bk$theta_idx]
  V <- exp(2 * tt[1]) * (exp(-outer(d$x, d$x, "-")^2 / (2 * exp(2 * tt[2]))))
  diag(V) <- diag(V) * (1 + 1e-12) + sigma(fr)[1L]^2
  R <- chol(V)
  one <- backsolve(R, rep(1, n), transpose = TRUE)
  ry <- backsolve(R, d$y, transpose = TRUE)
  xvx <- sum(one^2)
  bh <- sum(one * ry) / xvx
  ref <- 0.5 * sum((ry - bh * one)^2) + sum(log(diag(R))) +
    0.5 * log(xvx) + 0.5 * (n - 1) * log(2 * pi)
  expect_lt(abs(fr$opt$objective - ref), 1e-5 * abs(ref))
})

test_that("a Poisson field's objective is the Laplace of z, by hand", {
  set.seed(4401)
  n <- 40L
  x <- sort(round(stats::runif(n, 0, 6), 2))
  y <- stats::rpois(n, exp(12 + sin(x)))
  d <- data.frame(x = x, y = y)
  fit <- allow_warnings(frm(bf(y ~ gp(x)), data = d, family = poisson()),
                        c("Optimizer did not report", "gradient"))
  bk <- fit$frame$re_blocks[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  b0 <- fit$estimates$beta[[1L]]
  u <- sort(unique(x))
  J <- match(x, u)
  K <- exp(-outer(u, u, "-")^2 / (2 * exp(2 * th[2])))
  diag(K) <- diag(K) + 1e-12
  p <- bk$gp_order
  L <- exp(th[1]) * t(chol(K[p, p]))[order(p), , drop = FALSE]
  Lr <- L[J, , drop = FALSE]
  z <- fit$estimates$b[bk$b_idx]
  # ten Newton steps from frm()'s mode change nothing, and the value reads
  for (i in 1:10) {
    eta <- b0 + as.vector(Lr %*% z)
    g <- z - as.vector(crossprod(Lr, y - exp(eta)))
    H <- diag(length(z)) + crossprod(Lr * sqrt(exp(eta)))
    z <- z - solve(H, g)
  }
  eta <- b0 + as.vector(Lr %*% z)
  H <- diag(length(z)) + crossprod(Lr * sqrt(exp(eta)))
  ref <- 0.5 * sum(z^2) + 0.5 * length(z) * log(2 * pi) -
    sum(stats::dpois(y, exp(eta), log = TRUE)) +
    0.5 * as.numeric(determinant(H)$modulus) -
    0.5 * length(z) * log(2 * pi)
  expect_lt(abs(fit$opt$objective - ref), 1e-6 * abs(ref))
})

test_that("a central difference over theta reads the delta method", {
  set.seed(4402)
  n <- 40L
  x <- sort(round(stats::runif(n, 0, 6), 2))
  d <- data.frame(x = x, y = stats::rpois(n, exp(1.5 + sin(2 * x))))
  fit <- allow_warnings(frm(bf(y ~ gp(x)), data = d, family = poisson()),
                        c("Optimizer did not report", "gradient"))
  sm <- frmtmb:::smooth_b_idx(fit)
  ch <- frmtmb:::fd_eta_chain(fit, NULL, NULL, FALSE, TRUE, sm)
  f <- function(g) frmtmb:::fitted_point(g, NULL, NA)
  fd <- as.vector(frmtmb:::fit_fd_se(fit, f, b_idx = sm, b_chain = ch$chain))
  de <- frm_linpred(fit, type = "response", se.fit = TRUE)$se.fit
  expect_lt(max(abs(fd / de - 1)), 1e-3)
})

test_that("d L / d theta reproduces d K / d theta", {
  # A central difference at fixed z reads the factor's rounding: replaced by
  # one, fitted standard errors moved by 2e-2 at noise 1e-4, n 200
  # (dev/ncgp-findings.md). dL L' + L dL' = dK needs no reference but the kernel
  set.seed(3)
  n <- 200L
  x <- stats::runif(n, 0, 6)
  d <- data.frame(x = x, y = sin(2 * x) + 1e-4 * stats::rnorm(n))
  fit <- gpnc_fit(d)
  bk <- fit$frame$re_blocks[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  L <- frmtmb:::gp_factor(th, bk)
  m <- nrow(L)
  U <- vapply(seq_len(m), function(i) {
    frmtmb:::gp_field_dtheta(replace(numeric(m), i, 1), th, bk)[, 2L]
  }, numeric(m))
  dL <- L %*% U
  K0 <- exp(-bk$aux_D2[[1L]] / (2 * exp(2 * th[2L])))
  dK <- K0 * bk$aux_D2[[1L]] / exp(2 * th[2L])
  expect_lt(max(abs(tcrossprod(dL, L) + tcrossprod(L, dL) - dK)),
            1e-6 * max(abs(dK)))
})

test_that("b first gives the plain tape's gradient on a sorted hand tape", {
  # TMB's compressed Hessian reverses wrong with another inner parameter ahead of
  # the field: 0.232 against 0.0717 here, and 0.182 in a plain TMB C++ template
  # (dev/review-ncgp-tmbcpp.R, dev/upstream-bugs.md). Farthest-first takes that
  # to 8e-8 and does not remove it, so this tape is over sorted positions
  set.seed(20261102)
  n <- 40L
  x <- round(stats::runif(n, 0, 6), 3)
  y <- sin(2 * x) + 1e-2 * stats::rnorm(n)
  u <- sort(unique(x))
  J <- match(x, u)
  D2 <- outer(u, u, "-")^2
  stub <- list(re_blocks = list(list(covstruct = "gp")))
  f <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    L <- t(chol(exp(-D2 / (2 * exp(2 * p$lrho))) + diag(1e-12, length(u))))
    eta <- p$mu + exp(p$lsd) * as.vector(L %*% p$b)[J]
    -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) -
      sum(RTMB::dnorm(y, eta, exp(p$ls), log = TRUE))
  }
  pars <- list(mu = 0, ls = log(0.01), lsd = 1.08, lrho = 0.3,
               b = numeric(length(u)))
  at <- c(ls = log(0.01), lsd = 1.08396, lrho = 0.27)
  plain <- RTMB::MakeADFun(f, pars, random = c("mu", "b"), silent = TRUE)
  g0 <- plain$gr(at)
  v0 <- plain$fn(at)
  fast <- frmtmb:::gp_tape_scope(stub, RTMB::MakeADFun(
    f, frmtmb:::gp_tape_order(stub, pars), random = c("mu", "b"),
    silent = TRUE))
  # the outer parameters are named the same whichever order the list is in
  expect_identical(names(fast$par), names(plain$par))
  expect_lt(abs(fast$fn(at) - v0), 1e-8 * abs(v0))
  expect_lt(max(abs(fast$gr(at) - g0)), 1e-5 * max(abs(g0)))
})

test_that("a lost sdgp leaves a field prediction its standard error", {
  # the sd runs to zero on pure noise; the field's value does not move along
  # that direction. Holding z read all 60 rows as lost, NaN, where a centered
  # field gave them finite
  set.seed(4404)
  n <- 60L
  d <- data.frame(x = round(stats::runif(n, 0, 6), 3),
                  y = stats::rbinom(n, 1L, 0.5))
  fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = bernoulli()))
  # the construction holds: the direction was lost
  expect_gt(NCOL(frmtmb:::get_joint_cov(fit)$null), 0L)
  se <- suppressWarnings(frm_linpred(fit, se.fit = TRUE))$se.fit
  expect_true(all(is.finite(se)))
})

test_that("the compressed Hessian gives the plain tape's gradient", {
  # TMB's compressed reverse sweep is wrong with the intercept ahead of the
  # field in the parameter list, value right and gradient not
  # (dev/ncgp-remlarm2.R), so the tape puts b first. This pins it on every kind
  # of inner vector.
  set.seed(20261106)
  n <- 30L
  d <- data.frame(x = sort(round(stats::runif(n, 0, 6), 3)),
                  z = round(stats::runif(n, 0, 4), 2),
                  g = factor(rep(1:5, length.out = n)))
  d$y <- sin(2 * d$x) + 0.05 * stats::rnorm(n)
  d$cnt <- stats::rpois(n, exp(1 + sin(2 * d$x)))
  d$w <- d$z
  d$w[c(3L, 9L, 20L)] <- NA
  mk <- function(fit, comp) {
    old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
    on.exit(TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB"))
    TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "RTMB")
    fr <- fit$frame
    rnd <- unique(names(fit$obj$env$par)[fit$obj$env$random])
    o <- RTMB::MakeADFun(frmtmb::build_objective(fr),
                         frmtmb::gp_tape_order(fr, fit$estimates),
                         random = rnd, map = fr$map, silent = TRUE)
    force(o$env$spHess)
    o
  }
  specs <- list(list(bf(y ~ gp(x)), gaussian()),
                list(bf(y ~ gp(x) + (1 | g)), gaussian()),
                list(bf(y ~ gp(x) + gp(z)), gaussian()),
                list(bf(cnt ~ gp(x)), poisson()),
                list(bf(y ~ gp(x)), gaussian(), REML = TRUE),
                list(bf(y ~ z + gp(x) + (1 | g)), gaussian(), REML = TRUE),
                list(bf(y ~ mi(w) + gp(x)) + bf(w | mi() ~ 1), gaussian()))
  for (s in specs) {
    args <- c(list(s[[1L]], data = d, family = s[[2L]]), s[-(1:2)])
    fit <- suppressWarnings(do.call(frm, args))
    p <- fit$opt$par + 0.05
    g1 <- mk(fit, 1L)$gr(p)
    g0 <- mk(fit, 0L)$gr(p)
    expect_lt(max(abs(g1 - g0)), 1e-6 * max(abs(g0)))
  }
  # the frm() tape under REML, its gradient against its value at low noise
  dl <- transform(d, y = sin(2 * x) + 0.01 * (y - sin(2 * x)) / 0.05)
  fr <- suppressWarnings(frm(bf(y ~ gp(x)), data = dl, family = gaussian(),
                             REML = TRUE))
  o <- fr$obj
  p <- fr$opt$par + 0.01
  h <- 1e-4
  fd <- vapply(seq_along(p), function(i) {
    e <- replace(numeric(length(p)), i, h)
    (o$fn(p + e) - o$fn(p - e)) / (2 * h)
  }, 0)
  expect_lt(max(abs(o$gr(p) - fd)), 1e-3 * max(abs(fd)))
})

test_that("a kriging row keeps its digits at a condition of 1e13", {
  # positions on a grid at a long length scale. Ks K^-1 f, with K^-1 on a
  # near-singular kernel, against the same row through the factor alone
  set.seed(17)
  n <- 120L
  xg <- round(stats::runif(n, 0, 10), 1)
  pos <- sort(unique(xg))
  u <- drop(crossprod(chol(1.5^2 * exp(-outer(pos, pos, "-")^2 / 8) +
                             diag(1e-8, length(pos))),
                      stats::rnorm(length(pos))))
  d <- data.frame(y = 2 + u[match(xg, pos)] + stats::rnorm(n, 0, 0.5), x = xg)
  fit <- gpnc_fit(d)
  r <- gpnc_read(fit)
  po <- fit$frame$re_blocks[[1L]]$gp_order
  R <- chol((exp(-outer(pos, pos, "-")^2 / (2 * r$rho^2)) +
               diag(1e-12, length(pos)))[po, po])
  xs <- seq(0.05, 9.55, by = 0.5)
  P <- t(backsolve(R, t(exp(-outer(xs, pos[po], "-")^2 / (2 * r$rho^2))),
                   transpose = TRUE))
  want <- fit$estimates$beta[[1L]] + r$sdgp * drop(P %*% fit$estimates$b)
  got <- frm_linpred(fit, newdata = data.frame(x = xs))
  expect_lt(max(abs(got - want)), 1e-9 * max(abs(want)))
})

test_that("a jitter does not reach the tape", {
  # is.numeric() is TRUE on an advector, and a jitter on the tape is another
  # model, silently: a tape at 1e-16 read exactly as the tape at 1e-12
  # (dev/reviews/2026-11-02-gpgap.md, RF1). A failing start is an error here
  d <- gpnc_low(1e-1, n = 30L)
  fit <- gpnc_fit(d)
  bk <- fit$frame$re_blocks[[1L]]
  # a long length scale, where the factor needs the nugget
  th <- c(fit$estimates$theta[bk$theta_idx][1L], log(1e3))
  val <- function(nug) {
    local_mocked_bindings(gp_nugget = nug, .package = "frmtmb")
    tryCatch({
      tp <- RTMB::MakeTape(function(p) {
        L <- frmtmb:::gp_factor_ordered(p, bk)
        sum(log(diag(L)))
      }, th)
      tp(th)
    }, error = function(e) NA_real_)
  }
  v12 <- val(1e-12)
  v16 <- val(1e-16)
  expect_true(is.finite(v12))
  expect_false(isTRUE(all.equal(v12, v16, tolerance = 1e-8)))
})
