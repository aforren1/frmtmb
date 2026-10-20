# The outer Hessian the covariance is built from, at a covariate's own
# scale: the AD Hessian on a fit without random effects, central
# differences with steps scaled to the design columns with them
# (fit_outer_hessian(), outer_fd_steps(); lane adhess,
# dev/adhess-findings.md). 0.71.0 differenced the gradient with an
# absolute step of 1e-3, wrong wherever that step moves a predictor by
# much and the gradient is not linear along it, and nothing warned.

ah_warns <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

ah_cov_err <- function(V, ref) {
  s <- sqrt(diag(ref))
  c(se = max(abs(sqrt(diag(V)) / s - 1)),
    cor = max(abs(V - ref) / outer(s, s)))
}

# What rounding alone leaves between two inverses of the same
# information in double: the unit roundoff times the condition number
# of the correlation matrix, which governs the relative error of SEs
# and correlations, with the house's factor of 1e3
ah_round_tol <- function(ref) {
  s <- sqrt(diag(ref))
  1e3 * .Machine$double.eps * kappa(ref / outer(s, s), exact = TRUE)
}

# The error 0.71.0's Hessian had at the same point: optimHess()'s
# absolute step of 1e-3 (fd_hessian() at its default), so that a bound
# is a ratio to what the old route measures on this very fit
ah_old_err <- function(fit, at, ref) {
  ns <- asNamespace("frmtmb")
  st <- ns$obj_state_save(fit$obj)
  on.exit(ns$obj_state_restore(fit$obj, st))
  ah_cov_err(solve(ns$fd_hessian(at, fit$obj$gr)$H), ref)
}

test_that("a poisson GLM's covariance is the inverse of its information", {
  # 0.71.0: SEs 6.2 percent off, the correlation 0.12, no warning
  set.seed(7)
  d <- data.frame(x = stats::rnorm(300, 500, 100))
  d$y <- stats::rpois(300, exp(0.3 + 0.004 * (d$x - 500)))
  r <- ah_warns(frm(y ~ x, family = poisson(), data = d))
  expect_length(r$warnings, 0L)
  b <- fixef(r$value)[, "Estimate"]
  X <- cbind(1, d$x)
  mu <- exp(drop(X %*% b))
  ref <- solve(crossprod(X, X * mu))
  e <- ah_cov_err(unname(vcov(r$value)), ref)
  expect_lt(e[["se"]], ah_round_tol(ref))
  expect_lt(e[["cor"]], ah_round_tol(ref))
})

test_that("a bernoulli GLM on three covariate scales matches its information", {
  # covariates of sd 1, 1e3 and 1e-2: 0.71.0 was 6.9 percent off
  # (dev/adhess-glm.R, design d07)
  set.seed(702)
  n <- 500
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n, 0, 1e3),
                  x3 = stats::rnorm(n, 0, 1e-2))
  X <- cbind(1, d$x1, d$x2, d$x3)
  eta <- drop(scale(X[, -1]) %*% rep(0.3, 3))
  d$y <- stats::rbinom(n, 1, stats::plogis(0.2 + eta))
  r <- ah_warns(frm(y ~ x1 + x2 + x3, family = bernoulli(), data = d))
  expect_length(r$warnings, 0L)
  b <- fixef(r$value)[, "Estimate"]
  p <- stats::plogis(drop(X %*% b))
  ref <- solve(crossprod(X, X * (p * (1 - p))))
  e <- ah_cov_err(unname(vcov(r$value)), ref)
  expect_lt(e[["se"]], ah_round_tol(ref))
  expect_lt(e[["cor"]], ah_round_tol(ref))
})

test_that("a nonlinear body on a shifted covariate gets its own curvature", {
  # y ~ a exp(b x), x on 20 + [0, 1]: 0.71.0's correlation was 0.133
  # off and its SEs 0.87 percent. The reference is the observed
  # information in closed form at the estimate.
  set.seed(1)
  d <- data.frame(x = 20 + stats::runif(200))
  d$y <- 2 * exp(0.8 * (d$x - 20)) + stats::rnorm(200, 0, 0.3)
  r <- ah_warns(frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d,
                    start = list(beta = c(2 * exp(-16), 0.8))))
  expect_length(r$warnings, 0L)
  fit <- r$value
  p <- unname(fit$opt$par)
  e <- exp(p[2] * d$x)
  res <- d$y - p[1] * e
  s2 <- exp(2 * p[3])
  fa <- e
  fb <- p[1] * d$x * e
  H <- matrix(0, 3, 3)
  H[1, 1] <- sum(fa^2)
  H[1, 2] <- H[2, 1] <- sum(fa * fb - res * d$x * e)
  H[2, 2] <- sum(fb^2 - res * p[1] * d$x^2 * e)
  H[1, 3] <- H[3, 1] <- -2 * sum(res * fa)
  H[2, 3] <- H[3, 2] <- -2 * sum(res * fb)
  H[3, 3] <- 2 * sum(res^2)
  ref <- solve(H / s2)
  err <- ah_cov_err(unname(vcov(fit, full = TRUE)), ref)
  # the information sums terms of opposite sign (the residual
  # curvature), so its rounding is far above that of a GLM's; the bound
  # is a ratio to the error of 0.71.0's differences on this fit
  old <- ah_old_err(fit, fit$opt$par, ref)
  expect_lt(err[["se"]], 1e-3 * old[["se"]])
  expect_lt(err[["cor"]], 1e-3 * old[["cor"]])
})

test_that("without random effects the outer Hessian is the AD one", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(100))
  d$y <- stats::rpois(100, exp(0.5 + 0.3 * d$x))
  fit <- frm(y ~ x, family = poisson(), data = d)
  ns <- asNamespace("frmtmb")
  h <- ns$fit_outer_hessian(fit)
  expect_true(isTRUE(h$exact))
  st <- ns$obj_state_save(fit$obj)
  Hx <- fit$obj$he(h$at)
  ns$obj_state_restore(fit$obj, st)
  expect_identical(unname(h$H), unname((Hx + t(Hx)) / 2))
  # and sdreport() inverted exactly that matrix
  expect_identical(unname(vcov(fit, full = TRUE)), unname(solve(h$H)))
})

ah_pois_re <- function(seed, xs = FALSE) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::rnorm(n, 500, 100),
                  g = factor(rep(1:20, length.out = n)))
  u <- stats::rnorm(20, 0, 0.3)
  d$y <- stats::rpois(n, exp(0.3 + 0.004 * (d$x - 500) + u[d$g]))
  d$xs <- (d$x - 500) / 100
  d
}

test_that("a GLMM's covariance does not depend on the covariate's units", {
  # x ~ N(500, 100) with (1 | g): 0.71.0's step of 1e-3 moved the
  # predictor by 0.5 per row, the SEs were 0.37 percent off and the SE
  # check called the group sd flat. The same model on (x - 500) / 100 is
  # an exact reparameterization, so its covariance, mapped back, is the
  # reference (dev/adhess-restep.R measured both within 3e-7 of a
  # Richardson Jacobian of the AD gradient).
  d <- ah_pois_re(13)
  r <- ah_warns(frm(y ~ x + (1 | g), family = poisson(), data = d))
  expect_length(r$warnings, 0L)
  rs <- ah_warns(frm(y ~ xs + (1 | g), family = poisson(), data = d))
  expect_length(rs$warnings, 0L)
  expect_length(sdr_of(r$value)$se_lost, 0L)
  # b0 = c0 - 5 c1, b1 = c1 / 100, theta = theta
  A <- rbind(c(1, -5, 0), c(0, 0.01, 0), c(0, 0, 1))
  ref <- A %*% unname(vcov(rs$value, full = TRUE)) %*% t(A)
  e <- ah_cov_err(unname(vcov(r$value, full = TRUE)), ref)
  at <- asNamespace("frmtmb")$fit_outer_hessian(r$value)$at
  old <- ah_old_err(r$value, at, ref)
  expect_lt(e[["se"]], 1e-2 * old[["se"]])
  expect_lt(e[["cor"]], 1e-2 * old[["cor"]])
})

test_that("the steps follow the design's columns and nothing else", {
  ns <- asNamespace("frmtmb")
  d <- ah_pois_re(13)
  fit <- frm(y ~ x + (1 | g), family = poisson(), data = d)
  pn <- names(fit$opt$par)
  h <- ns$outer_fd_steps(fit$frame, pn)
  expect_identical(h[[1]], 1e-3)
  expect_equal(h[[2]], 1e-3 / sqrt(mean(d$x^2)),
               tolerance = 64 * .Machine$double.eps)
  expect_identical(h[[3]], 1e-3)
  # in the coordinates par / u the column's scale is its own times u,
  # and a parameter outside every design keeps 1e-3 whatever its unit
  hu <- ns$outer_fd_steps(fit$frame, pn, u = c(1, 1e-2, 50))
  expect_identical(hu[[1]], 1e-3)
  expect_equal(hu[[2]], 1e-3 / (1e-2 * sqrt(mean(d$x^2))),
               tolerance = 64 * .Machine$double.eps)
  expect_identical(hu[[3]], 1e-3)
})

test_that("a GLMM whose design is unit scale keeps optimHess()'s steps", {
  # the guard's absent case: no design column is wider than one unit, so
  # the Hessian is the one 0.71.0 built, bit for bit
  set.seed(5)
  n <- 200
  d <- data.frame(x = stats::rnorm(n, 0, 0.5),
                  g = factor(rep(1:20, length.out = n)))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x + stats::rnorm(20, 0, 0.3)[d$g]))
  fit <- frm(y ~ x + (1 | g), family = poisson(), data = d)
  ns <- asNamespace("frmtmb")
  expect_true(all(ns$outer_fd_steps(fit$frame, names(fit$opt$par)) ==
                    1e-3))
  h <- ns$fit_outer_hessian(fit)
  st <- ns$obj_state_save(fit$obj)
  a <- stats::optimHess(h$at, fit$obj$fn, fit$obj$gr)
  ns$obj_state_restore(fit$obj, st)
  expect_identical(unname(h$H), unname(a))
})

test_that("a kinked log density shows no standard error of 0", {
  # asym_laplace() is piecewise linear in its location, so its AD
  # Hessian is [[0, B], [B', C]]: solve() succeeds, and the block for
  # sigma and quantile is exactly 0. These fits stop short (the gradient
  # check), which keeps sdreport()'s inverse, so the SEs were 0 and the
  # intervals of zero width (lane adhess review, finding 1). A variance
  # that is not positive is NaN, with its reason.
  set.seed(2)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
  ok <- c("did not report convergence", "Large maximum absolute gradient")
  ns <- asNamespace("frmtmb")
  for (f in list(bf(y ~ x), bf(y ~ x, quantile = 0.5))) {
    fit <- allow_warnings(frm(f, family = asym_laplace(), data = d), ok)
    v <- diag(allow_warnings(vcov(fit, full = TRUE), ok))
    expect_false(any(v <= 0, na.rm = TRUE))
    ci <- allow_warnings(confint(fit), ok)
    expect_false(any(ci[, "upr"] - ci[, "lwr"] <= 0, na.rm = TRUE))
    # every parameter without a variance has a reason
    lost <- sdr_of(fit)$se_lost
    nm <- ns$outer_par_names(fit)
    expect_true(all(nm[is.na(v)] %in% names(lost)))
  }
})

ah_quad_data <- function() {
  set.seed(31)
  n <- 300
  d <- data.frame(x = stats::rnorm(n, 500, 100),
                  g = factor(rep(1:30, length.out = n)))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * (d$x - 500) / 100 +
                               stats::rnorm(30, 0, 0.5)[d$g]))
  d$xs <- d$x / 100
  d
}

test_that("an importance fit's covariance does not depend on x's units", {
  # importance = and quadrature = TRUE fits have no inner random effects
  # in TMB's sense, so they take the AD Hessian too. On this design
  # 0.71.0's importance-fit intercept SE was 0.498 against 0.270 on x /
  # 100, 85 percent off (lane adhess review, finding 2). x / 100 is an
  # exact reparameterization, so its covariance, mapped back, is the
  # reference.
  d <- ah_quad_data()
  ns <- asNamespace("frmtmb")
  fit <- frm(y ~ x + (1 | g), family = poisson(), data = d,
             importance = 200L)
  fs <- frm(y ~ xs + (1 | g), family = poisson(), data = d,
            importance = 200L)
  h <- ns$fit_outer_hessian(fit)
  expect_true(isTRUE(h$exact))
  A <- diag(c(1, 0.01, 1))
  ref <- A %*% unname(vcov(fs, full = TRUE)) %*% t(A)
  e <- ah_cov_err(unname(vcov(fit, full = TRUE)), ref)
  old <- ah_old_err(fit, h$at, ref)
  expect_lt(e[["se"]], 1e-2 * old[["se"]])
  expect_lt(e[["cor"]], 1e-2 * old[["cor"]])
  q <- frm(y ~ x + (1 | g), family = poisson(), data = d,
           quadrature = TRUE)
  expect_true(isTRUE(ns$fit_outer_hessian(q)$exact))
})
