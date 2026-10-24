# cox(), moved from frmtmb's tests/testthat/test-famgaps.R with the family
# (v0.35 added it there, with categorical() and von_mises()).

# ------------------------------------------------------------------ cox

test_that("the M-spline and I-spline bases match splines2", {
  skip_if_not_installed("splines2")
  set.seed(31)
  y <- rexp(80) + 0.1
  sp <- bhaz_spec(y, df = 5, degree = 3, intercept = TRUE)
  M <- mspline_design(y, sp)
  I <- ispline_design(y, sp)
  M2 <- splines2::mSpline(y, df = 5, degree = 3, intercept = TRUE,
                          Boundary.knots = sp$boundary)
  I2 <- splines2::iSpline(y, df = 5, degree = 3, intercept = TRUE,
                          Boundary.knots = sp$boundary)
  expect_lt(max(abs(M - as.matrix(M2))), 1e-12)
  expect_lt(max(abs(I - as.matrix(I2))), 1e-12)
  expect_equal(sp$internal, unname(attr(M2, "knots")))
  # the same without the boundary basis function
  sp0 <- bhaz_spec(y, df = 5, degree = 3, intercept = FALSE)
  M0 <- mspline_design(y, sp0)
  M02 <- splines2::mSpline(y, df = 5, degree = 3, intercept = FALSE,
                           Boundary.knots = sp0$boundary)
  expect_lt(max(abs(M0 - as.matrix(M02))), 1e-12)
})

test_that("the I-spline basis is the integral of the M-spline basis", {
  set.seed(32)
  y <- rexp(60) + 0.2
  sp <- bhaz_spec(y, df = 5, degree = 3, intercept = TRUE)
  xt <- stats::median(y)
  for (j in 1:5) {
    num <- stats::integrate(function(u) {
      mspline_design(u, sp)[, j]
    }, sp$boundary[1], xt, subdivisions = 2000)$value
    expect_equal(num, ispline_design(xt, sp)[, j],
                 tolerance = 1e-8)
  }
  # each M-spline integrates to one over its whole support, which is
  # what makes a simplex of weights a hazard rather than an arbitrary
  # positive combination
  Ifull <- ispline_design(sp$boundary[2], sp)
  expect_equal(as.numeric(Ifull), rep(1, 5), tolerance = 1e-10)
})

sim_cox_data <- function(seed = 33, n = 400, shape = 1.4, beta = 0.8) {
  set.seed(seed)
  x <- rnorm(n)
  lam <- exp(-0.5 + beta * x)
  tt <- stats::rweibull(n, shape = shape, scale = lam^(-1 / shape))
  ct <- stats::rexp(n, 0.3)
  data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct),
             ev = as.numeric(tt <= ct), x = x)
}

test_that("cox() matches a hand-rolled M-spline PH likelihood exactly", {
  dd <- sim_cox_data()
  fit <- frm(bf(time | cens(cens) ~ x), family = cox(), data = dd)

  # the same likelihood built from scratch: the baseline hazard is the
  # M-spline basis over a softmax simplex, the cumulative baseline is
  # the I-spline basis over the same weights, an event contributes
  # log h + log S and a right-censored row log S alone
  sp <- bhaz_spec(dd$time, df = 5, degree = 3, intercept = TRUE)
  Zb <- mspline_design(dd$time, sp)
  Zc <- ispline_design(dd$time, sp)
  ev <- dd$ev
  X <- cbind(1, dd$x)
  nll_ref <- function(p) {
    "c" <- RTMB::ADoverload("c")
    e <- exp(c(0, p$s))
    s <- e / sum(e)
    mu <- exp(as.vector(X %*% p$b))
    bh <- as.vector(Zb %*% s)
    cb <- as.vector(Zc %*% s)
    -sum(ev * (log(bh) + log(mu)) - cb * mu)
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), s = rep(0, 4)),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(rel.tol = 1e-12, eval.max = 2000,
                               iter.max = 2000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  expect_vector_equal(fixef_by_dpar(fit)$mu, opt$par[1:2], tol = 1e-4)
  e <- exp(c(0, opt$par[3:6]))
  expect_vector_equal(unname(cox_baseline(fit)), e / sum(e), tol = 1e-4)
})

test_that("cox() agrees with survival::coxph on the hazard ratio", {
  skip_if_not_installed("survival")
  dd <- sim_cox_data(seed = 34, n = 600)
  fit <- frm(bf(time | cens(cens) ~ x), family = cox(), data = dd)
  ref <- survival::coxph(survival::Surv(time, ev) ~ x, data = dd)
  # coxph leaves the baseline fully nonparametric while this family
  # spends five spline weights on it, so the agreement is close but not
  # exact; the hand-rolled test above is the exact one
  expect_equal(unname(fixef_by_dpar(fit)$mu["x"]), unname(coef(ref)["x"]),
               tolerance = 2e-2)
  # a bigger basis moves toward coxph rather than away from it
  f9 <- frm(bf(time | cens(cens) ~ x), family = cox(df = 9), data = dd)
  expect_equal(unname(fixef_by_dpar(f9)$mu["x"]), unname(coef(ref)["x"]),
               tolerance = 2e-2)
})

test_that("cox() with no censoring is a plain flexible parametric PH fit", {
  set.seed(35)
  n <- 400
  x <- rnorm(n)
  dd <- data.frame(x = x,
                   time = stats::rexp(n, exp(-0.4 + 0.6 * x)))
  fit <- frm(bf(time ~ x), family = cox(), data = dd)
  expect_equal(unname(fixef_by_dpar(fit)$mu["x"]), 0.6, tolerance = 0.15)
  s <- cox_baseline(fit)
  expect_equal(sum(s), 1, tolerance = 1e-12)
  expect_true(all(s > 0))
})

test_that("cox() frailty models come out of the Laplace approximation", {
  set.seed(36)
  n_g <- 60
  n_i <- 8
  n <- n_g * n_i
  g <- factor(rep(seq_len(n_g), each = n_i))
  u <- rnorm(n_g, 0, 0.8)
  x <- rnorm(n)
  lam <- exp(-0.5 + 0.7 * x + u[as.integer(g)])
  tt <- stats::rexp(n, lam)
  ct <- stats::rexp(n, 0.2)
  dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct),
                   ev = as.numeric(tt <= ct), x = x, g = g)
  # a baseline-hazard coordinate at its limit has no standard error
  # (lane setier); the estimates are the point here
  fit <- allow_warnings(frm(bf(time | cens(cens) ~ x + (1 | g)),
                            family = cox(), data = dd),
                        "Standard errors are not available")
  sd_hat <- sqrt(as.numeric(varcorr_matrices(fit)[[1L]]))
  expect_equal(sd_hat, 0.8, tolerance = 0.3)
  expect_equal(unname(fixef_by_dpar(fit)$mu["x"]), 0.7, tolerance = 0.2)
  # the frailty buys likelihood over the same model without it. ML puts
  # two of this baseline's five weights on the simplex boundary, so the
  # Hessian is singular in their softmax directions and the optimizer
  # says so; the gradient is still zero and the coefficients are still
  # at the optimum, which is what the checks below assert.
  f0 <- suppressWarnings(
    frm(bf(time | cens(cens) ~ x), family = cox(), data = dd))
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(f0)) + 2)
  expect_lt(max(abs(f0$obj$gr(f0$opt$par))), 1e-3)
  expect_lt(min(cox_baseline(f0)), 1e-8)   # a weight went to the boundary
  expect_equal(sum(cox_baseline(f0)), 1, tolerance = 1e-12)

  skip_if_not_installed("survival")
  ref <- survival::coxph(
    survival::Surv(time, ev) ~ x + frailty(g, distribution = "gaussian"),
    data = dd)
  expect_equal(unname(fixef_by_dpar(fit)$mu["x"]), unname(coef(ref)["x"]),
               tolerance = 0.1)
})

test_that("cox() handles left and interval censoring, and refuses a mean", {
  dd <- sim_cox_data(seed = 37, n = 300)
  # left censoring: a row known only to have failed before its time
  dl <- dd
  dl$cens <- ifelse(seq_len(nrow(dl)) %% 7L == 0L, -1, dl$cens)
  fl <- frm(bf(time | cens(cens) ~ x), family = cox(), data = dl)
  expect_true(is.finite(as.numeric(logLik(fl))))
  # interval censoring needs the upper bound
  di <- dd
  di$hi <- di$time * 1.5
  di$cens <- ifelse(seq_len(nrow(di)) %% 5L == 0L, 2, di$cens)
  fi <- frm(bf(time | cens(cens, hi) ~ x), family = cox(), data = di)
  expect_true(is.finite(as.numeric(logLik(fi))))

  fit <- frm(bf(time | cens(cens) ~ x), family = cox(), data = dd)
  expect_error(fitted(fit), "no mean on the response scale")
  expect_error(frm_linpred(fit, type = "response"), "no mean on the response")
  expect_error(simulate(fit), "has no simulator")
  # the log hazard ratio is still there
  expect_length(frm_linpred(fit, type = "link"), nrow(dd))
})

test_that("cox() validation", {
  dd <- sim_cox_data(seed = 38, n = 120)
  expect_error(frm(bf(time ~ x), family = cox(df = 3), data = dd),
               "df must be at least degree")
  expect_error(frm(bf(time ~ x), family = cox(),
                   data = transform(dd, time = time - min(time))),
               "strictly positive")
  d2 <- data.frame(time = rep(c(1, 2, 3), 20), x = rnorm(60))
  expect_error(frm(bf(time ~ x), family = cox(), data = d2),
               "fewer distinct event times")
})

# ------------------------------------------ moved from other core files
#
# Each block below was an assertion about cox() in a core test file
# whose subject is wider than cox(); the core file keeps its subject on
# another family and says where this part went.

test_that("cox() is refused as a mixture component, as in brms", {
  # from test-brms-families.R
  expect_error(mixture(poisson, cox), "real and integer support")
})

test_that("cox() is the family whose simulate row is refused", {
  # from test-hurdle-nb-zoib.R
  expect_identical(frm_compat("cox", "simulate")$status, "refused")
})

test_that("cox() states why it has no simulator, at both entry points", {
  # from test-simulate-density.R: a refusal that gives the reason is the
  # difference between a gap and a decision, and cox()'s is a decision
  set.seed(12)
  dd <- data.frame(x = rep(c(0, 1), length.out = 40L))
  dd$y <- stats::rexp(40L, rate = exp(-0.2 + 0.5 * dd$x))
  expect_error(
    frm_simulate(y ~ x, dd, family = cox(),
                 newparams = list(b_Intercept = 0, b_x = 0.5),
                 nsim = 1L),
    "cumulative baseline hazard")
  expect_error(
    frm_simulate(y ~ x, dd, family = cox(),
                 newparams = list(b_Intercept = 0, b_x = 0.5),
                 nsim = 1L),
    "family 'cox' has no simulator", fixed = TRUE)
  fit <- frm(y ~ x, family = cox(), data = dd)
  expect_error(stats::simulate(fit, nsim = 1L),
               "family 'cox' has no simulator", fixed = TRUE)
})

test_that("cox() in a multivariate model stays refused", {
  # from test-mv-gaps.R: its baseline is an extra-parameter block the
  # post-fit methods read without a response in hand
  set.seed(4)
  n <- 40
  dd <- data.frame(t = rexp(n), cc = rbinom(n, 1, 0.3), x = rnorm(n),
                   y = rnorm(n))
  expect_error(frm(bf(t | cens(cc) ~ x) + cox() + bf(y ~ x) + gaussian(),
                   data = dd, dry_run = "frame"),
               "not supported in multivariate fits yet. The ordinal")
})

test_that("update() refits a cox(df = 6) fit with the same family", {
  # from test-update-pool.R, which checks that the family a constructor
  # call cannot say is stored as the object; here through update(), the
  # public route that reads it
  set.seed(13)
  dd <- data.frame(x = stats::rnorm(150), z = stats::rnorm(150))
  dd$time <- stats::rexp(150, exp(-0.3 + 0.6 * dd$x))
  f0 <- frm(time ~ x, family = cox(df = 6), data = dd)
  expect_length(cox_baseline(f0), 6L)
  up <- update(f0, . ~ . + z)
  ref <- frm(time ~ x + z, family = cox(df = 6), data = dd)
  expect_length(cox_baseline(up), 6L)
  expect_identical(logLik(up), logLik(ref))
})
