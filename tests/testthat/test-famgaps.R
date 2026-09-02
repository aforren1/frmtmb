# v0.35: the three brms families frmtmb was missing - categorical() on a
# bare factor, von_mises() for circular responses, and cox() with a
# flexible M-spline baseline hazard.

# ---------------------------------------------------------- categorical

sim_categorical <- function(seed = 11, n = 400) {
  set.seed(seed)
  x <- rnorm(n)
  E <- cbind(0, 0.5 + 0.8 * x, -0.3 + 0.4 * x)
  P <- exp(E) / rowSums(exp(E))
  y <- vapply(seq_len(n), function(i) {
    sample(c("a", "b", "c"), 1L, prob = P[i, ])
  }, "")
  data.frame(x = x, z = rnorm(n), y = factor(y))
}

test_that("categorical() matches nnet::multinom exactly", {
  skip_if_not_installed("nnet")
  dd <- sim_categorical()
  fit <- frm(bf(y ~ x), family = categorical(), data = dd)
  ref <- nnet::multinom(y ~ x, data = dd, trace = FALSE,
                        reltol = 1e-14, maxit = 1000)
  # both likelihoods are the plain categorical pmf: no multinomial
  # coefficient to reconcile, unlike the count-matrix spelling
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-6)
  cf <- coef(ref)
  expect_vector_equal(fixef(fit)$mub, cf["b", ], tol = 1e-4)
  expect_vector_equal(fixef(fit)$muc, cf["c", ], tol = 1e-4)
})

test_that("categorical() equals multinomial() on the one-hot response", {
  dd <- sim_categorical(seed = 12)
  fit <- frm(bf(y ~ x), family = categorical(), data = dd)
  Y <- stats::model.matrix(~ y - 1, dd)
  colnames(Y) <- levels(dd$y)
  d2 <- data.frame(x = dd$x)
  d2$Y <- Y
  fm <- frm(bf(Y ~ x), family = multinomial(K = 3), data = d2)
  # one trial per row, so the multinomial coefficient is log(1) = 0 and
  # the two log-likelihoods are the same number, not merely proportional
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(fm))), 1e-8)
  expect_vector_equal(unname(fixef(fit)$mub), unname(fixef(fm)$mu2),
                      tol = 1e-6)
  expect_vector_equal(unname(fixef(fit)$muc), unname(fixef(fm)$mu3),
                      tol = 1e-6)
})

test_that("categorical() dpars follow the brms mu<Level> spelling", {
  dd <- sim_categorical(seed = 13)
  fit <- frm(bf(y ~ x), family = categorical(), data = dd)
  expect_named(fixef(fit), c("mub", "muc"))
  # the first level is the reference and gets no predictor
  expect_false("mua" %in% names(fixef(fit)))
  # one category may take its own formula; the rest keep the main one
  f2 <- frm(bf(y ~ x, muc ~ z), family = categorical(), data = dd)
  expect_named(fixef(f2)$mub, c("(Intercept)", "x"))
  expect_named(fixef(f2)$muc, c("(Intercept)", "z"))
  # relevel and the reference moves with it
  dd$y2 <- factor(dd$y, levels = c("c", "a", "b"))
  f3 <- frm(bf(y2 ~ x), family = categorical(), data = dd)
  expect_named(fixef(f3), c("mua", "mub"))
  expect_equal(as.numeric(logLik(f3)), as.numeric(logLik(fit)),
               tolerance = 1e-6)
})

test_that("categorical() accepts K= and character responses", {
  dd <- sim_categorical(seed = 14)
  fit <- frm(bf(y ~ x), family = categorical(), data = dd)
  dd$yi <- as.integer(dd$y)
  fk <- frm(bf(yi ~ x), family = categorical(K = 3), data = dd)
  expect_named(fixef(fk), c("mu2", "mu3"))
  expect_equal(as.numeric(logLik(fk)), as.numeric(logLik(fit)),
               tolerance = 1e-6)
  fl <- frm(bf(y ~ x), family = categorical(levels = c("a", "b", "c")),
            data = dd)
  expect_equal(as.numeric(logLik(fl)), as.numeric(logLik(fit)),
               tolerance = 1e-6)
  dd$ych <- as.character(dd$y)
  expect_message(frm(bf(ych ~ x), family = categorical(), data = dd),
                 "read as a factor with levels a, b, c")
})

test_that("categorical() fitted/predict give the n x K probabilities", {
  dd <- sim_categorical(seed = 15)
  fit <- frm(bf(y ~ x), family = categorical(), data = dd)
  P <- fitted(fit)
  expect_equal(dim(P), c(nrow(dd), 3L))
  expect_equal(colnames(P), c("a", "b", "c"))
  expect_equal(unname(rowSums(P)), rep(1, nrow(dd)), tolerance = 1e-12)
  expect_equal(P, predict(fit, type = "response"))
  # the probabilities are the family's own softmax of the two predictors
  eb <- predict(fit, type = "link", dpar = "mub")
  ec <- predict(fit, type = "link", dpar = "muc")
  den <- 1 + exp(eb) + exp(ec)
  expect_equal(unname(P[, "a"]), unname(1 / den), tolerance = 1e-10)
  expect_equal(unname(P[, "c"]), unname(exp(ec) / den), tolerance = 1e-10)
  nd <- data.frame(x = c(-1, 0, 1))
  Pn <- predict(fit, newdata = nd, type = "response")
  expect_equal(dim(Pn), c(3L, 3L))
  expect_equal(unname(rowSums(Pn)), rep(1, 3), tolerance = 1e-12)
  expect_error(predict(fit, type = "response", se.fit = TRUE),
               "se.fit is not supported on the response scale for a")
})

test_that("categorical() simulates factor levels and refuses residuals", {
  dd <- sim_categorical(seed = 16)
  fit <- frm(bf(y ~ x), family = categorical(), data = dd)
  s <- simulate(fit, nsim = 2, seed = 7)
  expect_true(is.factor(s$sim_1))
  expect_false(is.ordered(s$sim_1))   # nominal, so no claimed order
  expect_equal(levels(s$sim_1), c("a", "b", "c"))
  # the draws follow the fitted probabilities
  P <- fitted(fit)
  expect_lt(max(abs(prop.table(table(s$sim_1)) - colMeans(P))), 0.06)
  expect_error(residuals(fit), "not defined for a categorical family")
})

test_that("categorical() takes random effects and marginaleffects", {
  dd <- sim_categorical(seed = 17, n = 300)
  dd$g <- factor(rep(seq_len(30), 10))
  set.seed(18)
  u <- rnorm(30, 0, 1.2)
  E <- cbind(0, 0.5 + 0.8 * dd$x + u[as.integer(dd$g)], -0.3 + 0.4 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$y <- factor(c("a", "b", "c")[vapply(seq_len(nrow(dd)), function(i) {
    sample.int(3L, 1L, prob = P[i, ])
  }, 1L)])
  fit <- frm(bf(y ~ x + (1 | g)), family = categorical(), data = dd)
  vc <- VarCorr(fit)
  expect_true(length(vc) >= 1L)
  expect_equal(dim(fitted(fit)), c(nrow(dd), 3L))

  skip_if_not_installed("marginaleffects")
  # a categorical outcome predicts a DISTRIBUTION per row, which
  # marginaleffects keys with a `group` column
  p <- marginaleffects::predictions(fit, newdata = dd[1:5, ])
  expect_true("group" %in% names(p))
  expect_setequal(unique(as.character(p$group)), c("a", "b", "c"))
  expect_equal(nrow(p) %% 3L, 0L)
  sl <- marginaleffects::avg_slopes(fit, variables = "x")
  expect_setequal(as.character(sl$group), c("a", "b", "c"))
  # the reference category's probability falls as x rises, and the
  # probabilities are a simplex, so the slopes sum to zero
  expect_lt(abs(sum(sl$estimate)), 1e-8)
})

test_that("categorical() validation", {
  dd <- sim_categorical(seed = 19)
  expect_error(categorical(levels = "a"), "at least two distinct")
  expect_error(categorical(K = 1), "at least two categories")
  expect_error(categorical(link = "probit"), "'logit' link only")
  expect_error(frm(bf(y ~ x), family = categorical(levels = c("a", "b")),
                   data = dd),
               "not among the family's categories")
  dd$const <- factor(rep("a", nrow(dd)), levels = c("a", "b"))
  expect_error(frm(bf(const ~ x), family = categorical(), data = dd),
               "only one value")
})

# ------------------------------------------------------------ von Mises

test_that("von_mises() matches a hand-rolled circular ML fit", {
  set.seed(21)
  n <- 500
  x <- rnorm(n)
  mu <- 2 * atan(0.3 + 0.9 * x)
  y <- frmtmb:::rvon_mises(n, mu, 3)
  dd <- data.frame(x = x, y = y)
  fit <- frm(bf(y ~ x), family = von_mises(), data = dd)

  # the same density written out by hand: kappa cos(y - mu) minus the
  # log normalizer 2 pi I0(kappa), with mu through the tan-half link
  nll_ref <- function(p) {
    m <- 2 * atan(p$b[1] + p$b[2] * x)
    k <- exp(p$lk)
    -sum(k * cos(y - m) - log(2 * pi) -
           (log(RTMB::besselI(k, 0, expon.scaled = TRUE)) + k))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), lk = 0), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(rel.tol = 1e-12, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  expect_vector_equal(fixef(fit)$mu, opt$par[1:2], tol = 1e-4)
  expect_vector_equal(fixef(fit)$kappa, opt$par[3], tol = 1e-4)
})

test_that("von_mises() matches circular::mle.vonmises on an intercept fit", {
  skip_if_not_installed("circular")
  set.seed(22)
  y <- frmtmb:::rvon_mises(600, 0.6, 2.5)
  dd <- data.frame(y = y)
  fit <- frm(bf(y ~ 1), family = von_mises(), data = dd)
  ref <- circular::mle.vonmises(circular::circular(y))
  # mu comes back through the tan-half link; circular's kappa uses a
  # different root finder, so it agrees to about three digits
  expect_equal(2 * atan(unname(fixef(fit)$mu)), as.numeric(ref$mu),
               tolerance = 1e-5)
  expect_equal(exp(unname(fixef(fit)$kappa)), as.numeric(ref$kappa),
               tolerance = 1e-2)
})

test_that("von_mises() takes a distributional kappa", {
  set.seed(23)
  n <- 700
  x <- rnorm(n)
  kap <- exp(0.6 + 0.8 * x)
  y <- frmtmb:::rvon_mises(n, rep(0.4, n), kap)
  dd <- data.frame(x = x, y = y)
  fit <- frm(bf(y ~ 1, kappa ~ x), family = von_mises(), data = dd)
  expect_equal(unname(fixef(fit)$kappa), c(0.6, 0.8), tolerance = 0.15)

  nll_ref <- function(p) {
    m <- 2 * atan(p$b0)
    k <- exp(p$lk[1] + p$lk[2] * x)
    -sum(k * cos(y - m) - log(2 * pi) -
           (log(RTMB::besselI(k, 0, expon.scaled = TRUE)) + k))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b0 = 0, lk = c(0, 0)),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(rel.tol = 1e-12, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("von_mises() surface: link, fitted, simulate, validation", {
  set.seed(24)
  y <- frmtmb:::rvon_mises(300, 1.0, 4)
  dd <- data.frame(y = y, x = rnorm(300))
  fit <- frm(bf(y ~ x), family = von_mises(), data = dd)
  expect_equal(family(fit)$links$mu$name, "tan_half")
  # tan_half maps the whole line onto the circle
  ft <- fitted(fit)
  expect_true(all(ft > -pi & ft <= pi))
  expect_equal(unname(ft), unname(2 * atan(predict(fit, type = "link"))),
               tolerance = 1e-10)
  s <- simulate(fit, nsim = 1, seed = 3)[[1]]
  expect_true(all(s > -pi & s <= pi))
  expect_error(frm(bf(y ~ x), family = von_mises(),
                   data = transform(dd, y = y + 10)),
               "angles in radians")
})

test_that("the von Mises simulator draws from the density it declares", {
  set.seed(25)
  y <- frmtmb:::rvon_mises(20000, 0.7, 3)
  # the resultant vector recovers the mean direction and, through
  # A(kappa) = I1/I0, the concentration
  expect_equal(atan2(mean(sin(y)), mean(cos(y))), 0.7, tolerance = 0.02)
  rbar <- sqrt(mean(sin(y))^2 + mean(cos(y))^2)
  expect_equal(rbar, besselI(3, 1) / besselI(3, 0), tolerance = 0.01)
  # kappa = 0 is the uniform distribution on the circle
  y0 <- frmtmb:::rvon_mises(5000, 0, 0)
  expect_lt(sqrt(mean(sin(y0))^2 + mean(cos(y0))^2), 0.05)
})

# ------------------------------------------------------------------ cox

test_that("the M-spline and I-spline bases match splines2", {
  skip_if_not_installed("splines2")
  set.seed(31)
  y <- rexp(80) + 0.1
  sp <- frmtmb:::bhaz_spec(y, df = 5, degree = 3, intercept = TRUE)
  M <- frmtmb:::mspline_design(y, sp)
  I <- frmtmb:::ispline_design(y, sp)
  M2 <- splines2::mSpline(y, df = 5, degree = 3, intercept = TRUE,
                          Boundary.knots = sp$boundary)
  I2 <- splines2::iSpline(y, df = 5, degree = 3, intercept = TRUE,
                          Boundary.knots = sp$boundary)
  expect_lt(max(abs(M - as.matrix(M2))), 1e-12)
  expect_lt(max(abs(I - as.matrix(I2))), 1e-12)
  expect_equal(sp$internal, unname(attr(M2, "knots")))
  # the same without the boundary basis function
  sp0 <- frmtmb:::bhaz_spec(y, df = 5, degree = 3, intercept = FALSE)
  M0 <- frmtmb:::mspline_design(y, sp0)
  M02 <- splines2::mSpline(y, df = 5, degree = 3, intercept = FALSE,
                           Boundary.knots = sp0$boundary)
  expect_lt(max(abs(M0 - as.matrix(M02))), 1e-12)
})

test_that("the I-spline basis is the integral of the M-spline basis", {
  set.seed(32)
  y <- rexp(60) + 0.2
  sp <- frmtmb:::bhaz_spec(y, df = 5, degree = 3, intercept = TRUE)
  xt <- stats::median(y)
  for (j in 1:5) {
    num <- stats::integrate(function(u) {
      frmtmb:::mspline_design(u, sp)[, j]
    }, sp$boundary[1], xt, subdivisions = 2000)$value
    expect_equal(num, frmtmb:::ispline_design(xt, sp)[, j],
                 tolerance = 1e-8)
  }
  # each M-spline integrates to one over its whole support, which is
  # what makes a simplex of weights a hazard rather than an arbitrary
  # positive combination
  Ifull <- frmtmb:::ispline_design(sp$boundary[2], sp)
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
  sp <- frmtmb:::bhaz_spec(dd$time, df = 5, degree = 3, intercept = TRUE)
  Zb <- frmtmb:::mspline_design(dd$time, sp)
  Zc <- frmtmb:::ispline_design(dd$time, sp)
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
  expect_vector_equal(fixef(fit)$mu, opt$par[1:2], tol = 1e-4)
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
  expect_equal(unname(fixef(fit)$mu["x"]), unname(coef(ref)["x"]),
               tolerance = 2e-2)
  # a bigger basis moves toward coxph rather than away from it
  f9 <- frm(bf(time | cens(cens) ~ x), family = cox(df = 9), data = dd)
  expect_equal(unname(fixef(f9)$mu["x"]), unname(coef(ref)["x"]),
               tolerance = 2e-2)
})

test_that("cox() with no censoring is a plain flexible parametric PH fit", {
  set.seed(35)
  n <- 400
  x <- rnorm(n)
  dd <- data.frame(x = x,
                   time = stats::rexp(n, exp(-0.4 + 0.6 * x)))
  fit <- frm(bf(time ~ x), family = cox(), data = dd)
  expect_equal(unname(fixef(fit)$mu["x"]), 0.6, tolerance = 0.15)
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
  fit <- frm(bf(time | cens(cens) ~ x + (1 | g)), family = cox(),
             data = dd)
  sd_hat <- sqrt(as.numeric(VarCorr(fit)[[1L]]))
  expect_equal(sd_hat, 0.8, tolerance = 0.3)
  expect_equal(unname(fixef(fit)$mu["x"]), 0.7, tolerance = 0.2)
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
  expect_equal(unname(fixef(fit)$mu["x"]), unname(coef(ref)["x"]),
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
  expect_error(predict(fit, type = "response"), "no mean on the response")
  expect_error(simulate(fit), "has no simulator")
  # the log hazard ratio is still there
  expect_length(predict(fit, type = "link"), nrow(dd))
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
