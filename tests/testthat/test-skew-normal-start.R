# skew_normal() has a stationary point at alpha = 0 whose information is
# singular, so a fit that stops there converges cleanly and reports a
# logLik tens of units below the maximum. Two things keep a fit off it:
# a start taken from the RESIDUAL skew rather than the response's, and a
# refit from both sides when the optimum lands on the point anyway.
#
# Every tolerance below is a ratio to a quantity the run itself
# measures, never an absolute.

# The design the defect was found on: a covariate skewed one way and a
# residual skewed the other, so the marginal and the conditional skew
# have opposite signs.
sn_stall_data <- function(seed, n = 200) {
  set.seed(seed)
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  data.frame(y = y, xs = xs)
}
sn_form <- function() bf(y ~ xs, sigma ~ 1, alpha ~ 1)
sn_skew <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3

test_that("the alpha start reads the residual skew, not the response's", {
  dd <- sn_stall_data(1)
  # the two disagree in sign on this design, which is the whole defect
  expect_lt(sn_skew(dd$y) * sn_skew(residuals(lm(y ~ xs, data = dd))), 0)
  u <- frm(sn_form(), family = skew_normal(), data = dd,
           dry_run = "objective")
  a0 <- u$estimates$betad[["alpha_(Intercept)"]]
  expect_gt(a0 * sn_skew(residuals(lm(y ~ xs, data = dd))), 0)
})

test_that("the sigma start is the residual sd, not the marginal sd", {
  dd <- sn_stall_data(1)
  r <- residuals(lm(y ~ xs, data = dd))
  u <- frm(sn_form(), family = skew_normal(), data = dd,
           dry_run = "objective")
  s0 <- u$estimates$betad[["sigma_(Intercept)"]]
  expect_equal(s0, log(stats::sd(r)), tolerance = 1e-8)
  # the marginal sd is the number that used to be used, and on this
  # design it is a different number, not a rounding of the same one
  expect_gt(abs(log(stats::sd(dd$y)) - s0), 0.1 * abs(s0 + 1))
})

test_that("skew_normal reaches the maximum where it used to stop at zero", {
  dd <- sn_stall_data(1)
  fit <- frm(sn_form(), family = skew_normal(), data = dd)
  # the reference is this run's own best optimum over both sides of the
  # stationary point, so nothing absolute is pinned
  ref <- max(vapply(c(2, -2), function(a) {
    as.numeric(logLik(frm(sn_form(), family = skew_normal(), data = dd,
                          start = list(betad = c(log(stats::sd(dd$y)), a)))))
  }, 0))
  ll <- as.numeric(logLik(fit))
  expect_lt(ref - ll, abs(ref) * 1e-8)
  # and it is nowhere near the stationary point
  expect_gt(abs(fit$opt$par[[4]]), 1)
})

test_that("a fit driven onto the stationary point is refit off it", {
  dd <- sn_stall_data(1)
  # the start the family used before: the sign of the RAW skew
  m3 <- sn_skew(dd$y)
  bad <- c(log(stats::sd(dd$y)), 2 * sign(m3) + 0.5 * m3)
  fit <- frm(sn_form(), family = skew_normal(), data = dd,
             start = list(betad = bad))
  esc <- fit$opt[["stationary_escape"]]
  expect_false(is.null(esc))
  expect_identical(unname(esc[["starts"]]), 2)
  expect_gt(esc[["gain"]], 1)
  ref <- as.numeric(logLik(frm(sn_form(), family = skew_normal(), data = dd)))
  expect_lt(ref - as.numeric(logLik(fit)), abs(ref) * 1e-8)
})

test_that("a fit that is not on the stationary point pays no refit", {
  dd <- sn_stall_data(1)
  fit <- frm(sn_form(), family = skew_normal(), data = dd)
  expect_null(fit$opt[["stationary_escape"]])
})

test_that("an intercept-only skew_normal keeps the marginal sd start", {
  set.seed(4)
  n <- 300
  dd <- data.frame(y = 2 + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5)
  u <- frm(bf(y ~ 1, sigma ~ 1, alpha ~ 1), family = skew_normal(),
           data = dd, dry_run = "objective")
  # with no covariate the residual IS the centred response, so the
  # start must not have moved
  expect_identical(u$estimates$betad[["sigma_(Intercept)"]],
                   log(stats::sd(dd$y)))
})

test_that("the escape also runs under REML and profile", {
  dd <- sn_stall_data(3)
  m3 <- sn_skew(dd$y)
  bad <- c(log(stats::sd(dd$y)), 2 * sign(m3) + 0.5 * m3)
  for (ctl in list(list(REML = TRUE),
                   list(control = frmtmb_control(profile = TRUE)))) {
    args <- c(list(sn_form(), family = skew_normal(), data = dd), ctl)
    ref <- as.numeric(logLik(do.call(frm, args)))
    args$start <- list(betad = bad)
    fit <- do.call(frm, args)
    expect_lt(ref - as.numeric(logLik(fit)), abs(ref) * 1e-8)
    # REML and profile drop beta from the outer vector, so read alpha
    # from the estimates rather than from a position in opt$par
    expect_gt(abs(fit$estimates$betad[["alpha_(Intercept)"]]), 1)
  }
})

test_that("a two-argument init_dpars function is still called with two", {
  # nothing already written may have to change: a custom family whose
  # initializer takes (y, aterms) must keep working
  fam <- custom_family(
    "two_arg", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      -0.5 * ((y - dpars[["mu"]]) / dpars[["sigma"]])^2 -
        log(dpars[["sigma"]]) - 0.5 * log(2 * pi)
    },
    init_dpars = list(mu = function(y, aterms) mean(y) + 7,
                      sigma = function(y, aterms) stats::sd(y)),
    post = list(mean_fn = function(dpars, aterms) dpars[["mu"]])
  )
  set.seed(5)
  dd <- data.frame(y = rnorm(100), x = rnorm(100))
  u <- frm(bf(y ~ x, sigma ~ 1), family = fam, data = dd,
           dry_run = "objective")
  expect_identical(u$estimates$beta[["(Intercept)"]], mean(dd$y) + 7)
})

# A family declaring post$stationary is the public half of this, so the
# machinery is exercised through one rather than only through
# skew_normal(). sn_norm_family() is an ordinary normal likelihood: it
# has no stationary point, so a declaration on it fires only because
# the window is made absurd, and the fit must come out unchanged.
sn_norm_family <- function(stationary = NULL) {
  custom_family(
    "sn_norm", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      -0.5 * ((y - dpars[["mu"]]) / dpars[["sigma"]])^2 -
        log(dpars[["sigma"]]) - 0.5 * log(2 * pi)
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) stats::sd(y)),
    post = c(list(mean_fn = function(dpars, aterms) dpars[["mu"]]),
             if (!is.null(stationary)) list(stationary = stationary))
  )
}

test_that("a refit that finds nothing leaves the fit exactly as it was", {
  set.seed(6)
  dd <- data.frame(y = rnorm(150, 3, 2), x = rnorm(150))
  fo <- bf(y ~ x, sigma ~ 1)
  plain <- frm(fo, family = sn_norm_family(), data = dd)
  # a window wide enough to fire on every fit, so the refit always runs
  always <- sn_norm_family(list(sigma = list(at = 0, tol = 1e6,
                                             from = c(1, -1))))
  forced <- frm(fo, family = always, data = dd)
  esc <- forced$opt[["stationary_escape"]]
  expect_false(is.null(esc))
  expect_identical(unname(esc[["starts"]]), 2)
  # it found nothing, and it must not have cost anything either
  expect_lt(abs(esc[["gain"]]), abs(plain$opt$objective) * 1e-8)
  expect_lt(abs(as.numeric(logLik(forced)) - as.numeric(logLik(plain))),
            abs(as.numeric(logLik(plain))) * 1e-8)
})

test_that("a malformed post$stationary is ignored, not obeyed", {
  set.seed(7)
  dd <- data.frame(y = rnorm(120, 1, 1.5), x = rnorm(120))
  fo <- bf(y ~ x, sigma ~ 1)
  ref <- frm(fo, family = sn_norm_family(), data = dd)
  bad <- list(
    character_from = list(sigma = list(at = 0, tol = 1e6, from = "two")),
    empty_from     = list(sigma = list(at = 0, tol = 1e6,
                                       from = numeric(0))),
    na_tol         = list(sigma = list(at = 0, tol = NA_real_,
                                       from = c(1, -1))),
    negative_tol   = list(sigma = list(at = 0, tol = -1, from = c(1, -1))),
    not_a_list     = list(sigma = "escape please"),
    unknown_dpar   = list(nosuch = list(at = 0, tol = 1e6, from = c(1, -1))),
    # ATOMIC declarations: `[[` on these raises "subscript out of
    # bounds" if the shape of the field is not checked before it is
    # indexed, which killed the fit outright rather than being ignored
    atomic_chr     = "alpha",
    atomic_num     = c(0, 1),
    atomic_named   = c(alpha = 0),
    atomic_lgl     = TRUE,
    atomic_na      = NA
  )
  for (nm in names(bad)) {
    fit <- frm(fo, family = sn_norm_family(bad[[nm]]), data = dd)
    expect_null(fit$opt[["stationary_escape"]], info = nm)
    expect_identical(fit$opt$objective, ref$opt$objective, info = nm)
  }
})

test_that("a dpar design with no intercept gets the start and the escape", {
  # `alpha ~ 0 + g` spans the same columns as `alpha ~ g`, so it has the
  # same optimum; writing the start into an "(Intercept)" that is not
  # there left every coefficient on the stationary point instead.
  set.seed(31)
  n <- 300
  g <- factor(rep(1:4, length.out = n))
  x <- rnorm(n)
  dd <- data.frame(
    x = x, g = g,
    y = 0.4 + 0.8 * x +
      1.2 * (0.9486833 * abs(rnorm(n)) + 0.3162278 * rnorm(n)))
  with_i <- frm(bf(y ~ x, sigma ~ 1, alpha ~ g), family = skew_normal(),
                data = dd)
  no_i <- frm(bf(y ~ x, sigma ~ 1, alpha ~ 0 + g), family = skew_normal(),
              data = dd)
  ll_i <- as.numeric(logLik(with_i))
  # same column space, so the same maximum, to the optimizer's tolerance
  expect_lt(abs(as.numeric(logLik(no_i)) - ll_i), abs(ll_i) * 1e-6)
  # and it is nowhere near the stationary point
  expect_gt(max(abs(no_i$estimates$betad[grep("^alpha",
                                              names(no_i$estimates$betad))])),
            1)
})

test_that("an intercept design places the start in the intercept alone", {
  # predictor_at() must return the exact intercept answer, or every
  # model with an intercept moves in its last bits for no reason
  X <- cbind(`(Intercept)` = rep(1, 5), x = c(-2, -1, 0, 1, 2))
  expect_identical(frmtmb:::predictor_at(X, 2.5), c(2.5, 0))
  # without an intercept the value spreads over the columns that can
  # carry it: cell means all take it
  G <- cbind(a = c(1, 1, 0, 0), b = c(0, 0, 1, 1))
  expect_equal(as.vector(G %*% frmtmb:::predictor_at(G, 1.5)),
               rep(1.5, 4), tolerance = 1e-10)
})

test_that("offset() is taken out of the residual the start reads", {
  # offset() is part of the mu predictor but is not in X, so a model
  # that enters its covariate as an offset used to keep the covariate's
  # skew in the "residual" and start alpha from the wrong side
  set.seed(1)
  n <- 200
  o <- -abs(rnorm(n)) * 3
  y <- o + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y, o = o)
  expect_lt(sn_skew(dd$y) * sn_skew(dd$y - dd$o), 0)   # the trap
  u <- frm(bf(y ~ offset(o), sigma ~ 1, alpha ~ 1), family = skew_normal(),
           data = dd, dry_run = "objective")
  a0 <- u$estimates$betad[["alpha_(Intercept)"]]
  expect_gt(a0 * sn_skew(dd$y - dd$o), 0)
  # and the fit needs no refit to reach the optimum
  fit <- frm(bf(y ~ offset(o), sigma ~ 1, alpha ~ 1), family = skew_normal(),
             data = dd)
  expect_null(fit$opt[["stationary_escape"]])
  ref <- max(vapply(c(2, -2), function(a) {
    as.numeric(logLik(frm(bf(y ~ offset(o), sigma ~ 1, alpha ~ 1),
                          family = skew_normal(), data = dd,
                          start = list(betad = c(0, a)))))
  }, 0))
  expect_lt(ref - as.numeric(logLik(fit)), abs(ref) * 1e-8)
})

test_that("a sparse() design still gets its initializer", {
  # predictor_at() judges a design by dim(), not by is.matrix(): a
  # sparse() predictor is a Matrix, and an is.matrix() gate dropped the
  # initializer of every sparse model silently. Caught by test-sparsex.R
  # rather than here, which is why it is pinned here too.
  set.seed(12)
  n <- 200
  dd <- data.frame(y = rnorm(n, 5, 2), g = factor(rep(1:8, length.out = n)))
  dense <- frm(bf(y ~ g), data = dd, dry_run = "objective")
  sparse <- frm(bf(y ~ g), data = dd,
                control = frmtmb_control(sparse_x = TRUE),
                dry_run = "objective")
  expect_equal(unname(dense$estimates$beta[1]),
               unname(sparse$estimates$beta[1]), tolerance = 1e-12)
  expect_equal(unname(dense$estimates$beta[1]), mean(dd$y),
               tolerance = 1e-12)
  # and the predictor helpers take a Matrix without converting callers
  X <- Matrix::sparse.model.matrix(~ g, dd)
  expect_identical(frmtmb:::predictor_at(X, 1.25)[1], 1.25)
  expect_true(all(frmtmb:::predictor_at(X, 1.25)[-1] == 0))
})

test_that("a dpar declaring nothing is skipped without an intercept", {
  # The start placement is confined to a dpar that DECLARES a
  # stationary point. Spreading a least-squares start over every
  # family's intercept-less design is right on the predictor scale and
  # ruinous on the parameter scale: on `y ~ 0 + xt` with xt at 1e-6 it
  # starts the coefficient near 1.6e6, where nlminb's relative step
  # test fires before the fit moves. glm() is the oracle.
  set.seed(23)
  n <- 250
  xt <- rnorm(n)
  dd <- data.frame(y = rpois(n, pmin(exp(0.5 + 0.4 * xt), 1e6)))
  for (s in c(1, 1e-2, 1e-4, 1e-6)) {
    dd$xs <- xt * s
    fit <- frm(bf(y ~ 0 + xs), family = poisson(), data = dd)
    ll <- as.numeric(logLik(fit))
    # The reference is the ZERO start, which is what this design got
    # before the placement existed: the assertion is that the start
    # cannot make the fit worse. Pinning glm() at every scale would
    # assert away a pre-existing weakness at 1e-6, where frmtmb sits
    # below glm on the unchanged build too (dev/test-backlog.md).
    zero <- as.numeric(logLik(frm(bf(y ~ 0 + xs), family = poisson(),
                                  data = dd, start = list(beta = 0))))
    expect_gt(ll, zero - abs(zero) * 1e-8, label = paste("scale", s))
    if (s >= 1e-4) {
      rl <- as.numeric(stats::logLik(
        stats::glm(y ~ 0 + xs, family = stats::poisson(), data = dd)))
      expect_lt(abs(ll - rl), abs(rl) * 1e-8, label = paste("glm", s))
    }
  }
  # and the start itself must be untouched: a poisson mu declares no
  # stationary point, so `~ 0 + g` keeps every coefficient at zero
  set.seed(24)
  d2 <- data.frame(y = rpois(100, 4),
                   g = factor(rep(1:4, length.out = 100)))
  u <- frm(bf(y ~ 0 + g), family = poisson(), data = d2,
           dry_run = "objective")
  expect_true(all(u$estimates$beta == 0))
})

test_that("a declaring dpar IS placed across an intercept-less design", {
  # the complement of the test above, so the confinement cannot pass
  # by disabling the feature outright
  set.seed(23)
  n <- 250
  x <- rnorm(n)
  dd <- data.frame(
    x = x, g = factor(rep(1:4, length.out = n)),
    y = 0.4 + 0.8 * x +
      1.2 * (0.9486833 * abs(rnorm(n)) + 0.3162278 * rnorm(n)))
  u <- frm(bf(y ~ x, sigma ~ 1, alpha ~ 0 + g), family = skew_normal(),
           data = dd, dry_run = "objective")
  a <- u$estimates$betad[grep("^alpha", names(u$estimates$betad))]
  expect_length(a, 4L)
  expect_true(all(abs(a) > 1))        # every cell carries the value
  expect_lt(diff(range(a)), 1e-8)     # and they carry the same one
  # sigma declares nothing, so its intercept-less design stays at zero
  u2 <- frm(bf(y ~ x, sigma ~ 0 + g, alpha ~ 1), family = skew_normal(),
            data = dd, dry_run = "objective")
  s <- u2$estimates$betad[grep("^sigma", names(u2$estimates$betad))]
  expect_true(all(s == 0))
})

test_that("a non-finite `from` is dropped at the guard, not downstream", {
  set.seed(6)
  dd <- data.frame(y = rnorm(150, 3, 2), x = rnorm(150))
  fo <- bf(y ~ x, sigma ~ 1)
  one <- frm(fo, data = dd,
             family = sn_norm_family(list(sigma = list(at = 0, tol = 1e6,
                                                       from = c(Inf, -1)))))
  esc <- one$opt[["stationary_escape"]]
  expect_false(is.null(esc))
  expect_identical(unname(esc[["starts"]]), 1)   # Inf dropped, -1 kept
  none <- frm(fo, data = dd,
              family = sn_norm_family(list(sigma = list(at = 0, tol = 1e6,
                                                        from = c(Inf, NaN)))))
  expect_null(none$opt[["stationary_escape"]])
})
