# fitted() on a new row is the mean over the row's draws: a me() value, a
# latent residual, and the effects of a grouping level the fit did not see,
# as brms's posterior_epred() draws each of them. Each test was seen to fail
# on 0.79.1 (dev/readers-findings.md has the output).

dm_lp <- function(f, nd, dp) {
  as.numeric(frm_linpred(f, newdata = nd, dpar = dp, type = "link"))
}

test_that("the product Gauss-Hermite rule integrates", {
  for (s in c(1, 3)) {
    m <- me_gh_moments(function(z) exp(s * z[[1L]]), 1L, 41L)
    expect_lt(abs(m$m1 / exp(s^2 / 2) - 1), 1e-8)
  }
  m2 <- me_gh_moments(function(z) exp(0.7 * z[[1L]] - 0.4 * z[[2L]]), 2L, 21L)
  expect_lt(abs(m2$m1 / exp((0.49 + 0.16) / 2) - 1), 1e-10)
})

test_that("a rule that does not settle past two directions says so", {
  # exp(s z1 + 0.5 (z2 + z3)) on 11 nodes a direction: its sd is 7 percent
  # low at s = 2.5 and 3e-7 at s = 1
  for (s in c(1, 2.5)) {
    val <- function(z) exp(s * z[[1L]] + 0.5 * (z[[2L]] + z[[3L]]))
    mo <- me_gh_moments(val, 3L, 11L)
    w <- tryCatch(me_rule_check(val, 3L, 11L, mo, "y"),
                  warning = function(w) w)
    if (s == 1) expect_null(w) else expect_s3_class(w, "frmtmb_draw_rule")
  }
})

test_that("the per-row Cholesky factor is chol() row by row", {
  # the cross term is live from the third direction; dropping it moved a
  # three-direction mean by 2.9e-3 of itself, inside Monte Carlo noise
  set.seed(3301)
  n <- 6
  D <- 3
  S <- lapply(seq_len(n), function(i) {
    A <- matrix(stats::rnorm(D * D), D)
    crossprod(A) + diag(0.1, D)
  })
  C <- lapply(seq_len(D), function(i) lapply(seq_len(i), function(j) {
    vapply(S, function(s) s[i, j], 0)
  }))
  L <- chol_rows(C, D)
  for (k in seq_len(n)) {
    R <- t(chol(S[[k]]))
    for (i in seq_len(D)) for (j in seq_len(i)) {
      expect_lt(abs(L[[i]][[j]][k] - R[i, j]), 1e-12 * max(abs(R)))
    }
  }
})

test_that("a new level's mean is over its effects, not at their mean", {
  set.seed(4101)
  J <- 40
  n <- 400
  dd <- data.frame(x = stats::runif(n, -1, 1),
                   g = factor(sample(seq_len(J), n, TRUE)))
  u <- stats::rnorm(J)
  dd$y <- stats::rpois(n, exp(0.5 + 0.3 * dd$x + u[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd)
  s <- VarCorr(fit)$g$sd[1, 1]
  nd <- data.frame(x = c(-0.5, 0, 0.8), g = "new")
  eta <- frm_lp_basis(fit, newdata = nd, allow_new_levels = TRUE)$eta
  f <- fitted(fit, newdata = nd, allow_new_levels = TRUE)
  # 0.79.1 read exp(eta), low by exp(s^2 / 2) = 1.647 here
  expect_lt(max(abs(f[, "Estimate"] / exp(eta + s^2 / 2) - 1)), 1e-6)
  # the effects' own spread, lognormal's sd, is in Est.Error
  sd_u <- sqrt((exp(s^2) - 1) * exp(2 * eta + s^2))
  expect_true(all(f[, "Est.Error"] >= sd_u))
  expect_lt(max(f[, "Est.Error"] / sd_u), 1.1)
  expect_true(all(f[, "Q2.5"] < f[, "Estimate"] &
                    f[, "Estimate"] < f[, "Q97.5"]))
  # a level the fit saw, re_formula = NA, and a seen level drawn for the new
  # one take no mean over effects
  seen <- transform(nd, g = dd$g[1])
  expect_equal(unname(fitted(fit, newdata = seen)[, "Estimate"]),
               unname(as.numeric(frm_linpred(fit, newdata = seen,
                                             type = "response"))),
               tolerance = 1e-12)
  expect_equal(unname(fitted(fit, newdata = nd, re_formula = NA,
                             allow_new_levels = TRUE)[, "Estimate"]),
               unname(exp(eta)), tolerance = 1e-12)
  set.seed(1)
  old <- fitted(fit, newdata = nd, allow_new_levels = TRUE,
                sample_new_levels = "old_levels")[, "Estimate"]
  lev_vals <- vapply(levels(dd$g), function(l) {
    as.numeric(frm_linpred(fit, newdata = transform(nd[1, ], g = l),
                           type = "response"))
  }, 0)
  expect_true(any(abs(lev_vals / old[1] - 1) < 1e-10))
})

test_that("a new level of a correlated slope reads z' Sigma z", {
  set.seed(4102)
  J <- 40
  n <- 400
  dd <- data.frame(x = stats::runif(n, -1, 1),
                   g = factor(sample(seq_len(J), n, TRUE)))
  L <- t(chol(matrix(c(0.64, 0.2, 0.2, 0.36), 2)))
  u <- t(L %*% matrix(stats::rnorm(2 * J), 2))
  dd$y <- stats::rpois(n, exp(0.5 + 0.3 * dd$x + u[dd$g, 1] +
                                u[dd$g, 2] * dd$x))
  fit <- frm(bf(y ~ x + (1 + x | g)), family = poisson(), data = dd)
  nd <- data.frame(x = c(-0.5, 0, 0.8), g = "new")
  # extra_var is the new level's z' Sigma z on each row
  lb <- frm_lp_basis(fit, newdata = nd, allow_new_levels = TRUE)
  vc <- VarCorr(fit)$g
  S <- diag(vc$sd[, 1]) %*% vc$cor[, 1, ] %*% diag(vc$sd[, 1])
  Z <- cbind(1, nd$x)
  expect_lt(max(abs(lb$extra_var / rowSums((Z %*% S) * Z) - 1)), 1e-8)
  f <- fitted(fit, newdata = nd, allow_new_levels = TRUE)
  expect_lt(max(abs(f[, "Estimate"] / exp(lb$eta + lb$extra_var / 2) - 1)),
            1e-6)
})

test_that("a new level's mean on a logit link is the integral over it", {
  set.seed(4103)
  J <- 40
  n <- 400
  dd <- data.frame(x = stats::runif(n, -1, 1),
                   g = factor(sample(seq_len(J), n, TRUE)))
  u <- stats::rnorm(J, 0, 1.5)
  dd$y <- stats::rbinom(n, 1, stats::plogis(0.2 + 0.5 * dd$x + u[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = dd)
  s <- VarCorr(fit)$g$sd[1, 1]
  nd <- data.frame(x = c(-0.5, 0.8), g = "new")
  eta <- frm_lp_basis(fit, newdata = nd, allow_new_levels = TRUE)$eta
  ref <- vapply(eta, function(e) {
    stats::integrate(function(v) stats::plogis(e + v) * stats::dnorm(v, 0, s),
                     -Inf, Inf, rel.tol = 1e-10)$value
  }, 0)
  ref2 <- vapply(eta, function(e) {
    stats::integrate(function(v) {
      stats::plogis(e + v)^2 * stats::dnorm(v, 0, s)
    }, -Inf, Inf, rel.tol = 1e-10)$value
  }, 0)
  f <- fitted(fit, newdata = nd, allow_new_levels = TRUE)
  expect_lt(max(abs(f[, "Estimate"] / ref - 1)), 1e-6)
  sd_u <- sqrt(ref2 - ref^2)
  expect_true(all(f[, "Est.Error"] >= sd_u))
  expect_lt(max(f[, "Est.Error"] / sd_u), 1.2)
})

dm_latent_data <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = 1:25, s = factor(1:8))
  d$x <- stats::rnorm(nrow(d))
  d$zo <- stats::rnorm(nrow(d))
  d$sz <- 0.5
  e <- stats::ave(stats::rnorm(nrow(d)), d$s, FUN = function(z) {
    stats::filter(z, 0.5, "recursive")
  })
  d$e <- 0.4 * e
  d$zl <- d$zo + stats::rnorm(nrow(d), 0, 0.5)
  d
}

test_that("a latent and one me() on two predictors are two draws", {
  d <- dm_latent_data(9201)
  lam <- exp(0.4 + 0.3 * d$x + 0.4 * d$zl + d$e)
  cnt <- stats::rpois(nrow(d), lam)
  while (any(cnt == 0)) {
    i <- cnt == 0
    cnt[i] <- stats::rpois(sum(i), lam[i])
  }
  d$y <- ifelse(stats::runif(nrow(d)) < stats::plogis(-0.5 + 0.8 * d$zl), 0L,
                cnt)
  fit <- suppressWarnings(frm(
    bf(y ~ x + me(zo, sz) + ar(time = t, gr = s, cov = TRUE),
       hu ~ me(zo, sz)), data = d, family = hurdle_poisson()))
  nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1), zo = c(0, 1.5),
                   sz = 0.5)
  sl <- latres_row_sd(fit, "y", NULL)
  em <- dm_lp(fit, nd, "mu")
  eh <- dm_lp(fit, nd, "hu")
  bm <- dm_lp(fit, transform(nd, zo = zo + 1), "mu") - em
  bh <- dm_lp(fit, transform(nd, zo = zo + 1), "hu") - eh
  mean_at <- function(i, p) {
    inner <- function(z) {
      vapply(z, function(zz) {
        stats::integrate(function(u) {
          l <- exp(em[i] + bm[i] * 0.5 * zz + sl * u)
          ((1 - stats::plogis(eh[i] + bh[i] * 0.5 * zz)) * l /
             (1 - exp(-l)))^p * stats::dnorm(u)
        }, -10, 10, rel.tol = 1e-10)$value
      }, 0) * stats::dnorm(z)
    }
    stats::integrate(inner, -10, 10, rel.tol = 1e-9)$value
  }
  ref <- vapply(1:2, mean_at, 0, p = 1)
  ref2 <- vapply(1:2, mean_at, 0, p = 2)
  mm <- fitted_me_moments(fit, nd, NULL, NULL, NULL, FALSE)
  # 0.79.1 moved hu by its regression on mu: the draw's sd read 0.733
  # against 0.925 on row 2
  expect_lt(max(abs(mm$mean / ref - 1)), 1e-5)
  expect_lt(max(abs(sqrt(mm$var) / sqrt(ref2 - ref^2) - 1)), 1e-4)
  f <- fitted(fit, newdata = nd)
  expect_equal(unname(f[, "Estimate"]), mm$mean, tolerance = 1e-12)
  expect_true(all(f[, "Q2.5"] < f[, "Estimate"] &
                    f[, "Estimate"] < f[, "Q97.5"]))
})

test_that("a lognormal whose sigma reads a draw reads the quadrature mean", {
  # the exact mean is infinite from the draw's far tail; 0.79.1's lane rule
  # read Inf, and a mean of draws reads a finite value at these sizes
  set.seed(7102)
  n <- 800
  J <- 30
  d <- data.frame(x = stats::rnorm(n), g = factor(sample(J, n, TRUE)))
  d$y <- exp(0.5 + 0.4 * d$x + stats::rnorm(n, 0, exp(
    -0.5 + stats::rnorm(J, 0, 0.3)[d$g])))
  fit <- frm(bf(y ~ x, sigma ~ 1 + (1 | g)), family = lognormal(), data = d)
  nd <- data.frame(x = 0, g = "new")
  f <- allow_warnings(fitted(fit, newdata = nd, allow_new_levels = TRUE),
                      "is infinite", require = "is infinite")
  cnd <- tryCatch(fitted(fit, newdata = nd, allow_new_levels = TRUE),
                  warning = function(w) w)
  expect_s3_class(cnd, "frmtmb_mean_not_finite")
  s <- VarCorr(fit)$g$sd[1, 1]
  em <- as.numeric(frm_linpred(fit, newdata = nd, re_formula = NA))
  es <- as.numeric(frm_linpred(fit, newdata = nd, re_formula = NA,
                               dpar = "sigma"))
  set.seed(1)
  mc <- mean(exp(em + exp(2 * (es + s * stats::rnorm(4e6))) / 2))
  expect_lt(abs(f[1, "Estimate"] / mc - 1), 1e-3)
  expect_true(all(is.finite(f)))
  # a sweep of the draw's size on sigma: finite, said, and a mean of draws
  set.seed(7101)
  d <- data.frame(x = stats::rnorm(n), l = stats::rnorm(n), sz = 0.3)
  d$zo <- d$l + stats::rnorm(n, 0, 0.3)
  d$y <- exp(0.5 + 0.4 * d$x + stats::rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
  f1 <- frm(bf(y ~ x, sigma ~ me(zo, sz)), family = lognormal(), data = d)
  for (sz in c(0.01, 0.1, 0.2, 0.3)) {
    nd <- data.frame(x = 0, zo = 1, sz = sz)
    v <- allow_warnings(fitted(f1, newdata = nd), "is infinite",
                        require = "is infinite")
    em <- dm_lp(f1, nd, "mu")
    es <- dm_lp(f1, nd, "sigma")
    bs <- dm_lp(f1, transform(nd, zo = zo + 1), "sigma") - es
    set.seed(2)
    u <- bs * sz * stats::rnorm(1e6)
    mc <- exp(em + exp(2 * (es + u)) / 2)
    expect_lt(abs(v[1, "Estimate"] / mean(mc) - 1), 2e-3, label = sz)
    mm <- allow_warnings(fitted_me_moments(f1, nd, NULL, NULL, NULL, FALSE),
                         "is infinite")
    expect_lt(abs(sqrt(mm$var) / stats::sd(mc) - 1), 0.05, label = sz)
  }
  # me() on mu alone: the mean is finite and is the closed form
  set.seed(9203)
  n <- 600
  d <- data.frame(zo = stats::rnorm(n), sz = 0.5)
  zl <- d$zo + stats::rnorm(n, 0, 0.5)
  d$y <- exp(0.5 + 0.4 * zl + stats::rnorm(n, 0, 0.5))
  nd <- data.frame(zo = c(0, 1.5), sz = 0.5)
  fit2 <- frm(bf(y ~ me(zo, sz)), data = d, family = lognormal())
  em <- dm_lp(fit2, nd, "mu")
  b <- dm_lp(fit2, transform(nd, zo = zo + 1), "mu") - em
  sig <- as.numeric(frm_linpred(fit2, newdata = nd, dpar = "sigma",
                                type = "response"))
  f2 <- fitted(fit2, newdata = nd)
  expect_lt(max(abs(f2[, "Estimate"] /
                      exp(em + (b * 0.5)^2 / 2 + sig^2 / 2) - 1)), 1e-12)
})

test_that("a lognormal component's or a body's sigma reads the same way", {
  # the rule keyed on a dpar named sigma missed both, and read 2.5e207 and
  # 4.7e82 with no warning
  n <- 800
  set.seed(7201)
  d <- data.frame(x = stats::rnorm(n), l = stats::rnorm(n), sz = 0.3)
  d$zo <- d$l + stats::rnorm(n, 0, 0.3)
  comp <- stats::runif(n) < 0.4
  d$y <- ifelse(comp, exp(1.5 + stats::rnorm(n, 0, exp(-0.7 + 0.5 * d$l))),
                exp(-0.5 + stats::rnorm(n, 0, 0.4)))
  fm <- frm(bf(y ~ 1, sigma1 ~ me(zo, sz)),
            family = mixture(lognormal(), lognormal()), data = d)
  nd <- data.frame(zo = 1, sz = c(3, 10))
  f <- allow_warnings(fitted(fm, newdata = nd), "is infinite",
                      require = "is infinite")
  expect_true(all(is.finite(f)))
  # the brms-style mean at the small draw: component 1's integral, at the
  # estimates, plus component 2's mean
  em <- dm_lp(fm, nd[1, ], "mu1")
  es <- dm_lp(fm, nd[1, ], "sigma1")
  bs <- dm_lp(fm, transform(nd[1, ], zo = zo + 1), "sigma1") - es
  th1 <- as.numeric(frm_linpred(fm, newdata = nd[1, ], dpar = "theta1",
                                type = "response"))
  rest <- (1 - th1) * exp(dm_lp(fm, nd[1, ], "mu2") + as.numeric(
    frm_linpred(fm, newdata = nd[1, ], dpar = "sigma2", type = "response"))^2 /
      2)
  set.seed(3)
  mc <- th1 * mean(exp(em + exp(2 * (es + bs * 3 * stats::rnorm(1e6))) / 2)) +
    rest
  expect_lt(abs(f[1, "Estimate"] / mc - 1), 1e-3)
  # sigma reached through a body
  set.seed(7202)
  d <- data.frame(x = stats::rnorm(n), l = stats::rnorm(n), sz = 0.3)
  d$zo <- d$l + stats::rnorm(n, 0, 0.3)
  d$y <- exp(0.5 + 0.4 * d$x + stats::rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
  fs <- frm(bf(y ~ x) + nlf(sigma ~ exp(s1)) + lf(s1 ~ me(zo, sz)),
            family = lognormal(), data = d)
  nd <- data.frame(x = 0, zo = 1, sz = c(0.3, 0.8))
  w <- tryCatch(fitted(fs, newdata = nd), warning = function(w) w)
  expect_s3_class(w, "frmtmb_mean_not_finite")
  f <- allow_warnings(fitted(fs, newdata = nd), "is infinite")
  # at the data's own sz the value is finite and near the plug-in; at 0.8,
  # sigma = exp(exp(s1)) overflows inside 5.19 sds, and the row reads NaN
  # under the warning, not a number
  expect_true(all(is.finite(f[1, ])))
  expect_lt(f[1, "Estimate"], 10 * exp(dm_lp(fs, nd[1, ], "mu")))
  expect_true(is.finite(f[2, "Estimate"]) || is.nan(f[2, "Estimate"]))
})

test_that("the warning on a sigma draw is said once per call", {
  set.seed(7101)
  n <- 800
  d <- data.frame(x = stats::rnorm(n), l = stats::rnorm(n), sz = 0.3)
  d$zo <- d$l + stats::rnorm(n, 0, 0.3)
  d$y <- exp(0.5 + 0.4 * d$x + stats::rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
  fit <- frm(bf(y ~ x, sigma ~ me(zo, sz)), family = lognormal(), data = d)
  seen <- character()
  withCallingHandlers(
    conditional_effects(fit, effects = "x",
                        int_conditions = list(x = c(-1, 0, 1)),
                        conditions = data.frame(zo = c(0, 1), sz = 0.3)),
    frmtmb_mean_not_finite = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(seen, 1L)
  expect_match(seen, "move that mean by at most")
})

test_that("a draw whose unit step lands on a period is kept", {
  # y ~ sin(a)^2 + c at a = 0 and a draw sd of pi: steps of one sd read
  # sin(pi)^2 = sin(0)^2, and the draw, which doubles the mean, was dropped
  set.seed(7208)
  n <- 800
  d <- data.frame(l1 = stats::rnorm(n), l2 = stats::rnorm(n), s1 = 0.2,
                  s2 = 0.2)
  d$x1 <- d$l1 + stats::rnorm(n, 0, 0.2)
  d$x2 <- d$l2 + stats::rnorm(n, 0, 0.2)
  d$y <- sin(0.3 + 0.8 * d$l1)^2 + (0.5 + 0.5 * d$l2) +
    stats::rnorm(n, 0, 0.1)
  fit <- frm(bf(y ~ sin(a)^2 + c, a ~ me(x1, s1), c ~ me(x2, s2), nl = TRUE),
             family = gaussian(), data = d)
  q <- data.frame(x1 = 0, x2 = 0, s1 = 1, s2 = 1)
  ba <- dm_lp(fit, transform(q, x1 = 1), "a") - dm_lp(fit, q, "a")
  a0 <- dm_lp(fit, q, "a")
  # at 10 pi, steps of 0.5, 1 and 1.7 sds all land on zeros of the sine
  for (s in c(pi, 10 * pi)) {
    nd <- data.frame(x1 = -a0 / ba, x2 = 0, s1 = s / abs(ba), s2 = 0.5)
    loads <- me_row_loads(fit, fit$spec$responses[[1L]], nd, TRUE, FALSE)
    val <- me_node_value(fit, loads, function(f) {
      as.matrix(frm_linpred(f, newdata = nd, type = "response"))
    })
    kept <- me_keep_dirs(loads, val, val(0), "y")
    expect_length(kept$dirs, 2L)
    # and fitted() reads the mean, or says the rule cannot
    said <- FALSE
    f <- withCallingHandlers(fitted(fit, newdata = nd),
                             frmtmb_mean_unsettled = function(w) {
                               said <<- TRUE
                               invokeRestart("muffleWarning")
                             })
    exact <- (1 - exp(-2 * s^2)) / 2 + dm_lp(fit, nd, "c")
    expect_true(said || abs(f[, "Estimate"] / exact - 1) < 1e-3, label = s)
  }
  # a row that does not settle leaves the row beside it on its own rule
  nd <- data.frame(x1 = -a0 / ba, x2 = 0, s1 = c(pi, 10 * pi) / abs(ba),
                   s2 = 0.5)
  f <- allow_warnings(fitted(fit, newdata = nd), "does not settle")
  exact <- (1 - exp(-2 * pi^2)) / 2 + dm_lp(fit, nd[1L, ], "c")
  expect_lt(abs(f[1L, "Estimate"] / exact - 1), 1e-3)
})

test_that("a log-linear mean is finite at any draw size", {
  # the first round of this rule read Inf from a link sd of 6.6
  set.seed(7103)
  n <- 800
  d <- data.frame(l = stats::rnorm(n), sx = 0.3)
  d$xo <- d$l + stats::rnorm(n, 0, 0.3)
  d$y <- stats::rpois(n, exp(0.2 + d$l))
  fit <- frm(bf(y ~ me(xo, sx)), family = poisson(), data = d)
  b <- dm_lp(fit, data.frame(xo = 1, sx = 1), "mu") -
    dm_lp(fit, data.frame(xo = 0, sx = 1), "mu")
  for (s in c(1, 4, 6.6, 8)) {
    nd <- data.frame(xo = 0, sx = s / b)
    e <- dm_lp(fit, nd, "mu")
    f <- fitted(fit, newdata = nd)
    expect_lt(abs(f[, "Estimate"] / exp(e + s^2 / 2) - 1), 1e-10, label = s)
    sdx <- sqrt(expm1(s^2) * exp(2 * e + s^2))
    expect_gte(f[, "Est.Error"] / sdx, 1 - 1e-10, label = s)
    expect_lt(f[, "Est.Error"] / sdx, 1.01, label = s)
  }
})

test_that("se__ of conditional_effects() is on the displayed scale", {
  set.seed(7971)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  eta <- 0.4 + 0.7 * d$x
  d$cnt <- stats::rpois(n, exp(eta))
  d$pos <- exp(eta + stats::rnorm(n, 0, exp(-0.5 + 0.3 * d$x)))
  nd <- data.frame(x = c(-1, 0, 1))
  for (sp in list(list(bf(cnt ~ x), poisson(), NULL),
                  list(bf(pos ~ x, sigma ~ x), lognormal(), NULL),
                  list(bf(pos ~ x, sigma ~ x), lognormal(), "sigma"))) {
    fit <- frm(sp[[1]], family = sp[[2]], data = d)
    ce <- conditional_effects(fit, effects = "x", dpar = sp[[3]],
                              int_conditions = list(x = nd$x))[[1L]]
    # brms's se__ is the MAD of its draws of the displayed value; 0.79.1
    # read the link scale, 0.027 where the response scale reads 0.045
    want <- frm_linpred(fit, newdata = nd, dpar = sp[[3]], type = "response",
                        se.fit = TRUE)$se.fit
    expect_lt(max(abs(ce$se__ / want - 1)), 1e-8,
              label = paste(sp[[2]]$family, sp[[3]]))
    if (is.null(sp[[3]])) {
      expect_lt(max(abs(ce$se__ / fitted(fit, newdata = nd)[, "Est.Error"] -
                          1)), 1e-8)
    }
  }
})

test_that("a nonlinear body's Est.Error is the delta method through it", {
  set.seed(8101)
  d <- data.frame(x = stats::runif(200, 0, 3))
  d$y <- 3 * (1 - exp(-0.8 * d$x)) + stats::rnorm(200, 0, 0.2)
  fit <- suppressMessages(frm(
    bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE),
    family = gaussian(), data = d,
    prior = c(prior(normal(3, 2), nlpar = "b1"),
              prior(normal(1, 1), nlpar = "b2"))))
  nd <- data.frame(x = c(0.5, 2))
  f <- fitted(fit, newdata = nd)
  # 0.79.1 returned NA here with no message
  expect_true(all(is.finite(f[, "Est.Error"])))
  b <- fixef(fit)[, "Estimate"]
  G <- cbind(1 - exp(-b[2] * nd$x), b[1] * nd$x * exp(-b[2] * nd$x))
  V <- vcov(fit)[1:2, 1:2]
  expect_lt(max(abs(f[, "Est.Error"] / sqrt(rowSums((G %*% V) * G)) - 1)),
            1e-6)
  # against the spread of the body over the coefficients' own law
  set.seed(11)
  B <- MASS::mvrnorm(20000, b, V)
  sims <- vapply(nd$x, function(x) B[, 1] * (1 - exp(-B[, 2] * x)),
                 numeric(20000))
  ratio <- f[, "Est.Error"] / apply(sims, 2, stats::sd)
  expect_lt(max(abs(ratio - 1)), 4 / sqrt(2 * 20000) + 0.02)
})

test_that("conditional_effects() reads a body whose parameter reads me()", {
  set.seed(5101)
  d <- data.frame(xo = stats::rnorm(200), sx = 0.3, z = stats::runif(200))
  d$y <- 2 * exp(0.5 * (d$xo + stats::rnorm(200, 0, 0.3))) * d$z +
    stats::rnorm(200, 0, 0.5)
  fit <- frm(bf(y ~ a * exp(b * z), a ~ me(xo, sx), b ~ 1, nl = TRUE),
             family = gaussian(), data = d)
  # 0.79.1 stopped: "Its chain rule through a nonlinear body has not been
  # measured"
  ce <- conditional_effects(fit, effects = "z",
                            int_conditions = list(z = c(0.2, 0.8)))[[1L]]
  expect_true(all(is.finite(ce$se__)))
  nd <- data.frame(xo = ce$xo, sx = ce$sx, z = ce$z)
  lb <- frm_lp_basis(fit, newdata = nd)
  bb <- dm_lp(fit, nd, "b")
  ba <- dm_lp(fit, transform(nd, xo = xo + 1), "a") - dm_lp(fit, nd, "a")
  # a enters linearly, times exp(b z): its draw's variance at first order
  expect_lt(max(abs(lb$extra_var / (exp(bb * nd$z) * ba * nd$sx)^2 - 1)),
            1e-8)
  expect_error(frm_lp_basis(fit, newdata = nd, extra_cov = TRUE),
               "across rows")
})

test_that("a body's band under me() is fitted()'s, and a draw's own", {
  # the first order read se__ 1.82 where the draw's sd is 3.93, and a band
  # below 0 for a positive mean
  set.seed(8202)
  n <- 600
  d <- data.frame(l = stats::rnorm(n), sx = 0.3)
  d$xo <- d$l + stats::rnorm(n, 0, 0.3)
  d$y <- exp(0.2 + 0.8 * d$l) + stats::rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ exp(a), a ~ me(xo, sx), nl = TRUE), family = gaussian(),
             data = d)
  ba <- dm_lp(fit, data.frame(xo = 1, sx = 1), "a") -
    dm_lp(fit, data.frame(xo = 0, sx = 1), "a")
  set.seed(41)
  z <- stats::rnorm(1e6)
  for (s in c(0.5, 1)) {
    nd <- data.frame(xo = 0.5, sx = s / ba)
    f <- fitted(fit, newdata = nd)
    ce <- conditional_effects(fit, effects = "xo",
                              int_conditions = list(xo = 0.5),
                              conditions = data.frame(sx = s / ba))[[1L]]
    expect_equal(ce$se__, unname(f[, "Est.Error"]), tolerance = 1e-10)
    expect_equal(ce$lower__, unname(f[, "Q2.5"]), tolerance = 1e-10)
    expect_equal(ce$upper__, unname(f[, "Q97.5"]), tolerance = 1e-10)
    mc <- exp(dm_lp(fit, nd, "a") + s * z)
    q <- stats::quantile(mc, c(0.025, 0.975))
    expect_lt(abs(f[, "Estimate"] / mean(mc) - 1), 4 * stats::sd(mc) /
                sqrt(length(mc)) / mean(mc))
    expect_lt(abs(ce$se__ / stats::sd(mc) - 1), 0.02)
    expect_lt(max(abs(c(ce$lower__, ce$upper__) / q - 1)), 0.02)
  }
})

test_that("a body under allow_new_levels answers, or says why not", {
  # Est.Error was NA with no message under the flag, with no group term too
  set.seed(8201)
  n <- 300
  d <- data.frame(x = stats::runif(n, 0, 3), g = factor(rep(1:12, 25)))
  d$y <- (3 + stats::rnorm(12, 0, 0.3)[d$g]) * exp(-0.8 * d$x) +
    stats::rnorm(n, 0, 0.1)
  fit <- suppressMessages(frm(
    bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
    family = gaussian(), data = d,
    prior = c(prior(normal(3, 2), nlpar = "a"),
              prior(normal(1, 1), nlpar = "b"))))
  seen <- data.frame(x = c(0.5, 2), g = "1")
  expect_equal(fitted(fit, newdata = seen, allow_new_levels = TRUE),
               fitted(fit, newdata = seen))
  new <- transform(seen, g = "new")
  f <- fitted(fit, newdata = new, allow_new_levels = TRUE)
  expect_true(all(is.finite(f[, "Est.Error"])))
  expect_true(all(f[, "Est.Error"] > fitted(fit, newdata = seen)[,
                                                                 "Est.Error"]))
  expect_error(frm_linpred(fit, newdata = new, allow_new_levels = TRUE,
                           se.fit = TRUE), class = "frmtmb_nl_new_level")
})

test_that("past four draws the ones that move the mean are kept", {
  # five me() predictors of a two-gaussian mixture: by dpar order the two
  # sigmas, which the mean does not read, displaced theta2
  set.seed(6202)
  n <- 3000
  d <- data.frame(matrix(stats::rnorm(5 * n), n))
  names(d) <- paste0("l", 1:5)
  for (k in 1:5) {
    d[[paste0("x", k)]] <- d[[paste0("l", k)]] + stats::rnorm(n, 0, 0.3)
    d[[paste0("s", k)]] <- 0.3
  }
  comp2 <- stats::runif(n) < stats::plogis(0.2 + 0.8 * d$l5)
  d$y <- ifelse(comp2,
                3 + 0.8 * d$l2 + stats::rnorm(n, 0, exp(-0.3 + 0.3 * d$l4)),
                -2 + 0.8 * d$l1 + stats::rnorm(n, 0, exp(-0.5 + 0.3 * d$l3)))
  fit <- frm(bf(y ~ 1, mu1 ~ me(x1, s1), mu2 ~ me(x2, s2),
                sigma1 ~ me(x3, s3), sigma2 ~ me(x4, s4),
                theta2 ~ me(x5, s5)),
             family = mixture(gaussian(), gaussian()), data = d)
  nd <- data.frame(x1 = 0, x2 = 0, x3 = 0, x4 = 0, x5 = 0.5,
                   s1 = 1, s2 = 1, s3 = 1, s4 = 1, s5 = 2)
  f <- fitted(fit, newdata = nd)
  e0 <- function(p, v = NULL) {
    q <- nd
    if (!is.null(v)) q[[v]] <- q[[v]] + 1
    dm_lp(fit, q, p)
  }
  b1 <- e0("mu1", "x1") - e0("mu1")
  b2 <- e0("mu2", "x2") - e0("mu2")
  b5 <- e0("theta2", "x5") - e0("theta2")
  set.seed(52)
  N <- 2e6
  m1 <- e0("mu1") + b1 * stats::rnorm(N)
  m2 <- e0("mu2") + b2 * stats::rnorm(N)
  t2 <- stats::plogis(e0("theta2") + b5 * 2 * stats::rnorm(N))
  mc <- (1 - t2) * m1 + t2 * m2
  expect_lt(abs(f[, "Estimate"] - mean(mc)) / (stats::sd(mc) / sqrt(N)), 4)
  expect_lt(abs(f[, "Est.Error"] / stats::sd(mc) - 1), 0.02)
})

test_that("a category band under me() is the draw's quantiles", {
  set.seed(7971)
  n <- 400
  d <- data.frame(lat = stats::rnorm(n), sx = 0.4)
  d$xo <- d$lat + stats::rnorm(n, 0, 0.4)
  eta <- 0.4 + 0.7 * d$lat
  d$ord <- as.integer(cut(eta + stats::rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                          labels = FALSE))
  fit <- frm(bf(ord ~ me(xo, sx)), family = cumulative(), data = d)
  nd <- data.frame(xo = 1, sx = 1.2)
  fv <- fitted(fit, newdata = nd)
  ce <- conditional_effects(fit, effects = "xo", categorical = TRUE,
                            int_conditions = list(xo = 1),
                            conditions = data.frame(sx = 1.2))[[1L]]
  # 0.79.1 read category 3 at 0.196 to 0.305, where the draw's quantiles at
  # the estimates are 0.098 to 0.247
  expect_equal(ce$lower__, unname(fv[1, "Q2.5", ]), tolerance = 1e-12)
  expect_equal(ce$upper__, unname(fv[1, "Q97.5", ]), tolerance = 1e-12)
  set.seed(11)
  P <- as.matrix(frm_linpred(fit, type = "response", newdata = data.frame(
    xo = 1 + 1.2 * stats::rnorm(20000), sx = 1.2)))
  q <- apply(P, 2, stats::quantile, c(0.025, 0.975))
  # the coefficients widen the draw's own quantiles a little, no more
  expect_true(all(ce$lower__ <= q[1, ] * 1.02 & ce$upper__ >= q[2, ] * 0.98))
  expect_lt(max(abs(ce$upper__ - ce$lower__) / (q[2, ] - q[1, ])), 1.2)
  # a row with no draw keeps the Wald band on the logit
  d$x <- d$lat
  fit0 <- frm(bf(ord ~ x), family = cumulative(), data = d)
  ce0 <- conditional_effects(fit0, effects = "x", categorical = TRUE,
                             int_conditions = list(x = 1))[[1L]]
  pk <- ce0$estimate__
  sl <- ce0$se__ / (pk * (1 - pk))
  expect_equal(ce0$lower__, stats::plogis(stats::qlogis(pk) -
                                            stats::qnorm(0.975) * sl),
               tolerance = 1e-12)
})
