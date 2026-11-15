# se() and mi() beside a residual correlation term of the ar() family,
# and the cluster sandwich over its series. brms 2.23.0 fits the first
# two: a known se adds se^2 to the diagonal of each group's covariance
# (normal_time_hom_se_lpdf()), and a missing response is a coordinate of
# the group's joint draw (Yl). The identity with brms's own program is in
# test-brms-likelihood.R; this file checks the densities against dense
# references built from the natural parameters alone.

acse_data <- function(seed = 3, G = 10, K = 5, ar = 0.5, tdf = NULL) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  d$s <- stats::runif(n, 0.2, 0.6)
  e <- as.vector(vapply(seq_len(G), function(i) {
    as.vector(stats::arima.sim(list(ar = ar), K, n.start = 100, sd = 0.7))
  }, numeric(K)))
  if (!is.null(tdf)) e <- e * rep(sqrt(tdf / stats::rchisq(G, tdf)), each = K)
  d$y <- 1 + 0.5 * d$x + e + stats::rnorm(n, 0, d$s)
  d
}

# the AR(1) covariance of one group at the natural parameters, plus the
# known se^2: sigma is the innovation sd, as in brms
acse_cov <- function(tt, sigma, ar, se) {
  sigma^2 * ar^abs(outer(tt, tt, "-")) / (1 - ar^2) + diag(se^2, length(tt))
}

acse_nat <- function(fit) {
  ac <- fit$frame$autocor[[1L]]
  th <- fit$estimates$thetaac[ac$theta_idx]
  list(beta = unname(fixef_by_dpar(fit)$mu),
       sigma = unname(exp(fixef_by_dpar(fit)$sigma[["(Intercept)"]])),
       ar = unname(autocor_natural(th, ac)[["ar[1]"]]))
}

test_that("se() beside ar(cov = TRUE) adds se^2 to each group's covariance", {
  skip_if_not_installed("mvtnorm")
  d <- acse_data()
  fit <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)),
             family = gaussian(), data = d)
  expect_identical(fit$opt$convergence, 0L)
  p <- acse_nat(fit)
  mu <- p$beta[1] + p$beta[2] * d$x
  ref <- sum(vapply(split(seq_len(nrow(d)), d$g), function(i) {
    mvtnorm::dmvnorm(d$y[i], mu[i], acse_cov(d$t[i], p$sigma, p$ar, d$s[i]),
                     log = TRUE)
  }, 0))
  ll <- as.numeric(logLik(fit))
  expect_lt(abs(ll - ref) / abs(ref), 1e-10)
})

test_that("se() beside a residual correlation term needs sigma = TRUE", {
  d <- acse_data()
  for (tm in c("ar(t, g, cov = TRUE)", "cosy(t, g)", "unstr(t, g)")) {
    expect_error(
      frm(bf(stats::as.formula(paste("y | se(s) ~ x +", tm))),
          family = gaussian(), data = d),
      "needs se(..., sigma = TRUE)", fixed = TRUE)
  }
  # still refused beside cov = FALSE and fcor(), as brms refuses both
  expect_error(frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g)),
                   family = gaussian(), data = d),
               "se() cannot be combined with cov = FALSE", fixed = TRUE)
  M <- diag(nrow(d))
  expect_error(frm(bf(y | se(s, sigma = TRUE) ~ x + fcor(M)),
                   family = gaussian(), data = d, data2 = list(M = M)),
               "se() cannot be combined", fixed = TRUE)
})

test_that("se() beside cosy() on student() is the multivariate t of D R D + se^2", {
  skip_if_not_installed("mvtnorm")
  d <- acse_data(seed = 5, tdf = 4)
  fit <- frm(bf(y | se(s, sigma = TRUE) ~ x + cosy(t, g)),
             family = student(), data = d)
  ac <- fit$frame$autocor[[1L]]
  rho <- unname(autocor_natural(fit$estimates$thetaac[ac$theta_idx], ac))
  # the dpars at the estimates, each through its own link
  dp <- eval_dpars(fit)[[1L]]
  sigma <- dp$sigma[1]
  nu <- dp$nu[1]
  mu <- rep_len(dp$mu, nrow(d))
  ref <- sum(vapply(split(seq_len(nrow(d)), d$g), function(i) {
    k <- length(i)
    S <- sigma^2 * (diag(1 - rho, k) + rho) + diag(d$s[i]^2, k)
    mvtnorm::dmvt(d$y[i], mu[i], S, df = nu, log = TRUE)
  }, 0))
  ll <- as.numeric(logLik(fit))
  expect_lt(abs(ll - ref) / abs(ref), 1e-9)
})

test_that("simulate() under se() draws the summed covariance", {
  d <- acse_data(seed = 4)
  fit <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)),
             family = gaussian(), data = d)
  p <- acse_nat(fit)
  ns <- 4000
  sims <- as.matrix(simulate(fit, nsim = ns, seed = 11))
  r <- sims - fitted(fit)[, "Estimate"]
  # per-row variance against sigma^2 / (1 - ar^2) + se^2, in units of its
  # own Monte Carlo error, sqrt(2 / (ns - 1)) of a normal variance
  v_ref <- p$sigma^2 / (1 - p$ar^2) + d$s^2
  z <- (apply(r, 1, stats::var) / v_ref - 1) / sqrt(2 / (ns - 1))
  expect_lt(max(abs(z)), stats::qnorm(1 - 0.001 / (2 * nrow(d))))
  # and the lag-one covariance within a group, which se() does not reach
  i1 <- which(d$t < max(d$t))
  c1 <- vapply(i1, function(i) stats::cov(r[i, ], r[i + 1L, ]), 0)
  c_ref <- p$ar * p$sigma^2 / (1 - p$ar^2)
  expect_lt(abs(mean(c1) / c_ref - 1), 0.1)
})

test_that("mi() beside ar(cov = TRUE) integrates a missing row out exactly", {
  skip_if_not_installed("mvtnorm")
  for (fam in list(gaussian(), student())) {
    d <- acse_data(seed = 6, tdf = if (fam$family == "student") 5)
    miss <- c(3, 12, 26, 49)
    d$y[miss] <- NA
    fit <- frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), family = fam,
               data = d)
    # the rows no density reads leave the objective
    expect_setequal(fit$frame$mi_map[[1L]]$unread, miss)
    # the margin of the joint draw over the rows a group keeps is the same
    # block over its remaining time points, so the fit is the one with the
    # rows dropped, whose lags are read off the global time levels; the
    # dropped rows leave interior gaps, which frm() says brms reads as no
    # gap (under mi() the rows stay, and it says nothing)
    expect_message(
      f0 <- frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = fam,
                data = d[-miss, ]),
      class = "frmtmb_autocor_gap")
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(f0))) /
                abs(as.numeric(logLik(f0))), 1e-10)
    if (fam$family != "gaussian") next
    # and the dense margin itself, at the fit's natural parameters
    p <- acse_nat(fit)
    mu <- p$beta[1] + p$beta[2] * d$x
    ok <- !is.na(d$y)
    ref <- sum(vapply(split(which(ok), d$g[ok]), function(i) {
      mvtnorm::dmvnorm(d$y[i], mu[i], acse_cov(d$t[i], p$sigma, p$ar,
                                               numeric(length(i))),
                       log = TRUE)
    }, 0))
    expect_lt(abs(as.numeric(logLik(fit)) - ref) / abs(ref), 1e-10)
  }
})

test_that("mi() beside fcor() integrates the missing values out", {
  d <- acse_data(seed = 7)
  n <- nrow(d)
  M <- 0.4^abs(outer(seq_len(n), seq_len(n), "-"))
  miss <- c(4, 17, 33)
  d$y[miss] <- NA
  fit <- frm(bf(y | mi() ~ x + fcor(M)), family = gaussian(), data = d,
             data2 = list(M = M))
  # no other density reads the values, so they leave the fit, which takes
  # the margin of the observed rows: the fit over them with M cut to them,
  # which is what brms's posterior is
  expect_setequal(fit$frame$mi_map$y$unread, miss)
  expect_length(fit$estimates$miss, 0L)
  f0 <- frm(bf(y ~ x + fcor(M0)), family = gaussian(), data = d[-miss, ],
            data2 = list(M0 = M[-miss, -miss]))
  l1 <- as.numeric(logLik(fit))
  l0 <- as.numeric(logLik(f0))
  expect_lt(abs(l1 - l0) / abs(l0), 1e-7)
  expect_lt(max(abs(fixef(fit)[, "Estimate"] - fixef(f0)[, "Estimate"]) /
                  fixef(f0)[, "Est.Error"]), 1e-4)
})

# a fit with nu held at `nu` and every other parameter at its maximum: one
# t vector does not pin nu, which runs to about 1e10 under ML. nu's
# internal parameter is the second betad, after sigma's, on these models
acse_nu_held <- function(nu, ...) {
  held <- function(par, fn, gr, lower, upper, control) {
    inu <- which(names(par) == "betad")[2L]
    fx <- function(q) {
      p <- numeric(length(par))
      p[-inu] <- q
      p[inu] <- log(nu - 1)
      p
    }
    r <- stats::nlminb(par[-inu], function(q) fn(fx(q)),
                       function(q) gr(fx(q))[-inu])
    p <- par
    p[] <- fx(r$par)
    list(par = p, objective = fn(p), convergence = 0L, message = "nu held")
  }
  # nu held off its maximum leaves its gradient, which frm() reports
  allow_warnings(
    frm(..., control = frmtmb_control(optimizer = held, restarts = 0,
                                      check_se = "ignore")),
    "Large maximum absolute gradient")
}

test_that("mi() beside fcor() on student() integrates out exactly", {
  skip_if_not_installed("mvtnorm")
  # the margin of a multivariate t over the observed rows is a multivariate
  # t; 0.81.0 held the values latent, and the Laplace value of their
  # integral was low by 1.98 at m = 20, nu = 3
  n <- 60
  M <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
  set.seed(3)
  d0 <- data.frame(x = stats::rnorm(n))
  d0$y <- 1 + 0.5 * d0$x +
    as.vector(t(chol(M)) %*% stats::rnorm(n)) * sqrt(4 / stats::rchisq(1, 4))
  for (m in c(3, 10, 20)) {
    d <- d0
    d$y[sort(sample(n, m))] <- NA
    ok <- !is.na(d$y)
    for (nu in c(3, 10)) {
      f <- acse_nu_held(nu, bf(y | mi() ~ x + fcor(M)), family = student(),
                        data = d, data2 = list(M = M))
      dp <- eval_dpars(f)[[1L]]
      expect_lt(abs(dp$nu[1] - nu) / nu, 1e-12)
      ex <- mvtnorm::dmvt(d$y[ok], rep_len(dp$mu, n)[ok],
                          dp$sigma[1]^2 * M[ok, ok], df = nu, log = TRUE)
      expect_lt(abs(as.numeric(logLik(f)) - ex) / abs(ex), 1e-9)
    }
  }
  # brms's default prior pins nu at a moderate value (28.11 here)
  set.seed(4)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x +
    as.vector(t(chol(M)) %*% stats::rnorm(n)) * sqrt(4 / stats::rchisq(1, 4))
  d$y[sort(sample(n, 20))] <- NA
  ok <- !is.na(d$y)
  f <- frm(bf(y | mi() ~ x + fcor(M)), family = student(), data = d,
           data2 = list(M = M), prior = set_prior("gamma(2, 0.1)",
                                                   class = "nu"))
  dp <- eval_dpars(f)[[1L]]
  ex <- mvtnorm::dmvt(d$y[ok], rep_len(dp$mu, n)[ok],
                      dp$sigma[1]^2 * M[ok, ok], df = dp$nu[1], log = TRUE)
  expect_lt(abs(as.numeric(logLik(f)) - ex) / abs(ex), 1e-9)
  # nothing is latent, so the Laplace check has nothing to measure
  expect_setequal(f$frame$mi_map$y$unread, which(!ok))
  expect_true(mi_laplace_exact(f$frame))
  expect_null(f$mi_laplace)
})

test_that("mi() beside fcor() on student() with more integrated is exact", {
  # Held latent, the joint t's Laplace error over the missing values moved
  # with the observed rows' quadratic form, so with a random effect, or
  # with beta under REML, the objective was 0.3 to 3.1 nats low, and
  # silent (mi_laplace_note() read tk = -0.93 on a prior that held nu at
  # 28). The model is the one fit on the observed rows with M cut to them,
  # so its objective is that fit's at every point
  n <- 60
  M <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
  mk <- function(re_sd) {
    set.seed(4)
    d <- data.frame(x = stats::rnorm(n), h = factor(rep(1:6, each = 10)))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, re_sd)[d$h] +
      as.vector(t(chol(M)) %*% stats::rnorm(n)) *
      sqrt(4 / stats::rchisq(1, 4))
    d$y[sort(sample(n, 20))] <- NA
    d
  }
  cases <- list(
    list(f = y ~ x + (1 | h) + fcor(M), f0 = y ~ x + (1 | h) + fcor(M0),
         sd = 0, reml = FALSE, profile = FALSE),
    list(f = y ~ x + (1 | h) + fcor(M), f0 = y ~ x + (1 | h) + fcor(M0),
         sd = 1.5, reml = FALSE, profile = TRUE),
    list(f = y ~ x + fcor(M), f0 = y ~ x + fcor(M0), sd = 0, reml = TRUE,
         profile = FALSE))
  for (cs in cases) {
    d <- mk(cs$sd)
    ok <- !is.na(d$y)
    ctl <- frmtmb_control(profile = cs$profile)
    f <- frm(bf(update(cs$f, y | mi() ~ .)), family = student(), data = d,
             data2 = list(M = M), REML = cs$reml, control = ctl)
    f0 <- frm(bf(cs$f0), family = student(), data = d[ok, ],
              data2 = list(M0 = M[ok, ok]), REML = cs$reml, control = ctl)
    expect_setequal(f$frame$mi_map$y$unread, which(!ok))
    expect_true(mi_laplace_exact(f$frame))
    inu <- grep("^nu", outer_par_labels(f))
    # nu at 3, 6 and 21, where the latent values' Laplace error was large
    for (k in 0:2) {
      p <- f$opt$par + 0.05 * k * seq_along(f$opt$par) / length(f$opt$par)
      p[inu] <- log(c(2, 5, 20)[k + 1L])
      a <- f$obj$fn(p)
      b <- f0$obj$fn(p)
      expect_lt(abs(a - b) / abs(b), 1e-10)
    }
  }
})

test_that("mi() beside fcor() with a random effect in sigma is exact", {
  n <- 60
  M <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
  # gaussian(): held latent, the values' Hessian moved with a random
  # effect in sigma, and the joint Laplace was 1.1 to 2.0 nats low (data
  # seed 9 below). With the random effect in mu only it was exact, and the
  # margin is too
  set.seed(9)
  d <- data.frame(x = stats::rnorm(n), h = factor(rep(1:6, each = 10)))
  d$y <- 1 + 0.5 * d$x + exp(0.6 * stats::rnorm(6)[d$h]) *
    as.vector(t(chol(M)) %*% stats::rnorm(n))
  d$y[sort(sample(n, 15))] <- NA
  ok <- !is.na(d$y)
  for (fam in list(gaussian(), student())) {
    f <- frm(bf(y | mi() ~ x + fcor(M), sigma ~ 1 + (1 | h)), family = fam,
             data = d, data2 = list(M = M))
    f0 <- frm(bf(y ~ x + fcor(M0), sigma ~ 1 + (1 | h)), family = fam,
              data = d[ok, ], data2 = list(M0 = M[ok, ok]))
    expect_setequal(f$frame$mi_map$y$unread, which(!ok))
    isd <- which(names(f$opt$par) == "theta")
    for (s in c(0.3, 0.6, 1.2)) {
      p <- f$opt$par
      p[isd] <- log(s)
      expect_lt(abs(f$obj$fn(p) - f0$obj$fn(p)) / abs(f0$obj$fn(p)), 1e-10)
    }
  }
  fm <- frm(bf(y | mi() ~ x + (1 | h) + fcor(M)), family = gaussian(),
            data = d, data2 = list(M = M))
  fm0 <- frm(bf(y ~ x + (1 | h) + fcor(M0)), family = gaussian(),
             data = d[ok, ], data2 = list(M0 = M[ok, ok]))
  for (k in 0:2) {
    p <- fm$opt$par + 0.05 * k * seq_along(fm$opt$par) / length(fm$opt$par)
    expect_lt(abs(fm$obj$fn(p) - fm0$obj$fn(p)) / abs(fm0$obj$fn(p)), 1e-10)
  }
})

test_that("a registered term's missing values stay latent where read", {
  # the margin is for values only the term's density reads: one that
  # another density reads stays latent, and the sampler holds every value
  n <- 40
  M <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
  set.seed(5)
  d <- data.frame(x = stats::rnorm(n), h = factor(rep(1:4, each = 10)))
  d$y <- 1 + 0.5 * d$x + as.vector(t(chol(M)) %*% stats::rnorm(n))
  miss <- c(3, 11, 25)
  d$y[miss] <- NA
  d$z <- 0.3 + 0.7 * ifelse(is.na(d$y), 1, d$y) + stats::rnorm(n, 0, 0.5)
  fr <- frm(bf(y | mi() ~ x + (1 | h) + fcor(M), family = student()) +
              bf(z ~ mi(y), family = gaussian()) + set_rescor(FALSE),
            data = d, data2 = list(M = M),
            prior = set_prior("gamma(2, 0.1)", class = "nu", resp = "y"),
            control = frmtmb_control(check_mi = "ignore"))
  expect_length(fr$frame$mi_map$y$unread, 0L)
  expect_setequal(fr$frame$mi_map$y$rows, miss)
  fs <- frm(bf(y | mi() ~ x + (1 | h) + fcor(M)), family = student(),
            data = d, data2 = list(M = M),
            prior = set_prior("gamma(2, 0.1)", class = "nu"))
  expect_setequal(fs$frame$mi_map$y$unread, miss)
  fsm <- mi_sampling(fs)
  expect_setequal(fsm$frame$mi_map$y$rows, miss)
  expect_length(fsm$frame$mi_map$y$unread, 0L)
  expect_length(fsm$estimates$miss, length(miss))
})

test_that("the cluster sandwich works over the series of ar(cov = TRUE)", {
  d <- acse_data(seed = 8, G = 12)
  d$h <- factor((as.integer(d$g) + 1L) %/% 2L)
  fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = gaussian(),
             data = d)
  S <- cluster_scores(fit, ~ g)
  # each cluster's score is the gradient of its group's own density, here
  # by central differences of a dense AR(1) normal at the internal
  # parameters (beta, log sigma, the pacf parameter of ar[1])
  X <- stats::model.matrix(~ x, d)
  par <- fit$opt$par
  expect_equal(names(par), c("beta", "beta", "betad", "thetaac"))
  # the internal beta is the natural one when the design is not centered;
  # the comparison is meaningful only then
  expect_lt(max(abs(unname(par[1:2]) - unname(fixef_by_dpar(fit)$mu))), 1e-8)
  ll_g <- function(p, i) {
    sigma <- exp(p[3])
    ar <- p[4] / sqrt(1 + p[4]^2)
    mvtnorm::dmvnorm(d$y[i], as.vector(X[i, ] %*% p[1:2]),
                     acse_cov(d$t[i], sigma, ar, numeric(length(i))),
                     log = TRUE)
  }
  h <- 1e-5
  ref <- t(vapply(split(seq_len(nrow(d)), d$g), function(i) {
    vapply(seq_along(par), function(j) {
      e <- replace(numeric(length(par)), j, h)
      (ll_g(par + e, i) - ll_g(par - e, i)) / (2 * h)
    }, 0)
  }, numeric(length(par))))
  expect_lt(max(abs(unname(S) - unname(ref))) / max(abs(ref)), 1e-6)
  V <- vcov_cluster(fit, ~ h, type = "CR1")
  expect_true(all(is.finite(diag(V))))
  # a clustering that splits a series is refused
  d$k <- factor(d$t %% 2)
  expect_error(vcov_cluster(fit, d$k), "crosses `cluster`", fixed = TRUE)
  # cov = FALSE and fcor(): the first works over the series, the second
  # is one density over all the rows and stays refused
  f1 <- frm(bf(y ~ x + ar(t, g)), family = gaussian(), data = d)
  S1 <- cluster_scores(f1, ~ g)
  expect_lt(max(abs(colSums(S1) + f1$obj$gr(f1$opt$par))) /
              max(abs(S1)), 1e-6)
  M <- diag(nrow(d))
  f2 <- frm(bf(y ~ x + fcor(M)), family = gaussian(), data = d,
            data2 = list(M = M))
  expect_error(vcov_cluster(f2, ~ g), "one joint density over all the rows",
               fixed = TRUE)
})

test_that("marginaleffects reads a residual correlation term's series", {
  skip_if_not_installed("marginaleffects")
  set.seed(2)
  dd <- data.frame(t = rep(1:20, 3), h = factor(rep(1:3, each = 20)),
                   x = stats::rnorm(60))
  dd$y <- dd$x + stats::rnorm(60)
  # marginaleffects pads newdata with a row per level of h, copied from
  # its first row, so the padded rows repeat a time of the rows asked
  # for in their group; under cov = FALSE they are predicted apart
  f0 <- frm(bf(y ~ x + ar(time = t, gr = h)), family = gaussian(),
            data = dd)
  nd <- dd[c(1, 2, 21), ]
  p <- marginaleffects::predictions(f0, newdata = nd)
  expect_identical(nrow(p), 3L)
  ref <- fitted(f0, newdata = nd)[, "Estimate"]
  expect_lt(max(abs(p$estimate - ref)) / max(abs(ref)), 1e-10)
  # a cluster formula over the term's group is the sandwich
  f1 <- frm(bf(y ~ x + ar(time = t, gr = h, cov = TRUE)), family = gaussian(),
            data = dd)
  a <- marginaleffects::avg_predictions(f1, vcov = ~ h)
  expect_true(is.finite(a$std.error) && a$std.error > 0)
  expect_error(marginaleffects::avg_predictions(f1, vcov = "HAC"),
               "orders the rows within each group only")
})

test_that("a group that skips a time level says that brms counts rows", {
  d <- acse_data(seed = 9)
  dg <- d[-c(3, 14), ]
  expect_message(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = gaussian(),
                     data = dg),
                 class = "frmtmb_autocor_gap")
  # the guard's absent case: no gap, a gap only cosy() would see, and a
  # gap that every group has (which warns instead)
  expect_no_message(frm(bf(y ~ x + ar(t, g, cov = TRUE)),
                        family = gaussian(), data = d),
                    class = "frmtmb_autocor_gap")
  expect_no_message(frm(bf(y ~ x + cosy(t, g)), family = gaussian(),
                        data = dg),
                    class = "frmtmb_autocor_gap")
  expect_no_message(frm(bf(y ~ x + ar(t, g)), family = gaussian(),
                        data = dg),
                    class = "frmtmb_autocor_gap")
  d3 <- d[d$t != 3, ]
  allow_warnings(
    expect_no_message(frm(bf(y ~ x + ar(t, g, cov = TRUE)),
                          family = gaussian(), data = d3),
                      class = "frmtmb_autocor_gap"),
    "not consecutive")
})
