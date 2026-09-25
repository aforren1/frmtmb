# brms's cov = FALSE ARMA: ar(), ma() and arma() without cov = TRUE.
#
# The reference is brms 2.23.0's own model block, transliterated below
# from brms::make_stancode() and run over the data in brms's order with
# brms's J_lag, built the way brms's data_ac() builds it. Nothing on the
# reference side calls frmtmb's recursion or its row bookkeeping. The
# classical reference is stats::arima(method = "CSS"), the conditional
# sum of squares, which is the same likelihood for one series.
# dev/arcov-validate.R runs the wider grid and records the numbers.

# for (n in 1:N) {
#   mu[n] += Err[n, 1:Kma] * ma;
#   err[n] = Y[n] - mu[n];
#   for (i in 1:J_lag[n]) Err[n + 1, i] = err[n + 1 - i];
#   mu[n] += Err[n, 1:Kar] * ar;
# }
cond_brms_mu <- function(mu, Y, J_lag, ar, ma) {
  N <- length(Y)
  max_lag <- max(length(ar), length(ma))
  Err <- matrix(0, N + 1, max_lag)
  err <- numeric(N)
  for (n in seq_len(N)) {
    if (length(ma)) mu[n] <- mu[n] + sum(Err[n, seq_along(ma)] * ma)
    err[n] <- Y[n] - mu[n]
    for (i in seq_len(J_lag[n])) Err[n + 1, i] <- err[n + 1 - i]
    if (length(ar)) mu[n] <- mu[n] + sum(Err[n, seq_along(ar)] * ar)
  }
  mu
}

# brms's data_ac(): after sorting by (gr, time), row n counts the rows
# n, n - 1, ..., n + 1 - max_lag that share the group of row n + 1
cond_brms_jlag <- function(g, max_lag) {
  N <- length(g)
  J <- integer(N)
  for (n in seq_len(N - 1L)) {
    ind <- n:max(1L, n + 1L - max_lag)
    J[n] <- sum(g[ind] == g[n + 1L])
  }
  J
}

# ragged groups with interior gaps, rows shuffled: brms counts lags in
# rows of the (gr, time) order, so neither the gaps nor the data order
# may matter
cond_data <- function(seed = 3, G = 8, K = 7) {
  set.seed(seed)
  d <- expand.grid(week = seq_len(K), subj = factor(seq_len(G)))
  d <- d[-c(3, 10, 11, 30, 44), ]
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(nrow(d))
  d$w <- stats::runif(nrow(d), 0.5, 2)
  d$cc <- sample(c(0, 0, 0, 1, -1), nrow(d), TRUE)
  d <- d[sample(nrow(d)), ]
  rownames(d) <- NULL
  d
}

# the fit's natural parameters at a full parameter vector, plus the
# brms-order view of the data
cond_ref_parts <- function(f, d, full) {
  pl <- lapply(split(unname(full), factor(names(full),
                                          unique(names(full)))), identity)
  ac <- f$frame$autocor[[1L]]
  th <- pl$thetaac[ac$theta_idx]
  ord <- order(d$subj, d$week)
  X <- stats::model.matrix(~ x, d)
  mu0 <- drop(X %*% pl$beta)
  if (!is.null(pl$b)) mu0 <- mu0 + pl$b[as.integer(d$subj)]
  list(pl = pl, ord = ord, mu0 = mu0[ord], y = d$y[ord],
       J = cond_brms_jlag(as.integer(d$subj)[ord], max(ac$p, ac$q)),
       ar = th[seq_len(ac$p)], ma = th[ac$p + seq_len(ac$q)])
}

cond_points <- function(f) {
  full0 <- f$obj$env$last.par.best
  set.seed(11)
  list(full0, full0 + stats::rnorm(length(full0), 0, 0.05))
}

cond_tol <- 1e3 * .Machine$double.eps

test_that("the likelihood is brms's recursion (gaussian, ragged, shuffled)", {
  d <- cond_data()
  for (tm in c("ar(week, subj)", "ma(week, subj)",
               "arma(week, subj, p = 2, q = 1)",
               "ar(week, subj, p = 3)")) {
    f <- frm(stats::as.formula(paste("y ~ x +", tm)), data = d,
             family = gaussian())
    for (full in cond_points(f)) {
      r <- cond_ref_parts(f, d, full)
      mu <- cond_brms_mu(r$mu0, r$y, r$J, r$ar, r$ma)
      ref <- sum(stats::dnorm(r$y, mu, exp(r$pl$betad), log = TRUE))
      expect_equal(-f$obj$env$f(full), ref, tolerance = cond_tol,
                   info = tm)
    }
  }
})

test_that("J_lag is the one brms builds", {
  skip_if_not_installed("brms")
  d <- cond_data()
  sd <- brms::make_standata(y ~ x + arma(week, subj, p = 2, q = 1),
                            data = d)
  ord <- order(d$subj, d$week)
  expect_equal(as.integer(sd$J_lag),
               cond_brms_jlag(as.integer(d$subj)[ord], 2L))
  expect_equal(as.numeric(sd$Y), d$y[ord])
})

test_that("weights(), cens() and a random intercept keep brms's density", {
  d <- cond_data()
  fw <- frm(y | weights(w) ~ x + ma(week, subj), data = d,
            family = gaussian())
  fc <- frm(y | cens(cc) ~ x + arma(week, subj), data = d,
            family = gaussian())
  fr <- frm(y ~ x + (1 | subj) + arma(week, subj), data = d,
            family = gaussian())
  for (full in cond_points(fw)) {
    r <- cond_ref_parts(fw, d, full)
    mu <- cond_brms_mu(r$mu0, r$y, r$J, r$ar, r$ma)
    ref <- sum(d$w[r$ord] *
                 stats::dnorm(r$y, mu, exp(r$pl$betad), log = TRUE))
    expect_equal(-fw$obj$env$f(full), ref, tolerance = cond_tol)
  }
  for (full in cond_points(fc)) {
    r <- cond_ref_parts(fc, d, full)
    mu <- cond_brms_mu(r$mu0, r$y, r$J, r$ar, r$ma)
    s <- exp(r$pl$betad)
    cc <- d$cc[r$ord]
    ref <- sum(ifelse(
      cc == 0, stats::dnorm(r$y, mu, s, log = TRUE),
      ifelse(cc == 1,
             stats::pnorm(r$y, mu, s, lower.tail = FALSE, log.p = TRUE),
             stats::pnorm(r$y, mu, s, log.p = TRUE))))
    expect_equal(-fc$obj$env$f(full), ref, tolerance = cond_tol)
  }
  # the joint density at the modes: rows at the conditional mean plus
  # the 1-d block's dnorm(b, 0, exp(theta))
  for (full in cond_points(fr)) {
    r <- cond_ref_parts(fr, d, full)
    mu <- cond_brms_mu(r$mu0, r$y, r$J, r$ar, r$ma)
    ref <- sum(stats::dnorm(r$y, mu, exp(r$pl$betad), log = TRUE)) +
      sum(stats::dnorm(r$pl$b, 0, exp(r$pl$theta), log = TRUE))
    expect_equal(-fr$obj$env$f(full), ref, tolerance = cond_tol)
  }
})

test_that("one long series matches the conditional sum of squares", {
  set.seed(5)
  n <- 400
  s <- data.frame(t = seq_len(n),
                  y = 1 + as.numeric(stats::arima.sim(list(ma = 0.5), n)))
  f <- frm(y ~ 1 + ma(t), data = s, family = gaussian())
  a <- stats::arima(s$y, order = c(0, 0, 1), method = "CSS",
                    optim.control = list(reltol = 1e-15, maxit = 5000))
  # the objectives are one function: frmtmb's log-likelihood at arima's
  # estimates is the gaussian density of arima's own CSS residuals
  pr <- f$obj$par
  pr[names(pr) == "beta"] <- a$coef[["intercept"]]
  pr[names(pr) == "betad"] <- 0.5 * log(a$sigma2)
  pr[names(pr) == "thetaac"] <- a$coef[["ma1"]]
  expect_equal(-f$obj$fn(pr),
               sum(stats::dnorm(stats::residuals(a), 0, sqrt(a$sigma2),
                                log = TRUE)),
               tolerance = cond_tol)
  # and the two optimizers land on one optimum, to a small fraction of
  # the standard error the fit itself reports
  se <- sqrt(diag(vcov(f, full = TRUE)))
  est <- c(f$estimates$beta, f$estimates$thetaac)
  css <- c(a$coef[["intercept"]], a$coef[["ma1"]])
  expect_lt(max(abs(est - css) / se[c(1L, 3L)]), 1e-2)
})

test_that("families without a residual are refused in brms's words", {
  d <- cond_data()
  d$cnt <- stats::rpois(nrow(d), 3)
  expect_error(frm(cnt ~ x + ma(week, subj), data = d, family = poisson()),
               "Please set cov = TRUE when modeling MA structures")
  expect_error(frm(cnt ~ x + arma(week, subj), data = d, family = poisson()),
               "Please set cov = TRUE when modeling MA structures")
  # brms fits latent residuals for an AR term here; that is refused by
  # name, with the random-effect spelling of a latent AR process
  expect_error(frm(cnt ~ x + ar(week, subj), data = d, family = poisson()),
               "LATENT residuals.*ar1[(]factor[(]week[)]")
  d$s <- 0.3
  expect_error(frm(y | se(s) ~ x + ar(week, subj), data = d,
                   family = gaussian()),
               "Please set cov = TRUE in ARMA structures")
  # a lag longer than every group names the order, rather than fitting a
  # coefficient that never enters the likelihood
  expect_error(frm(y ~ x + ar(week, subj, p = 7), data = d,
                   family = gaussian()),
               "never reaches an earlier row")
})

test_that("fitted() is brms's one-step mean, residuals follow it", {
  d <- cond_data()
  f <- frm(y ~ x + arma(week, subj), data = d, family = gaussian())
  r <- cond_ref_parts(f, d, f$obj$env$last.par.best)
  mu <- cond_brms_mu(r$mu0, r$y, r$J, r$ar, r$ma)
  fv <- fitted(f)[, "Estimate"]
  expect_equal(unname(fv[r$ord]), unname(mu), tolerance = cond_tol)
  expect_equal(unname(residuals(f)[, "Estimate"]), d$y - fv,
               tolerance = cond_tol)
  expect_equal(as.vector(frm_linpred(f, type = "link")), unname(fv),
               tolerance = cond_tol)
  # newdata that carries the response is the same recursion over
  # newdata's own groups, so the training rows give fitted() back
  fn <- fitted(f, newdata = d)
  expect_equal(fn[, "Estimate"], fv, tolerance = cond_tol)
  expect_true(all(is.finite(fn[, "Est.Error"]) & fn[, "Est.Error"] > 0))
  expect_equal(fn[, "Est.Error"], fitted(f)[, "Est.Error"],
               tolerance = cond_tol)
  expect_error(fitted(f, newdata = transform(d, y = NULL)),
               "needs the observed response 'y'")
})

test_that("predict() draws around the one-step mean", {
  d <- cond_data()
  f <- frm(y ~ x + ar(week, subj), data = d, family = gaussian())
  nd <- 4000L
  set.seed(2)
  p <- predict(f, ndraws = nd, propagate_error = FALSE)
  fv <- fitted(f)[, "Estimate"]
  # a Monte Carlo mean, judged against its own Monte Carlo error
  z <- (p[, "Estimate"] - fv) / (p[, "Est.Error"] / sqrt(nd))
  expect_lt(max(abs(z)), 5)
  # the one-step predictive sd is sigma, not the larger marginal sd; a
  # draw sd has relative Monte Carlo error 1 / sqrt(2 nd), and the mean
  # over the rows divides that by sqrt(rows) once more
  zs <- (mean(p[, "Est.Error"]) - sigma(f)) /
    (sigma(f) / sqrt(2 * nd * nrow(d)))
  expect_lt(abs(zs), 5)
})

test_that("simulate() runs the recursion over its own draws", {
  d <- cond_data()
  f <- frm(y ~ x + arma(week, subj, p = 2, q = 1), data = d,
           family = gaussian())
  y1 <- simulate(f, nsim = 1, seed = 9)[[1L]]
  # the innovations of a replicate, recovered with brms's recursion run
  # over the replicate itself, are exactly the normal draws it was built
  # from: the same stream, in whatever order the positions drew them
  g <- f
  g$frame$y[[1L]] <- y1
  e <- (y1 - fitted(g)[, "Estimate"]) / sigma(f)
  set.seed(9)
  z <- stats::rnorm(nrow(d))
  expect_equal(sort(e), sort(z), tolerance = cond_tol)
  # a new group on newdata starts from an empty past
  nw <- data.frame(week = 1:5, subj = factor(99), x = 0)
  s2 <- simulate(f, nsim = 2, newdata = nw, seed = 1,
                 allow_new_levels = TRUE)
  expect_equal(dim(s2), c(5L, 2L))
  expect_true(all(is.finite(as.matrix(s2))))
})

test_that("the coefficients are unconstrained, named and bounded as brms's", {
  d <- cond_data()
  f <- frm(y ~ x + arma(week, subj, p = 2, q = 1), data = d,
           family = gaussian())
  cv <- confint_varcorr(f)
  expect_equal(cv$term, c("ar[1]", "ar[2]", "ma[1]"))
  expect_equal(cv$type, rep("raw", 3L))
  expect_equal(cv$estimate, as.numeric(f$estimates$thetaac),
               tolerance = cond_tol)
  expect_error(autocor_matrix(f), "defines no correlation matrix")
  expect_error(residuals(f, type = "osa"), "not available")
  # a prior on class "ar" is a density on the coefficient itself
  fp <- frm(y ~ x + ar(week, subj), data = d, family = gaussian(),
            prior = set_prior("normal(0.3, 0.001)", class = "ar"))
  f0 <- frm(y ~ x + ar(week, subj), data = d, family = gaussian())
  expect_lt(abs(fp$estimates$thetaac - 0.3),
            1e-2 * abs(f0$estimates$thetaac - 0.3))
  # and a box on an order-2 AR is a box on both coefficients, which the
  # transformed cov = TRUE parameterization cannot offer
  fb <- frm(y ~ x + ar(week, subj, p = 2), data = d, family = gaussian(),
            prior = set_prior("", class = "ar", lb = -0.05, ub = 0.05))
  expect_true(all(abs(fb$estimates$thetaac) <= 0.05))
})

test_that("conditional_effects() drops the term, as brms does", {
  d <- cond_data()
  f <- frm(y ~ x + ar(week, subj), data = d, family = gaussian())
  ce <- conditional_effects(f, "x")[[1L]]
  b <- f$estimates$beta
  expect_equal(ce$estimate__, unname(b[1] + b[2] * ce$x),
               tolerance = cond_tol)
})

test_that("rescor and student() fit with the term", {
  d <- cond_data()
  d$y2 <- d$y + stats::rnorm(nrow(d))
  fm <- frm(bf(y ~ x + ar(week, subj)) + bf(y2 ~ x + ar(week, subj)) +
              set_rescor(TRUE), data = d, family = gaussian())
  expect_length(fm$frame$autocor, 2L)
  expect_equal(dim(fitted(fm)), c(nrow(d), 4L, 2L))
  fs <- frm(bf(y ~ x + ma(week, subj), nu ~ x), data = d,
            family = student())
  expect_true(is.finite(logLik(fs)))
})
