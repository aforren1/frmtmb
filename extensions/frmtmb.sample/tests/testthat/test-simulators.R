# The simulator contract (frmtmb's R/families.R): one implementation per
# family, reached identically by simulate() on a fit,
# posterior_predict() on a draw, and frm_simulate() de novo.
#
# The file lives HERE, in the sampling package, because one of its three
# entry points does. The claim it makes is a cross-package one now -
# frmtmb's simulator, reached through frmtmb's simulate() and
# frm_simulate() and through this package's posterior_predict(), giving
# the same draw - and this is the only suite that can make it. Splitting
# it so that frmtmb kept a two-way version would have left two copies of
# 350 lines to drift apart, and the two-way claim is the weaker half.
#
# Two things are asserted for every family that gained a simulator here.
# DISTRIBUTIONAL: seeded moments or proportions match what the fitted
# parameters imply, loosely. CROSS-PATH: the three entry points agree.
# They agree exactly, not just in distribution, because the contract
# gives them one implementation and the same RNG call sequence - which
# is a stronger statement than moment agreement and the reason the
# unification is worth having.

# a frmtmb_draws whose every draw is the same parameter vector, so
# posterior_predict() runs its real code path (draws_fit_at, the
# context, the simulator) without a NUTS run in the test suite
ml_draws <- function(fit, R = 20L) {
  tpl <- fit$frame$par_template
  v <- numeric(0)
  for (cp in names(tpl)) {
    e <- fit$estimates[[cp]]
    if (cp == "betad" && length(fx <- fit$frame$betad_fixed_idx)) {
      e <- e[-fx]
    }
    v <- c(v, as.numeric(e))
  }
  structure(list(stanfit = NULL,
                 draws = matrix(v, R, length(v), byrow = TRUE),
                 fit = fit),
            class = "frmtmb_draws")
}

# the internal-spelling newparams that reproduce a fit's parameters:
# every optimized component plus whatever family-level extras the frame
# declares (ordinal thresholds, lca item profiles, mvn covariances)
np_of <- function(fit) {
  comps <- c("beta", "betad", "theta", "thetaac", "thetar",
             fit$frame$extra_names %||% character(0))
  out <- list()
  for (cp in comps) {
    if (!is.null(fit$estimates[[cp]])) {
      out[[cp]] <- as.numeric(fit$estimates[[cp]])
    }
  }
  out
}

## ---- hmm: the chain walk reaches all three entry points --------------

sim_hmm_data <- function(ng, tt, G, mu, sg, seed) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(ng), function(i) {
    s <- integer(tt)
    s[1L] <- 1L
    for (k in 2:tt) s[k] <- sample.int(2L, 1L, prob = G[s[k - 1L], ])
    data.frame(id = i, t = seq_len(tt), y = stats::rnorm(tt, mu[s], sg[s]))
  }))
}

test_that("an hmm() draw walks the chain from every entry point", {
  G2 <- matrix(c(0.9, 0.1, 0.25, 0.75), 2L, 2L, byrow = TRUE)
  tt <- 20L
  dd <- sim_hmm_data(20L, tt, G2, c(0, 4), c(0.5, 0.5), 4001)
  form <- bf(y ~ 1) + hmm(K = 2, gaussian(), time = t, group = id)
  fit <- frm(form, data = dd)

  e <- unlist(fixef(fit))
  mu <- c(e[["mu1.(Intercept)"]], e[["mu2.(Intercept)"]])
  sg <- exp(c(e[["sigma1.(Intercept)"]], e[["sigma2.(Intercept)"]]))
  G <- rbind(c(1, exp(e[["tr12.(Intercept)"]])),
             c(1, exp(e[["tr22.(Intercept)"]])))
  G <- G / rowSums(G)
  pi_ <- as.vector(solve(t(diag(2L) - G + 1), rep(1, 2L)))

  set.seed(7); s_fit <- as.matrix(simulate(fit, nsim = 20))
  set.seed(7); s_pp <- posterior_predict(ml_draws(fit, 20L))
  set.seed(7)
  s_new <- as.matrix(frm_simulate(form, dd, nsim = 20,
                                  newparams = np_of(fit)))

  # distributional: the marginal is the stationary mixture of the two
  # state-dependent normals
  expect_equal(mean(s_fit), sum(pi_ * mu), tolerance = 0.1)
  expect_equal(stats::var(as.vector(s_fit)),
               sum(pi_ * (sg^2 + mu^2)) - sum(pi_ * mu)^2,
               tolerance = 0.1)

  # the chain, not just the marginal: how often consecutive rows of a
  # sequence fall on opposite sides of the two means
  switch_rate <- function(y) {
    m <- matrix(y > mean(mu), nrow = tt)
    mean(m[-1L, ] != m[-tt, ])
  }
  expect_equal(mean(apply(s_fit, 2L, switch_rate)),
               1 - sum(pi_ * diag(G)), tolerance = 0.1)

  # cross-path: one implementation, one RNG sequence
  expect_equal(unname(s_pp), unname(t(s_fit)))
  expect_equal(unname(s_new), unname(s_fit))
})

## ---- group-level mixture ---------------------------------------------

test_that("mixture(groups =) draws one class per group everywhere", {
  set.seed(31)
  ng <- 50L
  m <- 8L
  cls <- stats::rbinom(ng, 1L, 0.4)
  g <- rep(seq_len(ng), each = m)
  dd <- data.frame(y = stats::rnorm(ng * m, c(-2, 2)[cls + 1L][g], 0.5),
                   g = factor(g))
  form <- bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g)
  fit <- frm(form, data = dd)

  set.seed(11); s_fit <- as.matrix(simulate(fit, nsim = 20))
  set.seed(11); s_pp <- posterior_predict(ml_draws(fit, 20L))
  set.seed(11)
  s_new <- as.matrix(frm_simulate(form, dd, nsim = 20,
                                  newparams = np_of(fit)))

  # the structural property the rowwise simulator cannot produce: every
  # row of a group comes from the SAME component
  pure <- function(y) {
    mm <- matrix(y > 0, nrow = m)
    mean(apply(mm, 2L, function(z) all(z) || !any(z)))
  }
  expect_gt(mean(apply(s_fit, 2L, pure)), 0.97)
  expect_gt(mean(apply(s_new, 2L, pure)), 0.97)
  # a rowwise draw of the same mixture would be pure only 2^-8 of the
  # time, so the check has real power
  expect_lt(2^-(m - 1L), 0.02)

  # the mixing weight itself: the share of GROUPS in the high component
  pi_hi <- unname(exp(fit$spec$responses[[1L]]$family$mix$log_pi(
    eval_dpars(fit)[[1L]])[[2L]][1L]))
  grp_hi <- function(y) mean(colMeans(matrix(y, nrow = m)) > 0)
  expect_equal(mean(apply(s_fit, 2L, grp_hi)), pi_hi, tolerance = 0.12)

  expect_equal(unname(s_pp), unname(t(s_fit)))
  expect_equal(unname(s_new), unname(s_fit))
})

## ---- mixture_mvn ------------------------------------------------------

test_that("mixture_mvn simulates its class covariances", {
  set.seed(46)
  n <- 300L
  cl <- stats::rbinom(n, 1L, 0.4)
  L1 <- t(chol(matrix(c(1, 0.5, 0.5, 1), 2L)))
  L2 <- t(chol(matrix(c(0.5, -0.2, -0.2, 0.8), 2L)))
  E <- matrix(stats::rnorm(2L * n), 2L)
  Y <- t(ifelse(matrix(cl == 1L, 2L, n, byrow = TRUE),
                c(0, 0) + L1 %*% E, c(3, 4) + L2 %*% E))
  dd <- data.frame(row = seq_len(n))
  dd$Y <- Y
  form <- bf(Y ~ 1) + mixture_mvn(K = 2, D = 2)
  fit <- frm(form, data = dd)

  set.seed(3); s_fit <- simulate(fit, nsim = 20)
  set.seed(3); s_pp <- posterior_predict(ml_draws(fit, 20L))
  set.seed(3)
  s_new <- frm_simulate(form, dd, nsim = 20, newparams = np_of(fit))

  # a draw is an n x D matrix, not a flattened vector
  expect_true(is.matrix(s_fit[[1L]]))
  expect_equal(dim(s_fit[[1L]]), c(n, 2L))
  expect_equal(dim(s_pp), c(20L, n, 2L))

  # distributional: the fitted two-component mixture reproduces the
  # data's own mean vector and covariance
  M <- do.call(rbind, as.list(s_fit))
  expect_equal(colMeans(M), colMeans(Y), tolerance = 0.15,
               ignore_attr = TRUE)
  expect_equal(as.vector(stats::cov(M)), as.vector(stats::cov(Y)),
               tolerance = 0.15)
  # the class covariances really differ, so the draw is not one normal
  expect_false(isTRUE(all.equal(
    as.matrix(fit$spec$responses[[1L]]$family$mix$sigma(
      fit$estimates[c("sigmaraw1", "sigmaraw2")], 1L)),
    as.matrix(fit$spec$responses[[1L]]$family$mix$sigma(
      fit$estimates[c("sigmaraw1", "sigmaraw2")], 2L)))))

  expect_equal(unname(s_pp[1L, , ]), unname(s_fit[[1L]]))
  expect_equal(unname(s_new[[1L]]), unname(s_fit[[1L]]))
})

## ---- residual correlation (frame-level, not family-level) ------------

test_that("an autocor residual is one group draw on every path", {
  set.seed(77)
  G <- 30L
  K <- 8L
  dd <- expand.grid(week = seq_len(K), subj = factor(seq_len(G)))
  dd <- dd[order(dd$subj, dd$week), ]
  dd$x <- stats::rnorm(nrow(dd))
  e <- as.vector(vapply(seq_len(G), function(i) {
    as.vector(stats::arima.sim(list(ar = 0.6), K, n.start = 200, sd = 1))
  }, numeric(K)))
  dd$y <- 1 + 0.5 * dd$x + e
  form <- bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian()
  fit <- frm(form, data = dd)

  set.seed(5); s_fit <- as.matrix(simulate(fit, nsim = 20))
  set.seed(5); s_pp <- posterior_predict(ml_draws(fit, 20L))
  set.seed(5)
  s_new <- as.matrix(frm_simulate(form, dd, nsim = 20,
                                  newparams = np_of(fit)))

  ac <- fit$frame$autocor[[1L]]
  phi <- unname(frmtmb:::autocor_natural(fit$estimates$thetaac[ac$theta_idx],
                                ac))[1L]
  lag1 <- function(y) {
    r <- matrix(y - as.vector(fitted(fit)), nrow = K)
    mean(vapply(seq_len(K - 1L), function(k) {
      stats::cor(r[k, ], r[k + 1L, ])
    }, numeric(1)))
  }
  # independent rows would give 0 here; the AR(1) gives phi
  expect_equal(mean(apply(s_fit, 2L, lag1)), phi, tolerance = 0.12)
  expect_equal(mean(apply(s_new, 2L, lag1)), phi, tolerance = 0.12)

  expect_equal(unname(s_pp), unname(t(s_fit)))
  expect_equal(unname(s_new), unname(s_fit))
})

## ---- lca: the extras-aware rowwise contract, all three paths ---------

test_that("lca() simulates item codes from every entry point", {
  set.seed(5)
  n <- 300L
  pr <- rbind(c(0.90, 0.85, 0.20, 0.75), c(0.20, 0.15, 0.85, 0.25))
  cl <- stats::rbinom(n, 1L, 0.35) + 1L
  Y <- matrix(0L, n, ncol(pr))
  for (j in seq_len(ncol(pr))) {
    Y[, j] <- 1L + stats::rbinom(n, 1L, pr[cl, j])
  }
  dd <- data.frame(row = seq_len(n))
  dd$Y <- Y
  form <- bf(Y ~ 1) + lca(K = 2)
  fit <- frm(form, data = dd)

  set.seed(13); s_fit <- simulate(fit, nsim = 10)
  set.seed(13); s_pp <- posterior_predict(ml_draws(fit, 10L))
  set.seed(13)
  s_new <- frm_simulate(form, dd, nsim = 10, newparams = np_of(fit))

  expect_equal(dim(s_fit[[1L]]), dim(Y))
  # per-item share of category 2, against the data
  M <- do.call(rbind, as.list(s_fit))
  expect_equal(colMeans(M == 2L), colMeans(Y == 2L), tolerance = 0.06)

  expect_equal(unname(s_pp[1L, , ]), unname(s_fit[[1L]]))
  expect_equal(unname(s_new[[1L]]), unname(s_fit[[1L]]))
})

## ---- the rowwise families still work on all three paths -------------

test_that("ordinal, categorical and multinomial keep all three paths", {
  set.seed(19)
  n <- 200L
  do_ <- data.frame(x = stats::rnorm(n))
  do_$y <- factor(cut(0.8 * do_$x + stats::rnorm(n),
                      c(-Inf, -0.5, 0.5, Inf), labels = c("a", "b", "c")),
                  ordered = TRUE)
  fo <- frm(bf(y ~ x) + cumulative(), data = do_)
  set.seed(2); so <- simulate(fo, nsim = 10)
  set.seed(2); po <- posterior_predict(ml_draws(fo, 10L))
  set.seed(2)
  no <- frm_simulate(bf(y ~ x) + cumulative(), do_, nsim = 10,
                     newparams = np_of(fo))
  # simulate() restores the response's own ordered levels, and so does
  # the de novo path since v0.36
  expect_s3_class(so[[1L]], "ordered")
  expect_s3_class(no[[1L]], "ordered")
  expect_equal(levels(no[[1L]]), levels(do_$y))
  expect_equal(as.integer(so[[1L]]), unname(po[1L, ]))
  expect_equal(prop.table(table(unlist(so))),
               prop.table(table(do_$y)), tolerance = 0.08)

  dc <- data.frame(x = stats::rnorm(n))
  dc$y <- factor(sample(c("a", "b", "c"), n, TRUE, c(0.5, 0.3, 0.2)))
  fc <- frm(bf(y ~ x) + categorical(), data = dc)
  set.seed(4); sc <- simulate(fc, nsim = 10)
  set.seed(4); pc <- posterior_predict(ml_draws(fc, 10L))
  # a nominal response comes back UNORDERED
  expect_s3_class(sc[[1L]], "factor")
  expect_false(is.ordered(sc[[1L]]))
  expect_equal(as.integer(sc[[1L]]), unname(pc[1L, ]))
  expect_equal(prop.table(table(unlist(sc))),
               prop.table(table(dc$y)), tolerance = 0.08)

  dm <- data.frame(x = stats::rnorm(80L), n = 10)
  dm$Y <- t(vapply(seq_len(80L),
                   function(i) stats::rmultinom(1L, 10L, c(0.3, 0.3, 0.4)),
                   numeric(3L)))
  fm <- frm(bf(Y | trials(n) ~ x) + multinomial(K = 3), data = dm)
  set.seed(6); sm <- simulate(fm, nsim = 10)
  set.seed(6); pm <- posterior_predict(ml_draws(fm, 10L))
  expect_equal(dim(pm), c(10L, 80L, 3L))
  expect_equal(unname(pm[1L, , ]), unname(sm[[1L]]))
  expect_equal(colMeans(do.call(rbind, as.list(sm))) / 10,
               colMeans(dm$Y) / 10, tolerance = 0.05,
               ignore_attr = TRUE)
})

## ---- refusals ---------------------------------------------------------

test_that("cox() refuses to simulate the same way at all three doors", {
  set.seed(23)
  dd <- data.frame(x = stats::rnorm(80L))
  dd$y <- stats::rexp(80L, exp(0.3 * dd$x))
  fit <- frm(bf(y ~ x) + cox(), data = dd)

  m <- vapply(list(
    function() simulate(fit),
    function() posterior_predict(ml_draws(fit, 2L)),
    function() frm_simulate(bf(y ~ x) + cox(), dd,
                            newparams = list(beta = c(0, 0.3)))
  ), function(f) tryCatch({ f(); "" }, error = conditionMessage), "")

  # unique per entry point (G5.2a), so a reported message names one line
  expect_length(unique(m), 3L)
  expect_match(m[1L], "^simulate\\(\\)")
  expect_match(m[2L], "^posterior_predict\\(\\)")
  expect_match(m[3L], "^frm_simulate\\(\\)")
  # and consistent: the same family-level reason on each
  expect_true(all(grepl("cumulative baseline hazard", m, fixed = TRUE)))
  expect_true(all(grepl("no simulator yet", m, fixed = TRUE)))
})

test_that("a structured draw refuses trunc() and newdata", {
  set.seed(31)
  ng <- 20L
  m <- 6L
  g <- rep(seq_len(ng), each = m)
  cls <- stats::rbinom(ng, 1L, 0.5)
  dd <- data.frame(y = stats::rnorm(ng * m, c(-2, 2)[cls + 1L][g], 0.5),
                   g = factor(g), lb = -10)
  fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
             data = dd)
  # the group structure indexes the training rows
  expect_error(
    posterior_predict(ml_draws(fit, 2L), newdata = dd[1:4, ]),
    "structured")

  # trunc() rejection resamples single rows, which a whole-group draw
  # cannot supply. Every structured model refuses trunc() at frame
  # assembly already - a mixture density carries no CDF, and an autocor
  # or hmm likelihood does not factorize by row - so the simulator's own
  # guard is a backstop for a future structured family that does carry
  # one. What the user meets is the upstream refusal.
  expect_error(
    frm(bf(y | trunc(lb = lb) ~ 1) +
          mixture(gaussian(), gaussian(), groups = ~g), data = dd),
    "CDF")
  d2 <- expand.grid(week = 1:5, subj = factor(1:8))
  d2 <- d2[order(d2$subj, d2$week), ]
  d2$lb <- -6
  d2$y <- stats::rnorm(nrow(d2))
  expect_error(
    frm(bf(y | trunc(lb = lb) ~ ar(week, subj, cov = TRUE)) + gaussian(),
        data = d2),
    "residual correlation term")
})

test_that("simulate() still refuses re.form and censored on an hmm", {
  G2 <- matrix(c(0.9, 0.1, 0.25, 0.75), 2L, 2L, byrow = TRUE)
  dd <- sim_hmm_data(10L, 10L, G2, c(0, 4), c(0.5, 0.5), 91)
  fit <- frm(bf(y ~ 1) + hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  expect_error(simulate(fit, re.form = NA), "re.form")
  expect_error(simulate(fit, censored = TRUE), "cens\\(\\)")
})
