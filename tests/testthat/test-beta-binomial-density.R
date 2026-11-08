# The beta-binomial density at a precision with no limit in the data: raw
# lgamma() differences there cancelled to -66752 and +327680 at code 0, against a
# true -548.982078 that no phi moves. The reference is the finite product, which
# has no large lgamma() anywhere.

bb_exact_row <- function(y, size, p, phi) {
  a <- p * phi
  b <- (1 - p) * phi
  lchoose(size, y) +
    (if (y > 0) sum(log(a + 0:(y - 1))) else 0) +
    (if (size - y > 0) sum(log(b + 0:(size - y - 1))) else 0) -
    (if (size > 0) sum(log(phi + 0:(size - 1))) else 0)
}

bb_exact <- function(y, size, p, phi) {
  vapply(seq_along(y), function(i) bb_exact_row(y[i], size[i], p[i], phi), 0)
}

bb_data <- function(s, n = 300, trials = 10L) {
  set.seed(s)
  d <- data.frame(x = round(rnorm(n), 3), ntr = trials)
  d$y <- rbinom(n, trials, plogis(0.3 + 0.5 * d$x))
  d
}

test_that("lgamma_count_shift() is the difference, short and long", {
  a <- c(0.3, 2, 1e3, 1e8, exp(39), exp(100))
  for (s in c(0, 1, 7, 24, 25, 200)) {
    got <- lgamma_count_shift(a, rep(s, length(a)))
    ref <- vapply(a, function(ai) {
      if (s == 0) 0 else sum(log(ai + 0:(s - 1)))
    }, 0)
    # relative to the sum's own size, and the rounding of a sum of s logs
    expect_true(all(abs(got - ref) <=
                      1e-13 * pmax(1, abs(ref)) * max(1, sqrt(s))),
                info = paste("s =", s))
  }
  # a scalar shape beside per-row shifts
  expect_length(lgamma_count_shift(3, c(0L, 1L, 2L)), 3L)
  expect_equal(lgamma_count_shift(3, c(0L, 1L, 2L)),
               c(0, log(3), log(3) + log(4)), tolerance = 1e-14)
  # and against lgamma() where lgamma() has its digits
  x <- c(0.5, 3, 40)
  expect_equal(lgamma_count_shift(x, c(3L, 30L, 300L)),
               lgamma(x + c(3, 30, 300)) - lgamma(x), tolerance = 1e-12)
})

test_that("beta_binomial()'s density holds its value at every phi", {
  d <- bb_data(9)
  fam <- beta_binomial()
  p <- plogis(0.2 + 0.4 * d$x)
  for (lp in c(2, 10, 20, 39, 45, 60, 100)) {
    phi <- exp(lp)
    got <- fam$lpdf(d$y, list(mu = p, phi = rep(phi, nrow(d))),
                    list(trials = d$ntr))
    ex <- bb_exact(d$y, d$ntr, p, phi)
    expect_lt(max(abs(got - ex) / pmax(1, abs(ex))), 1e-12)
  }
})

test_that("zero_inflated_beta_binomial()'s density holds its value too", {
  d <- bb_data(9)
  fam <- zero_inflated_beta_binomial()
  p <- plogis(0.2 + 0.4 * d$x)
  zi <- 0.1
  for (lp in c(2, 20, 39, 45, 60, 100)) {
    phi <- exp(lp)
    got <- fam$lpdf(d$y, list(mu = p, phi = rep(phi, nrow(d)),
                              zi = rep(zi, nrow(d))),
                    list(trials = d$ntr))
    ex <- log1p(-zi) + bb_exact(d$y, d$ntr, p, phi)
    i0 <- d$y == 0
    ex[i0] <- log(zi + (1 - zi) * exp(ex[i0] - log1p(-zi)))
    expect_lt(max(abs(got - ex) / pmax(1, abs(ex))), 1e-12)
  }
})

test_that("the density is unchanged where the old form had its digits", {
  d <- bb_data(3)
  fam <- beta_binomial()
  p <- plogis(0.2 + 0.4 * d$x)
  for (phi in c(0.5, 8, 200)) {
    got <- fam$lpdf(d$y, list(mu = p, phi = rep(phi, nrow(d))),
                    list(trials = d$ntr))
    old <- RTMBdist::dbetabinom(d$y, d$ntr, p * phi, (1 - p) * phi, log = TRUE)
    expect_lt(max(abs(got - old) / abs(old)), 1e-12)
  }
})

test_that("data with no overdispersion reach the binomial limit at code 0", {
  for (s in c(2, 4, 6, 9, 12)) {
    d <- bb_data(s)
    bi <- as.numeric(logLik(glm(cbind(y, ntr - y) ~ x, binomial, d)))
    fit <- frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d)
    ll <- as.numeric(logLik(fit))
    expect_identical(fit$opt$convergence, 0L)
    # seed 2 read 0.15 above this limit, and seed 9 read +94208
    expect_lt(abs(ll - bi) / abs(bi), 1e-7)
    p <- plogis(fit$estimates$beta[[1L]] + fit$estimates$beta[[2L]] * d$x)
    ex <- sum(bb_exact(d$y, d$ntr, p, exp(fit$estimates$betad[[1L]])))
    expect_lt(abs(ll - ex) / abs(ex), 1e-10)
    body <- allow_warnings(suppressMessages(
      frm(bf(y | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
          family = beta_binomial(), data = d)),
      allowed = c("Standard errors", "not identified"))
    expect_lt(abs(as.numeric(logLik(body)) - bi) / abs(bi), 1e-7)
  }
})

test_that("zero_inflated_beta_binomial() reaches the same limit", {
  for (s in c(2, 9)) {
    d <- bb_data(s)
    bi <- as.numeric(logLik(glm(cbind(y, ntr - y) ~ x, binomial, d)))
    fit <- allow_warnings(
      frm(bf(y | trials(ntr) ~ x), family = zero_inflated_beta_binomial(),
          data = d),
      allowed = c("Standard errors", "not identified", "gradient",
                  "did not report convergence"))
    # at most the limit plus what a zero-inflation on these data can add
    expect_gt(as.numeric(logLik(fit)), bi - 1e-6 * abs(bi))
    expect_lt(as.numeric(logLik(fit)), bi + 2)
  }
})

test_that("overdispersed fits and many trials keep their value", {
  set.seed(11)
  n <- 400
  d <- data.frame(x = rnorm(n), ntr = 200L)
  p <- plogis(0.3 + 0.5 * d$x)
  d$y <- rbinom(n, d$ntr, rbeta(n, p * 8, (1 - p) * 8))
  fit <- frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d)
  expect_identical(fit$opt$convergence, 0L)
  pf <- plogis(fit$estimates$beta[[1L]] + fit$estimates$beta[[2L]] * d$x)
  ph <- exp(fit$estimates$betad[[1L]])
  old <- sum(RTMBdist::dbetabinom(d$y, d$ntr, pf * ph, (1 - pf) * ph,
                                  log = TRUE))
  expect_lt(abs(as.numeric(logLik(fit)) - old) / abs(old), 1e-11)
  ex <- sum(bb_exact(d$y, d$ntr, pf, ph))
  expect_lt(abs(as.numeric(logLik(fit)) - ex) / abs(ex), 1e-11)
})

test_that("the tape's gradient in log phi is flat at the limit", {
  d <- bb_data(9)
  p <- plogis(0.2 + 0.4 * d$x)
  tp <- RTMB::MakeTape(function(lp) {
    phi <- exp(lp)
    sum(beta_binomial_lpmf(d$y, d$ntr, p * phi, (1 - p) * phi, phi))
  }, 39)
  v39 <- tp(39)
  v60 <- tp(60)
  expect_lt(abs(v39 - v60) / abs(v39), 1e-12)
  # against the slope at a phi where the density does move
  expect_lt(abs(as.numeric(tp$jacobian(39))),
            1e-6 * abs(as.numeric(tp$jacobian(2))))
})

test_that("the censored and truncated forms hold their value at every phi", {
  d <- bb_data(9, n = 30)
  fam <- beta_binomial()
  p <- plogis(0.2 + 0.4 * d$x)
  for (lp in c(2, 20, 39, 45)) {
    phi <- exp(lp)
    pair <- fam$cdf_log(d$y, list(mu = p, phi = rep(phi, nrow(d))),
                        list(trials = d$ntr))
    ex_row <- function(i, set) {
      lv <- vapply(set, function(k) bb_exact_row(k, d$ntr[i], p[i], phi), 0)
      m <- max(lv)
      m + log(sum(exp(lv - m)))
    }
    exF <- vapply(seq_len(nrow(d)), function(i) ex_row(i, 0:d$y[i]), 0)
    tops <- which(d$y < d$ntr)
    exS <- vapply(tops, function(i) ex_row(i, (d$y[i] + 1):d$ntr[i]), 0)
    expect_lt(max(abs(pair$lF - exF) / pmax(1, abs(exF))), 1e-11)
    expect_lt(max(abs(pair$lS[tops] - exS) / pmax(1, abs(exS))), 1e-11)
  }
})

test_that("the pair holds its edges and per-row trials", {
  fam <- beta_binomial()
  size <- c(1L, 5L, 30L, 30L, 5L, 1L, 0L)
  q <- c(0, 2, 29, 40, -1, 1, 0)
  p <- c(0.2, 0.5, 0.7, 0.4, 0.9, 0.3, 0.5)
  for (phi in c(3, exp(39))) {
    pair <- fam$cdf_log(q, list(mu = p, phi = rep(phi, 7L)),
                        list(trials = size))
    for (i in which(q >= 0 & q < size)) {
      lv <- vapply(0:size[i], function(k) bb_exact_row(k, size[i], p[i], phi), 0)
      inF <- 0:size[i] <= q[i]
      expect_lt(abs(pair$lF[i] - log(sum(exp(lv[inF])))) /
                  max(1, abs(pair$lF[i])), 1e-11)
      expect_lt(abs(pair$lS[i] - log(sum(exp(lv[!inF])))) /
                  max(1, abs(pair$lS[i])), 1e-11)
    }
    # outside the support the pair is exactly 0 and -Inf, in the right order
    expect_identical(pair$lF[4L], 0)
    expect_identical(pair$lS[4L], -Inf)
    expect_identical(pair$lF[5L], -Inf)
    expect_identical(pair$lS[5L], 0)
  }
  # a censored fit against the sum of its own masses
  d <- bb_data(5, n = 200)
  d$cens <- ifelse(d$y < d$ntr & seq_len(nrow(d)) %% 4L == 0L, "right", "none")
  fit <- frm(bf(y | trials(ntr) + cens(cens) ~ x), family = beta_binomial(),
             data = d)
  pf <- plogis(fit$estimates$beta[[1L]] + fit$estimates$beta[[2L]] * d$x)
  ph <- exp(fit$estimates$betad[[1L]])
  ll <- vapply(seq_len(nrow(d)), function(i) {
    if (d$cens[i] == "none") {
      return(bb_exact_row(d$y[i], d$ntr[i], pf[i], ph))
    }
    log(sum(exp(vapply((d$y[i] + 1L):d$ntr[i], function(k) {
      bb_exact_row(k, d$ntr[i], pf[i], ph)
    }, 0))))
  }, 0)
  expect_lt(abs(as.numeric(logLik(fit)) - sum(ll)) / abs(sum(ll)), 1e-10)
})

test_that("a discrete fit above 0 warns in its own class", {
  set.seed(1)
  d <- data.frame(y = rpois(50, 3))
  # a mass that forgot a constant: a positive check is all that could see it
  bad <- frmtmb_family(
    "bad_pois", dpars = "mu", links = list(mu = "log"), type = "discrete",
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      y * log(dpars$mu) - dpars$mu - lgamma(y + 1) + 10
    },
    init_dpars = list(mu = function(y, aterms) mean(y)))
  expect_warning(fit <- frm(y ~ 1, family = bad, data = d),
                 class = "frmtmb_discrete_loglik")
  expect_gt(as.numeric(logLik(fit)), 0)
  # the correct mass
  expect_no_warning(frm(y ~ 1, family = poisson(), data = d))
  good <- frmtmb_family(
    "ok_cont", dpars = "mu", links = list(mu = "identity"),
    type = "continuous",
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      -0.5 * ((y - dpars$mu) / 1e-3)^2 - log(1e-3) - 0.5 * log(2 * pi)
    },
    init_dpars = list(mu = function(y, aterms) mean(y)))
  dc <- data.frame(y = 1 + 1e-4 * rnorm(50))
  # a density above 1 is no fault
  expect_no_warning(frm(y ~ 1, family = good, data = dc))
  expect_silent(discrete_loglik_warning(
    list(REML = TRUE, prior = NULL, opt = list(objective = -5),
         frame = list(n_obs = 5L),
         spec = list(responses = list(list(family = list(type = "discrete")))))))
  # under a prior the check runs on logLik(): test-loglik-prior.R
  # a sum of near-zero log masses that passes 0 by its rounding
  expect_silent(discrete_loglik_warning(
    list(REML = FALSE, prior = NULL, opt = list(objective = -1e-12),
         frame = list(n_obs = 5L),
         spec = list(responses = list(list(family = list(type = "discrete")))))))
})

test_that("simulate() and the OSA refusal keep their routes", {
  d <- bb_data(9)
  fit <- frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d)
  sims <- simulate(fit, nsim = 2, seed = 1)
  expect_true(all(unlist(sims) >= 0 & unlist(sims) <= 10))
  # the mean of a draw at the fit's mean, within a few of its own se
  v <- unlist(sims, use.names = FALSE)
  mu <- rep(10 * plogis(fit$estimates$beta[[1L]] +
                          fit$estimates$beta[[2L]] * d$x), 2L)
  expect_lt(abs(mean(v - mu)) / (sd(v - mu) / sqrt(length(v))), 5)
  expect_error(beta_binomial_lpmf(structure(1, class = "osa"), 2, 1, 1, 2),
               "does not support OSA", fixed = TRUE)
  # RTMBdist 1.2.0's own route reads NaN with no warning, so the refusal must
  # hold on the real route and not depend on the installed RTMBdist
  expect_error(residuals(fit, type = "osa"), "does not support OSA",
               fixed = TRUE)
})
