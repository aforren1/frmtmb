# Addition terms on newdata. predict(newdata =) used to read only
# trunc(), rate() and thres(gr = ) from newdata and to hand the
# simulator an empty list otherwise, and binomial's simulator read a
# missing trials() as one trial: with newdata N = 1000 every draw was 0
# or 1, and nothing was raised (dev/ndaterm-findings.md). simulate() and
# fitted() read every term, so they are the references here.
#
# Every tolerance is a multiple of a Monte Carlo standard error that the
# run measures from its own draws.

nd_binom <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(7)
      n <- 200
      d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
      d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
      cache <<- list(d = d,
                     fit = frm(bf(y | trials(N) ~ x), family = binomial(),
                               data = d),
                     nd = data.frame(x = c(0, 0), N = c(1000, 2)))
    }
    cache
  }
})

# |a - b| within k standard errors of the difference of two means, each
# standard error measured from its own draws
expect_means_agree <- function(a, b, sd_a, sd_b, n_a, n_b, k = 6) {
  se <- sqrt(sd_a^2 / n_a + sd_b^2 / n_b)
  expect_true(all(abs(a - b) <= k * se),
              info = paste("a =", paste(format(a), collapse = ", "),
                           "b =", paste(format(b), collapse = ", "),
                           "k * se =", paste(format(k * se), collapse = ", ")))
}

test_that("predict(newdata, summary = FALSE) draws trials() from newdata", {
  cs <- nd_binom()
  set.seed(11)
  p <- predict(cs$fit, newdata = cs$nd, summary = FALSE, ndraws = 400)
  # a count out of N: at most N, and more than one where N is 1000
  # (0.71.0: the maximum of both columns was 1)
  expect_true(all(p[, 1] <= 1000) && all(p[, 2] <= 2))
  expect_gt(max(p[, 1]), 2)
  s <- as.matrix(simulate(cs$fit, nsim = 400, newdata = cs$nd, seed = 11))
  expect_means_agree(colMeans(p), rowMeans(s), apply(p, 2, stats::sd),
                     apply(s, 1, stats::sd), nrow(p), ncol(s))
})

test_that("predict(newdata) summarizes draws out of newdata's trials()", {
  cs <- nd_binom()
  set.seed(12)
  ps <- predict(cs$fit, newdata = cs$nd, ndraws = 400)
  f <- fitted(cs$fit, newdata = cs$nd)
  # Est.Error is the sd of the same draws, so it is their standard error
  # times sqrt(400); the expected count's own sd adds nothing that
  # size of draw does not already carry
  expect_means_agree(ps[, "Estimate"], f[, "Estimate"], ps[, "Est.Error"],
                     0, 400, 1)
  expect_true(all(ps[, "Q97.5"] <= c(1000, 2)))
})

test_that("a model without trials() in newdata is refused, as brms does", {
  cs <- nd_binom()
  # brms 2.23.0: "The following variables can neither be found in 'data'
  # nor in 'data2'" (dev/ndaterm-brms.R); 0.71.0 drew 0/1 counts
  expect_error(predict(cs$fit, newdata = data.frame(x = 0), ndraws = 5),
               "trials(N) could not be evaluated on newdata", fixed = TRUE)
})

test_that("a binomial-type family refuses a missing trials() by name", {
  for (nm in c("binomial", "zero_inflated_binomial", "beta_binomial",
               "zero_inflated_beta_binomial")) {
    ff <- frmtmb:::family_registry[[nm]]()
    dp <- list(mu = rep(0.4, 3), zi = rep(0.1, 3), phi = rep(5, 3))
    dp <- dp[names(dp) %in% ff[["dpars"]]]
    expect_error(ff[["sim"]](dp, list(), 3), "trials(), did not reach",
                 fixed = TRUE)
    expect_error(ff[["post"]][["mean_fn"]](dp, list()),
                 "trials(), did not reach", fixed = TRUE)
    # with the term the same calls answer
    expect_length(ff[["sim"]](dp, list(trials = rep(10, 3)), 3), 3)
  }
})

test_that("multinomial predict(newdata) draws newdata's trials()", {
  set.seed(7)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
  d$Y <- t(vapply(seq_len(n), function(i) {
    stats::rmultinom(1, d$N[i], c(0.5, 0.3, 0.2))
  }, numeric(3)))
  fit <- frm(bf(Y | trials(N) ~ x), family = multinomial(K = 3), data = d)
  set.seed(13)
  # 0.71.0 refused: "a multinomial fit needs trials() to know how many
  # draws each row gets"
  p <- predict(fit, newdata = data.frame(x = 0, N = c(1000, 2)),
               summary = FALSE, ndraws = 20)
  expect_identical(unique(apply(p, c(1, 2), sum)[, 1]), 1000)
  expect_identical(unique(apply(p, c(1, 2), sum)[, 2]), 2)
})

nd_se <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(7)
      n <- 200
      d <- data.frame(x = stats::rnorm(n), s = stats::runif(n, 0.2, 1))
      # residual sd beyond the known one, so sigma = TRUE has a sigma
      d$y <- 1 + d$x + stats::rnorm(n, 0, sqrt(d$s^2 + 0.3^2))
      cache <<- list(d = d, nd = data.frame(x = c(0, 0), s = c(0.05, 5)))
    }
    cache
  }
})

# The sd of predictive draws at a known se: the known sd, the estimated
# sigma (sigma = TRUE) and the expected value's own spread, which
# fitted() reports as Est.Error, in quadrature
expect_se_draws <- function(fit, nd, sigma) {
  set.seed(14)
  p <- predict(fit, newdata = nd, summary = FALSE, ndraws = 1000)
  f <- fitted(fit, newdata = nd)
  want <- sqrt(nd$s^2 + sigma^2 + f[, "Est.Error"]^2)
  got <- apply(p, 2, stats::sd)
  # the sd of a normal sample has relative standard error
  # 1 / sqrt(2 (n - 1)); sigma's own draw spread widens it a little
  rse <- 1 / sqrt(2 * (nrow(p) - 1))
  expect_true(all(abs(got / want - 1) <= 6 * rse),
              info = paste("got", paste(format(got), collapse = ", "),
                           "want", paste(format(want), collapse = ", ")))
}

test_that("predict(newdata) draws at newdata's se()", {
  cs <- nd_se()
  fit <- frm(bf(y | se(s) ~ x), family = gaussian(), data = cs$d)
  # 0.71.0: sd 0.029 in both columns, the expected value's spread alone
  expect_se_draws(fit, cs$nd, 0)
  fit2 <- frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(),
              data = cs$d)
  expect_se_draws(fit2, cs$nd, sigma(fit2))
})

test_that("se() and cov = FALSE ARMA are refused at the fit, as in brms", {
  # so the ARMA fill on newdata (arma_cond_fill_dpars()) has no se() to
  # read; it takes every term all the same
  set.seed(7)
  d <- data.frame(t = 1:60, s = stats::runif(60, 0.2, 1),
                  y = stats::rnorm(60))
  expect_error(frm(bf(y | se(s) ~ 1 + ar(time = t, p = 1, cov = FALSE)),
                   family = gaussian(), data = d),
               "se() cannot be combined with cov = FALSE", fixed = TRUE)
})

cust_fams <- function() {
  list(
    cbin = custom_family(
      "nd_cbin", dpars = "mu", links = list(mu = "logit"),
      lpdf = function(y, dpars, aterms) {
        RTMB::dbinom(y, aterms[["vint1"]], dpars[["mu"]], log = TRUE)
      },
      required_aterms = "vint1", type = "discrete",
      post = list(mean_fn = function(dpars, aterms) {
        dpars[["mu"]] * aterms[["vint1"]]
      }),
      sim = function(dpars, aterms, n) {
        stats::rbinom(n, aterms[["vint1"]], dpars[["mu"]])
      }),
    cnorm = custom_family(
      "nd_cnorm", dpars = "mu", links = list(mu = "identity"),
      lpdf = function(y, dpars, aterms) {
        RTMB::dnorm(y, dpars[["mu"]], aterms[["vreal1"]], log = TRUE)
      },
      required_aterms = "vreal1",
      post = list(mean_fn = function(dpars, aterms) dpars[["mu"]]),
      sim = function(dpars, aterms, n) {
        stats::rnorm(n, dpars[["mu"]], aterms[["vreal1"]])
      }))
}

test_that("a custom family's vint() reaches predict(newdata)", {
  fams <- cust_fams()
  cs <- nd_binom()
  fc <- frm(bf(y | vint(N) ~ x), family = fams$cbin, data = cs$d)
  set.seed(15)
  # 0.71.0: "invalid arguments", from rbinom() given a NULL size
  p <- predict(fc, newdata = cs$nd, summary = FALSE, ndraws = 400)
  s <- as.matrix(simulate(fc, nsim = 400, newdata = cs$nd, seed = 15))
  expect_true(all(p[, 1] <= 1000) && all(p[, 2] <= 2))
  expect_means_agree(colMeans(p), rowMeans(s), apply(p, 2, stats::sd),
                     apply(s, 1, stats::sd), nrow(p), ncol(s))
})

test_that("a custom family's vreal() reaches predict(newdata)", {
  ce <- nd_se()
  fr <- frm(bf(y | vreal(s) ~ x), family = cust_fams()$cnorm, data = ce$d)
  # 0.71.0: "invalid arguments", from rnorm() given a NULL sd
  expect_se_draws(fr, ce$nd, 0)
})

test_that("a declared required term that does not arrive is named", {
  fams <- cust_fams()
  dp <- list(mu = rep(0.4, 3))
  expect_error(frmtmb:::sim_response(fams$cbin, dp, list(), 3),
               "nd_cbin: `vint1` did not reach the family", fixed = TRUE)
  expect_error(frmtmb:::response_mean(fams$cnorm, dp, list()),
               "nd_cnorm: `vreal1` did not reach the family", fixed = TRUE)
  # the same calls with the term answer
  expect_length(frmtmb:::sim_response(fams$cbin, dp,
                                      list(vint1 = rep(4, 3)), 3), 3)
})

test_that("nothing new fires where no term is due", {
  # a family that reads no addition term, and a model without one
  set.seed(7)
  d <- data.frame(x = stats::rnorm(80))
  d$y <- stats::rpois(80, exp(0.5 + 0.3 * d$x))
  d$b <- stats::rbinom(80, 1, 0.4)
  d$g <- 1 + d$x + stats::rnorm(80)
  nd <- data.frame(x = c(-1, 1))
  for (fit in list(frm(bf(y ~ x), family = poisson(), data = d),
                   frm(bf(b ~ x), family = bernoulli(), data = d),
                   frm(bf(g ~ x), family = gaussian(), data = d))) {
    set.seed(16)
    expect_no_condition(p <- predict(fit, newdata = nd, ndraws = 50))
    expect_identical(dim(p), c(2L, 4L))
    expect_no_condition(predict(fit, summary = FALSE, ndraws = 5))
  }
  # a binomial model's own rows read the frame's trials(), as before
  cs <- nd_binom()
  set.seed(17)
  p <- predict(cs$fit, summary = FALSE, ndraws = 50)
  expect_true(all(sweep(p, 2, cs$d$N, "<=")))
})

# --- class-wise draws (fix round 1, review finding 1) -------------------
#
# A mixture draws each component on the rows of its class. It handed the
# component the whole term list, so the j-th row of a class read the
# trials() of the j-th row of the data: 25% of in-sample draws exceeded
# N on 0.71.0. Seed 7, the reviewer's construction
# (dev/ndaterm-rev-probe2.R, probe3.R).

nd_mix <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(7)
      n <- 200
      d <- data.frame(x = stats::rnorm(n), N = sample(c(2, 50), n, TRUE))
      cls <- stats::rbinom(n, 1, 0.5)
      d$y <- stats::rbinom(n, d$N, ifelse(cls == 1,
                                          stats::plogis(-1 + 0.5 * d$x),
                                          stats::plogis(1 + 0.5 * d$x)))
      fit <- frm(bf(y | trials(N) ~ 1, mu1 ~ x, mu2 ~ x),
                 family = mixture(binomial(), binomial()), data = d)
      cache <<- list(d = d, fit = fit)
    }
    cache
  }
})

# draws (rows) by observations (columns), N the trials of each column:
# none above N, and the mean of the rows with N == 2 within 6 standard
# errors of fitted()
expect_class_draws <- function(s, N, f) {
  expect_identical(sum(sweep(s, 2, N, ">")), 0L)
  small <- s[, N == 2, drop = FALSE]
  m <- mean(small)
  se <- stats::sd(as.numeric(small)) / sqrt(length(small))
  expect_lt(abs(m - mean(f[N == 2])), 6 * se)
}

test_that("a mixture draws each class at its own rows' trials()", {
  cs <- nd_mix()
  s <- t(as.matrix(simulate(cs$fit, nsim = 200, seed = 1)))
  # 0.71.0: share above N 0.25, mean 11.77 on the N = 2 rows where
  # fitted() says 0.98
  expect_class_draws(s, cs$d$N, fitted(cs$fit)[, "Estimate"])
  set.seed(41)
  p <- predict(cs$fit, summary = FALSE, ndraws = 200)
  expect_class_draws(p, cs$d$N, fitted(cs$fit)[, "Estimate"])
})

test_that("a mixture's predict(newdata) matches fitted() at newdata's N", {
  cs <- nd_mix()
  nd <- data.frame(x = c(0, 0), N = c(1000, 2))
  set.seed(42)
  p <- predict(cs$fit, newdata = nd, summary = FALSE, ndraws = 2000)
  f <- fitted(cs$fit, newdata = nd)
  expect_true(all(p[, 1] <= 1000) && all(p[, 2] <= 2))
  # the expected value's own spread is the parameter part of the draws,
  # which fitted()'s Est.Error reports
  se <- apply(p, 2, stats::sd) / sqrt(nrow(p))
  expect_true(all(abs(colMeans(p) - f[, "Estimate"]) <= 6 * se),
              info = paste(format(colMeans(p)), format(f[, "Estimate"])))
})

test_that("a grouped mixture draws each class at its own rows' trials()", {
  set.seed(7)
  ng <- 40
  m <- 6
  d <- data.frame(id = rep(seq_len(ng), each = m), x = stats::rnorm(ng * m),
                  N = sample(c(2, 50), ng * m, TRUE))
  cl <- rep(stats::rbinom(ng, 1, 0.5), each = m)
  d$y <- stats::rbinom(nrow(d), d$N, ifelse(cl == 1, 0.2, 0.8))
  fit <- frm(bf(y | trials(N) ~ 1),
             family = mixture(binomial(), binomial(), groups = ~id),
             data = d)
  s <- t(as.matrix(simulate(fit, nsim = 100, seed = 1)))
  # 0.71.0: share above N 0.21
  expect_class_draws(s, d$N, fitted(fit)[, "Estimate"])
})

test_that("predict() hands a mixture's theta to the simulator unreported", {
  # theta reports as the softmax probability; the simulator read that
  # probability as the predictor and took the softmax again, so every
  # draw away from the in-sample fast path (newdata, re_formula) used
  # the wrong mixing weights
  set.seed(7)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  k <- stats::rbinom(n, 1, 0.2)
  d$y <- ifelse(k == 1, 6, 0) + d$x + stats::rnorm(n, 0, 0.5)
  fit <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()),
             data = d)
  nd <- data.frame(x = c(-1, 0, 1))
  set.seed(43)
  p <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000)
  f <- fitted(fit, newdata = nd)
  se <- apply(p, 2, stats::sd) / sqrt(nrow(p))
  expect_true(all(abs(colMeans(p) - f[, "Estimate"]) <= 6 * se),
              info = paste(format(colMeans(p)), format(f[, "Estimate"])))
})

test_that("a factor addition term is coded by the fit's levels", {
  # vint() of a factor is its level codes, which newdata's own levels
  # used to set
  set.seed(7)
  n <- 120
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(rep(c("lo", "hi"), n / 2), c("lo", "hi")))
  d$y <- 1 + d$x + 2 * (d$g == "hi") + stats::rnorm(n, 0, 0.5)
  shift <- custom_family(
    "nd_shift", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]] + 2 * (aterms[["vint1"]] - 1),
                  dpars[["sigma"]], log = TRUE)
    },
    required_aterms = "vint1",
    post = list(mean_fn = function(dpars, aterms) {
      dpars[["mu"]] + 2 * (aterms[["vint1"]] - 1)
    }))
  fit <- frm(bf(y | vint(g) ~ x), family = shift, data = d)
  ref <- fitted(fit, newdata = data.frame(
    x = 0, g = factor(c("lo", "hi"), c("lo", "hi"))))[, "Estimate"]
  rel <- fitted(fit, newdata = data.frame(
    x = 0, g = factor(c("lo", "hi"), c("hi", "lo"))))[, "Estimate"]
  one <- unname(fitted(fit, newdata = data.frame(x = 0, g = "hi"))[,
    "Estimate"])
  # 0.71.0: the relevelled rows swapped, and "hi" alone was code 1
  expect_equal(rel, ref)
  expect_equal(one, ref[2])
  expect_error(fitted(fit, newdata = data.frame(x = 0, g = "mid")),
               "holds value(s) 'mid' the fit did not see", fixed = TRUE)
})

test_that("conditional_effects() holds a vint() count at a whole value", {
  # held at its mean, 9.985, the count drew NA with a warning per draw
  # and the band failed; brms holds trials, thres and vint variables at
  # round(median()), 10 here
  cs <- nd_binom()
  fc <- frm(bf(y | vint(N) ~ x), family = cust_fams()$cbin, data = cs$d)
  set.seed(44)
  expect_message(
    ce <- conditional_effects(fc, "x", method = "predict", ndraws = 51),
    paste0("holding the addition-term variable(s) N at ",
           round(stats::median(cs$d$N))), fixed = TRUE)
  expect_true(all(is.finite(ce[[1]]$estimate__)))
  expect_true(all(ce[[1]]$N == round(stats::median(cs$d$N))))
})
