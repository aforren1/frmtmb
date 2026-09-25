# simulate(newdata = ) and simulate()'s re_formula, read as predict()
# reads it (lane wt-simnewdata). The behaviors that change are each in a
# block of their own, so that on the base build an error in one cannot
# end a block before the assertions behind it run.

sn_data <- function(seed = 7, n_g = 12, n_per = 10, sd_g = 1, sigma = 0.5) {
  set.seed(seed)
  n <- n_g * n_per
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(rep(seq_len(n_g), each = n_per)),
                  h = factor(rep(seq_len(5), length.out = n)))
  d$y <- 0.5 + 0.8 * d$x + stats::rnorm(n_g, 0, sd_g)[d$g] +
    stats::rnorm(n, 0, sigma)
  d$cnt <- stats::rpois(n, exp(0.2 + 0.3 * d$x +
                                 stats::rnorm(n_g, 0, 0.4)[d$g]))
  d
}

sn_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- frm(bf(y ~ x + (1 | g)) + gaussian(), data = sn_data())
    }
    cached
  }
})

# the mean correlation, across replicates, of the row pairs in `pairs`
sn_pair_cor <- function(sims, pairs) {
  m <- as.matrix(sims)
  mean(vapply(seq_len(nrow(pairs)), function(k) {
    stats::cor(m[pairs[k, 1], ], m[pairs[k, 2], ])
  }, 0))
}

test_that("simulate(newdata = the fitted rows) draws what simulate() draws", {
  d <- sn_data()
  fit <- sn_fit()
  for (rf in list(NULL, NA)) {
    a <- simulate(fit, nsim = 4, seed = 11, re_formula = rf)
    b <- simulate(fit, nsim = 4, seed = 11, re_formula = rf, newdata = d)
    expect_equal(unname(as.matrix(b)), unname(as.matrix(a)),
                 info = format(rf))
  }
})

test_that("the newdata draw carries trials(), cs() and a count family", {
  set.seed(5)
  n <- 150
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:10, 15)))
  d$nt <- sample(5:15, n, replace = TRUE)
  d$yb <- stats::rbinom(n, d$nt, stats::plogis(0.2 + 0.6 * d$x))
  d$yo <- factor(cut(d$x + stats::rlogis(n), c(-Inf, -0.5, 0.5, Inf),
                     labels = FALSE), ordered = TRUE)
  fb <- frm(bf(yb | trials(nt) ~ x + (1 | g)) + binomial(), data = d)
  fo <- frm(bf(yo ~ cs(x)) + sratio(), data = d)
  for (f in list(fb, fo)) {
    a <- simulate(f, nsim = 3, seed = 2)
    b <- simulate(f, nsim = 3, seed = 2, newdata = d)
    expect_identical(lapply(b, as.integer), lapply(a, as.integer))
  }
  # the trials column is read from newdata, so it bounds the draws there
  nd <- d[1:20, ]
  nd$nt <- 1L
  s <- simulate(fb, nsim = 20, seed = 3, newdata = nd)
  expect_true(all(as.matrix(s) %in% c(0, 1)))
})

test_that("a newdata row that is new draws at its own covariates", {
  fit <- sn_fit()
  nd <- data.frame(x = c(-2, 0, 2), g = factor(c("1", "2", "3")))
  s <- simulate(fit, nsim = 2000, seed = 4, newdata = nd)
  mu <- as.vector(frm_linpred(fit, newdata = nd))
  se <- sigma(fit) / sqrt(2000)
  expect_true(all(abs(rowMeans(as.matrix(s)) - mu) < 5 * se))
  expect_identical(dim(s), c(3L, 2000L))
})

test_that("re_formula = NA at newdata shares one group draw across a level", {
  # 2000 replicates: the correlation below has a standard error of about
  # (1 - rho^2) / sqrt(2000), so zero and the implied value are many
  # standard errors apart. The feasibility probe read 0.027 against an
  # implied 0.116 at 200 replicates, which was noise.
  fit <- sn_fit()
  nd <- data.frame(x = c(0, 0.5, 0, 0.5), g = factor(c(1, 1, 2, 2)))
  s <- simulate(fit, nsim = 2000, seed = 8, re_formula = NA, newdata = nd)
  sg2 <- varcorr_matrices(fit)[[1]][1, 1]
  rho <- sg2 / (sg2 + sigma(fit)^2)
  se <- (1 - rho^2) / sqrt(2000)
  within <- sn_pair_cor(s, rbind(c(1, 2), c(3, 4)))
  across <- sn_pair_cor(s, rbind(c(1, 3), c(2, 4)))
  expect_lt(abs(within - rho), 5 * se)
  expect_lt(abs(across), 5 / sqrt(2000))
  # conditioning instead holds each level at its estimate: no shared draw
  s0 <- simulate(fit, nsim = 2000, seed = 8, newdata = nd)
  expect_lt(abs(sn_pair_cor(s0, rbind(c(1, 2), c(3, 4)))), 5 / sqrt(2000))
})

test_that("re_formula = ~1 and ~0 redraw the group effects, as NA does", {
  fit <- sn_fit()
  na <- simulate(fit, nsim = 3, seed = 21, re_formula = NA)
  expect_identical(simulate(fit, nsim = 3, seed = 21, re_formula = ~1), na)
  expect_identical(simulate(fit, nsim = 3, seed = 21, re_formula = ~0), na)
  expect_false(identical(simulate(fit, nsim = 3, seed = 21), na))
})

test_that("a formula naming every term conditions on them, as NULL does", {
  fit <- sn_fit()
  expect_identical(simulate(fit, nsim = 3, seed = 22,
                            re_formula = ~ (1 | g)),
                   simulate(fit, nsim = 3, seed = 22))
})

test_that("a partial re_formula keeps the named term and redraws the other", {
  d <- sn_data(seed = 9, n_g = 10, n_per = 20, sd_g = 1.2, sigma = 0.4)
  d$y <- d$y + stats::rnorm(5, 0, 1.2)[d$h]
  fit <- frm(bf(y ~ x + (1 | g) + (1 | h)), data = d)
  s <- simulate(fit, nsim = 2000, seed = 3, re_formula = ~ (1 | g))
  vc <- varcorr_matrices(fit)
  sh2 <- vc[["1 | h"]][1, 1]
  rho_h <- sh2 / (sh2 + sigma(fit)^2)
  se <- (1 - rho_h^2) / sqrt(2000)
  # rows 1 and 6 share h and differ in g; rows 1 and 2 share g and
  # differ in h
  expect_identical(as.integer(d$h[c(1, 6)]), c(1L, 1L))
  expect_identical(as.integer(d$g[c(1, 2)]), c(1L, 1L))
  expect_lt(abs(sn_pair_cor(s, rbind(c(1, 6))) - rho_h), 5 * se)
  expect_lt(abs(sn_pair_cor(s, rbind(c(1, 2)))), 5 / sqrt(2000))
  # and the kept term sits at its estimate: the row mean is the
  # prediction with g kept and h dropped
  mu <- as.vector(frm_linpred(fit, re_formula = ~ (1 | g)))
  sd_row <- sqrt(sh2 + sigma(fit)^2)
  expect_lt(max(abs(rowMeans(as.matrix(s)) - mu)), 5 * sd_row / sqrt(2000))
})

test_that("a dropped slope is drawn given the kept intercept", {
  # End to end, on the fitted rows so that the base build answers too:
  # three rows of one level under re_formula = ~ (1 | g) on a
  # (1 + x1 + x2 | g) fit. The kept intercept is at its estimate, and
  # the two slopes are drawn from their law GIVEN it, so the rows' mean
  # and covariance are those of the exact conditional Gaussian. Before
  # the change a partial formula conditioned on every column, which
  # gives covariance sigma^2 I; drawing the slopes from their MARGINAL
  # law would give a mean and covariance many standard errors away too,
  # and the last two assertions check that this instrument can tell.
  set.seed(31)
  n_g <- 30
  d <- data.frame(x1 = stats::rnorm(n_g * 25), x2 = stats::rnorm(n_g * 25),
                  g = factor(rep(seq_len(n_g), 25)))
  Vt <- matrix(c(1, .6, -.5, .6, .8, -.2, -.5, -.2, .6), 3)
  U <- matrix(stats::rnorm(n_g * 3), n_g) %*% chol(Vt)
  d$y <- 1 + d$x1 - d$x2 + U[d$g, 1] + U[d$g, 2] * d$x1 +
    U[d$g, 3] * d$x2 + stats::rnorm(nrow(d), 0, 0.5)
  fit <- frm(bf(y ~ x1 + x2 + (1 + x1 + x2 | g)), data = d)
  V <- varcorr_matrices(fit)[[1]]
  rows <- which(d$g == "4")[1:3]
  b4 <- as.numeric(ranef(fit)[[1]]["4", ])
  Z <- cbind(1, d$x1[rows], d$x2[rows])
  A <- V[2:3, 1, drop = FALSE] / V[1, 1]
  m <- as.vector(A * b4[1])
  C <- V[2:3, 2:3] - A %*% V[1, 2:3, drop = FALSE]
  fx <- as.vector(Z %*% fixef(fit, flatten = TRUE)[c("(Intercept)", "x1",
                                                      "x2")])
  mu <- fx + Z[, 1] * b4[1] + as.vector(Z[, 2:3] %*% m)
  S <- Z[, 2:3] %*% C %*% t(Z[, 2:3]) + diag(sigma(fit)^2, 3)
  R <- 4000
  s <- as.matrix(simulate(fit, nsim = R, seed = 7,
                          re_formula = ~ (1 | g)))[rows, ]
  E <- stats::cov(t(s))
  # standard errors of a sample mean and of a sample covariance under
  # normality, from the law being tested
  se_mu <- sqrt(diag(S) / R)
  se_S <- sqrt((outer(diag(S), diag(S)) + S^2) / R)
  expect_lt(max(abs((rowMeans(s) - mu) / se_mu)), 5)
  expect_lt(max(abs((E - S) / se_S)), 5)
  mu_marg <- fx + Z[, 1] * b4[1]
  S_marg <- Z[, 2:3] %*% V[2:3, 2:3] %*% t(Z[, 2:3]) +
    diag(sigma(fit)^2, 3)
  expect_gt(max(abs((rowMeans(s) - mu_marg) / se_mu)), 5)
  expect_gt(max(abs((E - S_marg) / se_S)), 5)
})

test_that("simulate() refuses a re_formula term the fit does not have", {
  expect_error(simulate(sn_fit(), re_formula = ~ (1 | nosuch)),
               "(1 | nosuch)", fixed = TRUE, class = "frmtmb_error")
})

test_that("simulate() refuses a re_formula that is not NULL, NA or a formula", {
  expect_error(simulate(sn_fit(), re_formula = "g"), "re_formula",
               class = "frmtmb_error")
})

test_that("simulate() refuses a partial re_formula beside a factor smooth", {
  # two named terms, so that ~ (1 | g) keeps SOME of them; with one
  # named term it would keep all of them and mean NULL, as in predict()
  set.seed(4)
  d <- data.frame(x = stats::runif(200), f = factor(rep(1:4, 50)),
                  g = factor(rep(1:10, 20)), h = factor(rep(1:8, 25)))
  d$y <- sin(3 * d$x) + stats::rnorm(4, 0, 0.3)[d$f] +
    stats::rnorm(10, 0, 0.5)[d$g] + stats::rnorm(8, 0, 0.5)[d$h] +
    stats::rnorm(200, 0, 0.3)
  fit <- suppressWarnings(frm(bf(y ~ s(x, f, bs = "fs", k = 5) + (1 | g) +
                                   (1 | h)), data = d))
  expect_error(simulate(fit, re_formula = ~ (1 | g)), "cannot name",
               class = "frmtmb_error")
})

test_that("re_formula = NA keeps a population smooth curve", {
  # 0.62.0 redrew the smooth's penalized coefficients from the smoothing
  # prior: the draws spread with an sd of 2.13 around a curve whose
  # residual sd is 0.28 (dev/simnewdata-log/probe2-base.txt)
  set.seed(21)
  d <- data.frame(x = stats::runif(200))
  d$y <- 2 * sin(2 * pi * d$x) + stats::rnorm(200, 0, 0.3)
  fit <- frm(bf(y ~ s(x)), data = d)
  expect_identical(simulate(fit, nsim = 5, seed = 1, re_formula = NA),
                   simulate(fit, nsim = 5, seed = 1))
})

test_that("an unseen level needs allow_new_levels under a kept term", {
  fit <- sn_fit()
  nd <- data.frame(x = 0, g = factor("new"))
  expect_error(simulate(fit, newdata = nd), "New levels",
               class = "frmtmb_error")
})

test_that("allow_new_levels draws an unseen level's effect", {
  fit <- sn_fit()
  nd <- data.frame(x = c(0, 0), g = factor(c("new", "new")))
  s <- simulate(fit, nsim = 2000, seed = 5, newdata = nd,
                allow_new_levels = TRUE)
  sg2 <- varcorr_matrices(fit)[[1]][1, 1]
  v <- sg2 + sigma(fit)^2
  # both rows share the one unseen level's draw
  rho <- sg2 / v
  expect_lt(abs(sn_pair_cor(s, rbind(c(1, 2))) - rho),
            5 * (1 - rho^2) / sqrt(2000))
  expect_lt(abs(stats::var(unlist(s[1, ])) / v - 1), 5 * sqrt(2 / 1999))
})

test_that("an unseen level is one more fresh level under re_formula = NA", {
  fit <- sn_fit()
  nd <- data.frame(x = 0, g = factor("new"))
  s <- simulate(fit, nsim = 5, seed = 6, re_formula = NA, newdata = nd)
  expect_identical(dim(s), c(1L, 5L))
})

test_that("a redrawn term needs its grouping column in newdata", {
  expect_error(simulate(sn_fit(), re_formula = NA,
                        newdata = data.frame(x = 0)),
               "no column `g`", class = "frmtmb_error")
})

test_that("simulate(newdata = ) refuses a group-level mixture", {
  set.seed(4)
  dg <- data.frame(g = factor(rep(seq_len(30), each = 4)))
  cls <- rep(c(1, 2), each = 15)
  dg$y <- stats::rnorm(nrow(dg), c(0, 4)[cls[as.integer(dg$g)]], 1)
  fg <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
            data = dg)
  expect_error(simulate(fg, newdata = dg[1:8, ]), "indexes the rows",
               class = "frmtmb_error")
})

test_that("a residual correlation is rebuilt on the newdata rows", {
  set.seed(8)
  n_g <- 40
  tt <- 1:6
  d <- expand.grid(time = tt, g = factor(seq_len(n_g)))
  e <- unlist(lapply(seq_len(n_g), function(i) {
    as.vector(stats::arima.sim(list(ar = 0.7), length(tt)))
  }))
  d$y <- 1 + e
  fit <- frm(bf(y ~ 1 + ar(time, gr = g, cov = TRUE)), data = d)
  R <- frmtmb:::autocor_matrix(fit)
  nd <- data.frame(time = c(2, 3, 2), g = factor(c("a", "a", "b")))
  s <- simulate(fit, nsim = 2000, seed = 9, newdata = nd)
  r12 <- R[2, 3]
  expect_lt(abs(sn_pair_cor(s, rbind(c(1, 2))) - r12),
            5 * (1 - r12^2) / sqrt(2000))
  expect_lt(abs(sn_pair_cor(s, rbind(c(1, 3)))), 5 / sqrt(2000))
})

test_that("a residual correlation refuses a time the fit never saw", {
  set.seed(8)
  d <- expand.grid(time = 1:5, g = factor(1:20))
  d$y <- stats::rnorm(nrow(d))
  fit <- frm(bf(y ~ 1 + ar(time, gr = g, cov = TRUE)), data = d)
  expect_error(simulate(fit, newdata = data.frame(time = 9, g = "a")),
               "outside that set", class = "frmtmb_error")
})

test_that("censored = TRUE applies the fitted window at newdata", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(200))
  ylat <- 1 + d$x + stats::rnorm(200)
  d$cen <- as.integer(ylat > 2)
  d$y <- pmin(ylat, 2)
  fit <- frm(bf(y | cens(cen) ~ x), data = d)
  nd <- data.frame(x = c(2, 3), cen = 0L)
  s <- simulate(fit, nsim = 200, seed = 1, newdata = nd, censored = TRUE)
  expect_true(all(as.matrix(s) <= 2))
  expect_true(any(as.matrix(s) == 2))
})

test_that("simulate() refuses a newdata that is not a data frame", {
  expect_error(simulate(sn_fit(), newdata = list(x = 1)), "data frame",
               class = "frmtmb_error")
})

test_that("frm_bootstrap() redraws a smooth by default, and NULL holds it", {
  # The whole-model bootstrap of 0.62.0, kept by the user's decision of
  # 2026-09-24 although simulate(re_formula = NA) now holds a population
  # smooth: every replicate redraws the smooth's penalized coefficients,
  # so the refitted curve at one x spreads far more than the conditional
  # bootstrap's (re_formula = NULL), which redraws the noise alone.
  # Measured 1.80 against 0.054 at 20 replicates on this design
  # (dev/simnewdata-log/boot-test-*.txt); the pre-fix lane build gave
  # the two the same spread.
  set.seed(21)
  d <- data.frame(x = stats::runif(200))
  d$y <- 2 * sin(2 * pi * d$x) + stats::rnorm(200, 0, 0.3)
  fit <- frm(bf(y ~ s(x)), data = d)
  nd <- data.frame(x = 0.25)
  FUN <- function(f) frm_linpred(f, newdata = nd)
  whole <- frm_bootstrap(fit, FUN = FUN, nsim = 20, seed = 1)
  cond <- frm_bootstrap(fit, FUN = FUN, nsim = 20, seed = 1,
                        re_formula = NULL)
  expect_gt(stats::sd(whole$t) / stats::sd(cond$t), 5)
})

test_that("a newdata row with a missing covariate draws NA, as predict()", {
  # predict() leaves such a row NA; the family's generator gave NaN and
  # a base "NAs produced" warning per replicate
  nd <- data.frame(x = c(NA, 1), g = factor(c("1", "2")))
  expect_warning(s <- simulate(sn_fit(), nsim = 3, seed = 1, newdata = nd),
                 "not finite at the estimates", class = "frmtmb_warning")
  expect_true(all(is.na(unlist(s[1, ]))))
  expect_false(anyNA(unlist(s[2, ])))
})

test_that("a grouping factor read from the environment must fit newdata", {
  # predict() reads a grouping variable that is not a column of newdata
  # from the formula's environment, as model.frame() does, and dies on
  # "non-conformable arrays" when its length is wrong
  # (dev/simnewdata-log/edge-predict-base.txt); simulate() names it.
  # Fitted here so that this block's environment is the formula's.
  fit <- frm(bf(y ~ x + (1 | g)), data = sn_data())
  g <- factor(rep("3", 5))
  expect_error(simulate(fit, nsim = 1, seed = 1, re_formula = NA,
                        newdata = data.frame(x = c(0, 1))),
               "has 5 value(s) for newdata's 2 rows", fixed = TRUE,
               class = "frmtmb_error")
})
