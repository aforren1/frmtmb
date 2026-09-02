# Priors and bounds for sampling/fitting, and gr(cov=) known-covariance
# random effects.

sim_lmm <- function(seed = 301, n = 150, ng = 15) {
  set.seed(seed)
  dd <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng),
                                                length.out = n)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(ng, 0, 0.7)[dd$g], 1)
  dd
}

test_that("prior-augmented objective equals nll + neg log prior", {
  dd <- sim_lmm()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  pr <- list(x = prior_normal(0, 1), theta = prior_normal(0, 2))
  ri <- frmtmb:::resolve_prior_input(fit, pr)
  obj2 <- frmtmb:::prior_augmented_obj(fit, ri$entries)
  est <- fit$opt$par
  nlp_manual <-
    -stats::dnorm(fixef(fit)$mu[["x"]], 0, 1, log = TRUE) -
    stats::dnorm(fit$estimates$theta[1], 0, 2, log = TRUE)
  expect_lt(abs(obj2$fn(obj2$par) -
                  (-as.numeric(logLik(fit)) + nlp_manual)), 1e-8)
})

test_that("prior name resolution: classes, coefficients, errors", {
  dd <- sim_lmm()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  e <- frmtmb:::resolve_priors(fit, list(beta = prior_normal(0, 5)))
  expect_identical(e[[1]]$comp, "beta")
  expect_length(e[[1]]$idx, 2)
  expect_error(frmtmb:::resolve_priors(fit, list(zzz = prior_normal())),
               "Unknown parameter")
  expect_error(frmtmb:::resolve_priors(fit, list(x = 5)),
               "prior object")
})

test_that("a tight prior pulls the posterior toward it", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  dd <- sim_lmm(seed = 302)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- suppressWarnings(
    frm_sample(fit, chains = 1, iter = 600, refresh = 0, seed = 1,
               priors = list(x = prior_normal(0, 0.01)))
  )
  m <- as.matrix(ds)
  # shrunk to ~0: the posterior sd itself proves the prior bit, and the
  # mean is judged against that sd rather than a platform-fragile number
  expect_lt(stats::sd(m[, "x"]), 0.05)
  expect_lt(abs(mean(m[, "x"])), 5 * stats::sd(m[, "x"]) + 1e-8)
  expect_gt(fixef(fit)$mu[["x"]], 0.3)          # ML untouched
})

test_that("hard bounds constrain the ML fit", {
  dd <- sim_lmm(seed = 303)
  fit0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  est0 <- fixef(fit0)$mu[["x"]]
  expect_lt(est0, 1)
  fitb <- suppressWarnings(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
        lower = c(x = 1))
  )
  expect_equal(fixef(fitb)$mu[["x"]], 1, tolerance = 1e-6)
  expect_lt(as.numeric(logLik(fitb)), as.numeric(logLik(fit0)))
  expect_error(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
                   lower = c(zzz = 0)), "Unknown parameter")
})

test_that("gr(cov=) matches a hand-rolled correlated-intercepts reference", {
  set.seed(304)
  ng <- 30
  # a random PSD correlation matrix over the levels
  R <- crossprod(matrix(rnorm(ng * ng), ng)) / ng
  A <- stats::cov2cor(R)
  dimnames(A) <- list(as.character(seq_len(ng)),
                      as.character(seq_len(ng)))
  b_true <- drop(crossprod(chol(A), rnorm(ng))) * 0.8
  n <- 600
  g <- factor(rep(seq_len(ng), each = n / ng))
  x <- rnorm(n)
  dd <- data.frame(y = 1 + 0.5 * x + b_true[g] + rnorm(n, 0, 0.6),
                   x = x, g = g)

  fit <- frm(bf(y ~ x + (1 | gr(g, cov = A))) + gaussian(), data = dd)
  expect_identical(fit$frame$re_blocks[[1]]$covstruct, "gr_cov")

  yv <- dd$y; xv <- dd$x; gi <- as.integer(dd$g)
  Ad <- unname(A)
  nll_ref <- function(p) {
    Sigma <- exp(2 * p$lsd) * Ad
    nll <- -sum(RTMB::dmvnorm(p$u, 0, Sigma, log = TRUE))
    mu <- p$b[1] + p$b[2] * xv + p$u[gi]
    nll - sum(RTMB::dnorm(yv, mu, exp(p$ls), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(b = c(0, 0), ls = 0, lsd = 0,
                              u = numeric(ng)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)

  # sd recovered, methods work
  sd_hat <- sqrt(VarCorr(fit)[[1]][1, 1])
  expect_lt(abs(sd_hat - 0.8), 0.4)
  expect_identical(dim(ranef(fit)[[1]]), c(30L, 1L))
  cv <- confint_varcorr(fit)
  expect_true(all(cv$lwr < cv$estimate & cv$estimate < cv$upr))
  expect_equal(predict(fit, newdata = dd), predict(fit),
               tolerance = 1e-8)

  # validations
  A2 <- A; dimnames(A2) <- NULL
  expect_error(frm(bf(y ~ x + (1 | gr(g, cov = A2))) + gaussian(),
                   data = dd), "dimnames")
})